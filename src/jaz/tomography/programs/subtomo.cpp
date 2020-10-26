#include "subtomo.h"
#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/projection/point_insertion.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/time.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <iostream>

using namespace gravis;


void SubtomoProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Catalogue .tbl or .star file");
		tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		boxSize = textToInteger(parser.getOption("--b", "Binned projection box size", "100"));
		cropSize = textToInteger(parser.getOption("--crop", "Output box size", "-1"));
		binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));
		write_multiplicity = parser.checkOption("--multi", "Write out multiplicity volumes");
		SNR = textToDouble(parser.getOption("--SNR", "Assumed signal-to-noise ratio (negative means use a heuristic)", "-1"));
		
		do_cone_weight = parser.checkOption("--cone_weight", "Weight down a double cone along Z");
		const double alpha = 0.5 * textToDouble(parser.getOption("--cone_angle", "Opening angle of the cone in degrees", "10"));
		cone_slope = sin(DEG2RAD(alpha));
		cone_sig0 = textToDouble(parser.getOption("--cone_sig0", "Cone width at Z = 0", "2"));	

		do_circle_crop = !parser.checkOption("--no_circle_crop", "Do not crop 2D images to a circle");
		do_narrow_circle_crop = parser.checkOption("--narrow_circle_crop", "Crop CTF-modulated 2D images to a circle that fits the (smaller) output box size");
		do_gridding_precorrection = parser.checkOption("--grid_precorr", "Perform gridding pre-correction on 2D images");

		taper = textToDouble(parser.getOption("--taper", "Taper against the sphere by this number of pixels", "5"));
		env_sigma = textToDouble(parser.getOption("--env", "Sigma of a Gaussian envelope applied before cropping", "-1"));

		do_whiten = parser.checkOption("--whiten", "Whiten the noise by flattening the power spectrum");
		do_center = !parser.checkOption("--no_center", "Do not subtract the mean from the voxel values");

		flip_value = !parser.checkOption("--no_ic", "Do not invert contrast (keep particles dark)");
		write_combined = !parser.checkOption("--no_comb", "Do not write the concatenated CTF-multiplicity image");
		write_ctf = parser.checkOption("--ctf", "Write 3D CTFs");
		write_divided = parser.checkOption("--div", "Write CTF-corrected subtomograms");
		write_normalised = parser.checkOption("--nrm", "Write multiplicity-normalised subtomograms");
		
		motFn = parser.getOption("--mot", "Particle trajectories", "");
		
		diag = parser.checkOption("--diag", "Write out diagnostic information");
		do_sum_all = parser.checkOption("--sum", "Sum up all subtomograms (for debugging)");

		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		outTag = parser.getOption("--o", "Output filename pattern");
		
		Log::readParams(parser);

		parser.checkForErrors();

	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	if (do_gridding_precorrection)
	{
		do_narrow_circle_crop = true;
	}

	ZIO::makeOutputDir(outTag);

	{
		std::ofstream ofs(outTag+"_note.txt");

		ofs << "Command:\n\n";

		for (int i = 0; i < argc; i++)
		{
			ofs << argv[i] << ' ';
		}

		ofs << '\n';
	}
}

void SubtomoProgram::run()
{
	TomogramSet tomogramSet(tomoSetFn);

	ParticleSet dataSet(particlesFn, motFn);
	std::vector<std::vector<ParticleIndex> > particles = dataSet.splitByTomogram(tomogramSet);
	
	if (cropSize < 0) cropSize = boxSize;
	
	bool do_ctf = true;
	
	const long int tc = particles.size();
	const long int s2D = boxSize;
	const long int sh2D = s2D / 2 + 1;
	
	const long int s3D = cropSize;
	const long int sh3D = s3D / 2 + 1;
	
	const long int s02D = (int)(binning * s2D + 0.5);
	
	const double relative_box_scale = cropSize / (double) boxSize;

	ParticleSet copy = dataSet;

	for (int t = 0; t < tc; t++)
	{
		const int pc = particles[t].size();

		if (pc == 0) continue;

		for (int p = 0; p < pc; p++)
		{
			const ParticleIndex part_id = particles[t][p];

			const int opticsGroup = copy.getOpticsGroup(part_id);
			const double originalPixelSize = copy.getOriginalPixelSize(opticsGroup);

			std::string outData = outTag + "/" + dataSet.getName(part_id) + "_data.mrc";
			std::string outWeight = outTag + "/" + dataSet.getName(part_id) + "_weights.mrc";

			copy.setImageFileNames(outData, outWeight, part_id);
			
			const d3Vector offset_A = copy.getParticleOffset(part_id);
			const d3Vector coord_0 = copy.getParticleCoord(part_id);
			
			const d3Vector coord_1 = coord_0 - offset_A / originalPixelSize;
			
			copy.setParticleOffset(part_id, d3Vector(0,0,0));
			copy.setParticleCoord(part_id, coord_1);
		}
	}

	for (int og = 0; og < copy.numberOfOpticsGroups(); og++)
	{
		const double ps_img = copy.optTable.getDouble(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, og);
		const double ps_out = binning * ps_img;

		copy.optTable.setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, 1, og);
		copy.optTable.setValue(EMDL_IMAGE_DIMENSIONALITY, 3, og);
		copy.optTable.setValue(EMDL_TOMO_SUBTOMOGRAM_BINNING, binning, og);
		copy.optTable.setValue(EMDL_IMAGE_PIXEL_SIZE, ps_out, og);
		copy.optTable.setValue(EMDL_IMAGE_SIZE, cropSize, og);
	}

	copy.write(outTag + "_particles.star");


	BufferedImage<float> sum_data, sum_weights;

	if (do_sum_all)
	{
		sum_data.resize(s3D,s3D,s3D);
		sum_data.fill(0.0);

		sum_weights.resize(sh3D,s3D,s3D);
		sum_weights.fill(0.0);
	}

	AberrationsCache aberrationsCache(dataSet.optTable, s2D);

	
	for (int t = 0; t < tc; t++)
	{
		const int pc = particles[t].size();
		if (pc == 0) continue;
		
		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));
		Log::print("Loading");
				
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		
		const int fc = tomogram.frameCount;
		
		dataSet.checkTrajectoryLengths(particles[t][0], pc, fc, "subtomo");
		
		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s2D, binning);
		BufferedImage<float> noiseWeights;

		if (do_whiten)
		{
			noiseWeights = tomogram.computeNoiseWeight(s2D, binning);
		}

		const int inner_thread_num = 1;
		const int outer_thread_num = num_threads / inner_thread_num;
		
		// @TODO: define input and output pixel sizes!
		
		const double binnedPixelSize = tomogram.optics.pixelSize * binning;
		
		Log::beginProgress("Backprojecting subtomograms", (int)ceil(pc/(double)outer_thread_num));
		
		#pragma omp parallel for num_threads(outer_thread_num)	
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();
			
			if (th == 0)
			{
				Log::updateProgress(p);
			}
						
			const ParticleIndex part_id = particles[t][p];
			
			const d3Vector pos = dataSet.getPosition(part_id);
			const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(
						part_id, fc, tomogram.optics.pixelSize);
			
			std::vector<d4Matrix> projCut(fc), projPart(fc);
			
			BufferedImage<fComplex> particleStack = BufferedImage<fComplex>(sh2D,s2D,fc);
			BufferedImage<float> weightStack(sh2D,s2D,fc);

			const bool circleCrop = false;
			
			TomoExtraction::extractAt3D_Fourier(
					tomogram.stack, s02D, binning, tomogram.projectionMatrices, traj,
					particleStack, projCut, inner_thread_num, circleCrop);
			
			
			if (!do_ctf) weightStack.fill(1.f);


			const int og = dataSet.getOpticsGroup(part_id);

			const BufferedImage<double>* gammaOffset =
				aberrationsCache.hasSymmetrical? &aberrationsCache.symmetrical[og] : 0;
			
			for (int f = 0; f < fc; f++)
			{
				projPart[f] = projCut[f] * d4Matrix(dataSet.getSubtomogramMatrix(part_id));
				
				if (do_ctf)
				{
					CTF ctf = tomogram.getCtf(f, pos);
					BufferedImage<float> ctfImg(sh2D, s2D);
					ctf.draw(s2D, s2D, binnedPixelSize, gammaOffset, &ctfImg(0,0,0));
					
					const float sign = flip_value? -1.f : 1.f;

					for (int y = 0; y < s2D;  y++)
					for (int x = 0; x < sh2D; x++)
					{
						const double c = ctfImg(x,y) * doseWeights(x,y,f);
						
						particleStack(x,y,f) *= sign * c;
						weightStack(x,y,f) = c * c;
					}

					if (aberrationsCache.hasAntisymmetrical)
					{
						if (aberrationsCache.phaseShift[og].ydim != s2D)
						{
							REPORT_ERROR_STR(
								"subtomo: wrong cached phase-shift size. Box size: "
								<< s2D << ", cache size: " << aberrationsCache.phaseShift[og].ydim);
						}

						for (int y = 0; y < s2D;  y++)
						for (int x = 0; x < sh2D; x++)
						{
							const fComplex r = aberrationsCache.phaseShift[og](x,y);
							const fComplex z = particleStack(x,y,f);

							particleStack(x,y,f).real = z.real * r.real + z.imag * r.imag;
							particleStack(x,y,f).imag = z.imag * r.real - z.real * r.imag;
						}
					}
				}
			}

			if (do_whiten)
			{
				particleStack *= noiseWeights;
				weightStack *= noiseWeights;
			}
						
			const int boundary = (boxSize - cropSize) / 2;

			if (do_gridding_precorrection || do_circle_crop)
			{
				BufferedImage<float> particlesRS;
				
				particlesRS = NewStackHelper::inverseFourierTransformStack(particleStack);

				if (do_circle_crop)
				{
					const double crop_boundary = do_narrow_circle_crop? boundary : 0.0;
					TomoExtraction::cropCircle(particlesRS, crop_boundary, 5, num_threads);
				}

				if (do_gridding_precorrection)
				{
					TomoExtraction::griddingPreCorrect(particlesRS, boundary, num_threads);
				}
				
				particleStack = NewStackHelper::FourierTransformStack(particlesRS);
			}
			
			BufferedImage<fComplex> dataImgFS(sh3D,s3D,s3D);
			dataImgFS.fill(fComplex(0.0, 0.0));
			
			BufferedImage<float> ctfImgFS(sh3D,s3D,s3D), psfImgFS(sh3D,s3D,s3D), 
					dataImgRS(s3D,s3D,s3D), dataImgDivRS(s3D,s3D,s3D),
					multiImageFS(sh3D,s3D,s3D);
			
			ctfImgFS.fill(0.0);
			psfImgFS.fill(0.0);
			dataImgRS.fill(0.0);
			dataImgDivRS.fill(0.0);
			
			for (int f = 0; f < fc; f++)
			{
				FourierBackprojection::backprojectSlice_forward_with_multiplicity(
					particleStack.getSliceRef(f),
					weightStack.getSliceRef(f),
					projPart[f] * relative_box_scale,
					dataImgFS,
					ctfImgFS,
					multiImageFS);
			}
			
			Centering::shiftInSitu(dataImgFS);
			FFT::inverseFourierTransform(dataImgFS, dataImgRS, FFT::Both);
			
			
			if (do_cone_weight)
			{ 
				FFT::FourierTransform(dataImgRS, dataImgFS);

				d3Matrix R = dataSet.getMatrix3x3(part_id);
				
				for (int z = 0; z < s3D;  z++)
				for (int y = 0; y < s3D;  y++)
				for (int x = 0; x < sh3D; x++)
				{
					const d3Vector p0(
						x,
						y < s3D/2? y : y - s3D,
						z < s3D/2? z : z - s3D);
					
					const d3Vector p = R * p0;
					
					const double rho = sqrt(p.x*p.x + p.y*p.y);
					const double t = rho / (std::abs(p.z) * cone_slope + cone_sig0);
					
					const double m = 1.0 - exp(-0.5*t*t);
					
					dataImgFS(x,y,z) *= m;
					psfImgFS(x,y,z) *= m;
					multiImageFS(x,y,z) *= m;
				}
				
				FFT::inverseFourierTransform(dataImgFS, dataImgRS);
			}

			std::string outData = outTag + "/" + dataSet.getName(part_id) + "_data.mrc";
			std::string outWeight = outTag + "/" + dataSet.getName(part_id) + "_weights.mrc";
			std::string outCTF = outTag + "/" + dataSet.getName(part_id) + "_CTF2.mrc";
			std::string outDiv = outTag + "/" + dataSet.getName(part_id) + "_div.mrc";
			std::string outMulti = outTag + "/" + dataSet.getName(part_id) + "_multi.mrc";
			std::string outNrm = outTag + "/" + dataSet.getName(part_id) + "_data_nrm.mrc";
			std::string outWeightNrm = outTag + "/" + dataSet.getName(part_id) + "_CTF2_nrm.mrc";

			// What if we didn't? The 2D image is already tapered.
			Reconstruction::taper(dataImgRS, taper, do_center, inner_thread_num);
			
			if (do_sum_all)
			{
				sum_data += dataImgRS;
				sum_weights += ctfImgFS;
			}
			
			
			dataImgRS.write(outData, binnedPixelSize);
			
			if (write_combined)
			{
				BufferedImage<float> ctfAndMultiplicity(sh3D,s3D,2*s3D);
				ctfAndMultiplicity.getSlabRef(0,s3D).copyFrom(ctfImgFS);
				ctfAndMultiplicity.getSlabRef(s3D,s3D).copyFrom(multiImageFS);
				
				ctfAndMultiplicity.write(outWeight, 1.0 / binnedPixelSize);
			}
			
			if (write_ctf)
			{
				Centering::fftwHalfToHumanFull(ctfImgFS).write(outCTF, 1.0 / binnedPixelSize);
			}
			
			if (write_multiplicity)
			{
				Centering::fftwHalfToHumanFull(multiImageFS).write(outMulti, 1.0 / binnedPixelSize);
			}
			
			if (write_normalised)
			{
				BufferedImage<float> ctfImgFSnrm = ctfImgFS;
				BufferedImage<fComplex> dataImgCorrFS;	
				
				FFT::FourierTransform(dataImgRS, dataImgCorrFS, FFT::Both);
				
				for (long int i = 0; i < ctfImgFSnrm.getSize(); i++)
				{
					const float n = multiImageFS[i];
					ctfImgFSnrm[i] = n > 0.f? ctfImgFS[i] / n : 0.f;
					dataImgCorrFS[i] = n > 0.f? dataImgCorrFS[i] / n : fComplex(0.f,0.f);
				}
				
				FFT::inverseFourierTransform(dataImgCorrFS, dataImgDivRS, FFT::Both);
				
				dataImgDivRS.write(outNrm, binnedPixelSize);
				Centering::fftwHalfToHumanFull(ctfImgFSnrm).write(outWeightNrm, 1.0 / binnedPixelSize);
			}
			
			if (write_divided)
			{
				if (SNR > 0.0)
				{
					Reconstruction::ctfCorrect3D_Wiener(
						dataImgRS, ctfImgFS, dataImgDivRS,
						1.0 / SNR, inner_thread_num);
				}
				else
				{
					Reconstruction::ctfCorrect3D_heuristic(
						dataImgRS, ctfImgFS, dataImgDivRS,
						0.001, inner_thread_num);
				}
				
				Reconstruction::taper(dataImgDivRS, taper, do_center, inner_thread_num);
				dataImgDivRS.write(outDiv, binnedPixelSize);
			}
		}
		
		Log::endProgress();
		
		Log::endSection(); // tomogram
	}

	if (do_sum_all)
	{
		sum_data.write(outTag + "_sum_data.mrc");
		sum_weights.write(outTag + "_sum_weight.mrc");

		BufferedImage<float> dataImgDivRS(s3D,s3D,s3D);
		dataImgDivRS.fill(0.0);

		if (SNR > 0.0)
		{
			Reconstruction::ctfCorrect3D_Wiener(
				sum_data, sum_weights, dataImgDivRS,
				1.0 / SNR, num_threads);
		}
		else
		{
			Reconstruction::ctfCorrect3D_heuristic(
				sum_data, sum_weights, dataImgDivRS,
				0.001, num_threads);
		}

		dataImgDivRS.write(outTag + "_sum_div.mrc");
	}
}

BufferedImage<float> SubtomoProgram::cropAndTaper(const BufferedImage<float>& imgFS, int boundary, int num_threads) const
{
	BufferedImage<fComplex> ctfImgFS_complex = FFT::toComplex(imgFS);
	
	BufferedImage<float> ctfImgRS;
	FFT::inverseFourierTransform(ctfImgFS_complex, ctfImgRS);
	
	ctfImgRS = Centering::fftwFullToHumanFull(ctfImgRS);
	ctfImgRS = Padding::unpadCenter3D_full(ctfImgRS, boundary);
	
	Reconstruction::GaussEnvelope(ctfImgRS, env_sigma, do_center, num_threads);				
	Reconstruction::taper(ctfImgRS, taper, do_center, num_threads);
	
	BufferedImage<float> ctfImgRS_cent = Centering::humanFullToFftwFull(ctfImgRS);
	
	FFT::FourierTransform(ctfImgRS_cent, ctfImgFS_complex);
	
	return FFT::toReal(ctfImgFS_complex);
}

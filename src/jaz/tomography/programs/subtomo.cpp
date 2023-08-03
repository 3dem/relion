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
#include <src/jaz/math/Euler_angles_relion.h>
#include <mpi.h>
#include <iostream>

using namespace gravis;


SubtomoProgram::SubtomoProgram()
: do_not_write_any(false),
  run_from_MPI(false)
{

}

void SubtomoProgram::readBasicParameters(IOParser& parser)
{
	optimisationSet.read(
		parser,
		true,           // optimisation set
		true,   true,   // particles
		true,   true,   // tomograms
		true,   false,  // trajectories
		false,  false,  // manifolds
		false,  false); // reference

	int gen_section = parser.addSection("Reconstruction options");

	boxSize = textToInteger(parser.getOption("--b", "Binned projection box size"));
	cropSize = textToInteger(parser.getOption("--crop", "Output box size", "-1"));
	binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));
	write_multiplicity = parser.checkOption("--multi", "Write out multiplicity volumes");
	SNR = textToDouble(parser.getOption("--SNR", "Assumed signal-to-noise ratio (negative means use a heuristic)", "-1"));

	do_cone_weight = parser.checkOption("--cone_weight", "Weight down a double cone along Z");
	const double alpha = 0.5 * textToDouble(parser.getOption("--cone_angle", "Opening angle of the cone in degrees", "10"));
	cone_slope = sin(DEG2RAD(alpha));
	cone_sig0 = textToDouble(parser.getOption("--cone_sig0", "Cone width at Z = 0", "2"));

	do_circle_crop = !parser.checkOption("--no_circle_crop", "Do not crop 2D images to a circle");
	do_circle_precrop = parser.checkOption("--circle_precrop", "Crop 2D images to the large circle (--b) prior to CTF modulation");
	do_narrow_circle_crop = true;
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

    apply_offsets = !parser.checkOption("--dont_apply_offsets", "By default, rlnOriginX/Y/ZAngst are combined with rlnCoordinateX/Y/Z to construct the particles in their refined translations. Use this argument to skip that.");
    apply_orientations = parser.checkOption("--apply_orientations", "rlnAngle<Rot/Tilt/Psi> are combined with rlnTomoSubtomogram<Rot/Tilt/Psi> to construct the particles in their refined orientations. This will also apply translations!");
    if (apply_orientations) apply_offsets = true;

	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Only process undone subtomograms");

	write_float16  = parser.checkOption("--float16", "Write in half-precision 16 bit floating point numbers (MRC mode 12), instead of 32 bit (MRC mode 0).");


	diag = parser.checkOption("--diag", "Write out diagnostic information");

	num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));

	freqCutoffFract = textToDouble(parser.getOption("--cutoff_fract", "Ignore shells for which the dose weight falls below this value", "0.01"));

	outDir = parser.getOption("--o", "Output filename pattern");

	run_from_GUI = is_under_pipeline_control();
}

void SubtomoProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;

	parser.setCommandLine(argc, argv);

	readBasicParameters(parser);

	do_sum_all = parser.checkOption("--sum", "Sum up all subtomograms (for debugging)");
	do_not_write_any = parser.checkOption("--no_writing", "Do not write out any files, only a sum");

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	if (do_gridding_precorrection)
	{
		do_narrow_circle_crop = true;
	}

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);
}

void SubtomoProgram::run()
{
	TomogramSet tomogramSet(optimisationSet.tomograms, true);

	ParticleSet particleSet(optimisationSet.particles, optimisationSet.trajectories, true);
	std::vector<std::vector<ParticleIndex> > particles = particleSet.splitByTomogram(tomogramSet, true);
	
	if (cropSize < 0) cropSize = boxSize;
	
	bool do_ctf = true;

	const long int s2D = boxSize;
	
	const long int s3D = cropSize;
	const long int sh3D = s3D / 2 + 1;
	
	const long int s02D = (int)(binning * s2D + 0.5);
	
	const double relative_box_scale = cropSize / (double) boxSize;
	const double binned_pixel_size = binning * particleSet.getOriginalPixelSize(0);

	initialise(particleSet, particles, tomogramSet);

	BufferedImage<float> sum_data, sum_weights;

	if (do_sum_all)
	{
		sum_data.resize(s3D,s3D,s3D);
		sum_data.fill(0.0);

		sum_weights.resize(sh3D,s3D,s3D);
		sum_weights.fill(0.0);
	}

	AberrationsCache aberrationsCache(particleSet.optTable, s2D, binned_pixel_size);


	std::vector<int> tomoIndices = ParticleSet::enumerate(particles);

	processTomograms(
		tomoIndices,
		tomogramSet,
		particleSet,
		particles,
		aberrationsCache,
		s02D,
		s2D,
		s3D,
		relative_box_scale,
		do_ctf,
		1,
		sum_data,
		sum_weights);


	if (do_sum_all)
	{
		const double pixel_size = binning * tomogramSet.getPixelSize(0);

		sum_data.write(outDir + "sum_data.mrc", pixel_size);
		Centering::fftwHalfToHumanFull(sum_weights).write(outDir + "sum_weight.mrc", pixel_size);

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

		dataImgDivRS.write(outDir + "sum_div.mrc", pixel_size);
	}
}

void SubtomoProgram::initialise(
		const ParticleSet& particleSet,
		const std::vector<std::vector<ParticleIndex>>& particles,
		const TomogramSet& tomogramSet)
{
	const int tc = tomogramSet.size();

	int firstGoodTomo = 0;

	for (int t = 0; t < tc; t++)
	{
		if (particles[t].size() > 0)
		{
			firstGoodTomo = t;
			break;
		}
	}

	const std::string firstName = particleSet.getName(particles[firstGoodTomo][0]);

	directoriesPerTomogram = firstName.find_first_of('/') == std::string::npos;

	if (directoriesPerTomogram)
	{
		Log::print("No slashes found in first particle name: creating subdirectories for each tomogram");
	}
	else
	{
		Log::print("Slash found in first particle name: not creating subdirectories for each tomogram");
	}

	for (int t = 0; t < tc; t++)
	{
		if (particles[t].size() > 0)
		{
			ZIO::ensureParentDir(getOutputFilename(
				particles[t][0], t, particleSet, tomogramSet));
		}
	}

	writeParticleSet(particleSet, particles, tomogramSet);
}

std::string SubtomoProgram::getOutputFilename(
		ParticleIndex p,
		int tomogramIndex,
		const ParticleSet& particleSet,
		const TomogramSet& tomogramSet)
{
	if (directoriesPerTomogram)
	{
		return outDir + "Subtomograms/"
				+ tomogramSet.getTomogramName(tomogramIndex) + "/"
				+ particleSet.getName(p);
	}
	else
	{
		return outDir + "Subtomograms/" + particleSet.getName(p);
	}
}

void SubtomoProgram::writeParticleSet(
		const ParticleSet& particleSet,
		const std::vector<std::vector<ParticleIndex>>& particles,
		const TomogramSet& tomogramSet)
{
	const int tc = particles.size();

	ParticleSet copy = particleSet;
	copy.clearParticles();

	int particles_removed = 0;

	for (int t = 0; t < tc; t++)
	{
		const int pc = particles[t].size();

		if (pc == 0) continue;

		for (int p = 0; p < pc; p++)
		{
			const ParticleIndex part_id = particles[t][p];

			Tomogram tomogram = tomogramSet.loadTomogram(t, false);

			const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
						part_id, tomogram.frameCount, tomogram.optics.pixelSize, !apply_offsets);

			if (tomogram.isVisibleAtAll(traj, boxSize / 2.0))
			{
				const ParticleIndex new_id = copy.addParticle(particleSet, part_id);

				const int opticsGroup = particleSet.getOpticsGroup(part_id);
				const double originalPixelSize = particleSet.getOriginalPixelSize(opticsGroup);

				const std::string filenameRoot = getOutputFilename(
					part_id, t, particleSet, tomogramSet);

                std::string outData = filenameRoot + "_data.mrc";
                std::string outWeight = filenameRoot + "_weights.mrc";

                copy.setImageFileNames(outData, outWeight, new_id);

                if (apply_offsets)
                {
                    const d3Matrix A_subtomogram = particleSet.getSubtomogramMatrix(part_id);
                    const d3Vector pos = particleSet.getParticleCoord(part_id) - (A_subtomogram * particleSet.getParticleOffset(part_id)) / originalPixelSize;
                    copy.setParticleOffset(new_id, d3Vector(0,0,0));
                    copy.setParticleCoord(new_id, pos);
                }

				if (apply_orientations)
				{
                    d3Matrix A = particleSet.getMatrix3x3(part_id);
					const gravis::d3Vector ang = Euler::matrixToAngles(A);

                    copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_ROT, RAD2DEG(ang[0]), new_id.value);
                    copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_TILT, RAD2DEG(ang[1]), new_id.value);
                    copy.partTable.setValue(EMDL_TOMO_SUBTOMOGRAM_PSI, RAD2DEG(ang[2]), new_id.value);

                    copy.partTable.setValue(EMDL_ORIENT_ROT, 0.0, new_id.value);
                    copy.partTable.setValue(EMDL_ORIENT_TILT, 0.0, new_id.value);
                    copy.partTable.setValue(EMDL_ORIENT_PSI, 0.0, new_id.value);
				}
			}
			else
			{
				particles_removed++;
			}
		}
	}

	if (particles_removed == 1)
	{
		Log::warn("One particle was removed because it was too close to the edge in all images.");
	}
	else if (particles_removed > 1)
	{
		Log::warn(ZIO::itoa(particles_removed)+" particles were removed because they were too close to the edge in all images.");
	}

	for (int og = 0; og < copy.numberOfOpticsGroups(); og++)
	{
		const double ps_img = copy.optTable.getDouble(EMDL_TOMO_TILT_SERIES_PIXEL_SIZE, og);
		const double ps_out = binning * ps_img;

		copy.optTable.setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, true, og);
		copy.optTable.setValue(EMDL_IMAGE_DIMENSIONALITY, 3, og);
		copy.optTable.setValue(EMDL_TOMO_SUBTOMOGRAM_BINNING, binning, og);
		copy.optTable.setValue(EMDL_IMAGE_PIXEL_SIZE, ps_out, og);
		copy.optTable.setValue(EMDL_IMAGE_SIZE, cropSize, og);
	}

	copy.write(outDir + "particles.star");

	if (copy.hasMotion && particles_removed > 0)
	{
		copy.writeTrajectories(outDir + "motion.star");
		optimisationSet.trajectories = outDir + "motion.star";
	}

	optimisationSet.particles = outDir + "particles.star";
	optimisationSet.write(outDir + "optimisation_set.star");
}

void SubtomoProgram::processTomograms(
		const std::vector<int>& tomoIndices,
		const TomogramSet& tomogramSet,
		const ParticleSet& particleSet,
		const std::vector<std::vector<ParticleIndex>>& particles,
		const AberrationsCache& aberrationsCache,
		long int s02D,
		long int s2D,
		long int s3D,
		double relative_box_scale,
		bool do_ctf,
		int verbosity,
		BufferedImage<float>& sum_data,
		BufferedImage<float>& sum_weights )
{
	const int tc = tomoIndices.size();
	const int sh2D = s2D / 2 + 1;
	const int sh3D = s3D / 2 + 1;

	for (int tt = 0; tt < tc; tt++)
	{
		const int t = tomoIndices[tt];

		const int pc = particles[t].size();
		if (pc == 0) continue;

		if (run_from_GUI && pipeline_control_check_abort_job())
		{
			if (run_from_MPI)
			{
				MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_ABORTED);
				exit(RELION_EXIT_ABORTED);
			}
			else
			{
				exit(RELION_EXIT_ABORTED);
			}
		}

		if (verbosity > 0)
		{
			Log::beginSection("Tomogram " + ZIO::itoa(tt+1) + " / " + ZIO::itoa(tc));
			Log::print("Loading");
		}

		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		tomogram.validateParticleOptics(particles[t], particleSet);

		const int fc = tomogram.frameCount;

		particleSet.checkTrajectoryLengths(particles[t], fc, "subtomo");

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s2D, binning);
		BufferedImage<float> noiseWeights;

		if (do_whiten)
		{
			noiseWeights = tomogram.computeNoiseWeight(s2D, binning);
		}

		BufferedImage<int> xRanges = tomogram.findDoseXRanges(doseWeights, freqCutoffFract);

		const int inner_thread_num = 1;
		const int outer_thread_num = num_threads / inner_thread_num;

		// @TODO: define input and output pixel sizes!

		const double binnedPixelSize = tomogram.optics.pixelSize * binning;

		if (verbosity > 0)
		{
			Log::beginProgress(
				"Backprojecting subtomograms",
				(int)ceil(pc/(double)outer_thread_num));
		}

		omp_lock_t writelock;
		if (do_sum_all) omp_init_lock(&writelock);

		#pragma omp parallel for num_threads(outer_thread_num)
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();

			if (verbosity > 0 && th == 0)
			{
				Log::updateProgress(p);
			}

			const ParticleIndex part_id = particles[t][p];

			const std::string filenameRoot = getOutputFilename(
				part_id, t, particleSet, tomogramSet);

			std::string outData = filenameRoot + "_data.mrc";
			std::string outWeight = filenameRoot + "_weights.mrc";
			std::string outCTF = filenameRoot + "_CTF2.mrc";
			std::string outDiv = filenameRoot + "_div.mrc";
			std::string outMulti = filenameRoot + "_multi.mrc";
			std::string outNrm = filenameRoot + "_data_nrm.mrc";
			std::string outWeightNrm = filenameRoot + "_CTF2_nrm.mrc";

			if (only_do_unfinished && ZIO::fileExists(outData))
			{
				continue;
			}
			
			const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
						part_id, fc, tomogram.optics.pixelSize, !apply_offsets);

			if (!tomogram.isVisibleAtAll(traj, s2D / 2.0))
			{
				continue;
			}

			const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s2D / 2.0);

			std::vector<d4Matrix> projCut(fc), projPart(fc);

			BufferedImage<fComplex> particleStack = BufferedImage<fComplex>(sh2D,s2D,fc);
			BufferedImage<float> weightStack(sh2D,s2D,fc);

			TomoExtraction::extractAt3D_Fourier(
					tomogram.stack, s02D, binning, tomogram, traj, isVisible,
					particleStack, projCut, inner_thread_num, do_circle_precrop);

			if (!do_ctf) weightStack.fill(1.f);


			const int og = particleSet.getOpticsGroup(part_id);

			const BufferedImage<double>* gammaOffset =
				aberrationsCache.hasSymmetrical? &aberrationsCache.symmetrical[og] : 0;

			for (int f = 0; f < fc; f++)
			{
				if (!isVisible[f]) continue;

				d3Matrix A;

				if (apply_orientations)
				{
					A = particleSet.getMatrix3x3(part_id);
				}
				else
				{
					A = particleSet.getSubtomogramMatrix(part_id);
				}

				projPart[f] = projCut[f] * d4Matrix(A);

				if (do_ctf)
				{
                    const d3Vector pos = (apply_offsets) ? particleSet.getPosition(part_id) : particleSet.getParticleCoord(part_id);

                    CTF ctf = tomogram.getCtf(f, pos);
					BufferedImage<float> ctfImg(sh2D, s2D);
					ctf.draw(s2D, s2D, binnedPixelSize, gammaOffset, &ctfImg(0,0,0));

					const float sign = flip_value? -1.f : 1.f;

					for (int y = 0; y < s2D;  y++)
					for (int x = 0; x < xRanges(y,f); x++)
					{
						const double c = ctfImg(x,y) * doseWeights(x,y,f);

						particleStack(x,y,f) *= sign * c;
						weightStack(x,y,f) = c * c;
					}
				}
			}

			aberrationsCache.correctObservations(particleStack, og);

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

			BufferedImage<float> ctfImgFS(sh3D,s3D,s3D),
					dataImgRS(s3D,s3D,s3D), dataImgDivRS(s3D,s3D,s3D),
					multiImageFS(sh3D,s3D,s3D);

			ctfImgFS.fill(0.0);
			dataImgRS.fill(0.0);
			dataImgDivRS.fill(0.0);

			for (int f = 0; f < fc; f++)
			{
				if (isVisible[f])
				{
					FourierBackprojection::backprojectSlice_forward_with_multiplicity(
						&xRanges(0,f),
						particleStack.getSliceRef(f),
						weightStack.getSliceRef(f),
						projPart[f] * relative_box_scale,
						dataImgFS,
						ctfImgFS,
						multiImageFS);
				}
			}

			Centering::shiftInSitu(dataImgFS);

			// correct FT scale after the implicit cropping:

			if (s3D != s2D)
			{
				dataImgFS *= (float) sqrt(s2D / (double) s3D);
			}

			FFT::inverseFourierTransform(dataImgFS, dataImgRS, FFT::Both);

			if (do_cone_weight)
			{
				FFT::FourierTransform(dataImgRS, dataImgFS);

				d3Matrix R = particleSet.getMatrix3x3(part_id);

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
					ctfImgFS(x,y,z) *= m;
					multiImageFS(x,y,z) *= m; // apply to both multiplicity and weight?
				}

				FFT::inverseFourierTransform(dataImgFS, dataImgRS);
			}

			// What if we didn't? The 2D image is already tapered.
			//Reconstruction::taper(dataImgRS, taper, do_center, inner_thread_num);

			if (do_sum_all)
			{
				omp_set_lock(&writelock);

				sum_data += dataImgRS;
				sum_weights += ctfImgFS;

				omp_unset_lock(&writelock);
			}

			if (do_not_write_any) continue;


			dataImgRS.write(outData, binnedPixelSize, write_float16);

			if (write_combined)
			{
				BufferedImage<float> ctfAndMultiplicity(sh3D,s3D,2*s3D);
				ctfAndMultiplicity.getSlabRef(0,s3D).copyFrom(ctfImgFS);
				ctfAndMultiplicity.getSlabRef(s3D,s3D).copyFrom(multiImageFS);

				ctfAndMultiplicity.write(outWeight, 1.0 / binnedPixelSize, write_float16);
			}

			if (write_ctf)
			{
				Centering::fftwHalfToHumanFull(ctfImgFS).write(outCTF, 1.0 / binnedPixelSize, write_float16);
			}

			if (write_multiplicity)
			{
				Centering::fftwHalfToHumanFull(multiImageFS).write(outMulti, 1.0 / binnedPixelSize, write_float16);
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

				dataImgDivRS.write(outNrm, binnedPixelSize, write_float16);
				Centering::fftwHalfToHumanFull(ctfImgFSnrm).write(outWeightNrm, 1.0 / binnedPixelSize, write_float16);
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
				dataImgDivRS.write(outDiv, binnedPixelSize, write_float16);
			}
		}

		if (verbosity > 0)
		{
			Log::endProgress();
			Log::endSection(); // tomogram
		}
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

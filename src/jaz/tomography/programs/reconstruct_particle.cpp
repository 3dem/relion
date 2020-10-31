#include "reconstruct_particle.h"
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/image/symmetry.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/time.h>
#include <iostream>


using namespace gravis;


void ReconstructParticleProgram::readParameters(int argc, char *argv[])
{
	readBasicParameters(argc, argv);

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);
}

void ReconstructParticleProgram::readBasicParameters(int argc, char *argv[])
{
	IOParser parser;
	
	try
	{
		parser.setCommandLine(argc, argv);

		optimisationSet.read(
			parser,
			true,           // optimisation set
			true,   true,   // particles
			true,   true,   // tomograms
			true,   false,  // trajectories
			false,  false,  // manifolds
			false,  false); // reference

		int gen_section = parser.addSection("Reconstruction options");

		boxSize = textToInteger(parser.getOption("--b", "Box size"));
		cropSize = textToInteger(parser.getOption("--crop", "Size of (additionally output) cropped image", "-1"));

		do_whiten = parser.checkOption("--whiten", "Whiten the noise by flattening the power spectrum");

		binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));
		taper = textToDouble(parser.getOption("--taper", "Taper against the sphere by this number of pixels (only if cropping)", "10"));
		SNR = textToDouble(parser.getOption("--SNR", "Assumed signal-to-noise ratio (negative means use a heuristic)", "-1"));
		symmName = parser.getOption("--sym", "Symmetry group", "C1");
				
		max_mem_GB = textToInteger(parser.getOption("--mem", "Max. amount of memory to use for accumulation (--j_out will be reduced)", "-1"));
				
		explicit_gridding = parser.checkOption("--xg", "Perform gridding correction using a measured spreading function");
		diag = parser.checkOption("--diag", "Write out diagnostic information");
		no_subpix_off = parser.checkOption("--nso", "No subpixel offset (debugging)");
		
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		inner_threads = textToInteger(parser.getOption("--j_in", "Number of inner threads (slower, needs less memory)", "3"));
		outer_threads = textToInteger(parser.getOption("--j_out", "Number of outer threads (faster, needs more memory)", "2"));
		
		no_reconstruction = parser.checkOption("--no_recon", "Do not reconstruct the volume, only backproject (for benchmarking purposes)");

		outDir = parser.getOption("--o", "Output directory");
		
		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
}

void ReconstructParticleProgram::run()
{
	Log::beginSection("Initialising");

	TomogramSet tomoSet(optimisationSet.tomograms);
	ParticleSet particleSet(optimisationSet.particles, optimisationSet.trajectories);

	std::vector<std::vector<ParticleIndex>> particles = particleSet.splitByTomogram(tomoSet);
	
	const int tc = particles.size();
	const int s = boxSize;
	const int sh = s/2 + 1;
	
	const int s02D = (int)(binning * s + 0.5);
	
	const bool flip_value = true;
	const bool do_ctf = true;

	Tomogram tomo0 = tomoSet.loadTomogram(0, false);
	const double binnedOutPixelSize = tomo0.optics.pixelSize * binning;

	
	const long int voxelNum = (long int) sh * (long int) s * (long int) s;
	
	if (max_mem_GB > 0)
	{
		const double GB_per_thread = 
		   2.0 * voxelNum * 3.0 * sizeof(double)   // two halves  *  box size  *  (data (x2) + ctf)
		    / (1024.0 * 1024.0 * 1024.0);          // in GB
		
		const double maxThreads = max_mem_GB / GB_per_thread;
		
		if (maxThreads < outer_threads)
		{
			int lastOuterThreads = outer_threads;
			outer_threads = (int) maxThreads;
		
			Log::print("Outer thread number reduced from " + ZIO::itoa(lastOuterThreads) + 
					  " to " + ZIO::itoa(outer_threads) + " due to memory constraints (--mem).");
		}
	}
	
	const int outCount = 2 * outer_threads;
		
	Log::print("Memory required for accumulation: " + ZIO::itoa(
			(3.0 * sizeof(double) * (long int) outCount * (double)voxelNum)
			  / (1024.0 * 1024.0 * 1024.0)
			) + " GB");
	
	std::vector<BufferedImage<double>> ctfImgFS(outCount), psfImgFS(outCount);
	std::vector<BufferedImage<dComplex>> dataImgFS(outCount);
	
	for (int i = 0; i < outCount; i++)
	{	
		dataImgFS[i] = BufferedImage<dComplex>(sh,s,s);
		ctfImgFS[i] = BufferedImage<double>(sh,s,s),
		psfImgFS[i] = BufferedImage<double>(sh,s,s);
		
		dataImgFS[i].fill(dComplex(0.0, 0.0));
		ctfImgFS[i].fill(0.0);
		psfImgFS[i].fill(0.0);
	}

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize);

	Log::endSection();
	


	processTomograms(
		0, tc-1, tomoSet, particleSet, particles, aberrationsCache,
		dataImgFS, ctfImgFS, psfImgFS, binnedOutPixelSize,
		s02D, do_ctf, flip_value, 1);



	if (outCount > 2)
	{		
		Log::print("Merging volumes");
	
		for (int i = 2; i < outCount; i++)
		{
			dataImgFS[i%2] += dataImgFS[i];
			ctfImgFS[i%2] += ctfImgFS[i];
			
			if (explicit_gridding)
			{
				psfImgFS[i%2] += psfImgFS[i];
			}
		}
	}

	
	if (no_reconstruction) return;
	

	if (symmName != "C1")
	{
		Log::print("Applying symmetry");
		
		for (int half = 0; half < 2; half++)
		{
			dataImgFS[half] = Symmetry::symmetrise_FS_complex(
						dataImgFS[half], symmName, num_threads);
			
			ctfImgFS[half] = Symmetry::symmetrise_FS_real(
						ctfImgFS[half], symmName, num_threads);
			
			if (explicit_gridding)
			{
				psfImgFS[half] = Symmetry::symmetrise_FS_real(
						psfImgFS[half], symmName, num_threads);
			}
		}
	}

	finalise(dataImgFS, ctfImgFS, psfImgFS, binnedOutPixelSize);
}

void ReconstructParticleProgram::processTomograms(
	int first_t,
	int last_t,
	const TomogramSet& tomoSet,
	const ParticleSet& particleSet,
	const std::vector<std::vector<ParticleIndex>>& particles,
	const AberrationsCache& aberrationsCache,
	std::vector<BufferedImage<dComplex>>& dataImgFS,
	std::vector<BufferedImage<double>>& ctfImgFS,
	std::vector<BufferedImage<double>>& psfImgFS,
	const double binnedOutPixelSize,
	int s02D,
	bool do_ctf,
	bool flip_value,
	int verbosity)
{
	const int s = dataImgFS[0].ydim;
	const int sh = s/2 + 1;
	const int tc = last_t - first_t + 1;

	for (int t = first_t; t <= last_t; t++)
	{
		const int pc = particles[t].size();

		if (pc == 0) continue;

		if (verbosity > 0)
		{
			Log::beginSection("Tomogram " + ZIO::itoa(t+1-first_t) + " / " + ZIO::itoa(tc));
			Log::print("Loading");
		}

		Tomogram tomogram = tomoSet.loadTomogram(t, true);

		const int fc = tomogram.frameCount;

		particleSet.checkTrajectoryLengths(particles[t][0], pc, fc, "backproject");

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, binning);
		BufferedImage<float> noiseWeights;

		if (do_whiten)
		{
			noiseWeights = tomogram.computeNoiseWeight(s, binning);
		}

		const double binnedPixelSize = tomogram.optics.pixelSize * binning;

		std::vector<BufferedImage<float>> weightStack(outer_threads, BufferedImage<float>(sh,s,fc));
		std::vector<BufferedImage<fComplex>> particleStack(outer_threads, BufferedImage<fComplex>(sh,s,fc));

		if (!do_ctf)
		{
			for (int i = 0; i < outer_threads; i++)
			{
				weightStack[i].fill(1.f);
			}
		}

		Log::beginProgress("Backprojecting", (int)ceil(pc/(double)outer_threads));

		#pragma omp parallel for num_threads(outer_threads)
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();

			if (th == 0)
			{
				Log::updateProgress(p);
			}

			const ParticleIndex part_id = particles[t][p];

			const d3Vector pos = particleSet.getPosition(part_id);
			const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
						part_id, fc, tomogram.optics.pixelSize);
			std::vector<d4Matrix> projCut(fc), projPart(fc);


			TomoExtraction::extractAt3D_Fourier(
					tomogram.stack, s02D, binning, tomogram.projectionMatrices, traj,
					particleStack[th], projCut, inner_threads, true);

			const d4Matrix particleToTomo = particleSet.getMatrix4x4(part_id, s,s,s);

			const int halfSet = particleSet.getHalfSet(part_id);

			const int og = particleSet.getOpticsGroup(part_id);

			const BufferedImage<double>* gammaOffset =
				aberrationsCache.hasSymmetrical? &aberrationsCache.symmetrical[og] : 0;

			for (int f = 0; f < fc; f++)
			{
				const double scaleRatio = binnedOutPixelSize / binnedPixelSize;
				projPart[f] = scaleRatio * projCut[f] * particleToTomo;

				if (do_ctf)
				{
					CTF ctf = tomogram.getCtf(f, pos);
					BufferedImage<float> ctfImg(sh,s);
					ctf.draw(s, s, binnedPixelSize, gammaOffset, &ctfImg(0,0,0));

					const float scale = flip_value? -1.f : 1.f;

					for (int y = 0; y < s;  y++)
					for (int x = 0; x < sh; x++)
					{
						const float c = scale * ctfImg(x,y) * doseWeights(x,y,f);

						particleStack[th](x,y,f) *= c;
						weightStack[th](x,y,f) = c * c;
					}
				}
			}

			if (aberrationsCache.hasAntisymmetrical)
			{
				aberrationsCache.correctObservations(particleStack[th], og);
			}

			if (do_whiten)
			{
				particleStack[th] *= noiseWeights;
				weightStack[th] *= noiseWeights;
			}

			if (explicit_gridding)
			{
				FourierBackprojection::backprojectStack_backward(
					particleStack[th], weightStack[th], projPart,
					dataImgFS[2*th + halfSet],
					psfImgFS[2*th + halfSet],
					ctfImgFS[2*th + halfSet],
					inner_threads);
			}
			else
			{
				for (int f = 0; f < fc; f++)
				{
					FourierBackprojection::backprojectSlice_backward(
						particleStack[th].getSliceRef(f),
						weightStack[th].getSliceRef(f),
						projPart[f],
						dataImgFS[2*th + halfSet],
						ctfImgFS[2*th + halfSet],
						inner_threads);
				}
			}
		}

		Log::endProgress();
		Log::endSection();
	}
}

void ReconstructParticleProgram::finalise(
	std::vector<BufferedImage<dComplex>>& dataImgFS,
	std::vector<BufferedImage<double>>& ctfImgFS,
	std::vector<BufferedImage<double>>& psfImgFS,
	const double binnedOutPixelSize)
{
	const int s = dataImgFS[0].ydim;

	std::vector<BufferedImage<double>> dataImgRS(2), dataImgDivRS(2);

	BufferedImage<dComplex> dataImgFS_both = dataImgFS[0] + dataImgFS[1];
	BufferedImage<double> ctfImgFS_both = ctfImgFS[0] + ctfImgFS[1];

	BufferedImage<double> psfImgFS_both;

	if (explicit_gridding)
	{
		psfImgFS_both = psfImgFS[0] + psfImgFS[1];
	}


	Log::beginSection("Reconstructing");

	for (int half = 0; half < 2; half++)
	{
		Log::print("Half " + ZIO::itoa(half));

		dataImgRS[half] = BufferedImage<double>(s,s,s);
		dataImgDivRS[half] = BufferedImage<double>(s,s,s);

		reconstruct(
			dataImgRS[half], dataImgDivRS[half], ctfImgFS[half],
			&psfImgFS[half], dataImgFS[half]);

		writeOutput(
			dataImgDivRS[half], dataImgRS[half], ctfImgFS[half],
			"half"+ZIO::itoa(half+1), binnedOutPixelSize);
	}

	reconstruct(
		dataImgRS[0], dataImgDivRS[0], ctfImgFS_both,
		&psfImgFS_both, dataImgFS_both);

	writeOutput(
		dataImgDivRS[0], dataImgRS[0], ctfImgFS[0],
			"merged", binnedOutPixelSize);

	Log::endSection();
}

void ReconstructParticleProgram::reconstruct(
		BufferedImage<double>& dataImgRS,
		BufferedImage<double>& dataImgDivRS,
		BufferedImage<double>& ctfImgFS,
		BufferedImage<double>* psfImgFS,
		BufferedImage<dComplex>& dataImgFS)
{
	if (explicit_gridding)
	{
		Reconstruction::griddingCorrect3D(
			dataImgFS, *psfImgFS, dataImgRS, true, num_threads);
	}
	else
	{
		Reconstruction::griddingCorrect3D_sinc2(
			dataImgFS, dataImgRS, true, num_threads);
	}

	if (SNR > 0.0)
	{
		Reconstruction::ctfCorrect3D_Wiener(
			dataImgRS, ctfImgFS, dataImgDivRS,
			1.0 / SNR, num_threads);
	}
	else
	{
		Reconstruction::ctfCorrect3D_heuristic(
			dataImgRS, ctfImgFS, dataImgDivRS,
			0.001, num_threads);
	}
}


void ReconstructParticleProgram::writeOutput(
		const BufferedImage<double>& corrected,
		const BufferedImage<double>& data,
		const BufferedImage<double>& weight,
		const std::string& tag,
		double pixelSize)
{
	data.write(outDir+"data_"+tag+".mrc", pixelSize);

	Centering::fftwHalfToHumanFull(weight).write(outDir+"weight_"+tag+".mrc", pixelSize);

	if (cropSize > 0 && cropSize < boxSize)
	{
		corrected.write(outDir+tag+"_full.mrc", pixelSize);

		BufferedImage<double> cropped = Padding::unpadCenter3D_full(
					corrected, (boxSize - cropSize)/2);

		Reconstruction::taper(cropped, taper, true, num_threads);

		cropped.write(outDir+tag+".mrc", pixelSize);
	}
	else
	{
		BufferedImage<double> tapered = corrected;
		Reconstruction::taper(tapered, taper, true, num_threads);

		tapered.write(outDir+tag+".mrc", pixelSize);
	}
}

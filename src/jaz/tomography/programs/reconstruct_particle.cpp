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
#include <mpi.h>
#include <iostream>


using namespace gravis;


ReconstructParticleProgram::ReconstructParticleProgram()
: run_from_MPI(false)
{

}

void ReconstructParticleProgram::readParameters(int argc, char *argv[])
{
	readBasicParameters(argc, argv);

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);

	ZIO::makeDir(outDir + "temp");
	tmpOutRoot = outDir + "temp/sum_";

	run_from_GUI = is_under_pipeline_control();
}

void ReconstructParticleProgram::readBasicParameters(int argc, char *argv[])
{
	IOParser parser;

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

	nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
	helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
	helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

	max_mem_GB = textToInteger(parser.getOption("--mem", "Max. amount of memory (in GB) to use for accumulation (--j_out will be reduced)", "-1"));

	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Only process undone subtomograms");
	no_backup = parser.checkOption("--no_backup", "Do not make backups (makes it impossible to use --only_do_unfinished)");

	do_circle_crop = !parser.checkOption("--no_circle_crop", "Do not crop 2D images to a circle prior to insertion");

	num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
	inner_threads = textToInteger(parser.getOption("--j_in", "Number of inner threads (slower, needs less memory)", "3"));
	outer_threads = textToInteger(parser.getOption("--j_out", "Number of outer threads (faster, needs more memory)", "2"));

	no_reconstruction = parser.checkOption("--no_recon", "Do not reconstruct the volume, only backproject (for benchmarking purposes)");
	freqCutoffFract = textToDouble(parser.getOption("--cutoff_fract", "Ignore shells for which the dose weight falls below this value", "0.01"));

	outDir = parser.getOption("--o", "Output directory");

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
}

void ReconstructParticleProgram::run()
{
	Log::beginSection("Initialising");

	TomogramSet tomoSet(optimisationSet.tomograms, true);
	ParticleSet particleSet(optimisationSet.particles, optimisationSet.trajectories, true);

	std::vector<std::vector<ParticleIndex>> particles = particleSet.splitByTomogram(tomoSet, true);
	
	const int tc = particles.size();
	const int s = boxSize;
	const int sh = s/2 + 1;
	
	const int s02D = (int)(binning * s + 0.5);
	
	const bool flip_value = true;
	const bool do_ctf = true;

	Tomogram tomo0 = tomoSet.loadTomogram(0, false);
	const double binnedOutPixelSize = tomo0.optics.pixelSize * binning;

	
	const long int voxelNum = (long int) sh * (long int) s * (long int) s;

	const double GB_per_thread =
			2.0 * voxelNum * 3.0 * sizeof(double)   // two halves  *  box size  *  (data (x2) + ctf)
			/ (1024.0 * 1024.0 * 1024.0);           // in GB

	if (max_mem_GB > 0)
	{
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
		
	Log::print("Memory required for accumulation: " + ZIO::itoa(GB_per_thread  * outer_threads) + " GB");
	
	std::vector<BufferedImage<double>> ctfImgFS(outCount);
	std::vector<BufferedImage<dComplex>> dataImgFS(outCount);
	
	for (int i = 0; i < outCount; i++)
	{	
		dataImgFS[i] = BufferedImage<dComplex>(sh,s,s);
		ctfImgFS[i] = BufferedImage<double>(sh,s,s);
		
		dataImgFS[i].fill(dComplex(0.0, 0.0));
		ctfImgFS[i].fill(0.0);
	}

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize, binnedOutPixelSize);

	Log::endSection();
	

	std::vector<int> tomoIndices = ParticleSet::enumerate(particles);

	processTomograms(
		tomoIndices, tomoSet, particleSet, particles, aberrationsCache,
		dataImgFS, ctfImgFS, binnedOutPixelSize,
		s02D, do_ctf, flip_value, 1, true);

	if (no_reconstruction) return;

	finalise(dataImgFS, ctfImgFS, binnedOutPixelSize);

	// Delete temporary files
	// No error checking - do not bother the user if it fails
	if (system(("rm -rf "+ tmpOutRoot + "*.mrc").c_str()))
	{
		Log::warn("Deleting temporary files in folder "+tmpOutRoot+" failed.");
	}
}

void ReconstructParticleProgram::processTomograms(
	const std::vector<int>& tomoIndices,
	const TomogramSet& tomoSet,
	const ParticleSet& particleSet,
	const std::vector<std::vector<ParticleIndex>>& particles,
	const AberrationsCache& aberrationsCache,
	std::vector<BufferedImage<dComplex>>& dataImgFS,
	std::vector<BufferedImage<double>>& ctfImgFS,
	const double binnedOutPixelSize,
	int s02D,
	bool do_ctf,
	bool flip_value,
	int verbosity,
	bool per_tomogram_progress)
{
	const int s = dataImgFS[0].ydim;
	const int sh = s/2 + 1;
	const int tc = tomoIndices.size();

	if (verbosity > 0 && !per_tomogram_progress)
	{
		int total_particles_on_first_thread = 0;

		for (int tt = 0; tt < tc; tt++)
		{
			const int t = tomoIndices[tt];
			const int pc_all = particles[t].size();
			const int pc_th0 = (int)ceil(pc_all/(double)outer_threads);

			total_particles_on_first_thread += pc_th0;
		}

		Log::beginProgress("Backprojecting", total_particles_on_first_thread);
	}

	int particles_in_previous_tomograms = 0;


	int ttIni = 0, ttPrevious = -1;

	if (only_do_unfinished)
	{
		for (int tt = tc-1; tt > -1; tt--)
		{
			if (ZIO::fileExists(tmpOutRoot + ZIO::itoa(tt) + "_data_half1.mrc"))
			{
				ttPrevious = tt;
				ttIni = tt + 1;
				break;
			}
		}

		if (ttIni > 0)
		{
			//Read temporary files
			std::vector<BufferedImage<double>> tmpDataImg(2);

			std::string tmpOutRootTT = tmpOutRoot + ZIO::itoa(ttIni-1);

			for (int half = 0; half < 2; half++)
			{
				tmpDataImg[half].read(tmpOutRootTT + "_data_half" + ZIO::itoa(half) + ".mrc");
				ctfImgFS[half].read(tmpOutRootTT + "_ctf_half" + ZIO::itoa(half) + ".mrc");

				for (int z = 0; z < s;  z++)
				for (int y = 0; y < s;  y++)
				for (int x = 0; x < sh; x++)
				{
					dataImgFS[half](x,y,z) = fComplex(tmpDataImg[half](x,y,z), tmpDataImg[half](x,y,z+s));
				}
			}
		}
	}

	for (int tt = ttIni; tt < tc; tt++)
	{
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

		const int t = tomoIndices[tt];
		const int pc = particles[t].size();

		if (pc == 0) continue;

		if (verbosity > 0)
		{
			if (per_tomogram_progress)
			{
				Log::beginSection("Tomogram " + ZIO::itoa(tt+1) + " / " + ZIO::itoa(tc));
				Log::print("Loading");
			}
		}

		Tomogram tomogram = tomoSet.loadTomogram(t, true);
		tomogram.validateParticleOptics(particles[t], particleSet);

		const int fc = tomogram.frameCount;

		particleSet.checkTrajectoryLengths(particles[t], fc, "reconstruct_particle");

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, binning);

		BufferedImage<int> xRanges = tomogram.findDoseXRanges(doseWeights, freqCutoffFract);

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

		if (verbosity > 0 && per_tomogram_progress)
		{
			Log::beginProgress("Backprojecting", (int)ceil(pc/(double)outer_threads));
		}

		#pragma omp parallel for num_threads(outer_threads)
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();

			if (th == 0 && verbosity > 0)
			{
				if (per_tomogram_progress)
				{
					Log::updateProgress(p);
				}
				else
				{
					Log::updateProgress(particles_in_previous_tomograms + p);
				}
			}

			const ParticleIndex part_id = particles[t][p];

			const d3Vector pos = particleSet.getPosition(part_id);
			const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
						part_id, fc, tomogram.optics.pixelSize);
			std::vector<d4Matrix> projCut(fc), projPart(fc);

			const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s/2.0);

			const bool circle_crop = do_circle_crop;

			TomoExtraction::extractAt3D_Fourier(
					tomogram.stack, s02D, binning, tomogram, traj, isVisible,
					particleStack[th], projCut, inner_threads, circle_crop);


			const d4Matrix particleToTomo = particleSet.getMatrix4x4(part_id, s,s,s);

			const int halfSet = particleSet.getHalfSet(part_id);

			const int og = particleSet.getOpticsGroup(part_id);

			const BufferedImage<double>* gammaOffset =
				aberrationsCache.hasSymmetrical? &aberrationsCache.symmetrical[og] : 0;

			for (int f = 0; f < fc; f++)
			{
				if (!isVisible[f]) continue;
				
				const double scaleRatio = binnedOutPixelSize / binnedPixelSize;
				projPart[f] = scaleRatio * projCut[f] * particleToTomo;

				if (do_ctf)
				{
					CTF ctf = tomogram.getCtf(f, pos);
					BufferedImage<float> ctfImg(sh,s);
					ctf.draw(s, s, binnedPixelSize, gammaOffset, &ctfImg(0,0,0));

					const float scale = flip_value? -1.f : 1.f;

					for (int y = 0; y < s;  y++)
					{
						for (int x = 0; x < xRanges(y,f); x++)
						{
							const float c = scale * ctfImg(x,y) * doseWeights(x,y,f);

							particleStack[th](x,y,f) *= c;
							weightStack[th](x,y,f) = c * c;
						}
						for (int x = xRanges(y,f); x < sh; x++)
						{

							particleStack[th](x,y,f) = fComplex(0.f, 0.f);
							weightStack[th](x,y,f) = 0.f;
						}
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

			for (int f = 0; f < fc; f++)
			{
				if (isVisible[f])
				{
					FourierBackprojection::backprojectSlice_backward(
						xRanges(0,f),
						particleStack[th].getSliceRef(f),
						weightStack[th].getSliceRef(f),
						projPart[f],
						dataImgFS[2*th + halfSet],
						ctfImgFS[2*th + halfSet],
						inner_threads);
				}
			}

		} // particles

		if (!no_backup)
		{
			for (int i = 2; i < dataImgFS.size(); i++)
			{
				dataImgFS[i%2] += dataImgFS[i];
				ctfImgFS[i%2]  += ctfImgFS[i];
			}

			for (int i = 2; i < dataImgFS.size(); i++)
			{
				dataImgFS[i].fill(dComplex(0.0, 0.0));
				ctfImgFS[i].fill(0.0);
			}

			//Save temporary files

			for (int half = 0; half < 2; half++)
			{
				BufferedImage<double> tmpDataImg(sh, s, s*2);

				for (int z = 0; z < s;  z++)
				for (int y = 0; y < s;  y++)
				for (int x = 0; x < sh; x++)
				{
					fComplex pv = dataImgFS[half](x,y,z);
					tmpDataImg(x,y,z) = pv.real;
					tmpDataImg(x,y,z+s) = pv.imag;
				}

				std::string tmpOutRootTT = tmpOutRoot + ZIO::itoa(tt);
				tmpDataImg.write(tmpOutRootTT + "_data_half" + ZIO::itoa(half) + ".mrc");
				ctfImgFS[half].write(tmpOutRootTT + "_ctf_half" + ZIO::itoa(half) + ".mrc");
			}

			// Delete temporary files from previous tomogram
			// Intentionally no error checking
			if (ttPrevious > -1)
			{
				if (system(("rm -rf "+ tmpOutRoot  + ZIO::itoa(ttPrevious) + "*.mrc").c_str()))
					std::cerr << "WARNING: deleting temporary files " <<
					tmpOutRoot  + ZIO::itoa(ttPrevious) + "*.mrc failed." << std::endl;
			}

			ttPrevious = tt;
		}

		if (verbosity > 0 && per_tomogram_progress)
		{
			Log::endProgress();
			Log::endSection();
		}


		particles_in_previous_tomograms += (int)ceil(pc/(double)outer_threads);

	} // tomograms

	if (no_backup)
	{
		for (int i = 2; i < dataImgFS.size(); i++)
		{
			dataImgFS[i%2] += dataImgFS[i];
			ctfImgFS[i%2]  += ctfImgFS[i];
		}
	}

	if (verbosity > 0 && !per_tomogram_progress)
	{
		Log::endProgress();
	}
}

void ReconstructParticleProgram::finalise(
	std::vector<BufferedImage<dComplex>>& dataImgFS,
	std::vector<BufferedImage<double>>& ctfImgFS,
	const double binnedOutPixelSize)
{
	const int s = dataImgFS[0].ydim;

	symmetrise(dataImgFS, ctfImgFS, binnedOutPixelSize);

	std::vector<BufferedImage<double>> dataImgRS(2), dataImgDivRS(2);

	BufferedImage<dComplex> dataImgFS_both = dataImgFS[0] + dataImgFS[1];
	BufferedImage<double> ctfImgFS_both = ctfImgFS[0] + ctfImgFS[1];


	Log::beginSection("Reconstructing");

	for (int half = 0; half < 2; half++)
	{
		Log::print("Half " + ZIO::itoa(half));

		dataImgRS[half] = BufferedImage<double>(s,s,s);
		dataImgDivRS[half] = BufferedImage<double>(s,s,s);

		reconstruct(
			dataImgRS[half], dataImgDivRS[half], ctfImgFS[half],
			dataImgFS[half]);

		writeOutput(
			dataImgDivRS[half], dataImgRS[half], ctfImgFS[half],
			"half"+ZIO::itoa(half+1), binnedOutPixelSize);
	}

	reconstruct(
		dataImgRS[0], dataImgDivRS[0], ctfImgFS_both,
		dataImgFS_both);

	writeOutput(
		dataImgDivRS[0], dataImgRS[0], ctfImgFS[0],
			"merged", binnedOutPixelSize);

	optimisationSet.refMap1 = outDir + "half1.mrc";
	optimisationSet.refMap2 = outDir + "half2.mrc";
	optimisationSet.refFSC = "";
	optimisationSet.write(outDir + "optimisation_set.star");

	Log::endSection();
}

void ReconstructParticleProgram::symmetrise(
		std::vector<BufferedImage<dComplex> >& dataImgFS,
		std::vector<BufferedImage<double> >& ctfImgFS,
		double binnedOutPixelSize)
{
	if (symmName != "C1" && symmName != "c1")
	{
		std::vector<gravis::d4Matrix> symmetryMatrices;

		if (nr_helical_asu == 1)
		{
			Log::print("Applying point-group symmetries");

			symmetryMatrices = Symmetry::getPointGroupMatrices(symmName);
		}
		else
		{
			Log::print("Applying helical symmetries");

			symmetryMatrices = Symmetry::getHelicalSymmetryMatrices(
						nr_helical_asu, helical_twist, helical_rise/binnedOutPixelSize);
		}

		for (int half = 0; half < 2; half++)
		{
			dataImgFS[half] = Symmetry::symmetrise_FS_complex(
						dataImgFS[half], symmetryMatrices, num_threads);

			ctfImgFS[half] = Symmetry::symmetrise_FS_real(
						ctfImgFS[half], symmetryMatrices, num_threads);
		}
	}
}

void ReconstructParticleProgram::reconstruct(
		BufferedImage<double>& dataImgRS,
		BufferedImage<double>& dataImgDivRS,
		BufferedImage<double>& ctfImgFS,
		BufferedImage<dComplex>& dataImgFS)
{
	Reconstruction::griddingCorrect3D_sinc2(
			dataImgFS, dataImgRS, true, num_threads);

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

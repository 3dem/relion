#include "frame_alignment_program.h"
#include <src/ctf.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/motion/proto_alignment.h>
#include <src/jaz/tomography/projection_IO.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/math/fcc.h>
#include <src/time.h>
#include <omp.h>


using namespace gravis;


FrameAlignmentProgram::FrameAlignmentProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{	
	IOParser parser;
	parser.setCommandLine(argc, argv);
	readParams(parser);
}

void FrameAlignmentProgram::readParams(IOParser &parser)
{
	try
	{
		_readParams(parser);
				
		int defocus_section = parser.addSection("Alignment options");
		
		shiftOnly = parser.checkOption("--shift_only", "Only apply an optimal rigid shift to each frame (no iterative optimisation)");
		range = textToInteger(parser.getOption("--r", "Max. shift allowed [Pixels]", "20"));
		padding = textToDouble(parser.getOption("--pad", "Apply Fourier padding to the cross-correlation images", "1"));
		whiten = !parser.checkOption("--no_whiten", "Do not whiten the noise spectra");
		whiten_abs = parser.checkOption("--whiten_abs", "Divide by the square root of the power spectrum");
		hiPass_px = textToDouble(parser.getOption("--hp", "High-pass filter the cross-correlation images by this sigma", "-1"));
		sig2RampPower = textToDouble(parser.getOption("--srp", "Noise variance is divided by k^this during whitening", "0"));
		const_particles = parser.checkOption("--const_p", "Keep the particle positions constant");
		const_angles = parser.checkOption("--const_a", "Keep the frame angles constant");
		const_shifts = parser.checkOption("--const_s", "Keep the frame shifts constant");
		num_iters = textToInteger(parser.getOption("--it", "Max. number of iterations", "10000"));
		
		Log::readParams(parser);
		
		if (parser.checkForErrors()) std::exit(-1);
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
}

void FrameAlignmentProgram::run()
{
	Log::beginSection("Initialising");

		initialise();

		AberrationsCache aberrationsCache(particleSet.optTable, boxSize);

	Log::endSection();


	std::vector<int> tomoIndices = ParticleSet::enumerate(particles);

	processTomograms(tomoIndices, aberrationsCache, 1, true);
	

	finalise();
}

void FrameAlignmentProgram::initialise()
{
	RefinementProgram::init();

	ZIO::makeDir(outDir + "/Trajectories");

	const int tc = particles.size();

	Log::beginSection("Configuration");
	Log::printBinaryChoice("Frame angles: ", const_angles, "static", "variable");
	Log::printBinaryChoice("Frame shifts: ", const_shifts, "static", "variable");
	Log::printBinaryChoice("Particle positions: ", const_particles, "static", "variable");
	Log::endSection();

	ZIO::ensureParentDir(getTempFilenameRoot(""));
}

void FrameAlignmentProgram::finalise()
{
	const int tc = particles.size();

	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;

		readTempData(t);
	}

	tomogramSet.write(outDir + "tomograms.star");
	optimisationSet.tomograms = outDir+"tomograms.star";

	if (!shiftOnly)
	{
		particleSet.write(outDir + "particles.star");
		optimisationSet.particles = outDir+"particles.star";
	}

	optimisationSet.write(outDir+"optimisation_set.star");
}

void FrameAlignmentProgram::processTomograms(
		const std::vector<int> &tomoIndices,
		const AberrationsCache &aberrationsCache,
		int verbosity,
		bool per_tomogram_progress)
{
	const int ttc = tomoIndices.size();

	if (verbosity > 0 && !per_tomogram_progress)
	{
		Log::beginProgress("Processing tomograms", ttc);
	}

	for (int tt = 0; tt < ttc; tt++)
	{
		const int t = tomoIndices[tt];

		if (verbosity > 0 && !per_tomogram_progress)
		{
			Log::updateProgress(tt);
		}

		const std::string temp_filename_root = getTempFilenameRoot(
					tomogramSet.getTomogramName(t));

		if (only_do_unfinished &&
				ZIO::fileExists(temp_filename_root + "_projections.star"))
		{
			continue;
		}

		int pc = particles[t].size();
		if (pc == 0) continue;

		if (verbosity > 0 && per_tomogram_progress)
		{
			Log::beginSection("Tomogram " + ZIO::itoa(tt+1) + " / " + ZIO::itoa(ttc));
			Log::print("Loading");
		}

		Tomogram tomogram = tomogramSet.loadTomogram(t, true);

		const int fc = tomogram.frameCount;
		std::vector<d4Matrix> newProjections = tomogram.projectionMatrices;

		std::string tag = ZIO::itoa(t);
		std::string diagPrefix = outDir + "diag_" + tag;


		BufferedImage<float> frqWeight = computeFrequencyWeights(
			tomogram, whiten, sig2RampPower, hiPass_px, true, num_threads);

		std::vector<int> dummySeq(fc);

		for (int f = 0; f < fc; f++)
		{
			dummySeq[f] = f;
		}


		std::vector<BufferedImage<double>> CCs = Prediction::computeCroppedCCs(
				particleSet, particles[t], tomogram, aberrationsCache,
				referenceMap, frqWeight, dummySeq,
				range, true, num_threads, padding);


		ProtoAlignment protoAlignment(
				CCs, tomogram.projectionMatrices, particleSet, particles[t], referenceMap.image_FS,
				const_particles, const_angles, const_shifts, range,
				tomogram.centre, 0, num_threads, padding);

		BufferedImage<double> FCC3, FCC1;

		if (diag)
		{
			FCC3 = FCC::compute3(
				particleSet, particles[t], tomogram, referenceMap.image_FS,
				true, num_threads);

			FCC3.write(diagPrefix + "_FCC3_initial.mrc");
			FCC1 = FCC::divide(FCC3);
			FCC1.write(diagPrefix + "_FCC_initial.mrc");
		}

		std::vector<d3Vector> newPositions;

		if (shiftOnly)
		{
			const int diam = (int)(2*range*padding);

			BufferedImage<float> CCsum(diam, diam, fc);

			CCsum.fill(0.f);

			if (verbosity > 0 && per_tomogram_progress)
			{
				Log::beginProgress("Adding up cross-correlations", pc);
			}

			for (int p = 0; p < pc; p++)
			{
				if (verbosity > 0 && per_tomogram_progress)
				{
					Log::updateProgress(p);
				}

				CCsum += protoAlignment.CCs[p];
			}

			if (verbosity > 0 && per_tomogram_progress)
			{
				Log::endProgress();
			}

			CCsum.write(outDir + "CCsum_" + tag + ".mrc");

			d2Vector origin(padding*range, padding*range);

			std::ofstream frameShifts(outDir + "frame_shifts_" + tag + ".txt");

			for (int f = 0; f < fc; f++)
			{
				d2Vector opt = (Interpolation::quadraticMaxXY(CCsum.getSliceRef(f)) - origin)/padding;

				frameShifts << f << " " << opt.x << " " << opt.y << std::endl;

				newProjections[f](0,3) += opt.x;
				newProjections[f](1,3) += opt.y;
			}

			newPositions.resize(pc);

			for (int p = 0; p < pc; p++)
			{
				newPositions[p] = particleSet.getPosition(particles[t][p]);
			}
		}
		else
		{
			std::vector<double> initial(protoAlignment.getParamCount(), 0.0);

			if (verbosity > 0 && per_tomogram_progress)
			{
				Log::beginProgress("Performing optimisation", num_iters);
			}

			std::vector<double> opt = LBFGS::optimize(
				initial, protoAlignment, 1, num_iters, 1e-6, 1e-4);

			if (verbosity > 0 && per_tomogram_progress)
			{
				Log::endProgress();
			}

			newProjections = protoAlignment.getProjections(opt);
			newPositions = protoAlignment.getParticlePositions(opt);
		}

		writeTempData(newProjections, newPositions, t);

		if (diag)
		{
			tomogram.projectionMatrices = newProjections;

			FCC3 = FCC::compute3(
				particleSet, particles[t], tomogram, referenceMap.image_FS,
				true, num_threads);

			FCC3.write(diagPrefix + "_FCC3_final.mrc");
			BufferedImage<double> FCC1_new = FCC::divide(FCC3);
			FCC1_new.write(diagPrefix + "_FCC_final.mrc");

			(FCC1_new - FCC1).write(diagPrefix + "_FCC_delta.mrc");
		}

		Log::endSection();
	}

	if (verbosity > 0 && !per_tomogram_progress)
	{
		Log::endProgress();
	}
}

std::string FrameAlignmentProgram::getTempFilenameRoot(
		const std::string &tomogram_name)
{
	return outDir + "temp/" + tomogram_name;
}

void FrameAlignmentProgram::writeTempData(
		const std::vector<d4Matrix> &proj,
		const std::vector<d3Vector> &pos,
		int t)
{
	const int pc = particles[t].size();
	const int fc = tomogramSet.getFrameCount(t);

	const std::string tomoName = tomogramSet.getTomogramName(t);
	const std::string temp_filename_root = getTempFilenameRoot(tomoName);


	MetaDataTable temp_positions;

	for (int p = 0; p < pc; p++)
	{
		temp_positions.addObject();

		temp_positions.setValue(EMDL_IMAGE_COORD_X, pos[p].x - 1, p);
		temp_positions.setValue(EMDL_IMAGE_COORD_Y, pos[p].y - 1, p);
		temp_positions.setValue(EMDL_IMAGE_COORD_Z, pos[p].z - 1, p);
	}

	temp_positions.write(temp_filename_root + "_positions.star");


	MetaDataTable temp_projections;

	for (int f = 0; f < fc; f++)
	{
		const d4Matrix& P = proj[f];

		temp_projections.addObject();

		temp_projections.setValue(EMDL_TOMO_PROJECTION_X, std::vector<double>{P(0,0), P(0,1), P(0,2), P(0,3)}, f);
		temp_projections.setValue(EMDL_TOMO_PROJECTION_Y, std::vector<double>{P(1,0), P(1,1), P(1,2), P(1,3)}, f);
		temp_projections.setValue(EMDL_TOMO_PROJECTION_Z, std::vector<double>{P(2,0), P(2,1), P(2,2), P(2,3)}, f);
		temp_projections.setValue(EMDL_TOMO_PROJECTION_W, std::vector<double>{P(3,0), P(3,1), P(3,2), P(3,3)}, f);
	}

	temp_projections.write(temp_filename_root + "_projections.star");
}

void FrameAlignmentProgram::readTempData(int t)
{
	const int pc = particles[t].size();

	if (pc == 0) return;

	const int fc = tomogramSet.getFrameCount(t);

	const std::string tomoName = tomogramSet.getTomogramName(t);
	const std::string temp_filename_root = getTempFilenameRoot(tomoName);


	MetaDataTable temp_positions;
	temp_positions.read(temp_filename_root + "_positions.star");

	for (int p = 0; p < pc; p++)
	{
		d3Vector imgCoord;

		imgCoord.x = temp_positions.getDouble(EMDL_IMAGE_COORD_X, p);
		imgCoord.y = temp_positions.getDouble(EMDL_IMAGE_COORD_Y, p);
		imgCoord.z = temp_positions.getDouble(EMDL_IMAGE_COORD_Z, p);

		particleSet.partTable.setValue(EMDL_IMAGE_COORD_X, imgCoord.x, particles[t][p].value);
		particleSet.partTable.setValue(EMDL_IMAGE_COORD_Y, imgCoord.y, particles[t][p].value);
		particleSet.partTable.setValue(EMDL_IMAGE_COORD_Z, imgCoord.z, particles[t][p].value);

		particleSet.partTable.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, 0.0, particles[t][p].value);
		particleSet.partTable.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, 0.0, particles[t][p].value);
		particleSet.partTable.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, 0.0, particles[t][p].value);
	}


	MetaDataTable temp_projections;
	temp_projections.read(temp_filename_root + "_projections.star");

	for (int f = 0; f < fc; f++)
	{
		const std::vector<double> X = temp_projections.getDoubleVector(EMDL_TOMO_PROJECTION_X, f);
		const std::vector<double> Y = temp_projections.getDoubleVector(EMDL_TOMO_PROJECTION_Y, f);
		const std::vector<double> Z = temp_projections.getDoubleVector(EMDL_TOMO_PROJECTION_Z, f);
		const std::vector<double> W = temp_projections.getDoubleVector(EMDL_TOMO_PROJECTION_W, f);

		const d4Matrix P(
				X[0], X[1], X[2], X[3],
				Y[0], Y[1], Y[2], Y[3],
				Z[0], Z[1], Z[2], Z[3],
				W[0], W[1], W[2], W[3] );

		tomogramSet.setProjection(t, f, P);
	}
}

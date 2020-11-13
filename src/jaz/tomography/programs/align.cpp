#include "align.h"
#include <src/ctf.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/motion/motion_fit.h>
#include <src/jaz/tomography/motion/proto_alignment.h>
#include <src/jaz/tomography/motion/trajectory_set.h>
#include <src/jaz/tomography/projection_IO.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/index_sort.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/math/fcc.h>
#include <src/jaz/util/log.h>
#include <omp.h>

using namespace gravis;


AlignProgram::AlignProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	readParams(parser);
}

void AlignProgram::readParams(IOParser &parser)
{
	try
	{
		_readParams(parser);


		int alignment_section = parser.addSection("General alignment options");

		shiftOnly = parser.checkOption("--shift_only", "Only apply an optimal rigid shift to each frame (no iterative optimisation)");
		range = textToInteger(parser.getOption("--r", "Max. particle shift allowed [Pixels]", "20"));
		mfSettings.constParticles = parser.checkOption("--const_p", "Keep the particle positions constant");
		mfSettings.constAngles = parser.checkOption("--const_a", "Keep the frame angles constant");
		mfSettings.constShifts = parser.checkOption("--const_s", "Keep the frame shifts constant");
		num_iters = textToInteger(parser.getOption("--it", "Max. number of iterations", "1000"));


		int motion_section = parser.addSection("Motion estimation options");

		do_motion = parser.checkOption("--motion", "Estimate particle motion (expensive)");

		motParams.sig_vel = textToDouble(parser.getOption("--s_vel", "Velocity sigma [Å/dose]", "0.5"));
		motParams.sig_div = textToDouble(parser.getOption("--s_div", "Divergence sigma [Å]", "5000.0"));

		mfSettings.params_scaled_by_dose = !parser.checkOption("--abs_params", "Do not scale the sigmas by the dose");

		mfSettings.sqExpKernel = parser.checkOption("--sq_exp_ker", "Use a square-exponential kernel instead of an exponential one");
		mfSettings.maxEDs = textToInteger(parser.getOption("--max_ed", "Maximum number of eigendeformations", "-1"));


		int expert_section = parser.addSection("Expert options");

		padding = textToDouble(parser.getOption("--pad", "Apply Fourier padding to the cross-correlation images", "1"));
		whiten = !parser.checkOption("--no_whiten", "Do not not whiten the noise spectra");
		whiten_abs = parser.checkOption("--whiten_abs", "Divide by the square root of the power spectrum");
		hiPass_px = textToDouble(parser.getOption("--hp", "High-pass filter the cross-correlation images by this sigma", "-1"));
		sig2RampPower = textToDouble(parser.getOption("--srp", "Noise variance is divided by k^this during whitening", "0"));
		
		Log::readParams(parser);
		
		if (parser.checkForErrors()) std::exit(-1);
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	if (shiftOnly && do_motion)
	{
		REPORT_ERROR("The options --shift_only and --motion are mutually exclusive");
	}
}

void AlignProgram::run()
{
	Log::beginSection("Initialising");

		initialise();

		AberrationsCache aberrationsCache(particleSet.optTable, boxSize);
	
	Log::endSection();


	std::vector<int> tomoIndices = ParticleSet::enumerate(particles);

	processTomograms(tomoIndices, aberrationsCache, 1, true);

	finalise();
}

void AlignProgram::initialise()
{
	RefinementProgram::init();

	const int tc = particles.size();

	Log::beginSection("Configuration");
	Log::printBinaryChoice("Particle motion: ", do_motion, "considered", "not considered");
	Log::printBinaryChoice("Frame angles: ", mfSettings.constAngles, "constant", "variable");
	Log::printBinaryChoice("Frame shifts: ", mfSettings.constShifts, "constant", "variable");
	Log::printBinaryChoice("Particle positions: ", mfSettings.constParticles, "constant", "variable");
	Log::endSection();

	if (do_motion)
	{
		ZIO::makeDir(outDir + "/Trajectories");

		int tpc = particleSet.getTotalParticleNumber();

		if (particleSet.motionTrajectories.size() != tpc)
		{
			particleSet.motionTrajectories.resize(tpc);
		}

		for (int t = 0; t < tc; t++)
		{
			int pc = particles[t].size();
			if (pc == 0) continue;

			const int fc = tomogramSet.getFrameCount(t);

			for (int p = 0; p < pc; p++)
			{
				if (particleSet.motionTrajectories[particles[t][p].value].shifts_Ang.size() != fc)
				{
					particleSet.motionTrajectories[particles[t][p].value] = Trajectory(fc);
				}
			}
		}

		particleSet.hasMotion = true;
	}

	ZIO::ensureParentDir(getTempFilenameRoot(""));
}

void AlignProgram::finalise()
{
	const int tc = particles.size();

	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;

		readTempData(t);
	}

	if (do_motion)
	{
		Trajectory::write(
			particleSet.motionTrajectories, particleSet,
			particles, outDir + "motion.star");

		optimisationSet.trajectories = outDir+"motion.star";
	}

	if (!shiftOnly)
	{
		particleSet.write(outDir + "particles.star");
		optimisationSet.particles = outDir+"particles.star";
	}

	tomogramSet.write(outDir + "tomograms.star");
	optimisationSet.tomograms = outDir+"tomograms.star";

	optimisationSet.write(outDir+"optimisation_set.star");
}

void AlignProgram::processTomograms(
		const std::vector<int>& tomoIndices,
		const AberrationsCache& aberrationsCache,
		int verbosity,
		bool per_tomogram_progress)
{
	const int ttc = tomoIndices.size();

	if (verbosity > 0 && !per_tomogram_progress)
	{
		Log::beginProgress("Processing tomograms", ttc * num_iters);
	}

	for (int tt = 0; tt < ttc; tt++)
	{
		const int t = tomoIndices[tt];

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


		std::string tag = ZIO::itoa(t);
		std::string diagPrefix = outDir + "diag_" + tag;


		BufferedImage<float> frqWeight = computeFrequencyWeights(
			tomogram, whiten, sig2RampPower, hiPass_px, true, num_threads);

		if (diag)
		{
			frqWeight.write(diagPrefix + "_frq_weight.mrc");
		}


		std::vector<d4Matrix> projByTime(fc);

		for (int f = 0; f < fc; f++)
		{
			projByTime[f] = tomogram.projectionMatrices[tomogram.frameSequence[f]];
		}

		// motion estimation requires the CCs to be given in chronological order

		std::vector<int> frameSequence(fc);

		if (do_motion)
		{
			frameSequence = tomogram.frameSequence;
		}
		else
		{
			for (int f = 0; f < fc; f++)
			{
				frameSequence[f] = f;
			}
		}

		std::vector<BufferedImage<double>> CCs = Prediction::computeCroppedCCs(
				particleSet, particles[t], tomogram, aberrationsCache,
				referenceMap, frqWeight, frameSequence,
				range, true, num_threads, padding);


		BufferedImage<double> FCC;

		if (diag)
		{
			BufferedImage<double> FCC3 = FCC::compute3(
				particleSet, particles[t], tomogram, referenceMap.image_FS,
				true, num_threads);

			FCC3.write(diagPrefix + "_FCC3_initial.mrc");

			FCC = FCC::divide(FCC3);
			FCC.write(diagPrefix + "_FCC_initial.mrc");
		}


		const int progress_bar_offset = per_tomogram_progress? 0 : tt * num_iters;

		if (do_motion)
		{

			MotionFit motionFit(
				CCs, projByTime, particleSet, particles[t], referenceMap.image_FS,
				motParams, mfSettings, tomogram.centre,
				tomogram.getFrameDose(), tomogram.optics.pixelSize, padding,
				progress_bar_offset, num_threads);

			if (diag)
			{
				std::ofstream evDat(diagPrefix + "_deformation_eigenvalues.dat");

				for (int i = 0; i < motionFit.deformationLambda.size(); i++)
				{
					evDat << i << ' ' << motionFit.deformationLambda[i] << '\n';
				}
			}

			std::vector<double> initial(motionFit.getParamCount(), 0.0);


			if (verbosity > 0 && per_tomogram_progress)
			{
				Log::beginProgress("Performing optimisation", num_iters);
			}

			std::vector<double> opt = LBFGS::optimize(
				initial, motionFit, 1, num_iters, 1e-3, 1e-4);


			if (verbosity > 0 && per_tomogram_progress)
			{
				Log::endProgress();
			}

			std::vector<d4Matrix> projections = motionFit.getProjections(opt, tomogram.frameSequence);
			std::vector<d3Vector> positions = motionFit.getParticlePositions(opt);
			std::vector<Trajectory> trajectories = motionFit.exportTrajectories(
						opt, particleSet, tomogram.frameSequence);


			writeTempData(&trajectories, projections, positions, t);


			Mesh mesh8 = motionFit.visualiseTrajectories(opt, 8.0);
			mesh8.writePly(outDir + "Trajectories/" + tomogram.name + "_x8.ply");

			Mesh mesh1 = motionFit.visualiseTrajectories(opt, 1.0);
			mesh1.writePly(outDir + "Trajectories/" + tomogram.name + "_x1.ply");


			// Update the particle set in case an FCC is to be evaluated

			for (int p = 0; p < pc; p++)
			{
				const ParticleIndex pp = particles[t][p];

				particleSet.motionTrajectories[pp.value] = trajectories[p];
				particleSet.moveParticleTo(pp, positions[p]);
			}

			tomogram.projectionMatrices = projections;
		}
		else
		{
			std::vector<d4Matrix> projections = tomogram.projectionMatrices;
			std::vector<d3Vector> positions;

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

					CCsum += CCs[p];
				}

				if (verbosity > 0 && per_tomogram_progress)
				{
					Log::endProgress();
				}

				if (diag) CCsum.write(outDir + "CCsum_" + tag + ".mrc");

				d2Vector origin(padding*range, padding*range);

				for (int f = 0; f < fc; f++)
				{
					d2Vector opt = (Interpolation::quadraticMaxXY(CCsum.getSliceRef(f)) - origin)/padding;

					projections[f](0,3) += opt.x;
					projections[f](1,3) += opt.y;
				}

				if (!per_tomogram_progress)
				{
					Log::updateProgress(progress_bar_offset + num_iters);
				}

				positions.resize(pc);

				for (int p = 0; p < pc; p++)
				{
					positions[p] = particleSet.getPosition(particles[t][p]);
				}
			}
			else
			{
				ProtoAlignment protoAlignment(
					CCs, tomogram.projectionMatrices, particleSet, particles[t],
					referenceMap.image_FS,
					mfSettings.constParticles,
					mfSettings.constAngles,
					mfSettings.constShifts,
					range,
					tomogram.centre, progress_bar_offset,
					num_threads, padding);

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

				projections = protoAlignment.getProjections(opt);
				positions = protoAlignment.getParticlePositions(opt);
			}

			writeTempData(0, projections, positions, t);


			// Update the particle set in case an FCC is to be evaluated

			for (int p = 0; p < pc; p++)
			{
				const ParticleIndex pp = particles[t][p];
				particleSet.moveParticleTo(pp, positions[p]);
			}

			tomogram.projectionMatrices = projections;
		}

		if (diag)
		{

			BufferedImage<double> FCC3 = FCC::compute3(
				particleSet, particles[t], tomogram, referenceMap.image_FS,
				true, num_threads);

			FCC3.write(diagPrefix + "_FCC3_final.mrc");
			BufferedImage<double> FCC_new = FCC::divide(FCC3);
			FCC_new.write(diagPrefix + "_FCC_final.mrc");

			(FCC_new - FCC).write(diagPrefix + "_FCC_delta.mrc");
		}

		if (verbosity > 0 && per_tomogram_progress)
		{
			Log::endSection();
		}
	}

	if (verbosity > 0 && !per_tomogram_progress)
	{
		Log::endProgress();
	}
}


std::string AlignProgram::getTempFilenameRoot(const std::string& tomogram_name)
{
	return outDir + "temp/" + tomogram_name;
}

void AlignProgram::writeTempData(
		const std::vector<Trajectory>* traj,
		const std::vector<d4Matrix>& proj,
		const std::vector<d3Vector>& pos,
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


	if (do_motion && traj != 0)
	{
		Trajectory::write(*traj, particleSet, {particles[t]}, temp_filename_root + "_motion.star");
	}


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

void AlignProgram::readTempData(int t)
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


	if (do_motion)
	{
		std::ifstream ifs(temp_filename_root + "_motion.star");

		std::vector<MetaDataTable> mdts = MetaDataTable::readAll(ifs, pc+1);

		for (int p = 0; p < pc; p++)
		{
			MetaDataTable& mdt = mdts[p+1];

			int fc = mdt.numberOfObjects();

			d3Vector shift;

			for (int f = 0; f < fc; f++)
			{
				mdt.getValueSafely(EMDL_ORIENT_ORIGIN_X_ANGSTROM, shift.x, f);
				mdt.getValueSafely(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, shift.y, f);
				mdt.getValueSafely(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, shift.z, f);

				particleSet.motionTrajectories[particles[t][p].value].shifts_Ang[f] = shift;
			}
		}
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

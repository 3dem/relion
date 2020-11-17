#include "polish.h"
#include <src/ctf.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/motion/motion_fit.h>
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


PolishProgram::PolishProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	readParams(parser);
}

void PolishProgram::readParams(IOParser &parser)
{
	try
	{
		_readParams(parser);
				
		int defocus_section = parser.addSection("Alignment options");
		
		motParams.sig_vel = textToDouble(parser.getOption("--s_vel", "Velocity sigma [Å/dose]", "0.5"));
		motParams.sig_div = textToDouble(parser.getOption("--s_div", "Divergence sigma [Å]", "5000.0"));
		
		mfSettings.params_scaled_by_dose = !parser.checkOption("--abs_params", "Do not scale the sigmas by the dose");
		
		mfSettings.sqExpKernel = parser.checkOption("--sq_exp_ker", "Use a square-exponential kernel instead of an exponential one");
		mfSettings.maxEDs = textToInteger(parser.getOption("--max_ed", "Maximum number of eigendeformations", "-1"));
		
		range = textToInteger(parser.getOption("--r", "Max. shift allowed [Pixels]", "20"));
		padding = textToDouble(parser.getOption("--pad", "Apply Fourier padding to the cross-correlation images", "1"));
		whiten = !parser.checkOption("--no_whiten", "Do not not whiten the noise spectra");
		whiten_abs = parser.checkOption("--whiten_abs", "Divide by the square root of the power spectrum");
		hiPass_px = textToDouble(parser.getOption("--hp", "High-pass filter the cross-correlation images by this sigma", "-1"));
		sig2RampPower = textToDouble(parser.getOption("--srp", "Noise variance is divided by k^this during whitening", "0"));
		mfSettings.constParticles = parser.checkOption("--const_p", "Keep the particle positions constant");
		mfSettings.constAngles = parser.checkOption("--const_a", "Keep the frame angles constant");
		mfSettings.constShifts = parser.checkOption("--const_s", "Keep the frame shifts constant");
		num_iters = textToInteger(parser.getOption("--it", "Max. number of iterations", "10000"));
		
		outputShiftedCCs = parser.checkOption("--diag_CC", "Output shifted CCs (expensive)");
		
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

void PolishProgram::run()
{
	Log::beginSection("Initialising");

		initialise();

		AberrationsCache aberrationsCache(particleSet.optTable, boxSize);
	
	Log::endSection();


	std::vector<int> tomoIndices = ParticleSet::enumerate(particles);

	processTomograms(tomoIndices, aberrationsCache, 1, true);

	finalise();
}

void PolishProgram::initialise()
{
	RefinementProgram::init();
	int res = system(("mkdir -p " + outDir + "/Trajectories").c_str());

	const int tc = particles.size();

	Log::beginSection("Configuration");
	Log::printBinaryChoice("Frame angles: ", mfSettings.constAngles, "static", "variable");
	Log::printBinaryChoice("Frame shifts: ", mfSettings.constShifts, "static", "variable");
	Log::printBinaryChoice("Particle positions: ", mfSettings.constParticles, "static", "variable");
	Log::endSection();


	int tpc = particleSet.getTotalParticleNumber();

	if (particleSet.motionTrajectories.size() != tpc)
	{
		particleSet.motionTrajectories.resize(tpc);
	}

	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;

		Tomogram tomogram = tomogramSet.loadTomogram(t, false);
		const int fc = tomogram.frameCount;

		for (int p = 0; p < pc; p++)
		{
			if (particleSet.motionTrajectories[particles[t][p].value].shifts_Ang.size() != fc)
			{
				particleSet.motionTrajectories[particles[t][p].value] = Trajectory(fc);
			}
		}
	}

	particleSet.hasMotion = true;

	ZIO::ensureParentDir(getTempFilenameRoot(""));
}

void PolishProgram::finalise()
{
	const int tc = particles.size();

	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;

		readTempData(t);
	}

	Trajectory::write(
		particleSet.motionTrajectories, particleSet,
		particles, outDir + "motion.star");

	tomogramSet.write(outDir + "tomograms.star");
	particleSet.write(outDir + "particles.star");

	optimisationSet.particles = outDir+"particles.star";
	optimisationSet.tomograms = outDir+"tomograms.star";
	optimisationSet.trajectories = outDir+"motion.star";

	optimisationSet.write(outDir+"optimisation_set.star");
}

void PolishProgram::processTomograms(
		const std::vector<int>& tomoIndices,
		const AberrationsCache& aberrationsCache,
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


		std::vector<BufferedImage<double>> CCs = Prediction::computeCroppedCCs(
				particleSet, particles[t], tomogram, aberrationsCache,
				referenceMap, frqWeight, tomogram.frameSequence,
				range, true, num_threads, padding);


		MotionFit motionFit(
				CCs, projByTime, particleSet, particles[t], referenceMap.image_FS,
				motParams, mfSettings, tomogram.centre,
				tomogram.getFrameDose(), tomogram.optics.pixelSize, padding, 0, num_threads);


		BufferedImage<double> FCC3, FCC1, specCC;

		if (diag)
		{
			{
				std::ofstream evDat(diagPrefix + "_deformation_eigenvalues.dat");

				for (int i = 0; i < motionFit.deformationLambda.size(); i++)
				{
					evDat << i << ' ' << motionFit.deformationLambda[i] << '\n';
				}
			}

			FCC3 = FCC::compute3(
				particleSet, particles[t], tomogram, referenceMap.image_FS,
				true, num_threads);

			FCC3.write(diagPrefix + "_FCC3_initial.mrc");
			FCC1 = FCC::divide(FCC3);
			FCC1.write(diagPrefix + "_FCC_initial.mrc");

			if (outputShiftedCCs)
			{
				const int diam = CCs[0].xdim;

				BufferedImage<float> CCsum(diam, diam, fc);

				CCsum.fill(0.f);

				for (int p = 0; p < pc; p++)
				{
					CCsum += CCs[p];
				}

				CCsum.write(diagPrefix + "_CC_sum_" + tag + "_initial.mrc");


				const int d = CCsum.xdim;
				const int dh = d/2 + 1;

				specCC = BufferedImage<double>(dh,fc);
				specCC.fill(0.0);

				BufferedImage<fComplex> CCsumFS;

				for (int f = 0; f < fc; f++)
				{
					BufferedImage<float> CCsum_f = CCsum.getSliceRef(f);
					FFT::FourierTransform(CCsum_f, CCsumFS, FFT::Both);

					for (int y = 0; y < d; y++)
					for (int x = 0; x < dh; x++)
					{
						const double yy = y < d/2? y : y - d;
						const double rd = sqrt(x*x + yy*yy);
						const int r = (int) rd;

						const double mod = (1 - 2 * (x % 2)) * (1 - 2 * (y % 2));

						if (r < dh) specCC(r,f) += mod * CCsumFS(x,y).real;
					}
				}

				specCC.write(diagPrefix + "_specCC_initial.mrc");
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

		writeTempData(trajectories, projections, positions, t);



		Mesh mesh8 = motionFit.visualiseTrajectories(opt, 8.0);
		mesh8.writePly(outDir + "Trajectories/" + tomogram.name + "_x8.ply");

		Mesh mesh1 = motionFit.visualiseTrajectories(opt, 1.0);
		mesh1.writePly(outDir + "Trajectories/" + tomogram.name + "_x1.ply");

		if (diag)
		{
			for (int p = 0; p < pc; p++)
			{
				const ParticleIndex pp = particles[t][p];

				particleSet.motionTrajectories[pp.value] = trajectories[p];
				particleSet.moveParticleTo(pp, positions[p]);
			}

			tomogram.projectionMatrices = projections;

			FCC3 = FCC::compute3(
				particleSet, particles[t], tomogram, referenceMap.image_FS,
				true, num_threads);

			FCC3.write(diagPrefix + "_FCC3_final.mrc");
			BufferedImage<double> FCC1_new = FCC::divide(FCC3);
			FCC1_new.write(diagPrefix + "_FCC_final.mrc");

			(FCC1_new - FCC1).write(diagPrefix + "_FCC_delta.mrc");
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



std::string PolishProgram::getTempFilenameRoot(const std::string& tomogram_name)
{
	return outDir + "temp/" + tomogram_name;
}

void PolishProgram::writeTempData(
		const std::vector<Trajectory>& traj,
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


	Trajectory::write(traj, particleSet, {particles[t]}, temp_filename_root + "_motion.star");


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

void PolishProgram::readTempData(int t)
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

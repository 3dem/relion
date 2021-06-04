#include "align.h"
#include <src/ctf.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/motion/motion_fit.h>
#include <src/jaz/tomography/motion/proto_alignment.h>
#include <src/jaz/tomography/motion/trajectory_set.h>
#include <src/jaz/tomography/motion/modular_alignment/modular_alignment.h>
#include <src/jaz/tomography/motion/modular_alignment/GP_motion_model.h>
#include <src/jaz/tomography/motion/modular_alignment/no_motion_model.h>
#include <src/jaz/tomography/motion/modular_alignment/no_2D_deformation_model.h>
#include <src/jaz/tomography/motion/modular_alignment/spline_2D_deformation_model.h>
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
#include <iomanip>
#include <mpi.h>
#include <omp.h>

using namespace gravis;


AlignProgram::AlignProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv), debug(true)
{
}

void AlignProgram::run()
{
	parseInput();

	Log::beginSection("Initialising");

		initialise();

		AberrationsCache aberrationsCache(particleSet.optTable, boxSize);
	
	Log::endSection();


	std::vector<int> tomoIndices = ParticleSet::enumerate(particles);

	processTomograms(tomoIndices, aberrationsCache, true);

	finalise();
}

void AlignProgram::parseInput()
{
	IOParser parser;
	parser.setCommandLine(argc, argv);

	_readParams(parser);

	int alignment_section = parser.addSection("General alignment options");

	shiftOnly = parser.checkOption("--shift_only", "Only apply an optimal rigid shift to each frame (no iterative optimisation)");
	range = textToInteger(parser.getOption("--r", "Max. particle shift allowed [Pixels]", "20"));
	alignmentSettings.constParticles = parser.checkOption("--const_p", "Keep the particle positions constant");
	alignmentSettings.constAngles = parser.checkOption("--const_a", "Keep the frame angles constant");
	alignmentSettings.constShifts = parser.checkOption("--const_s", "Keep the frame shifts constant");
	alignmentSettings.rangeRegulariser = textToDouble(parser.getOption("--range_reg", "Value of the range regulariser", "0.0"));
	do_anisotropy = parser.checkOption("--aniso", "Assume an anisotropic projection model");
	per_tilt_anisotropy = parser.checkOption("--per_tilt_aniso", "Fit independent view anisotropy for each tilt image");
	num_iters = textToInteger(parser.getOption("--it", "Max. number of iterations", "50000"));


	int motion_section = parser.addSection("Motion estimation options");

	do_motion = parser.checkOption("--motion", "Estimate particle motion (expensive)");

	motionParameters.sig_vel = textToDouble(parser.getOption("--s_vel", "Velocity sigma [Å/dose]", "0.5"));
	motionParameters.sig_div = textToDouble(parser.getOption("--s_div", "Divergence sigma [Å]", "5000.0"));

	motionParameters.params_scaled_by_dose = !parser.checkOption("--abs_params", "Do not scale the sigmas by the dose");

	motionParameters.sqExpKernel = parser.checkOption("--sq_exp_ker", "Use a square-exponential kernel instead of an exponential one");
	motionParameters.maxEDs = textToInteger(parser.getOption("--max_ed", "Maximum number of eigendeformations", "-1"));


	int deformation_section = parser.addSection("Deformation estimation options");

	do_deformation = parser.checkOption("--deformation", "Estimate 2D deformations");
	deformationParameters.grid_width = textToInteger(parser.getOption("--def_w", "Number of horizontal sampling points for the deformation grid", "3"));
	deformationParameters.grid_height = textToInteger(parser.getOption("--def_h", "Number of vertical sampling points for the deformation grid", "3"));

	alignmentSettings.perFrame2DDeformation = parser.checkOption("--per_frame_deformation", "Model separate 2D deformations for all tilts");

	int expert_section = parser.addSection("Expert options");

	padding = textToDouble(parser.getOption("--pad", "Apply Fourier padding to the cross-correlation images", "1"));
	whiten = !parser.checkOption("--no_whiten", "Do not not whiten the noise spectra");
	whiten_abs = parser.checkOption("--whiten_abs", "Divide by the square root of the power spectrum");
	hiPass_px = textToDouble(parser.getOption("--hp", "High-pass filter the cross-correlation images by this sigma", "-1"));
	sig2RampPower = textToDouble(parser.getOption("--srp", "Noise variance is divided by k^this during whitening", "0"));

	Log::readParams(parser);

	if (shiftOnly && do_motion)
	{
		parser.reportError("ERROR: The options --shift_only and --motion are mutually exclusive");
	}

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
}

void AlignProgram::initialise()
{
	RefinementProgram::init();

	const int tc = particles.size();

	if (verbosity > 0)
	{
		Log::beginSection("Configuration");
		Log::printBinaryChoice("Particle motion: ", do_motion, "considered", "not considered");
		Log::printBinaryChoice("Frame angles: ", alignmentSettings.constAngles, "constant", "variable");
		Log::printBinaryChoice("Frame shifts: ", alignmentSettings.constShifts, "constant", "variable");
		Log::printBinaryChoice("Particle positions: ", alignmentSettings.constParticles, "constant", "variable");
		Log::endSection();
	}

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

	mergeLogFiles();

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


		BufferedImage<float> freqWeight = computeFrequencyWeights(
			tomogram, whiten, sig2RampPower, hiPass_px, false, num_threads);

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(boxSize, 1);

		if (diag)
		{
			freqWeight.write(diagPrefix + "_frq_weight.mrc");
		}


		std::vector<d4Matrix> projByTime(fc);

		for (int f = 0; f < fc; f++)
		{
			projByTime[f] = tomogram.projectionMatrices[tomogram.frameSequence[f]];
		}

		std::vector<BufferedImage<double>> CCs = Prediction::computeCroppedCCs(
				particleSet, particles[t], tomogram, aberrationsCache,
				referenceMap, freqWeight, doseWeights, tomogram.frameSequence,
				range, true, num_threads, padding, Prediction::OwnHalf,
				per_tomogram_progress && verbosity > 0);

		
		const int progress_bar_offset = per_tomogram_progress? 0 : tt * num_iters;

		if (do_motion)
		{
			GPMotionModel motionModel(
				particleSet, particles[t], tomogram,
				motionParameters,
				per_tomogram_progress && verbosity > 0);

			performAlignment(
				motionModel, CCs, projByTime, tomogram,
				t, progress_bar_offset, per_tomogram_progress);
		}
		else
		{

			if (!shiftOnly)
			{
				NoMotionModel noMotionModel;
				
				performAlignment(
					noMotionModel, CCs, projByTime, tomogram,
					t, progress_bar_offset, per_tomogram_progress);
			}
			else
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

				d2Vector origin(padding*range + 4, padding*range + 4);
				
				std::vector<d4Matrix> projections = tomogram.projectionMatrices;
				
				for (int f = 0; f < fc; f++)
				{
					d2Vector opt = (Interpolation::quadraticMaxXY(CCsum.getSliceRef(f)) - origin)/padding;
					
					const int ff = tomogram.frameSequence[f];

					projections[ff](0,3) += opt.x;
					projections[ff](1,3) += opt.y;
				}

				if (verbosity > 0 && !per_tomogram_progress)
				{
					Log::updateProgress(progress_bar_offset + num_iters);
				}
				
				std::vector<d3Vector> positions(pc);
				
				for (int p = 0; p < pc; p++)
				{
					positions[p] = particleSet.getPosition(particles[t][p]);
				}
				
				writeTempAlignmentData(projections, positions, t);
			}
			
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

void AlignProgram::writeTempAlignmentData(
		const std::vector<d4Matrix>& proj, 
		const std::vector<d3Vector>& pos, int t)
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

void AlignProgram::writeTempMotionData(
		const std::vector<Trajectory>& traj, 
		int t)
{
	const std::string tomoName = tomogramSet.getTomogramName(t);
	const std::string temp_filename_root = getTempFilenameRoot(tomoName);

	Trajectory::write(traj, particleSet, {particles[t]}, temp_filename_root + "_motion.star");
}

void AlignProgram::writeTempDeformationData(
		const std::vector<std::vector<double>>& def, 
		int t)
{
	const std::string tomoName = tomogramSet.getTomogramName(t);
	const std::string temp_filename_root = getTempFilenameRoot(tomoName);
	
	MetaDataTable temp_deformations;
		
	const int fc = def.size();

	for (int f = 0; f < fc; f++)
	{
		temp_deformations.addObject();
		temp_deformations.setValue(EMDL_TOMO_DEFORMATION_COEFFICIENTS, def[f], f);
	}

	temp_deformations.write(temp_filename_root + "_deformations.star");
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
	
	
	if (do_deformation)
	{
		MetaDataTable temp_deformations;
		temp_deformations.read(temp_filename_root + "_deformations.star");
		
		const i2Vector gridSize(deformationParameters.grid_width, deformationParameters.grid_height);
		
		std::vector<std::vector<double>> coeffs(fc);
		
		for (int f = 0; f < fc; f++)
		{
			coeffs[f] = temp_deformations.getDoubleVector(EMDL_TOMO_DEFORMATION_COEFFICIENTS, f);
		}
		
		tomogramSet.setDeformation(t, gridSize, coeffs);
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

void AlignProgram::mergeLogFiles()
{
	const int tc = tomogramSet.size();

	std::vector<FileName> eps_files;
	const std::vector<std::string> plot_names {"XY", "XZ", "YZ"};

	for (int t = 0; t < tc; t++)
	{
		const std::string tomo_name = tomogramSet.getTomogramName(t);

		if (do_motion)
		{
			for (int dim = 0; dim < 3; dim++)
			{
				const std::string fn = getTempFilenameRoot(tomo_name)
						+ "_tracks_" + plot_names[dim] + ".eps";

				if (ZIO::fileExists(fn))
				{
					eps_files.push_back(fn);
				}
			}
		}
	}

	if (eps_files.size() > 0)
	{
		joinMultipleEPSIntoSinglePDF(outDir + "logfile.pdf", eps_files);
	}
}

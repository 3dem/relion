#include "ctf_refinement.h"
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tilt_geometry.h>
#include <src/jaz/math/Zernike_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/optics/magnification_helper.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <iostream>
#include <src/time.h>
#include <mpi.h>

#define TIMING 0


using namespace gravis;
using namespace aberration;


CtfRefinementProgram::CtfRefinementProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
}

void CtfRefinementProgram::run()
{
	parseInput();

	Log::beginSection("Initialising");

		RefinementProgram::init();

		initTempDirectories();

		AberrationsCache aberrationsCache(particleSet.optTable, boxSize, particleSet.getOriginalPixelSize(0));

	Log::endSection();


	std::vector<int> tomoIndices = ParticleSet::enumerate(particles);

	processTomograms(tomoIndices, aberrationsCache, true);


	finalise();
}

void CtfRefinementProgram::parseInput()
{
	IOParser parser;
	parser.setCommandLine(argc, argv);

	_readParams(parser);


	int defocus_section = parser.addSection("Defocus refinement options");

	do_refine_defocus = parser.checkOption("--do_defocus", "Refine the (astigmatic) defocus.");
	do_regularise_defocus = parser.checkOption("--do_reg_defocus", "Regularise defocus estimation (i.e. require tilts to have similar defoci).");
	lambda_reg = textToDouble(parser.getOption("--lambda", "Defocus regularisation scale", "0.001"));
	if (lambda_reg < 0.0) lambda_reg = 0.0;

	minDelta = textToDouble(parser.getOption("--d0", "Min. defocus offset to test [Å]", "-3000"));
	maxDelta = textToDouble(parser.getOption("--d1", "Max. defocus offset to test [Å]", "3000"));
	deltaSteps = textToInteger(parser.getOption("--ds", "Number of defocus steps in-between", "100"));

	k_min_Ang = textToDouble(parser.getOption("--kmin", "Lowest spatial frequency to consider [Å]", "30"));
	do_reset_to_common = parser.checkOption("--reset_to_common", "Reset the CTFs of all tilts to a common one prior to local refinement.");


	int scale_section = parser.addSection("Scale estimation options");

	do_refine_scale = parser.checkOption("--do_scale", "Refine the contrast scale");
	bool per_frame_scale = parser.checkOption("--per_frame_scale", "Estimate the scale per frame (no Lambert fit)");
	bool per_tomogram_scale = parser.checkOption("--per_tomogram_scale", "Estimate the scale per tomogram (luminance may become unstable)");

	if (per_frame_scale && per_tomogram_scale)
	{
		parser.reportError("ERROR: The options --per_frame_scale and --per_tomogram_scale are mutually exclusive");
	}

	do_fit_Lambert_per_tomo = do_refine_scale && per_tomogram_scale;
	do_fit_Lambert_globally = do_refine_scale && !per_frame_scale && !per_tomogram_scale;


	int aberr_section = parser.addSection("Aberration refinement options");

	do_even_aberrations = parser.checkOption("--do_even_aberrations", "Refine even higher-order aberrations");
	do_odd_aberrations = parser.checkOption("--do_odd_aberrations", "Refine odd higher-order aberrations");
	do_refine_aberrations = do_odd_aberrations || do_even_aberrations;

	n_even = textToInteger(parser.getOption("--ne", "Maximal N for even aberrations", "4"));
	n_odd = textToInteger(parser.getOption("--no", "Maximal N for odd aberrations", "3"));

	int expert_section = parser.addSection("Expert options");

	min_frame = textToInteger(parser.getOption("--min_frame", "First frame to consider", "0"));
	max_frame = textToInteger(parser.getOption("--max_frame", "Last frame to consider", "-1"));
	freqCutoffFract = textToDouble(parser.getOption("--cutoff_fract", "Ignore shells for which the relative dose or frequency weight falls below this fraction of the average", "0.02"));

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
}

void CtfRefinementProgram::initTempDirectories()
{
	if (do_refine_defocus)
	{
		ZIO::ensureParentDir(getDefocusTempFilenameRoot(""));
	}

	if (do_refine_scale)
	{
		ZIO::ensureParentDir(getScaleTempFilenameRoot(""));
	}

	if (do_refine_aberrations)
	{
		if (do_even_aberrations)
		{
			ZIO::ensureParentDir(getEvenAberrationsTempFilename("", 0));
		}

		if (do_odd_aberrations)
		{
			ZIO::ensureParentDir(getOddAberrationsTempFilename("", 0));
		}
	}
}

void CtfRefinementProgram::processTomograms(
		const std::vector<int>& tomoIndices,
		const AberrationsCache& aberrationsCache,
		bool per_tomogram_progress)
{
	const int ttc = tomoIndices.size();

	if (verbosity > 0 && !per_tomogram_progress)
	{
		Log::beginProgress("Processing tomograms", ttc);
	}

	for (int tt = 0; tt < ttc; tt++)
	{
		if (verbosity > 0 && !per_tomogram_progress)
		{
			Log::updateProgress(tt);
		}

		abortIfNeeded();

		const int t = tomoIndices[tt];

		int pc = particles[t].size();
		if (pc == 0) continue;

		const std::string tomogram_name = tomogramSet.getTomogramName(t);
		const int gc = particleSet.numberOfOpticsGroups();

		if (only_do_unfinished &&
				defocusAlreadyDone(tomogram_name) &&
				scaleAlreadyDone(tomogram_name) &&
				aberrationsAlreadyDone(tomogram_name, gc) )
		{
			continue;
		}

		if (verbosity > 0 && per_tomogram_progress)
		{
			Log::beginSection(
				"Tomogram " + ZIO::itoa(tt + 1) + " / " + ZIO::itoa(ttc));

			Log::print("Loading");
		}

		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		tomogram.validateParticleOptics(particles[t], particleSet);

		const int fc = tomogram.frameCount;

		const double k_min_px = boxSize * tomogram.optics.pixelSize / k_min_Ang;


		particleSet.checkTrajectoryLengths(
				particles[t], fc, "CtfRefinementProgram::run");

		BufferedImage<float> freqWeights = computeFrequencyWeights(
			tomogram, true, 0.0, 0.0, false, num_threads);

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(boxSize, 1);


		BufferedImage<int> xRanges = findXRanges(freqWeights, doseWeights, freqCutoffFract);

		const int item_verbosity = per_tomogram_progress? verbosity : 0;

		abortIfNeeded();

		if (do_refine_defocus)
		{
			refineDefocus(
				t, tomogram, aberrationsCache, freqWeights, doseWeights, xRanges,
				k_min_px, item_verbosity);

			abortIfNeeded();
		}


		if (do_refine_scale)
		{
			updateScale(
				t, tomogram, aberrationsCache, freqWeights, doseWeights,
				item_verbosity);

			abortIfNeeded();
		}


		if (do_refine_aberrations)
		{
			updateAberrations(
				t, tomogram, aberrationsCache, freqWeights, doseWeights, xRanges,
				item_verbosity);

			abortIfNeeded();
		}


		if (verbosity > 0 && per_tomogram_progress)
		{
			Log::endSection();
		}

	} // all tomograms

	if (verbosity > 0 && !per_tomogram_progress)
	{
		Log::endProgress();
	}
}

void CtfRefinementProgram::finalise()
{
	if (!do_fit_Lambert_globally && !do_refine_aberrations)
	{
		Log::print("Merging temporary data");
	}
	else
	{
		Log::print("Merging temporary data and performing global fits");
	}

	if (do_refine_defocus)
	{
		collectDefocus();
	}

	if (do_refine_scale)
	{
		if (do_fit_Lambert_globally)
		{
			fitGlobalScale();
		}
		else
		{
			collectScale();
		}
	}

	if (do_refine_aberrations)
	{
		Tomogram tomogram = tomogramSet.loadTomogram(0, false);
		const double k_min_px = boxSize * tomogram.optics.pixelSize / k_min_Ang;

		fitAberrations(k_min_px);
	}

	Log::print("Writing output data");

	if (do_refine_aberrations)
	{
		particleSet.write(outDir + "particles.star");
		optimisationSet.particles = outDir + "particles.star";
	}

	if (do_refine_defocus || do_refine_scale)
	{
		tomogramSet.write(outDir + "tomograms.star");
		optimisationSet.tomograms = outDir + "tomograms.star";
	}

	optimisationSet.write(outDir + "optimisation_set.star");

	mergeLogFiles();

	std::cout << std::endl;
}



void CtfRefinementProgram::refineDefocus(
		int t,
		Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights,
		const BufferedImage<int>& xRanges,
		double k_min_px,
		int verbosity)
{
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int pc = particles[t].size();

	EvenData evenZero({0.0, 0.0, 0.0, 0.0, 0.0});
	OddData oddZero({0.0, dComplex(0.0, 0.0)});

	if (only_do_unfinished && defocusAlreadyDone(tomogram.name))
	{
		return;
	}

	if (verbosity > 0)
	{
		Log::beginSection("Refining defocus");
	}

	BufferedImage<EvenData> evenData(sh,s,fc);
	evenData.fill(evenZero);

	BufferedImage<OddData> oddData(sh,s,fc);
	oddData.fill(oddZero);

	if (do_regularise_defocus && do_reset_to_common)
	{
		// temporarily set all CTFs to that of the (chronologically) first frame:

		std::vector<int> chronoOrder = IndexSort<double>::sortIndices(tomogram.cumulativeDose);

		for (int f = 1; f < fc; f++)
		{
			tomogram.centralCTFs[chronoOrder[f]] = tomogram.centralCTFs[chronoOrder[0]];
		}
	}

	if (verbosity > 0)
	{
		Log::beginProgress("Accumulating defocus evidence", fc);
	}

	const int f0 = min_frame;
	const int f1 = max_frame > 0? max_frame : fc - 1;

	for (int f = f0; f <= f1; f++)
	{
		if (verbosity > 0)
		{
			Log::updateProgress(f);
		}

		std::vector<BufferedImage<EvenData>> evenData_thread(num_threads);
		std::vector<BufferedImage<OddData>> oddData_thread(num_threads);

		for (int th = 0; th < num_threads; th++)
		{
			evenData_thread[th] = BufferedImage<EvenData>(sh,s);
			evenData_thread[th].fill(evenZero);

			oddData_thread[th] = BufferedImage<OddData>(sh,s);
			oddData_thread[th].fill(oddZero);
		}

		#pragma omp parallel for num_threads(num_threads)
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();

			AberrationFit::considerParticle(
				particles[t][p], tomogram, referenceMap, particleSet,
				aberrationsCache, true, freqWeights, doseWeights, xRanges,
				f, f,
				evenData_thread[th], oddData_thread[th]);
		}

		for (int th = 0; th < num_threads; th++)
		{
			evenData.getSliceRef(f) += evenData_thread[th];
			oddData.getSliceRef(f)  += oddData_thread[th];
		}
	}

	if (verbosity > 0)
	{
		Log::endProgress();

		Log::print("Fitting");
	}

	const BufferedImage<double> dataTerm = evaluateDefocusRange(
			evenData, tomogram.optics.pixelSize, tomogram.centralCTFs,
			minDelta, maxDelta, deltaSteps, k_min_px);

	if (diag)
	{
		dataTerm.write(outDir + tomogram.name + "_dataTerm.mrc");
	}


	std::vector<d3Vector> astigmatism(fc);

	if (do_regularise_defocus)
	{
		int best_di = deltaSteps / 2;
		double minCost = std::numeric_limits<double>::max();

		for (int di = 0; di < deltaSteps; di++)
		{
			double dataTermSum = 0.0;

			for (int f = 0; f < fc; f++)
			{
				dataTermSum += dataTerm(f,di);
			}

			if (dataTermSum < minCost)
			{
				minCost = dataTermSum;
				best_di = di;
			}
		}

		const double deltaStep = (maxDelta - minDelta) / (double) (deltaSteps - 1);
		const double bestDeltaZ = minDelta + best_di * deltaStep;

		if (verbosity > 0)
		{
			Log::print("Refining astigmatic defocus");
		}

		EvenSolution solution = AberrationFit::solveEven(evenData);

		astigmatism = findMultiAstigmatism(
			solution, tomogram.centralCTFs, bestDeltaZ, tomogram.optics.pixelSize,
			lambda_reg, k_min_px);
	}
	else
	{
		for (int f = f0; f <= f1; f++)
		{
			int best_di = deltaSteps / 2;
			double minCost = std::numeric_limits<double>::max();

			for (int di = 0; di < deltaSteps; di++)
			{
				if (dataTerm(f,di) < minCost)
				{
					minCost = dataTerm(f,di);
					best_di = di;
				}
			}

			const double deltaStep = (maxDelta - minDelta) / (double) (deltaSteps - 1);
			const double bestDeltaZ = minDelta + best_di * deltaStep;

			RawImage<EvenData> evenDataSlice = evenData.getSliceRef(f);
			EvenSolution solutionSlice = AberrationFit::solveEven(evenDataSlice);

			astigmatism[f] = findAstigmatism(
				solutionSlice, tomogram.centralCTFs[f], bestDeltaZ,
				tomogram.optics.pixelSize,
				0.02, k_min_px);
		}
	}

	MetaDataTable tempTable;

	for (int f = 0; f < fc; f++)
	{
		CTF ctf = tomogram.centralCTFs[f];

		if (f >= f0 && f <= f1)
		{
			ctf.DeltafU = astigmatism[f][0];
			ctf.DeltafV = astigmatism[f][1];
			ctf.azimuthal_angle = astigmatism[f][2];
		}

		tempTable.addObject();
		tempTable.setValue(EMDL_CTF_DEFOCUSU, ctf.DeltafU, f);
		tempTable.setValue(EMDL_CTF_DEFOCUSV, ctf.DeltafV, f);
		tempTable.setValue(EMDL_CTF_DEFOCUS_ANGLE, ctf.azimuthal_angle, f);

		// Also store the new defoci in the tomogram set, so they can be
		// used for consecutive fits on the same node.
		// Note: this does not work when a run has been resumed.

		tomogramSet.setCtf(t, f, ctf);
		tomogram.centralCTFs[f] = ctf;
	}

	tempTable.write(getDefocusTempFilenameRoot(tomogram.name) + ".star");

	writeDefocusEps(tempTable, tomogram.name);

	if (verbosity > 0)
	{
		Log::endSection();
	}
}

void CtfRefinementProgram::updateScale(
		int t,
		Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights,
		int verbosity)
{
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int pc = particles[t].size();

	const int f0 = min_frame;
	const int f1 = max_frame > 0? max_frame : fc - 1;

	if (only_do_unfinished && scaleAlreadyDone(tomogram.name))
	{
		return;
	}

	if (verbosity > 0)
	{
		Log::beginSection("Refining scale");
	}

	std::vector<double> sum_prdObs_f(fc, 0.0);
	std::vector<double> sum_prdSqr_f(fc, 0.0);

	if (verbosity > 0)
	{
		Log::beginProgress("Accumulating scale evidence", pc);
	}

	for (int p = 0; p < pc; p++)
	{
		if (verbosity > 0)
		{
			Log::updateProgress(p);
		}

		const ParticleIndex part_id = particles[t][p];

		const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
					part_id, fc, tomogram.optics.pixelSize);
		
		const std::vector<bool> isVisible = tomogram.determineVisiblity(traj, s/2.0);

		#pragma omp parallel for num_threads(num_threads)
		for (int f = f0; f <= f1; f++)
		{
			if (!isVisible[f]) continue;
			
			d4Matrix projCut;

			BufferedImage<tComplex<float>> observation(sh,s);

			TomoExtraction::extractFrameAt3D_Fourier(
				tomogram.stack, f, s, 1.0, tomogram, traj[f],
				observation, projCut, 1, true);

			CTF ctf = tomogram.getCtf(f, particleSet.getPosition(part_id));

			const RawImage<float> doseSlice = doseWeights.getConstSliceRef(f);

			BufferedImage<fComplex> prediction = Prediction::predictModulated(
				part_id, particleSet, tomogram.projectionMatrices[f], s,
				ctf, tomogram.optics.pixelSize, aberrationsCache,
				referenceMap.image_FS,
				Prediction::OwnHalf,
				Prediction::AmplitudeModulated,
				&doseSlice,
				Prediction::CtfUnscaled);

			for (int y = 0; y < sh; y++)
			for (int x = 0; x < s;  x++)
			{
				const double xx = x;
				const double yy = y < s/2? y : y - s;
				const double r = sqrt(xx*xx + yy*yy);

				const fComplex obs = -observation(x,y);
				const fComplex prd =  prediction(x,y);

				const int ri = (int) r;

				if (ri < sh)
				{
					sum_prdObs_f[f] += freqWeights(x,y) * (prd.real * obs.real + prd.imag * obs.imag);
					sum_prdSqr_f[f] += freqWeights(x,y) * (prd.real * prd.real + prd.imag * prd.imag);
				}
			}

		} // all frames

	} // all particles

	if (verbosity > 0)
	{
		Log::endProgress();
	}


	std::vector<double> per_frame_scale(fc);

	for (int f = f0; f <= f1; f++)
	{
		per_frame_scale[f] = sum_prdObs_f[f] / sum_prdSqr_f[f];
	}

	if (diag)
	{
		std::ofstream scaleFile(outDir + tomogram.name + "_raw_scale.dat");

		for (int f = f0; f <= f1; f++)
		{
			scaleFile << f << ' ' << per_frame_scale[f] << '\n';
		}

		scaleFile.flush();
	}

	if (do_fit_Lambert_globally)
	{
		MetaDataTable tempFile;

		for (int f = 0; f < fc; f++)
		{
			tempFile.addObject();
			tempFile.setValue(EMDL_TOMO_TEMP_PRED_TIMES_OBS, sum_prdObs_f[f], f);
			tempFile.setValue(EMDL_TOMO_TEMP_PRED_SQUARED, sum_prdSqr_f[f], f);
		}

		const std::string tempFilename = getScaleTempFilenameRoot(tomogram.name) + ".star";
		tempFile.write(tempFilename);
	}
	else if (do_fit_Lambert_per_tomo)
	{
		double max_scale = 0.0;
		int max_scale_f = 0;

		for (int f = f0; f <= f1; f++)
		{
			if (per_frame_scale[f] > max_scale)
			{
				max_scale = per_frame_scale[f];
				max_scale_f = f;
			}
		}

		d3Vector max_scale_view;

		for (int i = 0; i < 3; i++)
		{
			max_scale_view[i] = tomogram.projectionMatrices[max_scale_f](2,i);
		}

		LambertFit blf(
			tomogram.projectionMatrices,
			sum_prdObs_f,
			sum_prdSqr_f);

		std::vector<double> initial {
			2.0 * max_scale,
			0.5,
			atan2(max_scale_view.dot(blf.tilt_p), max_scale_view.dot(blf.tilt_q))};

		std::vector<double> opt = NelderMead::optimize(initial, blf, 0.01, 0.0001, 10000);


		MetaDataTable tempGlobalTable, tempPerFrameTable;

		for (int f = 0; f < fc; f++)
		{
			CTF ctf = tomogram.centralCTFs[f];

			if (f >= f0 && f <= f1)
			{
				ctf.scale = blf.getScale(f, opt);
			}

			tempPerFrameTable.addObject();
			tempPerFrameTable.setValue(EMDL_CTF_SCALEFACTOR, ctf.scale, f);

			tomogramSet.setCtf(t, f, ctf);
			tomogram.centralCTFs[f] = ctf;
		}

		const double rel_thickness = -log(opt[1]);
		d3Vector ice_normal = sin(opt[2]) * blf.tilt_p + cos(opt[2]) * blf.tilt_q;
		ice_normal.normalize();

		tempGlobalTable.addObject();

		tempGlobalTable.setValue(EMDL_TOMO_RELATIVE_LUMINANCE, opt[0], 0);
		tempGlobalTable.setValue(EMDL_TOMO_RELATIVE_ICE_THICKNESS, rel_thickness, 0);
		tempGlobalTable.setValue(EMDL_TOMO_ICE_NORMAL_X, ice_normal.x, 0);
		tempGlobalTable.setValue(EMDL_TOMO_ICE_NORMAL_Y, ice_normal.y, 0);
		tempGlobalTable.setValue(EMDL_TOMO_ICE_NORMAL_Z, ice_normal.z, 0);

		{
			const std::string tempFilename = getScaleTempFilenameRoot(tomogram.name) + ".star";
			std::ofstream ofs(tempFilename);

			if (!ofs)
			{
				REPORT_ERROR("TomogramSet::write: unable to write to "+tempFilename);
			}

			tempGlobalTable.write(ofs);
			tempPerFrameTable.write(ofs);
		}

		// Also store the new scale in the tomogram set, so it can be
		// used for consecutive fits on the same MPI node.
		// Note: this does not work when a run has been resumed.

		tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_LUMINANCE, opt[0], t);
		tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_ICE_THICKNESS, rel_thickness, t);
		tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_X, ice_normal.x, t);
		tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Y, ice_normal.y, t);
		tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Z, ice_normal.z, t);


		if (diag)
		{
			std::ofstream scaleFile(outDir + tomogram.name + "_fitted_scale.dat");

			for (int f = f0; f <= f1; f++)
			{
				const double est = blf.getScale(f, opt);

				scaleFile << f << ' ' << est << '\n';
			}

			scaleFile.flush();
		}
	}
	else // !do_fit_Lambert_per_tomo && !do_fit_Lambert_globally
	{
		MetaDataTable tempPerFrameTable;

		for (int f = 0; f < fc; f++)
		{
			CTF ctf = tomogram.centralCTFs[f];

			if (f >= f0 && f <= f1)
			{
				ctf.scale = per_frame_scale[f];
			}

			tempPerFrameTable.addObject();
			tempPerFrameTable.setValue(EMDL_CTF_SCALEFACTOR, ctf.scale, f);

			// Also store the new scale in the tomogram set, so it can be
			// used for consecutive fits on the same MPI node.
			// Note: this does not work when a run has been resumed.

			tomogramSet.setCtf(t, f, ctf);
			tomogram.centralCTFs[f] = ctf;
		}

		tempPerFrameTable.write(getScaleTempFilenameRoot(tomogram.name) + ".star");
	}

	if (verbosity > 0)
	{
		Log::endSection();
	}
}

void CtfRefinementProgram::updateAberrations(
		int t,
		const Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights,
		const BufferedImage<int>& xRanges,
		int verbosity)
{
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int pc = particles[t].size();
	const int gc = particleSet.numberOfOpticsGroups();

	const int f0 = min_frame;
	const int f1 = max_frame > 0? max_frame : fc - 1;


	if (only_do_unfinished && aberrationsAlreadyDone(tomogram.name, gc))
	{
		return;
	}

	if (verbosity > 0)
	{
		Log::beginSection("Updating aberrations");
	}


	std::vector<BufferedImage<EvenData>> evenData_perGroup(gc);
	std::vector<BufferedImage<OddData>>  oddData_perGroup(gc);

	std::vector<std::vector<BufferedImage<EvenData>>> evenData_perGroup_perThread(gc);
	std::vector<std::vector<BufferedImage<OddData>>> oddData_perGroup_perThread(gc);

	EvenData evenZero({0.0, 0.0, 0.0, 0.0, 0.0});
	OddData oddZero({0.0, dComplex(0.0, 0.0)});

	for (int g = 0; g < gc; g++)
	{
		evenData_perGroup[g] = BufferedImage<EvenData>(sh,s);
		evenData_perGroup[g].fill(evenZero);

		oddData_perGroup[g] = BufferedImage<OddData>(sh,s);
		oddData_perGroup[g].fill(oddZero);

		evenData_perGroup_perThread[g] = std::vector<BufferedImage<EvenData>>(num_threads);
		oddData_perGroup_perThread[g]  = std::vector<BufferedImage<OddData>>(num_threads);
	}


	for (int g = 0; g < gc; g++)
	{
		for (int th = 0; th < num_threads; th++)
		{
			evenData_perGroup_perThread[g][th] = BufferedImage<EvenData>(sh,s);
			evenData_perGroup_perThread[g][th].fill(evenZero);

			oddData_perGroup_perThread[g][th] = BufferedImage<OddData>(sh,s);
			oddData_perGroup_perThread[g][th].fill(oddZero);
		}
	}


	if (verbosity > 0)
	{
		Log::beginProgress("Accumulating aberrations evidence", pc/num_threads);
	}

	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();

		if (th == 0 && verbosity > 0)
		{
			Log::updateProgress(p);
		}

		const int g = particleSet.getOpticsGroup(particles[t][p]);

		AberrationFit::considerParticle(
			particles[t][p], tomogram, referenceMap, particleSet,
			aberrationsCache, true, freqWeights, doseWeights, xRanges,
			f0, f1,
			evenData_perGroup_perThread[g][th],
			oddData_perGroup_perThread[g][th]);
	}

	if (verbosity > 0)
	{
		Log::endProgress();
	}


	for (int g = 0; g < gc; g++)
	{
		for (int th = 0; th < num_threads; th++)
		{
			evenData_perGroup[g] += evenData_perGroup_perThread[g][th];
			oddData_perGroup[g]  += oddData_perGroup_perThread[g][th];
		}

		EvenData::write(
			evenData_perGroup[g],
			getEvenAberrationsTempFilename(tomogram.name, g));

		OddData::write(
			oddData_perGroup[g],
			getOddAberrationsTempFilename(tomogram.name, g));
	}

	if (verbosity > 0)
	{
		Log::endSection();
	}
}



void CtfRefinementProgram::collectDefocus()
{
	const int tc = tomogramSet.size();

	for (int t = 0; t < tc; t++)
	{
		Tomogram tomogram = tomogramSet.loadTomogram(t, false);

		const std::string tempFilename = getDefocusTempFilenameRoot(tomogram.name) + ".star";

		if (!ZIO::fileExists(tempFilename)) continue;

		MetaDataTable tempTable;
		tempTable.read(tempFilename);

		const int fc = tomogram.frameCount;

		for (int f = 0; f < fc; f++)
		{
			CTF& ctf = tomogram.centralCTFs[f];

			ctf.DeltafU = tempTable.getDouble(EMDL_CTF_DEFOCUSU, f);
			ctf.DeltafV = tempTable.getDouble(EMDL_CTF_DEFOCUSV, f);
			ctf.azimuthal_angle = tempTable.getDouble(EMDL_CTF_DEFOCUS_ANGLE, f);

			ctf.initialise();

			tomogramSet.setCtf(t, f, ctf);
		}
	}
}

void CtfRefinementProgram::fitGlobalScale()
{
	std::vector<int> tomoIndices = ParticleSet::enumerateNonEmpty(particles);
	const int tc = tomoIndices.size();

	std::vector<std::vector<double>> all_sum_prdObs(tc), all_sum_prdSqr(tc);
	std::vector<std::vector<d4Matrix>> all_proj_matrices(tc);
	std::vector<double> all_fract_doses(tc);

	std::vector<bool> tomo_good(tc, false);


	for (int t = 0; t < tc; t++)
	{
		Tomogram tomogram = tomogramSet.loadTomogram(tomoIndices[t], false);

		const std::string tempFilename = getScaleTempFilenameRoot(tomogram.name) + ".star";
		MetaDataTable tempFile;

		if (!ZIO::fileExists(tempFilename))
		{
			Log::warn("Temporary file "+tempFilename+" not found.");
			continue;
		}

		tempFile.read(tempFilename);

		const int fc = tomogram.frameCount;

		if (tempFile.numberOfObjects() != fc)
		{
			REPORT_ERROR("CtfRefinementProgram::fitGlobalScale: temporary file "+tempFilename+" is corrupted.");
		}

		all_sum_prdObs[t] = std::vector<double>(fc);
		all_sum_prdSqr[t] = std::vector<double>(fc);

		for (int f = 0; f < fc; f++)
		{
			all_sum_prdObs[t][f] = tempFile.getDouble(EMDL_TOMO_TEMP_PRED_TIMES_OBS, f);
			all_sum_prdSqr[t][f] = tempFile.getDouble(EMDL_TOMO_TEMP_PRED_SQUARED, f);
		}

		all_proj_matrices[t] = tomogram.projectionMatrices;
		all_fract_doses[t] = tomogram.getFrameDose();

		tomo_good[t] = true;
	}

	std::vector<double> initial(2*tc + 1);

	double avg_max_scale = 0.0;
	double nonempty_tomos = 0.0;

	for (int t = 0; t < tc; t++)
	{
		const int fc = all_sum_prdObs[t].size();

		if (fc > 0)
		{
			std::vector<double> per_frame_scale(fc);

			for (int f = 0; f < fc; f++)
			{
				per_frame_scale[f] = all_sum_prdObs[t][f] / all_sum_prdSqr[t][f];
			}

			double max_scale = 0.0;
			int max_scale_f = 0;

			for (int f = 0; f < fc; f++)
			{
				if (per_frame_scale[f] > max_scale)
				{
					max_scale = per_frame_scale[f];
					max_scale_f = f;
				}
			}

			d3Vector max_scale_view;

			for (int i = 0; i < 3; i++)
			{
				max_scale_view[i] = all_proj_matrices[t][max_scale_f](2,i);
			}

			d3Matrix w2t = TiltGeometry::worldToTiltSpace(all_proj_matrices[t]);

			d3Vector tilt_p = d3Vector(w2t(0,0), w2t(0,1), w2t(0,2));
			d3Vector tilt_q = d3Vector(w2t(1,0), w2t(1,1), w2t(1,2));

			initial[2*t + 1] = 0.5;
			initial[2*t + 2] = atan2(max_scale_view.dot(tilt_p), max_scale_view.dot(tilt_q));

			avg_max_scale += max_scale;
			nonempty_tomos += 1.0;
		}
		else
		{
			initial[2*t + 1] = 0.5;
			initial[2*t + 2] = 0.0;
		}
	}

	avg_max_scale /= nonempty_tomos;

	initial[0] = 2.0 * avg_max_scale;

	MultiLambertFit mlf(
			all_proj_matrices, all_sum_prdObs,
			all_sum_prdSqr,
			all_fract_doses);

	std::vector<double> opt = LBFGS::optimize(initial, mlf, 0, 1000, 1e-5, 1e-4);


	for (int t = 0; t < tc; t++)
	{
		if (!tomo_good[t]) continue;

		const int tt = tomoIndices[t];

		Tomogram tomogram = tomogramSet.loadTomogram(tt, false);

		const int fc = tomogram.frameCount;

		const double lum    = opt[0];
		const double kappa  = opt[2*t + 1];
		const double phi    = opt[2*t + 2];

		const bool failed = kappa <= 0.0;

		if (failed)
		{
			Log::warn("Scale estimation for "+tomogram.name+" has failed.");
		}
		else
		{
			for (int f = 0; f < fc; f++)
			{
				const double est = mlf.getScale(t, f, opt);

				CTF ctf = tomogram.centralCTFs[f];

				ctf.scale = est;

				tomogramSet.setCtf(tt, f, ctf);
				tomogram.centralCTFs[f] = ctf;

				const double rel_thickness = -log(kappa);
				d3Vector ice_normal = sin(phi) * mlf.tilt_p[t] + cos(phi) * mlf.tilt_q[t];
				ice_normal.normalize();

				tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_LUMINANCE, lum, tt);
				tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_ICE_THICKNESS, rel_thickness, tt);
				tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_X, ice_normal.x, tt);
				tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Y, ice_normal.y, tt);
				tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Z, ice_normal.z, tt);
			}

			writeScaleEps(tomogramSet.tomogramTables[tt], tomogram.name);
		}
	}
}

void CtfRefinementProgram::collectScale()
{
	const int tc = tomogramSet.size();

	for (int t = 0; t < tc; t++)
	{
		Tomogram tomogram = tomogramSet.loadTomogram(t, false);
		std::string tempFilename = getScaleTempFilenameRoot(tomogram.name) + ".star";

		if (!ZIO::fileExists(tempFilename))
		{
			continue;
		}

		std::ifstream ifs(tempFilename);

		if (!ifs)
		{
			REPORT_ERROR_STR("CtfRefinementProgram::collectScale: Unable to read " << tempFilename);
		}

		const int fc = tomogram.frameCount;

		if (do_fit_Lambert_per_tomo)
		{
			std::vector<MetaDataTable> bothTables = MetaDataTable::readAll(ifs, 2);

			MetaDataTable& tempGlobalTable = bothTables[0];
			MetaDataTable& tempPerFrameTable = bothTables[1];

			tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_LUMINANCE,
				   tempGlobalTable.getDouble(EMDL_TOMO_RELATIVE_LUMINANCE, 0), t);

			tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_ICE_THICKNESS,
				   tempGlobalTable.getDouble(EMDL_TOMO_RELATIVE_ICE_THICKNESS, 0), t);

			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_X,
				   tempGlobalTable.getDouble(EMDL_TOMO_ICE_NORMAL_X, 0), t);

			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Y,
				   tempGlobalTable.getDouble(EMDL_TOMO_ICE_NORMAL_Y, 0), t);

			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Z,
				   tempGlobalTable.getDouble(EMDL_TOMO_ICE_NORMAL_Z, 0), t);


			for (int f = 0; f < fc; f++)
			{
				CTF& ctf = tomogram.centralCTFs[f];
				ctf.scale = tempPerFrameTable.getDouble(EMDL_CTF_SCALEFACTOR, f);
				ctf.initialise();

				tomogramSet.setCtf(t, f, ctf);
			}
		}
		else
		{
			MetaDataTable tempPerFrameTable;
			tempPerFrameTable.read(getScaleTempFilenameRoot(tomogram.name) + ".star");

			for (int f = 0; f < fc; f++)
			{
				CTF& ctf = tomogram.centralCTFs[f];
				ctf.scale = tempPerFrameTable.getDouble(EMDL_CTF_SCALEFACTOR, f);
				ctf.initialise();

				tomogramSet.setCtf(t, f, ctf);
			}
		}

		writeScaleEps(tomogramSet.tomogramTables[t], tomogram.name);
	}
}

void CtfRefinementProgram::fitAberrations(int k_min_px)
{
	const int s  = boxSize;
	const int sh = s/2 + 1;
	const int gc = particleSet.numberOfOpticsGroups();
	const int tc = tomogramSet.size();

	Tomogram tomogram0 = tomogramSet.loadTomogram(0, false);
	const double pixelSize = tomogram0.optics.pixelSize;

	for (int g = 0; g < gc; g++)
	{
		if (do_even_aberrations)
		{
			EvenData evenZero({0.0, 0.0, 0.0, 0.0, 0.0});
			BufferedImage<EvenData> even_data_sum(sh,s);
			even_data_sum.fill(evenZero);

			for (int t = 0; t < tc; t++)
			{
				const std::string fn = getEvenAberrationsTempFilename(tomogramSet.getTomogramName(t), g);

				if (ZIO::fileExists(fn + "_by.mrc"))
				{
					BufferedImage<EvenData> even = EvenData::read(fn);
					even_data_sum += even;
				}
			}

			for (int yy = 0; yy < s;  yy++)
			for (int xx = 0; xx < sh; xx++)
			{
				const double x = xx;
				const double y = yy < s/2? yy : yy - s;
				const double r = sqrt(x*x + y*y);

				if (r < k_min_px)
				{
					even_data_sum(xx,yy) *= 0.0;
				}
			}

			std::vector<double> initialEven(Zernike::numberOfEvenCoeffs(n_even), 0.0);

			std::vector<double> evenCoeffs = AberrationFit::solveAndFitEven(
				even_data_sum, n_even, initialEven,
				pixelSize, outDir + "group_" + ZIO::itoa(g+1) + "_", true);

			if (particleSet.optTable.containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS))
			{
				const std::vector<double> evenCoeffs0 = particleSet.optTable.getDoubleVector(
					EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, g);

				for (int i = 0; i < evenCoeffs.size(); i++)
				{
					if (i < evenCoeffs0.size())
					{
						evenCoeffs[i] += evenCoeffs0[i];
					}
				}
			}

			particleSet.optTable.setValue(
				EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, evenCoeffs, g);
		}

		if (do_odd_aberrations)
		{
			OddData oddZero({0.0, dComplex(0.0, 0.0)});
			BufferedImage<OddData> odd_data_sum(sh,s);
			odd_data_sum.fill(oddZero);

			for (int t = 0; t < tc; t++)
			{
				const std::string fn = getOddAberrationsTempFilename(tomogramSet.getTomogramName(t), g);

				if (ZIO::fileExists(fn + "_b_real.mrc"))
				{
					BufferedImage<OddData> odd = OddData::read(fn);
					odd_data_sum += odd;
				}
			}

			for (int yy = 0; yy < s;  yy++)
			for (int xx = 0; xx < sh; xx++)
			{
				const double x = xx;
				const double y = yy < s/2? yy : yy - s;
				const double r = sqrt(x*x + y*y);

				if (r < k_min_px)
				{
					odd_data_sum(xx,yy) *= 0.0;
				}
			}

			std::vector<double> initialOdd(Zernike::numberOfOddCoeffs(n_odd), 0.0);

			std::vector<double> oddCoeffs = AberrationFit::solveAndFitOdd(
				odd_data_sum, n_odd, initialOdd,
				pixelSize, outDir + "group_" + ZIO::itoa(g+1) + "_", true);

			if (particleSet.optTable.containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS))
			{
				const std::vector<double> oddCoeffs0 = particleSet.optTable.getDoubleVector(
					EMDL_IMAGE_ODD_ZERNIKE_COEFFS, g);

				for (int i = 0; i < oddCoeffs.size(); i++)
				{
					if (i < oddCoeffs0.size())
					{
						oddCoeffs[i] += oddCoeffs0[i];
					}
				}
			}

			particleSet.optTable.setValue(
				EMDL_IMAGE_ODD_ZERNIKE_COEFFS, oddCoeffs, g);
		}
	}
}



BufferedImage<double> CtfRefinementProgram::evaluateDefocusRange(
		const BufferedImage<EvenData>& evenData,
		double pixelSize,
		const std::vector<CTF>& ctfs,
		double minDefocus,
		double maxDefocus,
		int steps,
		double k_min_px)
{
	const int s  = evenData.ydim;
	const int sh = evenData.xdim;
	const int fc = evenData.zdim;

	const double as = s * pixelSize;
	const double eps = 1e-30;
	const double deltaStep = (maxDefocus - minDefocus) / (double) (steps - 1);

	const double k_min_sq = k_min_px * k_min_px;


	BufferedImage<double> out(fc,steps);


	for (int f = 0; f < fc; f++)
	{
		CTF ctfz = ctfs[f];

		for (int di = 0; di < steps; di++)
		{
			const double deltaZ = minDefocus + di * deltaStep;

			ctfz.DeltafU = ctfs[f].DeltafU + deltaZ;
			ctfz.DeltafV = ctfs[f].DeltafV + deltaZ;

			ctfz.initialise();

			double cost = 0.0;

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				const double xx = x;
				const double yy = y < s/2? y : y - s;
				const double r2 = xx * xx + yy * yy;

				if (r2 < s*s/4 && r2 >= k_min_sq)
				{
					EvenData d = evenData(x,y,f);

					d2Vector b(d.bx, d.by);
					d2Matrix A(d.Axx, d.Axy, d.Axy, d.Ayy);

					const double det = A(0,0) * A(1,1) - A(0,1) * A(1,0);

					if (std::abs(det) > eps)
					{
						d2Matrix Ai = A;
						Ai.invert();

						const d2Vector opt = Ai * b;

						const double gamma_0 = ctfs[f].getLowOrderGamma(xx/as, yy/as);
						const double gamma_z = ctfz.getLowOrderGamma(xx/as, yy/as);
						const double delta = gamma_z - gamma_0;

						const d2Vector dx = d2Vector(cos(delta), sin(delta)) - opt;

						cost += dx.dot(A * dx);
					}
				}
			}

			out(f,di) = cost / (s*s);
		}
	}

	return out;
}

gravis::d3Vector CtfRefinementProgram::findAstigmatism(
		const EvenSolution& solution,
		const CTF& referenceCtf,
		double initialDeltaZ,
		double pixelSize,
		double initialStep,
		double k_min_px)
{
	const int s = solution.optimum.ydim;
	const int sh = solution.optimum.xdim;
	const double K1 = PI * referenceCtf.lambda;

	BufferedImage<double> astigBasis(sh,s,3);

	const double as = s * pixelSize;

	for (int yi = 0; yi < s;  yi++)
	for (int xi = 0; xi < sh; xi++)
	{
		const double xf = xi;
		const double yf = (yi < s/2)? yi : yi - s;

		const double r2 = xf * xf + yf * yf;

		if (r2 > k_min_px * k_min_px)
		{
			const double xx = xf / as;
			const double yy = yf / as;

			astigBasis(xi,yi,0) = xx * xx + yy * yy;
			astigBasis(xi,yi,1) = xx * xx - yy * yy;
			astigBasis(xi,yi,2) = 2.0 * xx * yy;
		}
		else
		{
			astigBasis(xi,yi,0) = 0.0;
			astigBasis(xi,yi,1) = 0.0;
			astigBasis(xi,yi,2) = 0.0;
		}
	}

	ZernikeHelper::AnisoBasisOptimisation problem(
				solution.optimum, solution.weight, astigBasis, false);

	std::vector<double> nmOpt = NelderMead::optimize(
		{-initialDeltaZ * K1, 0.0, 0.0},
		problem, initialStep, 0.000001, 2000, 1.0, 2.0, 0.5, 0.5, false);

	const double dz = nmOpt[0];
	const double a1 = nmOpt[1];
	const double a2 = nmOpt[2];

	d2Matrix A_delta((dz+a1)/ K1,      a2 / K1,
						 a2 / K1,  (dz-a1)/ K1);

	d2Matrix A_ref(referenceCtf.getAxx(), referenceCtf.getAxy(),
				   referenceCtf.getAxy(), referenceCtf.getAyy());

	d2Matrix A_total = A_ref + A_delta;

	RFLOAT defocusU, defocusV, angleDeg;
	MagnificationHelper::matrixToPolar(
			A_total, defocusU, defocusV, angleDeg);

	return d3Vector(-defocusU, -defocusV, angleDeg);
}

std::vector<d3Vector> CtfRefinementProgram::findMultiAstigmatism(
		const EvenSolution& solution,
		const std::vector<CTF>& referenceCtfs,
		double initialDeltaZ,
		double pixelSize,
		double lambda_reg,
		double k_min_px)
{
	const int s  = solution.optimum.ydim;
	const int sh = solution.optimum.xdim;
	const int fc = solution.optimum.zdim;

	BufferedImage<double> astigBasis(sh,s,3);

	const double as = s * pixelSize;

	for (int yi = 0; yi < s;  yi++)
	for (int xi = 0; xi < sh; xi++)
	{
		const double xf = xi;
		const double yf = (yi < s/2)? yi : yi - s;

		const double r2 = xf * xf + yf * yf;

		if (r2 > k_min_px * k_min_px)
		{
			const double xx = xf / as;
			const double yy = yf / as;

			astigBasis(xi,yi,0) = xx * xx + yy * yy;
			astigBasis(xi,yi,1) = xx * xx - yy * yy;
			astigBasis(xi,yi,2) = 2.0 * xx * yy;
		}
		else
		{
			astigBasis(xi,yi,0) = 0.0;
			astigBasis(xi,yi,1) = 0.0;
			astigBasis(xi,yi,2) = 0.0;
		}
	}

	ZernikeHelper::MultiAnisoBasisOptimisation problem(
			solution.optimum, solution.weight, astigBasis, lambda_reg, false);

	std::vector<double> initial(3 * fc + 3, 0.0);

	for (int f = 0; f < fc+1; f++)
	{
		initial[3*f] = -initialDeltaZ * PI * referenceCtfs[0].lambda;
	}

	std::vector<double> nmOpt = LBFGS::optimize(initial, problem, false, 300, 1e-7, 1e-6);

	std::vector<d3Vector> out(fc);

	for (int f = 0; f < fc; f++)
	{
		const double dz = nmOpt[3*f + 3];
		const double a1 = nmOpt[3*f + 4];
		const double a2 = nmOpt[3*f + 5];

		const double K1 = PI * referenceCtfs[f].lambda;

		d2Matrix A_delta((dz+a1)/ K1,      a2 / K1,
							 a2 / K1,  (dz-a1)/ K1);

		d2Matrix A_ref(referenceCtfs[f].getAxx(), referenceCtfs[f].getAxy(),
					   referenceCtfs[f].getAxy(), referenceCtfs[f].getAyy());

		d2Matrix A_total = A_ref + A_delta;

		RFLOAT defocusU, defocusV, angleDeg;

		MagnificationHelper::matrixToPolar(
				A_total, defocusU, defocusV, angleDeg);

		out[f] = d3Vector(-defocusU, -defocusV, angleDeg);
	}

	return out;
}



std::string CtfRefinementProgram::getDefocusTempFilenameRoot(const std::string &tomogram_name)
{
	return outDir + "temp/defocus/" + tomogram_name;
}

std::string CtfRefinementProgram::getScaleTempFilenameRoot(const std::string& tomogram_name)
{
	return outDir + "temp/scale/" + tomogram_name;
}

std::string CtfRefinementProgram::getEvenAberrationsTempFilename(
		const std::string& tomogram_name, int opticsGroup)
{
	return outDir + "temp/aberrations/" + tomogram_name
			+ "_group_" + ZIO::itoa(opticsGroup) + "_even";
}

std::string CtfRefinementProgram::getOddAberrationsTempFilename(
		const std::string& tomogram_name, int opticsGroup)
{
	return outDir + "temp/aberrations/" + tomogram_name
			+ "_group_" + ZIO::itoa(opticsGroup) + "_odd";
}


bool CtfRefinementProgram::defocusAlreadyDone(const std::string &tomogram_name)
{
	return !do_refine_defocus || ZIO::fileExists(getDefocusTempFilenameRoot(tomogram_name) + ".star");
}

bool CtfRefinementProgram::scaleAlreadyDone(const std::string &tomogram_name)
{
	return !do_refine_scale || ZIO::fileExists(getScaleTempFilenameRoot(tomogram_name) + ".star");
}

bool CtfRefinementProgram::aberrationsAlreadyDone(
		const std::string &tomogram_name, int group_count)
{
	return
		(!do_even_aberrations ||
			ZIO::fileExists(getEvenAberrationsTempFilename(
					tomogram_name, group_count-1) + "_by.mrc") ) &&
		(!do_odd_aberrations  ||
			ZIO::fileExists(getOddAberrationsTempFilename(
								tomogram_name, group_count-1) + "_b_imag.mrc") );
}

void CtfRefinementProgram::writeDefocusEps(const MetaDataTable& table, const std::string& tomo_name)
{
	const std::string root_name = getDefocusTempFilenameRoot(tomo_name);

	CPlot2D plot2D(tomo_name + " Defocus");

	plot2D.SetXAxisSize(600);
	plot2D.SetYAxisSize(600);
	plot2D.SetDrawLegend(true);
	plot2D.SetFlipY(true);

	const int fc = table.numberOfObjects();

	for (int dim = 0; dim < 2; dim++)
	{
		CDataSet dataSet;

		dataSet.SetDrawMarker(false);
		dataSet.SetDrawLine(true);
		dataSet.SetLineWidth(1.0);
		dataSet.SetMarkerSize(10);
		dataSet.SetDatasetColor(0.1, 0.1, 0.1);

		for (int f = 0; f < fc; f++)
		{
			double delta_Z = 0.0;

			if (dim == 0)
			{
				table.getValue(EMDL_CTF_DEFOCUSU, delta_Z, f);
			}
			else
			{
				table.getValue(EMDL_CTF_DEFOCUSV, delta_Z, f);
			}

			dataSet.AddDataPoint(CDataPoint(f, delta_Z));
		}

		plot2D.AddDataSet(dataSet);
	}

	plot2D.SetXAxisTitle("frame");
	plot2D.SetYAxisTitle("defocus");
	plot2D.OutputPostScriptPlot(root_name + ".eps");


	std::ofstream defocusFile(root_name + ".dat");

	for (int f = 0; f < fc; f++)
	{
		defocusFile << f << ' ' << table.getDouble(EMDL_CTF_DEFOCUSU, f) << '\n';
	}

	defocusFile << '\n';

	for (int f = 0; f < fc; f++)
	{
		defocusFile << f << ' ' << table.getDouble(EMDL_CTF_DEFOCUSV, f) << '\n';
	}

	defocusFile.flush();
}

void CtfRefinementProgram::writeScaleEps(
		const MetaDataTable& table,
		const std::string& tomo_name)
{
	const std::string root_name = getScaleTempFilenameRoot(tomo_name);

	CPlot2D plot2D(tomo_name + " Contrast Scale");

	plot2D.SetXAxisSize(600);
	plot2D.SetYAxisSize(600);
	plot2D.SetDrawLegend(true);
	plot2D.SetFlipY(false);

	const int fc = table.numberOfObjects();

	CDataSet dataSet;

	dataSet.SetDrawMarker(false);
	dataSet.SetDrawLine(true);
	dataSet.SetLineWidth(1.0);
	dataSet.SetMarkerSize(10);
	dataSet.SetDatasetColor(0.1, 0.1, 0.1);

	for (int f = 0; f < fc; f++)
	{
		dataSet.AddDataPoint(CDataPoint(f, table.getDouble(EMDL_CTF_SCALEFACTOR, f)));
	}

	plot2D.AddDataSet(dataSet);

	plot2D.SetXAxisTitle("frame");
	plot2D.SetYAxisTitle("scale");
	plot2D.OutputPostScriptPlot(root_name + ".eps");


	std::ofstream scaleFile(root_name + ".dat");

	for (int f = 0; f < fc; f++)
	{
		scaleFile << f << ' ' << table.getDouble(EMDL_CTF_SCALEFACTOR, f) << '\n';
	}

	scaleFile.flush();
}

void CtfRefinementProgram::mergeLogFiles()
{
	const int tc = tomogramSet.size();

	std::vector<FileName> fn_eps(0);

	for (int t = 0; t < tc; t++)
	{
		const std::string tomo_name = tomogramSet.getTomogramName(t);

		const std::string fn_defocus = getDefocusTempFilenameRoot(tomo_name) + ".eps";

		if (do_refine_defocus && ZIO::fileExists(fn_defocus))
		{
			fn_eps.push_back(fn_defocus);
		}

		const std::string fn_scale = getScaleTempFilenameRoot(tomo_name) + ".eps";

		if (do_refine_scale && ZIO::fileExists(fn_scale))
		{
			fn_eps.push_back(fn_scale);
		}
	}

	if (fn_eps.size() > 0)
	{
		joinMultipleEPSIntoSinglePDF(outDir + "logfile.pdf", fn_eps);
	}
}

void CtfRefinementProgram::abortIfNeeded()
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
}



LambertFit::LambertFit(
	const std::vector<d4Matrix> &projections,
	const std::vector<double> &sum_prdObs,
	const std::vector<double> &sum_prdSqr)
:
  projections(projections),
  sum_prdObs(sum_prdObs),
  sum_prdSqr(sum_prdSqr)
{
	const int fc = projections.size();
	view_dir = std::vector<d3Vector>(fc);

	for (int f = 0; f < fc; f++)
	{
		for (int i = 0; i < 3; i++)
		{
			view_dir[f][i] = projections[f](2,i);
		}

		view_dir[f].normalize();
	}

	d3Matrix w2t = TiltGeometry::worldToTiltSpace(projections);

	tilt_p = d3Vector(w2t(0,0), w2t(0,1), w2t(0,2));
	tilt_q = d3Vector(w2t(1,0), w2t(1,1), w2t(1,2));
}

double LambertFit::f(const std::vector<double> &x, void *tempStorage) const
{
	const double a      = x[0];
	const double kappa  = x[1];
	const double phi    = x[2];

	d3Vector n = cos(phi) * tilt_q + sin(phi) * tilt_p;
	n.normalize();

	const int fc = projections.size();

	double out = 0.0;

	for (int f = 0; f < fc; f++)
	{
		const double cos_f = n.dot(view_dir[f]);
		const double est = a * pow(kappa, 1.0 / std::abs(cos_f));
		const double d = sum_prdSqr[f] * est - sum_prdObs[f];

		out += d * d;
	}

	return out;
}

double LambertFit::gradAndValue(const std::vector<double> &x, std::vector<double> &gradDest) const
{
	const double a      = x[0];
	const double kappa  = x[1];
	const double phi    = x[2];

	d3Vector n = cos(phi) * tilt_q + sin(phi) * tilt_p;
	n.normalize();

	const int fc = projections.size();

	double out = 0.0;

	for (int i = 0; i < 3; i++)
	{
		gradDest[i] = 0.0;
	}

	for (int f = 0; f < fc; f++)
	{
		const double cos_f = n.dot(view_dir[f]);
		const double est = a * pow(kappa, 1.0 / std::abs(cos_f));
		const double d = sum_prdSqr[f] * est - sum_prdObs[f];

		out += d * d;

		const double dout_dd = 2.0 * d;
		const double dd_dest = sum_prdSqr[f];
		const double dest_da = pow(kappa, 1.0 / std::abs(cos_f));
		const double dest_dkappa = a * pow(kappa, 1.0 / std::abs(cos_f)) / (kappa * cos_f);
		const d3Vector dn_dphi = cos(phi) * tilt_p - sin(phi) * tilt_q;
		const double dcosf_dphi = view_dir[f].dot(dn_dphi);
		const double dest_dphi = a * pow(kappa, 1.0 / std::abs(cos_f)) * log(kappa) * dcosf_dphi / (-cos_f * cos_f);

		gradDest[0] += dout_dd * dd_dest * dest_da;
		gradDest[1] += dout_dd * dd_dest * dest_dkappa;
		gradDest[2] += dout_dd * dd_dest * dest_dphi;
	}

	return out;
}

double LambertFit::getScale(int f, const std::vector<double> &x)
{
	const double a      = x[0];
	const double kappa  = x[1];
	const double phi    = x[2];

	d3Vector n = cos(phi) * tilt_q + sin(phi) * tilt_p;
	n.normalize();

	const double cos_f = n.dot(view_dir[f]);

	return a * pow(kappa, 1.0 / std::abs(cos_f));
}


MultiLambertFit::MultiLambertFit(
	const std::vector<std::vector<d4Matrix>>& projections,
	const std::vector<std::vector<double>>& sum_prdObs,
	const std::vector<std::vector<double>>& sum_prdSqr,
	const std::vector<double>& fractional_dose)
:
  projections(projections),
  sum_prdObs(sum_prdObs),
  sum_prdSqr(sum_prdSqr),
  fractional_dose(fractional_dose)
{
	const int tc = projections.size();

	view_dir = std::vector<std::vector<d3Vector>>(tc);
	tilt_p = std::vector<d3Vector>(tc);
	tilt_q = std::vector<d3Vector>(tc);

	for (int t = 0; t < tc; t++)
	{
		const int fc = projections[t].size();

		if (fc == 0) continue;

		view_dir[t].resize(fc);

		for (int f = 0; f < fc; f++)
		{
			for (int i = 0; i < 3; i++)
			{
				view_dir[t][f][i] = projections[t][f](2,i);
			}

			view_dir[t][f].normalize();
		}

		d3Matrix w2t = TiltGeometry::worldToTiltSpace(projections[t]);

		tilt_p[t] = d3Vector(w2t(0,0), w2t(0,1), w2t(0,2));
		tilt_q[t] = d3Vector(w2t(1,0), w2t(1,1), w2t(1,2));
	}
}

double MultiLambertFit::gradAndValue(const std::vector<double> &x, std::vector<double> &gradDest) const
{
	// @TODO: add fractional dose!

	const double a0 = x[0];
	const int tc = projections.size();

	const double kappa_min = 1e-3;

	double out = 0.0;

	for (int i = 0; i < 2*tc + 1; i++)
	{
		gradDest[i] = 0.0;
	}

	for (int t = 0; t < tc; t++)
	{
		const double kappa  = x[2*t + 1] > kappa_min? x[2*t + 1] : kappa_min;
		const double phi    = x[2*t + 2];

		d3Vector n = cos(phi) * tilt_q[t] + sin(phi) * tilt_p[t];
		n.normalize();

		const double a = a0 * fractional_dose[t];

		const int fc = projections[t].size();

		for (int f = 0; f < fc; f++)
		{
			const double cos_f = n.dot(view_dir[t][f]);
			const double est = a * pow(kappa, 1.0 / std::abs(cos_f));
			const double d = sum_prdSqr[t][f] * est - sum_prdObs[t][f];

			out += d * d;

			const double dout_dd = 2.0 * d;
			const double dd_dest = sum_prdSqr[t][f];
			const double dest_da = pow(kappa, 1.0 / std::abs(cos_f));
			const double dest_dkappa = a * pow(kappa, 1.0 / std::abs(cos_f)) / (kappa * cos_f);
			const d3Vector dn_dphi = cos(phi) * tilt_p[t] - sin(phi) * tilt_q[t];
			const double dcosf_dphi = view_dir[t][f].dot(dn_dphi);
			const double dest_dphi = a * pow(kappa, 1.0 / std::abs(cos_f)) * log(kappa) * dcosf_dphi / (-cos_f * cos_f);

			gradDest[0] += dout_dd * dd_dest * dest_da;

			gradDest[2*t + 1] += dout_dd * dd_dest * dest_dkappa;
			gradDest[2*t + 2] += dout_dd * dd_dest * dest_dphi;
		}
	}

	return out;
}

double MultiLambertFit::getScale(int t, int f, const std::vector<double> &x)
{
	const double a = x[0] * fractional_dose[t];

	const double kappa  = x[2*t + 1];
	const double phi    = x[2*t + 2];

	d3Vector n = cos(phi) * tilt_q[t] + sin(phi) * tilt_p[t];
	n.normalize();

	const double cos_f = n.dot(view_dir[t][f]);

	return a * pow(kappa, 1.0 / std::abs(cos_f));
}

void MultiLambertFit::report(int iteration, double cost, const std::vector<double> &x) const
{
	std::cout.precision(16);

	std::cout << iteration << ": " << cost << std::endl;
}

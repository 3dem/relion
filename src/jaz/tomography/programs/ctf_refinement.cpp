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

#define TIMING 0


using namespace gravis;
using namespace aberration;


CtfRefinementProgram::CtfRefinementProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	
	try
	{
		_readParams(parser);
		

		int defocus_section = parser.addSection("Defocus refinement options");

		do_refine_defocus = parser.checkOption("--do_defocus", "Refine the (astigmatic) defocus.");
		lambda_reg = textToDouble(parser.getOption("--lambda", "Defocus regularisation scale", "0.1"));

		minDelta = textToDouble(parser.getOption("--d0", "Min. defocus offset to test [Å]", "-3000"));
		maxDelta = textToDouble(parser.getOption("--d1", "Max. defocus offset to test [Å]", "3000"));
		deltaSteps = textToInteger(parser.getOption("--ds", "Number of defocus steps in-between", "100"));


		int scale_section = parser.addSection("Scale estimation options");

		do_refine_scale = parser.checkOption("--do_scale", "Refine the contrast scale");
		bool per_frame_scale = parser.checkOption("--per_frame_scale", "Estimate the scale per frame (no Beer-Lambert fit)");
		bool per_tomogram_scale = parser.checkOption("--per_tomogram_scale", "Estimate the scale per tomogram (luminance becomes unstable)");

		if (per_frame_scale && per_tomogram_scale)
		{
			parser.reportError("The options --per_frame_scale and --per_tomogram_scale are mutually exclusive");
		}

		do_fit_Beer_Lambert_per_tomo = do_refine_scale && per_tomogram_scale;
		do_fit_Beer_Lambert_globally = do_refine_scale && !per_frame_scale && !per_tomogram_scale;


		int aberr_section = parser.addSection("Aberration refinement options");

		do_refine_aberrations = parser.checkOption("--do_aberrations", "Refine higher-order aberrations");
		do_even_aberrations = !parser.checkOption("--no_even_aberrations", "Do not refine even aberrations");
		do_odd_aberrations = !parser.checkOption("--no_odd_aberrations", "Do not refine odd aberrations");
		n_even = textToInteger(parser.getOption("--ne", "Maximal N for even aberrations", "4"));
		n_odd = textToInteger(parser.getOption("--no", "Maximal N for odd aberrations", "3"));

		
		Log::readParams(parser);
		
		if (parser.checkForErrors())
		{
			std::exit(-2);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
}

void CtfRefinementProgram::run()
{
	Log::beginSection("Initialising");
	
	RefinementProgram::init();

	if (do_fit_Beer_Lambert_globally)
	{
		ZIO::makeDir(outDir + "temp/scale");
	}

	if (do_refine_aberrations)
	{
		ZIO::makeDir(outDir + "temp/aberrations");
	}

	const int tc = particles.size();

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize);
	
	Log::endSection();


	processTomograms(0, tc-1, aberrationsCache, 1);


	if (do_fit_Beer_Lambert_globally)
	{
		fitGlobalScale();
	}

	if (do_refine_aberrations)
	{
		fitAberrations();

		particleSet.write(outDir+"particles.star");
		optimisationSet.particles = outDir+"particles.star";
	}

	if (do_refine_defocus || do_refine_scale)
	{
		tomogramSet.write(outDir+"tomograms.star");
		optimisationSet.tomograms = outDir+"tomograms.star";
	}

	optimisationSet.write(outDir+"optimisation_set.star");
}

void CtfRefinementProgram::processTomograms(
		int first_t,
		int last_t,
		const AberrationsCache& aberrationsCache,
		int verbosity)
{
	for (int t = first_t; t <= last_t; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;

		if (verbosity > 0)
		{
			Log::beginSection(
					"Tomogram " + ZIO::itoa(t - first_t + 1)
					+ " / " + ZIO::itoa(last_t - first_t + 1));

			Log::print("Loading");
		}

		Tomogram tomogram = tomogramSet.loadTomogram(t, true);

		const int fc = tomogram.frameCount;


		particleSet.checkTrajectoryLengths(
				particles[t][0], pc, fc, "CtfRefinementProgram::run");

		BufferedImage<float> freqWeights = computeFrequencyWeights(
			tomogram, true, 0.0, 0.0, false, num_threads);

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(boxSize,1);


		if (do_refine_defocus)
		{
			refineDefocus(t, tomogram, aberrationsCache, freqWeights, doseWeights);
		}


		if (do_refine_scale)
		{
			updateScale(t, tomogram, aberrationsCache, freqWeights, doseWeights);
		}


		if (do_refine_aberrations)
		{
			updateAberrations(t, tomogram, aberrationsCache, freqWeights, doseWeights);
		}

		if (verbosity > 0)
		{
			Log::endSection();
		}

	} // all tomograms
}

void CtfRefinementProgram::refineDefocus(
		int t,
		Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights)
{
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int pc = particles[t].size();

	EvenData evenZero({0.0, 0.0, 0.0, 0.0, 0.0});
	OddData oddZero({0.0, dComplex(0.0, 0.0)});


	Log::beginSection("Refining defocus");

	BufferedImage<EvenData> evenData(sh,s,fc);
	evenData.fill(evenZero);

	BufferedImage<OddData> oddData(sh,s,fc);
	oddData.fill(oddZero);

	// temporarily set all CTFs to that of the (chronologically) first frame:

	std::vector<int> chronoOrder = IndexSort<double>::sortIndices(tomogram.cumulativeDose);

	for (int f = 1; f < fc; f++)
	{
		tomogram.centralCTFs[chronoOrder[f]] = tomogram.centralCTFs[chronoOrder[0]];
	}



	Log::beginProgress("Accumulating defocus evidence", pc * fc / num_threads);

	for (int f = 0; f < fc; f++)
	{
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

			if (th == 0)
			{
				Log::updateProgress(pc * f / num_threads + p);
			}

			AberrationFit::considerParticle(
				particles[t][p], tomogram, referenceMap, particleSet,
				aberrationsCache, true, freqWeights, doseWeights,
				f, f,
				evenData_thread[th], oddData_thread[th]);
		}

		for (int th = 0; th < num_threads; th++)
		{
			evenData.getSliceRef(f) += evenData_thread[th];
			oddData.getSliceRef(f)  += oddData_thread[th];
		}
	}

	Log::endProgress();

	Log::print("Fitting");


	const BufferedImage<double> dataTerm = evaluateDefocusRange(
			evenData, tomogram.optics.pixelSize, tomogram.centralCTFs,
			minDelta, maxDelta, deltaSteps);

	if (diag)
	{
		dataTerm.write(outDir + tomogram.name + "_dataTerm.mrc");
	}


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

	Log::print("Refining astigmatic defocus");

	EvenSolution solution = AberrationFit::solveEven(evenData);

	std::vector<d3Vector> astig = findMultiAstigmatism(
		solution, tomogram.centralCTFs, bestDeltaZ, tomogram.optics.pixelSize, lambda_reg);

	for (int f = 0; f < fc; f++)
	{
		CTF ctf0 = tomogram.centralCTFs[f];
		CTF ctf_dz = ctf0;

		ctf_dz.DeltafU = ctf0.DeltafU + bestDeltaZ;
		ctf_dz.DeltafV = ctf0.DeltafV + bestDeltaZ;

		CTF ctf1 = ctf_dz;

		ctf1.DeltafU = astig[f][0];
		ctf1.DeltafV = astig[f][1];
		ctf1.azimuthal_angle = astig[f][2];

		tomogramSet.setCtf(t, f, ctf1);
		tomogram.centralCTFs[f] = ctf1;
	}

	if (diag)
	{
		std::ofstream meanDefocus(outDir + tomogram.name + "_mean_defocus.dat");

		for (int f = 0; f < fc; f++)
		{
			meanDefocus << f << ' ' << (astig[f][0] + astig[f][1]) / 2.0 << '\n';
		}

		meanDefocus.flush();
	}

	Log::endSection();
}

void CtfRefinementProgram::updateScale(
		int t,
		Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights)
{
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int pc = particles[t].size();

	Log::beginSection("Refining scale");

	std::vector<double> sum_prdObs_f(fc, 0.0);
	std::vector<double> sum_prdSqr_f(fc, 0.0);

	Log::beginProgress("Accumulating scale evidence", pc);

	for (int p = 0; p < pc; p++)
	{
		Log::updateProgress(p);

		const ParticleIndex part_id = particles[t][p];

		const std::vector<d3Vector> traj = particleSet.getTrajectoryInPixels(
					part_id, fc, tomogram.optics.pixelSize);

		#pragma omp parallel for num_threads(num_threads)
		for (int f = 0; f < fc; f++)
		{
			d4Matrix projCut;

			BufferedImage<tComplex<float>> observation(sh,s);

			TomoExtraction::extractFrameAt3D_Fourier(
				tomogram.stack, f, s, 1.0, tomogram.projectionMatrices[f], traj[f],
				observation, projCut, 1, true);

			CTF ctf = tomogram.getCtf(f, particleSet.getPosition(part_id));

			BufferedImage<fComplex> prediction = Prediction::predictModulated(
				part_id, particleSet, tomogram.projectionMatrices[f], s,
				ctf, tomogram.optics.pixelSize, aberrationsCache,
				referenceMap.image_FS,
				Prediction::OwnHalf,
				Prediction::AmplitudeModulated);

			for (int y = 0; y < sh; y++)
			for (int x = 0; x < s;  x++)
			{
				const double xx = x;
				const double yy = y < s/2? y : y - s;
				const double r = sqrt(xx*xx + yy*yy);

				const fComplex obs = -observation(x,y);
				const fComplex prd =  doseWeights(x,y,f) * prediction(x,y) / ctf.scale;

				const int ri = (int) r;

				if (ri < sh)
				{
					sum_prdObs_f[f] += freqWeights(x,y) * (prd.real * obs.real + prd.imag * obs.imag);
					sum_prdSqr_f[f] += freqWeights(x,y) * (prd.real * prd.real + prd.imag * prd.imag);
				}
			}

		} // all frames

	} // all particles

	Log::endProgress();


	std::vector<double> per_frame_scale(fc);

	for (int f = 0; f < fc; f++)
	{
		per_frame_scale[f] = sum_prdObs_f[f] / sum_prdSqr_f[f];
	}

	if (diag)
	{
		std::ofstream scaleFile(outDir + tomogram.name + "_raw_scale.dat");

		for (int f = 0; f < fc; f++)
		{
			scaleFile << f << ' ' << per_frame_scale[f] << '\n';
		}

		scaleFile.flush();
	}

	if (do_fit_Beer_Lambert_globally)
	{
		std::ofstream sum_prdObs_file(outDir + "temp/scale/" + tomogram.name + "_sum_predObs.dat");
		std::ofstream sum_prdSqr_file(outDir + "temp/scale/" + tomogram.name + "_sum_predSqr.dat");

		for (int f = 0; f < fc; f++)
		{
			sum_prdObs_file << sum_prdObs_f[f] << ' ';
			sum_prdSqr_file << sum_prdSqr_f[f] << ' ';
		}

		sum_prdObs_file << '\n';
		sum_prdSqr_file << '\n';
	}
	else if (do_fit_Beer_Lambert_per_tomo)
	{
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
			max_scale_view[i] = tomogram.projectionMatrices[max_scale_f](2,i);
		}

		BeerLambertFit blf(
			tomogram.projectionMatrices,
			sum_prdObs_f,
			sum_prdSqr_f);

		std::vector<double> initial {
			2.0 * max_scale,
			0.5,
			atan2(max_scale_view.dot(blf.tilt_p), max_scale_view.dot(blf.tilt_q))};

		std::vector<double> opt = NelderMead::optimize(initial, blf, 0.01, 0.0001, 10000);

		for (int f = 0; f < fc; f++)
		{
			const double est = blf.getScale(f, opt);

			CTF ctf = tomogram.centralCTFs[f];

			ctf.scale = est;

			tomogramSet.setCtf(t, f, ctf);
			tomogram.centralCTFs[f] = ctf;

			const double rel_thickness = -log(opt[1]);
			d3Vector ice_normal = sin(opt[2]) * blf.tilt_p + cos(opt[2]) * blf.tilt_q;
			ice_normal.normalize();

			tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_LUMINANCE, opt[0], t);
			tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_ICE_THICKNESS, rel_thickness, t);
			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_X, ice_normal.x, t);
			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Y, ice_normal.y, t);
			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Z, ice_normal.z, t);
		}

		if (diag)
		{
			std::ofstream scaleFile(outDir + tomogram.name + "_fitted_scale.dat");

			for (int f = 0; f < fc; f++)
			{
				const double est = blf.getScale(f, opt);

				scaleFile << f << ' ' << est << '\n';
			}

			scaleFile.flush();
		}
	}
	else // !do_fit_Beer_Lambert_per_tomo && !do_fit_Beer_Lambert_globally
	{
		for (int f = 0; f < fc; f++)
		{
			CTF ctf = tomogram.centralCTFs[f];

			ctf.scale = per_frame_scale[f];

			tomogramSet.setCtf(t, f, ctf);
			tomogram.centralCTFs[f] = ctf;
		}
	}

	Log::endSection();
}

void CtfRefinementProgram::updateAberrations(
		int t,
		const Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		const BufferedImage<float>& freqWeights,
		const BufferedImage<float>& doseWeights)
{
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const int pc = particles[t].size();
	const int gc = particleSet.numberOfOpticsGroups();


	Log::beginSection("Updating aberrations");


	std::vector<BufferedImage<EvenData>> evenData_perGroup(gc);
	std::vector<BufferedImage<OddData>>  oddData_perGroup(gc);

	std::vector<std::vector<BufferedImage<EvenData>>> evenData_perGroup_perThread(num_threads);
	std::vector<std::vector<BufferedImage<OddData>>> oddData_perGroup_perThread(num_threads);

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


	Log::beginProgress("Accumulating aberrations evidence", pc/num_threads);

	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();

		if (th == 0)
		{
			Log::updateProgress(p/num_threads);
		}

		const int g = particleSet.getOpticsGroup(particles[t][p]);

		AberrationFit::considerParticle(
			particles[t][p], tomogram, referenceMap, particleSet,
			aberrationsCache, true, freqWeights, doseWeights,
			0, fc - 1,
			evenData_perGroup_perThread[g][th],
			oddData_perGroup_perThread[g][th]);
	}

	Log::endProgress();


	for (int g = 0; g < gc; g++)
	{
		for (int th = 0; th < num_threads; th++)
		{
			evenData_perGroup[g] += evenData_perGroup_perThread[g][th];
			oddData_perGroup[g]  += oddData_perGroup_perThread[g][th];
		}

		EvenData::write(
			evenData_perGroup[g], outDir + "temp/aberrations/" +
			tomogram.name + "_group_" + ZIO::itoa(g) + "_even");

		OddData::write(
			oddData_perGroup[g], outDir + "temp/aberrations/" +
			tomogram.name + "_group_" + ZIO::itoa(g) + "_odd");
	}

	Log::endSection();
}

void CtfRefinementProgram::fitAberrations()
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
				Tomogram tomogram = tomogramSet.loadTomogram(t, false);

				BufferedImage<EvenData> even = EvenData::read(
						outDir + "temp/aberrations/" +
						tomogram.name + "_group_" + ZIO::itoa(g) + "_even");

				even_data_sum += even;
			}

			std::vector<double> initialEven(Zernike::numberOfEvenCoeffs(n_even), 0.0);

			std::vector<double> evenCoeffs = AberrationFit::solveAndFitEven(
				even_data_sum, n_even, initialEven,
				pixelSize, outDir + "group_" + ZIO::itoa(g+1) + "_", true);

			if (particleSet.optTable.labelExists(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS))
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
				Tomogram tomogram = tomogramSet.loadTomogram(t, false);

				BufferedImage<OddData> odd = OddData::read(
						outDir + "temp/aberrations/" +
						tomogram.name + "_group_" + ZIO::itoa(g) + "_odd");

				odd_data_sum += odd;
			}

			std::vector<double> initialOdd(Zernike::numberOfOddCoeffs(n_odd), 0.0);

			std::vector<double> oddCoeffs = AberrationFit::solveAndFitOdd(
				odd_data_sum, n_odd, initialOdd,
				pixelSize, outDir + "group_" + ZIO::itoa(g+1) + "_", true);

			if (particleSet.optTable.labelExists(EMDL_IMAGE_ODD_ZERNIKE_COEFFS))
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

void CtfRefinementProgram::fitGlobalScale()
{
	const int tc = tomogramSet.size();

	std::vector<std::vector<double>> all_sum_prdObs(tc), all_sum_prdSqr(tc);
	std::vector<std::vector<d4Matrix>> all_proj_matrices(tc);
	std::vector<double> all_fract_doses(tc);

	for (int t = 0; t < tc; t++)
	{
		Tomogram tomogram = tomogramSet.loadTomogram(t, false);

		std::ifstream sum_prdObs_file(outDir + "temp/scale/" + tomogram.name + "_sum_predObs.dat");
		std::ifstream sum_prdSqr_file(outDir + "temp/scale/" + tomogram.name + "_sum_predSqr.dat");

		if (sum_prdObs_file && sum_prdSqr_file)
		{
			const int fc = tomogram.frameCount;

			all_sum_prdObs[t] = std::vector<double>(fc);
			all_sum_prdSqr[t] = std::vector<double>(fc);

			for (int f = 0; f < fc; f++)
			{
				sum_prdObs_file >> all_sum_prdObs[t][f];
				sum_prdSqr_file >> all_sum_prdSqr[t][f];
			}

			all_proj_matrices[t] = tomogram.projectionMatrices;
		}

		all_fract_doses[t] = tomogram.getFrameDose();
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

	MultiBeerLambertFit mblf(
			all_proj_matrices, all_sum_prdObs,
			all_sum_prdSqr,
			all_fract_doses);

	std::vector<double> opt = LBFGS::optimize(initial, mblf, 0, 1000, 1e-5, 1e-4);

	for (int t = 0; t < tc; t++)
	{
		Tomogram tomogram = tomogramSet.loadTomogram(t, false);

		const int fc = tomogram.frameCount;

		const double lum    = opt[0];
		const double kappa  = opt[2*t + 1];
		const double phi    = opt[2*t + 2];

		for (int f = 0; f < fc; f++)
		{
			const double est = mblf.getScale(t, f, opt);

			CTF ctf = tomogram.centralCTFs[f];

			ctf.scale = est;

			tomogramSet.setCtf(t, f, ctf);
			tomogram.centralCTFs[f] = ctf;

			const double rel_thickness = -log(kappa);
			d3Vector ice_normal = sin(phi) * mblf.tilt_p[t] + cos(phi) * mblf.tilt_q[t];
			ice_normal.normalize();

			tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_LUMINANCE, lum, t);
			tomogramSet.globalTable.setValue(EMDL_TOMO_RELATIVE_ICE_THICKNESS, rel_thickness, t);
			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_X, ice_normal.x, t);
			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Y, ice_normal.y, t);
			tomogramSet.globalTable.setValue(EMDL_TOMO_ICE_NORMAL_Z, ice_normal.z, t);
		}

		if (diag)
		{
			std::ofstream scaleFile(outDir + tomogram.name + "_multifitted_scale.dat");

			for (int f = 0; f < fc; f++)
			{
				const double est = mblf.getScale(t, f, opt);

				scaleFile << f << ' ' << est << '\n';
			}

			scaleFile.flush();
		}
	}
}

BufferedImage<double> CtfRefinementProgram::evaluateDefocusRange(
		const BufferedImage<EvenData>& evenData,
		double pixelSize,
		const std::vector<CTF>& ctfs,
		double minDefocus,
		double maxDefocus,
		int steps)
{
	const int s  = evenData.ydim;
	const int sh = evenData.xdim;
	const int fc = evenData.zdim;

	const double as = s * pixelSize;
	const double eps = 1e-30;
	const double deltaStep = (maxDefocus - minDefocus) / (double) (steps - 1);


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

				if (r2 < s*s/4)
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
		double initialStep)
{
	const int s = solution.optimum.ydim;
	const int sh = solution.optimum.xdim;
	const double K1 = PI * referenceCtf.lambda;

	BufferedImage<double> astigBasis(sh,s,3);

	const double as = s * pixelSize;

	for (int yi = 0; yi < s;  yi++)
	for (int xi = 0; xi < sh; xi++)
	{
		const double xx = xi/as;
		const double yy = (yi < s/2)? yi/as : (yi - s)/as;

		astigBasis(xi,yi,0) = xx * xx + yy * yy;
		astigBasis(xi,yi,1) = xx * xx - yy * yy;
		astigBasis(xi,yi,2) = 2.0 * xx * yy;
	}

	ZernikeHelper::AnisoBasisOptimisation problem(
				solution.optimum, solution.weight, astigBasis, false);

	std::vector<double> nmOpt = NelderMead::optimize(
		{-initialDeltaZ * K1, 0.0, 0.0},
		problem, initialStep, 0.000001, 2000, 1.0, 2.0, 0.5, 0.5, false);

	std::cout.precision(16);

	for (int i = 0; i < nmOpt.size(); i++)
	{
		std::cout << nmOpt[i] << "  ";
	}

	std::cout << "\n";

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
		double lambda_reg)
{
	const int s  = solution.optimum.ydim;
	const int sh = solution.optimum.xdim;
	const int fc = solution.optimum.zdim;

	BufferedImage<double> astigBasis(sh,s,3);

	const double as = s * pixelSize;

	for (int yi = 0; yi < s;  yi++)
	for (int xi = 0; xi < sh; xi++)
	{
		const double xx = xi/as;
		const double yy = (yi < s/2)? yi/as : (yi - s)/as;

		astigBasis(xi,yi,0) = xx * xx + yy * yy;
		astigBasis(xi,yi,1) = xx * xx - yy * yy;
		astigBasis(xi,yi,2) = 2.0 * xx * yy;
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



BeerLambertFit::BeerLambertFit(
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

double BeerLambertFit::f(const std::vector<double> &x, void *tempStorage) const
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

double BeerLambertFit::gradAndValue(const std::vector<double> &x, std::vector<double> &gradDest) const
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

double BeerLambertFit::getScale(int f, const std::vector<double> &x)
{
	const double a      = x[0];
	const double kappa  = x[1];
	const double phi    = x[2];

	d3Vector n = cos(phi) * tilt_q + sin(phi) * tilt_p;
	n.normalize();

	const double cos_f = n.dot(view_dir[f]);

	return a * pow(kappa, 1.0 / std::abs(cos_f));
}


MultiBeerLambertFit::MultiBeerLambertFit(
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

double MultiBeerLambertFit::gradAndValue(const std::vector<double> &x, std::vector<double> &gradDest) const
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

double MultiBeerLambertFit::getScale(int t, int f, const std::vector<double> &x)
{
	const double a = x[0] * fractional_dose[t];

	const double kappa  = x[2*t + 1];
	const double phi    = x[2*t + 2];

	d3Vector n = cos(phi) * tilt_q[t] + sin(phi) * tilt_p[t];
	n.normalize();

	const double cos_f = n.dot(view_dir[t][f]);

	return a * pow(kappa, 1.0 / std::abs(cos_f));
}

void MultiBeerLambertFit::report(int iteration, double cost, const std::vector<double> &x) const
{
	std::cout.precision(16);

	std::cout << iteration << ": " << cost << std::endl;
}

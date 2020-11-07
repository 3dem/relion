#include "ctf_refinement.h"
#include "aberration_fit.h"
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
#include <src/jaz/optics/damage.h>
#include <src/jaz/math/Zernike_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
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
		
		int def_section = parser.addSection("Defocus refinement options");

		do_refine_defocus = !parser.checkOption("--no_defocus", "Do not refine the (astigmatic) defocus.");
		do_refine_defocus = !parser.checkOption("--no_aberrations", "Do not refine higher-order aberrations");
		do_refine_defocus = !parser.checkOption("--no_scale", "Do not refine the contrast scale");

		lambda_reg = textToDouble(parser.getOption("--lambda", "Regularisation scale", "0.1"));
		
		minDelta = textToDouble(parser.getOption("--d0", "Min. defocus offset to test [Å]", "-3000"));
		maxDelta = textToDouble(parser.getOption("--d1", "Max. defocus offset to test [Å]", "3000"));
		deltaSteps = textToInteger(parser.getOption("--ds", "Number of defocus steps in-between", "100"));
		
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
		
	const int tc = particles.size();
	const bool flip_value = true;

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize);
	
	Log::endSection();


	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;
		
		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));
		Log::print("Loading");
		
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);

		const int fc = tomogram.frameCount;
		
		particleSet.checkTrajectoryLengths(
				particles[t][0], pc, fc, "CtfRefinementProgram::run");

		BufferedImage<float> freqWeights = computeFrequencyWeights(
			tomogram, true, 0.0, 0.0, true, num_threads);

		BufferedImage<float> doseWeight = tomogram.computeDoseWeight(boxSize,1);


		const int s = referenceMap.image_real[0].xdim;
		const int sh = s/2 + 1;



		const EvenData evenZero({0.0, 0.0, 0.0, 0.0, 0.0});
		const OddData oddZero({0.0, dComplex(0.0, 0.0)});

		if (do_refine_defocus)
		{

			BufferedImage<EvenData> evenData(sh,s,fc);
			evenData.fill(evenZero);

			BufferedImage<OddData> oddData(sh,s,fc);
			oddData.fill(oddZero);

			std::vector<int> chronoOrder = IndexSort<double>::sortIndices(tomogram.cumulativeDose);

			for (int f = 1; f < fc; f++)
			{
				tomogram.centralCTFs[chronoOrder[f]] = tomogram.centralCTFs[chronoOrder[0]];
			}



			Log::beginProgress("Accumulating evidence", pc * fc / num_threads);

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
						aberrationsCache, flip_value, freqWeights, doseWeight,
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

			/*{
				EvenData::write(evenData, outDir+"DEBUG");
			}*/

			/*{
				evenData = EvenData::read(outDir+"DEBUG");
			}*/

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
			}

			if (diag)
			{
				std::ofstream meanDefocus(outDir + tomogram.name + "_mean_defocus.dat");

				for (int f = 0; f < fc; f++)
				{
					meanDefocus << f << ' ' << (astig[f][0] + astig[f][1]) / 2.0 << '\n';
				}
			}

		}

		Log::endSection();
		
	} // all tomograms


	tomogramSet.write(outDir+"tomograms.star");

	optimisationSet.tomograms = outDir+"tomograms.star";
	optimisationSet.write(outDir+"optimisation_set.star");
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

	/*std::vector<double> nmOpt = NelderMead::optimize(
		initial, problem, initialStep, 0.000001, 2000, 1.0, 2.0, 0.5, 0.5, false);*/

	std::vector<double> nmOpt = LBFGS::optimize(initial, problem, false, 300, 1e-7, 1e-6);

	std::vector<d3Vector> out(fc);

	//std::cout.precision(16);

	for (int f = 0; f < fc; f++)
	{
		/*for (int i = 0; i < 3; i++)
		{
			std::cout << nmOpt[3*f + i + 3] << "  ";
		}

		std::cout << "\n";*/

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

GlobalDefocusFit::GlobalDefocusFit(
		const BufferedImage<double>& dataTerm,
		double defocusStep,
		double lambda)
:
	dataTerm(dataTerm),
	defocusStep(defocusStep),
	lambda(lambda)
{
}

double GlobalDefocusFit::f(const std::vector<double>& x, void *tempStorage) const
{
	const double meanDefocus = x[0];
	const int fc = dataTerm.xdim;

	double sum = 0;

	for (int ff = 0; ff < fc; ff++)
	{
		const double dataCost = Interpolation::linearXY_clip(dataTerm, ff, x[ff+1]);
		const double delta = defocusStep * (x[ff+1] - meanDefocus);


	}
}

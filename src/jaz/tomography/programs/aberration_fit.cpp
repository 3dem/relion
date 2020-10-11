#include "aberration_fit.h"
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h> 
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/math/Zernike_helper.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/single_particle/ctf/magnification_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <omp.h> 

using namespace gravis;


AberrationFitProgram::AberrationFitProgram(int argc, char *argv[])
:	RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	readParams(parser);
}

void AberrationFitProgram::readParams(IOParser &parser)
{
	try
	{
		_readParams(parser);
				
		int defocus_section = parser.addSection("Alignment options");
		
		bool perTomogram = parser.checkOption("--per_tomo", "Estimate the aberrations per tomogram instead of globally");
		
		granularity = perTomogram? PerTomogram : Global;
		n_even = textToInteger(parser.getOption("--ne", "Maximal even N", "4"));
		n_odd = textToInteger(parser.getOption("--no", "Maximal odd N", "3"));
		
		Log::readParams(parser);
		
		if (parser.checkForErrors()) 
		{
			parser.writeUsage(std::cout);
			std::exit(-1);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
}

void AberrationFitProgram::run()
{
	Log::beginSection("Initialising");
	
	RefinementProgram::init();
		
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int tc = particles.size();
	const bool flip_value = true;
	
	Log::printBinaryChoice("Estimating the aberrations ", granularity == PerTomogram, 
						   "per tomogram", "globally");
	
	Log::endSection();
		
	
	
	std::vector<BufferedImage<EvenData>> evenData_thread(num_threads);
	std::vector<BufferedImage<OddData>> oddData_thread(num_threads);
	
	EvenData evenZero({0.0, 0.0, 0.0, 0.0, 0.0});
	OddData oddZero({0.0, dComplex(0.0, 0.0)});
	
	BufferedImage<EvenData> evenData(sh,s);
	evenData.fill(evenZero);
	
	BufferedImage<OddData> oddData(sh,s);
	oddData.fill(oddZero);
	
	// split tomograms by pixel size!
	
	double lastPixelSize = 0.0;
	
	
	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;
		
		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));		
		Log::print("Loading");
		
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		lastPixelSize = tomogram.optics.pixelSize;
		
		std::string diagPrefix = outDir + "diag_" + ZIO::itoa(t);
		
		
		BufferedImage<float> frqWeight = computeFrequencyWeights(
			tomogram, true, 1.0, 0.0, num_threads);
		
		if (diag)
		{
			frqWeight.write(diagPrefix + "_noise_weight.mrc");
			referenceMap.freqWeight.write(diagPrefix + "_FSC_weight.mrc");
		}
				
		for (int th = 0; th < num_threads; th++)
		{
			evenData_thread[th] = BufferedImage<EvenData>(sh,s);
			evenData_thread[th].fill(evenZero);
			
			oddData_thread[th] = BufferedImage<OddData>(sh,s);
			oddData_thread[th].fill(oddZero);
		}
				
		Log::beginProgress("Accumulating evidence", pc/num_threads);
		
		#pragma omp parallel for num_threads(num_threads)		
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();
			
			if (th == 0)
			{
				Log::updateProgress(p);
			}
			
			considerParticle(
				particles[t][p], tomogram, referenceMap, dataSet, flip_value, frqWeight,
				0, -1,
				evenData_thread[th], oddData_thread[th]);
		}
		
		Log::endProgress();
		
		for (int th = 0; th < num_threads; th++)
		{
			evenData += evenData_thread[th];
			oddData += oddData_thread[th];
		}
				
		if (granularity == PerTomogram)
		{
			std::vector<double> initialEven(Zernike::numberOfEvenCoeffs(n_even), 0.0);
			std::vector<double> initialOdd(Zernike::numberOfEvenCoeffs(n_odd), 0.0);

			solveAndFitEven(
				evenData, n_even, initialEven,
				tomogram.optics.pixelSize, outDir + ZIO::itoa(t) + "_", true);

			solveAndFitOdd(
				oddData, n_odd, initialOdd,
				tomogram.optics.pixelSize, outDir + ZIO::itoa(t) + "_", true);
			
			evenData.fill(evenZero);
			oddData.fill(oddZero);
		}
		
		Log::endSection();
	}
	
	if (granularity == Global)
	{
		std::vector<double> initialEven(Zernike::numberOfEvenCoeffs(n_even), 0.0);
		std::vector<double> initialOdd(Zernike::numberOfEvenCoeffs(n_odd), 0.0);

		solveAndFitEven(
			evenData, n_even, initialEven,
			lastPixelSize, outDir, true);

		solveAndFitOdd(
			oddData, n_odd, initialOdd,
			lastPixelSize, outDir, true);
	}
}

void AberrationFitProgram :: considerParticle(
		int part_id,
		const Tomogram& tomogram, 
		const TomoReferenceMap& referenceMap, 
		const ParticleSet& dataSet,
		bool flip_value,
		const BufferedImage<float>& frqWeight,
		int f0, int f1,
		BufferedImage<EvenData>& even_out, 
		BufferedImage<OddData>& odd_out)
{
	const int s = referenceMap.image_real[0].xdim;
	const int sh = s/2 + 1;
	const int fc = tomogram.frameCount;
	const double pix2ang = 1.0 / ((double)s * tomogram.optics.pixelSize);
	
	if (f1 < 0) f1 = fc - 1;
	
	
	const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(
				part_id, fc, tomogram.optics.pixelSize);
	
	d4Matrix projCut;
	
	
	BufferedImage<fComplex> observation(sh,s);
	
	for (int f = f0; f <= f1; f++)
	{
		TomoExtraction::extractFrameAt3D_Fourier(
				tomogram.stack, f, s, 1.0, tomogram.projectionMatrices[f], traj[f],
				observation, projCut, 1, false, true);
		
		CTF ctf = tomogram.getCtf(f, dataSet.getPosition(part_id));

		BufferedImage<fComplex> prediction = Prediction::predictModulated(
				part_id, dataSet, projCut, s, 
				ctf,
				tomogram.optics.pixelSize,
				referenceMap.image_FS, 
				Prediction::OwnHalf,
				Prediction::Unmodulated);
		
		const float scale = flip_value? -1.f : 1.f;
		
		observation(0,0) = fComplex(0.f, 0.f);
		prediction(0,0) = fComplex(0.f, 0.f);

		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double x_ang = pix2ang * x;
			const double y_ang = pix2ang * (y < s/2? y : y - s);

			const double gamma = ctf.getLowOrderGamma(x_ang, y_ang);
			const double cg = cos(gamma);
			const double sg = sin(gamma);
			const double c = -sg;

			fComplex zobs = observation(x,y);
			fComplex zprd = scale * prediction(x,y);
			
			const double zz = zobs.real * zprd.real + zobs.imag * zprd.imag;
			const double zq = zobs.imag * zprd.real - zobs.real * zprd.imag;
			const double nr = zprd.norm();
			const double wg = frqWeight(x,y,f);
			
			
			// NOTE: the prediction contains neither phase nor amp modulation!
			// @TODO: when phase shifts are supported, consider them here.
			
			
			EvenData& ed = even_out(x,y);
			
			ed.Axx += wg * nr * sg * sg;
			ed.Axy += wg * nr * cg * sg;
			ed.Ayy += wg * nr * cg * cg;
	
			ed.bx -= wg * zz * sg;
			ed.by -= wg * zz * cg;
			
			
			OddData& od = odd_out(x,y);
			
			od.a += wg * c * c * nr;
			
			od.b.real += wg * c * zz;
			od.b.imag += wg * c * zq;
		}
	}
}

AberrationFitProgram::EvenSolution AberrationFitProgram::solveEven(
		const BufferedImage<EvenData>& data)
{
	const double eps = 1e-30;
	const int s  = data.ydim;
	const int sh = data.xdim;

	EvenSolution out;

	out.optimum = BufferedImage<dComplex>(sh,s);
	out.phaseShift = BufferedImage<double>(sh,s);
	out.weight = BufferedImage<Tensor2x2<double>>(sh,s);

	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		EvenData d = data(x,y);

		d2Vector b(d.bx, d.by);
		d2Matrix A(d.Axx, d.Axy, d.Axy, d.Ayy);

		const double det = A(0,0) * A(1,1) - A(0,1) * A(1,0);
		const double absDet = std::abs(det);

		if (absDet > eps)
		{
			d2Matrix Ai = A;
			Ai.invert();

			const d2Vector opt = Ai * b;

			out.optimum(x,y) = dComplex(opt.x, opt.y);
			out.phaseShift(x,y) = std::abs(opt.x) > 0.0? atan2(opt.y, opt.x) : 0.0;
			out.weight(x,y) = Tensor2x2<double>(d.Axx, d.Axy, d.Ayy);
		}
		else
		{
			out.optimum(x,y) = dComplex(0.0, 0.0);
			out.phaseShift(x,y) = 0.0;
			out.weight(x,y) = Tensor2x2<double>(0.0, 0.0, 0.0);
		}
	}

	return out;
}

std::vector<double> AberrationFitProgram::fitEven(
		const EvenSolution& solution,
		int n_bands,
		const std::vector<double>& initialCoeffs,
		double pixelSize,
		const std::string& prefix,
		bool writeImages)
{
	const d2Matrix mag(1.0, 0.0, 0.0, 1.0);

	const int cc = Zernike::numberOfEvenCoeffs(n_bands);

	if (initialCoeffs.size() != cc)
	{
		REPORT_ERROR_STR(
			"AberrationFitProgram::solveEven: " << initialCoeffs.size() <<
			" initial coefficient provided, but " << cc << " are required.");
	}
	
	/*if (writeImages)
	{
		Centering::fftwHalfToHumanFull(phaseShift).write(prefix + "even_phase_per-pixel.mrc");
	}*/

	BufferedImage<double> nonlinearFit;

	/*if (writeImages)
	{
		Centering::fftwHalfToHumanFull(nonlinearFit).write(prefix + "even_phase_nonlinear-fit.mrc");
	}*/

	std::vector<double> coeffs = ZernikeHelper::optimiseEvenZernike(
		solution.optimum,
		solution.weight,
		pixelSize, 
		mag,
		n_bands-1,
		initialCoeffs,
		&nonlinearFit);
	
	/*if (writeImages)
	{
		Centering::fftwHalfToHumanFull(nonlinearFit).write(prefix + "even_phase_nonlinear-fit.mrc");
	}*/
	
	// @TODO: write coefficients to a file
	
	return coeffs;
}

gravis::d3Vector AberrationFitProgram::findAstigmatism(
		const AberrationFitProgram::EvenSolution& solution,
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

BufferedImage<double> AberrationFitProgram::plotAstigmatism(
		const AberrationFitProgram::EvenSolution& solution,
		const CTF& referenceCtf,
		double initialDeltaZ,
		double range,
		double pixelSize,
		int size)
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
		astigBasis(xi,yi,2) = xx * yy;
	}

	ZernikeHelper::AnisoBasisOptimisation problem(
				solution.optimum, solution.weight, astigBasis, false);

	std::vector<double> globalOpt(3,0.0);
	double minCost = std::numeric_limits<double>::max();

	void* tempStorage = problem.allocateTempStorage();

	const int steps = size-1;
	const double mid = steps/2;
	const double scale = 2.0 * K1 * range;

	BufferedImage<double> out(size, size, 1);

	for (int a1i = 0; a1i < size; a1i++)
	for (int a2i = 0; a2i < size; a2i++)
	{
		const double dz = -initialDeltaZ * K1;
		const double a1 = scale * (a1i - mid);
		const double a2 = scale * (a2i - mid);

		std::vector<double> params = {dz, a1, a2};

		const double cost = problem.f(params, tempStorage);

		if (cost < minCost)
		{
			minCost = cost;
			globalOpt = params;
		}

		out(a1i,a2i) = cost;
	}

	return out;
}

std::vector<double> AberrationFitProgram::solveAndFitEven(
		const BufferedImage<AberrationFitProgram::EvenData> &data,
		int n_bands,
		const std::vector<double> &initialCoeffs,
		double pixelSize,
		const std::string& prefix,
		bool writeImages)
{
	EvenSolution solution = solveEven(data);
	return fitEven(solution, n_bands, initialCoeffs, pixelSize, prefix, writeImages);
}

double AberrationFitProgram::findDefocus(
		const BufferedImage<EvenData>& evenData,
		double pixelSize,
		const CTF& ctf0,
		double minDefocus,
		double maxDefocus,
		int steps)
{
	const int s = evenData.ydim;
	const int sh = evenData.xdim;
	const double as = s * pixelSize;
	const double eps = 1e-30;
	const double deltaStep = (maxDefocus - minDefocus) / (double) (steps - 1);


	CTF ctfz = ctf0;

	double best_deltaZ = minDefocus;
	double minCost = std::numeric_limits<double>::max();


	for (int di = 0; di < steps; di++)
	{
		const double deltaZ = minDefocus + di * deltaStep;

		ctfz.DeltafU = ctf0.DeltafU + deltaZ;
		ctfz.DeltafV = ctf0.DeltafV + deltaZ;

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
				AberrationFitProgram::EvenData d = evenData(x,y);

				d2Vector b(d.bx, d.by);
				d2Matrix A(d.Axx, d.Axy, d.Axy, d.Ayy);

				const double det = A(0,0) * A(1,1) - A(0,1) * A(1,0);

				if (std::abs(det) > eps)
				{
					d2Matrix Ai = A;
					Ai.invert();

					const d2Vector opt = Ai * b;

					const double gamma_0 = ctf0.getLowOrderGamma(xx/as, yy/as);
					const double gamma_z = ctfz.getLowOrderGamma(xx/as, yy/as);
					const double delta = gamma_z - gamma_0;

					const d2Vector dx = d2Vector(cos(delta), sin(delta)) - opt;

					cost += dx.dot(A * dx);
				}
			}
		}

		if (cost < minCost)
		{
			minCost = cost;
			best_deltaZ = deltaZ;
		}
	}

	return best_deltaZ;
}

AberrationFitProgram::OddSolution AberrationFitProgram::solveOdd(
		const BufferedImage<OddData>& data)
{
	const int s  = data.ydim;
	const int sh = data.xdim;

	OddSolution out;

	out.optimum = BufferedImage<dComplex>(sh,s);
	out.phaseShift = BufferedImage<double> (sh,s);
	out.weight = BufferedImage<double>(sh,s);

	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		OddData d = data(x,y);

		if (d.a > 0.0)
		{
			out.optimum(x,y) = d.b / d.a;
			out.phaseShift(x,y) = d.b.arg();
			out.weight(x,y) = d.a;
		}
		else
		{
			out.optimum(x,y) = dComplex(0.0, 0.0);
			out.phaseShift(x,y) = 0.0;
			out.weight(x,y) = 0.0;
		}
	}

	return out;
}

std::vector<double> AberrationFitProgram::fitOdd(
		const OddSolution& solution,
		int n_bands,
		const std::vector<double>& initialCoeffs,
		double pixelSize,
		const std::string& prefix,
		bool writeImages)
{
	const d2Matrix mag(1.0, 0.0, 0.0, 1.0);

	const int cc = Zernike::numberOfOddCoeffs(n_bands - 1);

	if (initialCoeffs.size() != cc)
	{
		REPORT_ERROR_STR(
			"AberrationFitProgram::solveOdd: " << initialCoeffs.size() <<
			" initial coefficient provided, but " << cc << " are required.");
	}

	/*if (writeImages)
	{
		Centering::fftwHalfAntisymmetricalToHumanFull(phaseShift).write(
					prefix + "odd_phase_per-pixel.mrc");
	}*/
	
	BufferedImage<double> nonlinearFit;
	
	std::vector<double> coeffs = ZernikeHelper::optimiseOddZernike(
		solution.optimum,
		solution.weight,
		pixelSize, 
		mag,
		n_bands-1,
		initialCoeffs,
		&nonlinearFit);
	
	/*if (writeImages)
	{
		Centering::fftwHalfAntisymmetricalToHumanFull(nonlinearFit).write(
					prefix + "odd_phase_nonlinear-fit.mrc");
	}*/
	
	return coeffs;
}

std::vector<double> AberrationFitProgram::solveAndFitOdd(
		const BufferedImage<AberrationFitProgram::OddData> &data,
		int n_bands,
		const std::vector<double> &initialCoeffs,
		double pixelSize,
		const std::string& prefix,
		bool writeImages)
{

	OddSolution solution = solveOdd(data);
	return fitOdd(solution, n_bands, initialCoeffs, pixelSize, prefix, writeImages);
}

AberrationFitProgram::EvenData &AberrationFitProgram::EvenData::operator+=(const EvenData& d)
{
	Axx += d.Axx;
	Axy += d.Axy;
	Ayy += d.Ayy;
	bx += d.bx;
	by += d.by;
}

AberrationFitProgram::OddData &AberrationFitProgram::OddData::operator+=(const OddData& d)
{
	a += d.a;
	b += d.b;
}

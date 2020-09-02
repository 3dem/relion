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
			solveEven(evenData, n_even, tomogram.optics.pixelSize, outDir + ZIO::itoa(t) + "_", true);
			solveOdd(oddData, n_odd, tomogram.optics.pixelSize, outDir + ZIO::itoa(t) + "_", true);
			
			evenData.fill(evenZero);
			oddData.fill(oddZero);
		}
		
		Log::endSection();
	}
	
	if (granularity == Global)
	{
		solveEven(evenData, n_even, lastPixelSize, outDir, true);
		solveOdd(oddData, n_odd, lastPixelSize, outDir, true);
	}
}

void AberrationFitProgram :: considerParticle(
		int part_id,
		const Tomogram& tomogram, 
		const TomoReferenceMap& referenceMap, 
		const ParticleSet* dataSet,
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
	
	
	const std::vector<d3Vector> traj = dataSet->getTrajectoryInPix(
				part_id, fc, tomogram.optics.pixelSize);
	
	d4Matrix projCut;					
	
	
	BufferedImage<fComplex> observation(sh,s);
	
	for (int f = f0; f <= f1; f++)
	{
		TomoExtraction::extractFrameAt3D_Fourier(
				tomogram.stack, f, s, 1.0, tomogram.proj[f], traj[f],
				observation, projCut, 1, false, true);
					
		BufferedImage<fComplex> prediction = Prediction::predictFS(
				part_id, dataSet, projCut, s, 
				tomogram.centralCTFs[f], tomogram.centre,
				tomogram.handedness, tomogram.optics.pixelSize,
				referenceMap.image_FS, 
				Prediction::OppositeHalf,
				Prediction::Unmodulated);
		
		CTF ctf = TomoCtfHelper::adaptToParticle(
				tomogram.centralCTFs[f], tomogram.proj[f], traj[f], tomogram.centre, 
				tomogram.handedness, tomogram.optics.pixelSize);
		
		const float scale = flip_value? -1.f : 1.f;
		
		observation(0,0) = fComplex(0.f, 0.f);
		prediction(0,0) = fComplex(0.f, 0.f);
		
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			const double x_ang = pix2ang * x;
			const double y_ang = pix2ang * (y < sh? y : y - s);
			
	
			const double gamma = ctf.getGamma(x_ang, y_ang);
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

std::vector<double> AberrationFitProgram::solveEven(
		const BufferedImage<EvenData>& data,
		int n_bands,
		double pixelSize, 
		std::string prefix,
		bool writeImages)
{
	const double eps = 1e-30;
	const int s  = data.ydim;
	const int sh = data.xdim;
	const d2Matrix mag(1.0, 0.0, 0.0, 1.0);
	
	BufferedImage<dComplex> optimum(sh,s);
	BufferedImage<double> phaseShift(sh,s);
	BufferedImage<double> initialWeight(sh,s);
	BufferedImage<Tensor2x2<double>> weight(sh,s);
	
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
			
			optimum(x,y) = dComplex(opt.x, opt.y);
			phaseShift(x,y) = std::abs(opt.x) > 0.0? atan2(opt.y, opt.x) : 0.0;
			initialWeight(x,y) = sqrt(absDet);
			weight(x,y) = Tensor2x2<double>(d.Axx, d.Axy, d.Ayy);
		}
		else
		{
			optimum(x,y) = dComplex(0.0, 0.0);
			phaseShift(x,y) = 0.0;
			initialWeight(x,y) = 0.0;
			weight(x,y) = Tensor2x2<double>(0.0, 0.0, 0.0);
		}
	}
	
	if (writeImages)
	{
		Centering::fftwHalfToHumanFull(phaseShift).write(prefix + "even_phase_per-pixel.mrc");
	}
		
	BufferedImage<double> linearFit, nonlinearFit;
			
	std::vector<double> coeffs0 = ZernikeHelper::fitEvenZernike(
		phaseShift, 
		initialWeight, 
		pixelSize, 
		mag, 
		n_bands, 
		&linearFit);
	
	if (writeImages)
	{
		Centering::fftwHalfToHumanFull(linearFit).write(prefix + "even_phase_linear-fit.mrc");
	}
	
	std::vector<double> coeffs = ZernikeHelper::optimiseEvenZernike(
		optimum, 
		weight, 
		pixelSize, 
		mag,
		n_bands, 
		coeffs0, 
		&nonlinearFit);
	
	if (writeImages)
	{
		Centering::fftwHalfToHumanFull(nonlinearFit).write(prefix + "even_phase_nonlinear-fit.mrc");
	}
	
	// @TODO: write coefficients to a file
	
	return coeffs;
}

std::vector<double> AberrationFitProgram::solveOdd(
		const BufferedImage<OddData>& data,
		int n_bands,
		double pixelSize, 
		std::string prefix,
		bool writeImages)
{
	const int s  = data.ydim;
	const int sh = data.xdim;
	const d2Matrix mag(1.0, 0.0, 0.0, 1.0);
	
	BufferedImage<dComplex> optimum(sh,s);
	BufferedImage<double> phaseShift(sh,s);
	BufferedImage<double> weight(sh,s);
	
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		OddData d = data(x,y);
		
		if (d.a > 0.0)
		{
			optimum(x,y) = d.b / d.a;
			phaseShift(x,y) = d.b.arg();
			weight(x,y) = d.a;
		}
		else
		{
			optimum(x,y) = dComplex(0.0, 0.0);
			phaseShift(x,y) = 0.0;
			weight(x,y) = 0.0;
		}
	}
	
	if (writeImages)
	{
		Centering::fftwHalfAntisymmetricalToHumanFull(phaseShift).write(
					prefix + "odd_phase_per-pixel.mrc");
	}
	
	BufferedImage<double> linearFit, nonlinearFit;	
			
	std::vector<double> coeffs0 = ZernikeHelper::fitOddZernike(
		optimum, 
		weight, 
		pixelSize, 
		mag, 
		n_bands, 
		&linearFit);
	
	if (writeImages)
	{
		Centering::fftwHalfAntisymmetricalToHumanFull(linearFit).write(
					prefix + "odd_phase_linear-fit.mrc");
	}
	
	std::vector<double> coeffs = ZernikeHelper::optimiseOddZernike(
		optimum, 
		weight, 
		pixelSize, 
		mag,
		n_bands, 
		coeffs0, 
		&nonlinearFit);
	
	if (writeImages)
	{
		Centering::fftwHalfAntisymmetricalToHumanFull(nonlinearFit).write(
					prefix + "odd_phase_nonlinear-fit.mrc");
	}
	
	return coeffs;
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

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
	
	RefinementProgram::init();
		
	const int tc = particles.size();
	const bool flip_value = true;
	
	Log::beginSection("Configuration");
	Log::printBinaryChoice("Frame angles: ", mfSettings.constAngles, "static", "variable");
	Log::printBinaryChoice("Frame shifts: ", mfSettings.constShifts, "static", "variable");
	Log::printBinaryChoice("Particle positions: ", mfSettings.constParticles, "static", "variable");
	Log::endSection();
	
	
	int tpc = dataSet.getTotalParticleNumber();
		
	if (dataSet.motionTrajectories.size() != tpc)
	{
		dataSet.motionTrajectories.resize(tpc);
	}
			
	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;
		
		Tomogram tomogram = tomogramSet.loadTomogram(t, false);
		const int fc = tomogram.frameCount;
		
		for (int p = 0; p < pc; p++)
		{
			if (dataSet.motionTrajectories[particles[t][p]].shifts_Ang.size() != fc)
			{		
				dataSet.motionTrajectories[particles[t][p]] = Trajectory(fc);
			}
		}
	}
	
	dataSet.hasMotion = true;
	
	Log::endSection();
		
	
	TomogramSet tomogramSetOut = tomogramSet;
	
	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;
		
		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));		
		Log::print("Loading");
		
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		
		const int fc = tomogram.frameCount;
		std::vector<d4Matrix> projTomoCorr = tomogram.projectionMatrices;
		
		
		std::string tag = ZIO::itoa(t);
		std::string diagPrefix = outDir + "diag_" + tag;
		
		
		BufferedImage<float> frqWeight = computeFrequencyWeights(
			tomogram, whiten, sig2RampPower, hiPass_px, num_threads);
		
		if (diag)
		{
			frqWeight.write(diagPrefix + "_frq_weight.mrc");
		}
		
		
		std::vector<d4Matrix> projByTime(fc);
		
		for (int f = 0; f < fc; f++)
		{
			projByTime[f] = tomogram.projectionMatrices[tomogram.frameSequence[f]];
		}
		
		const double t0 = omp_get_wtime();
		
		
		std::vector<BufferedImage<double>> CCs = Prediction::computeCroppedCCs(
				dataSet, particles[t], tomogram, referenceMap, frqWeight, tomogram.frameSequence,
				range, flip_value, num_threads, padding);
					
		
		MotionFit motionFit(
				CCs, projByTime, dataSet, particles[t], referenceMap.image_FS, 
				motParams, mfSettings, tomogram.centre, 
				tomogram.getFrameDose(), tomogram.optics.pixelSize, padding, num_threads);
		
		
		
		
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
				dataSet, particles[t], tomogram, referenceMap.image_FS,
				flip_value, num_threads);
			
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
		
		Log::beginProgress("Performing optimisation", num_iters);
				
		const double t1 = omp_get_wtime();
		
		std::vector<double> opt = LBFGS::optimize(
			initial, motionFit, 1, num_iters, 1e-3, 1e-4);
		
		const double t2 = omp_get_wtime();
		
		Log::endProgress();
						
		if (timing)
		{
			Log::beginSection("Time");
			Log::print("setup:        " + ZIO::itoa(t1-t0) + " sec");
			Log::print("optimisation: " + ZIO::itoa(t2-t1) + " sec");
			Log::endSection();
		}
		
		projTomoCorr = motionFit.getProjections(opt, tomogram.frameSequence);
		motionFit.shiftParticles(opt, dataSet);
		motionFit.exportTrajectories(opt, dataSet, tomogram.frameSequence);
		
		tomogramSetOut.setProjections(t, projTomoCorr);
				
		
		Mesh mesh8 = motionFit.visualiseTrajectories(opt, 8.0);
		mesh8.writePly(outDir + "tracks_" + tag + "_x8.ply");
		
		Mesh mesh1 = motionFit.visualiseTrajectories(opt, 1.0);
		mesh1.writePly(outDir + "tracks_" + tag + "_x1.ply");
		
		if (diag)
		{
			tomogram.projectionMatrices = projTomoCorr;
			
			FCC3 = FCC::compute3(
				dataSet, particles[t], tomogram, referenceMap.image_FS,
				flip_value, num_threads);
			
			FCC3.write(diagPrefix + "_FCC3_final.mrc");
			BufferedImage<double> FCC1_new = FCC::divide(FCC3);
			FCC1_new.write(diagPrefix + "_FCC_final.mrc");
			
			(FCC1_new - FCC1).write(diagPrefix + "_FCC_delta.mrc");
			
			/*if (outputShiftedCCs)
			{
				Log::print("Computing shifted cross-correlations");
				
				std::vector<Image<double>> shiftedCCs = motionFit.drawShiftedCCs(opt);
				
				const int diam = shiftedCCs[0].xdim;
				
				Image<double> CCsum(diam, diam, fc);
						
				CCsum.fill(0.f);
				
				for (int p = 0; p < pc; p++)
				{
					CCsum += shiftedCCs[p];
				}				
				
				CCsum.write(diagPrefix + "_CC_sum_" + tag + "_final.mrc");
				
				Image<double> specCC0 = specCC;
				
				const int d = CCsum.xdim;
				const int dh = d/2 + 1;
				
				specCC = Image<double>(dh,fc);
				specCC.fill(0.0);
					
				Image<fComplex> CCsumFS;
				
				for (int f = 0; f < fc; f++)
				{
					Image<float> CCsum_f = CCsum.getSliceRef(f);
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
				
				specCC.write(diagPrefix + "_specCC_final.mrc");				
				(specCC - specCC0).write(diagPrefix + "_specCC_delta.mrc");
			}*/
		}
		
		Log::endSection();
	}
	
	Trajectory::write(dataSet.motionTrajectories, outDir + "motion.star");
	tomogramSetOut.write(outDir + "tomograms.star");
	dataSet.write(outDir + "particles.star");
}

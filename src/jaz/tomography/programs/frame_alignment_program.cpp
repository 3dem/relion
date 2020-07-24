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
	RefinementProgram::init();
		
	const int tc = particles.size();
	const bool flip_value = true;
	
	Log::beginSection("Configuration");
	Log::printBinaryChoice("Frame angles: ", const_angles, "static", "variable");
	Log::printBinaryChoice("Frame shifts: ", const_shifts, "static", "variable");
	Log::printBinaryChoice("Particle positions: ", const_particles, "static", "variable");
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
		
		std::vector<int> dummySeq(fc);
		
		for (int f = 0; f < fc; f++)
		{
			dummySeq[f] = f;
		}
		
		
		std::vector<BufferedImage<double>> CCs = Prediction::computeCroppedCCs(
				dataSet, particles[t], tomogram, referenceMap, frqWeight, dummySeq,
				range, flip_value, num_threads, padding);
		
		
		if (first_frame > 0 || last_frame >= 0)
		{
			#pragma omp parallel for num_threads(num_threads)		
			for (int p = 0; p < pc; p++)
			{
				for (int ff = 0; ff < fc; ff++)
				{
					int f = tomogram.frameSequence[ff];
					
					if (ff < first_frame || ff > last_frame)
					{
						CCs[p].getSliceRef(f) *= 0.f;
					}
				}
			}
		}
		
		ProtoAlignment protoAlignment(
				CCs, tomogram.projectionMatrices, dataSet, particles[t], referenceMap.image_FS, 
				const_particles, const_angles, const_shifts, range,
				tomogram.centre, num_threads, padding);
		
		BufferedImage<double> FCC3, FCC1;
		
		if (diag)
		{
			FCC3 = FCC::compute3(
				dataSet, particles[t], tomogram, referenceMap.image_FS,
				flip_value, num_threads);
			
			FCC3.write(diagPrefix + "_FCC3_initial.mrc");
			FCC1 = FCC::divide(FCC3);
			FCC1.write(diagPrefix + "_FCC_initial.mrc");
		}
		
		if (shiftOnly)
		{
			const int diam = (int)(2*range*padding);
			
			BufferedImage<float> CCsum(diam, diam, fc);
			
			CCsum.fill(0.f);
			
			Log::beginProgress("Adding up cross-correlations", pc);
			
			for (int p = 0; p < pc; p++)
			{
				Log::updateProgress(p);
				CCsum += protoAlignment.CCs[p];
			}	
			
			Log::endProgress();
						
			CCsum.write(outDir + "CCsum_" + tag + ".mrc");			
			
			d2Vector origin(padding*range, padding*range);
			
			std::ofstream frameShifts(outDir + "frame_shifts_" + tag + ".txt");
			
			for (int f = 0; f < fc; f++)
			{
				d2Vector opt = (Interpolation::quadraticMaxXY(CCsum.getSliceRef(f)) - origin)/padding;
				
				frameShifts << f << " " << opt.x << " " << opt.y << std::endl;
				
				projTomoCorr[f](0,3) += opt.x;
				projTomoCorr[f](1,3) += opt.y;
			}
		}
		else
		{
			std::vector<double> initial(protoAlignment.getParamCount(), 0.0);
			
			Log::beginProgress("Performing optimisation", num_iters);
			
			std::vector<double> opt = LBFGS::optimize(
				initial, protoAlignment, 1, num_iters, 1e-6, 1e-4);
			
			Log::endProgress();

			
			projTomoCorr = protoAlignment.getProjections(opt);
			protoAlignment.shiftParticles(opt, particles[t], dataSet);
			
		}
		
		tomogramSetOut.setProjections(t, projTomoCorr);
		
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
		}
		
		Log::endSection();
	}
	
	tomogramSetOut.write(outDir + "tomograms.star");
	
	if (!shiftOnly)
	{
		dataSet->write(outDir + "particles.star");
	}
}

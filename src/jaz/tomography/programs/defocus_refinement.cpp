#include "defocus_refinement.h"
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
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <iostream>
#include <src/jaz/math/Zernike_helper.h>
#include <src/time.h>

#define TIMING 0


using namespace gravis;

DefocusRefinementProgram::DefocusRefinementProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	
	try
	{
		_readParams(parser);
		
		int def_section = parser.addSection("Defocus refinement options");
				
		//scanDefocus = parser.checkOption("--scan_defocus", "Perform brute-force defocus scan");
		//refineFast = !parser.checkOption("--scan_only", "Perform only a brute-force scan");
		scanDefocus = true;
		refineFast = false;
		max_particles = textToInteger(parser.getOption("--max", "Max. number of particles to consider per tomogram", "-1"));
		group_count = textToInteger(parser.getOption("--g", "Number of independent groups", "10"));
		sigma_input = textToDouble(parser.getOption("--sig0", "Std. dev. of initial defoci (negative to turn off regularisation)", "-1"));
		regularise = sigma_input > 0.0;
		
		minDelta = textToDouble(parser.getOption("--d0", "Min. defocus offset to test [Å]", "-100"));
		maxDelta = textToDouble(parser.getOption("--d1", "Max. defocus offset to test [Å]", "100"));
		deltaSteps = textToInteger(parser.getOption("--ds", "Number of defocus steps in-between", "100"));
		
		clearAstigmatism = parser.checkOption("--ca", "Clear the current astigmatism estimate");
		
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

void DefocusRefinementProgram::run()
{
	Log::beginSection("Initialising");
	
	RefinementProgram::init();
		
	const int tc = particles.size();
	const bool flip_value = true;	
	
	Log::endSection();
		
	for (int t = 0; t < tc; t++)
	{
		int pc0 = particles[t].size();
		if (pc0 == 0) continue;
		
		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));		
		Log::print("Loading");
		
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		
		pc0 = group_count * (pc0 / group_count);		
		
		const int usedParticleCount = (max_particles > 0 && pc0 > max_particles)? max_particles : pc0;
		const int fc = tomogram.frameCount;
		
		dataSet->checkTrajectoryLengths(
				particles[t][0], usedParticleCount, fc, "DefocusRefinementProgram::run");
		
		
		if (clearAstigmatism)
		{
			for (int f = 0; f < fc; f++)
			{
				CTF& ctf = tomogram.centralCTFs[f];
				const double z0 = (ctf.DeltafU + ctf.DeltafV) / 2.0;
				
				ctf.DeltafU = z0;
				ctf.DeltafV = z0;
			}
		}
				
		const int f0 = first_frame;
		const int f1 = (last_frame > 0 && last_frame < fc)? last_frame : fc-1;
		
		std::vector<double> defocusOffset(fc, 0.0), offsetStdDev(fc, 0.0);
		
		/*
		  TODO: find best handedness from current defocus over all frames
		*/
		
		BufferedImage<float> freqWeights = computeFrequencyWeights(
			tomogram, true, 0.0, 0.0, num_threads);
		
					
		for (int f = f0; f <= f1; f++)
		{
			Log::beginSection("Frame " + ZIO::itoa(f+1));
			
			if (scanDefocus)
			{
				DefocusFit defocus = findDefocus(
					f, minDelta, maxDelta, deltaSteps,  
					group_count, sigma_input,	
					dataSet, particles[t], usedParticleCount, 
					tomogram, referenceMap.image_FS, 
					freqWeights, flip_value, num_threads);
				
				defocusOffset[f] = defocus.value;
				offsetStdDev[f] = defocus.stdDev;
				
				if (diag)
				{
					std::ofstream costOutByGroup(
						outDir+"cost_tomo_"+ZIO::itoa(t)+"_frame_"+ZIO::itoa(f)+"_by_group.dat");
					
					for (int g = 0; g < group_count; g++)
					{
						for (int di = 0; di < deltaSteps; di++)
						{
							costOutByGroup << defocus.offsets[di] << " " << defocus.costByGroup[g][di] << '\n';
						}
						
						costOutByGroup << '\n';
					}
						
					std::ofstream costOut(
						outDir+"cost_tomo_"+ZIO::itoa(t)+"_frame_"+ZIO::itoa(f)+".dat");
					
					for (int di = 0; di < deltaSteps; di++)
					{
						costOut << defocus.offsets[di] << " " << defocus.totalCost[di] << '\n';
					}
					
					costOut << '\n';
				}
	
				
				CTF& ctf = tomogram.centralCTFs[f];
				
				const double s02 = sigma_input * sigma_input;
				const double sf2 = offsetStdDev[f] * offsetStdDev[f];
				
				const double deltaZ = regularise?
					s02 * defocusOffset[f] / (sf2 + s02) :
					defocusOffset[f];

				ctf.DeltafU += deltaZ;
				ctf.DeltafV += deltaZ;
				
				ctf.initialise();
				
				tomogramSet.setCtf(t,f,ctf);
			}

			// @TODO: fix
			if (refineFast)
			{
				const int s = referenceMap.image_real[0].xdim;
				const int sh = s/2 + 1;
				const int pc = particles[t].size();
				
				std::vector<BufferedImage<AberrationFitProgram::EvenData>> evenData_thread(num_threads);
				std::vector<BufferedImage<AberrationFitProgram::OddData>> oddData_thread(num_threads);
				
				const AberrationFitProgram::EvenData evenZero({0.0, 0.0, 0.0, 0.0, 0.0});
				const AberrationFitProgram::OddData oddZero({0.0, dComplex(0.0, 0.0)});
				
				BufferedImage<AberrationFitProgram::EvenData> evenData(sh,s);
				evenData.fill(evenZero);
				
				BufferedImage<AberrationFitProgram::OddData> oddData(sh,s);
				oddData.fill(oddZero);
				
				for (int th = 0; th < num_threads; th++)
				{
					evenData_thread[th] = BufferedImage<AberrationFitProgram::EvenData>(sh,s);
					evenData_thread[th].fill(evenZero);
					
					oddData_thread[th] = BufferedImage<AberrationFitProgram::OddData>(sh,s);
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
					
					AberrationFitProgram::considerParticle(
						particles[t][p], tomogram, referenceMap, dataSet, flip_value, freqWeights,
						f, f, 
						evenData_thread[th], oddData_thread[th]);
				}
				
				Log::endProgress();
				
				for (int th = 0; th < num_threads; th++)
				{
					evenData += evenData_thread[th];
					oddData += oddData_thread[th];
				}

				{
					std::vector<double> coeffs = AberrationFitProgram::solveEven(
							evenData, 2, tomogram.optics.pixelSize, 
							outDir+"diag_"+ZIO::itoa(t)+","+ZIO::itoa(f)+"_", diag);
					
					CTF ctf = tomogram.centralCTFs[f];
					
					//evenData.fill(evenZero);
					
					std::vector<double> coeffs0 = ZernikeHelper::convertSymmetrical(ctf);
					
					const int cc = coeffs0.size() < coeffs.size()? coeffs0.size() : coeffs.size();
					
					for (int i = 0; i < cc; i++)
					{
						coeffs0[i] += coeffs[i];
					}
					
					ZernikeHelper::OldCtfBasis newCtf = ZernikeHelper::convertSymmetrical(
							coeffs0, tomogram.optics.voltage);
					
					ctf.DeltafU = newCtf.defocusU;
					ctf.DeltafV = newCtf.defocusV;
					ctf.azimuthal_angle = newCtf.astigAzimuth_deg;
					
					// Q0, Cs and phase_shift are currently not being handled per particle.
					
					tomogramSet.setCtf(t, f, ctf);
				}
			}
			
			Log::endSection();
			
		} // all frames

		Log::endSection();
		
	} // all tomograms
	
	tomogramSet.write(outDir+"tomograms.star");
}

BufferedImage<double> DefocusRefinementProgram::computeOffsetCost(
		int f, 
		double z0, double z1, int steps, 
		const ParticleSet* dataSet,
		std::vector<int>& particles, int max_particles,
		const Tomogram& tomogram,
		std::vector<BufferedImage<fComplex>>& referenceFS,
		const BufferedImage<float>& freqWeights,
		bool flip_value, 
		double handedness,
		int num_threads)
{
	const double deltaStep = (z1 - z0) / (double) (steps - 1);
	const int fc = tomogram.stack.zdim;
	const int pc = max_particles;
	const double pixelSize = tomogram.optics.pixelSize;
	const int s = referenceFS[0].ydim;
	const int sh = s/2 + 1;
	const double sh2 = sh * (double) sh;
	
	
	std::vector<BufferedImage<double>> 
			CC(num_threads);
	
	for (int th = 0; th < num_threads; th++)
	{
		CC[th] = BufferedImage<double>(pc,steps);
	}
	
	#if TIMING			
		Timer timer;
		int time_extract = timer.setNew("extraction");
		int time_fwd_proj = timer.setNew("fwd. projection");
		int time_CTFdraw = timer.setNew("CTF drawing");
		int time_NCCcalc = timer.setNew("NCC calculation");				
	#endif
	
		
	Log::beginProgress(
		"Evaluating defocus offsets from "
		+ZIO::itoa(z0)+"Å to "+ZIO::itoa(z1)+"Å in "
		+ZIO::itoa(steps)+" steps of "+ZIO::itoa(deltaStep)+"Å", 
		pc / num_threads);
					
	#pragma omp parallel for num_threads(num_threads)		
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();
		
		if (th == 0) Log::updateProgress(p);
		
		const int part_id = particles[p];	
		
		#if TIMING
			if (th==0) timer.tic(time_extract);
		#endif						
		
		const d3Vector pos = dataSet->getPosition(part_id);
		const std::vector<d3Vector> traj = dataSet->getTrajectoryInPixels(part_id, fc, pixelSize);
		
		d4Matrix projCut;	
		
		BufferedImage<fComplex> observation(sh,s);
		
		TomoExtraction::extractFrameAt3D_Fourier(
				tomogram.stack, f, s, 1.0, tomogram.projectionMatrices[f], traj[f],
				observation, projCut, 1, false, true);
		
		#if TIMING
			if (th==0) timer.toc(time_extract);
		#endif
				
		
		#if TIMING
			if (th==0) timer.tic(time_fwd_proj);
		#endif
								
		BufferedImage<fComplex> prediction = Prediction::predictFS(
				part_id, dataSet, projCut, s, referenceFS);
						
		#if TIMING
			if (th==0) timer.toc(time_fwd_proj);
		#endif
			
		CTF ctf_part_0 = tomogram.getCtf(f, pos);		
		CTF ctf_part = ctf_part_0;
		
		BufferedImage<float> CTFimage(sh,s);
								
		for (int di = 0; di < steps; di++)
		{
			const double deltaZ = z0 + di * deltaStep; 
			
			#if TIMING
				if (th==0) timer.tic(time_CTFdraw);
			#endif
					
			ctf_part.DeltafU = ctf_part_0.DeltafU + deltaZ;
			ctf_part.DeltafV = ctf_part_0.DeltafV + deltaZ;
			
			ctf_part.initialise();
			
			ctf_part.draw_fast(s, s, pixelSize, &CTFimage[0]);
			
			#if TIMING
				if (th==0) timer.toc(time_CTFdraw);
			#endif
			
				
			const float scale = flip_value? -1.f : 1.f;
			
			double CCp = 0.0;
			
			#if TIMING
				if (th==0) timer.tic(time_NCCcalc);
			#endif
			
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				const float c = scale * CTFimage(x,y);
				const float wg = freqWeights(x,y,f);
				
				const double xx = x;
				const double yy = y < sh/2? y : y - s;
				const double r2 = xx * xx + yy * yy;
				
				if (r2 < sh2)
				{
					const fComplex zp = c * wg * prediction(x,y);
					const fComplex zv = observation(x,y);
					
					CCp += zp.real * zv.real + zp.imag * zv.imag;
				}
			}
			
			#if TIMING
				if (th==0) timer.toc(time_NCCcalc);
			#endif
			
			CC[th](p, di) = CCp;
		}
	}
	
	Log::endProgress();
	
	
	BufferedImage<double> CC_out(pc,steps);
	
	CC_out.fill(0.0);
	
	for (int th = 0; th < num_threads; th++)
	{
		CC_out += CC[th];
	}
				
	#if TIMING
		timer.printTimes(true);
	#endif
		
	return CC_out;
}

DefocusRefinementProgram::DefocusFit DefocusRefinementProgram::findDefocus(
		int f,  
		double minDelta, 
		double maxDelta,
		int steps, int group_count, double sigma_input,
		const ParticleSet* dataSet,
		std::vector<int>& particles, int max_particles,
		const Tomogram& tomogram,
		std::vector<BufferedImage<fComplex>>& referenceFS,
		const BufferedImage<float>& freqWeights,
		bool flip_value, 
		int num_threads
		)
{
	const int pc = max_particles;
	const double deltaStep = (maxDelta - minDelta) / (double) (steps - 1);
	const bool regularise = sigma_input > 0.0;
	
	DefocusFit out;
	
	Log::beginSection("Estimating coarse defoci");
	
	BufferedImage<double> CC = computeOffsetCost(
		f, minDelta, maxDelta, steps, 
		dataSet, particles, pc, tomogram,
		referenceFS, freqWeights,   
		flip_value, tomogram.handedness, num_threads);
		
	out.totalCost = std::vector<double>(steps, 0.0);
	out.offsets.resize(steps);
	out.costByGroup.resize(group_count);
	
	for (int g = 0; g < group_count; g++)
	{
		out.costByGroup[g] = std::vector<double>(steps, 0.0);
	}
	
	for (int di = 0; di < steps; di++)
	{
		std::vector<double> ccByGroup(group_count, 0.0);
		double totalCC = 0.0;
				
		for (int p = 0; p < pc; p++)
		{
			const int g = p % group_count;
			
			totalCC += CC(p, di);
			ccByGroup[g] += CC(p, di);
		}
		
		for (int g = 0; g < group_count; g++)
		{
			out.costByGroup[g][di] = -ccByGroup[g];
		}
		
		out.offsets[di] = minDelta + di * deltaStep;
		out.totalCost[di] = -totalCC;		
	}
	
	
	std::vector<double> bestDeltaZbyGroup(group_count);
	
	for (int g = 0; g < group_count; g++)
	{		
		double minCost = std::numeric_limits<double>::max();
		
		for (int di = 0; di < steps; di++)
		{
			const double deltaZ = minDelta + di * deltaStep;
			
			if (out.costByGroup[g][di] < minCost)
			{
				minCost = out.costByGroup[g][di];
				bestDeltaZbyGroup[g] = deltaZ;
			}
		}
	}
	
	double bestDeltaZ(0);
	
	{
		double minCost = std::numeric_limits<double>::max();
		
		for (int di = 0; di < steps; di++)
		{
			const double deltaZ = minDelta + di * deltaStep;
			
			if (out.totalCost[di] < minCost)
			{
				minCost = out.totalCost[di];
				bestDeltaZ = deltaZ;
			}
		}
	}
	
	double var(0.0);
	int outliers(0), outliersLeft(0), outliersRight(0);
	
	for (int g = 0; g < group_count; g++)
	{
		const double d = bestDeltaZbyGroup[g] - bestDeltaZ;
		var += d*d;
		
		if (bestDeltaZbyGroup[g] == minDelta)
		{
			outliersLeft++;
			outliers++;
		}
		
		if (bestDeltaZbyGroup[g] == maxDelta)
		{
			outliersRight++;
			outliers++;
		}
	}
	
	var /= (group_count - 1);
	
	const int group_size = pc / group_count;
	const double varOpt = var * group_size / pc;
	const double stdDevOpt = sqrt(varOpt);
	
	
	std::string statusStr;
	
	if (regularise)
	{
		const double s02 = sigma_input * sigma_input;
		const double sf2 = varOpt;
		const double deltaZ_reg = s02 * bestDeltaZ / (sf2 + s02);
	
		statusStr = ZIO::itoa(bestDeltaZ) + " +- " + ZIO::itoa(sqrt(varOpt))
			  + "; after regularisation: " + ZIO::itoa(deltaZ_reg) + "  (";
	}
	else
	{
		statusStr = ZIO::itoa(bestDeltaZ) + " +- " + ZIO::itoa(sqrt(varOpt))
			  + " (";
	}
	
	
	out.value = bestDeltaZ;
	out.stdDev = stdDevOpt;
	
	
	if (bestDeltaZ > 0.0) statusStr = "+"+statusStr;
	
	if (outliers == 0)
	{
		statusStr += "no outliers";
	}
	else if (outliers == 1)
	{
		statusStr += "1 outlier";
	}
	else 
	{
		statusStr += ZIO::itoa(outliers) + " outliers";
	}
	
	statusStr += ")";
	
	Log::print("Best offset:  " + statusStr);
	
	
	if (outliers == group_count && (outliersLeft == outliers || outliersRight == outliers))
	{
		Log::warn("Warning: all defoci have converged to the same side of the search range.");
		Log::extend("The odds of this happening randomly are 1 in " 
				  + ZIO::itoa( (long int) (pow(2.0, group_count)) ) + ".");
		Log::extend("The search range might not be wide enough.");			 
	}
	
	Log::endSection();
	
	if (stdDevOpt < 2.0 * deltaStep)
	{
		Log::beginSection("Refining result");
		
		const double step0 = stdDevOpt / 2.0;
		
		double bestDeltaZ_fine = bestDeltaZ;	
		
		const double diam = 2.0 * deltaStep;
		const double minStep = 2.0;
		
		int steps_fine;
		
		if (step0 > minStep)
		{
			steps_fine = (int) std::round(diam / step0) + 1;
		}
		else
		{
			steps_fine = (int) std::round(diam / minStep) + 1;
		}
		
		if (steps_fine < 5) steps_fine = 5;
		
		const double delta_step_fine = diam / (double) (steps_fine - 1);
						
		BufferedImage<double> CC_fine = computeOffsetCost(
				f, 
				bestDeltaZ - deltaStep, 
				bestDeltaZ + deltaStep, 
				steps_fine, 
				dataSet, particles, pc, tomogram,
				referenceFS, freqWeights,   
				flip_value, tomogram.handedness, num_threads);
						
		
		double minCost_fine = std::numeric_limits<double>::max();
		
		for (int di = 0; di < steps_fine; di++)
		{
			const double deltaZ = bestDeltaZ - deltaStep + di * delta_step_fine;
			
			double totalCC = 0.0;
					
			for (int p = 0; p < pc; p++)
			{
				const int g = p % group_count;
				
				totalCC += CC_fine(p, di);
			}
			
			const double cost = -totalCC;
			
			if (cost < minCost_fine)
			{
				minCost_fine = cost;
				bestDeltaZ_fine = deltaZ;
			}
		}
		
		
		if (regularise)
		{
			const double s02 = sigma_input * sigma_input;
			const double sf2 = varOpt;
			
			const double deltaZ = s02 * bestDeltaZ_fine / (sf2 + s02);
			
			Log::print("Refined:  " + ZIO::itoa(bestDeltaZ_fine) 
					  + "; after regularisation: " + ZIO::itoa(deltaZ));
		}
		else
		{
			Log::print("Refined:  " + ZIO::itoa(bestDeltaZ_fine));	
		}
		
		out.value = bestDeltaZ_fine;
		
		Log::endSection();
	}
	
	return out;
}

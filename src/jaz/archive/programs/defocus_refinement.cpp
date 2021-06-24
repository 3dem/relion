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
#include <src/jaz/math/Zernike_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optics/magnification_helper.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <iostream>
#include <src/time.h>

#define TIMING 0


using namespace gravis;
using namespace aberration;


DefocusRefinementProgram::DefocusRefinementProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	
	try
	{
		_readParams(parser);
		
		int def_section = parser.addSection("Defocus refinement options");

		do_scanDefocus = !parser.checkOption("--no_scan", "Skip accelerated defocus scan");
		do_slowScan = parser.checkOption("--slow_scan", "Perform a slow, brute-force defocus scan instead");
		do_refineFast = !parser.checkOption("--slow_scan_only", "Only perform a brute-force scan");

		do_defocus = !parser.checkOption("--no_defocus", "Do not refine the defocus itself, only its slope");

		do_refineAstigmatism = !parser.checkOption("--no_astigmatism", "Do not refine the astigmatism");
		do_plotAstigmatism = parser.checkOption("--plot_astigmatism", "Plot the astigmatism cost function");
		do_slopeFit = parser.checkOption("--fit_slope", "Fit the slope of defocus over depth");
		do_perTomogramSlope = parser.checkOption("--tomo_slope", "Fit the defocus slope separately for each tomogram (not recommended)");
		max_slope_dose = textToDouble(parser.getOption("--max_slope_dose", "Maximum dose for a frame to be considered for defocus slope estimation", "20"));

		max_particles = textToInteger(parser.getOption("--max", "Max. number of particles to consider per tomogram", "-1"));
		group_count = textToInteger(parser.getOption("--g", "Number of independent groups", "10"));
		sigma_input = textToDouble(parser.getOption("--sig0", "Std. dev. of initial defoci (negative to turn off regularisation)", "-1"));
		do_regularise = sigma_input > 0.0;
		
		minDelta = textToDouble(parser.getOption("--d0", "Min. defocus offset to test [Å]", "-300"));
		maxDelta = textToDouble(parser.getOption("--d1", "Max. defocus offset to test [Å]", "300"));
		deltaSteps = textToInteger(parser.getOption("--ds", "Number of defocus steps in-between", "100"));
		
		do_clearAstigmatism = parser.checkOption("--ca", "Clear the current astigmatism estimate");
		
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

	if (do_defocus)
	{
		if (!do_refineFast)
		{
			do_slowScan = true;
		}
	}
	else
	{
		if (!do_slopeFit)
		{
			REPORT_ERROR("You need to either refine the defoci or their slopes");
		}
	}
}

void DefocusRefinementProgram::run()
{
	Log::beginSection("Initialising");
	
	RefinementProgram::init();
		
	const int tc = particles.size();
	const bool flip_value = true;

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize);
	
	Log::endSection();

	const double min_slope = 0.95;
	const double max_slope = 1.05;
	const int slope_steps = 2;

	std::vector<d3Vector> globalSlopeCost(slope_steps, d3Vector(0.0, 0.0, 0.0));

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
		
		particleSet.checkTrajectoryLengths(
				particles[t][0], usedParticleCount, fc, "DefocusRefinementProgram::run");
		
		
		if (do_clearAstigmatism)
		{
			for (int f = 0; f < fc; f++)
			{
				CTF& ctf = tomogram.centralCTFs[f];
				const double z0 = (ctf.DeltafU + ctf.DeltafV) / 2.0;
				
				ctf.DeltafU = z0;
				ctf.DeltafV = z0;
				ctf.initialise();
			}
		}
				
		const int first_frame = specified_first_frame;
		const int last_frame = (specified_last_frame > 0 && specified_last_frame < fc)? specified_last_frame : fc-1;
		
		std::vector<double> defocusOffset(fc, 0.0), offsetStdDev(fc, 0.0);
		
		/*
		  TODO: find best handedness from current defocus over all frames
		*/
		
		BufferedImage<float> freqWeights = computeFrequencyWeights(
			tomogram, true, 0.0, 0.0, true, num_threads);

		BufferedImage<float> doseWeight = tomogram.computeDoseWeight(boxSize,1);

		if (do_defocus)
		{
			for (int f = first_frame; f <= last_frame; f++)
			{
				Log::beginSection("Frame " + ZIO::itoa(f+1));

				if (do_slowScan)
				{
					DefocusFit defocus = findDefocus(
						f, minDelta, maxDelta, deltaSteps,
						group_count, sigma_input,
						particleSet, particles[t], usedParticleCount,
						tomogram, aberrationsCache, referenceMap.image_FS,
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

						costOut.precision(12);

						for (int di = 0; di < deltaSteps; di++)
						{
							costOut << defocus.offsets[di] << " " << defocus.totalCost[di] << '\n';
						}

						costOut << '\n';
					}


					CTF& ctf = tomogram.centralCTFs[f];

					const double s02 = sigma_input * sigma_input;
					const double sf2 = offsetStdDev[f] * offsetStdDev[f];

					const double deltaZ = do_regularise?
						s02 * defocusOffset[f] / (sf2 + s02) :
						defocusOffset[f];

					ctf.DeltafU += deltaZ;
					ctf.DeltafV += deltaZ;

					ctf.initialise();

					tomogramSet.setCtf(t,f,ctf);
				}

				if (do_refineFast)
				{
					const int s = referenceMap.image_real[0].xdim;
					const int sh = s/2 + 1;
					const int pc = particles[t].size();

					std::vector<BufferedImage<EvenData>> evenData_thread(num_threads);
					std::vector<BufferedImage<OddData>> oddData_thread(num_threads);

					const EvenData evenZero({0.0, 0.0, 0.0, 0.0, 0.0});
					const OddData oddZero({0.0, dComplex(0.0, 0.0)});

					BufferedImage<EvenData> evenData(sh,s);
					evenData.fill(evenZero);

					BufferedImage<OddData> oddData(sh,s);
					oddData.fill(oddZero);

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

						AberrationFit::considerParticle(
							particles[t][p], tomogram, referenceMap, particleSet,
							aberrationsCache, flip_value, freqWeights, doseWeight,
							f, f,
							evenData_thread[th], oddData_thread[th]);
					}

					Log::endProgress();


					for (int th = 0; th < num_threads; th++)
					{
						evenData += evenData_thread[th];
						oddData += oddData_thread[th];
					}

					CTF ctf0 = tomogram.centralCTFs[f];
					CTF ctf_dz = ctf0;

					double bestDeltaZ = 0;

					if (do_scanDefocus)
					{
						Log::print("Scanning for optimal defocus");

						bestDeltaZ = scanForDefocus(
								evenData, tomogram.optics.pixelSize, tomogram.centralCTFs[f],
								minDelta, maxDelta, deltaSteps);

						ctf_dz.DeltafU = ctf0.DeltafU + bestDeltaZ;
						ctf_dz.DeltafV = ctf0.DeltafV + bestDeltaZ;
					}

					CTF ctf1 = ctf_dz;

					if (do_refineAstigmatism)
					{
						Log::print("Refining astigmatism");

						EvenSolution solution = AberrationFit::solveEven(evenData);

						const double pixelSize = tomogram.optics.pixelSize;

						d3Vector astig = findAstigmatism(
							solution, ctf0, bestDeltaZ, pixelSize, 1.0);

						if (do_plotAstigmatism)
						{
							BufferedImage<double> astigPlot = plotAstigmatism(
										solution, ctf0, bestDeltaZ, 100.0, pixelSize, 32);

							astigPlot.write(
								outDir+"astig_cost_"+ZIO::itoa(t)+"_"+ZIO::itoa(f)+".mrc");
						}

						ctf1.DeltafU = astig[0];
						ctf1.DeltafV = astig[1];
						ctf1.azimuthal_angle = astig[2];
					}

					tomogramSet.setCtf(t, f, ctf1);

				}

				Log::endSection();

			} // all frames
		} // do_slow_scan || do_refine_fast

		if (do_slopeFit)
		{
			Log::beginSection("Fitting defocus slope");

			std::vector<d3Vector> tomogramSlopeCost = computeSlopeCost(
				max_slope_dose, min_slope, max_slope, slope_steps,
				particleSet, particles[t], pc0, tomogram, aberrationsCache,
				referenceMap.image_FS, freqWeights,
				flip_value, num_threads);

			if (diag)
			{
				writeSlopeCost(tomogramSlopeCost, outDir+"t_"+ZIO::itoa(t)+"_slope.dat");
			}

			if (do_perTomogramSlope)
			{
				const double slope = min_slope -
						tomogramSlopeCost[0][2] * (max_slope - min_slope)
						/ (tomogramSlopeCost[1][2] - tomogramSlopeCost[0][2]);

				tomogramSet.setDefocusSlope(t, slope);
			}
			else
			{
				for (int i = 0; i < tomogramSlopeCost.size(); i++)
				{
					globalSlopeCost[i][0]  = tomogramSlopeCost[i][0];
					globalSlopeCost[i][1] += tomogramSlopeCost[i][1];
					globalSlopeCost[i][2] += tomogramSlopeCost[i][2];
				}
			}

			Log::endSection();
		}

		Log::endSection();
		
	} // all tomograms

	if (do_slopeFit)
	{
		if (diag)
		{
			writeSlopeCost(globalSlopeCost, outDir+"slope.dat");
		}

		if (!do_perTomogramSlope)
		{
			for (int t = 0; t < tc; t++)
			{
				const double slope = min_slope -
						globalSlopeCost[0][2] * (max_slope - min_slope)
						/ (globalSlopeCost[1][2] - globalSlopeCost[0][2]);

				tomogramSet.setDefocusSlope(t, slope);
			}
		}
	}

	tomogramSet.write(outDir+"tomograms.star");

	optimisationSet.tomograms = outDir+"tomograms.star";
	optimisationSet.write(outDir+"optimisation_set.star");
}

void DefocusRefinementProgram::writeSlopeCost(
		const std::vector<d3Vector>& cost,
		const std::string& filename)
{
	double minVal = cost[0][1];

	for (int i = 0; i < cost.size(); i++)
	{
		if (cost[i][1] < minVal)
		{
			minVal = cost[i][1];
		}
	}

	std::ofstream slopeFile(filename);

	for (int i = 0; i < cost.size(); i++)
	{
		slopeFile.precision(12);
		slopeFile << cost[i][0] << ' ' << (cost[i][1] - minVal) << '\n';
	}

	slopeFile << '\n';

	for (int i = 0; i < cost.size(); i++)
	{
		slopeFile.precision(12);
		slopeFile << cost[i][0] << ' ' << cost[i][2] << '\n';
	}
}

BufferedImage<double> DefocusRefinementProgram::computeOffsetCost(
		int f,
		double z0, double z1, int steps,
		const ParticleSet& dataSet,
		std::vector<ParticleIndex>& particles, int max_particles,
		const Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		std::vector<BufferedImage<fComplex>>& referenceFS,
		const BufferedImage<float>& freqWeights,
		bool flip_value,
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

		const ParticleIndex part_id = particles[p];

		#if TIMING
			if (th==0) timer.tic(time_extract);
		#endif

		const d3Vector pos = dataSet.getPosition(part_id);
		const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(part_id, fc, pixelSize);

		d4Matrix projCut;

		BufferedImage<fComplex> observation(sh,s);

		TomoExtraction::extractFrameAt3D_Fourier(
				tomogram.stack, f, s, 1.0, tomogram, traj[f],
				observation, projCut, 1, true);

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

		const int og = dataSet.getOpticsGroup(part_id);

		const BufferedImage<double>* gammaOffset =
				aberrationsCache.hasSymmetrical? &aberrationsCache.symmetrical[og] : 0;

		for (int di = 0; di < steps; di++)
		{
			const double deltaZ = z0 + di * deltaStep;

			#if TIMING
				if (th==0) timer.tic(time_CTFdraw);
			#endif

			ctf_part.DeltafU = ctf_part_0.DeltafU + deltaZ;
			ctf_part.DeltafV = ctf_part_0.DeltafV + deltaZ;

			ctf_part.initialise();

			ctf_part.draw_fast(s, s, pixelSize, gammaOffset, &CTFimage[0]);

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
				const double yy = y < s/2? y : y - s;
				const double r2 = xx * xx + yy * yy;

				if (r2 < sh2)
				{
					const fComplex zp = c * prediction(x,y);
					const fComplex zo = observation(x,y);

					CCp += wg * (zp - zo).norm();
				}
			}

			#if TIMING
				if (th==0) timer.toc(time_NCCcalc);
			#endif

			CC[th](p, di) = CCp / (s * s);
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


std::vector<d3Vector> DefocusRefinementProgram::computeSlopeCost(
		double max_dose,
		double m0, double m1, int steps,
		const ParticleSet& dataSet,
		std::vector<ParticleIndex>& particles, int max_particles,
		const Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
		std::vector<BufferedImage<fComplex>>& referenceFS,
		const BufferedImage<float>& freqWeights,
		bool flip_value,
		int num_threads)
{
	const double deltaStep = (m1 - m0) / (double) (steps - 1);
	const int fc = tomogram.stack.zdim;
	const int pc = max_particles;
	const double pixelSize = tomogram.optics.pixelSize;
	const double ba = boxSize * pixelSize;
	const int s = referenceFS[0].ydim;
	const int sh = s/2 + 1;
	const double sh2 = sh * (double) sh;


	const d3Vector centre_of_mass = tomogram.computeCentreOfMass(dataSet, particles);


	std::vector<std::vector<double>>
			cost_per_thread(num_threads),
			slope_per_thread(num_threads);

	for (int th = 0; th < num_threads; th++)
	{
		cost_per_thread[th] = std::vector<double>(steps, 0.0);
		slope_per_thread[th] = std::vector<double>(steps, 0.0);
	}

	std::vector<int> good_frames;

	for (int f = 0; f < fc; f++)
	{
		if (tomogram.cumulativeDose[f] < max_dose)
		{
			good_frames.push_back(f);
		}
	}


	Log::beginProgress(
		"Evaluating defocus slopes from "
		+ZIO::itoa(m0)+" to "+ZIO::itoa(m1)+" using "+ZIO::itoa(good_frames.size())+" frames",
		good_frames.size() * pc / num_threads);


	for (int ff = 0; ff < good_frames.size(); ff++)
	{
		const int f = good_frames[ff];

		#pragma omp parallel for num_threads(num_threads)
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();

			if (th == 0) Log::updateProgress(ff * pc / num_threads + p);

			const ParticleIndex part_id = particles[p];

			const d3Vector pos = dataSet.getPosition(part_id);
			const std::vector<d3Vector> traj = dataSet.getTrajectoryInPixels(
						part_id, fc, pixelSize);

			d4Matrix projCut;

			BufferedImage<fComplex> observation(sh,s);

			TomoExtraction::extractFrameAt3D_Fourier(
					tomogram.stack, f, s, 1.0,
					tomogram, traj[f],
					observation, projCut, 1, true);

			BufferedImage<fComplex> prediction = Prediction::predictFS(
					part_id, dataSet, projCut, s, referenceFS);

			double dz0 = tomogram.getDepthOffset(f, pos);

			CTF ctf0 = tomogram.centralCTFs[f];

			CTF ctf_part = ctf0;
			ctf_part.initialise();

			const int og = dataSet.getOpticsGroup(part_id);

			const BufferedImage<double>* gammaOffset =
					aberrationsCache.hasSymmetrical? &aberrationsCache.symmetrical[og] : 0;

			const double avg_offset = tomogram.getDepthOffset(f, centre_of_mass);

			BufferedImage<float> CTFimage(sh,s);

			for (int di = 0; di < steps; di++)
			{
				const double m = m0 + di * deltaStep;
				const double deltaZ = avg_offset + m * (dz0 - avg_offset);
				const double deltaF = tomogram.handedness * tomogram.optics.pixelSize * deltaZ;
				const double ddeltaF_dm = tomogram.handedness * tomogram.optics.pixelSize * (dz0 - avg_offset);

				ctf_part.DeltafU = ctf0.DeltafU + deltaF;
				ctf_part.DeltafV = ctf0.DeltafV + deltaF;

				ctf_part.initialise();

				std::vector<double> K_ctf = ctf_part.getK();
				ctf_part.draw_fast(s, s, pixelSize, gammaOffset, &CTFimage[0]);


				const float scale = flip_value? -1.f : 1.f;

				double cost_p = 0.0;
				double slope_p = 0.0;

				for (int y = 0; y < s;  y++)
				for (int x = 0; x < sh; x++)
				{
					const double xp = x;
					const double yp = y < s/2? y : y - s;

					const double xa = xp / ba;
					const double ya = yp / ba;

					const double gamma = ctf_part.getLowOrderGamma(xa,ya);

					const float c = -scale * sin(gamma);

					const float wg = freqWeights(x,y,f);

					const fComplex pred = prediction(x,y);

					const fComplex dF = c * pred - observation(x,y);

					const float dgamma_ddeltaF = -K_ctf[1] * (xa*xa + ya*ya);
					const float dgamma_dmag = dgamma_ddeltaF * ddeltaF_dm;

					const float dc_dmag = -scale * cos(gamma) * dgamma_dmag;

					cost_p += wg * dF.norm();

					const fComplex dFp_dmag = dc_dmag * pred;

					slope_p += 2.0 * wg * (dF.real * dFp_dmag.real + dF.imag * dFp_dmag.imag);
				}

				cost_per_thread[th][di] += cost_p / (s * s);
				slope_per_thread[th][di] += slope_p / (s * s);
			}
		}
	}

	Log::endProgress();

	std::vector<d3Vector> out(steps, d3Vector(0.0, 0.0, 0.0));

	for (int di = 0; di < steps; di++)
	{
		out[di][0] = m0 + di * deltaStep;

		for (int th = 0; th < num_threads; th++)
		{
			out[di][1] += cost_per_thread[th][di];
			out[di][2] += slope_per_thread[th][di];
		}
	}

	return out;
}

DefocusRefinementProgram::DefocusFit DefocusRefinementProgram::findDefocus(
		int f,  
		double minDelta, 
		double maxDelta,
		int steps, int group_count, double sigma_input,
		const ParticleSet& dataSet,
		std::vector<ParticleIndex>& particles, int max_particles,
		const Tomogram& tomogram,
		const AberrationsCache& aberrationsCache,
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
	
	BufferedImage<double> cost = computeOffsetCost(
		f, minDelta, maxDelta, steps, 
		dataSet, particles, pc, tomogram, aberrationsCache,
		referenceFS, freqWeights,
		flip_value, num_threads);
		
	out.totalCost = std::vector<double>(steps, 0.0);
	out.offsets.resize(steps);
	out.costByGroup.resize(group_count);
	
	for (int g = 0; g < group_count; g++)
	{
		out.costByGroup[g] = std::vector<double>(steps, 0.0);
	}
	
	for (int di = 0; di < steps; di++)
	{
		std::vector<double> costByGroup(group_count, 0.0);
		double totalCost = 0.0;

		for (int p = 0; p < pc; p++)
		{
			const int g = p % group_count;

			totalCost += cost(p, di);
			costByGroup[g] += cost(p, di);
		}
		
		for (int g = 0; g < group_count; g++)
		{
			out.costByGroup[g][di] = costByGroup[g];
		}
		
		out.offsets[di] = minDelta + di * deltaStep;
		out.totalCost[di] = totalCost;
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
						
		BufferedImage<double> cost_fine = computeOffsetCost(
				f, 
				bestDeltaZ - deltaStep, 
				bestDeltaZ + deltaStep, 
				steps_fine, 
				dataSet, particles, pc, tomogram, aberrationsCache,
				referenceFS, freqWeights,   
				flip_value, num_threads);
						
		
		double minCost_fine = std::numeric_limits<double>::max();
		
		for (int di = 0; di < steps_fine; di++)
		{
			const double deltaZ = bestDeltaZ - deltaStep + di * delta_step_fine;
			
			double totalCost = 0.0;
					
			for (int p = 0; p < pc; p++)
			{
				const int g = p % group_count;
				
				totalCost += cost_fine(p, di);
			}
			
			const double cost = totalCost;
			
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

double DefocusRefinementProgram::scanForDefocus(
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
				EvenData d = evenData(x,y);

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

gravis::d3Vector DefocusRefinementProgram::findAstigmatism(
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

BufferedImage<double> DefocusRefinementProgram::plotAstigmatism(
		const EvenSolution& solution,
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

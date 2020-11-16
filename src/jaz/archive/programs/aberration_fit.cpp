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
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/optics/aberration_fit.h>
#include <omp.h> 

using namespace gravis;
using namespace aberration;


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
				
		int defocus_section = parser.addSection("Aberration fit options");

		do_even = !parser.checkOption("--no_symm", "Do not fit symmetrical aberrations");
		do_odd = !parser.checkOption("--no_antisymm", "Do not fit antisymmetrical aberrations");

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
	if (!do_even && !do_odd)
	{
		// This might happen if the program is called from a script.
		// Let's make sure a particles.star is still written out.

		Log::warn("Estimating neither symmetrical nor antisymmetrical aberrations: there is nothing to estimate.");

		particleSet.write(outDir+"particles.star");

		return;
	}

	if (do_even && do_odd)
	{
		Log::print("Estimating both symmetrical and antisymmetrical aberrations");
	}
	else if (do_even)
	{
		Log::print("Estimating symmetrical aberrations only");
	}
	else // do_odd
	{
		Log::print("Estimating antisymmetrical aberrations only");
	}

	Log::beginSection("Initialising");
	
	RefinementProgram::init();
		
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int tc = particles.size();
	const int gc = particleSet.numberOfOpticsGroups();
	const bool flip_value = true;
	
	Log::endSection();

	// check for pixel size consistency!

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

	AberrationsCache aberrationsCache(particleSet.optTable, boxSize);

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

		const int fc = tomogram.frameCount;
		const int first_frame = specified_first_frame;
		const int last_frame = (specified_last_frame > 0 && specified_last_frame < fc)? specified_last_frame : fc-1;


		BufferedImage<float> frqWeight = computeFrequencyWeights(
			tomogram, true, 1.0, 0.0, true, num_threads);

		BufferedImage<float> doseWeight = tomogram.computeDoseWeight(s,1);
		
		if (diag)
		{
			frqWeight.write(diagPrefix + "_noise_weight.mrc");
			referenceMap.freqWeight.write(diagPrefix + "_FSC_weight.mrc");
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
				
		Log::beginProgress("Accumulating evidence", pc/num_threads);
		
		#pragma omp parallel for num_threads(num_threads)
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();
			
			if (th == 0)
			{
				Log::updateProgress(p);
			}

			const int g = particleSet.getOpticsGroup(particles[t][p]);
			
			AberrationFit::considerParticle(
				particles[t][p], tomogram, referenceMap, particleSet,
				aberrationsCache, flip_value, frqWeight, doseWeight,
				first_frame, last_frame,
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
		}
		
		Log::endSection();
	}

	for (int g = 0; g < gc; g++)
	{
		if (do_even)
		{
			std::vector<double> initialEven(Zernike::numberOfEvenCoeffs(n_even), 0.0);

			std::vector<double> evenCoeffs = AberrationFit::solveAndFitEven(
				evenData_perGroup[g], n_even, initialEven,
				lastPixelSize, outDir + ZIO::itoa(g+1) + "_", true);

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

		if (do_odd)
		{
			std::vector<double> initialOdd(Zernike::numberOfOddCoeffs(n_odd), 0.0);

			std::vector<double> oddCoeffs = AberrationFit::solveAndFitOdd(
				oddData_perGroup[g], n_odd, initialOdd,
				lastPixelSize, outDir + ZIO::itoa(g+1) + "_", true);

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

	particleSet.write(outDir+"particles.star");

	optimisationSet.particles = outDir+"particles.star";
	optimisationSet.write(outDir+"optimisation_set.star");
}

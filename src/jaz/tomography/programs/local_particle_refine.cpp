#include "local_particle_refine.h"
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
#include <src/jaz/optics/tomo_mag_fit.h>
#include <src/jaz/tomography/projection/fwd_projection.h>
#include <src/jaz/tomography/local_particle_refinement.h>
#include <src/jaz/math/Tait_Bryan_angles.h>

#include <src/jaz/math/Euler_angles_relion.h>
#include <src/euler.h>


#include <omp.h>

using namespace gravis;
using namespace aberration;


LocalParticleRefineProgram::LocalParticleRefineProgram(int argc, char *argv[])
:	RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	readParams(parser);
}

void LocalParticleRefineProgram::readParams(IOParser &parser)
{
	try
	{
		_readParams(parser);

		int defocus_section = parser.addSection("Alignment options");


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

void LocalParticleRefineProgram::run()
{
	Log::beginSection("Initialising");

	RefinementProgram::init();

	const int s = boxSize;
	const int sh = s/2 + 1;
	const int tc = particles.size();
	const int gc = dataSet.numberOfOpticsGroups();
	const bool flip_value = true;

	AberrationsCache aberrationsCache(dataSet.optTable, boxSize);

	Log::endSection();

	for (int t = 0; t < 1; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;

		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));
		Log::print("Loading");

		Tomogram tomogram = tomogramSet.loadTomogram(t, true);

		const int fc = tomogram.frameCount;

		BufferedImage<float> freqWeights = computeFrequencyWeights(
			tomogram, true, 0.0, 0.0, num_threads);

		freqWeights.write("DEBUG_freqWeights.mrc");

		dataSet.checkTrajectoryLengths(
				particles[t][0], pc, fc, "LocalParticleRefineProgram::run");

		const int first_frame = specified_first_frame;
		const int last_frame = (specified_last_frame > 0 && specified_last_frame < fc)? specified_last_frame : fc-1;


		for (int p = 0; p < 15; p++)
		{
			LocalParticleRefinement refinement(
					particles[t][p], dataSet, tomogram, referenceMap,
					freqWeights, aberrationsCache, false);

			const std::vector<double> initial {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

			const std::vector<double> optimal = NelderMead::optimize(
						initial, refinement, 2, 0.001, 300, 1, 2, 0.5, 0.5, false);

			{
				for (int j = 0; j < 3; j++)
				{
					std::cout << optimal[j] << "  ";
				}

				std::cout << " : ";

				for (int j = 3; j < 6; j++)
				{
					std::cout << optimal[j] << "  ";
				}

				std::cout << "   ->  ";
			}

			std::cout << "f = " << refinement.f(optimal, 0) << std::endl;
		}

		Log::endSection();
	}
}


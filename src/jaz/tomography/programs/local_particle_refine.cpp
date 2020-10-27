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
#include <src/jaz/optimization/lbfgs.h>
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

	const int tc = particles.size();

	AberrationsCache aberrationsCache(dataSet.optTable, boxSize);

	Log::endSection();

	for (int t = 0; t < tc; t++)
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

		const int pc_max = 36;
		const int data_pad = 256;

		std::vector<double> results(pc * data_pad);


		Log::beginProgress("Aligning particles", pc/num_threads);

		#pragma omp parallel for num_threads(num_threads)
		for (int p = 0; p < pc; p++)
		{
			const int th = omp_get_thread_num();

			if (th == 0)
			{
				Log::updateProgress(p);
			}

			LocalParticleRefinement refinement(
					particles[t][p], dataSet, tomogram, referenceMap,
					freqWeights, aberrationsCache, false);

			const std::vector<double> initial {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

			/*const std::vector<double> optimal = NelderMead::optimize(
						initial, refinement, 2, 0.001, 300, 1, 2, 0.5, 0.5, false);*/

			const std::vector<double> optimal = LBFGS::optimize(
						initial, refinement, false, 300, 1e-5, 1e-4);

			for (int i = 0; i < 6; i++)
			{
				results[p * data_pad + i] = optimal[i];
			}
		}

		Log::endProgress();


		for (int p = 0; p < pc; p++)
		{
			std::vector<double> opt(6);

			for (int i = 0; i < 6; i++)
			{
				opt[i] = results[p * data_pad + i];
			}

			LocalParticleRefinement::applyChange(
				opt, dataSet, particles[t][p], tomogram.optics.pixelSize);
		}

		Log::endSection();
	}

	dataSet.write(outDir+"particles.star");
}


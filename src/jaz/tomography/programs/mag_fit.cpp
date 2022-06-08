#include "mag_fit.h"
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

#include <src/jaz/math/Euler_angles_relion.h>
#include <src/euler.h>


#include <omp.h>

using namespace gravis;
using namespace aberration;


MagFitProgram::MagFitProgram(int argc, char *argv[])
:	RefinementProgram(argc, argv)
{
}

void MagFitProgram::readParams()
{
	IOParser parser;
	parser.setCommandLine(argc, argv);

	_readParams(parser);

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
}

void MagFitProgram::run()
{
	readParams();

	Log::beginSection("Initialising");

		RefinementProgram::init();

		const int s = boxSize;
		const int tc = particles.size();
		const int gc = particleSet.numberOfOpticsGroups();

	Log::endSection();

	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;

		Log::beginSection("Tomogram " + ZIO::itoa(t+1) + " / " + ZIO::itoa(tc));
		Log::print("Loading");

		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		tomogram.validateParticleOptics(particles[t], particleSet);

		const int fc = tomogram.frameCount;

		BufferedImage<float> freqWeights = computeFrequencyWeights(
			tomogram, true, 0.0, 0.0, false, num_threads);

		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, 1.0);

		particleSet.checkTrajectoryLengths(
				particles[t], fc, "MagFitProgram::run");

		TomoAnisoMagFit anisoFit(
			particles[t],
			tomogram,
			particleSet,
			referenceMap,
			freqWeights,
			doseWeights,
			boxSize,
			0,
			fc-1,
			num_threads);

		Log::print("Estimating magnification matrix");

		BufferedImage<Equation2x2> equations = anisoFit.computeEquations();
		d2Matrix magMatrix = MagnificationHelper::solveLinearly(equations);

		std::cout << magMatrix << std::endl;


		Log::endSection();
	}
}


#include "fcc_computation.h"
#include <src/ctf.h>
#include <src/jaz/optics/damage.h>
#include <src/jaz/math/fcc.h>
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
#include <omp.h>


using namespace gravis;


FccProgram::FccProgram(int argc, char *argv[])
	: RefinementProgram(argc, argv)
{
}

void FccProgram::readParams()
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

void FccProgram::run()
{
	readParams();

	RefinementProgram::init();
	
	const int tc = particles.size();
	const int s = boxSize;
	const int sh = s/2 + 1;
	const bool flip_value = true;

	
	std::vector<double> pixelSizes(tc, 0.0);
			
	
	for (int t = 0; t < tc; t++)
	{
		int pc = particles[t].size();
		if (pc == 0) continue;
		
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		pixelSizes[t] = tomogram.optics.pixelSize;

		std::string tag = outDir + tomogram.name;

		BufferedImage<double> FCC3 = FCC::compute3(
				particleSet, particles[t], tomogram, referenceMap.image_FS,
				flip_value, num_threads);
			
		FCC3.write(tag + "_FCC3.mrc");
		FCC::divide(FCC3).write(tag + "_FCC.mrc");

		{
			const int fc = tomogram.frameCount;
			BufferedImage<double> scaleFactor(sh,fc);

			for (int f = 0; f < fc; f++)
			{
				scaleFactor(0,f) = 0.0;

				for (int x = 1; x < sh; x++)
				{
					scaleFactor(x,f) = FCC3(x,f,0) / FCC3(x,f,1);
				}
			}

			scaleFactor.write(tag + "_scaleFactor.mrc");
		}
	}
}

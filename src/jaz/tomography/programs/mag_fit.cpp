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
#include <omp.h>

using namespace gravis;
using namespace aberration;


MagFitProgram::MagFitProgram(int argc, char *argv[])
:	RefinementProgram(argc, argv)
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	readParams(parser);
}

void MagFitProgram::readParams(IOParser &parser)
{
	try
	{
		_readParams(parser);

		int defocus_section = parser.addSection("Alignment options");


		initial_step = textToDouble(parser.getOption("--ins", "Initial step (in %)", "2"));


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

void MagFitProgram::run()
{
	Log::beginSection("Initialising");

	RefinementProgram::init();

	const int s = boxSize;
	const int sh = s/2 + 1;
	const int tc = particles.size();
	const int gc = dataSet.numberOfOpticsGroups();
	const bool flip_value = true;

	Log::endSection();


}


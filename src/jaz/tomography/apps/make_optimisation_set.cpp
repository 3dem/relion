#include <src/args.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/optimisation_set.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	try
	{
		IOParser parser;

		OptimisationSet os;
		std::string outFn;

		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		os.read(
			parser,
			true,          // optimisation set
			true,   false,  // particles
			true,   false,  // tomograms
			true,   false,  // trajectories
			true,   false,  // manifolds
			true,   false); // reference

		outFn = parser.getOption("--o", "Output file name");

		if (parser.checkForErrors())
		{
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}

		os.write(outFn);
	}
	catch (RelionError XE)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

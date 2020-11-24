#include <src/args.h>
#include <src/jaz/tomography/programs/sample_manifold.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	try
	{
		SampleManifoldProgram sp;

		sp.readParameters(argc, argv);
		sp.run();
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

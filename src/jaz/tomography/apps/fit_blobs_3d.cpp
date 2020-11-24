#include <src/args.h>
#include <src/jaz/tomography/programs/fit_blobs_3d.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	try
	{
		FitBlobs3DProgram program;

		program.readParameters(argc, argv);
		program.run();
	}
	catch (RelionError XE)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

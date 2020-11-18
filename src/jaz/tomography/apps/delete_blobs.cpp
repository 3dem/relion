#include <src/args.h>
#include <src/jaz/tomography/programs/delete_blobs.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	DeleteBlobsProgram program;

	try
	{
		program.readParameters(argc, argv);
		program.run();
	}
	catch (RelionError XE)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

#include <src/args.h>
#include <src/jaz/tomography/programs/subtomo.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	try
	{
		SubtomoProgram program;

		program.readParameters(argc, argv);
		program.run();
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;

	return 0;
}

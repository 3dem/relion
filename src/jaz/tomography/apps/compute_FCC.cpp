#include <src/args.h>
#include <src/jaz/tomography/programs/fcc_computation.h>

int main(int argc, char *argv[])
{

	try
	{
		FccProgram program(argc, argv);

		program.run();
	}
	catch (RelionError XE)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

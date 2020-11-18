#include <src/args.h>
#include <src/jaz/tomography/programs/bfactor_fit.h>

int main(int argc, char *argv[])
{

	try
	{
		BfactorFitProgram program(argc, argv);

		program.run();
	}
	catch (RelionError XE)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

#include <src/args.h>
#include <src/jaz/tomography/programs/mag_fit.h>


int main(int argc, char *argv[])
{
	try
	{
		MagFitProgram program(argc, argv);
		program.run();
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

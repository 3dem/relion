#include <src/args.h>
#include <src/jaz/tomography/programs/reconstruct_tomogram.h>


int main(int argc, char *argv[])
{
	try
	{
		TomoBackprojectProgram program;

		program.readParameters(argc, argv);
		program.run();
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

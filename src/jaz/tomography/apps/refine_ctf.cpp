#include <src/args.h>
#include <src/jaz/tomography/programs/ctf_refinement.h>


int main(int argc, char *argv[])
{
	try
	{
		CtfRefinementProgram program(argc, argv);
		program.run();
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

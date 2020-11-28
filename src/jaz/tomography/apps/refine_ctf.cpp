#include <src/args.h>
#include <src/jaz/tomography/programs/ctf_refinement.h>


int main(int argc, char *argv[])
{
	CtfRefinementProgram program(argc, argv);

	try
	{
		program.run();
	}
	catch (RelionError e)
	{
		std::cerr << e << std::endl;

		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

#include <src/args.h>
#include <src/jaz/tomography/programs/align.h>


int main(int argc, char *argv[])
{
	AlignProgram ap(argc, argv);

	try
	{
		ap.run();
	}
	catch (RelionError e)
	{
		std::cerr << e << std::endl;

		return RELION_EXIT_FAILURE;
	}
	
	return RELION_EXIT_SUCCESS;
}

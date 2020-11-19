#include <src/args.h>
#include <src/jaz/tomography/programs/align.h>


int main(int argc, char *argv[])
{
	try
	{
		AlignProgram ap(argc, argv);
		ap.run();
	}
	catch (RelionError XE)
	{
		return RELION_EXIT_FAILURE;
	}
	
	return RELION_EXIT_SUCCESS;
}

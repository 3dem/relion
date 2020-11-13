#include <src/args.h>
#include <src/jaz/tomography/programs/align.h>


int main(int argc, char *argv[])
{
	AlignProgram ap(argc, argv);
	
	ap.run();
	
	return 0;
}

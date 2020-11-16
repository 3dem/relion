#include <src/args.h>
#include <src/jaz/tomography/programs/aberration_fit.h>


int main(int argc, char *argv[])
{
	AberrationFitProgram afp(argc, argv);
	
	afp.run();
	
	return 0;
}

#include <src/args.h>
#include <src/jaz/tomo_programs/bfactor_fit.h>

int main(int argc, char *argv[])
{
	BfactorFitProgram bffp(argc, argv);
	
	bffp.run();
	
	return 0;
}

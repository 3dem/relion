#include <src/args.h>
#include <src/jaz/tomography/programs/ctf_refinement.h>


int main(int argc, char *argv[])
{
	CtfRefinementProgram crp(argc, argv);
	
	crp.run();
	
	return 0;
}

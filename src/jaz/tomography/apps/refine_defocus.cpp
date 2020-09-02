#include <src/args.h>
#include <src/jaz/tomography/programs/defocus_refinement.h>


int main(int argc, char *argv[])
{
	DefocusRefinementProgram drp(argc, argv);
	
	drp.run();
	
	return 0;
}

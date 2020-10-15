#include <src/args.h>
#include <src/jaz/tomography/programs/reconstruct_tomogram.h>


int main(int argc, char *argv[])
{
	TomoBackprojectProgram tbp;
	
	tbp.readParameters(argc, argv);	
	tbp.run();
	
	return 0;
}

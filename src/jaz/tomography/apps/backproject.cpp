#include <src/args.h>
#include <src/jaz/tomography/programs/backproject.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	BackprojectProgram bp;
	
	bp.readParameters(argc, argv);	
	bp.run();
	
	return 0;
}

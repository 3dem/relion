#include <src/args.h>
#include <src/jaz/single_particle/programs/delete_blobs_2d.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	DeleteBlobs2DProgram program;
	
	program.readParameters(argc, argv);	
	program.run();
	
	return 0;
}

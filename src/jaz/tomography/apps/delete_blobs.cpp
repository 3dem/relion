#include <src/args.h>
#include <src/jaz/tomography/programs/delete_blobs.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	DeleteBlobsProgram program;
	
	program.readParameters(argc, argv);	
	program.run();
	
	return 0;
}

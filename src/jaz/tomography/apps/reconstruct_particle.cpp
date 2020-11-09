#include <src/args.h>
#include <src/jaz/tomography/programs/reconstruct_particle.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	ReconstructParticleProgram bp;
	
	bp.readParameters(argc, argv);	
	bp.run();
	
	return 0;
}

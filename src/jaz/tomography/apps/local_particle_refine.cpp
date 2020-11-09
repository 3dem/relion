#include <src/args.h>
#include <src/jaz/tomography/programs/local_particle_refine.h>


int main(int argc, char *argv[])
{
	LocalParticleRefineProgram lprp(argc, argv);

	lprp.run();

	return 0;
}

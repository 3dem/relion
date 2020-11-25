#include <src/args.h>
#include <src/jaz/tomography/programs/local_particle_refine.h>


int main(int argc, char *argv[])
{
	try
	{
		LocalParticleRefineProgram program(argc, argv);
		program.run();
	}
	catch (RelionError XE)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

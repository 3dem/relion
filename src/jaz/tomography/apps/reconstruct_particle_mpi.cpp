#include <src/args.h>
#include <mpi.h>
#include <src/jaz/tomography/programs/reconstruct_particle_mpi.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	ReconstructParticleProgramMpi prm;

	try
	{
		prm.readParameters(argc, argv);
		prm.run();
	}
	catch (RelionError XE)
	{
		MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_FAILURE);
		return RELION_EXIT_FAILURE;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return RELION_EXIT_SUCCESS;
}

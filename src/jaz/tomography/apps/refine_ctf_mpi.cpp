#include <src/args.h>
#include <mpi.h>
#include <src/jaz/tomography/programs/ctf_refinement_mpi.h>


int main(int argc, char *argv[])
{
	CtfRefinementProgramMpi crp(argc, argv);

	try
	{
		crp.run();
	}
	catch (RelionError XE)
	{
		if (crp.rank == 0)
			std::cerr << XE;
		MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_FAILURE);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return RELION_EXIT_SUCCESS;
}

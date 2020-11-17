#include <src/args.h>
#include <mpi.h>
#include <src/jaz/tomography/programs/align_mpi.h>


int main(int argc, char *argv[])
{
	AlignProgramMpi ap(argc, argv);

	try
	{
		ap.run();
	}
	catch (RelionError XE)
	{
		if (ap.rank == 0)
			std::cerr << XE;
		MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_FAILURE);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return RELION_EXIT_SUCCESS;
}

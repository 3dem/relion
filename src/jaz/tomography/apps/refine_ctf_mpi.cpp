#include <src/args.h>
#include <mpi.h>
#include <src/jaz/tomography/programs/ctf_refinement_mpi.h>


int main(int argc, char *argv[])
{
	CtfRefinementProgramMpi program(argc, argv);

	try
	{
		if (program.nodeCount < 2)
		{
			REPORT_ERROR("refine_ctf_mpi: this program needs to be run with at least two MPI processes!");
		}

		program.run();
	}
	catch (RelionError e)
	{
		std::cerr << e << std::endl;

		MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_FAILURE);

		return RELION_EXIT_FAILURE;
	}

	MPI_Barrier(MPI_COMM_WORLD);

	return RELION_EXIT_SUCCESS;
}

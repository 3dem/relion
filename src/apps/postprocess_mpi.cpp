/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
#include <src/postprocessing.h>

int main(int argc, char *argv[])
{
	Postprocessing prm;

	try
	{
		int rank, size;

		MPI_Init(&argc, &argv);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		// Handle errors
		MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

		prm.read(argc, argv);

		// Don't put any output to screen for mpi followers
		prm.verb = (rank == 0) ? 1 : 0;

		if (prm.do_locres)
			prm.run_locres(rank, size);
		else
			prm.run();
	}
	catch (RelionError XE)
	{
		//prm.usage();
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

        MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return RELION_EXIT_SUCCESS;
}

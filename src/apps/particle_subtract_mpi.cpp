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

#include <src/args.h>
#include "../particle_subtractor.h"
#include "../mpi.h"

int main(int argc, char *argv[])
{
	ParticleSubtractor prm;

	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	// Handle errors
	MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);
	MPI_Status status;

	try
	{
		prm.read(argc, argv);

		prm.initialise(rank, size);

		if (prm.fn_revert != "")
			REPORT_ERROR("You cannot use MPI for reverting subtraction.");

		prm.run();

		if (prm.do_ssnr)
		{
			MultidimArray<RFLOAT> Maux(prm.sum_S2);
			MPI_Allreduce(MULTIDIM_ARRAY(prm.sum_S2), MULTIDIM_ARRAY(Maux),
					MULTIDIM_SIZE(prm.sum_S2), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			prm.sum_S2 = Maux;
			MPI_Allreduce(MULTIDIM_ARRAY(prm.sum_N2), MULTIDIM_ARRAY(Maux),
					MULTIDIM_SIZE(prm.sum_N2), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			prm.sum_N2 = Maux;
			MPI_Allreduce(MULTIDIM_ARRAY(prm.sum_count), MULTIDIM_ARRAY(Maux),
					MULTIDIM_SIZE(prm.sum_count), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			prm.sum_count=Maux;
		}

		prm.saveStarFile(rank);

		MPI_Barrier(MPI_COMM_WORLD);

		if (rank == 0) prm.combineStarFile(rank);

		MPI_Barrier(MPI_COMM_WORLD);
	}
	catch (RelionError XE)
	{
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

        MPI_Barrier(MPI_COMM_WORLD);

	MPI_Finalize();

	return RELION_EXIT_SUCCESS;
}

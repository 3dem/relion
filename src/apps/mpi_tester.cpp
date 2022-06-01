/***************************************************************************
 *
 * Author: "Takanori Nakane"
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
#include "src/mpi.h"
#include "src/multidim_array.h"

// This simple program tests MPI communication by sending
// blocks of 1 MB to 600 MB.

int main(int argc, char *argv[])
{
	MpiNode node(argc, argv);

	const int max_mb = 600;

	MultidimArray<RFLOAT> buf(max_mb * 1024 * 1024);
	for (int i = 1; i < max_mb; i++)
	{
		printf("i = %d\n", i);
		if (node.rank == 1)
		{
			node.relion_MPI_Send(MULTIDIM_ARRAY(buf), i * 1024 * 1024, MY_MPI_DOUBLE, 0, MPITAG_PACK, MPI_COMM_WORLD);
		}
		else if (node.rank == 0)
		{
			MPI_Status status;
			node.relion_MPI_Recv(MULTIDIM_ARRAY(buf), i * 1024 * 1024, MY_MPI_DOUBLE, 1, MPITAG_PACK, MPI_COMM_WORLD, status);
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	return 0;
}

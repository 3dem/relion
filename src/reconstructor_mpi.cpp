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
#include "src/reconstructor_mpi.h"

void ReconstructorMpi::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	// First read in non-parallelisation-dependent variables
	Reconstructor::read(argc, argv);

	// Don't put any output to screen for mpi followers
	verb = (node->isLeader()) ? verb : 0;

	// Possibly also read parallelisation-dependent variables here

	if (node->size < 2)
		REPORT_ERROR("ReconstductMpi::read ERROR: this program needs to be run with at least two MPI processes!");

	// Print out MPI info
	printMpiNodesMachineNames(*node);

}

void ReconstructorMpi::run()
{


	if (fn_debug != "")
	{
		Reconstructor::readDebugArrays();
	}
	else
	{
		Reconstructor::initialise();
		Reconstructor::backproject(node->rank, node->size);

		MultidimArray<Complex> sumd(backprojector.data);
		MultidimArray<RFLOAT> sumw(backprojector.weight);
		MPI_Allreduce(MULTIDIM_ARRAY(backprojector.data), MULTIDIM_ARRAY(sumd), 2*MULTIDIM_SIZE(backprojector.data), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MULTIDIM_ARRAY(backprojector.weight), MULTIDIM_ARRAY(sumw), MULTIDIM_SIZE(backprojector.weight), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		if (node->isLeader())
		{
			backprojector.data = sumd;
			backprojector.weight = sumw;
		}

	}

	if (node->isLeader())
		reconstruct();

}

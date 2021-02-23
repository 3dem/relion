/***************************************************************************
 *
 * Author: "Jasenko Zivanov & Sjors H.W. Scheres"
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

#include "ctf_refiner_mpi.h"

void CtfRefinerMpi::read(int argc, char **argv)
{
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    CtfRefiner::read(argc, argv);

    // Don't put any output to screen for mpi followers
    verb = (node->isLeader()) ? verb : 0;

    // Possibly also read parallelisation-dependent variables here
	if (node->size < 2)
	{
		REPORT_ERROR_STR("ParticlePolisherMpi::read ERROR: this program needs to be run "
						 << "with at least two MPI processes!");
	}

    // Print out MPI info
	printMpiNodesMachineNames(*node);
}

void CtfRefinerMpi::run()
{
	// Parallel loop over micrographs

	long int total_nr_micrographs = unfinishedMdts.size();

	// Each node does part of the work
	long int my_first_micrograph, my_last_micrograph;
	divide_equally(total_nr_micrographs, node->size, node->rank, my_first_micrograph, my_last_micrograph);

	if (do_defocus_fit || do_bfac_fit || do_tilt_fit || do_aberr_fit || do_mag_fit)
    {
    	processSubsetMicrographs(my_first_micrograph, my_last_micrograph);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (node->isLeader())
    {
		finalise();
    }
}

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

#include "src/motion_fitter_mpi.h"

void MotionFitterMpi::read(int argc, char **argv)
{
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    MotionFitter::read(argc, argv);

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? verb : 0;

    // Possibly also read parallelisation-dependent variables here
	if (node->size < 2)
    {
		REPORT_ERROR("ParticlePolisherMpi::read ERROR: this program needs to be run with at least two MPI processes!");
    }

    if (node->isMaster() && (motionParamEstimator.estim2 || motionParamEstimator.estim3))
    {
        REPORT_ERROR("Parameter estimation is currently not supported in MPI mode.");
        return;
    }

    // Print out MPI info
	printMpiNodesMachineNames(*node);

}

void MotionFitterMpi::run()
{
    if (motionParamEstimator.estim2 || motionParamEstimator.estim3)
    {
        return;
    }

	// Parallel loop over micrographs

	if (mdts.size() > 0)
	{
		long int total_nr_micrographs = mdts.size();

		// Each node does part of the work
		long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
		divide_equally(total_nr_micrographs, node->size, node->rank, my_first_micrograph, my_last_micrograph);
		my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

		// The subsets will be used in openMPI parallelisation: instead of over g0->gc, they will be over smaller subsets
		processSubsetMicrographs(my_first_micrograph, my_last_micrograph);
	}

    MPI_Barrier(MPI_COMM_WORLD);

    if (doCombineFrames)
	{
		initialiseCombineFrames();
		if (mdts.size() > 0)
		{
			long int total_nr_micrographs = mdts.size();

			// Each node does part of the work
			long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
			divide_equally(total_nr_micrographs, node->size, node->rank, my_first_micrograph, my_last_micrograph);
			my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

			combineFramesSubsetMicrographs(my_first_micrograph, my_last_micrograph);
		}
	}


    MPI_Barrier(MPI_COMM_WORLD);

	if (node->isMaster())
		combineEPSAndSTARfiles();
}

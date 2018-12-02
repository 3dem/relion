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

#include "motion_refiner_mpi.h"

void MotionRefinerMpi::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	// First read in non-parallelisation-dependent variables
	MotionRefiner::read(argc, argv);

	// Don't put any output to screen for mpi slaves
	verb = (node->isMaster()) ? verb : 0;

	// Possibly also read parallelisation-dependent variables here
	if (node->size < 2)
	{
        REPORT_ERROR("ERROR: this program needs to be run with at least two MPI processes!");
	}

    if (node->isMaster() && (motionParamEstimator.anythingToDo()))
	{
        REPORT_ERROR("Parameter estimation is not supported in MPI mode.");
		return;
	}

	// Print out MPI info
	printMpiNodesMachineNames(*node);
}

void MotionRefinerMpi::run()
{
    if (estimateParams)
	{
        REPORT_ERROR("Parameter estimation is not supported in MPI mode.");
		return;
	}

	// Parallel loop over micrographs

    if (estimateMotion)
	{
        long int total_nr_micrographs = motionMdts.size();

		// Each node does part of the work
		long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
		divide_equally(total_nr_micrographs, node->size, node->rank, my_first_micrograph, my_last_micrograph);
		my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

        motionEstimator.process(motionMdts, my_first_micrograph, my_last_micrograph);
	}

	MPI_Barrier(MPI_COMM_WORLD);

    if (recombineFrames)
    {
        long int total_nr_micrographs = recombMdts.size();

        // Each node does part of the work
        long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
        divide_equally(total_nr_micrographs, node->size, node->rank,
                       my_first_micrograph, my_last_micrograph);
        my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;
		
		double k_out_A = reference.pixToAng(reference.k_out);
		
        frameRecombiner.init(
            allMdts, 
			verb, reference.s, fc, k_out_A, reference.angpix,
			nr_omp_threads, outPath, debug,
            &obsModel, &micrographHandler);

        frameRecombiner.process(recombMdts, my_first_micrograph, my_last_micrograph);
	}

	MPI_Barrier(MPI_COMM_WORLD);

    if (generateStar && node->isMaster())
	{
		combineEPSAndSTARfiles();
	}
}

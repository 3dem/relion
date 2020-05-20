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

void MotionRefinerMpi::runWithFccUpdate()
{
	if (estimateParams)
	{
		REPORT_ERROR("Parameter estimation is not supported in MPI mode.");
		return;
	}
	
	const int lastMgForFCC = lastMicrographForFCC();
	
	// Parallel loop over micrographs
	
	if (estimateMotion)
	{
		// Each node does part of the work
		long int my_first_micrograph_FCC, my_last_micrograph_FCC;
		
		divide_equally(lastMgForFCC, node->size, node->rank, 
					   my_first_micrograph_FCC, my_last_micrograph_FCC);
		
		motionEstimator.process(motionMdts, my_first_micrograph_FCC, my_last_micrograph_FCC, true);
	}
}

void MotionRefinerMpi::runWithRecombination()
{
	double k_out_A = reference.pixToAng(reference.k_out);
	
	const int lastMgForFCC = lastMicrographForFCC();
	const int firstMgWithoutFCC = lastMgForFCC + 1;
	const long int total_nr_micrographs = recombMdts.size();
	const bool anyWithoutFCC = firstMgWithoutFCC < motionMdts.size();
	
	
	// micrographs [0 ... lastMgForFCC] have already been aligned:
	long int my_first_micrograph_FCC, my_last_micrograph_FCC;
	
	divide_equally(
		lastMgForFCC, node->size, node->rank, 
		my_first_micrograph_FCC, my_last_micrograph_FCC);
	
	
	// [lastMgForFCC+1 ... total_nr_micrographs] have not:
	long int my_first_micrograph_recomb, my_last_micrograph_recomb;
	
	divide_equally(
		total_nr_micrographs - firstMgWithoutFCC, node->size, node->rank,
		my_first_micrograph_recomb, my_last_micrograph_recomb);
	
	my_first_micrograph_recomb += firstMgWithoutFCC;
	my_last_micrograph_recomb += firstMgWithoutFCC;
	
	
	if (recombineFrames)
	{
		frameRecombiner.init(
			allMdts, 
			verb, reference.s, fc, k_out_A, reference.angpix,
			nr_omp_threads, outPath, debug,
			&reference, &obsModel, &micrographHandler);
		
		frameRecombiner.process(recombMdts, my_first_micrograph_FCC, my_last_micrograph_FCC);
		
		
		if (anyWithoutFCC)
		{
			motionEstimator.setVerbosity(0);
			frameRecombiner.setVerbosity(0);
			
			const int left = motionMdts.size() - lastMgForFCC - 1;
			const int barstep = XMIPP_MAX(1, left/ 60);
			
			if (verb > 0)
			{
				std::cout << " + Aligning and combining frames for micrographs ... " << std::endl;
				init_progress_bar(left);
			}
			
			for (int m = my_first_micrograph_recomb; m <= my_last_micrograph_recomb; m++)
			{
				motionEstimator.process(motionMdts, m, m, false);
				frameRecombiner.process(recombMdts, m, m);
				
				const int nr_done = m - lastMgForFCC;
				
				if (verb > 0 && nr_done % barstep == 0)
				{
					progress_bar(nr_done);
				}
			}
			
			if (verb > 0)
			{
				progress_bar(left);
			}
		}
	}
	else if (anyWithoutFCC)
	{
		motionEstimator.process(
				motionMdts, my_first_micrograph_recomb, my_last_micrograph_recomb, false);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (generateStar && node->isMaster())
	{
		combineEPSAndSTARfiles();
	}
}

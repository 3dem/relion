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

#include "src/amyloid_finder_mpi.h"

void AmyloidFinderMpi::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	if (node->isLeader())
		PRINT_VERSION_INFO();

	// First read in non-parallelisation-dependent variables
	AmyloidFinder::read(argc, argv);

	// Don't put any output to screen for mpi followers
	if (!node->isLeader())
		verb = 0;

	// Possibly also read parallelisation-dependent variables here

    if (do_plot) REPORT_ERROR("ERROR: cannot use --plot in parallel execution!");
    
	// Print out MPI info
	printMpiNodesMachineNames(*node);
}

#if defined _CUDA_ENABLED
void AmyloidFinderMpi::deviceInitialise()
{
	int devCount;
	accGPUGetDeviceCount(&devCount);

	std::vector < std::vector < std::string > > allThreadIDs;
	untangleDeviceIDs(gpu_ids, allThreadIDs);

	// Sequential initialisation of GPUs on all ranks
	if (!std::isdigit(*gpu_ids.begin()))
		device_id = node->rank%devCount;
	else
		device_id = textToInteger((allThreadIDs[node->rank][0]).c_str());

	for (int follower = 0; follower < node->size; follower++)
	{
		if (follower == node->rank)
		{
			std::cout << " + Using GPU device: " << device_id << " on MPI node: " << node->rank << std::endl;
			std::cout.flush();
		}
		node->barrierWait();
	}
}
#endif

void AmyloidFinderMpi::run()
{
	// Each node does part of the work
	if (todo_micrographs_fom.size() > 0)
    {
        long int my_first_fom, my_last_fom;
        divide_equally(todo_micrographs_fom.size(), node->size, node->rank, my_first_fom, my_last_fom);
        runFOMBatch(my_first_fom, my_last_fom);
    }

    if (todo_micrographs_tracing.size() > 0)
    {
        long int my_first_tracing, my_last_tracing;
        divide_equally(todo_micrographs_tracing.size(), node->size, node->rank, my_first_tracing, my_last_tracing);
        runTracingBatch(my_first_tracing, my_last_tracing, node->rank);
    }

}

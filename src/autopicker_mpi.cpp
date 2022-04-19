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

#include "src/autopicker_mpi.h"

void AutoPickerMpi::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	if (node->isLeader())
		PRINT_VERSION_INFO();

	// First read in non-parallelisation-dependent variables
	AutoPicker::read(argc, argv);

	// Don't put any output to screen for mpi followers
	if (!node->isLeader())
		verb = 0;

	if (do_write_fom_maps && node->isLeader())
		std::cerr << "WARNING : --write_fom_maps is very heavy on disc I/O and is not advised in parallel execution. If possible, using --shrink 0 and lowpass makes I/O less significant." << std::endl;

	// Possibly also read parallelisation-dependent variables here

	// Print out MPI info
	printMpiNodesMachineNames(*node);
}

#ifdef _CUDA_ENABLED
void AutoPickerMpi::deviceInitialise()
{
	int devCount;
	cudaGetDeviceCount(&devCount);

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

void AutoPickerMpi::run()
{
	// Each node does part of the work
	long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
	divide_equally(fn_micrographs.size(), node->size, node->rank, my_first_micrograph, my_last_micrograph);
	my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

	int barstep;
	if (verb > 0)
	{
		std::cout << " Autopicking ..." << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs / 60);
	}

	FileName fn_olddir="";
	for (long int imic = my_first_micrograph; imic <= my_last_micrograph; imic++)
	{
		// Abort through the pipeline_control system
		if (pipeline_control_check_abort_job())
			MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_ABORTED);

		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		// Check new-style outputdirectory exists and make it if not!
		FileName fn_dir = getOutputRootName(fn_micrographs[imic]);
		fn_dir = fn_dir.beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			mktree(fn_dir);
			fn_olddir = fn_dir;
		}

		if (do_topaz_extract)
			autoPickTopazOneMicrograph(fn_micrographs[imic], node->rank);
		else if (do_LoG)
			autoPickLoGOneMicrograph(fn_micrographs[imic], imic);
		else
			autoPickOneMicrograph(fn_micrographs[imic], imic);
	}

	if (verb > 0)
		progress_bar(my_nr_micrographs);
}

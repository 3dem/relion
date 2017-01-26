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

    // First read in non-parallelisation-dependent variables
    AutoPicker::read(argc, argv);

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? 1 : 0;

    if (do_write_fom_maps && node->isMaster())
    	std::cerr << "WARNING : --write_fom_maps is very heavy on disc I/O and is not advised in parallel execution. If possible, using --shrink 0 and lowpass makes I/O less significant." << std::endl;

    // Possibly also read parallelisation-dependent variables here

    // Print out MPI info
	printMpiNodesMachineNames(*node);


}

#ifdef CUDA
int AutoPickerMpi::deviceInitialise()
{
	int devCount;
	cudaGetDeviceCount(&devCount);

	std::vector < std::vector < std::string > > allThreadIDs;
	untangleDeviceIDs(gpu_ids, allThreadIDs);

	// Sequential initialisation of GPUs on all ranks
	int dev_id;
	if (!std::isdigit(*gpu_ids.begin()))
		dev_id = node->rank%devCount;
	else
		dev_id = textToInteger((allThreadIDs[node->rank][0]).c_str());

    for (int slave = 0; slave < node->size; slave++)
    {
    	if (slave == node->rank)
    	{
    		std::cout << " + Using GPU device: " << dev_id << " on MPI node: " << node->rank << std::endl;
    		std::cout.flush();
    	}
    	node->barrierWait();
    }

	return(dev_id);

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
    	if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		// Check new-style outputdirectory exists and make it if not!
		FileName fn_dir = getOutputRootName(fn_micrographs[imic]);
		fn_dir = fn_dir.beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			int res = system(("mkdir -p " + fn_dir).c_str());
			fn_olddir = fn_dir;
		}

    	autoPickOneMicrograph(fn_micrographs[imic], imic);
	}
	if (verb > 0)
		progress_bar(my_nr_micrographs);


}

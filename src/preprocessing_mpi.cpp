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
#include "src/preprocessing_mpi.h"

void PreprocessingMpi::read(int argc, char **argv)
{
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    Preprocessing::read(argc, argv, node->rank);

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? 1 : 0;

    // Possibly also read parallelisation-dependent variables here

    // Print out MPI info
	printMpiNodesMachineNames(*node);


}


void PreprocessingMpi::runExtractParticles()
{

	// Each node does part of the work
	long int my_first_coord, my_last_coord, my_nr_coords;
	divide_equally(fn_coords.size(), node->size, node->rank, my_first_coord, my_last_coord);
	my_nr_coords = my_last_coord - my_first_coord + 1;

	int barstep;
	if (verb > 0)
	{
		std::cout << " Extracting particles from the micrographs ..." << std::endl;
		init_progress_bar(my_nr_coords);
		barstep = XMIPP_MAX(1, my_nr_coords / 60);

		// Make a Particles directory
		int res = system("mkdir -p Particles");
	}

	FileName fn_olddir = "";
	for (long int ipos = my_first_coord; ipos <= my_last_coord; ipos++)
    {
		FileName fn_dir = "Particles/" + fn_coords[ipos].beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			int res = system(("mkdir -p " + fn_dir).c_str());
			fn_olddir = fn_dir;
		}

		if (verb > 0 && ipos % barstep == 0)
			progress_bar(ipos);

    	extractParticlesFromFieldOfView(fn_coords[ipos]);
	}

	if (verb > 0)
		progress_bar(my_nr_coords);

}



void PreprocessingMpi::run()
{

	// Extract and operate on particles in parallel
	if (do_extract)
	{
		runExtractParticles();

		// Wait until all nodes have finished to make final star file
		MPI_Barrier(MPI_COMM_WORLD);
	}

	if (do_join_starfile && node->isMaster())
		Preprocessing::joinAllStarFiles();

	// The following has not been parallelised....
	if (fn_operate_in != "" && node->isMaster())
		Preprocessing::runOperateOnInputFile(fn_operate_in);

	if (verb > 0)
		std::cout << " Done!" <<std::endl;

}

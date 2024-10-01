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
#include "src/ctffind_runner_mpi.h"

void CtffindRunnerMpi::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	// First read in non-parallelisation-dependent variables
	CtffindRunner::read(argc, argv);

	// Don't put any output to screen for mpi followers
	verb = (node->isLeader()) ? 1 : 0;

	// Possibly also read parallelisation-dependent variables here

	// Print out MPI info
	printMpiNodesMachineNames(*node);
}

void CtffindRunnerMpi::run()
{
	if (!do_only_join_results)
	{
		// Each node does part of the work
		long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
		divide_equally(fn_micrographs.size(), node->size, node->rank, my_first_micrograph, my_last_micrograph);
		my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

		int barstep;
		if (verb > 0)
		{
            std::cout << " Estimating CTF parameters using Niko Grigorieff's CTFFIND ..." << std::endl;
			init_progress_bar(my_nr_micrographs);
			barstep = XMIPP_MAX(1, my_nr_micrographs / 60);
		}

		std::vector<std::string> allmicnames;
		for (long int imic = my_first_micrograph; imic <= my_last_micrograph; imic++)
		{

			// Abort through the pipeline_control system
			if (pipeline_control_check_abort_job())
				MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_ABORTED);

			// Get angpix and voltage from the optics groups:
			obsModel.opticsMdt.getValue(EMDL_CTF_CS, Cs, optics_group_micrographs[imic]-1);
			obsModel.opticsMdt.getValue(EMDL_CTF_VOLTAGE, Voltage, optics_group_micrographs[imic]-1);
			obsModel.opticsMdt.getValue(EMDL_CTF_Q0, AmplitudeConstrast, optics_group_micrographs[imic]-1);
            EMDLabel mylabel = (is_tomo) ? EMDL_TOMO_TILT_SERIES_PIXEL_SIZE : EMDL_MICROGRAPH_PIXEL_SIZE;
			obsModel.opticsMdt.getValue(mylabel, angpix, optics_group_micrographs[imic]-1);

			if (is_ctffind4)
			{
				executeCtffind4(imic);
			}
			else if (is_ctffind5)
			{
				executeCtffind5(imic);
			}
			else
			{
				executeCtffind3(imic);
			}

			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);

		}

		if (verb > 0)
			progress_bar(my_nr_micrographs);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Only the leader writes the joined result file
	if (node->isLeader())
	{
		joinCtffindResults();
	}
}

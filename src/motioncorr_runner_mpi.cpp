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
#include "src/motioncorr_runner_mpi.h"

void MotioncorrRunnerMpi::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	// First read in non-parallelisation-dependent variables
	MotioncorrRunner::read(argc, argv);

	// Don't put any output to screen for mpi followers
	verb = (node->isLeader()) ? 1 : 0;

	// Print out MPI info
	printMpiNodesMachineNames(*node);
}

void MotioncorrRunnerMpi::run()
{
	prepareGainReference(node->isLeader());
	MPI_Barrier(MPI_COMM_WORLD); // wait for the leader to write the gain reference

	// Each node does part of the work
	long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
	divide_equally(fn_micrographs.size(), node->size, node->rank, my_first_micrograph, my_last_micrograph);
	my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

	int barstep;
	if (verb > 0)
	{
		if (do_own)
			 std::cout << " Correcting beam-induced motions using our own implementation ..." << std::endl;
		else if (do_motioncor2)
			std::cout << " Correcting beam-induced motions using Shawn Zheng's MOTIONCOR2 ..." << std::endl;
		else
			REPORT_ERROR("Bug: by now it should be clear whether to use MotionCor2 or Unblur...");

		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs / 60);
	}

	for (long int imic = my_first_micrograph; imic <= my_last_micrograph; imic++)
	{
		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		// Abort through the pipeline_control system
		if (pipeline_control_check_abort_job())
			MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_ABORTED);

		Micrograph mic(fn_micrographs[imic], fn_gain_reference, bin_factor, eer_upsampling, eer_grouping);

		// Get angpix and voltage from the optics groups:
		obsModel.opticsMdt.getValue(EMDL_CTF_VOLTAGE, voltage, optics_group_micrographs[imic]-1);
		obsModel.opticsMdt.getValue(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, angpix, optics_group_micrographs[imic]-1);

		bool result;
		if (do_own)
			result = executeOwnMotionCorrection(mic);
		else if (do_motioncor2)
			result = executeMotioncor2(mic, node->rank);
		else
			REPORT_ERROR("Bug: by now it should be clear whether to use MotionCor2 or Unblur...");

		if (result) {
			saveModel(mic);
			plotShifts(fn_micrographs[imic], mic);
		}
	}
	if (verb > 0)
		progress_bar(my_nr_micrographs);

	MPI_Barrier(MPI_COMM_WORLD);

	// Only the leader writes the joined result file
	if (node->isLeader())
		generateLogFilePDFAndWriteStarFiles();

}

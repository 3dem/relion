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
#include "src/motioncorr_own_devolved.h"

void MotioncorrOwnDevolved::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	parser.setCommandLine(argc, argv);

	movie_path = parser.getOption("--in_movie", "Path to input movie");
	micrograph_path = parser.getOption("--out_mic", "Output micrograph path");
	fn_gain_reference = parser.getOption("--gainref","Location of MRC file with the gain reference to be applied","");
	bin_factor =  textToFloat(parser.getOption("--bin_factor", "Binning factor (can be non-integer)", "1"));
	n_threads = textToInteger(parser.getOption("--j", "Number of threads per movie (= process)", "1"));
	ps_size = textToInteger(parser.getOption("--ps_size", "Output size of power spectrum", "512"));
	if (ps_size % 2 != 0) REPORT_ERROR("--ps_size must be an even number.");
	first_frame_sum =  textToInteger(parser.getOption("--first_frame_sum", "First movie frame used in output sum (start at 1)", "1"));
	if (first_frame_sum < 1) first_frame_sum = 1;
	last_frame_sum =  textToInteger(parser.getOption("--last_frame_sum", "Last movie frame used in output sum (0 or negative: use all)", "-1"));
	eer_grouping = textToInteger(parser.getOption("--eer_grouping", "EER grouping", "40"));
	eer_upsampling = textToInteger(parser.getOption("--eer_upsampling", "EER upsampling (1 = 4K or 2 = 8K)", "1"));

	// Don't put any output to screen for mpi followers
	verb = (node->isLeader()) ? 1 : 0;

	// Print out MPI info
	printMpiNodesMachineNames(*node);
}

void MotioncorrOwnDevolved::run()
{
	prepareGainReference(node->isLeader());
	MPI_Barrier(MPI_COMM_WORLD); // wait for the leader to write the gain reference

	Micrograph mic(movie_path, fn_gain_reference, bin_factor, eer_upsampling, eer_grouping);
	bool result;
	result = executeOwnMotionCorrection(mic);
	if (result) saveModel(mic);

	MPI_Barrier(MPI_COMM_WORLD);

}

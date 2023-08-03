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

	MotioncorrRunner::read(argc, argv);

	// Don't put any output to screen for mpi followers
	verb = (node->isLeader()) ? 1 : 0;

	// Print out MPI info
	printMpiNodesMachineNames(*node);
}

void MotioncorrOwnDevolved::addClArgs()
{
	int path_section =  parser.addSection("In/out paths options");
	movie_path = parser.getOption("--in_movie", "Path to input movie");
	micrograph_path = parser.getOption("--out_mic", "Output micrograph path");
	MotioncorrRunner::addClArgs();
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

FileName MotioncorrOwnDevolved::getOutputFileNames(FileName fn_mic)
{
	return micrograph_path;
}

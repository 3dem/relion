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

#include "src/ctf_refinement_mpi.h"

void CtfRefinerMpi::read(int argc, char **argv)
{
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    CtfRefiner::read(argc, argv);

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? verb : 0;

    // Possibly also read parallelisation-dependent variables here
	if (node->size < 2)
		REPORT_ERROR("ParticlePolisherMpi::read ERROR: this program needs to be run with at least two MPI processes!");

    // Print out MPI info
	printMpiNodesMachineNames(*node);

}

void CtfRefinerMpi::run()
{
	// Parallel loop over micrographs


	long int total_nr_micrographs = gc - g0 + 1;

	// Each node does part of the work
	long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
	divide_equally(total_nr_micrographs, node->size, node->rank, my_first_micrograph, my_last_micrograph);
	my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

    // If there were multiple groups of micrographs, we could consider introducing a loop over those here...
    //for (int igroup = 0; igroup < nr_micrographs_groups; igroup++)
    //{

    Image<Complex> xyAccSum;
    Image<RFLOAT> wAccSum;
	if (do_tilt_fit && !precomputed)
    {
    	xyAccSum().initZeros(sh,s);
    	wAccSum().initZeros(sh,s);
    }
	else
	{
		xyAccSum = lastXY;
		wAccSum = lastW;
	}

    if (do_defocus_fit || (do_tilt_fit && !precomputed) )
    {
    	// The subsets will be used in openMPI parallelisation: instead of over g0->gc, they will be over smaller subsets
    	processSubsetMicrographs(my_first_micrograph, my_last_micrograph, xyAccSum, wAccSum);
    }

    MPI_Barrier(MPI_COMM_WORLD);

	if (do_tilt_fit)
	{
		Image<Complex> xyAccSum2;
		Image<RFLOAT> wAccSum2;
		xyAccSum2().initZeros(xyAccSum());
		wAccSum2().initZeros(wAccSum());
		MPI_Allreduce(MULTIDIM_ARRAY(xyAccSum()), MULTIDIM_ARRAY(xyAccSum2()), 2*MULTIDIM_SIZE(xyAccSum()), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		MPI_Allreduce(MULTIDIM_ARRAY(wAccSum()), MULTIDIM_ARRAY(wAccSum2()), MULTIDIM_SIZE(wAccSum()), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		fitBeamTiltFromSumsAllMicrographs(xyAccSum2, wAccSum2);
	}

	//} // end loop over igroup

	// TODO: design mechanism to set the defocus parameters of all individual particles from different MPI ranks...
	// TODO: will need to pass the values through MPI_Send....
    MetaDataTable mdtAll;
    mdtAll.reserve(mdt0.numberOfObjects());
	for (long g = g0; g <= gc; g++)
		mdtAll.append(mdts[g]);

	mdtAll.write(outPath + "particles.star");


}

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
#include "src/particle_sorter_mpi.h"

void ParticleSorterMpi::read(int argc, char **argv)
{
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    ParticleSorter::read(argc, argv);

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? 1 : 0;

    // Possibly also read parallelisation-dependent variables here

    // Print out MPI info
	printMpiNodesMachineNames(*node);


}
void ParticleSorterMpi::run()
{

	int total_nr_images = MDin.numberOfObjects();
	features.resize(total_nr_images, NR_FEATURES);

	// Each node does part of the work
	long int my_first_image, my_last_image, my_nr_images;
	divide_equally(total_nr_images, node->size, node->rank, my_first_image, my_last_image);
	my_nr_images = my_last_image - my_first_image + 1;

	int barstep;
	if (verb > 0)
	{
		std::cout << "Calculating sorting features for all input particles..." << std::endl;
		init_progress_bar(my_nr_images);
		barstep = XMIPP_MAX(1, my_nr_images/ 60);
	}

	long int ipart = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
	{

		if (ipart >= my_first_image && ipart <= my_last_image)
		{
			if (verb > 0 && ipart % barstep == 0)
				progress_bar(ipart);

			calculateFeaturesOneParticle(ipart);

		}
		ipart++;
	}

	if (verb > 0)
		progress_bar(my_nr_images);

	// Combine results from all nodes
	MultidimArray<RFLOAT> allnodes_features;
	allnodes_features.resize(features);
	MPI_Allreduce(MULTIDIM_ARRAY(features), MULTIDIM_ARRAY(allnodes_features), MULTIDIM_SIZE(features), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	features = allnodes_features;

	// Only the master writes out files
	if (verb > 0)
	{
		normaliseFeatures();

		writeFeatures();
	}

}

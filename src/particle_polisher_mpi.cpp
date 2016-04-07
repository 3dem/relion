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
#include "src/particle_polisher_mpi.h"

void ParticlePolisherMpi::read(int argc, char **argv)
{
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    ParticlePolisher::read(argc, argv);

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? verb : 0;

    // Possibly also read parallelisation-dependent variables here

	if (node->size < 2)
		REPORT_ERROR("ParticlePolisherMpi::read ERROR: this program needs to be run with at least two MPI processes!");

    // Print out MPI info
	printMpiNodesMachineNames(*node);


}

void ParticlePolisherMpi::generateMicrographList()
{
	// Only the master makes the list of micrographs and the STAR files of movie-particles in each of them
	if (node->isMaster())
		ParticlePolisher::generateMicrographList();

	MPI_Barrier(MPI_COMM_WORLD);

	// All slaves read the list in
	if (!node->isMaster())
	{
		MetaDataTable MDmics;
		MDmics.read(fn_out+"micrograph_list.star");

		fn_mics.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmics)
		{
			FileName fn_mic;
			MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			fn_mics.push_back(fn_mic);
		}
	}

}

// Fit the beam-induced translations for all average micrographs
void ParticlePolisherMpi::fitMovementsAllMicrographs()
{

	long int total_nr_micrographs = fn_mics.size();

	// Each node does part of the work
	long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
	divide_equally(total_nr_micrographs, node->size, node->rank, my_first_micrograph, my_last_micrograph);
	my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

	// Loop over all average micrographs
	int barstep;
	if (verb > 0)
	{
		std::cout << " + Fitting straight paths for beam-induced movements in all micrographs ... " << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
	}

	for (long int i = my_first_micrograph; i <= my_last_micrograph; i++)
	{
    	if (verb > 0 && i % barstep == 0)
			progress_bar(i);

		fitMovementsOneMicrograph(i);
	}

	// Wait until all micrographs have been done
	MPI_Barrier(MPI_COMM_WORLD);

	if (verb > 0)
	{
		progress_bar(my_nr_micrographs);
	}


}

void ParticlePolisherMpi::calculateAllSingleFrameReconstructionsAndBfactors()
{

	FileName fn_star = fn_out + "bfactors.star";
	if (only_do_unfinished && readStarFileBfactors(fn_star))
	{
		if (verb > 0)
			std::cout << " + " << fn_star << " already exists: skipping calculation average of per-frame B-factors." <<std::endl;
		return;
	}

	RFLOAT bfactor, offset, corr_coeff;

	// Loop over all frames (two halves for each frame!) to be included in the reconstruction
	// Each node does part of the work
	long int my_first_frame, my_last_frame, my_nr_frames;
	divide_equally(2*movie_frame_numbers.size(), node->size, node->rank, my_first_frame, my_last_frame);
	my_nr_frames = my_last_frame - my_first_frame + 1;

	if (verb > 0)
	{
		std::cout << " + Calculating per-frame reconstructions ... " << std::endl;
		init_progress_bar(my_nr_frames);
	}

	for (long int i = my_first_frame; i <= my_last_frame; i++)
	{

		int iframe = (i >= movie_frame_numbers.size()) ? i - movie_frame_numbers.size() : i;
		int ihalf = (i >= movie_frame_numbers.size()) ? 2 : 1;

		calculateSingleFrameReconstruction(iframe, ihalf);

		if (verb > 0)
			progress_bar(i - my_first_frame + 1);
	}

	if (verb > 0)
	{
		progress_bar(my_nr_frames);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	// Also calculate the average of all single-frames for both halves
    if (node->rank == 0)
    	calculateAverageAllSingleFrameReconstructions(1);
    else if (node->rank == 1)
    	calculateAverageAllSingleFrameReconstructions(2);

	// Wait until all reconstructions have been done, and calculate the B-factors per-frame
	MPI_Barrier(MPI_COMM_WORLD);

	calculateBfactorSingleFrameReconstruction(-1, bfactor, offset, corr_coeff); // FSC between the two averages, also reads mask

	MPI_Barrier(MPI_COMM_WORLD);

	// Loop over all frames (two halves for each frame!) to be included in the reconstruction
	// Each node does part of the work
	divide_equally(movie_frame_numbers.size(), node->size, node->rank, my_first_frame, my_last_frame);
	my_nr_frames = my_last_frame - my_first_frame + 1;

	if (verb > 0)
	{
		std::cout << " + Calculating per-frame B-factors ... " << std::endl;
		init_progress_bar(my_nr_frames);
	}

	for (long int iframe = my_first_frame; iframe <= my_last_frame; iframe++)
	{

		calculateBfactorSingleFrameReconstruction(iframe, bfactor, offset, corr_coeff);
		DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 0) = bfactor;
       	DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 1) = offset;
       	DIRECT_A1D_ELEM(perframe_bfactors, iframe * 3 + 2) = corr_coeff;

    	if (verb > 0)
    		progress_bar(iframe - my_first_frame + 1);
	}

	// Combine results from all nodes
	MultidimArray<RFLOAT> allnodes_perframe_bfactors;
	allnodes_perframe_bfactors.resize(perframe_bfactors);
	MPI_Allreduce(MULTIDIM_ARRAY(perframe_bfactors), MULTIDIM_ARRAY(allnodes_perframe_bfactors), MULTIDIM_SIZE(perframe_bfactors), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	perframe_bfactors = allnodes_perframe_bfactors;

	if (verb > 0)
	{
		progress_bar(my_nr_frames);
		writeStarFileBfactors(fn_star);

	    // Also write a STAR file with the relative contributions of each frame to all frequencies
	    fn_star = fn_out + "relweights.star";
	    writeStarFileRelativeWeights(fn_star);
	}


}

void ParticlePolisherMpi::polishParticlesAllMicrographs()
{

	if (only_do_unfinished && exists(fn_out + "shiny.star"))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_out << "shiny.star already exists: skipping polishing of the particles." << std::endl;
		return;
	}

	long int total_nr_micrographs = fn_mics.size();

	// Each node does part of the work
	long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
	divide_equally(total_nr_micrographs, node->size, node->rank, my_first_micrograph, my_last_micrograph);
	my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

	// Loop over all average micrographs
	int barstep;
	if (verb > 0)
	{
		std::cout << " + Write out polished particles for all micrographs ... " << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
	}

    for (long int i = my_first_micrograph; i <= my_last_micrograph; i++)
	{
    	if (verb > 0 && i % barstep == 0)
			progress_bar(i);

    	polishParticlesOneMicrograph(i);
	}

   	if (verb > 0)
   		progress_bar(my_nr_micrographs);

    if (node->isMaster())
    	writeStarFilePolishedParticles();

    MPI_Barrier(MPI_COMM_WORLD);

}

void ParticlePolisherMpi::reconstructShinyParticlesAndFscWeight(int ipass)
{
	if (verb > 0)
		std::cout << "+ Reconstructing two halves of shiny particles ..." << std::endl;

	// Re-read the shiny particles' metadatatable (ignore original particle names here...)
	Experiment exp_model;
	exp_model.read(fn_out + "shiny.star", true);

	 // Do the reconstructions for both halves
	if (node->rank == 0)
		reconstructShinyParticlesOneHalf(1, exp_model);
	else if (node->rank == 1)
		reconstructShinyParticlesOneHalf(2, exp_model);

	// Wait until both reconstructions have been done
	MPI_Barrier(MPI_COMM_WORLD);

	// Only the master performs the FSC-weighting
	FileName fn_post = (ipass == 1) ? "shiny_post" : "shiny_post2";
	if (node->rank == 0)
	{

		if (only_do_unfinished && exists(fn_out + fn_post + "_masked.mrc")
						       && exists(fn_out + fn_post + ".star") )
		{
			if (verb > 0)
				std::cout << std::endl << " + " << fn_out << fn_post << "_masked.mrc already exists: re-reading map into memory." << std::endl;

			if (verb > 0)
				std::cout << std::endl << " + " << fn_out << fn_post << ".star already exists: re-reading resolution from it." << std::endl;

			MetaDataTable MD;
			MD.read(fn_out + fn_post + ".star", "general");
			MD.getValue(EMDL_POSTPROCESS_FINAL_RESOLUTION, maxres_model);
		}
		else
		{
			// Re-read the two halves to calculate FSCs
			Postprocessing prm;

			prm.clear();
			prm.fn_in = fn_out + "shiny";
			prm.fn_out = fn_out + fn_post;
			prm.angpix = angpix;
			prm.do_auto_mask = false;
			prm.fn_mask = fn_mask;
			prm.do_auto_bfac = false;
			prm.do_fsc_weighting = true;
			prm.verb=0;
			prm.run();

			maxres_model = prm.global_resol;
		}

		std::cout << " Resolution of reconstructions from shiny particles: " << maxres_model << std::endl;
		std::cout << " But you probably want to re-run at least a 3D auto-refinement with the shiny particles." << std::endl;

	}

}

void ParticlePolisherMpi::run()
{
	// Fit straight lines through all beam-induced translations
	if (fitting_mode != NO_FIT)
		fitMovementsAllMicrographs();

	// Perform single-frame reconstructions to estimate dose-dependent B-factors
	if (do_weighting)
		calculateAllSingleFrameReconstructionsAndBfactors();

	// Make a logfile in pdf format
	if (node->isMaster())
		generateLogFilePDF();

	// Write out the polished particles
	polishParticlesAllMicrographs();

	// Now reconstruct with all polished particles: two independent halves, FSC-weighting of the sum of the two...
	reconstructShinyParticlesAndFscWeight(1);

	if (verb > 0)
		std::cout << " done!" << std::endl;

}

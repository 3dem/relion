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
#include "src/motioncorr_own_devolved_mpi.h"

void MotioncorrOwnDevolvedMpi::addClArgs()
{
	int path_section =  parser.addSection("In/out paths options");
	motion_correction_star_path = parser.getOption("--mc_star", "Path to star file containing motion correction model information from a previous run", "");
	MotioncorrRunner::addClArgs();
}

void MotioncorrOwnDevolvedMpi::run()
{
    prepareGainReference(node->isLeader());

    bool fromStarFile = false;

    if (motion_correction_star_path != "") 
    {
        MetaDataTable MDin;
        MDin.read(motion_correction_star_path, "micrographs");
        
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
        {
            FileName fn_stars;
            MDin.getValue(EMDL_MICROGRAPH_METADATA_NAME, fn_stars);

            fn_stars_all.push_back(fn_stars);
        }
    }
    
    
	MPI_Barrier(MPI_COMM_WORLD); // wait for the leader to write the gain reference

	// Each node does part of the work
	long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
	divide_equally(fn_micrographs.size(), node->size, node->rank, my_first_micrograph, my_last_micrograph);
	my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

    for (long int imic = my_first_micrograph; imic <= my_last_micrograph; imic++)
	{
        // Abort through the pipeline_control system
		if (pipeline_control_check_abort_job())
			MPI_Abort(MPI_COMM_WORLD, RELION_EXIT_ABORTED);

        Micrograph mic(fn_micrographs[imic], fn_gain_reference, bin_factor, eer_upsampling, eer_grouping);
        if (motion_correction_star_path != "") 
        {
            FileName shifts_star_path = fn_stars_all[imic];
    		Micrograph mic2(shifts_star_path, "", 1);

            mic = mic2;
	    	fromStarFile = true;
	    }
        
        if (pre_exposure_micrographs.size() == imic)
    		mic.pre_exposure = pre_exposure + pre_exposure_micrographs[imic];
		else
    		mic.pre_exposure = pre_exposure;

        // Get angpix and voltage from the optics groups:
		obsModel.opticsMdt.getValue(EMDL_CTF_VOLTAGE, voltage, optics_group_micrographs[imic]-1);
		obsModel.opticsMdt.getValue(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, angpix, optics_group_micrographs[imic]-1);

        bool result;
        result = executeOwnMotionCorrection(mic, fromStarFile);

        if (result && !fromStarFile)
			saveModel(mic);
	}

    MPI_Barrier(MPI_COMM_WORLD);

    // Only the leader writes the joined result file
    if (node->isLeader() && !fromStarFile)
        generateLogFilePDFAndWriteStarFiles();
}

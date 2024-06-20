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
#ifndef RELION_ALIGN_TILTSERIES_RUNNER_H
#define RELION_ALIGN_TILTSERIES_RUNNER_H


#include  <glob.h>
#include  <vector>
#include  <string>
#include  <stdlib.h>
#include  <unistd.h>
#include  <stdio.h>
#include <src/metadata_table.h>
#include <src/image.h>
#include <src/time.h>
#include <src/jaz/single_particle/obs_model.h>
#include "src/jaz/tomography/tomogram_set.h"

class AlignTiltseriesRunner
{
public:

    // I/O Parser
    IOParser parser;

    // Verbosity
    int verb;

    // Output rootname
    FileName fn_in, fn_out;

    // Filenames of all the tomograms to estimate the CTF from
    std::vector<long> idx_tomograms, idx_tomograms_all;

    // Information about tomography experiment
    TomogramSet tomogramSet;

    // Imod-wrapper executables
    FileName fn_imodwrapper_exe;

    // AreTomo executable
    FileName fn_aretomo_exe;

    // Use IMOD:fiducials
    bool do_imod_fiducials;

    // Nominal fiducial diameter (nm)
    RFLOAT fiducial_diam;

    // Use IMOD:patch-tracking
    bool do_imod_patchtrack;

    // Unbinned patch size (pixels)
    RFLOAT patch_size;

    // Patch overlap (percentage 0-100)
    RFLOAT patch_overlap;

    // Use AreTomo
    bool do_aretomo;

    // Resolution used in AreTomo alignment
    RFLOAT aretomo_resolution;

    // Alignment thickness in AreTomo
    RFLOAT aretomo_thickness;

    // Perform tilt angle correction in AreTomo
    bool do_aretomo_tiltcorrect;

    // Which GPU devices to use?
    int devCount;
    std::string gpu_ids;
    std::vector < std::vector < std::string > > allThreadIDs;

    // Additional gctf command line options
    std::string other_wrapper_args;

    // Continue an old run: only estimate CTF if logfile WITH Final Values line does not yet exist, otherwise skip the tomogram
    bool continue_old;

    // Process at most this number of unprocessed tomograms
    long do_at_most;

    // Do CTF estimation in aretomo?
    bool do_aretomo_ctf;

    // Do phase shift estimation in AreTomo?
    bool do_aretomo_phaseshift;

public:
    // Read command line arguments
    void read(int argc, char **argv, int rank = 0);

    // Print usage instructions
    void usage();

    // Initialise some stuff after reading
    void initialise(bool is_leader = true);

    // Execute all CTFFIND jobs to get CTF parameters
    void run();

    // Check STAR file for tomogram exists and has the correct labels
    bool checkResults(long idx_tomo);

    // Find the etomo directives file (.edf)
    bool checkEtomoDirectiveFile(long idx_tomo, FileName &filename);

    // Execute IMOD for a single tomogram
    void executeImodWrapper(long idx_tomo, int rank = 0);

    // Harvest all IMOD results into the single tomograms set starfile, and write it out
    void joinImodWrapperResults();

    // Execute AreTomo for a single tomogram
    void executeAreTomo(long idx_tomo, int rank = 0);

    // Read AreTomo results files and insert data into relion's MetaDataTable
    bool readAreTomoResults(long idx_tomo, std::string &error_message);

    // Harvest all AreTomo results into the single tomograms set starfile, and write it out
    void joinResults();


};


#endif //RELION_ALIGN_TILTSERIES_RUNNER_H

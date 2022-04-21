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

    // CTFFIND and Gctf executables and shell
    FileName fn_imodwrapper_exe, fn_shell;

    // Use IMOD:fiducials
    bool do_imod_fiducials;

    // Use IMOD:patch-tracking
    bool do_imod_patchtrack;

    // Nominal fiducial diameter (nm)
    RFLOAT fiducial_diam;

    // Unbinned patch size (pixels)
    int patch_size;

    // Patch overlap (percentage 0-100)
    RFLOAT patch_overlap;

    // Additional gctf command line options
    std::string other_wrapper_args;

    // Continue an old run: only estimate CTF if logfile WITH Final Values line does not yet exist, otherwise skip the tomogram
    bool continue_old;

    // Process at most this number of unprocessed tomograms
    long do_at_most;

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
    bool checkImodWrapperResults(long idx_tomo);

    // Execute CTFFIND for a single tomogram
    void executeImodWrapper(long idx_tomo);

    // Harvest all IMOD results into the single tomograms set starfile, and write it out
    void joinImodWrapperResults();

};


#endif //RELION_ALIGN_TILTSERIES_RUNNER_H

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

#ifndef CTFFIND_RUNNER_H_
#define CTFFIND_RUNNER_H_

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

class CtffindRunner
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Output rootname
	FileName fn_in, fn_out;

	// Estimate CTFs from rlnMicrographNameWithoutDoseWeighting instead of rlnMicrographName?
	bool do_use_without_doseweighting;

	// Filenames of all the micrographs to estimate the CTF from
	std::vector<FileName> fn_micrographs, fn_micrographs_ctf, fn_micrographs_all, fn_micrographs_ctf_all;

	// Optics groups for all micrographs
	std::vector<int> optics_group_micrographs, optics_group_micrographs_all;

	// Information about the optics groups
	ObservationModel obsModel;

    // Tilt movie index for each micrograph in a tilt serie (needed for converting back to tomographyExperiment)
    std::vector<RFLOAT> pre_exposure_micrographs, nominal_defocus_micrographs;

    // Is this a tomography experiment?
    bool is_tomo;

    // Information about tomography experiment
    TomogramSet tomogramSet;

    // Use this search range (in A) around the nominal defocus value
    RFLOAT localsearch_nominal_defocus_range;

    // Exponential factor to decrease maxres per unit dose
    RFLOAT bfactor_dose;

	// Dimension of squared area of the micrograph to use for CTF estimation
	int ctf_win;

	// CTFFIND executable and shell
	FileName fn_ctffind_exe, fn_shell;

	// Is this ctffind4?
	bool is_ctffind4;

	// Is this ctffind5?
	bool is_ctffind5;

	// Number of OMP threads for CTFFIND4
	int nr_threads;

	// Use pre-calculated power spectra
	bool use_given_ps;

	// Calculate Thon rings from movies?
	bool do_movie_thon_rings;

	// Movie rootname
	FileName movie_rootname;

	// Number of movie frames to average
	int avg_movie_frames;

	// Estimate phaseshift from a phase-plate?
	bool do_phaseshift;

	// Min, max and step phase-shift
	RFLOAT phase_min, phase_max, phase_step;

	// Continue an old run: only estimate CTF if logfile WITH Final Values line does not yet exist, otherwise skip the micrograph
	bool continue_old;

	// Process at most this number of unprocessed micrographs
	long do_at_most;

	////// CTFFIND parameters
	// Size of the box to calculate FFTw
	RFLOAT box_size;

	// Minimum and maximum resolution (in A) to be taken into account
	RFLOAT resol_min, resol_max;

	// Defocus search parameters (in A, positive is underfocus)
	RFLOAT min_defocus, max_defocus, step_defocus;

	// Amount of astigmatism (in A)
	RFLOAT amount_astigmatism;

	// Voltage (kV)
	RFLOAT Voltage;

	// Spherical aberration
	RFLOAT Cs;

	// Amplitude contrast (e.g. 0.07)
	RFLOAT AmplitudeConstrast;

	// Magnification
	RFLOAT Magnification;

	RFLOAT angpix;

	// Flag to only join results into a star file
	bool do_only_join_results;

	// Micrograph size (for Gctf check that all are equal size)
	int xdim, ydim;

	// Current working directory to make absolute-path symlinks
	std::string currdir;

	// Disable "Slower, more exhaustive search?" in CTFFIND 4.1.5-
	bool do_fast_search;

	// CTFFIND 5+ parameters
	// Determine sample thickness (CTFFIND 5)
	bool do_determine_thickness;

	// Determine sample tilt (CTFFIND 5)
	bool do_determine_tilt;

	// Dont use brute force 1D search? 
	bool no_brute_force1d;

	// Dont use 2D refinement?
	bool no_refine2d;

	// Low resolution limit for nodes
	RFLOAT node_lowres_limit;

	// High resolution limit for nodes
	RFLOAT node_highres_limit;

	// Use rounded square for nodes
	bool node_rounded_square;

	// Downweight nodes
	bool node_downweight;

	// Which GPU devices to use?
	std::string gpu_ids;
	std::vector < std::vector < std::string > > allThreadIDs;
	int devCount;

public:
	// Read command line arguments
	void read(int argc, char **argv, int rank = 0);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise(bool is_leader = true);

	// Execute all CTFFIND jobs to get CTF parameters
	void run();

	// Harvest all CTFFIND results into a single STAR file
	void joinCtffindResults();

	void getMySearchParameters(long int imic, RFLOAT &def_min, RFLOAT &def_max, RFLOAT &maxres);

    // Execute CTFFIND for a single micrograph
	void executeCtffind3(long int imic);

	// Execute CTFFIND4.1+ for a single micrograph
	void executeCtffind4(long int imic);

	// Execute CTFFIND5+ for a single micrograph
	void executeCtffind5(long int imic);

	// Get micrograph metadata
	bool getCtffindResults(FileName fn_mic, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
			RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep,
			RFLOAT &maxres, RFLOAT &valscore, RFLOAT &phaseshift, RFLOAT &icering, RFLOAT &samplethick, bool do_warn = true);
	bool getCtffind3Results(FileName fn_mic, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
			RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep,
			RFLOAT &maxres, RFLOAT &phaseshift, RFLOAT &valscore, bool do_warn = true);
	bool getCtffind4Results(FileName fn_mic, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
			RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep,
			RFLOAT &maxres, RFLOAT &phaseshift, RFLOAT &icering, bool do_warn = true);
	bool getCtffind5Results(FileName fn_mic, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
			RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep,
			RFLOAT &maxres, RFLOAT &phaseshift, RFLOAT &icering, RFLOAT &samplethick, bool do_warn = true);
};



#endif /* CTFFIND_RUNNER_H_ */

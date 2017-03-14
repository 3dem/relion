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

#ifndef MOTIONCORR_RUNNER_H_
#define MOTIONCORR_RUNNER_H_

#include  <glob.h>
#include  <vector>
#include  <string>
#include  <stdlib.h>
#include  <stdio.h>
#include <src/time.h>
#include "src/metadata_table.h"
#include "src/image.h"
#include <algorithm>

class MotioncorrRunner
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Output rootname
	FileName fn_in, fn_out, fn_movie;

	// Filenames of all the micrographs to run Motioncorr on
	std::vector<FileName> fn_micrographs, fn_ori_micrographs;

	// Use MOTIONCOR2 instead of UNBLUR?
	bool do_motioncor2;

	// First and last movie frames to use in alignment and written-out corrected average and movie (default: do all)
	int first_frame_ali, last_frame_ali, first_frame_sum, last_frame_sum;

	// Binning factor for binning inside MOTIONCORR/MOTIONCOR2
	double bin_factor;

	// B-factor for MOTIONCOR2
	double bfactor;

	// Also save the aligned movies?
	bool do_save_movies;

	// Additional arguments that need to be passed to MOTIONCORR
	FileName fn_other_motioncor2_args;

	// MOTIONCOR2 executable
	FileName fn_motioncor2_exe;

	// Voltage and dose per frame for MOTIONCOR2/UNBLUR dose-weighting
	bool do_dose_weighting;
	double voltage;
	double dose_per_frame;
	double pre_exposure;

	// Gain reference file
	FileName fn_gain_reference;

	// Defect file
	FileName fn_defect;

	// Archive directory
	FileName fn_archive;

	// Number of patches in X, Y direction for MOTIONCOR2
	int patch_x, patch_y;

	// How many frames to group in MOTIONCOR2
	int group;

	// Use Unblur/summovie instead of MOTIONCORR?
	bool do_unblur;

	// UNBLUR, SUMMOVIE executable
	FileName fn_unblur_exe, fn_summovie_exe;

	// Pixel size for UNBLUR
	double angpix;

	// Number of threads for unblur
	int nr_threads;

	// Continue an old run: only estimate CTF if logfile WITH Final Values line does not yet exist, otherwise skip the micrograph
	bool continue_old;

	// Output STAR file
	MetaDataTable MDavg, MDmov;

	// Which GPU devices to use?
	int devCount;
	std::string gpu_ids;
	std::vector < std::vector < std::string > > allThreadIDs;

public:
	// Read command line arguments
	void read(int argc, char **argv, int rank = 0);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	// Execute all MOTIONCORR jobs
	void run();

	// Given an input fn_mic filename, this function will determine the names of the output corrected image (fn_avg) and the corrected movie (fn_mov).
	void getOutputFileNames(FileName fn_mic, FileName &fn_avg, FileName &fn_mov);

    // Execute MOTIONCOR2 for a single micrograph
    bool executeMotioncor2(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts, int rank = 0);

    // Get the shifts from MOTIONCOR2
    void getShiftsMotioncor2(FileName fn_log, std::vector<float> &xshifts, std::vector<float> &yshifts);

	// Execute UNBLUR for a single micrograph
    bool executeUnblur(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts);

	// Get the shifts from UNBLUR
    void getShiftsUnblur(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts);

	// Plot the FRC curve from SUMMOVIE
	void plotFRC(FileName fn_frc);

	// Plot the shifts
    void plotShifts(FileName fn_eps, std::vector<float> &xshifts, std::vector<float> &yshifts);

	// Make a PDF file with all the shifts and write output STAR files
	void generateLogFilePDFAndWriteStarFiles();

	// Write out final STAR file
	void writeSTAR();

};


#endif /* MOTIONCORR_RUNNER_H_ */

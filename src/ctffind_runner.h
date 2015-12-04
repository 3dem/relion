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
#include "src/metadata_table.h"
#include "src/image.h"
#include <src/time.h>

class CtffindRunner
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Output rootname
	FileName fn_in, fn_out;

	// Filenames of all the micrographs to estimate the CTF from
	std::vector<FileName> fn_micrographs;

	// Dimension of squared area of the micrograph to use for CTF estimation
	int ctf_win;

	// CTFFIND and Gctf executables
	FileName fn_ctffind_exe, fn_gctf_exe;

	// use Kai Zhang's Gctf instead of CTFFIND?
	bool do_use_gctf;

	// When using Gctf, ignore CTFFIND parameters and use Gctf defaults instead?
	bool do_ignore_ctffind_params;

	// When using Gctf, use equi-phase averaging?
	bool do_EPA;

	// Continue an old run: only estimate CTF if logfile WITH Final Values line does not yet exist, otherwise skip the micrograph
	bool continue_old;

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

	// Detector pixel size (um)
	RFLOAT PixelSize;

	// For Gctf: directly provide angpix!
	RFLOAT angpix;

	// Flag to only join results into a star file
	bool do_only_join_results;

	// Micrograph size (for Gctf check that all are equal size)
	int xdim, ydim;

	// Current working directory to make absolute-path symlinks
	std::string currdir;

public:
	// Read command line arguments
	void read(int argc, char **argv, int rank = 0);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	// Execute all CTFFIND jobs to get CTF parameters
	void run();

	// Harvest all CTFFIND results into a single STAR file
	void joinCtffindResults();

	// Execute CTFFIND for a single micrograph
	void executeCtffind(long int imic);

	// Check micrograph size and add name to the list of micrographs to run Gctf on
	void addToGctfJobList(long int imic, std::string &allmicnames);

	// Execute Gctf for many micrographs
	void executeGctf(std::string &allmicnames);

	// function to remove everything before the UNIQDATE label in the input filename
	FileName getOutputFile(FileName fn_input);

};

// Get micrograph metadata
bool getCtffindResults(FileName fn_mic, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
		RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep, bool die_if_not_found = true);


#endif /* CTFFIND_RUNNER_H_ */

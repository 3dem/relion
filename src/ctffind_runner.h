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

	// CTFFIND executable
	FileName fn_ctffind_exe;

	// Continue an old run: only estimate CTF if logfile WITH Final Values line does not yet exist, otherwise skip the micrograph
	bool continue_old;

	////// CTFFIND parameters
	// Size of the box to calculate FFTw
	double box_size;

	// Minimum and maximum resolution (in A) to be taken into account
	double resol_min, resol_max;

	// Defocus search parameters (in A, positive is underfocus)
	double min_defocus, max_defocus, step_defocus;

	// Amount of astigmatism (in A)
	double amount_astigmatism;

	// Voltage (kV)
	double Voltage;

	// Spherical aberration
	double Cs;

	// Amplitude contrast (e.g. 0.07)
	double AmplitudeConstrast;

	// Magnification
	double Magnification;

	// Detector pixel size (um)
	double PixelSize;

	// Flag to only join results into a star file
	bool do_only_join_results;

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
	void executeCtffind(FileName fn_mic);

};

// Get micrograph metadata
bool getCtffindResults(FileName fn_mic, double &defU, double &defV, double &defAng, double &CC,
		double &HT, double &CS, double &AmpCnst, double &XMAG, double &DStep, bool die_if_not_found = true);


#endif /* CTFFIND_RUNNER_H_ */

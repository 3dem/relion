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
#include "src/motioncorr_own_devolved.h"

void MotioncorrOwnDevolved::read(int argc, char **argv)
{
	MotioncorrRunner::read(argc, argv);
}

void MotioncorrOwnDevolved::addClArgs()
{
	int path_section =  parser.addSection("In/out paths options");
	movie_path = parser.getOption("--in_movie", "Path to input movie");
	micrograph_path = parser.getOption("--out_mic", "Output micrograph path");
	motion_correction_star_path = parser.getOption("--mc_star", "Path to star file containing motion correction model information from a previous run", "");
	MotioncorrRunner::addClArgs();
}

void MotioncorrOwnDevolved::run()
{
	prepareGainReference(1);

	bool fromStarFile = false;
	Micrograph mic(movie_path, fn_gain_reference, bin_factor, eer_upsampling, eer_grouping);
	if (motion_correction_star_path != "") {
		Micrograph mic2(motion_correction_star_path, "", 1);
		mic = mic2;
		fromStarFile = true;
	}
	bool result;
	result = executeOwnMotionCorrection(mic, fromStarFile);
	if (result) saveModel(mic);
}

FileName MotioncorrOwnDevolved::getOutputFileNames(FileName fn_mic, bool continue_even_odd)
{
	return micrograph_path;
}

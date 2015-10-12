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


#ifndef PREPROCESSING_H_
#define PREPROCESSING_H_
#include  <glob.h>
#include  <vector>
#include  <string>
#include  <stdlib.h>
#include  <stdio.h>
#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/ctffind_runner.h"
#include <src/fftw.h>
#include <src/time.h>

class Preprocessing
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Output rootname
	FileName fn_in, fn_out;

	////////////////////////////////////// Extract particles from the micrographs
	// Perform particle extraction?
	bool do_extract;

	// Extract particles from movies instead of single micrographs
	bool do_movie_extract;

	// First frame to extract from movies
	int movie_first_frame;

	// Last frame to extract from movies
	int movie_last_frame;

	// Number of individual movie frames to average over
	int avg_n_frames;

	// Rootname to identify movies, e.g. mic001_movie.mrcs will be the movie of mic001.mrc if fn_movie="movie"
	FileName fn_movie;

	// Filenames (may include wildcards) for all coordinate files to be used for particle extraction
	FileName fns_coords_in;

	// Alternative niput: STAR file with all (selected) micrographs and their rootname for the picked coordinates files
	FileName fn_star_in, fn_pick_suffix;

	// Filenames of all the coordinate files to use for particle extraction
	std::vector<FileName> fn_coords;

	// Dimensionality of the micrographs (2 for normal micrographs, 3 for tomograms)
	int dimensionality;

	// Flag to project subtomograms along Z
	bool do_project_3d;

	// Box size to extract the particles in
	int extract_size;

	// Bias in picked coordinates in X and in Y direction (in pixels)
	RFLOAT extract_bias_x, extract_bias_y;

	////////////////////////////////////// Post-extraction image modifications
	// Perform re-scaling of extracted images
	bool do_rescale;
	int scale;

	// Perform re-windowing of extracted images
	bool do_rewindow;
	int window;

	// Perform normalization of the extract images
	bool do_normalise;

	// Perform contrast inversion of the extracted images
	bool do_invert_contrast;

	// Standard deviations to remove black and white dust
	RFLOAT white_dust_stddev, black_dust_stddev;

	// Radius of a circle in the extracted images outside of which one calculates background mean and stddev
	int bg_radius;

	// Use input stack to perform the image modifications
	FileName fn_operate_in;

	// Name of output stack (only when fn_operate in is given)
	FileName fn_operate_out;

	//////////////////////////////////// Output STAR file
	bool do_join_starfile;

public:
	// Read command line arguments
	void read(int argc, char **argv, int rank = 0);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	// General Running
	void run();

	// join all STAR files into one
	// This is done separate from runExtractParticles to allow particle extraction to be done in parallel...
	void joinAllStarFiles();

	// Extract particles from the micrographs
	void runExtractParticles();

	// Read coordinates from text files
	void readCoordinates(FileName fn_coord, MetaDataTable &MD);

	// For the given coordinate file, read the micrograph and/or movie and extract all particles
	void extractParticlesFromFieldOfView(FileName fn_coord);

	// Actually extract particles. This can be from one (average) micrgraph or from a single frame from a movie
	void extractParticlesFromOneFrame(MetaDataTable &MD,
			FileName fn_mic, int iframe, int n_frames, FileName fn_output_img_root, long int &my_current_nr_images, long int my_total_nr_images,
			RFLOAT &all_avg, RFLOAT &all_stddev, RFLOAT &all_minval, RFLOAT &all_maxval);

	// Perform per-image operations (e.g. normalise, rescaling, rewindowing and inverting contrast) on an input stack (or STAR file)
	void runOperateOnInputFile(FileName fn_perimage_in);

	// Here normalisation, windowing etc is performed on an individual image and it is written to disc
	void performPerImageOperations(Image<RFLOAT> &Ipart, FileName fn_output_img_root, int nframes, long int image_nr, long int nr_of_images,
			RFLOAT &all_avg, RFLOAT &all_stddev, RFLOAT &all_minval, RFLOAT &all_maxval);

	// Get micrograph name from the rootname
	// The rootname may have an additional string after the uniqye micrograph name
	// That way, multiple "families" of distinct particle types may be extracted from the same micrographs
	FileName getMicrographNameFromRootName(FileName fn_root);

	// The inverse of the function above
	FileName getRootNameFromMicrographName(FileName fn_mic);

};

#endif /* PREPROCESSING_H_ */

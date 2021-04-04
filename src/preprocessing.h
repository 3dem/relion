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
#include "src/ctf.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/ctffind_runner.h"
#include "src/helix.h"
#include <src/jaz/single_particle/obs_model.h>
#include <src/fftw.h>
#include <src/time.h>

class Preprocessing
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Name for directory of output Particle stacks and Particle STAR file
	FileName fn_part_dir, fn_part_star, fn_pick_star;

	// Write in float16 (MRC mode 12)?
	bool write_float16;

	// Does the input micrograph STAR file or the input data STAR file have CTF information?
	bool mic_star_has_ctf, data_star_has_ctf;

	// If input micrograph STAR file has CTF and so does input particle STAR (fn_data), then
	// by default the fn_data CTF parameters will be in the output particles.star file (as they could have been refined in ctf_refine)
	// Use the option below to let the micrographs_ctf.star CTF parameters be output. (This is useful for example, when CTF estimation has been re-run with a calibrated pixel size)
	bool keep_ctf_from_micrographs;

	/////////////////? Do phase flipping?
	bool do_phase_flip;
	bool do_premultiply_ctf;
	bool do_ctf_intact_first_peak;
	RFLOAT angpix, output_angpix, ref_angpix;

	////////////////// Extract particles from the micrographs
	// Perform particle extraction?
	bool do_extract;

	// Only extract particles when the STAR file for that micrograph doesn't exist yet
	bool only_extract_unfinished;

	// Skip gathering CTF information from the ctffind logfiles (e.g. when the info is already there from Gctf)?
	bool do_skip_ctf_logfiles;

	// STAR file with all (selected) micrographs, the suffix of the coordinates files, and the directory where the coordinate files are
	FileName fn_star_in, fn_coord_suffix, fn_coord_dir;

	// Map to go from micrograph name to coordinate file name
	std::map<FileName, FileName> micname2coordname;

	// Alternative to fn_coord_suffix: a list of all coordinate files
	FileName fn_coord_list;

	// STAR file with refined particle coordinates (to re-extract particles, for example with different binning)
	FileName fn_data;

	// How many micrographs are joined together in batches of movie-particles?
	int join_nr_mics;

	// Reset origin offsets in the input fn_data STAR file to zero
	bool do_reset_offsets;

	// Re-center particles according to rlnOriginX/Y in fn_data STAR file?
	bool do_recenter;

	// Coordinates (in pixel in the reference) to re-center re-extracted particles on
	RFLOAT recenter_x, recenter_y, recenter_z;

	// MetadataTable with all refined particle coordinates (given through fn_data)
	MetaDataTable MDimg;
	ObservationModel obsModelMic, obsModelPart;

	// Filenames of all the coordinate files to use for particle extraction
	std::vector<FileName> fn_coords;

	// Filenames of all the micrographs to use for particle extraction
	std::vector<FileName> fn_mics;

	// Metadata table with CTF information for all micrographs
	MetaDataTable MDmics;

	// Dimensionality of the micrographs (2 for normal micrographs, 3 for tomograms)
	int dimensionality;

	// Flag to project subtomograms along Z
	bool do_project_3d;

	// Box size to extract the particles in
	int extract_size;

	// Minimum threshold for autopickFigureOfMerit to extract particles
	RFLOAT extract_minimum_fom;

	// Box size to extract the particle for premultiplication with the CTF
	int premultiply_ctf_extract_size;

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

	// Subtract ramp instead of a level background in normalization
	bool do_ramp;

	// Perform contrast inversion of the extracted images
	bool do_invert_contrast;

	// Standard deviations to remove black and white dust
	RFLOAT white_dust_stddev, black_dust_stddev;

	// Radius of a circle in the extracted images outside of which one calculates background mean and stddev (in pixels)
	int bg_radius;

	// Extract helical segments
	bool do_extract_helix;

	// Outer diameter of helical tubes in Angstroms (for masks of helical segments)
	RFLOAT helical_tube_outer_diameter;

	// Extract helical segments from tube coordinates
	bool do_extract_helical_tubes;

	// Number of helical asymmetrical units
	int helical_nr_asu;

	// Helical rise in Angstroms
	RFLOAT helical_rise;

	// Add bimodal angular priors for helical segments
	bool helical_bimodal_angular_priors;

	// Cut helical tubes into segments?
	bool helical_cut_into_segments;

	// Use input stack to perform the image modifications
	FileName fn_operate_in;

	// Name of output stack (only when fn_operate in is given)
	FileName fn_operate_out;

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

	// Read helical coordinates from text files
	void readHelicalCoordinates(FileName fn_mic, FileName fn_coord, MetaDataTable &MD);

	// For the given coordinate file, read the micrograph and/or movie and extract all particles
	bool extractParticlesFromFieldOfView(FileName fn_mic, long int imic);

	// Actually extract particles. This can be from one micrgraph
	void extractParticlesFromOneMicrograph(MetaDataTable &MD,
			FileName fn_mic, int ipos, FileName fn_output_img_root, FileName fn_oristack,
			long int &my_current_nr_images, long int my_total_nr_images,
			RFLOAT &all_avg, RFLOAT &all_stddev, RFLOAT &all_minval, RFLOAT &all_maxval);

	// Perform per-image operations (e.g. normalise, rescaling, rewindowing and inverting contrast) on an input stack (or STAR file)
	void runOperateOnInputFile();

	// Here normalisation, windowing etc is performed on an individual image and it is written to disc
	// Jun24,2015 - Shaoda, extract helical segments
	void performPerImageOperations(
			Image<RFLOAT> &Ipart,
			FileName fn_output_img_root,
			long int image_nr,
			long int nr_of_images,
			RFLOAT tilt_deg,
			RFLOAT psi_deg,
			RFLOAT &all_avg,
			RFLOAT &all_stddev,
			RFLOAT &all_minval,
			RFLOAT &all_maxval);


	// Get the coordinate metadatatable from fn_data
	MetaDataTable getCoordinateMetaDataTable(FileName fn_mic);
	FileName getOutputFileNameRoot(FileName fn_mic);

};

#endif /* PREPROCESSING_H_ */

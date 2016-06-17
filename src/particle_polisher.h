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

#ifndef PARTICLE_POLISHER_H_
#define PARTICLE_POLISHER_H_
#include "src/image.h"
#include "src/metadata_table.h"
#include "src/exp_model.h"
#include "src/fftw.h"
#include "src/time.h"
#include "src/mask.h"
#include "src/funcs.h"
#include "src/backprojector.h"
#include "src/ctf.h"
#include "src/postprocessing.h"
#include "src/CPlot2D.h"

#define LINEAR_FIT 0
#define LOGARITHMIC_FIT 1
#define SQRT_FIT 2
#define NO_FIT 3

class ParticlePolisher
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Input & Output rootname
	FileName fn_in, fn_out, fn_sym, fn_mask;

	// list of individual micrographs
	std::vector<FileName> fn_mics;

	// Standard deviation for a Gaussian-weight on the distance between particles in the micrograph
	RFLOAT sigma_neighbour_distance;

	// Maximum resolution in pre-frame reconstructions
	RFLOAT perframe_highres;

	// Flag to indicate all calculations have to be repeated from scratch
	// if false, then intermediate files are re-read from disc and earlier calculations are skipped
	bool only_do_unfinished;

	// List of frame-numbers for each movie-particle in the Experiment
	std::vector<int> movie_frame_numbers;

	// Which fitting mode (lienar/logarithmic/nofit)
	int fitting_mode;

	// CTF stuff for the reconstructions
	bool do_ctf, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak;

	// Pixel size (for B-factors)
	RFLOAT angpix;

	// Original image size
	int ori_size;

	// Skip B-factor weighting
	bool do_weighting;

	// Minimum resolution (in Angstrom) for fitting of B-factor in Guinier plot
	RFLOAT fit_minres;

	// Width of a running average window for the single-frame reconstructions
	int frame_running_average;

	// Vector with the B-factors for all individual frames
	MultidimArray<RFLOAT> perframe_bfactors;

	// Fitted movement coordinates for all input images
	MultidimArray<RFLOAT> fitted_movements;

	// Image with the mask (used for relative weighting of each frame)
	Image<RFLOAT> Imask;

	// FSC curve of the masked, averages of all single-frame reconstructions
	MultidimArray<RFLOAT> fsc_average;

	// Metadatatable with the information from the polished particles
	MetaDataTable MDshiny;

	// Reference volume reconstructed from the initially-polished particles to be used for per-particle CTF-refinement and beamtilt-refinement
	//Projector PPrefvol_half1, PPrefvol_half2;

	// Normalise the polished particles?
	bool do_normalise;

	// Subtract a ramping background in the normalisation?
	bool do_ramp;

	// Radius of the background-circle for noise normalisation (in pixels)
	int bg_radius;

	// Sigma-levels for dust removal
	RFLOAT white_dust_stddev, black_dust_stddev;

	// Maximum useful resolution in the reconstruction
	RFLOAT maxres_model;

	// Maximum beam tilt to analyse, and step-size to sample in X and Y
	RFLOAT beamtilt_max, beamtilt_step;

	// Number of sampled beamtilts
	int nr_sampled_beam_tilts;

	// Names of the data sets to be separated in the beamtilt refinement
	std::vector<FileName> fn_beamtilt_groups;

	// Minimum resolution to take beamtilt into account
	RFLOAT minres_beamtilt;

	// Weighted squared-differences for all beamtilts
	MultidimArray<RFLOAT> diff2_beamtilt;

	// Weighted squared-differences for all defocus values
	MultidimArray<RFLOAT> defocus_shift_allmics;

	// Optimal beamtilts for each data set
	std::vector<Matrix1D<RFLOAT> > best_beamtilts;

	// Per-particle CTF optimisation
	RFLOAT defocus_shift_max, defocus_shift_step;

	// Sep24,2015 - Shaoda, Helical reconstruction
	bool is_helix;

	int helical_nr_asu;

	RFLOAT helical_twist, helical_rise;

	// Make output directories only if it doesn't exist yet
	FileName fn_olddir;

public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	// Generate a list of all individual micrographs
	void generateMicrographList();

	// Fit the beam-induced translations for all average micrographs
	void fitMovementsAllMicrographs();

	// Fit a function through th observed movements
	void fitMovementsOneMicrograph(long int imic);

	// Get the B-factor for all single-frame reconstruction
	void calculateAllSingleFrameReconstructionsAndBfactors();

	// Read/write of STAR file with all per-frame B-factors
	bool readStarFileBfactors(FileName fn_star);
	void writeStarFileBfactors(FileName fn_star);

	// Write out an additional STAR file with the resolution-dependent, relative weights per frame
	void writeStarFileRelativeWeights(FileName fn_star);

	// Get the B-factor for a single-frame reconstruction
	void calculateSingleFrameReconstruction(int iframe, int ihalf);

	// Run standard post-processing (only unmasked FSC  on the single-frame reconstruction.
	void postProcessSingleFrameReconstruction(int this_frame);

	// Calculate the B-factors for a single-frame reconstruction
	void calculateBfactorSingleFrameReconstruction(int this_frame, RFLOAT &bfactor, RFLOAT &offset, RFLOAT &corr_coeff);

	// Calculate the average of all single-frame rconstructions (for a given half)
	void calculateAverageAllSingleFrameReconstructions(int ihalf);

	// Change the name of the particle stack to the output directory
	void changeParticleStackName(FileName &fn_part);

	// Movie frame re-alignment for a single micrograph
	void polishParticlesOneMicrograph(long int imic);

	// Movie frame re-alignment for all micrographs
	void polishParticlesAllMicrographs();

	// Write out the resulting STAR files
	void writeStarFilePolishedParticles();

	// Reconstruct the two independent halves of the shiny particles
	void reconstructShinyParticlesAndFscWeight(int ipass);

	// Reconstruct one half of the shiny particles
	void reconstructShinyParticlesOneHalf(int ihalf, Experiment &exp_model);

	// Generate a single PDF file with motion tracks, B-factor plots etc
	void generateLogFilePDF();

	// General Running
	void run();


};



#endif /* PARTICLE_POLISHER_H_ */

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

#ifndef ML_MODEL_H_
#define ML_MODEL_H_
#include "src/projector.h"
#include "src/backprojector.h"
#include "src/metadata_table.h"
#include "src/exp_model.h"
#include "src/healpix_sampling.h"

#define ML_BLOB_ORDER 0
#define ML_BLOB_RADIUS 1.9
#define ML_BLOB_ALPHA 15

class MlModel
{
public:

	// Dimension of the references (2D or 3D)
	int ref_dim;

	// Dimension of the data (2D or 3D)
	int data_dim;

	// Original size of the images
	int ori_size;

	// Pixel size (in Angstrom)
	RFLOAT pixel_size;

	// Current size of the images to be used in the expectation
	int current_size;

	// Current resolution (in 1/Ang)
	RFLOAT current_resolution;

	// Number of classes
	int nr_classes;

	// Number of independent bodies for multi-body refinement
	int nr_bodies;

	// Number of image groups with separate sigma2_noise spectra
	int nr_groups;

	// Number of particles in each group
	std::vector<long int> nr_particles_group;

	// Number of directions (size of pdf_direction);
	int nr_directions;

	// Log-likelihood target value
	RFLOAT LL;

	// Padding factor
	int padding_factor;

	// Fourier space interpolator
	int interpolator;

	// Minimum number of shells to perform linear interpolation
	int r_min_nn;

	// Average Pmax of the normalised probability distributions
	RFLOAT ave_Pmax;

	// Average normalisation correction factor
	RFLOAT avg_norm_correction;

	// Variance in the origin offsets
	RFLOAT sigma2_offset;

	// Fudge factor to adjust estimated tau2_class spectra
	RFLOAT tau2_fudge_factor;

	// Vector with all reference images
	std::vector<MultidimArray<RFLOAT> > Iref;

	// Vector with masks for all bodies in multi-body refinement
	std::vector<MultidimArray<RFLOAT> > masks_bodies;

	// Vector with center-of-mass coordinates for all bodies in multi-body refinement
	std::vector<Matrix1D<RFLOAT> > com_bodies;

	// One projector for each class;
	std::vector<Projector > PPref;

	// One name for each group
	std::vector<FileName> group_names;

	// One noise spectrum for each group
	std::vector<MultidimArray<RFLOAT > > sigma2_noise;

	// One intensity scale for each group
	std::vector<RFLOAT> scale_correction;

	// One intensity B-factor for each group
	std::vector<RFLOAT> bfactor_correction;

	// Prior information: one restrained power_class spectrum for each class (inverse of right-hand side in Wiener-filter-like update formula)
	std::vector<MultidimArray<RFLOAT > > tau2_class;

	// Radial average of the estimated variance in the reconstruction (inverse of left-hand side in Wiener-filter-like update formula)
	std::vector<MultidimArray<RFLOAT > > sigma2_class;

	// FSC spectra between random halves of the data
	MultidimArray<RFLOAT > fsc_halves_class;

	// One likelihood vs prior ratio spectrum for each class
	std::vector<MultidimArray<RFLOAT > > data_vs_prior_class;

	// One value for each class
	std::vector<RFLOAT > pdf_class;

	// One array for each class
	std::vector<MultidimArray<RFLOAT> > pdf_direction;

	// Priors for offsets for each class (only in 2D)
	std::vector<Matrix1D<RFLOAT> > prior_offset_class;

	// Mode for orientational prior distributions
	int orientational_prior_mode;

	// Variance in rot angle for the orientational pdf
	RFLOAT sigma2_rot;

	// Variance in tilt angle for the orientational pdf
	RFLOAT sigma2_tilt;

	// Variance in psi angle for the orientational pdf
	RFLOAT sigma2_psi;

	// Estimated accuracy at which rotations can be assigned, one for each class
	std::vector<RFLOAT> acc_rot;

	// Estimated accuracy at which translations can be assigned, one for each class
	std::vector<RFLOAT> acc_trans;

	// Spectral contribution to orientability of individual particles, one for each class
	std::vector<MultidimArray<RFLOAT > > orientability_contrib;


public:

	// Constructor
	MlModel():
		ref_dim(0),
		data_dim(0),
		ori_size(0),
		pixel_size(0),
		current_size(0),
		current_resolution(0),
		nr_classes(0),
		nr_groups(0),
		nr_directions(0),
		LL(0),
		padding_factor(0),
		interpolator(0),
		r_min_nn(0),
		ave_Pmax(0),
		avg_norm_correction(0),
		sigma2_offset(0),
		tau2_fudge_factor(0),
		orientational_prior_mode(0),
		sigma2_rot(0),
		sigma2_tilt(0),
		sigma2_psi(0)
	{
		clear();
	}

	// Destructor
	~MlModel()
	{
		clear();
	}

	// Clear everything
	void clear()
	{
		Iref.clear();
		masks_bodies.clear();
		com_bodies.clear();
		PPref.clear();
		group_names.clear();
		sigma2_noise.clear();
		scale_correction.clear();
		bfactor_correction.clear();
		tau2_class.clear();
		fsc_halves_class.clear();
		sigma2_class.clear();
		data_vs_prior_class.clear();
		pdf_class.clear();
		pdf_direction.clear();
		nr_particles_group.clear();
		ref_dim = ori_size = nr_classes = nr_bodies = nr_groups = nr_directions = interpolator = r_min_nn = padding_factor = 0;
		ave_Pmax = avg_norm_correction = LL = sigma2_offset = tau2_fudge_factor = 0.;
		sigma2_rot = sigma2_tilt = sigma2_psi = 0.;
		acc_rot.clear();
		acc_trans.clear();
		orientability_contrib.clear();
	}

	// Initialise vectors with the right size
	void initialise();

	//Read a model from a file
	void read(FileName fn_in);

	// Write a model to disc
	void write(FileName fn_out, HealpixSampling &sampling, bool do_write_bild = true);

	//Read a tau-spectrum from a STAR file
	void readTauSpectrum(FileName fn_tau, int verb);

	// Read images from disc and initialise
	// Also set do_average_unaligned and do_generate_seeds flags
	void readImages(FileName fn_ref, int _ori_size, Experiment &_mydata,
			bool &do_average_unaligned, bool &do_generate_seeds, bool &refs_are_ctf_corrected);

	// Given the Experiment of the already expanded dataset of movieframes, expand the current MlModel to contain all movie frames
	// Make a new group for each unique rlnGroupName in the expanded Experiment, copying the values from the groups in the current MlModel
	// For that: remove "00000i@" as well as movie extension from the rlnGroupName in the expanded Experiment and compare with group_names in current MlModel
	void expandToMovieFrames(Experiment &moviedataexpand, int running_avg_side);

	RFLOAT getResolution(int ipix)	{ return (RFLOAT)ipix/(pixel_size * ori_size); }

	RFLOAT getResolutionAngstrom(int ipix)	{ return (ipix==0) ? 999. : (pixel_size * ori_size)/(RFLOAT)ipix; }

	int getPixelFromResolution(RFLOAT resol)	{ return (int)(resol * pixel_size * ori_size); }

	/** Initialise pdf_orient arrays to the given size
	* If the pdf_orient vectors were empty, resize them to the given size and initialise with an even distribution
	* If they were not empty, check that the new size is equal to the old one, and otherwise throw an exception
	* because one cannot use an old pdf_orient with size unequal to the new one
	*/
	void initialisePdfDirection(int newsize);

	/** Read in the binary masks provided by the user and then make a soft edge on those */
	void initialiseBodyMasks(FileName fn_masks, FileName fn_root_out);

	// Set FourierTransforms in Projector of each class
	// current_size will determine the size of the transform (in number of Fourier shells) to be held in the projector
	void setFourierTransformMaps(bool update_tau2_spectra, int nr_threads = 1);

	/* Initialises the radial average of the data-versus-prior ratio
	 */
	void initialiseDataVersusPrior(bool fix_tau);

};

class MlWsumModel: public MlModel
{
public:
	// One backprojector for CTF-corrected estimate of each class;
	std::vector<BackProjector > BPref;

	// Store the sum of the weights inside each group
	// That is the number of particles inside each group
	std::vector<RFLOAT> sumw_group;

	// For the refinement of group intensity scales and bfactors
	// For each group store weighted sums of experimental image times reference image as a function of resolution
	std::vector<MultidimArray<RFLOAT > > wsum_signal_product_spectra;

	// For each group store weighted sums of squared reference as a function of resolution
	std::vector<MultidimArray<RFLOAT > > wsum_reference_power_spectra;

	// Constructor
	MlWsumModel()
	{
		clear();
	}

	// Destructor
	~MlWsumModel()
	{
		clear();
	}

	// Clear everything
	void clear()
	{
		BPref.clear();
		sumw_group.clear();
		MlModel::clear();
	}

	// Initialise all weighted sums (according to size of corresponding model
	void initialise(MlModel &_model, FileName fn_sym = "c1");

	// Initialize all weighted sums to zero (with resizing the BPrefs to current_size)
	void initZeros();

	// Pack entire structure into one large MultidimArray<RFLOAT> for reading/writing to disc
	// To save memory, the model itself will be cleared after packing.
	void pack(MultidimArray<RFLOAT> &packed);

	// Fill the model again using unpack (this is the inverse operation from pack)
	void unpack(MultidimArray<RFLOAT> &packed);

	// Pack entire structure into one large MultidimArray<RFLOAT> for shipping over with MPI
	// To save memory, the model itself will be cleared after packing.
    // If the whole thing becomes bigger than 1Gb (see MAX_PACK_SIZE in ml_model.cpp), then break it up into pieces because MPI cannot handle very large messages
	// When broken up: nr_pieces > 1
	void pack(MultidimArray<RFLOAT> &packed, int &piece, int &nr_pieces, bool do_clear=true);

	// Fill the model again using unpack (this is the inverse operation from pack)
	void unpack(MultidimArray<RFLOAT> &packed, int piece, bool do_clear=true);

};

#endif /* ML_MODEL_H_ */

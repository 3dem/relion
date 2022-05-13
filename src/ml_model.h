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
#include "src/gradient_optimisation.h"

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

	// Number of image groups with separate scale corrections
	int nr_groups;

	// Number of optics groups for separate sigma2_noise spectra
	int nr_optics_groups;

	// Keep track of the first and/or second moment of the gradient
	bool do_grad;

	// Number of particles in each (micrograph) group
	std::vector<long int> nr_particles_per_group;

	// Number of directions (size of pdf_direction);
	long long int nr_directions;

	// Log-likelihood target value
	RFLOAT LL;

	// Padding factor
	RFLOAT padding_factor;

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

	// Vector with all gradient moments
	std::vector<MultidimArray<Complex> > Igrad1;
	std::vector<MultidimArray<Complex> > Igrad2;

	// Vector with masks for all bodies in multi-body refinement
	std::vector<MultidimArray<RFLOAT> > masks_bodies;

	// Vector with center-of-mass coordinates for all bodies in multi-body refinement
	std::vector<Matrix1D<RFLOAT> > com_bodies;

	// Vector with 2D matrices that pre-orient all bodies in multi-body refinement
	std::vector<Matrix2D<RFLOAT> > orient_bodies;

	// Vector with directions around which to rotate each body in multi-body refinement
	std::vector<Matrix1D<RFLOAT> > rotate_direction_bodies;

	// One projector for each class;
	std::vector<Projector > PPref;
	std::vector<bool> PPrefRank;

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

	// FSC spectra between random halves of the data (multiple ones for each body in multibody-refinement)
	std::vector< MultidimArray<RFLOAT > > fsc_halves_class;

	// One likelihood vs prior ratio spectrum for each class
	std::vector<MultidimArray<RFLOAT > > data_vs_prior_class;

	// One Fourier-coverage spectrum for each class
	std::vector<MultidimArray<RFLOAT > > fourier_coverage_class;

	// One value for each class
	std::vector<RFLOAT > pdf_class;
	std::vector<RFLOAT > class_age;

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

	// Stddev in tilt angle for the orientational pdf of each body
	std::vector<RFLOAT> sigma_tilt_bodies;

	// Stddev in psi angle for the orientational pdf of each body
	std::vector<RFLOAT> sigma_psi_bodies;

	// Stddev in offsets for the orientational pdf of each body
	std::vector<RFLOAT> sigma_offset_bodies;

	// Is this body kept fixed in refinement?
	std::vector<int> keep_fixed_bodies;

	// Maximum radius of mask (in Angstrom!)
	std::vector<int> max_radius_mask_bodies;

	// 2D Matrix with pointers to the PPrefs for overlapping bodies
	MultidimArray<int> pointer_body_overlap;

	std::vector<int> pointer_body_overlap_inv;

	// Estimated accuracy at which rotations can be assigned, one for each class
	std::vector<RFLOAT> acc_rot;

	// Estimated accuracy at which translations can be assigned, one for each class
	std::vector<RFLOAT> acc_trans;

	// The estimate resolution, one for each class
	std::vector<RFLOAT> estimated_resolution;

	// Fourier coverage up to the estimate resolution, one  for each class
	std::vector<RFLOAT> total_fourier_coverage;

	// Spectral contribution to orientability of individual particles, one for each class
	std::vector<MultidimArray<RFLOAT > > orientability_contrib;

	// Nov20,2015 - Shaoda, Helical refinement
	bool is_helix;

	// Number of helical asymmetrical units
	int helical_nr_asu;

	// Helical twist (in degrees)
	std::vector<RFLOAT> helical_twist;

	// Helical rise (in Angstroms)
	std::vector<RFLOAT> helical_rise;

	// Self-organizing map
	SomGraph som;
	int last_som_add_iter;

	// Search range of helical twist (in degrees)
	RFLOAT helical_twist_min,  helical_twist_max, helical_twist_inistep;

	// Search range of helical rise (in Angstroms)
	RFLOAT helical_rise_min,  helical_rise_max, helical_rise_inistep;

	// Normalize overlapping regions in multibody masks
	bool norm_body_mask_overlap;

	// Process data on GPU
	bool do_gpu;

	bool pseudo_halfsets;

	// Store filenames of references for Liyi's class feature program
	std::vector<FileName> ref_names;

public:

	// Constructor
	MlModel():
		ref_dim(0),
		data_dim(0),
		ori_size(0),
		pixel_size (0),
		current_size(0),
		current_resolution(0),
		nr_classes(0),
		nr_bodies(0),
		nr_groups(0),
		nr_optics_groups(0),
		nr_directions(0),
		LL(0),
		padding_factor(0.),
		interpolator(0),
		r_min_nn(0),
		ave_Pmax(0),
		avg_norm_correction(0),
		sigma2_offset(0),
		tau2_fudge_factor(0),
		orientational_prior_mode(0),
		sigma2_rot(0),
		sigma2_tilt(0),
		sigma2_psi(0),
		is_helix(0),
		helical_nr_asu(1),
		helical_twist_min(0),
		helical_twist_max(0),
		helical_twist_inistep(0),
		helical_rise_min(0),
		helical_rise_max(0),
		helical_rise_inistep(0),
		norm_body_mask_overlap(false),
		som(),
		last_som_add_iter(0),
		do_gpu(false),
		pseudo_halfsets(false)
	{
		clear();
	}

	// Destructor
	~MlModel()
	{
		clear();
	}

	/** Assignment operator
	 */
	MlModel& operator =(const MlModel &MD)
	{
		if (this != &MD)
		{
			clear();
			ref_dim = MD.ref_dim;
			data_dim = MD.data_dim;
			ori_size = MD.ori_size;
			pixel_size = MD.pixel_size;
			current_size = MD.current_size;
			current_resolution = MD.current_resolution;
			nr_classes = MD.nr_classes;
			nr_bodies = MD.nr_bodies;
			nr_groups = MD.nr_groups;
			nr_optics_groups = MD.nr_optics_groups;
			do_grad = MD.do_grad;
			pseudo_halfsets = MD.pseudo_halfsets;
			nr_directions = MD.nr_directions;
			LL = MD.LL;
			padding_factor = MD.padding_factor;
			interpolator = MD.interpolator;
			r_min_nn = MD.r_min_nn;
			ave_Pmax = MD.ave_Pmax;
			avg_norm_correction = MD.avg_norm_correction;
			sigma2_offset = MD.sigma2_offset;
			tau2_fudge_factor = MD.tau2_fudge_factor;
			orientational_prior_mode = MD.orientational_prior_mode;
			sigma2_rot = MD.sigma2_rot;
			sigma2_tilt = MD.sigma2_tilt;
			sigma2_psi = MD.sigma2_psi;
			is_helix = MD.is_helix;
			helical_nr_asu = MD.helical_nr_asu;
			helical_twist_min = MD.helical_twist_min;
			helical_twist_max = MD.helical_twist_max;
			helical_twist_inistep = MD.helical_twist_inistep;
			helical_rise_min = MD.helical_rise_min;
			helical_rise_max = MD.helical_rise_max;
			helical_rise_inistep= MD.helical_rise_inistep;
			Iref = MD.Iref;
			Igrad1 = MD.Igrad1;
			Igrad2 = MD.Igrad2;
			masks_bodies = MD.masks_bodies;
			com_bodies = MD.com_bodies;
			orient_bodies = MD.orient_bodies;
			sigma_tilt_bodies = MD.sigma_tilt_bodies;
			sigma_psi_bodies = MD.sigma_psi_bodies;
			sigma_offset_bodies = MD.sigma_offset_bodies;
			keep_fixed_bodies = MD.keep_fixed_bodies;
			max_radius_mask_bodies = MD.max_radius_mask_bodies;
			PPref = MD.PPref;
			PPrefRank = MD.PPrefRank;
			group_names = MD.group_names;
			sigma2_noise = MD.sigma2_noise;
			scale_correction = MD.scale_correction;
			bfactor_correction = MD.bfactor_correction;
			tau2_class = MD.tau2_class;
			sigma2_class = MD.sigma2_class;
			fsc_halves_class = MD.fsc_halves_class;
			data_vs_prior_class = MD.data_vs_prior_class;
			fourier_coverage_class = MD.fourier_coverage_class;
			pdf_class = MD.pdf_class;
			class_age = MD.class_age;
			pdf_direction = MD.pdf_direction;
			prior_offset_class = MD.prior_offset_class;
			nr_particles_per_group = MD.nr_particles_per_group;
			acc_rot = MD.acc_rot;
			acc_trans = MD.acc_trans;
			estimated_resolution = MD.estimated_resolution;
			total_fourier_coverage = MD.total_fourier_coverage;
			orientability_contrib = MD.orientability_contrib;
			helical_twist = MD.helical_twist;
			helical_rise = MD.helical_rise;
			do_gpu = MD.do_gpu;
			pseudo_halfsets = MD.pseudo_halfsets;
			ref_names = MD.ref_names;
	        }
        	return *this;
	}

	// Clear everything
	void clear()
	{
		Iref.clear();
		Igrad1.clear();
		Igrad2.clear();
		masks_bodies.clear();
		com_bodies.clear();
		orient_bodies.clear();
		sigma_tilt_bodies.clear();
		sigma_psi_bodies.clear();
		sigma_offset_bodies.clear();
		keep_fixed_bodies.clear();
		max_radius_mask_bodies.clear();
		PPref.clear();
		PPrefRank.clear();
		group_names.clear();
		sigma2_noise.clear();
		scale_correction.clear();
		bfactor_correction.clear();
		tau2_class.clear();
		fsc_halves_class.clear();
		sigma2_class.clear();
		data_vs_prior_class.clear();
		fourier_coverage_class.clear();
		prior_offset_class.clear();
		pdf_class.clear();
		class_age.clear();
		pdf_direction.clear();
		nr_particles_per_group.clear();
		ref_dim = data_dim = ori_size = nr_classes = nr_bodies = nr_groups = nr_directions = interpolator = r_min_nn;
		padding_factor = 0.;
		ave_Pmax = avg_norm_correction = LL = sigma2_offset = tau2_fudge_factor = 0.;
		sigma2_rot = sigma2_tilt = sigma2_psi = 0.;
		acc_rot.clear();
		acc_trans.clear();
		estimated_resolution.clear();
		total_fourier_coverage.clear();
		orientability_contrib.clear();
		helical_twist.clear();
		helical_rise.clear();
		ref_names.clear();
		do_grad=false;
		pseudo_halfsets=false;
	}

	// Initialise vectors with the right size
	void initialise(bool _do_grad = false, bool _pseudo_halfsets = false);

	//Read a model from a file
	void read(FileName fn_in, int nr_optics_groups_from_mydata, bool _do_grad=false, bool _pseudo_halfsets=false);

	// Write a model to disc
	void write(FileName fn_out, HealpixSampling &sampling,
			bool do_write_bild = true, bool do_only_write_images = false);

	//Read a tau-spectrum from a STAR file
	void readTauSpectrum(FileName fn_tau, int verb);

	// Read images from disc and initialise
	// Also set do_average_unaligned and do_generate_seeds flags
	void initialiseFromImages(FileName fn_ref, bool _is_3d_model, Experiment &_mydata,
			bool &do_average_unaligned, bool &do_generate_seeds, bool &refs_are_ctf_corrected,
			RFLOAT ref_angpix = -1., bool _do_grad = false, bool _pseudo_halfsets = false, bool do_trust_ref = false, bool verb = false);

	RFLOAT getResolution(int ipix)	{ return (RFLOAT)ipix/(pixel_size * ori_size); }

	RFLOAT getResolutionAngstrom(int ipix)	{ return (ipix==0) ? 999. : (pixel_size * ori_size)/(RFLOAT)ipix; }

	int getPixelFromResolution(RFLOAT resol)	{ return (int)ROUND(resol * pixel_size * ori_size); }

	/** Initialise pdf_orient arrays to the given size
	* If the pdf_orient vectors were empty, resize them to the given size and initialise with an even distribution
	* If they were not empty, check that the new size is equal to the old one, and otherwise throw an exception
	* because one cannot use an old pdf_orient with size unequal to the new one
	*/
	void initialisePdfDirection(long long int newsize);

	/** Read in the binary masks provided by the user and then make a soft edge on those */
	void initialiseBodies(FileName fn_masks, FileName fn_root_out, bool also_initialise_rest = false, int rank = 0);

	/** Write out a Bild file with the COMs and directions or rotation for each body */
	void writeBildFileBodies(FileName fn_bild);

	// Set FourierTransforms in Projector of each class
	// current_size will determine the size of the transform (in number of Fourier shells) to be held in the projector ( thisClass == -1  => do all classes this call)
	void setFourierTransformMaps(bool update_tau2_spectra, int nr_threads = 1, RFLOAT strict_lowres_exp = -1,
			   const MultidimArray<RFLOAT> *fourier_mask = NULL );

	// current_size will determine the size of the transform (in number of Fourier shells) to be held in the projector ( thisClass == -1  => do all classes this call)
	void setFourierTransformMaps(bool update_tau2_spectra, std::vector<bool> ListCheapSetup, int nr_threads = 1, RFLOAT strict_lowres_exp = -1);

	/* Initialises the radial average of the data-versus-prior ratio
	 */
	void initialiseDataVersusPrior(bool fix_tau);

	void initialiseHelicalParametersLists(RFLOAT _helical_twist, RFLOAT _helical_rise);

	void calculateTotalFourierCoverage();

	void reset_class(int class_idx, int to_class_idx = -1);
};

class MlWsumModel: public MlModel
{
public:
	// One backprojector for CTF-corrected estimate of each class;
	std::vector<BackProjector > BPref;

	// Store the sum of the weights inside each optics group
	// That is the number of particles inside each optics group
	std::vector<RFLOAT> sumw_group;

    // Resolution-dependent sum of CTF^2 for ctf_premultiplied correction of tau2 estimates
    std::vector<MultidimArray<RFLOAT> > sumw_ctf2;

    // Resolution-dependent sum of multiplicities for subtomogram averaging
    std::vector<MultidimArray<RFLOAT> > sumw_stMulti;

	// For the refinement of group intensity scales and bfactors
	// For each group store weighted sums of experimental image times reference image as a function of resolution
	std::vector<RFLOAT > wsum_signal_product;

	// For each group store weighted sums of squared reference as a function of resolution
	std::vector<RFLOAT > wsum_reference_power;

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
        sumw_stMulti.clear();
        sumw_ctf2.clear();
        wsum_signal_product.clear();
        wsum_reference_power.clear();
		MlModel::clear();
	}

	// Initialise all weighted sums (according to size of corresponding model
	void initialise(MlModel &_model, FileName fn_sym = "c1", bool asymmetric_padding = false, bool _skip_gridding = false, bool _pseudo_halfsets = false);

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

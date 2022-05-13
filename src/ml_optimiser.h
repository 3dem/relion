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

#ifndef ML_OPTIMISER_H_
#define ML_OPTIMISER_H_

#ifdef ALTCPU
	#include <tbb/enumerable_thread_specific.h>
	#include <complex>
#endif

#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include "src/ml_model.h"
#include "src/parallel.h"
#include "src/exp_model.h"
#include "src/ctf.h"
#include "src/time.h"
#include "src/mask.h"
#include "src/healpix_sampling.h"
#include "src/helix.h"
#include "src/local_symmetry.h"
#include "src/acc/settings.h"

#define ML_SIGNIFICANT_WEIGHT 1.e-8
#define METADATA_LINE_LENGTH METADATA_LINE_LENGTH_ALL

#define METADATA_ROT 0
#define METADATA_TILT 1
#define METADATA_PSI 2
#define METADATA_XOFF 3
#define METADATA_YOFF 4
#define METADATA_ZOFF 5
#define METADATA_CLASS 6
#define METADATA_DLL 7
#define METADATA_PMAX 8
#define METADATA_NR_SIGN 9
#define METADATA_NORM 10

#define METADATA_CTF_DEFOCUS_U 11
#define METADATA_CTF_DEFOCUS_V 12
#define METADATA_CTF_DEFOCUS_ANGLE 13
#define METADATA_CTF_BFACTOR 14
#define METADATA_CTF_KFACTOR 15
#define METADATA_CTF_PHASE_SHIFT 16

#define METADATA_ROT_PRIOR 17
#define METADATA_TILT_PRIOR 18
#define METADATA_PSI_PRIOR 19
#define METADATA_XOFF_PRIOR 20
#define METADATA_YOFF_PRIOR 21
#define METADATA_ZOFF_PRIOR 22
#define METADATA_PSI_PRIOR_FLIP_RATIO 23
#define METADATA_ROT_PRIOR_FLIP_RATIO 24 	// KThurber

#define METADATA_LINE_LENGTH_BEFORE_BODIES 25
#define METADATA_NR_BODY_PARAMS 6

#define DO_WRITE_DATA true
#define DONT_WRITE_DATA false
#define DO_WRITE_SAMPLING true
#define DONT_WRITE_SAMPLING false
#define DO_WRITE_MODEL true
#define DONT_WRITE_MODEL false
#define DO_WRITE_OPTIMISER true
#define DONT_WRITE_OPTIMISER false

#define WIDTH_FMASK_EDGE 2.
#define MAX_NR_ITER_WO_RESOL_GAIN 1
#define MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES 1
#define MAX_NR_ITER_WO_RESOL_GAIN_GRAD 4

// for profiling
//#define TIMING

class MlOptimiser;

class MlOptimiser
{
public:

	// For GPU-maps
	std::vector<int> cudaDevices;
	std::vector<int> cudaOptimiserDeviceMap;
	std::vector<void*> cudaOptimisers;
	std::vector<void*> accDataBundles;

#ifdef ALTCPU
	std::vector<void*> cpuOptimisers;

	// Container for TBB thread-local data
	typedef tbb::enumerable_thread_specific< void * > CpuOptimiserType;

	CpuOptimiserType   tbbCpuOptimiser;

	std::complex<XFLOAT> **mdlClassComplex __attribute__((aligned(64)));
#endif


	// I/O Parser
	IOParser parser;

	// Experimental metadata model
	Experiment mydata;

	// Current ML model
	MlModel mymodel;

	// Current weighted sums
	MlWsumModel wsum_model;

	// HEALPix sampling object for coarse sampling
	HealpixSampling sampling;

	// Filename for the experimental images
	FileName fn_data;

	// Output root filename
	FileName fn_out;

	// Filename for input reference images (stack, star or image)
	FileName fn_ref;

	// Filename for the subtomo optimiser set. Used to set fn_data, fn_ref and fn_mask
	FileName fn_OS;

	// Optimiser set table for subtomo
	MetaDataTable optimisationSet;

	// Generate a 3D model from 2D particles de novo?
	bool is_3d_model;

	//User-specified pixel size for the reference
	RFLOAT ref_angpix;

	// Filename for input tau2-spectrum
	FileName fn_tau;

	// Filename for input masks for multi-body refinement
	FileName fn_body_masks;

	// Initialise bodies for a new multi-body refinement
	bool do_initialise_bodies;

	// Use subtracted images for reconstructions in multi-body refinement?
	bool do_reconstruct_subtracted_bodies;

	// Precalculated rotation matrix for (0,90,0) rotation, and its transpose
	Matrix2D<RFLOAT> A_rot90, A_rot90T;

	// Flag to keep tau-spectrum constant
	bool fix_tau;

	// some parameters for debugging
	RFLOAT debug1, debug2, debug3;

	// Starting and finishing particles (for parallelisation)
	long int my_first_particle_id, my_last_particle_id;

	// Total number iterations and current iteration
	int iter, nr_iter, nr_iter_max;

	// Flag whether to split data from the beginning into two random halves
	bool do_split_random_halves;

	// Debug flag to process 1 half when doing do_split_random_halves, without using MPI
	int debug_split_random_half;

	// For safe-guarding the gold-standard separation
	int my_halfset;

	// resolution (in Angstrom) to join the two random halves
	RFLOAT low_resol_join_halves;

	// Flag to join random halves again
	bool do_join_random_halves;

	// Flag to always join random halves, this is a developmental option for testing of sub-optimal FSC-usage only!
	bool do_always_join_random_halves;

	// Flag whether to use different images for recnostruction and alignment
	bool do_use_reconstruct_images;

	// Flag whether to do CTF correction
	bool do_ctf_correction;

	// Do not correct CTFs until after the first peak
	bool intact_ctf_first_peak;

	// Calculate CTF in 2x larger window, and then rescale it to mimic PSF of boxing out small particles (this deals with CTF aliaising)
	bool do_ctf_padding;

	// Pnly perform phase-flipping CTF correction
	bool only_flip_phases;

	// Images have been CTF phase flipped al;ready
	bool ctf_phase_flipped;

	// Flag whether current references are ctf corrected
	bool refs_are_ctf_corrected;

	// Flag whether directions have changed upon continuing an old run
	bool directions_have_changed;

	// Flag whether to do image-wise intensity-scale correction
	bool do_norm_correction;

	// Flag whether to do group-wise intensity bfactor correction
	bool do_scale_correction;

	// Flag whether to use the auto-refine procedure
	bool do_auto_refine;

	// Flag whether to use auto-sampling (outside auto_refine)
	bool do_auto_sampling;

	// Flag whether to ignore changes in hidden variables in auto-refine (which makes it faster)
	bool auto_ignore_angle_changes;

	// Flag whether to proceed faster with resolution-based angular sampling updates in auto-refine
	bool auto_resolution_based_angles;

	// Flag whether to accelerate early iterations of 2D/3D classification with coarser samplings
	bool allow_coarser_samplings;

	// Force auto-refinement to converge
	bool do_force_converge;

	// Number of iterations without a resolution increase
	int nr_iter_wo_resol_gain;

	// Best resolution obtained thus far
	RFLOAT best_resol_thus_far;

	// Is the FSC still high at the resolution limit?
	bool has_high_fsc_at_limit;

	// How many iterations ago was there a big step in the increase of current_size
	int has_large_incr_size_iter_ago;

	// In auto-sampling mode, use local searches from this sampling rate on
	int autosampling_hporder_local_searches;

	// Smallest changes thus far in the optimal translational offsets, orientations and classes
	RFLOAT smallest_changes_optimal_offsets;
	RFLOAT smallest_changes_optimal_orientations;
	int smallest_changes_optimal_classes;

	// Number of iterations without a decrease in OffsetChanges
	int nr_iter_wo_large_hidden_variable_changes;

	// Strict resolution limits in the expectation step
	RFLOAT strict_highres_exp, strict_lowres_exp;

	// Flag to indicate to estimate angular accuracy until current_size (and not coarse_size) when restricting high-res limit in the expectation step
	// This is used for testing purposes only
	bool do_acc_currentsize_despite_highres_exp;

	// Global parameters to store accuracy on rot and trans
	RFLOAT acc_rot, acc_trans;

	// Flag to indicate to use all data out to Nyquist
	bool do_use_all_data;

	// Flag to indicate the refinement has converged
	bool has_converged;

	// Flag to indicate that angular sampling in auto-sampling has reached its limit
	bool has_fine_enough_angular_sampling;

	// Minimum and maximum angular sampling to achieve in auto-refinement (in degrees)
	RFLOAT minimum_angular_sampling, maximum_angular_sampling;

	// Flag to keep sigma2_offset fixed
	bool fix_sigma_offset;

	// Flag to keep sigma2_noise fixed
	bool fix_sigma_noise;

	//  Use images only up to a certain resolution in the expectation step (one for each optics_group)
	// image_coarse_size is for first pass, image_current_size is for second pass, image_full_size is original image size
	std::vector<int> image_coarse_size, image_current_size, image_full_size;

	// Use images only up to a certain resolution in the expectation step
	int max_coarse_size;

	// Particle diameter (in Ang)
	RFLOAT particle_diameter;

	// How many fourier shells should be included beyond the highest shell where evidenceVsPriorRatio < 1?
	int incr_size;

	// Minimum resolution to perform Bayesian estimates of the model
	int minres_map;

	// Abort when resolution reaches beyond this value
	float abort_at_resolution;

	// Flag to flatten solvent
	bool do_solvent;

	// Filename for a user-provided mask
	FileName fn_mask;

	// Correct gold-standard FSC for solvent mask?
	bool do_phase_random_fsc;

	// Filename for a user-provided second solvent mask
	// This solvent mask will have its own average density and may be useful for example to fill the interior of an icosahedral virus
	FileName fn_mask2;

	// Mask for regions to be low-pass filtered
	FileName fn_lowpass_mask;
	RFLOAT lowpass;

	// Width of the soft-edges of the circular masks
	int width_mask_edge;

	// Number of particles to be processed simultaneously
	int nr_pool;

	//////////////// Gradient optimisation
	// If current refinement is gradient based
	bool gradient_refine;

	// Id current iteration is gradient based
	bool do_grad;
	bool do_grad_next_iter;

	// Do pseudo half-sets to estimate noise
	bool grad_pseudo_halfsets;

	// Number of iterations at the end of a gradient refinement using Expectation-Maximization
	int grad_em_iters;

	// Number of iterations in the initial phase of refinement
	int grad_ini_iter;

	// Fraction of iterations in the initial phase of refinement
	RFLOAT grad_ini_frac;

	// Number of iterations in the final phase of refinement
	int grad_fin_iter;

	// Fraction of iterations in the initial phase of refinement
	RFLOAT grad_fin_frac;

	// Number of iterations between the initial and the final ones
	// (during which a linear transform from grad_ini_resol->grad_fin_resol and grad_ini_subset_size->grad_fin_subset_size will be done)
	int grad_inbetween_iter;

	// Size of the subsets used in the initial iterations
	int grad_ini_subset_size;

	//Size of the subsets used in the final iterations
	int grad_fin_subset_size;

	// Effective size of subsets
	int effective_setsize;

	// Adjusting the under-estimation of signal of low-resolution components during gradient optimization to this resolution
	RFLOAT grad_min_resol;

	// The resolution in the initial iterations
	RFLOAT grad_ini_resol; // in A

	// The resolution in the final iterations
	RFLOAT grad_fin_resol; // in A

	// Skip annealing of multiple reference in SGD
	// (by default refs are kept the same during grad_nr_iter_initial and then slowly annealed during grad_nr_iter_inbetween)
	bool do_grad_skip_anneal;

	// Momentum update parameter
	RFLOAT mu;

	// Step size of the gradient updates
	RFLOAT grad_stepsize;
	RFLOAT grad_current_stepsize;
	std::string grad_stepsize_scheme;

	// For automated adjusting of tau2 fudge
	std::string tau2_fudge_scheme;
	RFLOAT tau2_fudge_arg;

	// Self-organizing map
	bool do_init_blobs;
	bool do_som;
	bool is_som_iter;
	int som_starting_nodes;
	float som_connectivity;
	float som_neighbour_pull;
	float som_inactivity_threshold;

	float class_inactivity_threshold;

	// Size of the random subsets
	long int subset_size;

	// Automatic doubling of subset size during auto refinement
	// Usage: subset_size = grad_ini_subset * 2 ^ auto_subest_size_order
	int auto_subset_size_order;

	// Gradient auto refinement has converged
	bool grad_has_converged;

	// Suspended finer sampling for this many iteration
	int grad_suspended_finer_sampling_iter;

	// Suspended finer sampling with local searches for this many iteration
	int grad_suspended_local_searches_iter;

	// Every how many iterations should be written to disk when using subsets
	int write_every_grad_iter;

	// Use subsets like in cisTEM to speed up 2D/3D classification
	bool do_fast_subsets;

	// Available memory (in Gigabyte)
	size_t available_gpu_memory;
	size_t requested_free_gpu_memory;

	// Perform combination of weight through files written on disc
	bool combine_weights_thru_disc;

	// Only print metadata label definitions and exit
	bool do_print_metadata_labels;

	// Use parallel access to disc?
	bool do_parallel_disc_io;

	// Use gpu resources?
	bool do_gpu;
	bool anticipate_oom;

	// Use alternate cpu implementation
	bool do_cpu;

	// Which GPU devices to use?
	std::string gpu_ids;

	// Or preread all images into RAM on the leader node?
	bool do_preread_images;

	// Place on scratch disk to copy particle stacks temporarily
	FileName fn_scratch;

	// Amount of scratch space to be left free (in Gb)
	RFLOAT keep_free_scratch_Gb;

	// Re-use data on scratch dir, i.e. dont delete data already there and copy again
	bool do_reuse_scratch;

	// Don't delete scratch after finishing
	bool keep_scratch;

	// Print the symmetry transformation matrices
	bool do_print_symmetry_ops;

	/** Name of the multiple symmetry groups */
	std::vector<FileName> fn_multi_sym;

	/** Perform reconstruction outside of relion_refine, e.g. for learned priors */
	bool do_external_reconstruct;

	/* Flag whether to use the Adaptive approach as by Tagare et al (2010) J. Struc. Biol.
	 * where two passes through the integrations are made: a first one with a coarse angular sampling and
	 * smaller images and a second one with finer angular sampling and larger images
	 * Only the orientations that accumulate XXX% (adaptive_fraction) of the top weights from the first pass are oversampled in the second pass
	 * Only these oversampled orientations will then contribute to the weighted sums
	 * Note that this is a bit different from Tagare et al, in that we do not perform any B-spline interpolation....
	 *
	 * For (adaptive_oversampling==0): no adaptive approach is being made, there is only one pass
	 * For (adaptive_oversampling==1): the coarse sampling grid is 2x coarse than the fine one
	 * For (adaptive_oversampling==2): the coarse sampling grid is 4x coarse than the fine one
	 * etc..
	 *
	 * */
	int adaptive_oversampling;

	/* Fraction of the weights for which to consider the finer sampling
	 * The closer to one, the more orientations will be oversampled
	 * The default is 0.999.
	 */
	RFLOAT adaptive_fraction;

	// Seed for random number generator
	int random_seed;

	/* Flag to indicate orientational (i.e. rotational AND translational) searches will be skipped */
	bool do_skip_align;

	/* Flag to indicate rotational searches will be skipped */
	bool do_skip_rotate;

	/* Flag to indicate maximization step will be skipped: only data.star file will be written out */
	bool do_skip_maximization;

	/* Flag to only sample tilt angles (from -180->180), ignore rot angles: temporary fix for DNA-origami frames */
	bool do_only_sample_tilt;

	/* Flag to use bimodal prior distributions on psi (2D classification of helical segments) */
	bool do_bimodal_psi;

	/* Flag to center classes */
	bool do_center_classes;

	//////// Special stuff for the first iterations /////////////////

	// Skip marginalisation in first iteration and use signal cross-product instead of Gaussian
	bool do_firstiter_cc;

	/// Always perform cross-correlation instead of marginalization
	bool do_always_cc;

	// Initial low-pass filter for all references (in digital frequency)
	RFLOAT ini_high;

	// Flag whether to generate seeds
	// TODO: implement!
	bool do_generate_seeds;

	// Flag whether to calculate the average of the unaligned images (for 2D refinements)
	bool do_average_unaligned;

	// Flag whether to calculate initial sigma_noise spectra
	bool do_calculate_initial_sigma_noise;

	// Flag to switch off error message about normalisation
	bool dont_raise_norm_error;

	///////////// Helical symmetry /////////////
	// Flag whether to do helical refinement
	bool do_helical_refine;

	// Ignore helical symmetry?
	bool ignore_helical_symmetry;

	// Initial helical twist in degrees
	RFLOAT helical_twist_initial;

	// Initial helical rise in Angstroms
	RFLOAT helical_rise_initial;

	// N-number of an N-start helix (only for rotational priors, a la KThurber)
	int helical_nstart;

	// Only expand this amount of Z axis proportion to full when imposing real space helical symmetry
	RFLOAT helical_z_percentage;

	// Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)
	RFLOAT helical_tube_inner_diameter;

	// Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)
	RFLOAT helical_tube_outer_diameter;

	// Flag whether to do local refinement of helical parameters
	bool do_helical_symmetry_local_refinement;

	// Sigma of distance along the helical tracks (in Angstroms)
	RFLOAT helical_sigma_distance;

	// Keep helical tilt priors fixed (at 90 degrees) in global angular searches?
	bool helical_keep_tilt_prior_fixed;

	// Apply fourier_mask for helical refinements
	std::string helical_fourier_mask_resols;
	FileName fn_fourier_mask;
	MultidimArray<RFLOAT> helical_fourier_mask;

	///////// Hidden stuff, does not work with read/write: only via command-line ////////////////

	// Skip gridding in reconstruction
	bool skip_gridding;

	// Number of iterations for gridding preweighting reconstruction
	int gridding_nr_iter;

	// Flag whether to do group-wise B-factor correction or not
	bool do_bfactor;

	// Flag whether to use maximum a posteriori (MAP) estimation
	bool do_map;

	// For debugging/analysis: Name for initial sigma_noise sepctrum;
	FileName fn_sigma;

	// Multiplicative fdge factor for the sigma estimates
	RFLOAT sigma2_fudge;

	// Perform two reconstructions of random halves sequentially (to save memory for very big cases)
	bool do_sequential_halves_recons;

	// Do zero-filled soft-masks instead of noisy masks on experimental particles
	// This will increase SNR but introduce correlations that are not modelled...
	// Until now the best refinements have used the noisy mask, not the soft mask....
	bool do_zero_mask;

	// Prepare for automated 2D class average selection: nice to have access to unmasked reference
	bool do_write_unmasked_refs;

	/////////// Keep track of hidden variable changes ////////////////////////

	// Changes from one iteration to the next in the angles
	RFLOAT current_changes_optimal_orientations, sum_changes_optimal_orientations;

	// Changes from one iteration to the next in the translations
	RFLOAT current_changes_optimal_offsets, sum_changes_optimal_offsets;

	// Changes from one iteration to the next in the class assignments
	RFLOAT current_changes_optimal_classes, sum_changes_optimal_classes;

	// Just count how often the optimal changes are summed
	RFLOAT sum_changes_count;

	///////// Special stuff for subtomogram avg //////////

	// Do CTF3D files contain CTF^2 or not?
	bool ctf3d_squared;
	// Apply regular Probability and Noise estimation
	bool do_skip_subtomo_correction;
	// Have subtomogram files been multiplicity normalised during reconstruction?
	bool normalised_subtomos;
	// Threshold applied to multiplicity volume to mask edge voxels
	RFLOAT subtomo_multi_thr;

	/////////// Some internal stuff ////////////////////////

	// Array with pointers to the resolution of each point in a Fourier-space FFTW-like array (one for each optics_group)
	std::vector<MultidimArray<int> > Mresol_fine, Mresol_coarse;
	MultidimArray<int> Npix_per_shell;

	// Verbosity flag
	int verb;

	// Thread Managers for the expectation step: one for all (pooled) particles
	ThreadTaskDistributor *exp_ipart_ThreadTaskDistributor;

	// Number of threads to run in parallel
	int x_pool;
	int nr_threads;

	//for catching exceptions in threads
	RelionError * threadException;

	long int exp_my_first_part_id, exp_my_last_part_id;
	MultidimArray<RFLOAT> exp_metadata, exp_imagedata;
	std::string exp_fn_img, exp_fn_ctf, exp_fn_recimg;
	std::vector<MultidimArray<RFLOAT> > exp_imgs;
	std::vector<int> exp_random_class_some_particles;

	// Calculate translated images on-the-fly
	bool do_shifts_onthefly;
	std::vector< std::vector<MultidimArray<Complex> > > global_fftshifts_ab_coarse, global_fftshifts_ab_current, global_fftshifts_ab2_coarse, global_fftshifts_ab2_current;

	//TMP DEBUGGING
	MultidimArray<RFLOAT> DEBUGGING_COPY_exp_Mweight;

	//Use different padding factors for the projector (PF+0) and backprojector (PF+1)
	bool asymmetric_padding;

	//Maximum number of significant weights in coarse pass of expectation specified by user
	int maximum_significants_arg;

	//Maximum number of significant weights in coarse pass of expectation used in refinement
	// This can be set by user or automatically
	int maximum_significants;

	//Do regularization by denoising
	bool do_red;

	// Tabulated sine and cosine values (for 3D helical sub-tomogram averaging with on-the-fly shifts)
	TabSine tab_sin;
	TabCosine tab_cos;

	// Local symmetry information (STAR or plain-text format)
	FileName fn_local_symmetry;

	// Local symmetry - list of masks
	std::vector<FileName> fn_local_symmetry_masks;

	// Local symmetry - list of operators
	std::vector<std::vector<Matrix1D<RFLOAT> > > fn_local_symmetry_operators;

	//Maximum number of particles permitted to be drop, due to zero sum of weights, before exiting with an error (GPU only).
	int failsafe_threshold;

	// Trust box size and angpix for the input reference
	bool do_trust_ref_size;

	// Number of particles for initial sigma2_noise estimation per optics group
	int minimum_nr_particles_sigma2_noise;

	// Do not call makeGoodHelixForEachRef
    bool skip_realspace_helical_sym;

#ifdef TIMING
	Timer timer;
	int TIMING_DIFF_PROJ, TIMING_DIFF_SHIFT, TIMING_DIFF_DIFF2;
	int TIMING_WSUM_PROJ, TIMING_WSUM_BACKPROJ, TIMING_WSUM_DIFF2, TIMING_WSUM_SUMSHIFT;
	int TIMING_EXP, TIMING_MAX, TIMING_RECONS, TIMING_SOLVFLAT, TIMING_UPDATERES;
	int TIMING_EXP_1,TIMING_EXP_1a,TIMING_EXP_2,TIMING_EXP_3,TIMING_EXP_4,TIMING_EXP_4a,TIMING_EXP_4b,TIMING_EXP_4c,TIMING_EXP_4d,TIMING_EXP_5,TIMING_EXP_6,TIMING_EXP_7,TIMING_EXP_8,TIMING_EXP_9;
	int TIMING_ESP, TIMING_ESP_THR, TIMING_ESP_ONEPART, TIMING_ESP_ONEPARTN, TIMING_EXP_METADATA, TIMING_EXP_CHANGES;
	int TIMING_ESP_FT, TIMING_ESP_INI, TIMING_ESP_DIFF1, TIMING_ESP_DIFF2;
	int TIMING_ESP_DIFF2_A, TIMING_ESP_DIFF2_B, TIMING_ESP_DIFF2_C, TIMING_ESP_DIFF2_D, TIMING_ESP_DIFF2_E;
	int TIMING_ESP_PREC1, TIMING_ESP_PREC2, TIMING_ESP_PRECW, TIMING_WSUM_GETSHIFT, TIMING_DIFF2_GETSHIFT, TIMING_WSUM_SCALE, TIMING_WSUM_LOCALSUMS;
	int TIMING_ESP_WEIGHT1, TIMING_ESP_WEIGHT2, TIMING_WEIGHT_EXP, TIMING_WEIGHT_SORT, TIMING_ESP_WSUM;
	int TIMING_EXTRA1, TIMING_EXTRA2, TIMING_EXTRA3;
	int TIMING_EXP_SETUP, TIMING_ITER_HELICALREFINE;
	int TIMING_ITER_WRITE, TIMING_ITER_LOCALSYM;

	int RCT_1, RCT_2, RCT_3, RCT_4, RCT_5, RCT_6, RCT_7, RCT_8;
	int ES_A, ES_B, ES_C, ES_D, ES_E, ES_F;
	int ESS_1, ESS_2, ESS_3, ESS_4;
	int EINS_1, EINS_2, EINS_3, EINS_4, EINS_5,EINS_6, EINS_7, EINS_8, EINS_9;
	int EINS_10, EINS_11;
#endif

public:

	MlOptimiser():
            do_zero_mask(0),
            do_write_unmasked_refs(0),
            do_generate_seeds(0),
            sum_changes_count(0),
            current_changes_optimal_orientations(0),
            do_average_unaligned(0),
            sigma2_fudge(0),
            do_always_cc(0),
            best_resol_thus_far(0),
            debug1(0),
            incr_size(0),
            has_large_incr_size_iter_ago(0),
            nr_iter_wo_resol_gain(0),
            nr_pool(0),
            refs_are_ctf_corrected(0),
            has_high_fsc_at_limit(0),
            do_acc_currentsize_despite_highres_exp(0),
            low_resol_join_halves(0),
            do_auto_refine(0),
            do_auto_sampling(0),
            has_converged(0),
            only_flip_phases(0),
            gridding_nr_iter(0),
            do_use_reconstruct_images(0),
            fix_sigma_noise(0),
            current_changes_optimal_offsets(0),
            smallest_changes_optimal_classes(0),
            do_print_metadata_labels(0),
            adaptive_fraction(0),
            do_print_symmetry_ops(0),
            do_bfactor(0),
            do_use_all_data(0),
            minres_map(0),
            debug2(0),
            do_always_join_random_halves(0),
            my_first_particle_id(0),
            x_pool(1),
            nr_threads(0),
            do_shifts_onthefly(0),
            exp_ipart_ThreadTaskDistributor(0),
            do_parallel_disc_io(0),
            sum_changes_optimal_orientations(0),
            do_solvent(0),
            strict_highres_exp(0),
            sum_changes_optimal_classes(0),
            acc_trans(0),
            width_mask_edge(0),
            has_fine_enough_angular_sampling(0),
            sum_changes_optimal_offsets(0),
            do_scale_correction(0),
            ctf_phase_flipped(0),
            nr_iter_wo_large_hidden_variable_changes(0),
            adaptive_oversampling(0),
            nr_iter(0),
            intact_ctf_first_peak(0),
            do_join_random_halves(0),
            do_skip_align(0),
            do_calculate_initial_sigma_noise(0),
            fix_sigma_offset(0),
            do_firstiter_cc(0),
            do_bimodal_psi(0),
            do_center_classes(0),
            exp_my_last_part_id(0),
            particle_diameter(0),
            smallest_changes_optimal_orientations(0),
            verb(0),
            do_norm_correction(0),
            fix_tau(0),
            directions_have_changed(0),
            acc_rot(0),
            do_sequential_halves_recons(0),
            do_skip_rotate(0),
            current_changes_optimal_classes(0),
            do_skip_maximization(0),
            dont_raise_norm_error(0),
            do_map(0),
            combine_weights_thru_disc(0),
            smallest_changes_optimal_offsets(0),
            exp_my_first_part_id(0),
            iter(0),
            my_last_particle_id(0),
            ini_high(0),
            do_ctf_correction(0),
            max_coarse_size(0),
            autosampling_hporder_local_searches(0),
            do_split_random_halves(0),
            my_halfset(-1),
            debug_split_random_half(0),
            random_seed(0),
            do_gpu(0),
            anticipate_oom(0),
            do_helical_refine(0),
            do_preread_images(0),
            ignore_helical_symmetry(0),
            helical_twist_initial(0),
            helical_rise_initial(0),
            helical_z_percentage(0),
            helical_tube_inner_diameter(0),
            helical_tube_outer_diameter(0),
            do_helical_symmetry_local_refinement(0),
            helical_sigma_distance(0),
            helical_keep_tilt_prior_fixed(0),
            ctf3d_squared(0),
            do_skip_subtomo_correction(0),
            normalised_subtomos(0),
            subtomo_multi_thr(0),
		//directional_lowpass(0),
		asymmetric_padding(false),
            maximum_significants(-1),
            threadException(NULL),
            do_init_blobs(false),
            do_som(false),
            is_som_iter(false),
            som_starting_nodes(0),
            som_connectivity(0),
            som_inactivity_threshold(0),
            som_neighbour_pull(0),
            class_inactivity_threshold(0),
            gradient_refine(false),
            do_grad(false),
            grad_em_iters(0),
            grad_ini_iter(0),
            grad_ini_frac(0),
            grad_fin_iter(0),
            grad_fin_frac(0),
            grad_inbetween_iter(0),
            subset_size(0),
            grad_ini_subset_size(0),
            grad_fin_subset_size(0),
            grad_min_resol(0),
            grad_ini_resol(0),
            grad_fin_resol(0),
            do_grad_skip_anneal(false),
            mu(0),
            grad_stepsize(-1),
            grad_current_stepsize(0),
            auto_subset_size_order(0),
            grad_has_converged(false),
            grad_suspended_local_searches_iter(-1),
            grad_suspended_finer_sampling_iter(-1),
            grad_pseudo_halfsets(false),
            skip_realspace_helical_sym(false),
#ifdef ALTCPU
		mdlClassComplex(NULL),
#endif
            failsafe_threshold(40),
            do_trust_ref_size(0),
            minimum_nr_particles_sigma2_noise(1000)
	{
#ifdef ALTCPU
		tbbCpuOptimiser = CpuOptimiserType((void*)NULL);
#endif
	};

	/** ========================== I/O operations  =========================== */
	/// Print help message
	void usage();

	/// Interpret command line
	void read(int argc, char **argv, int rank = 0);

	/// Interpret command line for the initial start of a run
	void parseInitial(int argc, char **argv);

	/// Interpret command line for the re-start of a run
	void parseContinue(int argc, char **argv);

	/// Read from STAR file
	void read(FileName fn_in, int rank = 0, bool do_prevent_preread = false);

	// Write files to disc
	void write(bool do_write_sampling, bool do_write_data, bool do_write_optimiser, bool do_write_model,
			int random_subset = 0);

    /** ========================== Initialisation  =========================== */

	// Initialise the whole optimiser
	void initialise();

	// Check the mask is thr ight size
	void checkMask(FileName &_fn_mask, int solvent_nr, int rank);

	// Some general stuff that is shared between MPI and sequential code in the early stages of initialisation
	void initialiseGeneral(int rank = 0);

	// Randomise particle processing order and resize metadata array
	void initialiseWorkLoad();

	// Fetch initial noise spectra, from input file or calculate from subset of input images,
	// also calculate initial Iref images if fn_ref == "None"
	void initialiseSigma2Noise();

	// Initialise initial reference images
	void initialiseReferences();

	/* Calculates the sum of all individual power spectra and the average of all images for initial sigma_noise estimation
	 * The rank is passed so that if one splits the data into random halves one can know which random half to treat
	 */
	void calculateSumOfPowerSpectraAndAverageImage(MultidimArray<RFLOAT> &Mavg, bool myverb = true);

	/** Use the sum of the individual power spectra to calculate their average and set this in sigma2_noise
	 * Also subtract the power spectrum of the average images,
	 * and if (do_average_unaligned) then also set Mavg to all Iref
	 */
	void setSigmaNoiseEstimatesAndSetAverageImage(MultidimArray<RFLOAT> &Mavg);

	/* Perform an initial low-pass filtering of the references
	 * Note that because of the MAP estimation, this is not necessary inside the refinement
	 */
	void initialLowPassFilterReferences();

	/** ========================== EM-Iteration  ================================= */

	/* Launch threads and set up task distributors*/
	void iterateSetup();

	/* Delete threads and task distributors */
	void iterateWrapUp();

	/* Perform expectation-maximization iterations */
	void iterate();

	/* Expectation step: loop over all particles
	 */
	void expectation();

	/* Setup expectation step. We divide the heavy steps over mpi followers,
	 * so each call needs a list of which to skip heavy setup for. For
	 * these classes, only some formatting is done. Data is copied
	 * explicitly later.*/
	void expectationSetup(std::vector<bool> cheapSetup);

	void expectationSetup();

	/* Check whether everything fits into memory, possibly adjust nr_pool and setup thread task managers */
	void expectationSetupCheckMemory(int myverb = 1);

	/* For on-the-fly shifts, precalculates AB-matrices */
	void precalculateABMatrices();

	/* Perform the expectation integration over all k, phi and series elements for a number (some) of pooled particles
	 * The number of pooled particles is determined by max_nr_pool and some memory checks in expectationSetup()
	 */
	void expectationSomeParticles(long int my_first_particle, long int my_last_particle);

	/* Perform expectation step for some particles using threads */
	void doThreadExpectationSomeParticles(int thread_id);

	/* Perform the expectation integration over all k, phi and series elements for a given particle */
	void expectationOneParticle(long int part_id_sorted, int thread_id);

	/* Function to call symmetrise of BackProjector helical objects for each class or body
	 * Do rise and twist for all asymmetrical units in Fourier space
	 * */
	void symmetriseReconstructions();

	/* Apply local symmetry according to a list of masks and their operators
	 * */
	void applyLocalSymmetryForEachRef();

	/* Apply helical symmetry (twist and rise) to the central Z slices of real space references. (Do average for every particle)
	 * Then elongate and fill the perfect helix along Z axis.
	 * */
	void makeGoodHelixForEachRef();

    /*
     * Calculate the average CTF^2 from the summed weights in wsum_model and return boolean whether to
     * correct tau2 estimates by the average CTF^2 (which is needed when using ctf_premultiplied image
     * because the tau2 is estimated from the CTF-corrected map, whereas the sigma2_noise estimates are
     * from the CTF-affected differences)
     */
    bool setAverageCTF2(MultidimArray<RFLOAT> &avgctf2);

	/* Maximization step
	 * Updates the current model: reconstructs and updates all other model parameter
	 */
	void maximization();

	/* Perform the actual reconstructions
	 * This is officially part of the maximization, but it is separated because of parallelisation issues.
	 */
	void maximizationReconstructClass(int iclass);

	/* Update gradient related parameters, returns class index that should be skipped during SOM
	 */
	int maximizationGradientParameters();

	/* Updates all other model parameters (besides the reconstructions)
	 */
	void maximizationOtherParameters();


    /** ========================== Deeper functions ================================= */

	/* Apply a solvent flattening to a map
	 */
	void solventFlatten();

	/* Center classes based on their center-of-mass
	 * and also update the origin offsets in the _data.star file correspondingly
	 */
	void centerClasses();

	/* Updates the current resolution (from data_vs_prior array) and keeps track of best resolution thus far
	 *  and precalculates a 2D Fourier-space array with pointers to the resolution of each point in a FFTW-centered array
	 */
	void updateCurrentResolution();

	/* Update the current and coarse image size
	 *  and precalculates a 2D Fourier-space array with pointers to the resolution of each point in a FFTW-centered array
	 */
	void updateImageSizeAndResolutionPointers();

	/* Read image and its metadata from disc (threaded over all pooled particles)
	 */
	void getFourierTransformsAndCtfs(long int part_id, int ibody, int metadata_offset,
			std::vector<MultidimArray<Complex > > &exp_Fimg,
			std::vector<MultidimArray<Complex > > &exp_Fimg_nomask,
			std::vector<MultidimArray<RFLOAT> > &exp_Fctf,
			std::vector<Matrix1D<RFLOAT> > &exp_old_offset,
			std::vector<Matrix1D<RFLOAT> > &exp_prior,
			std::vector<MultidimArray<RFLOAT> > &exp_power_img,
			std::vector<RFLOAT> &exp_highres_Xi2_img,
			std::vector<int> &exp_pointer_dir_nonzeroprior,
			std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<RFLOAT> &exp_directions_prior,
			std::vector<RFLOAT> &exp_psi_prior,
			std::vector<MultidimArray<RFLOAT> > &exp_STweight);

	/* Store all shifted FourierTransforms in a vector
	 * also store precalculated 2D matrices with 1/sigma2_noise
	 */
	void precalculateShiftedImagesCtfsAndInvSigma2s(bool do_also_unmasked, bool is_for_store_wsums, long int my_particle,
			int exp_current_oversampling, int metadata_offset,
			int exp_itrans_min, int exp_itrans_max,
			std::vector<MultidimArray<Complex > > &exp_Fimg,
			std::vector<MultidimArray<Complex > > &exp_Fimg_nomask,
			std::vector<MultidimArray<RFLOAT> > &exp_Fctf,
			std::vector<std::vector<MultidimArray<Complex > > > &exp_local_Fimgs_shifted,
			std::vector<std::vector<MultidimArray<Complex > > > &exp_local_Fimgs_shifted_nomask,
			std::vector<MultidimArray<RFLOAT> > &exp_local_Fctf,
			std::vector<RFLOAT> &exp_local_sqrtXi2,
			std::vector<MultidimArray<RFLOAT> > &exp_local_Minvsigma2,
			std::vector<MultidimArray<RFLOAT> > &exp_STweight,
			std::vector<MultidimArray<RFLOAT> > &exp_local_STMulti);

	// Given exp_Mcoarse_significant, check for iorient whether any of the particles has any significant (coarsely sampled) translation
	bool isSignificantAnyImageAnyTranslation(long int iorient,
			int exp_itrans_min, int exp_itrans_max, MultidimArray<bool> &exp_Mcoarse_significant);

	// Get squared differences for all iclass, idir, ipsi and itrans...
	void getAllSquaredDifferences(long int part_id, int ibody,
			int exp_ipass, int exp_current_oversampling, int metadata_offset,
			int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
			int exp_itrans_min, int exp_itrans_max, int my_iclass_min, int my_iclass_max,
			std::vector<RFLOAT> &exp_min_diff2,
			std::vector<RFLOAT> &exp_highres_Xi2_img,
			std::vector<MultidimArray<Complex > > &exp_Fimg,
			std::vector<MultidimArray<RFLOAT> > &exp_Fctf,
			MultidimArray<RFLOAT> &exp_Mweight,
			MultidimArray<bool> &exp_Mcoarse_significant,
			std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior,
			std::vector<std::vector<MultidimArray<Complex > > > &exp_local_Fimgs_shifted,
			std::vector<MultidimArray<RFLOAT> > &exp_local_Minvsigma2,
			std::vector<MultidimArray<RFLOAT> > &exp_local_Fctf,
			std::vector<RFLOAT> &exp_local_sqrtXi,
			std::vector<MultidimArray<RFLOAT> > &exp_STweight);

	// Convert all squared difference terms to weights.
	// Also calculates exp_sum_weight and, for adaptive approach, also exp_significant_weight
	void convertAllSquaredDifferencesToWeights(long int part_id, int ibody, int exp_ipass,
			int exp_current_oversampling, int metadata_offset,
			int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
			int exp_itrans_min, int exp_itrans_max, int my_iclass_min, int my_iclass_max,
			MultidimArray<RFLOAT> &exp_Mweight, MultidimArray<bool> &exp_Mcoarse_significant,
			std::vector<RFLOAT> &exp_significant_weight, std::vector<RFLOAT> &exp_sum_weight,
			std::vector<Matrix1D<RFLOAT> > &exp_old_offset, std::vector<Matrix1D<RFLOAT> > &exp_prior, std::vector<RFLOAT> &exp_min_diff2,
			std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior);

	// Store all relevant weighted sums, also return optimal hidden variables, max_weight and dLL
	void storeWeightedSums(long int part_id, int ibody,
			int exp_current_oversampling, int metadata_offset,
			int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
			int exp_itrans_min, int exp_itrans_max, int my_iclass_min, int my_iclass_max,
			std::vector<RFLOAT> &exp_min_diff2,
			std::vector<RFLOAT> &exp_highres_Xi2_img,
			std::vector<MultidimArray<Complex > > &exp_Fimg,
			std::vector<MultidimArray<Complex > > &exp_Fimg_nomask,
			std::vector<MultidimArray<RFLOAT> > &exp_Fctf,
			std::vector<MultidimArray<RFLOAT> > &exp_power_img,
			std::vector<Matrix1D<RFLOAT> > &exp_old_offset,
			std::vector<Matrix1D<RFLOAT> > &exp_prior,
			MultidimArray<RFLOAT> &exp_Mweight,
			MultidimArray<bool> &exp_Mcoarse_significant,
			std::vector<RFLOAT> &exp_significant_weight,
			std::vector<RFLOAT> &exp_sum_weight,
			std::vector<RFLOAT> &exp_max_weight,
			std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior,
			std::vector<std::vector<MultidimArray<Complex > > > &exp_local_Fimgs_shifted,
			std::vector<std::vector<MultidimArray<Complex > > > &exp_local_Fimgs_shifted_nomask,
			std::vector<MultidimArray<RFLOAT> > &exp_local_Minvsigma2,
			std::vector<MultidimArray<RFLOAT> > &exp_local_Fctf,
			std::vector<RFLOAT> &exp_local_sqrtXi2,
			std::vector<MultidimArray<RFLOAT> > &exp_STweight);

	/** Monitor the changes in the optimal translations, orientations and class assignments for some particles */
	void monitorHiddenVariableChanges(long int my_first_part_id, long int my_last_part_id);

	// Updates the overall changes in the hidden variables and keeps track of nr_iter_wo_large_changes_in_hidden_variables
	void updateOverallChangesInHiddenVariables();

	// Calculate expected error in orientational assignments
	// Based on comparing projections of the model and see how many degrees apart gives rise to difference of power > 3*sigma^ of the noise
	void calculateExpectedAngularErrors(long int my_first_part_id, long int my_last_part_id);

	// Adjust angular sampling based on the expected angular accuracies for auto-refine procedure
	void updateAngularSampling(bool verb = true);

	// Adjust subset size in fast_subsets or Gradient algorithms
	void updateSubsetSize(bool verb = true);

	// Adjust step size for the gradient algorithms
	void updateStepSize();

	// Adjust tau2 fudge factor
	void updateTau2Fudge();

	// Check convergence for auto-refine procedure
	// Also print convergence information to screen for auto-refine procedure
	void checkConvergence(bool myverb = true);

	// Set metadata of a subset of particles to the experimental model
	void setMetaDataSubset(long int my_first_part_id, long int my_last_part_id);

	// Get metadata array of a subset of particles from the experimental model
	void getMetaAndImageDataSubset(long int my_first_part_id, long int my_last_part_id, bool do_also_imagedata = true);

	// Get the CTF (and Multiplicity weights where available) volumes from the stored files and correct them
	void get3DCTFAndMulti(MultidimArray<RFLOAT> &Ictf, MultidimArray<RFLOAT> &Fctf, MultidimArray<RFLOAT> &FstMulti,
			bool ctf_premultiplied);

	// Apply the Multiplicity weights to correct for Images and CTFs to properly estimate squared differences and averages
	void applySubtomoCorrection(MultidimArray<Complex > &Fimg, MultidimArray<Complex > &Fimg_nomask ,
								MultidimArray<RFLOAT> &Fctf, MultidimArray<RFLOAT> &FstMulti);

};

// Global call to threaded core of doThreadExpectationSomeParticles
void globalThreadExpectationSomeParticles(void *self, int thread_id);

#endif /* MAXLIK_H_ */

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

#include <pthread.h>
#include "src/ml_model.h"
#include "src/parallel.h"
#include "src/exp_model.h"
#include "src/ctf.h"
#include "src/time.h"
#include "src/mask.h"
#include "src/healpix_sampling.h"

#define ML_SIGNIFICANT_WEIGHT 1.e-8
#define METADATA_LINE_LENGTH METADATA_LINE_LENGTH_ALL

#define METADATA_LINE_LENGTH_NOCTF 11
#define METADATA_LINE_LENGTH_NOPRIOR 18
#define METADATA_LINE_LENGTH_NOTILT 26
#define METADATA_LINE_LENGTH_ALL 35

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
#define METADATA_CTF_VOLTAGE 14
#define METADATA_CTF_Q0 15
#define METADATA_CTF_CS 16
#define METADATA_CTF_BFAC 17

#define METADATA_ROT_PRIOR 18
#define METADATA_TILT_PRIOR 19
#define METADATA_PSI_PRIOR 20
#define METADATA_XOFF_PRIOR 21
#define METADATA_YOFF_PRIOR 22
#define METADATA_ZOFF_PRIOR 23

#define METADATA_BEAMTILT_X 24
#define METADATA_BEAMTILT_Y 25

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

// for profiling
//#define TIMING

class MlOptimiser;

class MlOptimiser
{
public:

	//TODO Fix void pointer
	void *cudaMlOptimiser; //CudaMlOptimiser

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

	// Filename for input tau2-spectrum
	FileName fn_tau;

	// Flag to keep tau-spectrum constant
	bool fix_tau;

    // some parameters for debugging
	double debug1, debug2;

	// Starting and finishing particles (for parallelisation)
    long int my_first_ori_particle_id, my_last_ori_particle_id;

	// Total number iterations and current iteration
	int iter, nr_iter;

	// Flag whether to split data from the beginning into two random halves
	bool do_split_random_halves;

	// resolution (in Angstrom) to join the two random halves
	double low_resol_join_halves;

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

	// Number of iterations without a resolution increase
	int nr_iter_wo_resol_gain;

	// Best resolution obtained thus far
	double best_resol_thus_far;

	// Is the FSC still high at the resolution limit?
	bool has_high_fsc_at_limit;

	// How many iterations ago was there a big step in the increase of current_size
	int has_large_incr_size_iter_ago;

	// In auto-sampling mode, use local searches from this sampling rate on
	int autosampling_hporder_local_searches;

	// Smallest changes thus far in the optimal translational offsets, orientations and classes
	double smallest_changes_optimal_offsets;
	double smallest_changes_optimal_orientations;
	int smallest_changes_optimal_classes;

	// Number of iterations without a decrease in OffsetChanges
	int nr_iter_wo_large_hidden_variable_changes;

	// Strict high-res limit in the expectation step
	double strict_highres_exp;

	// Flag to indicate to estimate angular accuracy until current_size (and not coarse_size) when restricting high-res limit in the expectation step
	// This is used for testing purposes only
	bool do_acc_currentsize_despite_highres_exp;

	// Global parameters to store accuracy on rot and trans
	double acc_rot, acc_trans;

	// Flag to indicate to use all data out to Nyquist
	bool do_use_all_data;

	// Flag to indicate the refinement has converged
	bool has_converged;

	// Flag to indicate that angular sampling in auto-sampling has reached its limit
	bool has_fine_enough_angular_sampling;

	// Flag to keep sigma2_offset fixed
	bool fix_sigma_offset;

	// Flag to keep sigma2_noise fixed
	bool fix_sigma_noise;

	//  Use images only up to a certain resolution in the expectation step
	int coarse_size;

	// Use images only up to a certain resolution in the expectation step
	int max_coarse_size;

	// Particle diameter (in Ang)
	double particle_diameter;

	// How many fourier shells should be included beyond the highest shell where evidenceVsPriorRatio < 1?
	int incr_size;

	// Minimum resolution to perform Bayesian estimates of the model
	int minres_map;

	// Flag to flatten solvent
	bool do_solvent;

	// Filename for a user-provided mask
	FileName fn_mask;

	// Filename for a user-provided second solvent mask
	// This solvent mask will have its own average density and may be useful for example to fill the interior of an icosahedral virus
	FileName fn_mask2;

	// Width of the soft-edges of the circular masks
	int width_mask_edge;

	// Number of particles to be processed simultaneously
	int nr_pool;

	// Available memory (in Gigabyte)
	double available_memory;

	// Perform combination of weight through files written on disc
	bool combine_weights_thru_disc;

	// Only print metadata label definitions and exit
	bool do_print_metadata_labels;

    // Use parallel access to disc?
    bool do_parallel_disc_io;

    // Use gpu resources?
    bool do_gpu;

	// Print the symmetry transformation matrices
	bool do_print_symmetry_ops;

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
	double adaptive_fraction;

	// Seed for random number generator
	int random_seed;

	/* Flag to indicate orientational (i.e. rotational AND translational) searches will be skipped */
	bool do_skip_align;

	/* Flag to indicate rotational searches will be skipped */
	bool do_skip_rotate;

	/* Flag to indicate maximization step will be skipped: only data.star file will be written out */
	bool do_skip_maximization;

	//////// Special stuff for the first iterations /////////////////

	// Skip marginalisation in first iteration and use signal cross-product instead of Gaussian
	bool do_firstiter_cc;

	/// Always perform cross-correlation instead of marginalization
	bool do_always_cc;

	// Initial low-pass filter for all references (in digital frequency)
	double ini_high;

	// Flag whether to generate seeds
	// TODO: implement!
	bool do_generate_seeds;

	// Flag whether to calculate the average of the unaligned images (for 2D refinements)
	bool do_average_unaligned;

	// Flag whether to calculate initial sigma_noise spectra
	bool do_calculate_initial_sigma_noise;

	// Flag to switch off error message about normalisation
	bool dont_raise_norm_error;

	///////// Re-align individual frames of movies /////////////

	// Flag whether to realign frames of movies
	bool do_realign_movies;

	// Starfile with the movie-frames
	FileName fn_data_movie;

	// How many individual frames contribute to the priors?
	int nr_frames_per_prior;

	// How wide are the running averages of the frames to use for alignment?
	// If set to 1, then running averages will be n-1, n, n+1, i.e. 3 frames wide
	int movie_frame_running_avg_side;

	///////// Hidden stuff, does not work with read/write: only via command-line ////////////////

	// Number of iterations for gridding preweighting reconstruction
	int gridding_nr_iter;

	// Flag whether to do group-wise B-factor correction or not
	bool do_bfactor;

	// Flag whether to use maximum a posteriori (MAP) estimation
	bool do_map;

	// For debugging/analysis: Name for initial sigma_noise sepctrum;
	FileName fn_sigma;

	// Multiplicative fdge factor for the sigma estimates
	double sigma2_fudge;

	// Perform two reconstructions of random halves sequentially (to save memory for very big cases)
	bool do_sequential_halves_recons;

	// Do zero-filled soft-masks instead of noisy masks on experimental particles
	// This will increase SNR but introduce correlations that are not modelled...
	// Until now the best refinements have used the noisy mask, not the soft mask....
	bool do_zero_mask;

	/////////// Keep track of hidden variable changes ////////////////////////

	// Changes from one iteration to the next in the angles
	double current_changes_optimal_orientations, sum_changes_optimal_orientations;

	// Changes from one iteration to the next in the translations
	double current_changes_optimal_offsets, sum_changes_optimal_offsets;

	// Changes from one iteration to the next in the class assignments
	double current_changes_optimal_classes, sum_changes_optimal_classes;

	// Just count how often the optimal changes are summed
	double sum_changes_count;

	/////////// Some internal stuff ////////////////////////

    // Array with pointers to the resolution of each point in a Fourier-space FFTW-like array
	MultidimArray<int> Mresol_fine, Mresol_coarse, Npix_per_shell;

	// Tabulated sin and cosine functions for shifts in Fourier space
	TabSine tab_sin;
	TabCosine tab_cos;

	// Verbosity flag
	int verb;

	// Thread Managers for the expectation step: one for all (pooled) particles
	ThreadTaskDistributor *exp_ipart_ThreadTaskDistributor;

	// Number of threads to run in parallel
	int nr_threads;

	long int exp_my_first_ori_particle, exp_my_last_ori_particle;
	MultidimArray<double> exp_metadata, exp_imagedata;
	std::string exp_fn_img, exp_fn_ctf, exp_fn_recimg;
	int exp_nr_images;

	// Calculate translated images on-the-fly
	bool do_shifts_onthefly;
	std::vector<MultidimArray<Complex> > global_fftshifts_ab_coarse, global_fftshifts_ab_current, global_fftshifts_ab2_coarse, global_fftshifts_ab2_current;

	/** Some global variables that are only for thread visibility */
	/// Taken from getAllSquaredDifferences
	/*
	std::vector<MultidimArray<Complex > > exp_Fimgs, exp_Fimgs_nomask, exp_local_Fimgs_shifted, exp_local_Fimgs_shifted_nomask;
	std::vector<MultidimArray<double> > exp_Fctfs, exp_local_Fctfs, exp_local_Minvsigma2s;
	Matrix2D<double> exp_R_mic;
	int exp_iseries, exp_iclass, exp_ipass, exp_iimage, exp_ipart, exp_current_image_size, exp_current_oversampling, exp_nr_ori_particles, exp_nr_particles, exp_nr_images;
	long int exp_nr_oversampled_rot, exp_nr_oversampled_trans, exp_nr_rot, exp_nr_dir, exp_nr_psi, exp_nr_trans;
	long int exp_part_id, exp_my_first_ori_particle, exp_my_last_ori_particle
	std::vector<int> exp_starting_image_no;
	std::vector<long int> exp_ipart_to_part_id, exp_ipart_to_ori_part_id, exp_ipart_to_ori_part_nframe, exp_iimg_to_ipart;
	std::vector<double> exp_highres_Xi2_imgs, exp_min_diff2, exp_local_sqrtXi2, exp_local_oldcc;
	MultidimArray<double> exp_Mweight;
	MultidimArray<bool> exp_Mcoarse_significant;
	// And from storeWeightedSums
	std::vector<double> exp_sum_weight, exp_significant_weight, exp_max_weight;
	std::vector<Matrix1D<double> > exp_old_offset, exp_prior;
	std::vector<double> exp_wsum_norm_correction;
	std::vector<MultidimArray<double> > exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA, exp_power_imgs;
	double exp_thisparticle_sumweight;
	*/
	//TMP DEBUGGING
	MultidimArray<double> DEBUGGING_COPY_exp_Mweight;

#ifdef TIMING
    Timer timer;
	int TIMING_DIFF_PROJ, TIMING_DIFF_SHIFT, TIMING_DIFF_DIFF2;
	int TIMING_WSUM_PROJ, TIMING_WSUM_BACKPROJ, TIMING_WSUM_DIFF2, TIMING_WSUM_SUMSHIFT;
	int TIMING_EXP, TIMING_MAX, TIMING_RECONS;
	int TIMING_ESP, TIMING_ESP_THR, TIMING_ESP_ONEPART, TIMING_ESP_ONEPARTN, TIMING_EXP_METADATA, TIMING_EXP_CHANGES;
	int TIMING_ESP_FT, TIMING_ESP_INI, TIMING_ESP_DIFF1, TIMING_ESP_DIFF2;
	int TIMING_ESP_DIFF2_A, TIMING_ESP_DIFF2_B, TIMING_ESP_DIFF2_C, TIMING_ESP_DIFF2_D, TIMING_ESP_DIFF2_E;
	int TIMING_ESP_PREC1, TIMING_ESP_PREC2, TIMING_ESP_PRECW, TIMING_WSUM_GETSHIFT, TIMING_DIFF2_GETSHIFT, TIMING_WSUM_SCALE, TIMING_WSUM_LOCALSUMS;
	int TIMING_ESP_WEIGHT1, TIMING_ESP_WEIGHT2, TIMING_WEIGHT_EXP, TIMING_WEIGHT_SORT, TIMING_ESP_WSUM;
	int TIMING_EXTRA1, TIMING_EXTRA2, TIMING_EXTRA3;
#endif

public:

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
	void read(FileName fn_in, int rank = 0);

	// Write files to disc
	void write(bool do_write_sampling, bool do_write_data, bool do_write_optimiser, bool do_write_model, int random_subset = 0);

    /** ========================== Initialisation  =========================== */

	// Initialise the whole optimiser
	void initialise();

	// Some general stuff that is shared between MPI and sequential code
	void initialiseGeneral(int rank = 0);

	// Randomise particle processing order and resize metadata array
	void initialiseWorkLoad();

	/* Calculates the sum of all individual power spectra and the average of all images for initial sigma_noise estimation
	 * The rank is passed so that if one splits the data into random halves one can know which random half to treat
	 */
	void calculateSumOfPowerSpectraAndAverageImage(MultidimArray<double> &Mavg, bool myverb = true);

	/** Use the sum of the individual power spectra to calculate their average and set this in sigma2_noise
	 * Also subtract the power spectrum of the average images,
	 * and if (do_average_unaligned) then also set Mavg to all Iref
	 */
	void setSigmaNoiseEstimatesAndSetAverageImage(MultidimArray<double> &Mavg);

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

	/* Setup expectation step */
	void expectationSetup();

	/* Check whether everything fits into memory, possibly adjust nr_pool and setup thread task managers */
	void expectationSetupCheckMemory(bool myverb = true);

	/* For on-the-fly shifts, precalculates AB-matrices */
	void precalculateABMatrices();

	/* Perform the expectation integration over all k, phi and series elements for a number (some) of pooled particles
	 * The number of pooled particles is determined by max_nr_pool and some memory checks in expectationSetup()
	 */
	void expectationSomeParticles(long int my_first_particle, long int my_last_particle);

	/* Perform expectation step for some particles using threads */
	void doThreadExpectationSomeParticles(int thread_id);

	/* Perform the expectation integration over all k, phi and series elements for a given particle */
	void expectationOneParticle(long int my_ori_particle, int thread_id);

	/* Maximization step
	 * Updates the current model: reconstructs and updates all other model parameter
	 */
	void maximization();

	/* Perform the actual reconstructions
	 * This is officially part of the maximization, but it is separated because of parallelisation issues.
	 */
	void maximizationReconstructClass(int iclass);

	/* Updates all other model parameters (besides the reconstructions)
	 */
	void maximizationOtherParameters();


    /** ========================== Deeper functions ================================= */

	/* Apply a solvent flattening to a map
	 */
	void solventFlatten();

	/* Updates the current resolution (from data_vs_prior array) and keeps track of best resolution thus far
	 *  and precalculates a 2D Fourier-space array with pointers to the resolution of each point in a FFTW-centered array
	 */
	void updateCurrentResolution();

	/* Update the current and coarse image size
	 *  and precalculates a 2D Fourier-space array with pointers to the resolution of each point in a FFTW-centered array
	 */
	void updateImageSizeAndResolutionPointers();

	/* From the vectors of Fourier transforms of the images, calculate running averages over the movie frames
	 */
	void calculateRunningAveragesOfMovieFrames(long int my_ori_particle,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<double> > &exp_power_imgs,
		std::vector<double> &exp_highres_Xi2_imgs);

	/* Read image and its metadata from disc (threaded over all pooled particles)
	 */
	void getFourierTransformsAndCtfs(long int my_ori_particle, int metadata_offset,
			std::vector<MultidimArray<Complex > > &exp_Fimgs,
			std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
			std::vector<MultidimArray<double> > &exp_Fctfs,
			std::vector<Matrix1D<double> > &exp_old_offset,
			std::vector<Matrix1D<double> > &exp_prior,
			std::vector<MultidimArray<double> > &exp_power_imgs,
			std::vector<double> &exp_highres_Xi2_imgs,
			std::vector<int> &exp_pointer_dir_nonzeroprior,
			std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<double> &exp_directions_prior,
			std::vector<double> &exp_psi_prior);

	/* Store all shifted FourierTransforms in a vector
	 * also store precalculated 2D matrices with 1/sigma2_noise
	 */
	void precalculateShiftedImagesCtfsAndInvSigma2s(bool do_also_unmasked, long int my_ori_particle,
			int exp_current_image_size, int exp_current_oversampling,
			int exp_itrans_min, int exp_itrans_max,
			std::vector<MultidimArray<Complex > > &exp_Fimgs,
			std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
			std::vector<MultidimArray<double> > &exp_Fctfs,
			std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
			std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
			std::vector<MultidimArray<double> > &exp_local_Fctfs,
			std::vector<double> &exp_local_sqrtXi2,
			std::vector<MultidimArray<double> > &exp_local_Minvsigma2s);

	// Given exp_Mcoarse_significant, check for iorient whether any of the particles has any significant (coarsely sampled) translation
	bool isSignificantAnyParticleAnyTranslation(long int iorient,
			int exp_itrans_min, int exp_itrans_max, MultidimArray<bool> &exp_Mcoarse_significant);

	// Get squared differences for all iclass, idir, ipsi and itrans...
	void getAllSquaredDifferences(long int my_ori_particle, int exp_current_image_size,
			int exp_ipass, int exp_current_oversampling, int metadata_offset,
			int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
			int exp_itrans_min, int exp_itrans_max, int my_iclass_min, int my_iclass_max,
			std::vector<double> &exp_min_diff2,
			std::vector<double> &exp_highres_Xi2_imgs,
			std::vector<MultidimArray<Complex > > &exp_Fimgs,
			std::vector<MultidimArray<double> > &exp_Fctfs,
			MultidimArray<double> &exp_Mweight,
			MultidimArray<bool> &exp_Mcoarse_significant,
			std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<double> &exp_directions_prior, std::vector<double> &exp_psi_prior,
			std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
			std::vector<MultidimArray<double> > &exp_local_Minvsigma2s,
			std::vector<MultidimArray<double> > &exp_local_Fctfs,
			std::vector<double> &exp_local_sqrtXi2);

	// Convert all squared difference terms to weights.
	// Also calculates exp_sum_weight and, for adaptive approach, also exp_significant_weight
	void convertAllSquaredDifferencesToWeights(long int my_ori_particle, int exp_ipass,
			int exp_current_oversampling, int metadata_offset,
			int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
			int exp_itrans_min, int exp_itrans_max, int my_iclass_min, int my_iclass_max,
			MultidimArray<double> &exp_Mweight, MultidimArray<bool> &exp_Mcoarse_significant,
			std::vector<double> &exp_significant_weight, std::vector<double> &exp_sum_weight,
			std::vector<Matrix1D<double> > &exp_old_offset, std::vector<Matrix1D<double> > &exp_prior,
			std::vector<double> &exp_min_diff2,
			std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<double> &exp_directions_prior, std::vector<double> &exp_psi_prior);

	// Store all relevant weighted sums, also return optimal hidden variables, max_weight and dLL
	void storeWeightedSums(long int my_ori_particle, int exp_current_image_size,
			int exp_current_oversampling, int metadata_offset,
			int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
			int exp_itrans_min, int exp_itrans_max, int my_iclass_min, int my_iclass_max,
			std::vector<double> &exp_min_diff2,
			std::vector<double> &exp_highres_Xi2_imgs,
			std::vector<MultidimArray<Complex > > &exp_Fimgs,
			std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
			std::vector<MultidimArray<double> > &exp_Fctfs,
			std::vector<MultidimArray<double> > &exp_power_imgs,
			std::vector<Matrix1D<double> > &exp_old_offset,
			std::vector<Matrix1D<double> > &exp_prior,
			MultidimArray<double> &exp_Mweight,
			MultidimArray<bool> &exp_Mcoarse_significant,
			std::vector<double> &exp_significant_weight,
			std::vector<double> &exp_sum_weight,
			std::vector<double> &exp_max_weight,
			std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
			std::vector<double> &exp_directions_prior, std::vector<double> &exp_psi_prior,
			std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
			std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
			std::vector<MultidimArray<double> > &exp_local_Minvsigma2s,
			std::vector<MultidimArray<double> > &exp_local_Fctfs,
			std::vector<double> &exp_local_sqrtXi2);

	/** Monitor the changes in the optimal translations, orientations and class assignments for some particles */
	void monitorHiddenVariableChanges(long int my_first_ori_particle, long int my_last_ori_particle);

	// Updates the overall changes in the hidden variables and keeps track of nr_iter_wo_large_changes_in_hidden_variables
	void updateOverallChangesInHiddenVariables();

	// Calculate expected error in orientational assignments
	// Based on comparing projections of the model and see how many degrees apart gives rise to difference of power > 3*sigma^ of the noise
	void calculateExpectedAngularErrors(long int my_first_ori_particle, long int my_last_ori_particle);

	// Adjust angular sampling based on the expected angular accuracies for auto-refine procedure
	void updateAngularSampling(bool verb = true);

	// Check convergence for auto-refine procedure
	void checkConvergence();

	// Print convergence information to screen for auto-refine procedure
	void printConvergenceStats();

	// Set metadata of a subset of particles to the experimental model
	void setMetaDataSubset(int first_ori_particle_id, int last_ori_particle_id);

	// Get metadata array of a subset of particles from the experimental model
	void getMetaAndImageDataSubset(int first_ori_particle_id, int last_ori_particle_id, bool do_also_imagedata = true);

};

// Global call to threaded core of doThreadExpectationSomeParticles
void globalThreadExpectationSomeParticles(ThreadArgument &thArg);

#endif /* MAXLIK_H_ */

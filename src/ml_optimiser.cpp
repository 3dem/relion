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

//#define DEBUG
//#define DEBUG_CHECKSIZES
//#define DEBUG_BODIES
//#define DEBUG_HELICAL_ORIENTATIONAL_SEARCH
//#define PRINT_GPU_MEM_INFO

#ifdef TIMING
		#define RCTIC(timer,label) (timer.tic(label))
		#define RCTOC(timer,label) (timer.toc(label))
#else
		#define RCTIC(timer,label)
    	#define RCTOC(timer,label)
#endif

#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <string>
#include <fstream>
#include "src/macros.h"
#include "src/error.h"
#include "src/ml_optimiser.h"
#ifdef CUDA
#include "src/gpu_utils/cuda_ml_optimiser.h"
#endif

#define NR_CLASS_MUTEXES 5

//Some global threads management variables
static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;
Barrier * global_barrier;
ThreadManager * global_ThreadManager;

/** ========================== Threaded parallelization of expectation === */

void globalThreadExpectationSomeParticles(ThreadArgument &thArg)
{
	MlOptimiser *MLO = (MlOptimiser*) thArg.workClass;

	try
	{
#ifdef CUDA
		if (MLO->do_gpu)
			((MlOptimiserCuda*) MLO->cudaOptimisers[thArg.thread_id])->doThreadExpectationSomeParticles(thArg.thread_id);
		else
#endif
		MLO->doThreadExpectationSomeParticles(thArg.thread_id);
	}
	catch (RelionError XE)
	{
		RelionError *gE = new RelionError(XE.msg, XE.file, XE.line);
		gE->msg = XE.msg;
		MLO->threadException = gE;
	}
}


/** ========================== I/O operations  =========================== */


void MlOptimiser::usage()
{
	parser.writeUsage(std::cout);
}

void MlOptimiser::read(int argc, char **argv, int rank)
{
//#define DEBUG_READ

	parser.setCommandLine(argc, argv);

	if (checkParameter(argc, argv, "--continue"))
	{
		// Do this before reading in the data.star file below!
		do_preread_images   = checkParameter(argc, argv, "--preread_images");
		do_parallel_disc_io = !checkParameter(argc, argv, "--no_parallel_disc_io");

		parser.addSection("Continue options");
		FileName fn_in = parser.getOption("--continue", "_optimiser.star file of the iteration after which to continue");
		// Read in previously calculated parameters
		if (fn_in != "")
			read(fn_in, rank);
		// And look for additional command-line options...
		parseContinue(argc, argv);
	}
	else
	{
		// Start a new run from scratch
		parseInitial(argc, argv);
	}
}

void MlOptimiser::parseContinue(int argc, char **argv)
{
#ifdef DEBUG
	std::cerr << "Entering parseContinue" << std::endl;
#endif

	int general_section = parser.addSection("General options");
	// Not all parameters are accessible here...
	FileName fn_out_new = parser.getOption("--o", "Output rootname", "OLD_ctX");
	if (fn_out_new == "OLD_ctX" || fn_out_new == fn_out )
		fn_out += "_ct" + integerToString(iter);
	else
		fn_out = fn_out_new;

	do_force_converge =  parser.checkOption("--force_converge", "Force an auto-refinement run to converge immediately upon continuation.");

	std::string fnt;
	fnt = parser.getOption("--iter", "Maximum number of iterations to perform", "OLD");
	if (fnt != "OLD")
		nr_iter = textToInteger(fnt);

	fnt = parser.getOption("--tau2_fudge", "Regularisation parameter (values higher than 1 give more weight to the data)", "OLD");
	if (fnt != "OLD")
		mymodel.tau2_fudge_factor = textToFloat(fnt);

	// Solvent flattening
	if (parser.checkOption("--flatten_solvent", "Switch on masking on the references?", "OLD"))
		do_solvent = true;

	// Check whether the mask has changed
	fnt = parser.getOption("--solvent_mask", "User-provided mask for the references", "OLD");
	if (fnt != "OLD")
		fn_mask = fnt;

	// Check whether the secondary mask has changed
	fnt = parser.getOption("--solvent_mask2", "User-provided secondary mask", "OLD");
	if (fnt != "OLD")
		fn_mask2 = fnt;

	// Check whether tau2-spectrum has changed
	fnt = parser.getOption("--tau", "STAR file with input tau2-spectrum (to be kept constant)", "OLD");
	if (fnt != "OLD")
		fn_tau = fnt;

	// Check whether particle diameter has changed
	fnt = parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms)", "OLD");
	if (fnt != "OLD")
		particle_diameter = textToFloat(fnt);

	// SGD stuff
	do_sgd = parser.checkOption("--sgd", "Perform stochastic gradient descent instead of default expectation-maximization");

	if (parser.checkOption("--stop_sgd", "Switch from stochastic gradient descent to expectation maximisation for the remainder of iterations."))
		do_sgd = false;

	fnt = parser.getOption("--max_subsets", "Stop SGD after processing this many subsets (possibly more than 1 iteration)", "OLD");
	if (fnt != "OLD")
		sgd_max_subsets = textToInteger(fnt);

	fnt = parser.getOption("--subset_size", "Size of the subsets for SGD", "OLD");
	if (fnt != "OLD")
		subset_size = textToInteger(fnt);

	fnt = parser.getOption("--sgd_stepsize", "Step size parameter for SGD updates", "OLD");
	if (fnt != "OLD")
		sgd_stepsize = textToInteger(fnt);

	fnt = parser.getOption("--write_subsets", "Write out model every so many subsets", "OLD");
	if (fnt != "OLD")
		write_every_subset = textToInteger(fnt);

	fnt = parser.getOption("--mu", "Momentum parameter for SGD updates", "OLD");
	if (fnt != "OLD")
		mu = textToFloat(fnt);


	do_join_random_halves = parser.checkOption("--join_random_halves", "Join previously split random halves again (typically to perform a final reconstruction).");

	// Re-align movie frames
	int movie_section = parser.addSection("Re-align movie frames");

	fn_data_movie = parser.getOption("--realign_movie_frames", "Input STAR file with the movie frames", "");

	do_movies_in_batches = parser.checkOption("--process_movies_in_batches", "Perform movie-processing one micrograph at a time (this reduces RAM costs for large data sets)");
	movie_identifier = parser.getOption("--movie_name", "Identified string for all movie files, e.g. (_movie.mrcs)", "movie");

	// TODO: add this to EMDL_OPTIMISER and read/write of optimiser.star
	nr_frames_per_prior = textToInteger(parser.getOption("--nr_frames_prior", "Number of movie frames to calculate running-average priors", "5"));

	// (integer-) divide running average width by 2 to have the side only
	movie_frame_running_avg_side = textToInteger(parser.getOption("--movie_frames_running_avg", "Number of movie frames in each running average", "3"));
	movie_frame_running_avg_side /= 2;

	// ORIENTATIONS
	int orientations_section = parser.addSection("Orientations");

	fnt = parser.getOption("--oversampling", "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)", "OLD");
	if (fnt != "OLD")
		adaptive_oversampling = textToInteger(fnt);

	// Check whether angular sampling has changed
	// Do not do this for auto_refine, but make sure to do this when realigning movies!
	if (!do_auto_refine || fn_data_movie != "")
	{
		directions_have_changed = false;
		fnt = parser.getOption("--healpix_order", "Healpix order for the angular sampling rate on the sphere (before oversampling): hp2=15deg, hp3=7.5deg, etc", "OLD");
		if (fnt != "OLD")
		{
			int _order = textToInteger(fnt);
			if (_order != sampling.healpix_order)
			{
				directions_have_changed = true;
				sampling.healpix_order = _order;
			}
		}

		fnt = parser.getOption("--psi_step", "Angular sampling (before oversampling) for the in-plane angle (default=10deg for 2D, hp sampling for 3D)", "OLD");
		if (fnt != "OLD")
			sampling.psi_step = textToFloat(fnt);

		fnt = parser.getOption("--offset_range", "Search range for origin offsets (in pixels)", "OLD");
		if (fnt != "OLD")
			sampling.offset_range = textToFloat(fnt);

		fnt = parser.getOption("--offset_step", "Sampling rate for origin offsets (in pixels)", "OLD");
		if (fnt != "OLD")
			sampling.offset_step = textToFloat(fnt);
	}

	fnt = parser.getOption("--auto_local_healpix_order", "Minimum healpix order (before oversampling) from which auto-refine procedure will use local searches", "OLD");
	if (fnt != "OLD")
		autosampling_hporder_local_searches = textToInteger(fnt);

	// Check whether the prior mode changes
	RFLOAT _sigma_rot, _sigma_tilt, _sigma_psi, _sigma_off;
	int _mode;
	fnt = parser.getOption("--sigma_ang", "Stddev on all three Euler angles for local angular searches (of +/- 3 stddev)", "OLD");
	if (fnt != "OLD")
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--sigma_rot", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "OLD");
	if (fnt != "OLD")
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_rot = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--sigma_tilt", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "OLD");
	if (fnt != "OLD")
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_tilt = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--sigma_psi", "Stddev on the in-plane angle for local angular searches (of +/- 3 stddev)", "OLD");
	if (fnt != "OLD")
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_psi = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--sigma_off", "Stddev. on the translations", "OLD");
	if (fnt != "OLD")
	{
		mymodel.sigma2_offset = textToFloat(fnt) * textToFloat(fnt);
	}
	fnt = parser.getOption("--helical_inner_diameter", "Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)", "OLD");
	if (fnt != "OLD")
	{
		helical_tube_inner_diameter = textToFloat(fnt);
	}
	fnt = parser.getOption("--helical_outer_diameter", "Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)", "OLD");
	if (fnt != "OLD")
	{
		helical_tube_outer_diameter = textToFloat(fnt);
	}

	if (parser.checkOption("--skip_align", "Skip orientational assignment (only classify)?"))
		do_skip_align = true;
	else
		do_skip_align = false; // do_skip_align should normally be false...

	if (parser.checkOption("--skip_rotate", "Skip rotational assignment (only translate and classify)?"))
		do_skip_rotate = true;
	else
		do_skip_rotate = false; // do_skip_rotate should normally be false...

	if (parser.checkOption("--bimodal_psi", "Do bimodal searches of psi angle?")) // Oct07,2015 - Shaoda, bimodal psi
		do_bimodal_psi = true;
	else
		do_bimodal_psi = false;



	do_skip_maximization = parser.checkOption("--skip_maximize", "Skip maximization step (only write out data.star file)?");

	int corrections_section = parser.addSection("Corrections");

	// Can also switch the following option OFF
	if (parser.checkOption("--scale", "Switch on intensity-scale corrections on image groups", "OLD"))
		do_scale_correction = true;
	if (parser.checkOption("--no_scale", "Switch off intensity-scale corrections on image groups", "OLD"))
		do_scale_correction = false;

	// Can also switch the following option OFF
	if (parser.checkOption("--norm", "Switch on normalisation-error correction","OLD"))
		do_norm_correction = true;
	if (parser.checkOption("--no_norm", "Switch off normalisation-error correction","OLD"))
		do_norm_correction = false;

	int computation_section = parser.addSection("Computation");

	x_pool = textToInteger(parser.getOption("--pool", "Number of images to pool for each thread task", "1"));
	nr_threads = textToInteger(parser.getOption("--j", "Number of threads to run in parallel (only useful on multi-core machines)", "1"));
	do_parallel_disc_io = !parser.checkOption("--no_parallel_disc_io", "Do NOT let parallel (MPI) processes access the disc simultaneously (use this option with NFS)");
	combine_weights_thru_disc = !parser.checkOption("--dont_combine_weights_via_disc", "Send the large arrays of summed weights through the MPI network, instead of writing large files to disc");
	do_shifts_onthefly = parser.checkOption("--onthefly_shifts", "Calculate shifted images on-the-fly, do not store precalculated ones in memory");
	do_preread_images  = parser.checkOption("--preread_images", "Use this to let the master process read all particles into memory. Be careful you have enough RAM for large data sets!");
	fn_scratch = parser.getOption("--scratch_dir", "If provided, particle stacks will be copied to this local scratch disk prior to refinement.", "");
	keep_free_scratch_Gb = textToInteger(parser.getOption("--keep_free_scratch", "Space available for copying particle stacks (in Gb)", "10"));
	do_reuse_scratch = parser.checkOption("--reuse_scratch", "Re-use data on scratchdir, instead of wiping it and re-copying all data.");


	failsafe_threshold = textToInteger(parser.getOption("--failsafe_threshold", "Maximum number of particles permitted to be drop, due to zero sum of weights, before exiting with an error (GPU only).", "40"));

	do_gpu = parser.checkOption("--gpu", "Use available gpu resources for some calculations");
	gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread","default");
#ifndef CUDA
	if(do_gpu)
	{
		std::cerr << "+ WARNING : Relion was compiled without CUDA of at least version 7.0 - you do NOT have support for GPUs" << std::endl;
		do_gpu = false;
	}
#endif
	double temp_reqSize = textToDouble(parser.getOption("--free_gpu_memory", "GPU device memory (in Mb) to leave free after allocation.", "0"));
	temp_reqSize *= 1000*1000;
	if(temp_reqSize<0)
		REPORT_ERROR("Invalid free_gpu_memory value.");
	else
		requested_free_gpu_memory =  temp_reqSize;

	do_phase_random_fsc = parser.checkOption("--solvent_correct_fsc", "Correct FSC curve for the effects of the solvent mask?");
	verb = textToInteger(parser.getOption("--verb", "Verbosity (1=normal, 0=silent)", "1"));

	int expert_section = parser.addSection("Expert options");

	fnt = parser.getOption("--strict_highres_exp", "Resolution limit (in Angstrom) to restrict probability calculations in the expectation step", "OLD");
	if (fnt != "OLD")
		strict_highres_exp = textToFloat(fnt);

	fnt = parser.getOption("--strict_highres_sgd", "Resolution limit (in Angstrom) to restrict probability calculations in SGD", "OLD");
	if (fnt != "OLD")
		strict_highres_sgd = textToFloat(fnt);

	// Debugging/analysis/hidden stuff
	do_map = !checkParameter(argc, argv, "--no_map");
	minres_map = textToInteger(getParameter(argc, argv, "--minres_map", "5"));
    gridding_nr_iter = textToInteger(getParameter(argc, argv, "--gridding_iter", "10"));
	debug1 = textToFloat(getParameter(argc, argv, "--debug1", "0."));
	debug2 = textToFloat(getParameter(argc, argv, "--debug2", "0."));
	debug3 = textToFloat(getParameter(argc, argv, "--debug3", "0."));
    do_bfactor = checkParameter(argc, argv, "--bfactor");
	// Read in initial sigmaNoise spectrum
	fn_sigma = getParameter(argc, argv, "--sigma","");
	sigma2_fudge = textToFloat(getParameter(argc, argv, "--sigma2_fudge", "1."));
	do_acc_currentsize_despite_highres_exp = checkParameter(argc, argv, "--accuracy_current_size");
	do_sequential_halves_recons  = checkParameter(argc, argv, "--sequential_halves_recons");
	do_always_join_random_halves = checkParameter(argc, argv, "--always_join_random_halves");
	do_use_all_data = checkParameter(argc, argv, "--use_all_data");
	do_always_cc  = checkParameter(argc, argv, "--always_cc");
	do_only_sample_tilt  = checkParameter(argc, argv, "--only_sample_tilt");
	minimum_angular_sampling = textToFloat(getParameter(argc, argv, "--minimum_angular_sampling", "0"));
	asymmetric_padding = parser.checkOption("--asymmetric_padding", "", "false", true);
	maximum_significants = textToInteger(parser.getOption("--maxsig", "", "0", true));
	skip_gridding = parser.checkOption("--skip_gridding", "", "false", true);

	do_print_metadata_labels = false;
	do_print_symmetry_ops = false;
#ifdef DEBUG
	std::cerr << "Leaving parseContinue" << std::endl;
#endif

}

void MlOptimiser::parseInitial(int argc, char **argv)
{
#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::parseInitial Entering "<<std::endl;
#endif

	// Read/initialise mymodel and sampling from a STAR file
    FileName fn_model = getParameter(argc, argv, "--model", "None");
	if (fn_model != "None")
	{
		mymodel.read(fn_model);
	}
	// Read in the sampling information from a _sampling.star file
    FileName fn_sampling = getParameter(argc, argv, "--sampling", "None");
	if (fn_sampling != "None")
	{
		sampling.read(fn_sampling);
	}

	// General optimiser I/O stuff
    int general_section = parser.addSection("General options");
    fn_data = parser.getOption("--i", "Input images (in a star-file or a stack)");
    fn_out = parser.getOption("--o", "Output rootname");
    nr_iter = textToInteger(parser.getOption("--iter", "Maximum number of iterations to perform", "50"));
    mymodel.pixel_size = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "-1"));
	mymodel.tau2_fudge_factor = textToFloat(parser.getOption("--tau2_fudge", "Regularisation parameter (values higher than 1 give more weight to the data)", "1"));
	mymodel.nr_classes = textToInteger(parser.getOption("--K", "Number of references to be refined", "1"));
    particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms)", "-1"));
	do_zero_mask = parser.checkOption("--zero_mask","Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)");
	do_solvent = parser.checkOption("--flatten_solvent", "Perform masking on the references as well?");
	fn_mask = parser.getOption("--solvent_mask", "User-provided mask for the references (default is to use spherical mask with particle_diameter)", "None");
	fn_mask2 = parser.getOption("--solvent_mask2", "User-provided secondary mask (with its own average density)", "None");
        //fn_body_masks = parser.getOption("--multibody_masks", "STAR file with binary masks for multi-body refinement", "");
        fn_body_masks = "";
    fn_tau = parser.getOption("--tau", "STAR file with input tau2-spectrum (to be kept constant)", "None");
	fn_local_symmetry = parser.getOption("--local_symmetry", "Local symmetry description file containing list of masks and their operators", "None");
    do_split_random_halves = parser.checkOption("--split_random_halves", "Refine two random halves of the data completely separately");
	low_resol_join_halves = textToFloat(parser.getOption("--low_resol_join_halves", "Resolution (in Angstrom) up to which the two random half-reconstructions will not be independent to prevent diverging orientations","-1"));

	// Initialisation
	int init_section = parser.addSection("Initialisation");
	fn_ref = parser.getOption("--ref", "Image, stack or star-file with the reference(s). (Compulsory for 3D refinement!)", "None");
	is_3d_model = parser.checkOption("--denovo_3dref", "Make an initial 3D model from randomly oriented 2D particles");
	mymodel.sigma2_offset = textToFloat(parser.getOption("--offset", "Initial estimated stddev for the origin offsets", "3"));
	mymodel.sigma2_offset *= mymodel.sigma2_offset;

	// Perform cross-product comparison at first iteration
	do_firstiter_cc = parser.checkOption("--firstiter_cc", "Perform CC-calculation in the first iteration (use this if references are not on the absolute intensity scale)");
	ini_high = textToFloat(parser.getOption("--ini_high", "Resolution (in Angstroms) to which to limit refinement in the first iteration ", "-1"));

	// Set the orientations
    int orientations_section = parser.addSection("Orientations");
	// Move these to sampling
	adaptive_oversampling = textToInteger(parser.getOption("--oversampling", "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)", "1"));
	sampling.healpix_order = textToInteger(parser.getOption("--healpix_order", "Healpix order for the angular sampling (before oversampling) on the (3D) sphere: hp2=15deg, hp3=7.5deg, etc", "2"));
	sampling.psi_step = textToFloat(parser.getOption("--psi_step", "Sampling rate (before oversampling) for the in-plane angle (default=10deg for 2D, hp sampling for 3D)", "-1"));
	sampling.limit_tilt = textToFloat(parser.getOption("--limit_tilt", "Limited tilt angle: positive for keeping side views, negative for keeping top views", "-91"));
	sampling.fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
	sampling.offset_range = textToFloat(parser.getOption("--offset_range", "Search range for origin offsets (in pixels)", "6"));
	sampling.offset_step = textToFloat(parser.getOption("--offset_step", "Sampling rate (before oversampling) for origin offsets (in pixels)", "2"));
	// Jun19,2015 - Shaoda, Helical refinement
	sampling.helical_offset_step = textToFloat(parser.getOption("--helical_offset_step", "Sampling rate (before oversampling) for offsets along helical axis (in pixels)", "-1"));
	sampling.perturbation_factor = textToFloat(parser.getOption("--perturb", "Perturbation factor for the angular sampling (0=no perturb; 0.5=perturb)", "0.5"));
	do_auto_refine = parser.checkOption("--auto_refine", "Perform 3D auto-refine procedure?");
	autosampling_hporder_local_searches = textToInteger(parser.getOption("--auto_local_healpix_order", "Minimum healpix order (before oversampling) from which autosampling procedure will use local searches", "4"));
	parser.setSection(orientations_section);
	RFLOAT _sigma_ang = textToFloat(parser.getOption("--sigma_ang", "Stddev on all three Euler angles for local angular searches (of +/- 3 stddev)", "-1"));
	RFLOAT _sigma_rot = textToFloat(parser.getOption("--sigma_rot", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "-1"));
	RFLOAT _sigma_tilt = textToFloat(parser.getOption("--sigma_tilt", "Stddev on the second Euler angle for local angular searches (of +/- 3 stddev)", "-1"));
	RFLOAT _sigma_psi = textToFloat(parser.getOption("--sigma_psi", "Stddev on the in-plane angle for local angular searches (of +/- 3 stddev)", "-1"));
	if (_sigma_ang > 0.)
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		// the sigma-values for the orientational prior are in model (and not in sampling) because one might like to estimate them
		// from the data by calculating weighted sums of all angular differences: therefore it needs to be in wsum_model and thus in mymodel.
		mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = _sigma_ang * _sigma_ang;
	}
	else if (_sigma_rot > 0. || _sigma_tilt > 0. || _sigma_psi > 0.)
	{
		mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
		mymodel.sigma2_rot  = (_sigma_rot > 0. ) ? _sigma_rot * _sigma_rot   : 0.;
		mymodel.sigma2_tilt = (_sigma_tilt > 0.) ? _sigma_tilt * _sigma_tilt : 0.;
		mymodel.sigma2_psi  = (_sigma_psi > 0. ) ? _sigma_psi * _sigma_psi   : 0.;
	}
	else
	{
		//default
		mymodel.orientational_prior_mode = NOPRIOR;
		mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 0.;
	}
	do_skip_align = parser.checkOption("--skip_align", "Skip orientational assignment (only classify)?");
	do_skip_rotate = parser.checkOption("--skip_rotate", "Skip rotational assignment (only translate and classify)?");
	do_bimodal_psi = parser.checkOption("--bimodal_psi", "Do bimodal searches of psi angle?"); // Oct07,2015 - Shaoda, bimodal psi
	do_skip_maximization = false;

	// Helical reconstruction
    int helical_section = parser.addSection("Helical reconstruction (in development...)");
    do_helical_refine = parser.checkOption("--helix", "Perform 3D classification or refinement for helices?");
    ignore_helical_symmetry = parser.checkOption("--ignore_helical_symmetry", "Ignore helical symmetry?");
    mymodel.helical_nr_asu = textToInteger(parser.getOption("--helical_nr_asu", "Number of new helical asymmetric units (asu) per box (1 means no helical symmetry is present)", "1"));
    helical_twist_initial = textToFloat(parser.getOption("--helical_twist_initial", "Helical twist (in degrees, positive values for right-handedness)", "0."));
    mymodel.helical_twist_min = textToFloat(parser.getOption("--helical_twist_min", "Minimum helical twist (in degrees, positive values for right-handedness)", "0."));
    mymodel.helical_twist_max = textToFloat(parser.getOption("--helical_twist_max", "Maximum helical twist (in degrees, positive values for right-handedness)", "0."));
    mymodel.helical_twist_inistep = textToFloat(parser.getOption("--helical_twist_inistep", "Initial step of helical twist search (in degrees)", "0."));
    helical_rise_initial = textToFloat(parser.getOption("--helical_rise_initial", "Helical rise (in Angstroms)", "0."));
    mymodel.helical_rise_min = textToFloat(parser.getOption("--helical_rise_min", "Minimum helical rise (in Angstroms)", "0."));
    mymodel.helical_rise_max = textToFloat(parser.getOption("--helical_rise_max", "Maximum helical rise (in Angstroms)", "0."));
    mymodel.helical_rise_inistep = textToFloat(parser.getOption("--helical_rise_inistep", "Initial step of helical rise search (in Angstroms)", "0."));
    helical_z_percentage = textToFloat(parser.getOption("--helical_z_percentage", "This box length along the center of Z axis contains good information of the helix. Important in imposing and refining symmetry", "0.3"));
    helical_tube_inner_diameter = textToFloat(parser.getOption("--helical_inner_diameter", "Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)", "-1."));
    helical_tube_outer_diameter = textToFloat(parser.getOption("--helical_outer_diameter", "Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)", "-1."));
    do_helical_symmetry_local_refinement = parser.checkOption("--helical_symmetry_search", "Perform local refinement of helical symmetry?");
    helical_sigma_distance = textToFloat(parser.getOption("--helical_sigma_distance", "Sigma of distance along the helical tracks", "-1."));
    helical_keep_tilt_prior_fixed = parser.checkOption("--helical_keep_tilt_prior_fixed", "Keep helical tilt priors fixed (at 90 degrees) in global angular searches?");
    if (ignore_helical_symmetry)
    {
    	mymodel.helical_nr_asu = 1; // IMPORTANT !
    	do_helical_symmetry_local_refinement = false;
    	helical_twist_initial = mymodel.helical_twist_min = mymodel.helical_twist_max = mymodel.helical_twist_inistep = 0.;
    	helical_rise_initial = mymodel.helical_rise_min = mymodel.helical_rise_max = mymodel.helical_rise_inistep = 0.;
    	helical_z_percentage = 0.;
    }
    mymodel.initialiseHelicalParametersLists(helical_twist_initial, helical_rise_initial);
    mymodel.is_helix = do_helical_refine;
    RFLOAT tmp_RFLOAT = 0.;
    if (mymodel.helical_rise_min > mymodel.helical_rise_max)
    	SWAP(mymodel.helical_rise_min, mymodel.helical_rise_max, tmp_RFLOAT);
    if (mymodel.helical_twist_min > mymodel.helical_twist_max)
    	SWAP(mymodel.helical_twist_min, mymodel.helical_twist_max, tmp_RFLOAT);

    // CTF, norm, scale, bfactor correction etc.
	int corrections_section = parser.addSection("Corrections");
	do_ctf_correction = parser.checkOption("--ctf", "Perform CTF correction?");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
	refs_are_ctf_corrected = parser.checkOption("--ctf_corrected_ref", "Have the input references been CTF-amplitude corrected?");
	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Have the data been CTF phase-flipped?");
	ctf_premultiplied = parser.checkOption("--ctf_multiplied", "Have the data been premultiplied with their CTF?");
	only_flip_phases = parser.checkOption("--only_flip_phases", "Only perform CTF phase-flipping? (default is full amplitude-correction)");
	do_norm_correction = parser.checkOption("--norm", "Perform normalisation-error correction?");
	do_scale_correction = parser.checkOption("--scale", "Perform intensity-scale corrections on image groups?");
	// Allow switching off norm and scale (which is on by default in the GUI)
	if (parser.checkOption("--no_norm", "Switch off normalisation-error correction?"))
		do_norm_correction = false;
	if (parser.checkOption("--no_scale", "Switch off intensity-scale corrections on image groups?"))
		do_scale_correction = false;

	// SGD stuff
	int sgd_section = parser.addSection("Stochastic Gradient Descent");
	do_sgd = parser.checkOption("--sgd", "Perform stochastic gradient descent instead of default expectation-maximization");
	mu = textToFloat(parser.getOption("--mu", "Momentum parameter for SGD updates", "0.9"));
	subset_size = textToInteger(parser.getOption("--subset_size", "Size of the subsets for SGD", "-1"));
	sgd_stepsize = textToFloat(parser.getOption("--sgd_stepsize", "Step size parameter for SGD updates", "0.5"));
	sgd_max_subsets = textToInteger(parser.getOption("--max_subsets", "Stop SGD after processing this many subsets (possibly more than 1 iteration)", "-1"));
	sgd_sigma2fudge_ini = textToFloat(parser.getOption("--sgd_sigma2fudge_initial", "Initial factor by which the noise variance will be multiplied for SGD (not used if halftime is negative)", "8"));
	sgd_sigma2fudge_halflife = textToInteger(parser.getOption("--sgd_sigma2fudge_halflife", "Initialise SGD with 8x higher noise-variance, and reduce with this half-life in # of particles (default is keep normal variance)", "-1"));
	write_every_subset = textToInteger(parser.getOption("--write_subsets", "Write out model every so many subsets (default is not writing any)", "-1"));
	strict_highres_sgd = textToFloat(parser.getOption("--strict_highres_sgd", "Resolution limit (in Angstrom) to restrict probability calculations in SGD", "20"));

	// Computation stuff
	// The number of threads is always read from the command line
	int computation_section = parser.addSection("Computation");
	x_pool = textToInteger(parser.getOption("--pool", "Number of images to pool for each thread task", "1"));
	nr_threads = textToInteger(parser.getOption("--j", "Number of threads to run in parallel (only useful on multi-core machines)", "1"));
	combine_weights_thru_disc = !parser.checkOption("--dont_combine_weights_via_disc", "Send the large arrays of summed weights through the MPI network, instead of writing large files to disc");
	do_shifts_onthefly = parser.checkOption("--onthefly_shifts", "Calculate shifted images on-the-fly, do not store precalculated ones in memory");
	do_parallel_disc_io = !parser.checkOption("--no_parallel_disc_io", "Do NOT let parallel (MPI) processes access the disc simultaneously (use this option with NFS)");
	do_preread_images  = parser.checkOption("--preread_images", "Use this to let the master process read all particles into memory. Be careful you have enough RAM for large data sets!");
	fn_scratch = parser.getOption("--scratch_dir", "If provided, particle stacks will be copied to this local scratch disk prior to refinement.", "");
	keep_free_scratch_Gb = textToInteger(parser.getOption("--keep_free_scratch", "Space available for copying particle stacks (in Gb)", "10"));
	do_reuse_scratch = parser.checkOption("--reuse_scratch", "Re-use data on scratchdir, instead of wiping it and re-copying all data.");

	do_gpu = parser.checkOption("--gpu", "Use available gpu resources for some calculations");
	gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread","default");
#ifndef CUDA
	if(do_gpu)
	{
		std::cerr << "+ WARNING : Relion was compiled without CUDA of at least version 7.0 - you do NOT have support for GPUs" << std::endl;
		do_gpu = false;
	}
#endif
	double temp_reqSize = textToDouble(parser.getOption("--free_gpu_memory", "GPU device memory (in Mb) to leave free after allocation.", "0"));
	temp_reqSize *= 1000*1000;
	if(temp_reqSize<0)
		REPORT_ERROR("Invalid free_gpu_memory value.");
	else
		requested_free_gpu_memory =  temp_reqSize;

	// Expert options
	int expert_section = parser.addSection("Expert options");
	mymodel.padding_factor = textToFloat(parser.getOption("--pad", "Oversampling factor for the Fourier transforms of the references", "2"));
	mymodel.interpolator = (parser.checkOption("--NN", "Perform nearest-neighbour instead of linear Fourier-space interpolation?")) ? NEAREST_NEIGHBOUR : TRILINEAR;
	mymodel.r_min_nn = textToInteger(parser.getOption("--r_min_nn", "Minimum number of Fourier shells to perform linear Fourier-space interpolation", "10"));
	verb = textToInteger(parser.getOption("--verb", "Verbosity (1=normal, 0=silent)", "1"));
	random_seed = textToInteger(parser.getOption("--random_seed", "Number for the random seed generator", "-1"));
	max_coarse_size = textToInteger(parser.getOption("--coarse_size", "Maximum image size for the first pass of the adaptive sampling approach", "-1"));
	adaptive_fraction = textToFloat(parser.getOption("--adaptive_fraction", "Fraction of the weights to be considered in the first pass of adaptive oversampling ", "0.999"));
	width_mask_edge = textToInteger(parser.getOption("--maskedge", "Width of the soft edge of the spherical mask (in pixels)", "5"));
	fix_sigma_noise = parser.checkOption("--fix_sigma_noise", "Fix the experimental noise spectra?");
	fix_sigma_offset = parser.checkOption("--fix_sigma_offset", "Fix the stddev in the origin offsets?");
	incr_size = textToInteger(parser.getOption("--incr_size", "Number of Fourier shells beyond the current resolution to be included in refinement", "10"));
	do_print_metadata_labels = parser.checkOption("--print_metadata_labels", "Print a table with definitions of all metadata labels, and exit");
	do_print_symmetry_ops = parser.checkOption("--print_symmetry_ops", "Print all symmetry transformation matrices, and exit");
	strict_highres_exp = textToFloat(parser.getOption("--strict_highres_exp", "Resolution limit (in Angstrom) to restrict probability calculations in the expectation step", "-1"));
	dont_raise_norm_error = parser.checkOption("--dont_check_norm", "Skip the check whether the images are normalised correctly");
	do_always_cc  = parser.checkOption("--always_cc", "Perform CC-calculation in all iterations (useful for faster denovo model generation?)");
	do_phase_random_fsc = parser.checkOption("--solvent_correct_fsc", "Correct FSC curve for the effects of the solvent mask?");
	do_skip_maximization = parser.checkOption("--skip_maximize", "Skip maximization step (only write out data.star file)?");
	failsafe_threshold = textToInteger(parser.getOption("--failsafe_threshold", "Maximum number of particles permitted to be handled by fail-safe mode, due to zero sum of weights, before exiting with an error (GPU only).", "40"));
	///////////////// Special stuff for first iteration (only accessible via CL, not through readSTAR ////////////////////

	// When reading from the CL: always start at iteration 1 and subset 1
	iter = 0;
	subset = 0;
	subset_start = 1;
    // When starting from CL: always calculate initial sigma_noise
    do_calculate_initial_sigma_noise = true;
    // Start average norm correction at 1!
    mymodel.avg_norm_correction = 1.;
    // Always initialise the PDF of the directions
    directions_have_changed = true;

    // Only reconstruct and join random halves are only available when continuing an old run
    do_join_random_halves = false;

    // For auto-sampling and convergence check
    nr_iter_wo_resol_gain = 0;
    nr_iter_wo_large_hidden_variable_changes = 0;
    current_changes_optimal_classes = 9999999;
    current_changes_optimal_offsets = 999.;
    current_changes_optimal_orientations = 999.;
    smallest_changes_optimal_classes = 9999999;
    smallest_changes_optimal_offsets = 999.;
    smallest_changes_optimal_orientations = 999.;
    acc_rot = acc_trans = 999.;

    best_resol_thus_far = 1./999.;
    has_converged = false;
    has_high_fsc_at_limit = false;
    has_large_incr_size_iter_ago = 0;

    // Never realign movies from the start
    do_realign_movies = false;

    // Debugging/analysis/hidden stuff
	do_map = !checkParameter(argc, argv, "--no_map");
	minres_map = textToInteger(getParameter(argc, argv, "--minres_map", "5"));
    do_bfactor = checkParameter(argc, argv, "--bfactor");
    gridding_nr_iter = textToInteger(getParameter(argc, argv, "--gridding_iter", "10"));
	debug1 = textToFloat(getParameter(argc, argv, "--debug1", "0"));
	debug2 = textToFloat(getParameter(argc, argv, "--debug2", "0"));
	debug3 = textToFloat(getParameter(argc, argv, "--debug3", "0"));
	// Read in initial sigmaNoise spectrum
	fn_sigma = getParameter(argc, argv, "--sigma","");
	do_calculate_initial_sigma_noise = (fn_sigma == "") ? true : false;
	sigma2_fudge = textToFloat(getParameter(argc, argv, "--sigma2_fudge", "1"));
	do_acc_currentsize_despite_highres_exp = checkParameter(argc, argv, "--accuracy_current_size");
	do_sequential_halves_recons  = checkParameter(argc, argv, "--sequential_halves_recons");
	do_always_join_random_halves = checkParameter(argc, argv, "--always_join_random_halves");
	do_use_all_data = checkParameter(argc, argv, "--use_all_data");
	do_only_sample_tilt  = checkParameter(argc, argv, "--only_sample_tilt");
	minimum_angular_sampling = textToFloat(getParameter(argc, argv, "--minimum_angular_sampling", "0"));
	asymmetric_padding = parser.checkOption("--asymmetric_padding", "", "false", true);
	maximum_significants = textToInteger(parser.getOption("--maxsig", "", "0", true));
	skip_gridding = parser.checkOption("--skip_gridding", "", "false", true);

#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::parseInitial Done"<<std::endl;
#endif

}


void MlOptimiser::read(FileName fn_in, int rank)
{
//#define DEBUG_READ
#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar entering ..."<<std::endl;
#endif

    // Open input file
    std::ifstream in(fn_in.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "MlOptimiser::readStar: File " + fn_in + " cannot be read." );

    MetaDataTable MD;

    // Read general stuff
    FileName fn_model, fn_model2, fn_sampling;
    MD.readStar(in, "optimiser_general");
    in.close();

    if (!MD.getValue(EMDL_OPTIMISER_OUTPUT_ROOTNAME, fn_out) ||
        !MD.getValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model) ||
		!MD.getValue(EMDL_OPTIMISER_DATA_STARFILE, fn_data) ||
		!MD.getValue(EMDL_OPTIMISER_SAMPLING_STARFILE, fn_sampling) ||
        !MD.getValue(EMDL_OPTIMISER_ITERATION_NO, iter) ||
        !MD.getValue(EMDL_OPTIMISER_NR_ITERATIONS, nr_iter) ||
        !MD.getValue(EMDL_OPTIMISER_DO_SPLIT_RANDOM_HALVES, do_split_random_halves) ||
        !MD.getValue(EMDL_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, low_resol_join_halves) ||
        !MD.getValue(EMDL_OPTIMISER_ADAPTIVE_OVERSAMPLING, adaptive_oversampling) ||
		!MD.getValue(EMDL_OPTIMISER_ADAPTIVE_FRACTION, adaptive_fraction) ||
		!MD.getValue(EMDL_OPTIMISER_RANDOM_SEED, random_seed) ||
		!MD.getValue(EMDL_OPTIMISER_PARTICLE_DIAMETER, particle_diameter) ||
		!MD.getValue(EMDL_OPTIMISER_WIDTH_MASK_EDGE, width_mask_edge) ||
		!MD.getValue(EMDL_OPTIMISER_DO_ZERO_MASK, do_zero_mask) ||
		!MD.getValue(EMDL_OPTIMISER_DO_SOLVENT_FLATTEN, do_solvent) ||
		!MD.getValue(EMDL_OPTIMISER_SOLVENT_MASK_NAME, fn_mask) ||
		!MD.getValue(EMDL_OPTIMISER_SOLVENT_MASK2_NAME, fn_mask2) ||
		!MD.getValue(EMDL_OPTIMISER_TAU_SPECTRUM_NAME, fn_tau) ||
		!MD.getValue(EMDL_OPTIMISER_COARSE_SIZE, coarse_size) ||
		!MD.getValue(EMDL_OPTIMISER_MAX_COARSE_SIZE, max_coarse_size) ||
		!MD.getValue(EMDL_OPTIMISER_HIGHRES_LIMIT_EXP, strict_highres_exp) ||
		!MD.getValue(EMDL_OPTIMISER_INCR_SIZE, incr_size) ||
		!MD.getValue(EMDL_OPTIMISER_DO_MAP, do_map) ||
		!MD.getValue(EMDL_OPTIMISER_DO_AUTO_REFINE, do_auto_refine) ||
		!MD.getValue(EMDL_OPTIMISER_AUTO_LOCAL_HP_ORDER, autosampling_hporder_local_searches) ||
	    !MD.getValue(EMDL_OPTIMISER_NR_ITER_WO_RESOL_GAIN, nr_iter_wo_resol_gain) ||
	    !MD.getValue(EMDL_OPTIMISER_BEST_RESOL_THUS_FAR, best_resol_thus_far) ||
	    !MD.getValue(EMDL_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, nr_iter_wo_large_hidden_variable_changes) ||
		!MD.getValue(EMDL_OPTIMISER_DO_SKIP_ALIGN, do_skip_align) ||
		//!MD.getValue(EMDL_OPTIMISER_DO_SKIP_ROTATE, do_skip_rotate) ||
	    !MD.getValue(EMDL_OPTIMISER_ACCURACY_ROT, acc_rot) ||
	    !MD.getValue(EMDL_OPTIMISER_ACCURACY_TRANS, acc_trans) ||
	    !MD.getValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, current_changes_optimal_orientations) ||
	    !MD.getValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, current_changes_optimal_offsets) ||
	    !MD.getValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_CLASSES, current_changes_optimal_classes) ||
	    !MD.getValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, smallest_changes_optimal_orientations) ||
	    !MD.getValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, smallest_changes_optimal_offsets) ||
	    !MD.getValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, smallest_changes_optimal_classes) ||
	    !MD.getValue(EMDL_OPTIMISER_HAS_CONVERGED, has_converged) ||
	    !MD.getValue(EMDL_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, has_high_fsc_at_limit) ||
	    !MD.getValue(EMDL_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, has_large_incr_size_iter_ago) ||
		!MD.getValue(EMDL_OPTIMISER_DO_CORRECT_NORM, do_norm_correction) ||
		!MD.getValue(EMDL_OPTIMISER_DO_CORRECT_SCALE, do_scale_correction) ||
		!MD.getValue(EMDL_OPTIMISER_DO_CORRECT_CTF, do_ctf_correction) ||
		!MD.getValue(EMDL_OPTIMISER_DO_REALIGN_MOVIES, do_realign_movies) ||
		!MD.getValue(EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, intact_ctf_first_peak) ||
		!MD.getValue(EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, ctf_phase_flipped) ||
		!MD.getValue(EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, only_flip_phases) ||
		!MD.getValue(EMDL_OPTIMISER_REFS_ARE_CTF_CORRECTED, refs_are_ctf_corrected) ||
		!MD.getValue(EMDL_OPTIMISER_FIX_SIGMA_NOISE, fix_sigma_noise) ||
		!MD.getValue(EMDL_OPTIMISER_FIX_SIGMA_OFFSET, fix_sigma_offset) ||
		!MD.getValue(EMDL_OPTIMISER_MAX_NR_POOL, nr_pool)  )
    	REPORT_ERROR("MlOptimiser::readStar: incorrect optimiser_general table");

    // Backward compatibility with RELION-1.4
    if (!MD.getValue(EMDL_OPTIMISER_LOCAL_SYMMETRY_FILENAME, fn_local_symmetry))
    	fn_local_symmetry = "None";
    if (!MD.getValue(EMDL_OPTIMISER_DO_HELICAL_REFINE, do_helical_refine))
    	do_helical_refine = false;
    if (!MD.getValue(EMDL_OPTIMISER_IGNORE_HELICAL_SYMMETRY, ignore_helical_symmetry))
    	ignore_helical_symmetry = false;
    if (!MD.getValue(EMDL_OPTIMISER_HELICAL_TWIST_INITIAL, helical_twist_initial))
    	helical_twist_initial = 0.;
    if (!MD.getValue(EMDL_OPTIMISER_HELICAL_RISE_INITIAL, helical_rise_initial))
    	helical_rise_initial = 0.;
	if (!MD.getValue(EMDL_OPTIMISER_HELICAL_Z_PERCENTAGE, helical_z_percentage))
		helical_z_percentage = 0.3;
	if (!MD.getValue(EMDL_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER, helical_tube_inner_diameter))
		helical_tube_inner_diameter = -1.;
	if (!MD.getValue(EMDL_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER, helical_tube_outer_diameter))
		helical_tube_outer_diameter = -1.;
    if (!MD.getValue(EMDL_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT, do_helical_symmetry_local_refinement))
    	do_helical_symmetry_local_refinement = false;
    if (!MD.getValue(EMDL_OPTIMISER_HELICAL_SIGMA_DISTANCE, helical_sigma_distance))
    	helical_sigma_distance = -1.;
    if (!MD.getValue(EMDL_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED, helical_keep_tilt_prior_fixed))
    	helical_keep_tilt_prior_fixed = false;
	if (!MD.getValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, ctf_premultiplied))
		ctf_premultiplied = false;
	if (!MD.getValue(EMDL_OPTIMISER_HIGHRES_LIMIT_SGD, strict_highres_sgd))
		strict_highres_sgd = -1.;
	if (!MD.getValue(EMDL_OPTIMISER_DO_SGD, do_sgd))
		do_sgd = false;
	if (!MD.getValue(EMDL_OPTIMISER_SGD_MU, mu))
		mu = 0.9;
	if (!MD.getValue(EMDL_OPTIMISER_SGD_SIGMA2FUDGE_INI, sgd_sigma2fudge_ini))
		sgd_sigma2fudge_ini = 8.;
	if (!MD.getValue(EMDL_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE, sgd_sigma2fudge_halflife))
		sgd_sigma2fudge_halflife = -1;
	if (!MD.getValue(EMDL_OPTIMISER_SGD_SUBSET_SIZE, subset_size))
		subset_size = -1;
	if (!MD.getValue(EMDL_OPTIMISER_SGD_SUBSET_START, subset_start))
		subset_start = 1;
	if (!MD.getValue(EMDL_OPTIMISER_SGD_STEPSIZE, sgd_stepsize))
		sgd_stepsize = 0.5;
	// The following line is only to avoid strange filenames when writing optimiser.star upon restarting
	subset = subset_start - 1;
	if (!MD.getValue(EMDL_OPTIMISER_SGD_WRITE_EVERY_SUBSET, write_every_subset))
		write_every_subset = -1;
	if (!MD.getValue(EMDL_OPTIMISER_SGD_MAX_SUBSETS, sgd_max_subsets))
		sgd_max_subsets = -1;

    if (do_split_random_halves &&
    		!MD.getValue(EMDL_OPTIMISER_MODEL_STARFILE2, fn_model2))
    	REPORT_ERROR("MlOptimiser::readStar: splitting data into two random halves, but rlnModelStarFile2 not found in optimiser_general table");

    // Initialise some stuff for first-iteration only (not relevant here...)
    do_calculate_initial_sigma_noise = false;
    do_average_unaligned = false;
    do_generate_seeds = false;
    do_firstiter_cc = false;
    ini_high = 0;

    // Initialise some of the other, hidden or debugging stuff
    minres_map = 5;
    do_bfactor = false;
    gridding_nr_iter = 10;
    debug1 = debug2 = debug3 = 0.;

    // Then read in sampling, mydata and mymodel stuff
    // If do_preread_images: when not do_parallel_disc_io: only the master reads all images into RAM; otherwise: everyone reads in images into RAM
#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar before data."<<std::endl;
#endif
    bool do_preread = (do_preread_images) ? (do_parallel_disc_io || rank == 0) : false;
    bool is_helical_segment = (do_helical_refine) || ((mymodel.ref_dim == 2) && (helical_tube_outer_diameter > 0.));
    mydata.read(fn_data, false, false, do_preread, is_helical_segment);

#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar before model."<<std::endl;
#endif
    if (do_split_random_halves)
	{
		if (rank % 2 == 1)
			mymodel.read(fn_model);
		else
			mymodel.read(fn_model2);
	}
	else
	{
		mymodel.read(fn_model);
	}
#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar before sampling."<<std::endl;
#endif
	sampling.read(fn_sampling);

#ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar done."<<std::endl;
#endif

}


void MlOptimiser::write(bool do_write_sampling, bool do_write_data, bool do_write_optimiser, bool do_write_model, int random_subset)
{

	FileName fn_root, fn_tmp, fn_model, fn_model2, fn_data, fn_sampling, fn_root2;
	std::ofstream  fh;
	if (iter > -1)
		fn_root.compose(fn_out+"_it", iter, "", 3);
	else
		fn_root = fn_out;
	// fn_root2 is used to write out the model and optimiser, and adds a subset number in SGD
	fn_root2 = fn_root;
	bool do_write_bild = !(do_skip_align || do_skip_rotate);
	int writeout_subset_start = 1;

	if (iter > 0 && nr_subsets > 1 && subset != nr_subsets)
	{
		if ((subset % write_every_subset) == 0)
		{
			fn_root2.compose(fn_root+"_sub", subset, "", 4);
			fn_root += "_sub";
			do_write_bild = false;
			writeout_subset_start = subset + 1;
		}
		else
		{
			do_write_sampling = false;
			do_write_data = false;
			do_write_optimiser = false;
			do_write_model = false;
		}
	}

	// First write "main" STAR file with all information from this run
	// Do this for random_subset==0 and random_subset==1
	if (do_write_optimiser && random_subset < 2)
	{
		fn_tmp = fn_root2+"_optimiser.star";
		fh.open((fn_tmp).c_str(), std::ios::out);
		if (!fh)
			REPORT_ERROR( (std::string)"MlOptimiser::write: Cannot write file: " + fn_tmp);

		// Write the command line as a comment in the header
		fh << "# RELION optimiser"<<std::endl;
		fh << "# ";
		parser.writeCommandLine(fh);

		if (do_split_random_halves && !do_join_random_halves)
		{
			fn_model  = fn_root2 + "_half1_model.star";
			fn_model2 = fn_root2 + "_half2_model.star";
		}
		else
		{
			fn_model = fn_root2 + "_model.star";
		}
		fn_data = fn_root + "_data.star";
		fn_sampling = fn_root + "_sampling.star";

		MetaDataTable MD;
		MD.setIsList(true);
		MD.setName("optimiser_general");
		MD.addObject();
		MD.setValue(EMDL_OPTIMISER_OUTPUT_ROOTNAME, fn_out);
		if (do_split_random_halves)
		{
			MD.setValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);
			MD.setValue(EMDL_OPTIMISER_MODEL_STARFILE2, fn_model2);
		}
		else
		{
			MD.setValue(EMDL_OPTIMISER_MODEL_STARFILE, fn_model);
		}
		MD.setValue(EMDL_OPTIMISER_DATA_STARFILE, fn_data);
		MD.setValue(EMDL_OPTIMISER_SAMPLING_STARFILE, fn_sampling);
		MD.setValue(EMDL_OPTIMISER_ITERATION_NO, iter);
		MD.setValue(EMDL_OPTIMISER_NR_ITERATIONS, nr_iter);
		MD.setValue(EMDL_OPTIMISER_DO_SPLIT_RANDOM_HALVES, do_split_random_halves);
		MD.setValue(EMDL_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, low_resol_join_halves);
		MD.setValue(EMDL_OPTIMISER_ADAPTIVE_OVERSAMPLING, adaptive_oversampling);
		MD.setValue(EMDL_OPTIMISER_ADAPTIVE_FRACTION, adaptive_fraction);
		MD.setValue(EMDL_OPTIMISER_RANDOM_SEED, random_seed);
		MD.setValue(EMDL_OPTIMISER_PARTICLE_DIAMETER, particle_diameter);
		MD.setValue(EMDL_OPTIMISER_WIDTH_MASK_EDGE, width_mask_edge);
		MD.setValue(EMDL_OPTIMISER_DO_ZERO_MASK, do_zero_mask);
		MD.setValue(EMDL_OPTIMISER_DO_SOLVENT_FLATTEN, do_solvent);
		MD.setValue(EMDL_OPTIMISER_SOLVENT_MASK_NAME, fn_mask);
		MD.setValue(EMDL_OPTIMISER_SOLVENT_MASK2_NAME, fn_mask2);
		MD.setValue(EMDL_OPTIMISER_TAU_SPECTRUM_NAME, fn_tau);
		MD.setValue(EMDL_OPTIMISER_COARSE_SIZE, coarse_size);
		MD.setValue(EMDL_OPTIMISER_MAX_COARSE_SIZE, max_coarse_size);
		MD.setValue(EMDL_OPTIMISER_HIGHRES_LIMIT_EXP, strict_highres_exp);
		MD.setValue(EMDL_OPTIMISER_INCR_SIZE, incr_size);
		MD.setValue(EMDL_OPTIMISER_DO_MAP, do_map);
		MD.setValue(EMDL_OPTIMISER_DO_SGD, do_sgd);
		MD.setValue(EMDL_OPTIMISER_SGD_MU, mu);
		MD.setValue(EMDL_OPTIMISER_SGD_SIGMA2FUDGE_INI, sgd_sigma2fudge_ini);
		MD.setValue(EMDL_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE, sgd_sigma2fudge_halflife);
		MD.setValue(EMDL_OPTIMISER_SGD_SUBSET_START, writeout_subset_start);
		MD.setValue(EMDL_OPTIMISER_SGD_SUBSET_SIZE, subset_size);
		MD.setValue(EMDL_OPTIMISER_SGD_WRITE_EVERY_SUBSET, write_every_subset);
		MD.setValue(EMDL_OPTIMISER_SGD_MAX_SUBSETS, sgd_max_subsets);
		MD.setValue(EMDL_OPTIMISER_SGD_STEPSIZE, sgd_stepsize);
		MD.setValue(EMDL_OPTIMISER_HIGHRES_LIMIT_SGD, strict_highres_sgd);
		MD.setValue(EMDL_OPTIMISER_DO_AUTO_REFINE, do_auto_refine);
		MD.setValue(EMDL_OPTIMISER_AUTO_LOCAL_HP_ORDER, autosampling_hporder_local_searches);
	    MD.setValue(EMDL_OPTIMISER_NR_ITER_WO_RESOL_GAIN, nr_iter_wo_resol_gain);
	    MD.setValue(EMDL_OPTIMISER_BEST_RESOL_THUS_FAR,best_resol_thus_far);
	    MD.setValue(EMDL_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, nr_iter_wo_large_hidden_variable_changes);
		MD.setValue(EMDL_OPTIMISER_DO_SKIP_ALIGN, do_skip_align);
		MD.setValue(EMDL_OPTIMISER_DO_SKIP_ROTATE, do_skip_rotate);
	    MD.setValue(EMDL_OPTIMISER_ACCURACY_ROT, acc_rot);
	    MD.setValue(EMDL_OPTIMISER_ACCURACY_TRANS, acc_trans);
	    MD.setValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, current_changes_optimal_orientations);
	    MD.setValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, current_changes_optimal_offsets);
	    MD.setValue(EMDL_OPTIMISER_CHANGES_OPTIMAL_CLASSES, current_changes_optimal_classes);
	    MD.setValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, smallest_changes_optimal_orientations);
	    MD.setValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, smallest_changes_optimal_offsets);
	    MD.setValue(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, smallest_changes_optimal_classes);
	    MD.setValue(EMDL_OPTIMISER_LOCAL_SYMMETRY_FILENAME, fn_local_symmetry);
	    MD.setValue(EMDL_OPTIMISER_DO_HELICAL_REFINE, do_helical_refine);
	    MD.setValue(EMDL_OPTIMISER_IGNORE_HELICAL_SYMMETRY, ignore_helical_symmetry);
	    MD.setValue(EMDL_OPTIMISER_HELICAL_TWIST_INITIAL, helical_twist_initial);
	    MD.setValue(EMDL_OPTIMISER_HELICAL_RISE_INITIAL, helical_rise_initial);
	    MD.setValue(EMDL_OPTIMISER_HELICAL_Z_PERCENTAGE, helical_z_percentage);
	    MD.setValue(EMDL_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER, helical_tube_inner_diameter);
	    MD.setValue(EMDL_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER, helical_tube_outer_diameter);
	    MD.setValue(EMDL_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT, do_helical_symmetry_local_refinement);
	    MD.setValue(EMDL_OPTIMISER_HELICAL_SIGMA_DISTANCE, helical_sigma_distance);
	    MD.setValue(EMDL_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED, helical_keep_tilt_prior_fixed);
	    MD.setValue(EMDL_OPTIMISER_HAS_CONVERGED, has_converged);
	    MD.setValue(EMDL_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, has_high_fsc_at_limit);
	    MD.setValue(EMDL_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, has_large_incr_size_iter_ago);
		MD.setValue(EMDL_OPTIMISER_DO_CORRECT_NORM, do_norm_correction);
		MD.setValue(EMDL_OPTIMISER_DO_CORRECT_SCALE, do_scale_correction);
		MD.setValue(EMDL_OPTIMISER_DO_CORRECT_CTF, do_ctf_correction);
		MD.setValue(EMDL_OPTIMISER_DO_REALIGN_MOVIES, do_realign_movies);
		MD.setValue(EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, intact_ctf_first_peak);
		MD.setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, ctf_phase_flipped);
		MD.setValue(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, ctf_premultiplied);
		MD.setValue(EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, only_flip_phases);
		MD.setValue(EMDL_OPTIMISER_REFS_ARE_CTF_CORRECTED, refs_are_ctf_corrected);
		MD.setValue(EMDL_OPTIMISER_FIX_SIGMA_NOISE, fix_sigma_noise);
		MD.setValue(EMDL_OPTIMISER_FIX_SIGMA_OFFSET, fix_sigma_offset);
		MD.setValue(EMDL_OPTIMISER_MAX_NR_POOL, nr_pool);

		MD.write(fh);
		fh.close();
	}

	// Then write the mymodel to file
	if (do_write_model)
	{
		if (do_split_random_halves && !do_join_random_halves)
			mymodel.write(fn_root2 + "_half" + integerToString(random_subset), sampling, do_write_bild);
		else
			mymodel.write(fn_root2, sampling, do_write_bild);
	}

	// And write the mydata to file
	if (do_write_data)
		mydata.write(fn_root);

	// And write the sampling object
	if (do_write_sampling)
		sampling.write(fn_root);

}

/** ========================== Initialisation  =========================== */

void MlOptimiser::initialise()
{
#ifdef DEBUG
    std::cerr<<"MlOptimiser::initialise Entering"<<std::endl;
#endif

	if (do_gpu)
	{
#ifdef CUDA
		int devCount;
		HANDLE_ERROR(cudaGetDeviceCount(&devCount));

		cudaDeviceProp deviceProp;
		int compatibleDevices(0);
		// Send device count seen by this slave
		HANDLE_ERROR(cudaGetDeviceCount(&devCount));
		for(int i=0; i<devCount; i++ )
                {
                	HANDLE_ERROR(cudaGetDeviceProperties(&deviceProp, i));
                        if(deviceProp.major>CUDA_CC_MAJOR)
                        	compatibleDevices+=1;
                       	else if(deviceProp.major==CUDA_CC_MAJOR && deviceProp.minor>=CUDA_CC_MINOR)
                        	compatibleDevices+=1;
                        //else
                        //std::cout << "Found a " << deviceProp.name << " GPU with compute-capability " << deviceProp.major << "." << deviceProp.minor << std::endl;
                }
		if(compatibleDevices==0)
			REPORT_ERROR("You have no GPUs compatible with RELION (CUDA-capable and compute-capability >= 3.5");
		else if(compatibleDevices!=devCount)
			std::cerr << "WARNING : at least one of your GPUs is not compatible with RELION (CUDA-capable and compute-capability >= 3.5)" << std::endl;

		std::vector < std::vector < std::string > > allThreadIDs;
		untangleDeviceIDs(gpu_ids, allThreadIDs);

		// Sequential initialisation of GPUs on all ranks
		bool fullAutomaticMapping(true);
		bool semiAutomaticMapping(true);
		if (allThreadIDs[0].size()==0 || (!std::isdigit(*gpu_ids.begin())) )
			std::cout << "gpu-ids not specified, threads will automatically be mapped to devices (incrementally)."<< std::endl;
		else
		{
			fullAutomaticMapping=false;
			if(allThreadIDs[0].size()!=nr_threads)
			{
				std::cout << " Will distribute threads over devices ";
				for (int j = 0; j < allThreadIDs[0].size(); j++)
					std::cout << " "  << allThreadIDs[0][j];
				std::cout  << std::endl;
			}
			else
			{
				semiAutomaticMapping = false;
				std::cout << " Using explicit indexing to assign devices ";
				for (int j = 0; j < allThreadIDs[0].size(); j++)
					std::cout << " "  << allThreadIDs[0][j];
				std::cout  << std::endl;
			}
		}

		for (int i = 0; i < nr_threads; i ++)
		{
			int dev_id;
			if (semiAutomaticMapping)
			{
				// Sjors: hack to make use of several cards; will only work if all MPI slaves are on the same node!
				// Bjorn: Better hack
				if (fullAutomaticMapping)
					dev_id = devCount*i / nr_threads;
				else
					dev_id = textToInteger(allThreadIDs[0][i % (allThreadIDs[0]).size()].c_str());
					//dev_id = textToInteger(allThreadIDs[0][(allThreadIDs[0].size()*i)/nr_threads].c_str());
			}
			else // not semiAutomatic => explicit
			{
					dev_id = textToInteger(allThreadIDs[0][i].c_str());
			}


			std::cout << " Thread " << i << " mapped to device " << dev_id << std::endl;

			//Only make a new bundle of not existing on device
			int bundleId(-1);

			for (int j = 0; j < cudaDevices.size(); j++)
				if (cudaDevices[j] == dev_id)
					bundleId = j;

			if (bundleId == -1)
			{
				bundleId = cudaDevices.size();
				cudaDevices.push_back(dev_id);
			}

			cudaOptimiserDeviceMap.push_back(bundleId);
		}
#else
        REPORT_ERROR("GPU usage requested, but RELION was compiled without CUDA support");
#endif
	}

    initialiseGeneral();

    initialiseWorkLoad();

	if (fn_sigma != "")
	{
		// Read in sigma_noise spetrum from file DEVELOPMENTAL!!! FOR DEBUGGING ONLY....
		MetaDataTable MDsigma;
		RFLOAT val;
		int idx;
		MDsigma.read(fn_sigma);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDsigma)
		{
			MDsigma.getValue(EMDL_SPECTRAL_IDX, idx);
			MDsigma.getValue(EMDL_MLMODEL_SIGMA2_NOISE, val);
			if (idx < XSIZE(mymodel.sigma2_noise[0]))
				mymodel.sigma2_noise[0](idx) = val;
		}
		if (idx < XSIZE(mymodel.sigma2_noise[0]) - 1)
		{
			if (verb > 0) std::cout<< " WARNING: provided sigma2_noise-spectrum has fewer entries ("<<idx+1<<") than needed ("<<XSIZE(mymodel.sigma2_noise[0])<<"). Set rest to zero..."<<std::endl;
		}
		// Use the same spectrum for all classes
		for (int igroup = 0; igroup< mymodel.nr_groups; igroup++)
			mymodel.sigma2_noise[igroup] =  mymodel.sigma2_noise[0];

	}
	else if (do_calculate_initial_sigma_noise || do_average_unaligned)
	{
		MultidimArray<RFLOAT> Mavg;

		// Calculate initial sigma noise model from power_class spectra of the individual images
		calculateSumOfPowerSpectraAndAverageImage(Mavg);

		// Set sigma2_noise and Iref from averaged poser spectra and Mavg
		setSigmaNoiseEstimatesAndSetAverageImage(Mavg);
	}

	// First low-pass filter the initial references
	if (iter == 0)
		initialLowPassFilterReferences();

	// Initialise the data_versus_prior ratio to get the initial current_size right
	if (iter == 0)
		mymodel.initialiseDataVersusPrior(fix_tau); // fix_tau was set in initialiseGeneral

	// Check minimum group size of 10 particles
	if (verb > 0)
	{
		bool do_warn = false;
		for (int igroup = 0; igroup< mymodel.nr_groups; igroup++)
		{
			if (mymodel.nr_particles_group[igroup] < 10)
			{
				std:: cout << "WARNING: There are only " << mymodel.nr_particles_group[igroup] << " particles in group " << igroup + 1 << std::endl;
				do_warn = true;
			}
		}
		if (do_warn)
		{
			std:: cout << "WARNING: You may want to consider joining some micrographs into larger groups to obtain more robust noise estimates. " << std::endl;
			std:: cout << "         You can do so by using the same rlnMicrographName label for particles from multiple different micrographs in the input STAR file. " << std::endl;
		}
	}

	// Write out initial mymodel
	if (!do_movies_in_batches)
		write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, 0);

	// Do this after writing out the model, so that still the random halves are written in separate files.
	if (do_realign_movies)
	{
		// Resolution seems to decrease again after 1 iteration. Therefore, just perform a single iteration until we figure out what exactly happens here...
		has_converged = true;
		// Then use join random halves
		do_join_random_halves = true;

		// If we skip the maximization step, then there is no use in using all data
		if (!do_skip_maximization)
		{
			// Use all data out to Nyquist because resolution gains may be substantial
			do_use_all_data = true;
		}
	}


#ifdef DEBUG
    std::cerr<<"MlOptimiser::initialise Done"<<std::endl;
#endif
}

void MlOptimiser::initialiseGeneral(int rank)
{

#ifdef DEBUG
	std::cerr << "Entering initialiseGeneral" << std::endl;
#endif

#ifdef TIMING
	//DIFFF = timer.setNew("difff");
	TIMING_EXP =           timer.setNew("expectation");
	TIMING_EXP_1 =           timer.setNew("expectation_1");
	TIMING_EXP_1a =           timer.setNew("expectation_1a");
	TIMING_EXP_2 =           timer.setNew("expectation_2");
	TIMING_EXP_3 =           timer.setNew("expectation_3");
	TIMING_EXP_4 =           timer.setNew("expectation_4");
	TIMING_EXP_4a =           timer.setNew("expectation_4a");
	TIMING_EXP_4b =           timer.setNew("expectation_4b");
	TIMING_EXP_4c =           timer.setNew("expectation_4c");
	TIMING_EXP_4d =           timer.setNew("expectation_4d");
	TIMING_EXP_5 =           timer.setNew("expectation_5");
	TIMING_EXP_6 =           timer.setNew("expectation_6");
	TIMING_EXP_7 =           timer.setNew("expectation_7");
	TIMING_EXP_8 =           timer.setNew("expectation_8");
	TIMING_EXP_9 =           timer.setNew("expectation_9");

	TIMING_EXP_METADATA =  timer.setNew(" - EXP: metadata shuffling");
	TIMING_EXP_CHANGES =   timer.setNew(" - EXP: monitor changes hidden variables");
	TIMING_MAX =           timer.setNew("maximization");
	TIMING_SOLVFLAT =      timer.setNew("flatten solvent");
	TIMING_UPDATERES =     timer.setNew("update resolution");
	TIMING_RECONS =        timer.setNew("reconstruction");
	TIMING_ESP =           timer.setNew("expectationSomeParticles");
	TIMING_ESP_THR =       timer.setNew("doThreadExpectationSomeParticles");
	TIMING_ESP_ONEPART =   timer.setNew("expectationOneParticle (thr0)");
	TIMING_ESP_ONEPARTN =  timer.setNew("expectationOneParticle (thrN)");
	TIMING_ESP_INI=        timer.setNew(" - EOP: initialise memory");
	TIMING_ESP_FT =        timer.setNew(" - EOP: getFourierTransforms");
	TIMING_ESP_PREC1 =     timer.setNew(" - EOP: precalcShifts1");
	TIMING_ESP_PREC2 =     timer.setNew(" - EOP: precalcShifts2");
	TIMING_ESP_DIFF1 =     timer.setNew(" - EOP: getAllSquaredDifferences1");
	TIMING_ESP_DIFF2 =     timer.setNew(" - EOP: getAllSquaredDifferences2");
	TIMING_ESP_DIFF2_A =   timer.setNew(" - EOP: getD2: A");
	TIMING_ESP_DIFF2_B =   timer.setNew(" - EOP: getD2: B");
	TIMING_ESP_DIFF2_C =   timer.setNew(" - EOP: getD2: C");
	TIMING_ESP_DIFF2_D =   timer.setNew(" - EOP: getD2: D");
	TIMING_ESP_DIFF2_E =   timer.setNew(" - EOP: getD2: E");
	TIMING_DIFF_PROJ =     timer.setNew(" -  - EOPdiff2: project");
	TIMING_DIFF_SHIFT =    timer.setNew(" -  - EOPdiff2: shift");
	TIMING_DIFF2_GETSHIFT =timer.setNew(" -  - EOPdiff2: get shifted img");
	TIMING_DIFF_DIFF2 =    timer.setNew(" -  - EOPdiff2: diff2");
	TIMING_ESP_WEIGHT1 =   timer.setNew(" - EOP: convertDiff2ToWeights1");
	TIMING_ESP_WEIGHT2 =   timer.setNew(" - EOP: convertDiff2ToWeights2");
	TIMING_WEIGHT_EXP =    timer.setNew(" -  - EOPweight: exp");
	TIMING_WEIGHT_SORT =   timer.setNew(" -  - EOPweight: sort");
	TIMING_ESP_WSUM =      timer.setNew(" - EOP: storeWeightedSums");
	TIMING_ESP_PRECW =     timer.setNew(" -  - EOPwsum: precalcShiftsW");
	TIMING_WSUM_PROJ =     timer.setNew(" -  - EOPwsum: project");
	TIMING_WSUM_GETSHIFT = timer.setNew(" -  - EOPwsum: get shifted img");
	TIMING_WSUM_DIFF2 =    timer.setNew(" -  - EOPwsum: diff2");
	TIMING_WSUM_LOCALSUMS =timer.setNew(" -  - EOPwsum: localsums");
	TIMING_WSUM_SUMSHIFT = timer.setNew(" -  - EOPwsum: shiftimg");
	TIMING_WSUM_BACKPROJ = timer.setNew(" -  - EOPwsum: backproject");

	TIMING_EXTRA1= timer.setNew(" -extra1");
	TIMING_EXTRA2= timer.setNew(" -extra2");
	TIMING_EXTRA3= timer.setNew(" -extra3");

	RCT_1 = timer.setNew(" RcT1_BPrefRecon ");
	RCT_2 = timer.setNew(" RcT2_BroadCast ");
	RCT_3 = timer.setNew(" RcT3_maxiOthers ");
	RCT_4 = timer.setNew(" RcT4_monitorHidden ");
	RCT_5 = timer.setNew(" RcT5_sums ");
	RCT_6 = timer.setNew(" RcT6_updatePdf ");
	RCT_7 = timer.setNew(" RcT7_updateNoise ");
	RCT_8 = timer.setNew(" RcT8_initials ");

#endif

	// Check for errors in the command-line option
	if (parser.checkForErrors(verb))
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

#ifdef RELION_SINGLE_PRECISION
        if (verb > 0)
            std::cout << " Running CPU instructions in single precision. Runs might not be exactly reproducible." << std::endl;
#else
        if (verb > 0)
            std::cout << " Running CPU instructions in double precision. " << std::endl;
#endif

    // Check if output directory exists
    FileName fn_dir = fn_out.beforeLastOf("/");
    if (!exists(fn_dir))
    	REPORT_ERROR("ERROR: output directory does not exist!");

    // Just die if trying to use GPUs and skipping alignments
    if (do_skip_align && do_gpu)
    	REPORT_ERROR("ERROR: you cannot use GPUs when skipping alignments");

	if (do_always_cc)
		do_calculate_initial_sigma_noise = false;

    if (do_print_metadata_labels)
	{
		if (verb > 0)
			EMDL::printDefinitions(std::cout);
		exit(0);
	}

	// Print symmetry operators to cout
	if (do_print_symmetry_ops)
	{
		if (verb > 0)
		{
			SymList SL;
			SL.writeDefinition(std::cout, sampling.symmetryGroup());
		}
		exit(0);
	}

	// If we are not continuing an old run, now read in the data and the reference images
	if (iter == 0)
	{

		// Read in the experimental image metadata
		// If do_preread_images: only the master reads all images into RAM
		bool do_preread = (do_preread_images) ? (do_parallel_disc_io || rank == 0) : false;
		if (do_realign_movies)
			do_preread = false; // as we will overwrite mydata.read with the movies anyway....
		bool is_helical_segment = (do_helical_refine) || ((mymodel.ref_dim == 2) && (helical_tube_outer_diameter > 0.));
		mydata.read(fn_data, true, false, do_preread, is_helical_segment); // true means ignore original particle name

		if (fn_body_masks != "")
		{
			// TODO: fix this!
			do_norm_correction = false;

			// Read in the masks in the input STAR file, add a soft-edge to them, and write them to disc with the standard name
			mymodel.initialiseBodyMasks(fn_body_masks, fn_out);

			// For multi-body refinement: expand the MetaDataTables with orientations for all bodies
			mydata.initialiseBodies(mymodel.nr_bodies);
		}
		else
		{
			mymodel.nr_bodies = 1;
			mydata.initialiseBodies(mymodel.nr_bodies);
		}

		// Also get original size of the images to pass to mymodel.read()
		int ori_size = -1;
		mydata.MDexp.getValue(EMDL_IMAGE_SIZE, ori_size);
		if (ori_size%2 != 0)
			REPORT_ERROR("This program only works with even values for the image dimensions!");

		// Read in the reference(s) and initialise mymodel
		int refdim = (fn_ref == "denovo") ? 3 : 2;
		mymodel.readImages(fn_ref, is_3d_model, ori_size, mydata,
				do_average_unaligned, do_generate_seeds, refs_are_ctf_corrected, do_sgd);

    	// Check consistency of EMDL_CTF_MAGNIFICATION and MEBL_CTF_DETECTOR_PIXEL_SIZE with mymodel.pixel_size
    	RFLOAT mag, dstep, first_angpix, my_angpix;
    	bool has_magn = false;
    	if (mydata.MDimg.containsLabel(EMDL_CTF_MAGNIFICATION) && mydata.MDimg.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
    	{
    		FOR_ALL_OBJECTS_IN_METADATA_TABLE(mydata.MDimg)
			{
    			mydata.MDimg.getValue(EMDL_CTF_MAGNIFICATION, mag);
    			mydata.MDimg.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
    			my_angpix = 10000. * dstep / mag;
    			if (!has_magn)
    			{
    				first_angpix = my_angpix;
    				has_magn = true;
    			}
    			else if (ABS(first_angpix - my_angpix) > 0.01)
    			{
    				std::cerr << " first_angpix= " << first_angpix << " my_angpix= " << my_angpix << " mag= " << mag  << " dstep= " << dstep << std::endl;
    				REPORT_ERROR("MlOptimiser::initialiseGeneral: ERROR inconsistent magnification and detector pixel sizes in images in input STAR file");
    			}
			}

    	}
    	if (has_magn && ABS(first_angpix - mymodel.pixel_size) > 0.01)
    	{
    		if (verb > 0 && mymodel.pixel_size > 0.)
    			std::cout << "MlOptimiser::initialiseGeneral: WARNING modifying pixel size from " << mymodel.pixel_size <<" to "<<first_angpix << " based on magnification information in the input STAR file" << std::endl;
    		mymodel.pixel_size = first_angpix;
    	}

    	if (!has_magn && mymodel.pixel_size < 0.)
    		REPORT_ERROR("ERROR: you did not specify the pixel size. Use the --angpix option to do so.");
	}
	// Expand movies if fn_data_movie is given AND we were not doing expanded movies already
	else if (fn_data_movie != "" && !do_realign_movies)
	{

		do_realign_movies = true;
		do_parallel_disc_io = false;
		nr_iter_wo_resol_gain = -1;
		nr_iter_wo_large_hidden_variable_changes = 0;
		smallest_changes_optimal_offsets = 999.;
		smallest_changes_optimal_orientations = 999.;
		current_changes_optimal_orientations = 999.;
		current_changes_optimal_offsets = 999.;

		// If we're realigning movie frames, then now read in the metadata of the movie frames and combine with the metadata of the average images
		if (verb > 0)
			std::cout << " Reading in pre-expanded data model for movie frames... " << std::endl;

		mydata.read(fn_data_movie, false, false, do_preread_images);

		// The group numbering might be different: re-assign groups based on group_names
		mymodel.reassignGroupsForMovies(mydata, movie_identifier);

		// Re-normalise the power spectra of the noise, as running averages of movie frames will have more noise than the average micrographs!
		int nframes;
		mydata.MDimg.getValue(EMDL_PARTICLE_NR_FRAMES, nframes, 0);
		int avg_frames;
		mydata.MDimg.getValue(EMDL_PARTICLE_NR_FRAMES_AVG, avg_frames, 0);

		// Also, set the width of the running averages inside the MDimg table of exp_model, for convenient hand-over to polishing
		int ravgwidth = 2 * movie_frame_running_avg_side + 1;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(mydata.MDimg)
		{
			mydata.MDimg.setValue(EMDL_PARTICLE_MOVIE_RUNNING_AVG, ravgwidth);
		}

		// Correct the input sigma2noise spectra by a factor of nframes
		for (int i=0; i < mymodel.nr_groups; i++)
		{
			mymodel.sigma2_noise[i] *= (RFLOAT)nframes/((RFLOAT)(avg_frames * 2 * movie_frame_running_avg_side + 1));
		}

		// Don't do norm correction for realignment of movies.
		do_norm_correction = false;

	}

	if (mymodel.nr_bodies > 1 && !do_auto_refine)
		REPORT_ERROR("ERROR: One cannot use multi-body refinement outside the auto-refine procedure!");

	if (mymodel.nr_classes > 1 && do_split_random_halves)
		REPORT_ERROR("ERROR: One cannot use --split_random_halves with more than 1 reference... You could first classify, and then refine each class separately using --random_halves.");

	if (do_join_random_halves && !do_split_random_halves)
		REPORT_ERROR("ERROR: cannot join random halves because they were not split in the previous run");

	// Local symmetry operators
	fn_local_symmetry_masks.clear();
	fn_local_symmetry_operators.clear();
	if (fn_local_symmetry != "None")
		readRelionFormatMasksAndOperators(fn_local_symmetry, fn_local_symmetry_masks, fn_local_symmetry_operators, mymodel.pixel_size, false);

	// Jun09, 2015 - Shaoda, Helical refinement
	if (do_helical_refine)
	{
		if (mymodel.ref_dim != 3)
			REPORT_ERROR("ERROR: cannot do 2D helical refinement!");

		if ( (do_realign_movies) && (!do_skip_rotate) )
			REPORT_ERROR("ERROR: can only realign movies without rotational searches for helices!");

		if (mymodel.nr_bodies != 1)
			REPORT_ERROR("ERROR: cannot do multi-body refinement for helices!");

		if ( (!do_gpu) && (do_shifts_onthefly) && (!ignore_helical_symmetry) && (verb > 0) )
		{
			std::cerr << " WARNING: On-the-fly shifts slow down helical reconstructions with CPUs considerably. "
					<< "Enable this option only if limited RAM causes trouble (e.g. too large segment boxes used or in 3D sub-tomogram averaging). "
					<< std::endl;
		}

		if ( (helical_keep_tilt_prior_fixed) && (!(helical_sigma_distance < 0.)) )
			REPORT_ERROR("ERROR: cannot keep tilt priors fixed while doing local averaging of helical segments along the same filaments!");

		// Set particle diameter to 90% the box size if user does not give this parameter
		if (particle_diameter < 0.)
		{
			if (verb > 0)
				std::cout << " Automatically set particle diameter to 90% the box size for 3D helical reconstruction." << std::endl;

			particle_diameter = ((RFLOAT)(mymodel.ori_size));
			if ( (((RFLOAT)(mymodel.ori_size)) * 0.05) < width_mask_edge)
				particle_diameter -= 2. * width_mask_edge;
			particle_diameter *= (0.90 * mymodel.pixel_size);
		}

		if (ignore_helical_symmetry)
			mymodel.helical_nr_asu = 1;

		if ( (!do_helical_symmetry_local_refinement) || (ignore_helical_symmetry) )
		{
			mymodel.helical_twist_min = mymodel.helical_twist_max = helical_twist_initial;
			mymodel.helical_rise_min = mymodel.helical_rise_max = helical_rise_initial;
		}

		checkParametersFor3DHelicalReconstruction(
				ignore_helical_symmetry,
				do_helical_symmetry_local_refinement,
				mymodel.helical_nr_asu,
				helical_rise_initial,
				mymodel.helical_rise_min,
				mymodel.helical_rise_max,
				helical_twist_initial,
				mymodel.helical_twist_min,
				mymodel.helical_twist_max,
				mymodel.ori_size,
				mymodel.pixel_size,
				helical_z_percentage,
				particle_diameter,
				helical_tube_inner_diameter,
				helical_tube_outer_diameter);
	}
	else
	{
		if (do_helical_symmetry_local_refinement)
			REPORT_ERROR("ERROR: cannot do local refinement of helical parameters for non-helical segments!");
	}
	if ( (mymodel.ref_dim == 2) && (helical_tube_outer_diameter > 0.) && (particle_diameter < helical_tube_outer_diameter) )
		REPORT_ERROR("ERROR: 2D classification: Helical tube diameter should be smaller than particle diameter!");

	if (do_always_join_random_halves)
		std::cout << " Joining half-reconstructions at each iteration: this is a developmental option to test sub-optimal FSC usage only! " << std::endl;

	// If fn_tau is provided, read in the tau spectrum
	fix_tau = false;
	if (fn_tau != "None")
	{
		fix_tau = true;
		mymodel.readTauSpectrum(fn_tau, verb);
	}

	if (do_auto_refine)
	{
		nr_iter = 999;
		has_fine_enough_angular_sampling = false;
		has_converged = false;

		if (mymodel.tau2_fudge_factor > 1. && verb > 0)
		{
			std::cerr << " WARNING: Using tau2_fudge of " <<mymodel.tau2_fudge_factor << ", which will lead to inflated resolution estimates during refinement." << std::endl;
			std::cerr << " WARNING: You have to run postprocessing afterwards to get a reliable, gold-standard resolution estimate! " << std::endl;
		}

		if (iter == 0 && sampling.healpix_order >= autosampling_hporder_local_searches)
		{
			mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
			sampling.orientational_prior_mode = PRIOR_ROTTILT_PSI;
			RFLOAT rottilt_step = sampling.getAngularSampling(adaptive_oversampling);
			mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 2. * 2. * rottilt_step * rottilt_step;
		}

		// If this is a normal (non-movie) auto-refinement run: check whether we had converged already
		// Jobs often fail in the last iteration, these lines below allow restarting of only the last iteration
		if (do_auto_refine && iter > 0 && !do_realign_movies)
		{
			if (do_force_converge)
			{
				has_converged = true;
				do_join_random_halves = true;
				// In the last iteration, include all data until Nyquist
				do_use_all_data = true;
			}
		}

	}

	// Initialise the sampling object (sets prior mode and fills translations and rotations inside sampling object)
	// May06,2015 - Shaoda & Sjors, initialise for helical translations
	bool do_local_searches = ((do_auto_refine) && (sampling.healpix_order >= autosampling_hporder_local_searches));
	sampling.initialise(mymodel.orientational_prior_mode, mymodel.ref_dim, (mymodel.data_dim == 3), do_gpu, (verb>0),
			do_local_searches, (do_helical_refine) && (!ignore_helical_symmetry), helical_rise_initial / mymodel.pixel_size, helical_twist_initial);

	// Now that sampling is initialised, also modify sigma2_rot for the helical refinement
        if (do_auto_refine && do_helical_refine && !ignore_helical_symmetry && iter == 0 && sampling.healpix_order >= autosampling_hporder_local_searches)
	{
		// Aug20,2015 - Shaoda, Helical refinement
		RFLOAT rottilt_step = sampling.getAngularSampling(adaptive_oversampling);
		mymodel.sigma2_rot = getHelicalSigma2Rot((helical_rise_initial / mymodel.pixel_size), helical_twist_initial, sampling.helical_offset_step, rottilt_step, mymodel.sigma2_rot);
	}

	// Default max_coarse_size is original size
	if (max_coarse_size < 0)
		max_coarse_size = mymodel.ori_size;

	if (particle_diameter < 0.)
    	particle_diameter = (mymodel.ori_size - width_mask_edge) * mymodel.pixel_size;

    // For do_average_unaligned, always use initial low_pass filter
    if (do_average_unaligned && ini_high < 0.)
    {
    	// By default, use 0.07 dig.freq. low-pass filter
    	// See S.H.W. Scheres (2010) Meth Enzym.
    	ini_high = 1./mymodel.getResolution(ROUND(0.07 * mymodel.ori_size));
    }

	// For skipped alignments
	// Also do not perturb this orientation, nor do oversampling or priors
	// Switch off on-the-fly shifts, as that will not work when skipping alignments! (it isn't necessary anyway in that case)
    if (do_skip_align || do_skip_rotate)
	{
		mymodel.orientational_prior_mode = NOPRIOR;
		sampling.orientational_prior_mode = NOPRIOR;
		adaptive_oversampling = 0;
		sampling.perturbation_factor = 0.;
		sampling.random_perturbation = 0.;
		sampling.addOneOrientation(0.,0.,0., true);
		directions_have_changed = true;
	}
    else if (do_only_sample_tilt)
    {
		std::cout << " Only sampling tilt, keep rot and psi fixed." << std::endl;

    	mymodel.orientational_prior_mode = NOPRIOR;
		sampling.orientational_prior_mode = NOPRIOR;
		adaptive_oversampling = 0;
		sampling.perturbation_factor = 0.;
		sampling.random_perturbation = 0.;

		sampling.directions_ipix.clear();
		sampling.rot_angles.clear();
		sampling.tilt_angles.clear();
		RFLOAT rot = 0., psi = 0., tilt;
		int ipix = 0;
		for (tilt = -180.; tilt < 180.; tilt+= sampling.getAngularSampling(), ipix++)
		{
			sampling.rot_angles.push_back(rot);
			sampling.tilt_angles.push_back(tilt);
			sampling.directions_ipix.push_back(ipix);
		}
		sampling.psi_angles.clear();
		sampling.psi_angles.push_back(psi);
		directions_have_changed = true;
    }
    if (do_skip_align)
	{
		RFLOAT dummy=0.;
		sampling.addOneTranslation(dummy, dummy, dummy, true);
		do_shifts_onthefly = false; // on-the-fly shifts are incompatible with do_skip_align!
	}
    if ( (do_bimodal_psi) && (mymodel.sigma2_psi > 0.) && (verb > 0) )
    	std::cout << " Using bimodal priors on the psi angle..." << std::endl;

	// Resize the pdf_direction arrays to the correct size and fill with an even distribution
	if (directions_have_changed)
		mymodel.initialisePdfDirection(sampling.NrDirections());

	// Initialise the wsum_model according to the mymodel
	wsum_model.initialise(mymodel, sampling.symmetryGroup(), asymmetric_padding, skip_gridding);

	// Initialise sums of hidden variable changes
	// In later iterations, this will be done in updateOverallChangesInHiddenVariables
	sum_changes_optimal_orientations = 0.;
	sum_changes_optimal_offsets = 0.;
	sum_changes_optimal_classes = 0.;
	sum_changes_count = 0.;

	if (mymodel.data_dim == 3)
	{

		if (do_gpu)
			std::cerr << "NOTE: Input 3D data on GPU is development feature: please report unexpected behaviour" << std::endl;

		// TODO: later do norm correction?!
		// Don't do norm correction for volume averaging at this stage....
		do_norm_correction = false;

		if (!((do_helical_refine) && (!ignore_helical_symmetry))) // For 3D helical sub-tomogram averaging, either is OK, so let the user decide
			do_shifts_onthefly = true; // save RAM for volume data (storing all shifted versions would take a lot!)

		if (do_skip_align)
			do_shifts_onthefly = false; // on-the-fly shifts are incompatible with do_skip_align!
		// getMetaAndImageData is not made for passing multiple volumes!
		do_parallel_disc_io = true;
	}
	// Tabulated sine and cosine values (for 2D helical segments / 3D helical sub-tomogram averaging with on-the-fly shifts)
	if ( (!do_gpu) && (do_helical_refine) && (!ignore_helical_symmetry) && (do_shifts_onthefly) )
	{
		tab_sin.initialise(100000);
		tab_cos.initialise(100000);
	}

	// Skip scale correction if there are no groups
	if (mymodel.nr_groups == 1 && !do_realign_movies)
		do_scale_correction = false;

	// Check for rlnReconstructImageName in the data.star file. If it is present, set do_use_reconstruct_images to true
	do_use_reconstruct_images = mydata.MDimg.containsLabel(EMDL_IMAGE_RECONSTRUCT_NAME);
	if (do_use_reconstruct_images && verb > 0)
		std::cout <<" Using rlnReconstructImageName from the input data.star file!" << std::endl;

	// For new thread-parallelization: each thread does 1 particle, so nr_pool=nr_threads
	nr_pool = x_pool*nr_threads;

	// Check image size model and data are the same
	int data_image_size;
	mydata.MDexp.getValue(EMDL_IMAGE_SIZE, data_image_size);
	if (XSIZE(mymodel.Iref[0]) != data_image_size)
		REPORT_ERROR("ERROR: reference and data image sizes are not the same!");

	// Make subsets?
	nr_subsets = 1;
	if (subset_size > 0)
	{
    	//do_norm_correction = false;
	    //do_scale_correction = false;
	    if (random_seed != 0)
    		mydata.randomiseOriginalParticlesOrder(random_seed);
	    nr_subsets = mydata.numberOfOriginalParticles() / subset_size;
	    if (write_every_subset < 0)
	    	write_every_subset = nr_subsets;
    }
	else
	{
	    mu = 0.;
	}

#ifdef DEBUG
	std::cerr << "Leaving initialiseGeneral" << std::endl;
#endif

}

void MlOptimiser::initialiseWorkLoad()
{

	// Note, this function is overloaded in ml_optimiser_mpi (where random_seed is only set by the master and then send to all slaves!)

	// Randomise the order of the particles
	if (random_seed == -1) random_seed = time(NULL);
    // Also randomize random-number-generator for perturbations on the angles
    init_random_generator(random_seed);

    divide_equally(mydata.numberOfOriginalParticles(), 1, 0, my_first_ori_particle_id, my_last_ori_particle_id);

    // Now copy particle stacks to scratch if needed
    if (fn_scratch != "" && !do_preread_images && !do_reuse_scratch)
    {
    	mydata.prepareScratchDirectory(fn_scratch);
    	bool also_do_ctfimage = (mymodel.data_dim == 3 && do_ctf_correction);
    	mydata.copyParticlesToScratch(1, true, also_do_ctfimage, keep_free_scratch_Gb);
    }

}

void MlOptimiser::calculateSumOfPowerSpectraAndAverageImage(MultidimArray<RFLOAT> &Mavg, bool myverb)
{

#ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::calculateSumOfPowerSpectraAndAverageImage Entering"<<std::endl;
#endif

    int barstep, my_nr_ori_particles = my_last_ori_particle_id - my_first_ori_particle_id + 1;
    if (my_nr_ori_particles < 1)
    {
    	// Master doesn't do anything here...
    	// But still set Mavg the right size for AllReduce later on
    	FileName fn_img;
    	mydata.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, 0);
    	Image<RFLOAT> img;
    	img.read(fn_img, false); // don't read data
    	Mavg.initZeros(img());
    	return;
    }

    if (myverb > 0)
	{
		std::cout << " Estimating initial noise spectra " << std::endl;
		init_progress_bar(my_nr_ori_particles);
		barstep = XMIPP_MAX(1, my_nr_ori_particles / 60);
	}

    // Only open stacks once and then read multiple images
	fImageHandler hFile;
	long int dump;
	FileName fn_open_stack="";

	// Note the loop over the particles (part_id) is MPI-parallelized
	int nr_ori_particles_done = 0;
	FileName fn_img, fn_stack;
	// For spectrum calculation: recycle the transformer (so do not call getSpectrum all the time)
	MultidimArray<Complex > Faux;
    FourierTransformer transformer;
	MetaDataTable MDimg;

	// Start reconstructions at ini_high or 0.07 digital frequencies....
	if (ini_high <= 0.)
		wsum_model.current_size = 1./mymodel.getResolution(ROUND(0.07 * mymodel.ori_size));
	else
		wsum_model.current_size  = mymodel.getPixelFromResolution(1./ini_high);
	wsum_model.initZeros();

	for (long int ori_part_id = my_first_ori_particle_id; ori_part_id <= my_last_ori_particle_id; ori_part_id++, nr_ori_particles_done++)
	{

		for (long int i = 0; i < mydata.ori_particles[ori_part_id].particles_id.size(); i++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[i];
			long int group_id = mydata.getGroupId(part_id);

			// May24,2015 - Shaoda & Sjors, Helical refinement
			RFLOAT psi_deg = 0., tilt_deg = 0.;
			bool is_helical_segment = (do_helical_refine) || ((mymodel.ref_dim == 2) && (helical_tube_outer_diameter > 0.));

			// Extract the relevant MetaDataTable row from MDimg
			MDimg = mydata.getMetaDataImage(part_id);

			if (!mydata.getImageNameOnScratch(part_id, fn_img))
				MDimg.getValue(EMDL_IMAGE_NAME, fn_img);

			// May24,2015 - Shaoda & Sjors, Helical refinement
			if (is_helical_segment)
			{
				if (!MDimg.getValue(EMDL_ORIENT_PSI_PRIOR, psi_deg))
				{
					if (!MDimg.getValue(EMDL_ORIENT_PSI, psi_deg))
						REPORT_ERROR("ml_optimiser.cpp::calculateSumOfPowerSpectraAndAverageImage: Psi priors of helical segments are missing!");
				}
				if (!MDimg.getValue(EMDL_ORIENT_TILT_PRIOR, tilt_deg))
				{
					if (!MDimg.getValue(EMDL_ORIENT_TILT, tilt_deg))
						REPORT_ERROR("ml_optimiser.cpp::calculateSumOfPowerSpectraAndAverageImage: Tilt priors of helical segments are missing!");
				}
			}

			// Read image from disc
            Image<RFLOAT> img;
			if (do_preread_images && do_parallel_disc_io)
            {
                img().reshape(mydata.particles[part_id].img);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mydata.particles[part_id].img)
				{
                	DIRECT_MULTIDIM_ELEM(img(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(mydata.particles[part_id].img, n);
				}
            }
            else
            {
				fn_img.decompose(dump, fn_stack);
				if (fn_stack != fn_open_stack)
				{
					hFile.openFile(fn_stack, WRITE_READONLY);
					fn_open_stack = fn_stack;
				}
				img.readFromOpenFile(fn_img, hFile, -1, false);
                img().setXmippOrigin();
            }

			// May24,2015 - Shaoda & Sjors, Helical refinement
			// Check that the average in the noise area is approximately zero and the stddev is one
			if (!dont_raise_norm_error)
			{
				// NEW METHOD
				RFLOAT sum, sum2, sphere_radius_pix, cyl_radius_pix;
				cyl_radius_pix = helical_tube_outer_diameter / (2. * mymodel.pixel_size);
				sphere_radius_pix = particle_diameter / (2. * mymodel.pixel_size);
				calculateBackgroundAvgStddev(img, sum, sum2, (int)(ROUND(sphere_radius_pix)), is_helical_segment, cyl_radius_pix, tilt_deg, psi_deg);

				// Average should be close to zero, i.e. max +/-50% of stddev...
				// Stddev should be close to one, i.e. larger than 0.5 and smaller than 2)
				if (ABS(sum/sum2) > 0.5 || sum2 < 0.5 || sum2 > 2.0)
				{
					std::cerr << " fn_img= " << fn_img << " bg_avg= " << sum << " bg_stddev= " << sum2 << std::flush;
					if (is_helical_segment)
						std::cerr << " tube_bg_radius= " << cyl_radius_pix << " psi_deg= " << psi_deg << " tilt_deg= " << tilt_deg << " (this is a particle from a helix)" << std::flush;
					else
						std::cerr << " bg_radius= " << sphere_radius_pix << std::flush;
					std::cerr << std::endl;
					REPORT_ERROR("ERROR: It appears that these images have not been normalised to an average background value of 0 and a stddev value of 1. \n \
							Note that the average and stddev values for the background are calculated: \n \
							(1) for single particles: outside a circle with the particle diameter \n \
							(2) for helical segments: outside a cylinder (tube) with the helical tube diameter \n \
							You can use the relion_preprocess program to normalise your images \n \
							If you are sure you have normalised the images correctly (also see the RELION Wiki), you can switch off this error message using the --dont_check_norm command line option");
				}
			}

			// Apply a similar softMask as below (assume zero translations)
			if (do_zero_mask)
			{
				// May24,2015 - Shaoda & Sjors, Helical refinement
				if (is_helical_segment)
				{
					softMaskOutsideMapForHelix(img(), psi_deg, tilt_deg, (particle_diameter / (2. * mymodel.pixel_size)),
							(helical_tube_outer_diameter / (2. * mymodel.pixel_size)), width_mask_edge);
				}
				else
					softMaskOutsideMap(img(), particle_diameter / (2. * mymodel.pixel_size), width_mask_edge);
			}

			// Keep track of the average image (only to correct power spectra, no longer for initial references!)
			if (part_id == mydata.ori_particles[my_first_ori_particle_id].particles_id[0])
				Mavg = img();
			else
				Mavg += img();

			// Calculate the power spectrum of this particle
			CenterFFT(img(), true);
   			transformer.FourierTransform(img(), Faux);
   			MultidimArray<RFLOAT> ind_spectrum, count;
   			ind_spectrum.initZeros(XSIZE(img()));
			count.initZeros(XSIZE(img()));
			// recycle the same transformer for all images
			transformer.FourierTransform(img(), Faux, false);
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
			{
				long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
				ind_spectrum(idx) += norm(dAkij(Faux, k, i, j));
				count(idx) += 1.;
			}
			ind_spectrum /= count;

			// Resize the power_class spectrum to the correct size and keep sum
			ind_spectrum.resize(wsum_model.sigma2_noise[0]); // Store sum of all groups in group 0
			wsum_model.sigma2_noise[group_id] += ind_spectrum;
			wsum_model.sumw_group[group_id] += 1.;

			if (fn_ref == "None")
			{

				MultidimArray<RFLOAT> Fctf, Fweight;
				MultidimArray<Complex > Fimg;

				// Make sure MPI and sequentialo behave exactly the same
				init_random_generator(random_seed + part_id);
				// Randomize the initial orientations for initial reference generation at this step....
				// TODO: this is not an even angular distribution....
				RFLOAT rot  = (mymodel.ref_dim == 2) ? 0. : rnd_unif() * 360.;
				RFLOAT tilt = (mymodel.ref_dim == 2) ? 0. :rnd_unif() * 180.;
				RFLOAT psi  = rnd_unif() * 360.;
				int iclass  = rnd_unif() * mymodel.nr_classes;
				Matrix2D<RFLOAT> A;
				Euler_angles2matrix(rot, tilt, psi, A, true);

				// Construct initial references from random subsets
	   			windowFourierTransform(Faux, Fimg, wsum_model.current_size);
				Fctf.resize(Fimg);
				Fctf.initConstant(1.);
				// Apply CTF if necessary (skip this for subtomograms!)
				if (do_ctf_correction && mymodel.data_dim != 3)
				{
					CTF ctf;
					ctf.read(MDimg, MDimg);
					ctf.getFftwImage(Fctf, mymodel.ori_size, mymodel.ori_size, mymodel.pixel_size,
						ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
					{
						DIRECT_MULTIDIM_ELEM(Fimg, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}
				(wsum_model.BPref[iclass]).set2DFourierTransform(Fimg, A, IS_NOT_INV, &Fctf);
			}

			// For sub-tomogram averaging: take effect of rotations on stddev of images into account here
			// TODO: CHECK THAT NOW THAT WE USE BPref TO CALCULATE AVERAGE: THIS IS NO LONGER NECESSARY?!!!
			//if (mymodel.data_dim == 3)
			//{
			//	RFLOAT stddev1 = img().computeStddev();
			//	selfApplyGeometry(img(), A, IS_INV, WRAP);
			//	RFLOAT stddev2 = img().computeStddev();
			//	// Correct for interpolation errors that drive down the average density...
			//	img() *= stddev1 / stddev2;
			//}


		} // end loop part_id (i)

		if (myverb > 0 && nr_ori_particles_done % barstep == 0)
			progress_bar(nr_ori_particles_done);

	} // end loop ori_part_id


	// Clean up the fftw object completely
	// This is something that needs to be done manually, as among multiple threads only one of them may actually do this
	transformer.cleanup();

	if (myverb > 0)
		progress_bar(my_nr_ori_particles);

#ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::calculateSumOfPowerSpectraAndAverageImage Leaving"<<std::endl;
#endif

}

void MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage(MultidimArray<RFLOAT> &Mavg)
{

#ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage Entering"<<std::endl;
#endif

	// First calculate average image
	RFLOAT total_sum = 0.;
    for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
    {
    	mymodel.nr_particles_group[igroup] = ROUND(wsum_model.sumw_group[igroup]);
    	total_sum += wsum_model.sumw_group[igroup];
    }
    Mavg /= total_sum;

	if (fn_ref == "None")
	{
		for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++)
		{

			MultidimArray<RFLOAT> dummy;
			(wsum_model.BPref[iclass]).reconstruct(mymodel.Iref[iclass], gridding_nr_iter, false,
					1., dummy, dummy, dummy, dummy, dummy);
			// 2D projection data were CTF-corrected, subtomograms were not
			refs_are_ctf_corrected = (mymodel.data_dim == 3) ? false : true;
		}
	}

	// Calculate sigma2_noise estimates as average of power class spectra, and subtract power spectrum of the average image from that
	if (do_calculate_initial_sigma_noise)
	{
		// Calculate power spectrum of the average image
		MultidimArray<RFLOAT> spect;
		getSpectrum(Mavg, spect, POWER_SPECTRUM);
		spect /= 2.; // because of 2-dimensionality of the complex plane
		spect.resize(mymodel.sigma2_noise[0]);

		for (int igroup = 0; igroup < wsum_model.nr_groups; igroup++)
		{
			// Factor 2 because of 2-dimensionality of the complex plane
			mymodel.sigma2_noise[igroup] = wsum_model.sigma2_noise[igroup] / ( 2. * wsum_model.sumw_group[igroup] );

			// Now subtract power spectrum of the average image from the average power spectrum of the individual images
			mymodel.sigma2_noise[igroup] -= spect;

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(spect)
			{
				// Remove any negative sigma2_noise values: replace by positive neighbouring value
				if (DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n) < 0. )
				{
					// First try the previous value
					if (n - 1 >= 0 && DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n - 1) > 0.)
					{
						DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n) = DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n - 1);
					}
					else
					{
						bool is_positive = false;
						int nn = n;
						while (!is_positive)
						{
							nn++;
							if (nn > XSIZE(mymodel.sigma2_noise[igroup]))
							{
								std::cerr << " igroup= " << igroup << " n= " << n << " mymodel.sigma2_noise[igroup]= " << mymodel.sigma2_noise[igroup] << std::endl;
								REPORT_ERROR("BUG! cannot find positive values in sigma2_noise spectrum");
							}
							if (DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], nn) > 0.)
							{
								is_positive = true;
								DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n) = DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], nn);
							}
						}
					}
				}
			}

		}
	}

#ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage Leaving"<<std::endl;
#endif

}

void MlOptimiser::initialLowPassFilterReferences()
{
	if (ini_high > 0.)
	{

		// Make a soft (raised cosine) filter in Fourier space to prevent artefacts in real-space
		// The raised cosine goes through 0.5 at the filter frequency and has a width of width_mask_edge fourier pixels
		RFLOAT radius = mymodel.ori_size * mymodel.pixel_size / ini_high;
		radius -= WIDTH_FMASK_EDGE / 2.;
		RFLOAT radius_p = radius + WIDTH_FMASK_EDGE;
		FourierTransformer transformer;
		MultidimArray<Complex > Faux;
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		{
			transformer.FourierTransform(mymodel.Iref[iclass], Faux);
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
			{
				RFLOAT r = sqrt((RFLOAT)(kp*kp + ip*ip + jp*jp));
				if (r < radius)
					continue;
				else if (r > radius_p)
					DIRECT_A3D_ELEM(Faux, k, i, j) = 0.;
				else
				{
					DIRECT_A3D_ELEM(Faux, k, i, j) *= 0.5 - 0.5 * cos(PI * (radius_p - r) / WIDTH_FMASK_EDGE);
				}
			}
			transformer.inverseFourierTransform(Faux, mymodel.Iref[iclass]);
		}

	}

}

/** ========================== EM-Iteration  ================================= */

void MlOptimiser::iterateSetup()
{

	// Make a barrier where all working threads wait
	global_barrier = new Barrier(nr_threads - 1);

    // Create threads to start working
	global_ThreadManager = new ThreadManager(nr_threads, this);

	// Set up the thread task distributors for the particles and the orientations (will be resized later on)
	exp_ipart_ThreadTaskDistributor = new ThreadTaskDistributor(nr_threads, 1);

}
void MlOptimiser::iterateWrapUp()
{

	// delete barrier, threads and task distributors
    delete global_barrier;
	delete global_ThreadManager;
    delete exp_ipart_ThreadTaskDistributor;

    // Delete volatile space on scratch
    mydata.deleteDataOnScratch();

}

void MlOptimiser::iterate()
{

	if (do_split_random_halves)
		REPORT_ERROR("ERROR: Cannot split data into random halves without using MPI!");


	// launch threads etc
	iterateSetup();

	// Update the current resolution and image sizes, and precalculate resolution pointers
	// The rest of the time this will be done after maximization and before writing output files,
	// so that current resolution is in the output files of the current iteration
	updateCurrentResolution();

	// If we're doing a restart from subsets, then do not increment the iteration number in the restart!
	if (subset > 0)
		iter--;

	bool has_already_reached_convergence = false;
	for (iter = iter + 1; iter <= nr_iter; iter++)
    {
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
		// May18,2015 - Shaoda & Sjors, Helical refinement (orientational searches)
		std::cerr << std::endl << std::endl;
		std::cerr << "MlOptimiser::iterate()" << std::endl;
		std::cerr << "iter = " << iter << ", do_helical_refine == " << std::flush;
		if(do_helical_refine)
			std::cerr << "true" << std::flush;
		else
			std::cerr << "false" << std::flush;
		std::cerr << std::endl;
#endif

		for (subset = subset_start; subset <= nr_subsets; subset++)
		{

#ifdef TIMING
			timer.tic(TIMING_EXP);
#endif

			if (do_auto_refine)
			{
				// Check whether we have converged by now
				// If we have, set do_join_random_halves and do_use_all_data for the next iteration
				checkConvergence();
			}

			expectation();

			int old_verb = verb;
			if (nr_subsets > 1) // be quiet
				verb = 0;

			// Sjors & Shaoda Apr 2015
			// This function does enforceHermitianSymmetry, applyHelicalSymmetry and applyPointGroupSymmetry sequentially.
			// First it enforces Hermitian symmetry to the back-projected Fourier 3D matrix.
			// Then helical symmetry is applied in Fourier space. It does rise and twist for all asymmetrical units in Fourier space.
			// Finally it applies point group symmetry (such as Cn, ...).
			// DEBUG
			if (verb > 0)
			{
				if ( (do_helical_refine) && (!ignore_helical_symmetry) )
				{
					if (mymodel.helical_nr_asu > 1)
						std::cout << " Applying helical symmetry from the last iteration for all asymmetrical units in Fourier space..." << std::endl;
					if ( (iter > 1) && (do_helical_symmetry_local_refinement) )
					{
						std::cout << " Refining helical symmetry in real space..." << std::endl;
						std::cout << " Applying refined helical symmetry in real space..." << std::endl;
					}
					else
						std::cout << " Applying helical symmetry from the last iteration in real space..." << std::endl;
				}
			}
			symmetriseReconstructions();

#ifdef TIMING
			timer.toc(TIMING_EXP);
			timer.tic(TIMING_MAX);
#endif

			if (do_skip_maximization)
			{
				// Only write data.star file and break from the iteration loop
				write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, 0);
				break;
			}

			maximization();

#ifdef TIMING
			timer.toc(TIMING_MAX);
#endif

			// Apply local symmetry according to a list of masks and their operators
			applyLocalSymmetryForEachRef();

			// Shaoda Jul26,2015
			// Helical symmetry refinement and imposition of real space helical symmetry.
			if (do_helical_refine)
			{
				if (!ignore_helical_symmetry)
					makeGoodHelixForEachRef();
				if ( (!do_skip_align) && (!do_skip_rotate) )
				{
					int nr_same_polarity = 0, nr_opposite_polarity = 0;
					RFLOAT opposite_percentage = 0.;
					bool do_auto_refine_local_searches = (do_auto_refine) && (sampling.healpix_order >= autosampling_hporder_local_searches);
					bool do_classification_local_searches = (!do_auto_refine) && (mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI)
							&& (mymodel.sigma2_rot > 0.) && (mymodel.sigma2_tilt > 0.) && (mymodel.sigma2_psi > 0.);
					bool do_local_angular_searches = (do_auto_refine_local_searches) || (do_classification_local_searches);

					if (helical_sigma_distance < 0.)
						updateAngularPriorsForHelicalReconstruction(mydata.MDimg, helical_keep_tilt_prior_fixed);
					else
					{
						updatePriorsForHelicalReconstruction(
								mydata.MDimg,
								nr_opposite_polarity,
								helical_sigma_distance * ((RFLOAT)(mymodel.ori_size)),
								(mymodel.data_dim == 3),
								do_auto_refine,
								do_local_angular_searches,
								mymodel.sigma2_rot,
								mymodel.sigma2_tilt,
								mymodel.sigma2_psi,
								mymodel.sigma2_offset,
								helical_keep_tilt_prior_fixed);

						nr_same_polarity = ((int)(mydata.MDimg.numberOfObjects())) - nr_opposite_polarity;
						opposite_percentage = (100.) * ((RFLOAT)(nr_opposite_polarity)) / ((RFLOAT)(mydata.MDimg.numberOfObjects()));
						if ( (verb > 0) && (!do_local_angular_searches) )
						{
							//std::cout << " DEBUG: auto_refine, healpix_order, min_for_local = " << do_auto_refine << ", " << sampling.healpix_order << ", " << autosampling_hporder_local_searches << std::endl;
							//std::cout << " DEBUG: orient_prior_mode = " << PRIOR_ROTTILT_PSI << ", sigma_ang2 = " << mymodel.sigma2_rot << ", " << mymodel.sigma2_tilt << ", " << mymodel.sigma2_psi << ", sigma_offset2 = " << mymodel.sigma2_offset << std::endl;
							std::cout << " Number of helical segments with psi angles similar/opposite to their priors: " << nr_same_polarity << " / " << nr_opposite_polarity << " (" << opposite_percentage << "%)" << std::endl;
						}
					}
				}
			}

			// Apply masks to the reference images
			// At the last iteration, do not mask the map for validation purposes
			if (do_solvent && !has_converged)
				solventFlatten();

			// Re-calculate the current resolution, do this before writing to get the correct values in the output files
			updateCurrentResolution();

			// Write output files
			write(DO_WRITE_SAMPLING, DO_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, 0);

			if (do_auto_refine && has_converged)
			{
				if (verb > 0)
				{
					std::cout << " Auto-refine: Refinement has converged, stopping now... " << std::endl;
					std::cout << " Auto-refine: + Final reconstruction from all particles is saved as: " <<  fn_out << "_class001.mrc" << std::endl;
					std::cout << " Auto-refine: + Final model parameters are stored in: " << fn_out << "_model.star" << std::endl;
					std::cout << " Auto-refine: + Final data parameters are stored in: " << fn_out << "_data.star" << std::endl;
                                        if (mymodel.tau2_fudge_factor > 1.)
                                        {
                                            std::cout << " Auto-refine: + SEVERE WARNING: Because you used a tau2_fudge of " << mymodel.tau2_fudge_factor << " your resolution during this refinement will be inflated!" << std::endl;
                                            std::cout << " Auto-refine: + SEVERE WARNING: You have to run a postprocessing on the unfil.mrc maps to get a gold-standard resolution estimate!"  << std::endl;
                                        }
                                        else if (do_phase_random_fsc)
                                            std::cout << " Auto-refine: + Final resolution (already with masking) is: " << 1./mymodel.current_resolution << std::endl;
                                        else
                                        {
                                            std::cout << " Auto-refine: + Final resolution (without masking) is: " << 1./mymodel.current_resolution << std::endl;
                                            std::cout << " Auto-refine: + But you may want to run relion_postprocess to mask the unfil.mrc maps and calculate a higher resolution FSC" << std::endl;
                                        }
					if (acc_rot > 10.)
					{
						std::cout << " Auto-refine: + SEVERE WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles!" << std::endl;
						std::cout << " Auto-refine: + SEVERE WARNING: This has been observed to lead to spurious FSC curves, so be VERY wary of inflated resolution estimates..." << std::endl;
						std::cout << " Auto-refine: + SEVERE WARNING: You most probably do NOT want to publish these results!" << std::endl;
						std::cout << " Auto-refine: + SEVERE WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
					}
				}
				break;
			}

			verb = old_verb;

			if (nr_subsets > 1 && sgd_max_subsets > 0 && subset > sgd_max_subsets)
				break; // break out of loop over the subsets

#ifdef TIMING
			if (verb > 0)
				timer.printTimes(false);
#endif
		} // end loop subsets

		// In the next iteration, start again from the first subset
		subset_start = 1;

		// Stop subsets after sgd_max_subsets has been reached
		if (nr_subsets > 1 && sgd_max_subsets > 0 && subset > sgd_max_subsets)
		{
			// Write out without a _sub in the name
			nr_subsets = 1;
			write(DO_WRITE_SAMPLING, DO_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, 0);


			if (do_sgd)
			{
				if (verb > 0)
					std::cout << " SGD has reached the maximum number of subsets, so stopping now..." << std::endl;
				// For initial model generation, just stop after the sgd_max_subsets has been reached
				break;
			}
			else
			{
				if (verb > 0)
					std::cout << " Run has reached the maximum number of subsets, continuing without subsets now..." << std::endl;
				// For subsets in 2D classification, now continue rest of iterations without subsets in the next iteration
				subset_size = -1;
			}
		}

    } // end loop iters

	// delete threads etc
	iterateWrapUp();
}

void MlOptimiser::expectation()
{

//#define DEBUG_EXP
#ifdef DEBUG_EXP
	std::cerr << "Entering expectation" << std::endl;
#endif

	// Initialise some stuff
	// A. Update current size (may have been changed to ori_size in autoAdjustAngularSampling) and resolution pointers
	updateImageSizeAndResolutionPointers();

	// B. Initialise Fouriertransform, set weights in wsum_model to zero, initialise AB-matrices for FFT-phase shifts, etc
	expectationSetup();

#ifdef DEBUG_EXP
	std::cerr << "Expectation: done setup" << std::endl;
#endif

	// C. Calculate expected minimum angular errors (only for 3D refinements)
	// And possibly update orientational sampling automatically
	if (!((iter==1 && do_firstiter_cc) || do_always_cc) && !(do_skip_align || do_only_sample_tilt) && subset == 1)
	{
		// Set the exp_metadata (but not the exp_imagedata which is not needed for calculateExpectedAngularErrors)
		int n_trials_acc = (mymodel.ref_dim==3 && mymodel.data_dim != 3) ? 100 : 10;
		n_trials_acc = XMIPP_MIN(n_trials_acc, mydata.numberOfOriginalParticles());
		getMetaAndImageDataSubset(0, n_trials_acc-1, false);
		calculateExpectedAngularErrors(0, n_trials_acc-1);
	}

	// D. Update the angular sampling (all nodes except master)
	if ( (do_auto_refine) && iter > 1 )
		updateAngularSampling();

	// E. Check whether everything fits into memory
	if (subset == 1)
		expectationSetupCheckMemory(verb);

	// F. Precalculate AB-matrices for on-the-fly shifts
	// Use tabulated sine and cosine values instead for 2D helical segments / 3D helical sub-tomogram averaging with on-the-fly shifts
	if ( (do_shifts_onthefly) && (subset == 1) && (!((do_helical_refine) && (!ignore_helical_symmetry))) )
		precalculateABMatrices();


#ifdef DEBUG_EXP
	std::cerr << "Expectation: done setupCheckMemory" << std::endl;
#endif

#ifdef CUDA
	/************************************************************************/
	//GPU memory setup

	if (do_gpu)
    {
		for (int i = 0; i < cudaDevices.size(); i ++)
		{
			MlDeviceBundle *b = new MlDeviceBundle(this);
			b->setDevice(cudaDevices[i]);
			b->setupFixedSizedObjects();
			cudaDeviceBundles.push_back((void*)b);
		}

		std::vector<unsigned> threadcountOnDevice(cudaDeviceBundles.size(),0);

		for (int i = 0; i < cudaOptimiserDeviceMap.size(); i ++)
		{
			std::stringstream didSs;
			didSs << "RRt" << i;
			MlOptimiserCuda *b = new MlOptimiserCuda(this, (MlDeviceBundle*) cudaDeviceBundles[cudaOptimiserDeviceMap[i]], didSs.str().c_str());
			b->resetData();
			cudaOptimisers.push_back((void*)b);
			threadcountOnDevice[cudaOptimiserDeviceMap[i]] ++;
		}

		int devCount;
		HANDLE_ERROR(cudaGetDeviceCount(&devCount));
		HANDLE_ERROR(cudaDeviceSynchronize());
		for (int i = 0; i < cudaDeviceBundles.size(); i ++)
		{
			if(((MlDeviceBundle*)cudaDeviceBundles[i])->device_id >= devCount || ((MlDeviceBundle*)cudaDeviceBundles[i])->device_id < 0 )
			{
				//std::cerr << " using device_id=" << ((MlDeviceBundle*)cudaDeviceBundles[i])->device_id << " (device no. " << ((MlDeviceBundle*)cudaDeviceBundles[i])->device_id+1 << ") which is not within the available device range" << devCount << std::endl;
				CRITICAL(ERR_GPUID);
			}
			else
				HANDLE_ERROR(cudaSetDevice(((MlDeviceBundle*)cudaDeviceBundles[i])->device_id));

			size_t free, total, allocationSize;
			HANDLE_ERROR(cudaMemGetInfo( &free, &total ));

			size_t required_free = requested_free_gpu_memory + GPU_THREAD_MEMORY_OVERHEAD_MB*1000*1000*threadcountOnDevice[i];

			if (free < required_free)
			{
				printf("WARNING: Ignoring required free GPU memory amount of %zu MB, due to space insufficiency.\n", required_free/1000000);
				allocationSize = (double)free *0.7;
			}
			else
				allocationSize = free - required_free;

			if (allocationSize < 200000000)
				printf("WARNING: The available space on the GPU after initialization (%zu MB) might be insufficient for the expectation step.\n", allocationSize/1000000);

#ifdef PRINT_GPU_MEM_INFO
			printf("INFO: Projector model size %dx%dx%d\n", (int)mymodel.PPref[0].data.xdim, (int)mymodel.PPref[0].data.ydim, (int)mymodel.PPref[0].data.zdim );
			printf("INFO: Free memory for Custom Allocator of device bundle %d is %d MB\n", i, (int) ( ((float)allocationSize)/1000000.0 ) );
#endif

			((MlDeviceBundle*)cudaDeviceBundles[i])->setupTunableSizedObjects(allocationSize);
		}
	}
#endif

	/************************************************************************/

	// Now perform real expectation over all particles
	// Use local parameters here, as also done in the same overloaded function in MlOptimiserMpi

	int old_verb = verb;
	long int prev_barstep = 0;
	int barstep = XMIPP_MAX(1, mydata.numberOfOriginalParticles() / 60);
	if (nr_subsets > 1)
		barstep = XMIPP_MIN(barstep, subset_size);

	long int my_subset_first_ori_particle, my_subset_last_ori_particle, nr_particles_todo;
	if (nr_subsets > 1)
	{
		divide_equally(mydata.numberOfOriginalParticles(), nr_subsets, subset-1, my_subset_first_ori_particle, my_subset_last_ori_particle);
		subset_size = nr_particles_todo = my_subset_last_ori_particle - my_subset_first_ori_particle + 1;
		if (verb > 0)
		{
			// SGD: progress bar over entire iteration
			if (subset == subset_start)
			{
				if (do_sgd)
					std::cout << " Stochastic Gradient Descent iteration " << iter << " of " << nr_iter << std::endl;
				else
					std::cout << " Incomplete expectation iteration " << iter << " of " << nr_iter << std::endl;
				long int barsize = (sgd_max_subsets > 0) ? sgd_max_subsets * subset_size : mydata.numberOfOriginalParticles();
				barsize = XMIPP_MIN(barsize, mydata.numberOfOriginalParticles());
				init_progress_bar(barsize);
			}
		}
	}
	else
	{
		my_subset_first_ori_particle = 0.;
		my_subset_last_ori_particle = mydata.numberOfOriginalParticles() - 1;
		subset_size = nr_particles_todo = mydata.numberOfOriginalParticles();
		if (verb > 0)
		{
			std::cout << " Expectation iteration " << iter;
			if (!do_auto_refine)
				std::cout << " of " << nr_iter;
			std::cout << std::endl;
			init_progress_bar(nr_particles_todo);
		}
	}
	long int nr_ori_particles_done = my_subset_first_ori_particle;

	long int nr_subset_particles_done = 0;
	while (nr_subset_particles_done < subset_size)
	{

#ifdef TIMING
		timer.tic(TIMING_EXP_METADATA);
#endif

		long int my_first_ori_particle = my_subset_first_ori_particle + nr_subset_particles_done;
		long int my_last_ori_particle = XMIPP_MIN(my_subset_last_ori_particle, my_first_ori_particle + nr_pool - 1);

		// Get the metadata for these particles
		getMetaAndImageDataSubset(my_first_ori_particle, my_last_ori_particle, !do_parallel_disc_io);

#ifdef TIMING
		timer.toc(TIMING_EXP_METADATA);
#endif

		// perform the actual expectation step on several particles
		expectationSomeParticles(my_first_ori_particle, my_last_ori_particle);

#ifdef TIMING
		timer.tic(TIMING_EXP_CHANGES);
#endif

		// Also monitor the changes in the optimal orientations and classes
		monitorHiddenVariableChanges(my_first_ori_particle, my_last_ori_particle);

#ifdef TIMING
		timer.toc(TIMING_EXP_CHANGES);
		timer.tic(TIMING_EXP_METADATA);
#endif

		// Set the metadata for these particles
		setMetaDataSubset(my_first_ori_particle, my_last_ori_particle);

#ifdef TIMING
		timer.toc(TIMING_EXP_METADATA);
#endif

		nr_subset_particles_done += my_last_ori_particle - my_first_ori_particle + 1;
		nr_ori_particles_done += my_last_ori_particle - my_first_ori_particle + 1;

		if (verb > 0 && nr_ori_particles_done - prev_barstep > barstep)
		{
			prev_barstep = nr_ori_particles_done;
			progress_bar(nr_ori_particles_done);
		}
	}

	if (subset_size < 0 && verb > 0)
		progress_bar(nr_particles_todo);

#ifdef CUDA
	if (do_gpu)
	{
		for (int i = 0; i < cudaDeviceBundles.size(); i ++)
		{
			MlDeviceBundle* b = ((MlDeviceBundle*)cudaDeviceBundles[i]);
			b->syncAllBackprojects();

			for (int j = 0; j < b->cudaProjectors.size(); j++)
			{
				unsigned long s = wsum_model.BPref[j].data.nzyxdim;
				XFLOAT *reals = new XFLOAT[s];
				XFLOAT *imags = new XFLOAT[s];
				XFLOAT *weights = new XFLOAT[s];

				b->cudaBackprojectors[j].getMdlData(reals, imags, weights);

				for (unsigned long n = 0; n < s; n++)
				{
					wsum_model.BPref[j].data.data[n].real += (RFLOAT) reals[n];
					wsum_model.BPref[j].data.data[n].imag += (RFLOAT) imags[n];
					wsum_model.BPref[j].weight.data[n] += (RFLOAT) weights[n];
				}

				delete [] reals;
				delete [] imags;
				delete [] weights;

				b->cudaProjectors[j].clear();
				b->cudaBackprojectors[j].clear();
				b->coarseProjectionPlans[j].clear();
			}
		}

		for (int i = 0; i < cudaOptimisers.size(); i ++)
			delete (MlOptimiserCuda*) cudaOptimisers[i];

		cudaOptimisers.clear();


		for (int i = 0; i < cudaDeviceBundles.size(); i ++)
		{

			((MlDeviceBundle*)cudaDeviceBundles[i])->allocator->syncReadyEvents();
			((MlDeviceBundle*)cudaDeviceBundles[i])->allocator->freeReadyAllocs();

#ifdef DEBUG_CUDA
			if (((MlDeviceBundle*) cudaDeviceBundles[i])->allocator->getNumberOfAllocs() != 0)
			{
				printf("DEBUG_ERROR: Non-zero allocation count encountered in custom allocator between iterations.\n");
				((MlDeviceBundle*) cudaDeviceBundles[i])->allocator->printState();
				fflush(stdout);
				CRITICAL(ERR_CANZ);
			}

#endif
		}

		for (int i = 0; i < cudaDeviceBundles.size(); i ++)
			delete (MlDeviceBundle*) cudaDeviceBundles[i];

		cudaDeviceBundles.clear();
	}
#endif

	// Set back verb
	verb = old_verb;

	// Clean up some memory
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		mymodel.PPref[iclass].data.clear();
#ifdef DEBUG_EXP
	std::cerr << "Expectation: done " << std::endl;
#endif

}


void MlOptimiser::expectationSetup()
{
#ifdef DEBUG
	std::cerr << "Entering expectationSetup" << std::endl;
#endif

	// Re-initialise the random seed, because with a noisy_mask, inside the previous iteration different timings of different MPI nodes may have given rise to different number of calls to ran1
	// Use the iteration number so that each iteration has a different random seed
	init_random_generator(random_seed + iter);

	// Reset the random perturbation for this sampling
	sampling.resetRandomlyPerturbedSampling();

    // Initialise Projectors and fill vector with power_spectra for all classes
	mymodel.setFourierTransformMaps(!fix_tau, nr_threads, do_gpu);

	// TMP for helices of Anthiony 12 july 2016
	if (debug1 > 0.)
	{
		for (int iclass = 0; iclass < mymodel.nr_classes*mymodel.nr_bodies; iclass++)
			mymodel.PPref[iclass].applyFourierMask((int)debug1, (int)debug2, debug3);
	}

	// Initialise all weighted sums to zero
	wsum_model.initZeros();

	// If we're doing SGD with gradual decrease of sigma2_fudge: calculate current fudge-factor here
	if (nr_subsets > 1 && sgd_sigma2fudge_halflife > 0)
	{
		// Subtract 1 from subset, as this is done BEFORE the current subset
		long int total_nr_subsets = ((iter - 1) * nr_subsets) + (subset - 1);
		long int NN = (total_nr_subsets * subset_size);
		RFLOAT f = (RFLOAT)NN / (RFLOAT)(NN + sgd_sigma2fudge_halflife);
		// new sigma2_fudge = f * 1.0 * (1 - f) * sgd_ini_sigma2fudge
		sigma2_fudge = f + (1. - f) * sgd_sigma2fudge_ini;
	}

}

void MlOptimiser::expectationSetupCheckMemory(int myverb)
{

#ifdef DEBUG_BODIES
	myverb = false;
#endif
	std::vector<int> pointer_dir_nonzeroprior, pointer_psi_nonzeroprior;
	std::vector<RFLOAT> directions_prior, psi_prior;
	if (mymodel.orientational_prior_mode != NOPRIOR)
	{
		// First select one random direction and psi-angle for selectOrientationsWithNonZeroPriorProbability
		// This is to get an idea how many non-zero probabilities there will be
		RFLOAT ran_rot, ran_tilt, ran_psi;
		int randir = (int)(rnd_unif() * sampling.NrDirections() );
		int ranpsi = (int)(rnd_unif() * sampling.NrPsiSamplings() );
		sampling.getDirection(randir, ran_rot, ran_tilt);
		sampling.getPsiAngle(ranpsi, ran_psi);
		// Calculate local searches for these angles
		// Jun04,2015 - Shaoda & Sjors, bimodal psi searches for helices
		if (do_helical_refine)
		{
			bool do_auto_refine_local_searches = (do_auto_refine) && (sampling.healpix_order >= autosampling_hporder_local_searches);
			bool do_classification_local_searches = (!do_auto_refine) && (mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI)
					&& (mymodel.sigma2_rot > 0.) && (mymodel.sigma2_tilt > 0.) && (mymodel.sigma2_psi > 0.);
			bool do_local_angular_searches = (do_auto_refine_local_searches) || (do_classification_local_searches);
			sampling.selectOrientationsWithNonZeroPriorProbabilityFor3DHelicalReconstruction(ran_rot, ran_tilt, ran_psi,
									sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
									pointer_dir_nonzeroprior, directions_prior, pointer_psi_nonzeroprior, psi_prior,
									do_local_angular_searches);
		}
		else
		{
			sampling.selectOrientationsWithNonZeroPriorProbability(ran_rot, ran_tilt, ran_psi,
									sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
									pointer_dir_nonzeroprior, directions_prior, pointer_psi_nonzeroprior, psi_prior,
									((do_bimodal_psi) && (mymodel.sigma2_psi > 0.)) );
		}
	}


	if (myverb > 0)
	{
		// Calculate number of sampled hidden variables:
		int nr_ang_steps = CEIL(PI * particle_diameter * mymodel.current_resolution);
		RFLOAT myresol_angstep = 360. / nr_ang_steps;
		std::cout << " CurrentResolution= " << 1./mymodel.current_resolution << " Angstroms, which requires orientationSampling of at least "<< myresol_angstep
				   <<" degrees for a particle of diameter "<< particle_diameter << " Angstroms"<< std::endl;
		for (int oversampling = 0; oversampling <= adaptive_oversampling; oversampling++)
		{
			std::cout << " Oversampling= " << oversampling << " NrHiddenVariableSamplingPoints= " << mymodel.nr_classes * sampling.NrSamplingPoints(oversampling, &pointer_dir_nonzeroprior, &pointer_psi_nonzeroprior) << std::endl;
			int nr_orient = (do_only_sample_tilt) ? sampling.NrDirections(oversampling, &pointer_dir_nonzeroprior) : sampling.NrDirections(oversampling, &pointer_dir_nonzeroprior) * sampling.NrPsiSamplings(oversampling, &pointer_psi_nonzeroprior);
			std::cout << " OrientationalSampling= " << sampling.getAngularSampling(oversampling) << " NrOrientations= "<< nr_orient <<std::endl;
			if ( (do_helical_refine) && (!ignore_helical_symmetry) )
				std::cout << " TranslationalSamplingAlongHelicalAxis= " << sampling.getHelicalTranslationalSampling(oversampling) << std::flush;
			std::cout << " TranslationalSampling= " << sampling.getTranslationalSampling(oversampling)
					<< " NrTranslations= "<< sampling.NrTranslationalSamplings(oversampling) << std::endl;
			std::cout << "=============================" << std::endl;
		}
	}


	if (myverb > 1)
	{
		// Check whether things will fit into memory
		// Each RFLOAT takes 8 bytes, and their are mymodel.nr_classes references, express in Gb
		RFLOAT Gb = sizeof(RFLOAT) / (1024. * 1024. * 1024.);
		// A. Calculate approximate size of the reference maps
		// Forward projector has complex data, backprojector has complex data and real weight
		RFLOAT mem_references = Gb * mymodel.nr_classes * (2 * MULTIDIM_SIZE((mymodel.PPref[0]).data) + 3 * MULTIDIM_SIZE((wsum_model.BPref[0]).data));
		// B. Weight vectors
		RFLOAT mem_pool = Gb * mymodel.nr_classes * sampling.NrSamplingPoints(adaptive_oversampling,
				&pointer_dir_nonzeroprior, &pointer_psi_nonzeroprior);
		// C. The original image data
		int nr_pix = (mymodel.data_dim == 2) ? mymodel.current_size * mymodel.current_size : mymodel.current_size * mymodel.current_size * mymodel.current_size;
		mem_pool += Gb * nr_pix;
		if (!do_shifts_onthefly)
		{
			// D. All precalculated shifted images as well (both masked and unmasked)
			mem_pool += Gb * nr_pix * 2 * sampling.NrTranslationalSamplings(adaptive_oversampling);
		}
		// Estimate the rest of the program at 0.1 Gb?
		RFLOAT mem_rest = 0.1; // This one does NOT scale with nr_pool
		// Use tabulated sine and cosine values instead for 2D helical segments / 3D helical sub-tomogram averaging with on-the-fly shifts
		if ( (do_shifts_onthefly) && (!((do_helical_refine) && (!ignore_helical_symmetry))) )
		{
			// E. Store all AB-matrices
			mem_rest += Gb * nr_pix * sampling.NrTranslationalSamplings(adaptive_oversampling);
		}

		RFLOAT total_mem_Gb_exp = mem_references + nr_pool * mem_pool + mem_rest;
		// Each reconstruction has to store 1 extra complex array (Fconv) and 4 extra RFLOAT arrays (Fweight, Fnewweight. vol_out and Mconv in convoluteBlobRealSpace),
		// in adddition to the RFLOAT weight-array and the complex data-array of the BPref
		// That makes a total of 2*2 + 5 = 9 * a RFLOAT array of size BPref
		RFLOAT total_mem_Gb_max = Gb * 9 * MULTIDIM_SIZE((wsum_model.BPref[0]).data);

		std::cout << " Estimated memory for expectation  step > " << total_mem_Gb_exp << " Gb."<<std::endl;
		std::cout << " Estimated memory for maximization step > " << total_mem_Gb_max << " Gb."<<std::endl;
	}

#ifdef DEBUG
	std::cerr << "Leaving expectationSetup" << std::endl;
#endif

}

void MlOptimiser::precalculateABMatrices()
{

	// Set the global AB-matrices for the FFT phase-shifted images
	global_fftshifts_ab_coarse.clear();
	global_fftshifts_ab_current.clear();
	global_fftshifts_ab2_coarse.clear();
	global_fftshifts_ab2_current.clear();
	MultidimArray<Complex> Fab_current, Fab_coarse;
	if (mymodel.data_dim == 3)
		Fab_current.resize(mymodel.current_size, mymodel.current_size, mymodel.current_size / 2 + 1);
	else
		Fab_current.resize(mymodel.current_size, mymodel.current_size / 2 + 1);
	long int exp_nr_trans = sampling.NrTranslationalSamplings();
	std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
	// Note that do_shifts_onthefly is incompatible with do_skip_align because of the loop below
	for (long int itrans = 0; itrans < exp_nr_trans; itrans++)
	{
		// First get the non-oversampled translations as defined by the sampling object
		// Feb01,2017 - Shaoda, obsolete, helical reconstuctions never call this function
		sampling.getTranslations(itrans, 0, oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
				(do_helical_refine) && (!ignore_helical_symmetry), helical_rise_initial / mymodel.pixel_size, helical_twist_initial); // need getTranslations to add random_perturbation

		// Precalculate AB-matrices
		RFLOAT tmp_zoff = (mymodel.data_dim == 2) ? (0.) : oversampled_translations_z[0];
		getAbMatricesForShiftImageInFourierTransform(Fab_current, Fab_current, (RFLOAT)mymodel.ori_size, oversampled_translations_x[0], oversampled_translations_y[0], tmp_zoff);

		windowFourierTransform(Fab_current, Fab_coarse, coarse_size);
		global_fftshifts_ab_coarse.push_back(Fab_coarse);
		if (adaptive_oversampling == 0)
		{
			global_fftshifts_ab_current.push_back(Fab_current);
		}
		else
		{
			// Then also loop over all its oversampled relatives
			// Then loop over all its oversampled relatives
			// Feb01,2017 - Shaoda, obsolete, helical reconstuctions never call this function
			sampling.getTranslations(itrans, adaptive_oversampling, oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
					(do_helical_refine) && (!ignore_helical_symmetry), helical_rise_initial / mymodel.pixel_size, helical_twist_initial);
			for (long int iover_trans = 0; iover_trans < oversampled_translations_x.size(); iover_trans++)
			{
				// Shift through phase-shifts in the Fourier transform
				// Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)

				RFLOAT tmp_zoff = (mymodel.data_dim == 2) ? (0.) : oversampled_translations_z[iover_trans];
				getAbMatricesForShiftImageInFourierTransform(Fab_current, Fab_current, (RFLOAT)mymodel.ori_size, oversampled_translations_x[iover_trans], oversampled_translations_y[iover_trans], tmp_zoff);

				global_fftshifts_ab2_current.push_back(Fab_current);
				if (strict_highres_exp > 0.)
				{
					windowFourierTransform(Fab_current, Fab_coarse, coarse_size);
					global_fftshifts_ab2_coarse.push_back(Fab_coarse);
				}
			}
		} // end else (adaptive_oversampling == 0)
	} // end loop itrans

#ifdef DEBUG_AB
	std::cerr << " global_fftshifts_ab_coarse.size()= " << global_fftshifts_ab_coarse.size() << " global_fftshifts_ab_current.size()= " << global_fftshifts_ab_current.size() << std::endl;
	std::cerr << " global_fftshifts_ab2_coarse.size()= " << global_fftshifts_ab2_coarse.size() << " global_fftshifts_ab2_current.size()= " << global_fftshifts_ab2_current.size() << std::endl;
#endif

}

void MlOptimiser::expectationSomeParticles(long int my_first_ori_particle, long int my_last_ori_particle)
{

#ifdef TIMING
	timer.tic(TIMING_ESP);
#endif

//#define DEBUG_EXPSOME
#ifdef DEBUG_EXPSOME
	std::cerr << "Entering expectationSomeParticles..." << std::endl;
#endif

    // Use global variables for thread visibility (before there were local ones for similar call in MPI version!)
	exp_my_first_ori_particle = my_first_ori_particle;
    exp_my_last_ori_particle = my_last_ori_particle;

	// Make sure random division is always the same with the same seed
    if (do_generate_seeds && ((do_firstiter_cc && iter == 2) || (!do_firstiter_cc && iter == 1)) )
	{
    	// calculate the random class for these SomeParticles
    	exp_random_class_some_particles.clear();
    	for (long int ori_part_id = my_first_ori_particle; ori_part_id <= my_last_ori_particle; ori_part_id++)
    	{
        	init_random_generator(random_seed + ori_part_id);
    		int random_class = rand() % mymodel.nr_classes;
    		exp_random_class_some_particles.push_back(random_class);
    	}
	}

    // Only open/close stacks once
    fImageHandler hFile;
	long int dump;
	FileName fn_img, fn_stack, fn_open_stack="";

	// Store total number of particle images in this bunch of SomeParticles, and set translations and orientations for skip_align/rotate
    exp_nr_images = 0;
    long int istop = 0;
    exp_imgs.clear();
    for (long int ori_part_id = my_first_ori_particle; ori_part_id <= my_last_ori_particle; ori_part_id++)
	{

		// If skipping alignment (not for movies) or rotations (for movies, but all frames (part_id) have the same rotation)
		// then store the old translation and orientation for each ori_particle (take first part_id!)
		// If we do local angular searches, get the previously assigned angles to center the prior
    	bool do_clear = (ori_part_id == my_first_ori_particle);
		if (do_skip_align || do_skip_rotate)
		{
			// Also set the rotations
			RFLOAT old_rot, old_tilt, old_psi;
			old_rot = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_ROT);
			old_tilt = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_TILT);
			old_psi = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_PSI);
			sampling.addOneOrientation(old_rot, old_tilt, old_psi, do_clear);
		}
		else if (do_only_sample_tilt)
		{
			if (do_clear) // only clear psi_angles for the first particle, as one psi-angle is stored for each particle!
				sampling.psi_angles.clear();
			RFLOAT old_psi = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_PSI);
			sampling.psi_angles.push_back(old_psi);
		}
		if (do_skip_align)
		{
			// Rounded translations will be applied to the image upon reading,
			// set the unique translation in the sampling object to the fractional difference
			RFLOAT my_old_offset_x, my_old_offset_y, my_old_offset_z;
			RFLOAT rounded_offset_x, rounded_offset_y, rounded_offset_z;
			RFLOAT rot_deg, tilt_deg, psi_deg;
			my_old_offset_x = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_XOFF);
			my_old_offset_y = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_YOFF);
			rounded_offset_x = my_old_offset_x - ROUND(my_old_offset_x);
			rounded_offset_y = my_old_offset_y - ROUND(my_old_offset_y);
			if (mymodel.data_dim == 3)
			{
				my_old_offset_z = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_ZOFF);
				rounded_offset_z = my_old_offset_z - ROUND(my_old_offset_z);
			}
			if (do_helical_refine)
			{
				rot_deg = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_ROT);
				tilt_deg = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_TILT);
				psi_deg = DIRECT_A2D_ELEM(exp_metadata, exp_nr_images, METADATA_PSI);
			}
			sampling.addOneTranslation(rounded_offset_x, rounded_offset_y, rounded_offset_z, do_clear, (do_helical_refine) && (!ignore_helical_symmetry), rot_deg, tilt_deg, psi_deg); // clear for first ori_particle
		}

		// Store total number of particle images in this bunch of SomeParticles
		exp_nr_images += mydata.ori_particles[ori_part_id].particles_id.size();

		// Sjors 7 March 2016 to prevent too high disk access... Read in all pooled images simultaneously
		// TODO: open and close stacks only once!
		// Don't do this for sub-tomograms to save RAM!
		if (do_parallel_disc_io && !do_preread_images && mymodel.data_dim != 3)
		{
			// Read in all images, only open/close common stacks once
			for (int ipart = 0; ipart < mydata.ori_particles[ori_part_id].particles_id.size(); ipart++, istop++)
			{

				long int part_id = mydata.ori_particles[ori_part_id].particles_id[ipart];

				// Read from disc
				// Get the filename
				if (!mydata.getImageNameOnScratch(part_id, fn_img))
				{
					std::istringstream split(exp_fn_img);
					for (int i = 0; i <= istop; i++)
						getline(split, fn_img);
				}

				// Only open again a new stackname
				fn_img.decompose(dump, fn_stack);
				if (fn_stack != fn_open_stack)
				{
					hFile.openFile(fn_stack, WRITE_READONLY);
					fn_open_stack = fn_stack;
				}
			    Image<RFLOAT> img;
				img.readFromOpenFile(fn_img, hFile, -1, false);
				img().setXmippOrigin();
				exp_imgs.push_back(img());

			} // end loop over all particles of this ori_particle

		} // end if do_parallel_disc_io

	} //end loop over ori_part_id


#ifdef DEBUG_EXPSOME
	std::cerr << " exp_my_first_ori_particle= " << exp_my_first_ori_particle << " exp_my_last_ori_particle= " << exp_my_last_ori_particle << std::endl;
	std::cerr << " exp_nr_images= " << exp_nr_images << std::endl;
#endif

	exp_ipart_ThreadTaskDistributor->resize(my_last_ori_particle - my_first_ori_particle + 1, 1);
	exp_ipart_ThreadTaskDistributor->reset();
    global_ThreadManager->run(globalThreadExpectationSomeParticles);

    if (threadException != NULL)
    	throw *threadException;

#ifdef TIMING
    timer.toc(TIMING_ESP);
#endif


}


void MlOptimiser::doThreadExpectationSomeParticles(int thread_id)
{

#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		timer.tic(TIMING_ESP_THR);
#endif

	size_t first_ipart = 0, last_ipart = 0;
	while (exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
//#define DEBUG_EXPSOMETHR
#ifdef DEBUG_EXPSOMETHR
		pthread_mutex_lock(&global_mutex);
		std::cerr << " thread_id= " << thread_id << " first_ipart= " << first_ipart << " last_ipart= " << last_ipart << std::endl;
		std::cerr << " exp_my_first_ori_particle= " << exp_my_first_ori_particle << " exp_my_last_ori_particle= " << exp_my_last_ori_particle << std::endl;
		pthread_mutex_unlock(&global_mutex);
#endif

		for (long int ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
#ifdef TIMING
			// Only time one thread
			if (thread_id == 0)
				timer.tic(TIMING_ESP_ONEPART);
			else if (thread_id == nr_threads -1)
				timer.tic(TIMING_ESP_ONEPARTN);
#endif
			expectationOneParticle(exp_my_first_ori_particle + ipart, thread_id);

#ifdef TIMING
			// Only time one thread
			if (thread_id == 0)
				timer.toc(TIMING_ESP_ONEPART);
			else if (thread_id == nr_threads -1)
				timer.toc(TIMING_ESP_ONEPARTN);
#endif

		}
	}

#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		timer.toc(TIMING_ESP_THR);
#endif

}


void MlOptimiser::expectationOneParticle(long int my_ori_particle, int thread_id)
{

#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
		timer.tic(TIMING_ESP_INI);
#endif

	// In the first iteration, multiple seeds will be generated
	// A single random class is selected for each pool of images, and one does not marginalise over the orientations
	// The optimal orientation is based on signal-product (rather than the signal-intensity sensitive Gaussian)
    // If do_firstiter_cc, then first perform a single iteration with K=1 and cross-correlation criteria, afterwards

    // Decide which classes to integrate over (for random class assignment in 1st iteration)
    int exp_iclass_min = 0;
    int exp_iclass_max = mymodel.nr_classes - 1;
    // low-pass filter again and generate the seeds
    if (do_generate_seeds)
    {
    	if (do_firstiter_cc && iter == 1)
    	{
    		// In first (CC) iter, use a single reference (and CC)
    		exp_iclass_min = exp_iclass_max = 0;
    	}
    	else if ( (do_firstiter_cc && iter == 2) || (!do_firstiter_cc && iter == 1))
		{
			// In second CC iter, or first iter without CC: generate the seeds
    		// Now select a single random class
    		// exp_part_id is already in randomized order (controlled by -seed)
    		// WARNING: USING SAME iclass_min AND iclass_max FOR SomeParticles!!
    		long int idx = my_ori_particle - exp_my_first_ori_particle;
    		if (idx >= exp_random_class_some_particles.size())
    			REPORT_ERROR("BUG: expectationOneParticle idx>random_class_some_particles.size()");
    		exp_iclass_min = exp_iclass_max = exp_random_class_some_particles[idx];
		}
    }


// This debug is a good one to step through the separate steps of the expectation to see where trouble lies....
//#define DEBUG_ESP_MEM
#ifdef DEBUG_ESP_MEM

	std::cerr << "Entering MlOptimiser::expectationOneParticle" << std::endl;
    std::cerr << " my_ori_particle= " << my_ori_particle << std::endl;
    if (thread_id==0)
	{
		char c;
		std::cerr << "Before getFourierTransformsAndCtfs, press any key to continue... " << std::endl;
		std::cin >> c;
	}
    global_barrier->wait();
#endif


    // Loop over all bodies of the multi-body refinement
    // Basically, subsequently align and store weighted sums for each body
    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++)
    {

		// Here define all kind of local arrays that will be needed
		std::vector<MultidimArray<Complex > > exp_Fimgs, exp_Fimgs_nomask, exp_local_Fimgs_shifted, exp_local_Fimgs_shifted_nomask;
		std::vector<MultidimArray<RFLOAT> > exp_Fctfs, exp_local_Fctfs, exp_local_Minvsigma2s;
		std::vector<int> exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior;
		std::vector<RFLOAT> exp_directions_prior, exp_psi_prior, exp_local_sqrtXi2;
		int exp_current_image_size, exp_current_oversampling;
		std::vector<RFLOAT> exp_highres_Xi2_imgs, exp_min_diff2;
		MultidimArray<RFLOAT> exp_Mweight;
		MultidimArray<bool> exp_Mcoarse_significant;
		// And from storeWeightedSums
		std::vector<RFLOAT> exp_sum_weight, exp_significant_weight, exp_max_weight;
		std::vector<Matrix1D<RFLOAT> > exp_old_offset, exp_prior;
		std::vector<RFLOAT> exp_wsum_norm_correction;
		std::vector<MultidimArray<RFLOAT> > exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA, exp_power_imgs;
		RFLOAT exp_thisparticle_sumweight;

		int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
		// Global exp_metadata array has metadata of all ori_particles. Where does my_ori_particle start?
		int metadata_offset = 0;
		for (long int iori = exp_my_first_ori_particle; iori <= exp_my_last_ori_particle; iori++)
		{
			if (iori == my_ori_particle)
				break;
			metadata_offset += mydata.ori_particles[iori].particles_id.size();
		}

		// Resize vectors for all particles
		exp_power_imgs.resize(exp_nr_particles);
		exp_highres_Xi2_imgs.resize(exp_nr_particles);
		exp_Fimgs.resize(exp_nr_particles);
		exp_Fimgs_nomask.resize(exp_nr_particles);
		exp_Fctfs.resize(exp_nr_particles);
		exp_old_offset.resize(exp_nr_particles);
		exp_prior.resize(exp_nr_particles);

		// Then calculate all Fourier Transforms
#ifdef TIMING
		if (my_ori_particle == exp_my_first_ori_particle)
		{
			timer.toc(TIMING_ESP_INI);
			timer.tic(TIMING_ESP_FT);
		}
#endif

		getFourierTransformsAndCtfs(my_ori_particle, ibody, metadata_offset, exp_Fimgs, exp_Fimgs_nomask, exp_Fctfs,
				exp_old_offset, exp_prior, exp_power_imgs, exp_highres_Xi2_imgs,
				exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior, exp_directions_prior, exp_psi_prior);

#ifdef TIMING
		if (my_ori_particle == exp_my_first_ori_particle)
			timer.toc(TIMING_ESP_FT);
#endif

#ifdef DEBUG_ESP_MEM
		if (thread_id==0)
		{
			char c;
			std::cerr << " exp_nr_particles= " << exp_nr_particles << " metadata_offset= " << metadata_offset << std::endl;
			std::cerr << "After getFourierTransformsAndCtfs, press any key to continue... " << std::endl;
			std::cin >> c;
		}
		global_barrier->wait();
#endif

		if (do_realign_movies && movie_frame_running_avg_side > 0)
		{
			calculateRunningAveragesOfMovieFrames(my_ori_particle, exp_Fimgs, exp_power_imgs, exp_highres_Xi2_imgs);
		}

		// To deal with skipped alignments/rotations
		int exp_itrans_min, exp_itrans_max, exp_idir_min, exp_idir_max, exp_ipsi_min, exp_ipsi_max;
		if (do_skip_align)
		{
			exp_itrans_min = exp_itrans_max = my_ori_particle - exp_my_first_ori_particle;
		}
		else
		{
			exp_itrans_min = 0;
			exp_itrans_max = sampling.NrTranslationalSamplings() - 1;
		}
		if (do_skip_align || do_skip_rotate)
		{
			exp_idir_min = exp_idir_max = exp_ipsi_min = exp_ipsi_max =
					my_ori_particle - exp_my_first_ori_particle;
		}
		else if (do_only_sample_tilt)
		{
			exp_idir_min = 0;
			exp_idir_max = sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior) - 1;
			exp_ipsi_min = exp_ipsi_max = my_ori_particle - exp_my_first_ori_particle;
		}
		else
		{
			exp_idir_min = exp_ipsi_min = 0;
			exp_idir_max = sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior) - 1;
			exp_ipsi_max = sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior ) - 1;
		}

		// Initialise significant weight to minus one, so that all coarse sampling points will be handled in the first pass
		exp_significant_weight.resize(exp_nr_particles, -1.);

		// Only perform a second pass when using adaptive oversampling
		int nr_sampling_passes = (adaptive_oversampling > 0) ? 2 : 1;

		// Pass twice through the sampling of the entire space of rot, tilt and psi
		// The first pass uses a coarser angular sampling and possibly smaller FFTs than the second pass.
		// Only those sampling points that contribute to the highest x% of the weights in the first pass are oversampled in the second pass
		// Only those sampling points will contribute to the weighted sums in the third loop below
		for (int exp_ipass = 0; exp_ipass < nr_sampling_passes; exp_ipass++)
		{

			if (strict_highres_exp > 0.)
				// Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
				exp_current_image_size = coarse_size;
			else if (adaptive_oversampling > 0)
				// Use smaller images in the first pass, larger ones in the second pass
				exp_current_image_size = (exp_ipass == 0) ? coarse_size : mymodel.current_size;
			else
				exp_current_image_size = mymodel.current_size;

			// Use coarse sampling in the first pass, oversampled one the second pass
			exp_current_oversampling = (exp_ipass == 0) ? 0 : adaptive_oversampling;

#ifdef DEBUG_ESP_MEM
			if (thread_id==0)
			{
				char c;
				std::cerr << " exp_current_image_size= " << exp_current_image_size << " exp_current_oversampling= " << exp_current_oversampling << " nr_sampling_passes= " << nr_sampling_passes << std::endl;
				std::cerr << "Before getAllSquaredDifferences, use top to see memory usage and then press any key to continue... " << std::endl;
				std::cin >> c;
			}
			global_barrier->wait();
#endif

			// Calculate the squared difference terms inside the Gaussian kernel for all hidden variables
			getAllSquaredDifferences(my_ori_particle, ibody, exp_current_image_size, exp_ipass, exp_current_oversampling,
					metadata_offset, exp_idir_min, exp_idir_max, exp_ipsi_min, exp_ipsi_max,
					exp_itrans_min, exp_itrans_max, exp_iclass_min, exp_iclass_max, exp_min_diff2, exp_highres_Xi2_imgs,
					exp_Fimgs, exp_Fctfs, exp_Mweight, exp_Mcoarse_significant,
					exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior, exp_directions_prior, exp_psi_prior,
					exp_local_Fimgs_shifted, exp_local_Minvsigma2s, exp_local_Fctfs, exp_local_sqrtXi2);


#ifdef DEBUG_ESP_MEM
			if (thread_id==0)
			{
				char c;
				std::cerr << "After getAllSquaredDifferences, use top to see memory usage and then press any key to continue... " << std::endl;
				std::cin >> c;
			}
			global_barrier->wait();
#endif

			// Now convert the squared difference terms to weights,
			// also calculate exp_sum_weight, and in case of adaptive oversampling also exp_significant_weight
			convertAllSquaredDifferencesToWeights(my_ori_particle, exp_ipass, exp_current_oversampling, metadata_offset,
					exp_idir_min, exp_idir_max, exp_ipsi_min, exp_ipsi_max,
					exp_itrans_min, exp_itrans_max, exp_iclass_min, exp_iclass_max,
					exp_Mweight, exp_Mcoarse_significant, exp_significant_weight,
					exp_sum_weight, exp_old_offset, exp_prior, exp_min_diff2,
					exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior, exp_directions_prior, exp_psi_prior);

#ifdef DEBUG_ESP_MEM
		if (thread_id==0)
		{
			char c;
			std::cerr << "After convertAllSquaredDifferencesToWeights, press any key to continue... " << std::endl;
			std::cin >> c;
		}
		global_barrier->wait();
#endif

		}// end loop over 2 exp_ipass iterations

#ifdef RELION_TESTING
		std::string mode;
		if (do_gpu)
		{
			mode="gpu";
		}
		else
		{
			mode="cpu";
		}
		std::cerr << " "<< std::endl;
		std::cerr << " finished running diffs in  " << mode << " mode."<< std::endl;
		Image<RFLOAT> tt;
		tt().resize(exp_current_image_size, exp_current_image_size);
		MultidimArray<Complex> Fimg1;
		Fimg1 = exp_local_Fimgs_shifted[0];
		FourierTransformer transformer;
		transformer.inverseFourierTransform(Fimg1, tt());
		CenterFFT(tt(),false);
		std::string fnm = mode + std::string("_out_shifted_image.mrc");
		tt.write(fnm);
		tt().resize(YSIZE(Mresol_coarse),XSIZE(Mresol_coarse));
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(tt())
		{
			DIRECT_MULTIDIM_ELEM(tt(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(Mresol_coarse, n);
		}
		fnm = mode + std::string("_out_mresol_coarse.mrc");
		tt.write(fnm);
		tt().resize(YSIZE(Mresol_fine),XSIZE(Mresol_fine));
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(tt())
		{
			DIRECT_MULTIDIM_ELEM(tt(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
		}
		fnm = mode + std::string("_out_mresol_fine.mrc");
		tt.write(fnm);
		tt().resize(YSIZE(exp_local_Fctfs[0]),XSIZE(exp_local_Fctfs[0]));
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fctfs[0])
		{
			DIRECT_MULTIDIM_ELEM(tt(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[0], n);
		}
		fnm = mode + std::string("_out_ctf.mrc");
		tt.write(fnm);
//		exp_local_Fctfs[0]), tt());
//		CenterFFT(tt(),false);
//		fnm = mode + std::string("_out_ctf.mrc");
//		tt.write(fnm);
		fnm = mode + std::string("_out_10k_weights.txt");
		char *text = &fnm[0];
		freopen(text,"w",stdout);
		for(int n=0; n<10000; n++)
		{
			printf("%4.8f \n",DIRECT_MULTIDIM_ELEM(exp_Mweight, n)); // << std::endl;
		}
		fclose(stdout);
//      	exit(0);
#endif
		// For the reconstruction step use mymodel.current_size!
		exp_current_image_size = mymodel.current_size;

#ifdef DEBUG_ESP_MEM
		if (thread_id==0)
		{
			char c;
			std::cerr << "Before storeWeightedSums, press any key to continue... " << std::endl;
			std::cin >> c;
		}
		global_barrier->wait();
#endif

		storeWeightedSums(my_ori_particle, ibody, exp_current_image_size, exp_current_oversampling, metadata_offset,
				exp_idir_min, exp_idir_max, exp_ipsi_min, exp_ipsi_max,
				exp_itrans_min, exp_itrans_max, exp_iclass_min, exp_iclass_max,
				exp_min_diff2, exp_highres_Xi2_imgs, exp_Fimgs, exp_Fimgs_nomask, exp_Fctfs,
				exp_power_imgs, exp_old_offset, exp_prior, exp_Mweight, exp_Mcoarse_significant,
				exp_significant_weight, exp_sum_weight, exp_max_weight,
				exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior, exp_directions_prior, exp_psi_prior,
				exp_local_Fimgs_shifted, exp_local_Fimgs_shifted_nomask, exp_local_Minvsigma2s, exp_local_Fctfs, exp_local_sqrtXi2);

#ifdef RELION_TESTING
//		std::string mode;
		if (do_gpu)
		{
			mode="gpu";
		}
		else
		{
			mode="cpu";
		}
		std::cerr << " "<< std::endl;
		std::cerr << " finished running diffs in  " << mode << " mode."<< std::endl;
		fnm = mode + std::string("_out_10k_weights_afterstore.txt");
		text = &fnm[0];
		freopen(text,"w",stdout);
		// Write the first 10k diffs to be sure
		for(int n=0; n<10000; n++)
		{
			//std::cout << DIRECT_MULTIDIM_ELEM(exp_Mweight, n) << std::endl;
			printf("%4.8f \n",DIRECT_MULTIDIM_ELEM(exp_Mweight, n));
		}
		//For tests we want to exit now
		//if(iter == 2)
		//	exit(0);

#endif

#ifdef DEBUG_ESP_MEM
		if (thread_id==0)
		{
			char c;
			std::cerr << "After storeWeightedSums, press any key to continue... " << std::endl;
			std::cin >> c;
		}
		global_barrier->wait();
#endif

    } // end for ibody


#ifdef DEBUG_EXPSINGLE
		std::cerr << "Leaving expectationOneParticle..." << std::endl;
#endif

}

void MlOptimiser::symmetriseReconstructions()
{
	for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++)
	{
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		{
			// either ibody or iclass can be larger than 0, never 2 at the same time!
			int ith_recons = (mymodel.nr_bodies > 1) ? ibody : iclass;

			if (mymodel.pdf_class[iclass] > 0.)
			{
				// Immediately after expectation process. Do rise and twist for all asymmetrical units in Fourier space
				// Also convert helical rise to pixels for BPref object
				wsum_model.BPref[ith_recons].symmetrise(mymodel.helical_nr_asu, mymodel.helical_twist[ith_recons], mymodel.helical_rise[ith_recons] / mymodel.pixel_size);
			}
		}
	}
	return;
}

void MlOptimiser::applyLocalSymmetryForEachRef()
{
	if ( (fn_local_symmetry_masks.size() < 1) || (fn_local_symmetry_operators.size() < 1) )
		return;

	if (verb > 0)
		std::cout << " Applying local symmetry in real space according to " << fn_local_symmetry_operators.size() << " operators..." << std::endl;

	for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++)
	{
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		{
			// either ibody or iclass can be larger than 0, never 2 at the same time!
			int ith_recons = (mymodel.nr_bodies > 1) ? ibody : iclass;
			applyLocalSymmetry(mymodel.Iref[ith_recons], fn_local_symmetry_masks, fn_local_symmetry_operators);
		}
	}
}

void MlOptimiser::makeGoodHelixForEachRef()
{
	if ( (!do_helical_refine) || (ignore_helical_symmetry) )
		return;

	for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++)
	{
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		{
			// either ibody or iclass can be larger than 0, never 2 at the same time!
			int ith_recons = (mymodel.nr_bodies > 1) ? ibody : iclass;
			if (mymodel.pdf_class[iclass] > 0.)
			{
				if ( (iter > 1) && (do_helical_symmetry_local_refinement) )
				{
					localSearchHelicalSymmetry(
							mymodel.Iref[ith_recons],
							mymodel.pixel_size,
							(particle_diameter / 2.),
							(helical_tube_inner_diameter / 2.),
							(helical_tube_outer_diameter / 2.),
							helical_z_percentage,
							mymodel.helical_rise_min,
							mymodel.helical_rise_max,
							mymodel.helical_rise_inistep,
							mymodel.helical_rise[iclass],
							mymodel.helical_twist_min,
							mymodel.helical_twist_max,
							mymodel.helical_twist_inistep,
							mymodel.helical_twist[iclass]);
				}
				imposeHelicalSymmetryInRealSpace(
						mymodel.Iref[ith_recons],
						mymodel.pixel_size,
						(particle_diameter / 2.),
						(helical_tube_inner_diameter / 2.),
						(helical_tube_outer_diameter / 2.),
						helical_z_percentage,
						mymodel.helical_rise[iclass],
						mymodel.helical_twist[iclass],
						width_mask_edge);
			}
		}
	}

	if (verb > 0)
	{
		outputHelicalSymmetryStatus(
				iter,
				helical_rise_initial,
				mymodel.helical_rise_min,
				mymodel.helical_rise_max,
				helical_twist_initial,
				mymodel.helical_twist_min,
				mymodel.helical_twist_max,
				do_helical_symmetry_local_refinement,
				mymodel.helical_rise,
				mymodel.helical_twist,
				0., 0., 0., 0., false,
				std::cout);
	}
	return;
}

void MlOptimiser::maximization()
{

	if (verb > 0)
	{
		std::cout << " Maximization ..." << std::endl;
		init_progress_bar(mymodel.nr_classes);
	}

	// First reconstruct the images for each class
	// multi-body refinement will never get here, as it is only 3D auto-refine and that requires MPI!
	for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++)
	{
		RCTIC(timer,RCT_1);
		if (mymodel.pdf_class[iclass] > 0. || mymodel.nr_bodies > 1 )
		{

			if ((wsum_model.BPref[iclass].weight).sum() > XMIPP_EQUAL_ACCURACY)
			{
				MultidimArray<RFLOAT> Iref_old;
				long int total_nr_subsets;
				RFLOAT total_mu_fraction, number_of_effective_particles, tau2_fudge;
				if(do_sgd)
				{
					Iref_old = mymodel.Iref[iclass];
					// Still regularise here. tau2 comes from the reconstruction, sum of sigma2 is only over a single subset
					// Gradually increase tau2_fudge to account for ever increasing number of effective particles in the reconstruction
					total_nr_subsets = ((iter - 1) * nr_subsets) + subset;
					total_mu_fraction = pow (mu, (RFLOAT)total_nr_subsets);
					number_of_effective_particles = (iter == 1) ? subset * subset_size : mydata.numberOfParticles();
					number_of_effective_particles *= (1. - total_mu_fraction);
					tau2_fudge = number_of_effective_particles * mymodel.tau2_fudge_factor / subset_size;
				}
				else
				{
					tau2_fudge = mymodel.tau2_fudge_factor;
				}

				(wsum_model.BPref[iclass]).reconstruct(mymodel.Iref[iclass], gridding_nr_iter, do_map,
								tau2_fudge, mymodel.tau2_class[iclass], mymodel.sigma2_class[iclass],
								mymodel.data_vs_prior_class[iclass], mymodel.fourier_coverage_class[iclass],
								mymodel.fsc_halves_class, wsum_model.pdf_class[iclass], false, false, nr_threads, minres_map, (iclass==0));

				if(do_sgd)
				{
					// Now update formula: dV_kl^(n) = (mu) * dV_kl^(n-1) + (1-mu)*step_size*G_kl^(n)
					// where G_kl^(n) is now in mymodel.Iref[iclass]!!!
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mymodel.Igrad[iclass])
						DIRECT_MULTIDIM_ELEM(mymodel.Igrad[iclass], n) = mu * DIRECT_MULTIDIM_ELEM(mymodel.Igrad[iclass], n) + (1. - mu) * sgd_stepsize * DIRECT_MULTIDIM_ELEM(mymodel.Iref[iclass], n);

					// update formula: V_kl^(n+1) = V_kl^(n) + dV_kl^(n)
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mymodel.Iref[iclass])
					{
						DIRECT_MULTIDIM_ELEM(mymodel.Iref[iclass], n) = DIRECT_MULTIDIM_ELEM(Iref_old, n) + DIRECT_MULTIDIM_ELEM(mymodel.Igrad[iclass], n);
					}

//#define DEBUG_SGD
#ifdef DEBUG_SGD
					FileName fn_tmp="grad_class"+integerToString(iclass)+".spi";
					Image<RFLOAT> It;
					It()=mymodel.Igrad[iclass];
					It.write(fn_tmp);
					fn_tmp="ref_class"+integerToString(iclass)+".spi";
					It()=mymodel.Iref[iclass];
					It.write(fn_tmp);
#endif
					// Enforce positivity?
					// Low-pass filter according to current resolution??
					// Some sort of regularisation may be necessary....?

				}
			}
		}
		else
		{
			// When not doing SGD, initialise to zero, but when doing SGD just keep the previous reference
			if (!do_sgd)
				mymodel.Iref[iclass].initZeros();
			// When doing SGD also re-initialise the gradient to zero
			if (do_sgd)
				mymodel.Igrad[iclass].initZeros();
		}
		RCTOC(timer,RCT_1);
		if (verb > 0)
			progress_bar(iclass);
	}

	RCTIC(timer,RCT_3);
	// Then perform the update of all other model parameters
	maximizationOtherParameters();
	RCTOC(timer,RCT_3);
	RCTIC(timer,RCT_4);
	// Keep track of changes in hidden variables
	updateOverallChangesInHiddenVariables();
	RCTOC(timer,RCT_4);
	if (verb > 0)
		progress_bar(mymodel.nr_classes);

}

void MlOptimiser::maximizationOtherParameters()
{
	// Note that reconstructions are done elsewhere!
#ifdef DEBUG
	std::cerr << "Entering maximizationOtherParameters" << std::endl;
#endif

	RCTIC(timer,RCT_5);
	// Calculate total sum of weights, and average CTF for each class (for SSNR estimation)
	RFLOAT sum_weight = 0.;
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		sum_weight += wsum_model.pdf_class[iclass];

	// Update average norm_correction
	if (do_norm_correction)
	{
		mymodel.avg_norm_correction *= mu;
		mymodel.avg_norm_correction += (1. - mu) * wsum_model.avg_norm_correction / sum_weight;
	}

	if (do_scale_correction && !((iter==1 && do_firstiter_cc) || do_always_cc) )
	{
		for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
		{
			mymodel.scale_correction[igroup] *= mu;
			RFLOAT sumXA = wsum_model.wsum_signal_product_spectra[igroup].sum();
			RFLOAT sumAA = wsum_model.wsum_reference_power_spectra[igroup].sum();
			if (sumAA > 0.)
				mymodel.scale_correction[igroup] += (1. - mu) * sumXA / sumAA;
			else
				mymodel.scale_correction[igroup] += (1. - mu);
		}

		// TODO! Avoid extremities in scale estimates, because they lead to catastrophic events and instabilities in refinement
		// Let's exclude anything bigger than 5x the median or smaller than 1/5 of the median...
		// Use the median instead of the mean, because it is much more robust to outliers.
		std::vector<RFLOAT> sorted = mymodel.scale_correction;
		std::sort(sorted.begin(), sorted.end());
		RFLOAT median = sorted[mymodel.nr_groups / 2];

		RFLOAT avg_scale_correction = 0., nr_part = 0.;
		for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
		{

			if (mymodel.scale_correction[igroup] > 5. * median)
				mymodel.scale_correction[igroup] = 5. * median;
			else if (mymodel.scale_correction[igroup] < median / 5.)
				mymodel.scale_correction[igroup] =  median / 5.;

			avg_scale_correction += (RFLOAT)(mymodel.nr_particles_group[igroup]) * mymodel.scale_correction[igroup];
			nr_part += (RFLOAT)(mymodel.nr_particles_group[igroup]);
		}

		// Constrain average scale_correction to one.
		avg_scale_correction /= nr_part;
		for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
		{
			mymodel.scale_correction[igroup] /= avg_scale_correction;
//#define DEBUG_UPDATE_SCALE
#ifdef DEBUG_UPDATE_SCALE
			if (verb > 0)
			{
				std::cerr<< "Group "<<igroup+1<<": scale_correction= "<<mymodel.scale_correction[igroup]<<std::endl;
				for (int i = 0; i < XSIZE(wsum_model.wsum_reference_power_spectra[igroup]); i++)
					if (wsum_model.wsum_reference_power_spectra[igroup](i)> 0.)
						std::cerr << " i= " << i << " XA= " << wsum_model.wsum_signal_product_spectra[igroup](i)
											<< " A2= " << wsum_model.wsum_reference_power_spectra[igroup](i)
											<< " XA/A2= " << wsum_model.wsum_signal_product_spectra[igroup](i)/wsum_model.wsum_reference_power_spectra[igroup](i) << std::endl;

			}
#endif
		}

	}
	RCTOC(timer,RCT_5);
	RCTIC(timer,RCT_6);
	// Update model.pdf_class vector (for each k)
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
	{

		// Update pdf_class (for SGD: update with taking mu into account! For non-SGD: mu equals zero)
		mymodel.pdf_class[iclass] *= mu;
		mymodel.pdf_class[iclass] += (1. - mu) * wsum_model.pdf_class[iclass] / sum_weight;

		// for 2D also update priors of translations for each class!
		if (mymodel.ref_dim == 2)
		{
			if (wsum_model.pdf_class[iclass] > 0.)
			{
				mymodel.prior_offset_class[iclass] *= mu;
				mymodel.prior_offset_class[iclass] += (1. - mu) * wsum_model.prior_offset_class[iclass] / sum_weight;
			}
			else
				mymodel.prior_offset_class[iclass].initZeros();
		}

		// Use sampling.NrDirections() to include all directions (also those with zero prior probability for any given image)
		if (!(do_skip_align || do_skip_rotate))
		{
			for (int idir = 0; idir < sampling.NrDirections(); idir++)
			{
				mymodel.pdf_direction[iclass](idir) *= mu;
				mymodel.pdf_direction[iclass](idir) += (1. - mu) * wsum_model.pdf_direction[iclass](idir) / sum_weight;
			}
		}
	}

	// Update sigma2_offset
	// Factor 2 because of the 2-dimensionality of the xy-plane
	if (!fix_sigma_offset)
	{
		mymodel.sigma2_offset *= mu;
		if (mymodel.data_dim == 3)
		{
			if ( (do_helical_refine) && (!ignore_helical_symmetry) )
				mymodel.sigma2_offset += (1. - mu) * (wsum_model.sigma2_offset) / (2. * sum_weight);
			else
				mymodel.sigma2_offset += (1. - mu) * (wsum_model.sigma2_offset) / (3. * sum_weight);
		}
		else
		{
			if ( (do_helical_refine) && (!ignore_helical_symmetry) )
				mymodel.sigma2_offset += (1. - mu) * (wsum_model.sigma2_offset) / (1. * sum_weight);
			else
				mymodel.sigma2_offset += (1. - mu) * (wsum_model.sigma2_offset) / (2. * sum_weight);
		}
	}

	// TODO: update estimates for sigma2_rot, sigma2_tilt and sigma2_psi!
	RCTOC(timer,RCT_6);
	RCTIC(timer,RCT_7);
	// Also refrain from updating sigma_noise after the first iteration with first_iter_cc!
	if (!fix_sigma_noise && !((iter == 1 && do_firstiter_cc) || do_always_cc) )
	{
		for (int igroup = 0; igroup < mymodel.nr_groups; igroup++)
		{
			float tsum = 0;
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mymodel.sigma2_noise[igroup])
				tsum += wsum_model.sigma2_noise[igroup].data[n];

//			if(tsum==0) //if nothing has been done for this group, use previous intr noise2_sigma
//				wsum_model.sigma2_noise[igroup].data = mymodel.sigma2_noise[igroup].data;
			if(tsum!=0)
			{
				// Factor 2 because of the 2-dimensionality of the complex-plane
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mymodel.sigma2_noise[igroup])
				{
					DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n) *= mu;
					DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n) +=
							(1. - mu) * DIRECT_MULTIDIM_ELEM(wsum_model.sigma2_noise[igroup], n ) /
								(2. * wsum_model.sumw_group[igroup] * DIRECT_MULTIDIM_ELEM(Npix_per_shell, n));
					// Watch out for all-zero sigma2 in case of CTF-premultiplication!
					if (ctf_premultiplied)
						DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n) = XMIPP_MAX(DIRECT_MULTIDIM_ELEM(mymodel.sigma2_noise[igroup], n), 1e-15);

				}
			}
		}
	}
	RCTIC(timer,RCT_7);
	RCTOC(timer,RCT_8);
	// After the first iteration the references are always CTF-corrected
    if (do_ctf_correction)
    	refs_are_ctf_corrected = true;

	// Some statistics to output
    mymodel.LL = 	wsum_model.LL;
	if ((iter==1 && do_firstiter_cc) || do_always_cc)
		mymodel.LL /= sum_weight; // this now stores the average ccf
	mymodel.ave_Pmax = wsum_model.ave_Pmax / sum_weight;

	// After the first, special iteration, apply low-pass filter of -ini_high again
	if (iter == 1 && do_firstiter_cc)
	{
		initialLowPassFilterReferences();
		if (ini_high > 0.)
		{
			// Adjust the tau2_class and data_vs_prior_class, because they were calculated on the unfiltered maps
			// This is merely a matter of having correct output in the model.star file (these values are not used in the calculations)
			RFLOAT radius = mymodel.ori_size * mymodel.pixel_size / ini_high;
			radius -= WIDTH_FMASK_EDGE / 2.;
			RFLOAT radius_p = radius + WIDTH_FMASK_EDGE;

			for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
			{
				for (int rr = 0; rr < XSIZE(mymodel.tau2_class[iclass]); rr++)
				{
					RFLOAT r = (RFLOAT)rr;
					if (r < radius)
						continue;
					else if (r > radius_p)
					{
						DIRECT_A1D_ELEM(mymodel.tau2_class[iclass], rr) = 0.;
						DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[iclass], rr) = 0.;
					}
					else
					{
						RFLOAT raisedcos = 0.5 - 0.5 * cos(PI * (radius_p - r) / WIDTH_FMASK_EDGE);
						DIRECT_A1D_ELEM(mymodel.tau2_class[iclass], rr) *= raisedcos * raisedcos;
						DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[iclass], rr) *= raisedcos * raisedcos;
					}
				}
			}
		}

		if (do_generate_seeds && mymodel.nr_classes > 1)
		{
			// In the first CC-iteration only a single reference was used
			// Now copy this one reference to all K references, for seed generation in the second iteration
			for (int iclass = 1; iclass < mymodel.nr_classes; iclass++)
			{
				mymodel.tau2_class[iclass] =  mymodel.tau2_class[0];
				mymodel.data_vs_prior_class[iclass] = mymodel.data_vs_prior_class[0];
				mymodel.pdf_class[iclass] = mymodel.pdf_class[0] / mymodel.nr_classes;
				mymodel.pdf_direction[iclass] = mymodel.pdf_direction[0];
				mymodel.Iref[iclass] = mymodel.Iref[0];
			}
			mymodel.pdf_class[0] /= mymodel.nr_classes;
		}

	}
	RCTOC(timer,RCT_8);
#ifdef DEBUG
	std::cerr << "Leaving maximizationOtherParameters" << std::endl;
#endif
}


void MlOptimiser::solventFlatten()
{
#ifdef DEBUG
	std::cerr << "Entering MlOptimiser::solventFlatten" << std::endl;
#endif
	// First read solvent mask from disc, or pre-calculate it
	Image<RFLOAT> Isolvent, Isolvent2;
    Isolvent().resize(mymodel.Iref[0]);
	Isolvent().setXmippOrigin();
	Isolvent().initZeros();
	if (fn_mask.contains("None"))
	{
		// Jun09,2015 - Shaoda, Helical refinement
		// Solvent flatten for helices has already been done in 'makeHelicalReferenceInRealSpace()'
		if (do_helical_refine)
		{
			if (ignore_helical_symmetry)
			{
				createCylindricalReference(
						Isolvent(),
						XSIZE(mymodel.Iref[0]),
						helical_tube_inner_diameter / mymodel.pixel_size,
						helical_tube_outer_diameter / mymodel.pixel_size,
						width_mask_edge);
			}
			else
			{
				FOR_ALL_ELEMENTS_IN_ARRAY3D(Isolvent())
				{
					A3D_ELEM(Isolvent(), k, i, j) = 1.;
				}
			}
		}
		else
		{
			RFLOAT radius = (particle_diameter / (2. * mymodel.pixel_size));
			RFLOAT radius_p = radius + width_mask_edge;
			FOR_ALL_ELEMENTS_IN_ARRAY3D(Isolvent())
			{
				RFLOAT r = sqrt((RFLOAT)(k * k + i * i + j * j));
				if (r < radius)
					A3D_ELEM(Isolvent(), k, i, j) = 1.;
				else if (r > radius_p)
					A3D_ELEM(Isolvent(), k, i, j) = 0.;
				else
					A3D_ELEM(Isolvent(), k, i, j) = 0.5 - 0.5 * cos(PI * (radius_p - r) / width_mask_edge );
			}
		}
	}
	else
	{
		Isolvent.read(fn_mask);
		Isolvent().setXmippOrigin();

		if (Isolvent().computeMin() < 0. || Isolvent().computeMax() > 1.)
			REPORT_ERROR("MlOptimiser::solventFlatten: ERROR solvent mask should contain values between 0 and 1 only...");
	}

	// Also read a second solvent mask if necessary
	if (!fn_mask2.contains("None"))
	{
		Isolvent2.read(fn_mask2);
		Isolvent2().setXmippOrigin();
		if (!Isolvent2().sameShape(Isolvent()))
			REPORT_ERROR("MlOptimiser::solventFlatten ERROR: second solvent mask is of incorrect size.");
	}

	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
	{

		// Then apply the expanded solvent mask to the map
		mymodel.Iref[iclass] *= Isolvent();

		// Apply a second solvent mask if necessary
		// This may for example be useful to set the interior of icosahedral viruses to a constant density value that is higher than the solvent
		// Invert the solvent mask, so that an input mask can be given where 1 is the masked area and 0 is protein....
		if (!fn_mask2.contains("None"))
			softMaskOutsideMap(mymodel.Iref[iclass], Isolvent2(), true);

	} // end for iclass
#ifdef DEBUG
	std::cerr << "Leaving MlOptimiser::solventFlatten" << std::endl;
#endif

}

void MlOptimiser::updateCurrentResolution()
{
//#define DEBUG
#ifdef DEBUG
	std::cerr << "Entering MlOptimiser::updateCurrentResolution" << std::endl;
#endif


	int maxres = 0;
	if (do_map)
	{
		// Set current resolution
		if (ini_high > 0. && (iter == 0 || (iter == 1 && do_firstiter_cc)))
		{
			maxres = ROUND(mymodel.ori_size * mymodel.pixel_size / ini_high);
		}
		else
		{
			// Calculate at which resolution shell the data_vs_prior drops below 1
			int ires;
			for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
			{
				for (ires = 1; ires < mymodel.ori_size/2; ires++)
				{
					if (DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[iclass], ires) < 1.)
						break;
				}
				// Subtract one shell to be back on the safe side
				ires--;
				if (ires > maxres)
					maxres = ires;
			}

			// Never allow smaller maxres than minres_map
			maxres = XMIPP_MAX(maxres, minres_map);
		}
	}
	else
	{
		// If we are not doing MAP-estimation, set maxres to Nyquist
		maxres = mymodel.ori_size/2;
	}
    RFLOAT newres = mymodel.getResolution(maxres);

    // Check whether resolution improved, if not increase nr_iter_wo_resol_gain
    //if (newres <= best_resol_thus_far)
    if (newres <= mymodel.current_resolution+0.0001) // Add 0.0001 to avoid problems due to rounding error
    	nr_iter_wo_resol_gain++;
    else
    	nr_iter_wo_resol_gain = 0;

    // Store best resolution thus far (but no longer do anything with it anymore...)
    if (newres > best_resol_thus_far)
    	best_resol_thus_far = newres;

    mymodel.current_resolution = newres;

}

void MlOptimiser::updateImageSizeAndResolutionPointers()
{

	// Increment the current_size
    // If we are far from convergence (in the initial stages of refinement) take steps of 25% the image size
    // Do this whenever the FSC at the current_size is larger than 0.2, but NOT when this is in combination with very low Pmax values,
    // in the latter case, over-marginalisation may lead to spuriously high FSCs (2 smoothed maps may look very similar at high-res: all zero!)
    //
	int maxres = mymodel.getPixelFromResolution(mymodel.current_resolution);
	if (mymodel.ave_Pmax > 0.1 && has_high_fsc_at_limit)
    {
		maxres += ROUND(0.25 * mymodel.ori_size / 2);
    }
	else
	{
		// If we are near our resolution limit, use incr_size (by default 10 shells)
		maxres += incr_size;
	}

	// In SGD, optionally restrict maximum resolution to a fixed value
	if (do_sgd && strict_highres_sgd > 0.)
	{
		maxres = XMIPP_MIN(maxres, mymodel.getPixelFromResolution(1./strict_highres_sgd) );
	}
    // Go back from resolution shells (i.e. radius) to image size, which are BTW always even...
	mymodel.current_size = maxres * 2;

	// If realigning movies: go all the way because resolution increase may be substantial
	if (do_use_all_data)
		mymodel.current_size = mymodel.ori_size;

	// current_size can never be larger than ori_size:
	mymodel.current_size = XMIPP_MIN(mymodel.current_size, mymodel.ori_size);

	// The current size is also used in wsum_model (in unpacking)
	wsum_model.current_size = mymodel.current_size;

	// Update coarse_size
	if (strict_highres_exp > 0.)
    {
    	// Strictly limit the coarse size to the one corresponding to strict_highres_exp
    	coarse_size = 2 * ROUND(mymodel.ori_size * mymodel.pixel_size / strict_highres_exp);
    }
    else if (adaptive_oversampling > 0.)
	{
    	// Dependency of coarse_size on the angular sampling used in the first pass
    	RFLOAT rotated_distance = (sampling.getAngularSampling() / 360.) * PI * particle_diameter;
		RFLOAT keepsafe_factor = (mymodel.ref_dim == 3) ? 1.2 : 1.5;
		RFLOAT coarse_resolution = rotated_distance / keepsafe_factor;
		// Note coarse_size should be even-valued!
		coarse_size = 2 * CEIL(mymodel.pixel_size * mymodel.ori_size / coarse_resolution);
		// Coarse size can never be larger than max_coarse_size
		coarse_size = XMIPP_MIN(max_coarse_size, coarse_size);
	}
	else
		coarse_size = mymodel.current_size;
    // Coarse_size can never become bigger than current_size
    coarse_size = XMIPP_MIN(mymodel.current_size, coarse_size);

	/// Also update the resolution pointers here

	// Calculate number of pixels per resolution shell
	Npix_per_shell.initZeros(mymodel.ori_size / 2 + 1);
	MultidimArray<RFLOAT> aux;
	if (mymodel.data_dim == 3)
		aux.resize(mymodel.ori_size, mymodel.ori_size, mymodel.ori_size / 2 + 1);
	else
		aux.resize(mymodel.ori_size, mymodel.ori_size / 2 + 1);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(aux)
	{
		int ires = ROUND(sqrt((RFLOAT)(kp*kp + ip*ip + jp*jp)));
		// TODO: better check for volume_refine, but the same still seems to hold... Half of the yz plane (either ip<0 or kp<0 is redundant at jp==0)
		// Exclude points beyond XSIZE(Npix_per_shell), and exclude half of the x=0 column that is stored twice in FFTW
		if (ires < mymodel.ori_size / 2 + 1 && !(jp==0 && ip < 0))
			Npix_per_shell(ires) += 1;
	}

	if (mymodel.data_dim == 3)
		Mresol_fine.resize(mymodel.current_size, mymodel.current_size, mymodel.current_size / 2 + 1);
	else
		Mresol_fine.resize(mymodel.current_size, mymodel.current_size / 2 + 1);
	Mresol_fine.initConstant(-1);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Mresol_fine)
	{
		int ires = ROUND(sqrt((RFLOAT)(kp*kp + ip*ip + jp*jp)));
		// TODO: better check for volume_refine, but the same still seems to hold... Half of the yz plane (either ip<0 or kp<0 is redundant at jp==0)
		// Exclude points beyond ires, and exclude and half (y<0) of the x=0 column that is stored twice in FFTW
		if (ires < mymodel.current_size / 2 + 1  && !(jp==0 && ip < 0))
		{
			DIRECT_A3D_ELEM(Mresol_fine, k, i, j) = ires;
		}
	}

	if (mymodel.data_dim == 3)
		Mresol_coarse.resize(coarse_size, coarse_size, coarse_size/ 2 + 1);
	else
		Mresol_coarse.resize(coarse_size, coarse_size/ 2 + 1);

	Mresol_coarse.initConstant(-1);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Mresol_coarse)
	{
		int ires = ROUND(sqrt((RFLOAT)(kp*kp + ip*ip + jp*jp)));
		// Exclude points beyond ires, and exclude and half (y<0) of the x=0 column that is stored twice in FFTW
		// exclude lowest-resolution points
		if (ires < coarse_size / 2 + 1 && !(jp==0 && ip < 0))
		{
			DIRECT_A3D_ELEM(Mresol_coarse, k, i, j) = ires;
		}
	}

//#define DEBUG_MRESOL
#ifdef DEBUG_MRESOL
	Image<RFLOAT> img;
	img().resize(YSIZE(Mresol_fine),XSIZE(Mresol_fine));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
	{
		DIRECT_MULTIDIM_ELEM(img(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
	}
	img.write("Mresol_fine.mrc");
	img().resize(YSIZE(Mresol_coarse),XSIZE(Mresol_coarse));
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
	{
		DIRECT_MULTIDIM_ELEM(img(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(Mresol_coarse, n);
	}
	img.write("Mresol_coarse.mrc");
#endif


#ifdef DEBUG
	std::cerr << " current_size= " << mymodel.current_size << " coarse_size= " << coarse_size << " current_resolution= " << mymodel.current_resolution << std::endl;
	std::cerr << "Leaving MlOptimiser::updateCurrentResolution" << std::endl;
#endif

}


void MlOptimiser::calculateRunningAveragesOfMovieFrames(long int my_ori_particle,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<RFLOAT> > &exp_power_imgs,
		std::vector<RFLOAT> &exp_highres_Xi2_imgs)
{

	std::vector<MultidimArray<Complex > > runavg_Fimgs;
	std::vector<int> count_runavg;
	MultidimArray<Complex > Fzero;
	Fzero.resize(exp_Fimgs[0]);
	Fzero.initZeros();
	runavg_Fimgs.resize(exp_Fimgs.size(), Fzero);
	count_runavg.resize(exp_Fimgs.size(), 0);

	// Calculate the running sums
	for (int iframe = 0; iframe < exp_Fimgs.size(); iframe++)
	{

		long int my_first_runavg_frame = XMIPP_MAX(0, iframe - movie_frame_running_avg_side);
		long int my_last_runavg_frame = XMIPP_MIN(exp_Fimgs.size() - 1, iframe + movie_frame_running_avg_side);

		// Run over all images again and see which ones to sum
		for (int iframe2 = my_first_runavg_frame; iframe2 <=  my_last_runavg_frame; iframe2++)
		{
			// Add to sum
			runavg_Fimgs[iframe] += exp_Fimgs[iframe2];
			count_runavg[iframe] += 1;
		}
	}

	// Calculate averages from sums and set back in exp_ vectors
	for (int iframe = 0; iframe < exp_Fimgs.size(); iframe++)
	{
		RFLOAT sum = (RFLOAT)count_runavg[iframe];
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_Fimgs[iframe])
		{
			DIRECT_MULTIDIM_ELEM(exp_Fimgs[iframe], n) = DIRECT_MULTIDIM_ELEM(runavg_Fimgs[iframe], n) / sum;
		}
		// Also lower the power of the images for the sigma2_noise and diff2 calculations beyond current_size....
		// sigma2_(a+b) = sigma2_(a) + sigma2_(b)
		// The actual values are lost, just hope the images obey statistics...
		exp_power_imgs[iframe] /= sum;
		exp_highres_Xi2_imgs[iframe] /= sum;
	}

}

void MlOptimiser::getFourierTransformsAndCtfs(long int my_ori_particle, int ibody, int metadata_offset,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
		std::vector<MultidimArray<RFLOAT> > &exp_Fctfs,
		std::vector<Matrix1D<RFLOAT> > &exp_old_offset,
		std::vector<Matrix1D<RFLOAT> > &exp_prior,
		std::vector<MultidimArray<RFLOAT> > &exp_power_imgs,
		std::vector<RFLOAT> &exp_highres_Xi2_imgs,
		std::vector<int> &exp_pointer_dir_nonzeroprior,
		std::vector<int> &exp_pointer_psi_nonzeroprior,
		std::vector<RFLOAT> &exp_directions_prior,
		std::vector<RFLOAT> &exp_psi_prior)
{

	FourierTransformer transformer;

	for (int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		FileName fn_img;
		Image<RFLOAT> img, rec_img;
		MultidimArray<Complex > Fimg, Faux;
		MultidimArray<RFLOAT> Fctf;

		// Get the right line in the exp_fn_img strings (also exp_fn_recimg and exp_fn_ctfs)
		int istop = 0;
		for (long int ii = exp_my_first_ori_particle; ii < my_ori_particle; ii++)
			istop += mydata.ori_particles[ii].particles_id.size();
		istop += ipart;

		// What is my particle_id?
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		// Which group do I belong?
		int group_id = mydata.getGroupId(part_id);

		// Get the norm_correction
		int icol_norm = (mymodel.nr_bodies == 1) ? METADATA_NORM : 6 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
		RFLOAT normcorr = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_norm);

		// Get the old offsets and the priors on the offsets
		Matrix1D<RFLOAT> my_old_offset(mymodel.data_dim), my_prior(mymodel.data_dim);
		int icol_xoff, icol_yoff, icol_zoff;
		if (mymodel.nr_bodies == 1)
		{
			XX(my_old_offset) = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_XOFF);
			XX(my_prior)      = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_XOFF_PRIOR);
			YY(my_old_offset) = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_YOFF);
			YY(my_prior)      = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_YOFF_PRIOR);
			if (mymodel.data_dim == 3)
			{
				ZZ(my_old_offset) = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ZOFF);
				ZZ(my_prior)      = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ZOFF_PRIOR);
			}
		}
		else
		{
			icol_xoff = 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
			icol_yoff = 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
			icol_zoff = 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;

			XX(my_old_offset) = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_XOFF) + DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_xoff);
			XX(my_prior)      = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_XOFF_PRIOR);
			YY(my_old_offset) = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_YOFF) + DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_yoff);
			YY(my_prior)      = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_YOFF_PRIOR);
			if (mymodel.data_dim == 3)
			{
				ZZ(my_old_offset) = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ZOFF) + DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_zoff);
				ZZ(my_prior)      = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ZOFF_PRIOR);
			}
		}

		// Uninitialised priors were set to 999.
		if (XX(my_prior) > 998.99 && XX(my_prior) < 999.01)
			XX(my_prior) = 0.;
		if (YY(my_prior) > 998.99 && YY(my_prior) < 999.01)
			YY(my_prior) = 0.;
		if (mymodel.data_dim == 3 && ZZ(my_prior) > 998.99 && ZZ(my_prior) < 999.01)
			ZZ(my_prior) = 0.;

		if (mymodel.orientational_prior_mode != NOPRIOR && !(do_skip_align || do_skip_rotate || do_only_sample_tilt))
		{
			// First try if there are some fixed prior angles
			// For multi-body refinements, have the prior angles set to zero
			RFLOAT prior_rot =  (mymodel.nr_bodies > 1 ) ? 0. : DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT_PRIOR);
			RFLOAT prior_tilt = (mymodel.nr_bodies > 1 ) ? 0. : DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT_PRIOR);
			RFLOAT prior_psi =  (mymodel.nr_bodies > 1 ) ? 0. : DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI_PRIOR);
			RFLOAT prior_psi_flip_ratio =  (mymodel.nr_bodies > 1 ) ? 0. : DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI_PRIOR_FLIP_RATIO);

			bool do_auto_refine_local_searches = (do_auto_refine) && (sampling.healpix_order >= autosampling_hporder_local_searches);
			bool do_classification_local_searches = (!do_auto_refine) && (mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI)
					&& (mymodel.sigma2_rot > 0.) && (mymodel.sigma2_tilt > 0.) && (mymodel.sigma2_psi > 0.);
			bool do_local_angular_searches = (do_auto_refine_local_searches) || (do_classification_local_searches);

			// If there were no defined priors (i.e. their values were 999.), then use the "normal" angles
			if (prior_rot > 998.99 && prior_rot < 999.01)
				prior_rot = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT);
			if (prior_tilt > 998.99 && prior_tilt < 999.01)
				prior_tilt = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT);
			if ( (do_helical_refine) && (helical_keep_tilt_prior_fixed) && (do_local_angular_searches) )
				prior_tilt = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT);
			if (prior_psi > 998.99 && prior_psi < 999.01)
				prior_psi = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI);
			if (prior_psi_flip_ratio > 998.99 && prior_psi_flip_ratio < 999.01)
				prior_psi_flip_ratio = 0.5;

			////////// How does this work now: each particle has a different sampling object?!!!
			// Select only those orientations that have non-zero prior probability
			// Jun04,2015 - Shaoda & Sjors, bimodal psi searches for helices
			if (do_helical_refine)
			{
				sampling.selectOrientationsWithNonZeroPriorProbabilityFor3DHelicalReconstruction(prior_rot, prior_tilt, prior_psi,
										sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
										exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior,
										do_local_angular_searches, prior_psi_flip_ratio);
			}
			else
			{
				sampling.selectOrientationsWithNonZeroPriorProbability(prior_rot, prior_tilt, prior_psi,
										sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
										exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior,
										((do_bimodal_psi) && (mymodel.sigma2_psi > 0.)) );
			}

			long int nr_orients = sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior) * sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
			if (nr_orients == 0)
			{
				std::cerr << " sampling.NrDirections()= " << sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior)
						<< " sampling.NrPsiSamplings()= " << sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior) << std::endl;
				REPORT_ERROR("Zero orientations fall within the local angular search. Increase the sigma-value(s) on the orientations!");
			}

		}

		// Get the image and recimg data
		if (do_parallel_disc_io)
		{
			// If all slaves had preread images into RAM: get those now
			if (do_preread_images)
			{
                img().reshape(mydata.particles[part_id].img);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mydata.particles[part_id].img)
				{
                	DIRECT_MULTIDIM_ELEM(img(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(mydata.particles[part_id].img, n);
				}
			}
			else
			{
//#define DEBUG_SIMULTANEOUS_READ
#ifdef DEBUG_SIMULTANEOUS_READ
				// Read from disc
				FileName fn_img;
				std::istringstream split(exp_fn_img);
				for (int i = 0; i <= istop; i++)
					getline(split, fn_img);

				img.read(fn_img);
				img().setXmippOrigin();

				// Check that this is the same as the image in exp_imgs vector
				Image<RFLOAT> diff;
				if (istop >= exp_imgs.size())
				{
					std::cerr << " istop= " <<istop << " exp_imgs.size()= "<<exp_imgs.size()<<std::endl;
					REPORT_ERROR("BUG: istop runs out of bounds!");
				}
				diff() = img() - exp_imgs[istop];
				if (diff().computeMax() > 1e-6)
				{
					std::cerr << "istop= " <<istop<<" fn_img=" << fn_img << " diff avg=" << diff().computeAvg()<<std::endl;
					diff.write("diff.spi");
					diff()= exp_imgs[istop];
					diff.write("preread.spi");
					img.write("img.spi");
					REPORT_ERROR("unequal pre-read images... BUG!");
				}
#else
				if (mymodel.data_dim == 3)
				{
					// Read sub-tomograms from disc in parallel (to save RAM in exp_imgs)
					FileName fn_img;
					if (!mydata.getImageNameOnScratch(part_id, fn_img))
					{
						std::istringstream split(exp_fn_img);
						for (int i = 0; i <= istop; i++)
							getline(split, fn_img);
					}
					img.read(fn_img);
					img().setXmippOrigin();
				}
				else
				{
					img() = exp_imgs[istop];
				}
#endif
			}
			if (has_converged && do_use_reconstruct_images)
			{
				FileName fn_recimg;
				std::istringstream split2(exp_fn_recimg);
				// Get the right line in the exp_fn_img string
				for (int i = 0; i <= istop; i++)
					getline(split2, fn_recimg);
				rec_img.read(fn_recimg);
				rec_img().setXmippOrigin();
			}
		}
		else
		{

			// Unpack the image from the imagedata
			if (mymodel.data_dim == 3)
			{
				img().resize(mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
				// Only allow a single image per call of this function!!! nr_pool needs to be set to 1!!!!
				// This will save memory, as we'll need to store all translated images in memory....
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
				{
					DIRECT_A3D_ELEM(img(), k, i, j) = DIRECT_A3D_ELEM(exp_imagedata, k, i, j);
				}
				img().setXmippOrigin();

				if (has_converged && do_use_reconstruct_images)
				{
					rec_img().resize(mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
					int offset = (do_ctf_correction) ? 2 * mymodel.ori_size : mymodel.ori_size;
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(rec_img())
					{
						DIRECT_A3D_ELEM(rec_img(), k, i, j) = DIRECT_A3D_ELEM(exp_imagedata, offset + k, i, j);
					}
					rec_img().setXmippOrigin();

				}

			}
			else
			{
				img().resize(mymodel.ori_size, mymodel.ori_size);
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img())
				{
					DIRECT_A2D_ELEM(img(), i, j) = DIRECT_A3D_ELEM(exp_imagedata, metadata_offset + ipart, i, j);
				}
				img().setXmippOrigin();
				if (has_converged && do_use_reconstruct_images)
				{

					////////////// TODO: think this through for no-threads here.....
					rec_img().resize(mymodel.ori_size, mymodel.ori_size);
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(rec_img())
					{
						DIRECT_A2D_ELEM(rec_img(), i, j) = DIRECT_A3D_ELEM(exp_imagedata, exp_nr_images + metadata_offset + ipart, i, j);
					}
					rec_img().setXmippOrigin();
				}
			}
		}

//#define DEBUG_SOFTMASK
#ifdef DEBUG_SOFTMASK
		Image<RFLOAT> tt;
		tt()=img();
		tt.write("Fimg_unmasked.spi");
		std::cerr << "written Fimg_unmasked.spi; press any key to continue..." << std::endl;
		char c;
		std::cin >> c;
#endif
		// Apply the norm_correction term
		if (do_norm_correction)
		{
//#define DEBUG_NORM
#ifdef DEBUG_NORM
			if (normcorr < 0.001 || normcorr > 1000. || mymodel.avg_norm_correction < 0.001 || mymodel.avg_norm_correction > 1000.)
			{
				std::cerr << " ** normcorr= " << normcorr << std::endl;
				std::cerr << " ** mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
				std::cerr << " ** fn_img= " << fn_img << " part_id= " << part_id << std::endl;
				std::cerr << " ** iseries= " << iseries << " ipart= " << ipart << " part_id= " << part_id << std::endl;
				int group_id = mydata.getGroupId(part_id);
				std::cerr << " ml_model.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << " group_id= " << group_id <<std::endl;
				std::cerr << " part_id= " << part_id << " iseries= " << iseries << std::endl;
				std::cerr << " img_id= " << img_id << std::endl;
				REPORT_ERROR("Very small or very big (avg) normcorr!");
			}
#endif
			img() *= mymodel.avg_norm_correction / normcorr;
		}

		// Helical reconstruction: calculate old_offset in the system of coordinates of the helix, i.e. parallel & perpendicular, depending on psi-angle!
		// For helices do NOT apply old_offset along the direction of the helix!!
		Matrix1D<RFLOAT> my_old_offset_helix_coords;
		RFLOAT rot_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT);
		RFLOAT tilt_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT);
		RFLOAT psi_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI);
		if ( (do_helical_refine) && (!ignore_helical_symmetry) )
		{
			// Calculate my_old_offset_helix_coords from my_old_offset and psi angle
			transformCartesianAndHelicalCoords(my_old_offset, my_old_offset_helix_coords, rot_deg, tilt_deg, psi_deg, CART_TO_HELICAL_COORDS);
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
			// May 18, 2015 - Shaoda & Sjors - Helical refinement (orientational searches)
			std::cerr << "MlOptimiser::getFourierTransformsAndCtfs()" << std::endl;
			std::cerr << " Transform old Cartesian offsets to helical ones..." << std::endl;
			if(my_old_offset.size() == 2)
			{
				std::cerr << "  psi_deg = " << psi_deg << " degrees" << std::endl;
				std::cerr << "  old_offset(x, y) = (" << XX(my_old_offset) << ", " << YY(my_old_offset) << ")" << std::endl;
				std::cerr << "  old_offset_helix(r, p) = (" << XX(my_old_offset_helix_coords) << ", " << YY(my_old_offset_helix_coords) << ")" << std::endl;
			}
			else
			{
				std::cerr << "  psi_deg = " << psi_deg << " degrees, tilt_deg = " << tilt_deg << " degrees"<< std::endl;
				std::cerr << "  old_offset(x, y, z) = (" << XX(my_old_offset) << ", " << YY(my_old_offset) << ", " << ZZ(my_old_offset) << ")" << std::endl;
				std::cerr << "  old_offset_helix(p1, p2, z) = (" << XX(my_old_offset_helix_coords) << ", " << YY(my_old_offset_helix_coords) << "," << ZZ(my_old_offset_helix_coords) << ")" << std::endl;
			}
#endif
			// We do NOT want to accumulate the offsets in the direction along the helix (which is X in the helical coordinate system!)
			// However, when doing helical local searches, we accumulate offsets
			// Do NOT accumulate offsets in 3D classification of helices
			if ( (mymodel.ref_dim == 3) && (!do_skip_align) && (!do_skip_rotate) )
			{
				// TODO: check whether the following lines make sense
				bool do_auto_refine_local_searches = (do_auto_refine) && (sampling.healpix_order >= autosampling_hporder_local_searches);
				bool do_classification_local_searches = (!do_auto_refine) && (mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI)
						&& (mymodel.sigma2_rot > 0.) && (mymodel.sigma2_tilt > 0.) && (mymodel.sigma2_psi > 0.);
				bool do_local_angular_searches = (do_auto_refine_local_searches) || (do_classification_local_searches);
				if (!do_local_angular_searches)
				{
					if (mymodel.data_dim == 2)
						XX(my_old_offset_helix_coords) = 0.;
					else if (mymodel.data_dim == 3)
						ZZ(my_old_offset_helix_coords) = 0.;
				}
			}
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
			std::cerr << " Set r (translation along helical axis) to zero..." << std::endl;
			if(my_old_offset.size() == 2)
				std::cerr << "  old_offset_helix(r, p) = (" << XX(my_old_offset_helix_coords) << ", " << YY(my_old_offset_helix_coords) << ")" << std::endl;
			else
				std::cerr << "  old_offset_helix(p1, p2, z) = (" << XX(my_old_offset_helix_coords) << ", " << YY(my_old_offset_helix_coords) << "," << ZZ(my_old_offset_helix_coords) << ")" << std::endl;
#endif
			// Now re-calculate the my_old_offset in the real (or image) system of coordinate (rotate -psi angle)
			transformCartesianAndHelicalCoords(my_old_offset_helix_coords, my_old_offset, rot_deg, tilt_deg, psi_deg, HELICAL_TO_CART_COORDS);
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
			std::cerr << " Transform helical offsets back to Cartesian ones..." << std::endl;
			if(my_old_offset.size() == 2)
				std::cerr << "  old_offset(x, y) = (" << XX(my_old_offset) << ", " << YY(my_old_offset) << ")" << std::endl;
			else
				std::cerr << "  old_offset(x, y, z) = (" << XX(my_old_offset) << ", " << YY(my_old_offset) << ", " << ZZ(my_old_offset) << ")" << std::endl;
#endif
		}

		my_old_offset.selfROUND();
		selfTranslate(img(), my_old_offset, DONT_WRAP);
		if (has_converged && do_use_reconstruct_images)
		{
			selfTranslate(rec_img(), my_old_offset, DONT_WRAP);
		}

#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
		if ( (do_helical_refine) && (!ignore_helical_symmetry) )
		{
			std::cerr << " Apply (rounded) old offsets (r = 0, p) & (psi, tilt) for helices..." << std::endl;
			if(my_old_offset.size() == 2)
				std::cerr << "  old_offset(x, y) = (" << XX(my_old_offset) << ", " << YY(my_old_offset) << ")" << std::endl;
			else
				std::cerr << "  old_offset(x, y, z) = (" << XX(my_old_offset) << ", " << YY(my_old_offset) << ", " << ZZ(my_old_offset) << ")" << std::endl;
			Image<RFLOAT> tt;
			tt = img;
			tt.write("selftranslated_helix.spi");
			tt.clear();
			std::cerr << " written selftranslated_helix.spi; press any key to continue..." << std::endl;
			std::string str;
			std::cin >> str;
		}
#endif

		if ( (do_helical_refine) && (!ignore_helical_symmetry) )
		{
			// Transform rounded Cartesian offsets to corresponding helical ones
			transformCartesianAndHelicalCoords(my_old_offset, my_old_offset_helix_coords, rot_deg, tilt_deg, psi_deg, CART_TO_HELICAL_COORDS);
			exp_old_offset[ipart] = my_old_offset_helix_coords;
		}
		else
		{
			exp_old_offset[ipart] = my_old_offset;  // Not doing helical refinement. Rounded Cartesian offsets are stored.
		}
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
		if ( (do_helical_refine) && (!ignore_helical_symmetry) )
		{
			std::cerr << " ipart = " << ipart << ", " << std::flush;
			if(exp_old_offset[ipart].size() == 2)
				std::cerr << "exp_old_offset = (" << XX(exp_old_offset[ipart]) << ", " << YY(exp_old_offset[ipart]) << ")" << std::endl;
			else
				std::cerr << "exp_old_offset = (" << XX(exp_old_offset[ipart]) << ", " << YY(exp_old_offset[ipart]) << ", " << ZZ(exp_old_offset[ipart]) << ")" << std::endl;
		}
#endif
		// Also store priors on translations
		exp_prior[ipart] = my_prior;

		// Always store FT of image without mask (to be used for the reconstruction)
		MultidimArray<RFLOAT> img_aux;
		img_aux = (has_converged && do_use_reconstruct_images) ? rec_img() : img();
		CenterFFT(img_aux, true);
		transformer.FourierTransform(img_aux, Faux);
		windowFourierTransform(Faux, Fimg, mymodel.current_size);

		// Here apply the beamtilt correction if necessary
		// This will only be used for reconstruction, not for alignment
		// But beamtilt only affects very high-resolution components anyway...
		//
		RFLOAT beamtilt_x = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_BEAMTILT_X);
		RFLOAT beamtilt_y = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_BEAMTILT_Y);
		RFLOAT Cs = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_CS);
		RFLOAT V = 1000. * DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_VOLTAGE);
		RFLOAT lambda = 12.2643247 / sqrt(V * (1. + V * 0.978466e-6));
		if (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.)
			selfApplyBeamTilt(Fimg, beamtilt_x, beamtilt_y, lambda, Cs, mymodel.pixel_size, mymodel.ori_size);
		exp_Fimgs_nomask.at(ipart) = Fimg;

		MultidimArray<RFLOAT> Mnoise;
		bool is_helical_segment = (do_helical_refine) || ((mymodel.ref_dim == 2) && (helical_tube_outer_diameter > 0.));
		if (!do_zero_mask)
		{
			// Make a noisy background image with the same spectrum as the sigma2_noise

			// Different MPI-distributed subsets may otherwise have different instances of the random noise below,
			// because work is on an on-demand basis and therefore variable with the timing of distinct nodes...
			// Have the seed based on the part_id, so that each particle has a different instant of the noise
			if (do_realign_movies)
				init_random_generator(random_seed + part_id);
			else
				init_random_generator(random_seed + my_ori_particle); // This only serves for exact reproducibility tests with 1.3-code...

			// If we're doing running averages, then the sigma2_noise was already adjusted for the running averages.
			// Undo this adjustment here in order to get the right noise in the individual frames
			MultidimArray<RFLOAT> power_noise = sigma2_fudge * mymodel.sigma2_noise[group_id];
			if (do_realign_movies)
				power_noise *= (2. * movie_frame_running_avg_side + 1.);

			// Create noisy image for outside the mask
			MultidimArray<Complex > Fnoise;
			Mnoise.resize(img());
			transformer.setReal(Mnoise);
			transformer.getFourierAlias(Fnoise);
			// Fill Fnoise with random numbers, use power spectrum of the noise for its variance
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fnoise)
			{
				int ires = ROUND( sqrt( (RFLOAT)(kp * kp + ip * ip + jp * jp) ) );
				if (ires >= 0 && ires < XSIZE(Fnoise))
				{
					RFLOAT sigma = sqrt(DIRECT_A1D_ELEM(power_noise, ires));
					DIRECT_A3D_ELEM(Fnoise, k, i, j).real = rnd_gaus(0., sigma);
					DIRECT_A3D_ELEM(Fnoise, k, i, j).imag = rnd_gaus(0., sigma);
				}
				else
				{
					DIRECT_A3D_ELEM(Fnoise, k, i, j) = 0.;
				}
			}
			// Back to real space Mnoise
			transformer.inverseFourierTransform();
			Mnoise.setXmippOrigin();

			// May24,2014 - Shaoda & Sjors, Helical refinement
			if (is_helical_segment)
			{
				softMaskOutsideMapForHelix(img(), psi_deg, tilt_deg, (particle_diameter / (2. * mymodel.pixel_size)),
						(helical_tube_outer_diameter / (2. * mymodel.pixel_size)), width_mask_edge, &Mnoise);
			}
			else
				softMaskOutsideMap(img(), particle_diameter / (2. * mymodel.pixel_size), (RFLOAT)width_mask_edge, &Mnoise);
		}
		else
		{
			// May24,2014 - Shaoda & Sjors, Helical refinement
			if (is_helical_segment)
			{
				softMaskOutsideMapForHelix(img(), psi_deg, tilt_deg, (particle_diameter / (2. * mymodel.pixel_size)),
						(helical_tube_outer_diameter / (2. * mymodel.pixel_size)), width_mask_edge);
			}
			else
				softMaskOutsideMap(img(), particle_diameter / (2. * mymodel.pixel_size), (RFLOAT)width_mask_edge);
		}
#ifdef DEBUG_SOFTMASK
		tt()=img();
		tt.write("Fimg_masked.spi");
		std::cerr << "written Fimg_masked.spi; press any key to continue..." << std::endl;
		std::cin >> c;
#endif

		// Inside Projector and Backprojector the origin of the Fourier Transform is centered!
		CenterFFT(img(), true);

		// Store the Fourier Transform of the image Fimg
		transformer.FourierTransform(img(), Faux);

		// Store the power_class spectrum of the whole image (to fill sigma2_noise between current_size and ori_size
		if (mymodel.current_size < mymodel.ori_size)
		{
			MultidimArray<RFLOAT> spectrum;
			spectrum.initZeros(mymodel.ori_size/2 + 1);
			RFLOAT highres_Xi2 = 0.;
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
			{
				int ires = ROUND( sqrt( (RFLOAT)(kp*kp + ip*ip + jp*jp) ) );
				// Skip Hermitian pairs in the x==0 column

				if (ires > 0 && ires < mymodel.ori_size/2 + 1 && !(jp==0 && ip < 0) )
				{
					RFLOAT normFaux = norm(DIRECT_A3D_ELEM(Faux, k, i, j));
					DIRECT_A1D_ELEM(spectrum, ires) += normFaux;
					// Store sumXi2 from current_size until ori_size
					if (ires >= mymodel.current_size/2 + 1)
						highres_Xi2 += normFaux;
				}
			}

			// Let's use .at() here instead of [] to check whether we go outside the vectors bounds
			exp_power_imgs.at(ipart) = spectrum;
			exp_highres_Xi2_imgs.at(ipart) = highres_Xi2;
		}
		else
		{
			exp_highres_Xi2_imgs.at(ipart) = 0.;
		}

		// We never need any resolutions higher than current_size
		// So resize the Fourier transforms
		windowFourierTransform(Faux, Fimg, mymodel.current_size);

		// Store Fimg
		exp_Fimgs.at(ipart) = Fimg;

		// Also store its CTF
		Fctf.resize(Fimg);

		// Now calculate the actual CTF
		if (do_ctf_correction)
		{
			if (mymodel.data_dim == 3)
			{
				Image<RFLOAT> Ictf;
				if (do_parallel_disc_io)
				{
					// Read CTF-image from disc
					FileName fn_ctf;
					if (!mydata.getImageNameOnScratch(part_id, fn_ctf, true))
					{
						std::istringstream split(exp_fn_ctf);
						// Get the right line in the exp_fn_img string
						for (int i = 0; i <= istop; i++)
							getline(split, fn_ctf);
					}
					Ictf.read(fn_ctf);
				}
				else
				{
					// Unpack the CTF-image from the exp_imagedata array
					Ictf().resize(mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Ictf())
					{
						DIRECT_A3D_ELEM(Ictf(), k, i, j) = DIRECT_A3D_ELEM(exp_imagedata, mymodel.ori_size + k, i, j);
					}
				}
				// Set the CTF-image in Fctf
				Ictf().setXmippOrigin();
				FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf)
				{
					// Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
					DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
				}
			}
			else
			{
				CTF ctf;
				ctf.setValues(DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_DEFOCUS_U),
							  DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_DEFOCUS_V),
							  DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_DEFOCUS_ANGLE),
							  DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_VOLTAGE),
							  DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_CS),
							  DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_Q0),
							  DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_BFAC),
							  1.,
							  DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CTF_PHASE_SHIFT));

				ctf.getFftwImage(Fctf, mymodel.ori_size, mymodel.ori_size, mymodel.pixel_size,
						ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
			}
//#define DEBUG_CTF_FFTW_IMAGE
#ifdef DEBUG_CTF_FFTW_IMAGE
			Image<RFLOAT> tt;
			tt()=Fctf;
			tt.write("relion_ctf.spi");
			std::cerr << "Written relion_ctf.spi, now exiting..." << std::endl;
			exit(1);
#endif
//#define DEBUG_GETCTF
#ifdef DEBUG_GETCTF
			std::cerr << " intact_ctf_first_peak= " << intact_ctf_first_peak << std::endl;
			ctf.write(std::cerr);
			Image<RFLOAT> tmp;
			tmp()=Fctf;
			tmp.write("Fctf.spi");
			tmp().resize(mymodel.ori_size, mymodel.ori_size);
			ctf.getCenteredImage(tmp(), mymodel.pixel_size, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
			tmp.write("Fctf_cen.spi");
			std::cerr << "Written Fctf.spi, Fctf_cen.spi. Press any key to continue..." << std::endl;
			char c;
			std::cin >> c;
#endif
		}
		else
		{
			Fctf.initConstant(1.);
		}

		// Store Fctf
		exp_Fctfs.at(ipart) = Fctf;

		// If we're doing multibody refinement, now subtract projections of the other bodies from both the masked and the unmasked particle
		if (mymodel.nr_bodies > 1)
		{
			// Aori is the matrix for the original angles in the 'consensus' refinement
			Matrix2D<RFLOAT> Aori;
			Euler_angles2matrix(DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT),
								DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT),
								DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI), Aori, false);

			for (int obody = 0; obody < mymodel.nr_bodies; obody++)
			{
				if (obody != ibody) // other body is not this body....
				{
					// Get the right metadata
					int ocol_rot  = 0 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_tilt = 1 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_psi  = 2 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_xoff = 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_yoff = 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_zoff = 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_norm = 6 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;

					Matrix2D<RFLOAT> Aresi,  Abody;
					// Aresi is the residual orientation for this obody
					Euler_angles2matrix(DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, ocol_rot),
										DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, ocol_tilt),
										DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, ocol_psi), Aresi, false);
					// The real orientation to be applied is the obody transformation applied after the original one
					Abody = Aresi * Aori;

					// Get the FT of the projection in the right direction
					MultidimArray<Complex> FTo;
					FTo.initZeros(Fimg);
					mymodel.PPref[obody].get2DFourierTransform(FTo, Abody, IS_NOT_INV);

#ifdef DEBUG_BODIES
					if (my_ori_particle == 0)
					{
						std::cerr << " A= " << A << std::endl;
						//std::cerr << " mymodel.com_bodies[obody]= " << mymodel.com_bodies[obody] << " comp= " << comp << std::endl;
						//FTo.resize(Fimg);
						//FTo.initZeros();
						//mymodel.PPref[0].get2DFourierTransform(FTo, A, IS_NOT_INV);
						windowFourierTransform(FTo, Faux, mymodel.ori_size);
						transformer.inverseFourierTransform(Faux, img());
						CenterFFT(img(), false);
						FileName fn_img = "FTo_ini.spi";
						fn_img = fn_img.insertBeforeExtension("_body" + integerToString(ibody+1));
						img.write(fn_img);
						std::cerr << "Written::: " << fn_img << std::endl;

						windowFourierTransform(exp_Fimgs_nomask.at(ipart), Faux, mymodel.ori_size);
						transformer.inverseFourierTransform(Faux, img());
						CenterFFT(img(), false);
						fn_img = "img_unmask";
						fn_img = fn_img.insertBeforeExtension("_body" + integerToString(ibody+1));
						img.write(fn_img);
						std::cerr << "Written" << fn_img << std::endl;
					}
					else
					{
						std::cerr << "sleeping!!" << std::endl;
						sleep(2000000);
						std::cerr << "waking!!!" << std::endl;
					}
#endif

					// Also apply the correct translation
					//21may2015
					Faux = FTo;
					RFLOAT xoff, yoff, zoff;
					xoff = XX(my_old_offset) - (DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_XOFF) + DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, ocol_xoff));
					yoff = YY(my_old_offset) - (DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_YOFF) + DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, ocol_yoff));
					zoff = (mymodel.data_dim == 2) ? 0. :
							ZZ(my_old_offset) - (DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ZOFF) + DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, ocol_zoff));


					//19may2015: now body is centered at its own com!
					// Put it back to where it was!
					Matrix1D<RFLOAT> comp(3);
					comp = Abody * (mymodel.com_bodies[obody]);
					xoff += XX(comp);
					yoff += YY(comp);
					zoff += ZZ(comp);

					shiftImageInFourierTransform(Faux, FTo, (RFLOAT)mymodel.ori_size, xoff, yoff, zoff);

#ifdef DEBUG_BODIES
					if (my_ori_particle == 0)
					{
						windowFourierTransform(FTo, Faux, mymodel.ori_size);
						transformer.inverseFourierTransform(Faux, img());
						CenterFFT(img(), false);
						FileName fn_img = "FTo.spi";
						fn_img = fn_img.insertBeforeExtension("_body" + integerToString(ibody+1));
						img.write(fn_img);
						std::cerr << "written " << fn_img << std::endl;
					}
#endif

					// Apply the CTF to this reference projection
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FTo)
					{
						DIRECT_MULTIDIM_ELEM(FTo, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}

					// Subtract the other-body FT from the current image FT
					exp_Fimgs.at(ipart) -= FTo;
					exp_Fimgs_nomask.at(ipart)  -= FTo;
				}
			}

			if (do_parallel_disc_io)
			{
//#define DEBUG_WRITE_BODY_DIFFS
#ifdef DEBUG_WRITE_BODY_DIFFS


				// Write out the corrected particles. Assume the mrcs files exist
				FileName fn_img;
				std::istringstream split(exp_fn_img);
				for (int i = 0; i <= istop; i++)
					getline(split, fn_img);

				fn_img = fn_img.insertBeforeExtension("_body" + integerToString(ibody+1));
				FileName fn_tmp;
				long ii;
				fn_img.decompose(ii, fn_tmp);
				windowFourierTransform(exp_Fimgs_nomask.at(ipart), Faux, mymodel.ori_size);
				transformer.inverseFourierTransform(Faux, img());
				CenterFFT(img(), false);
				img.write(fn_tmp, ii-1, true, WRITE_REPLACE);
#endif
#ifdef DEBUG_BODIES
				if (my_ori_particle==0)
				{
					img.write("diff_body"+ integerToString(ibody+1)+".spi");
					std::cerr << " written fn_img= " << fn_img << " and " << "diff_body"<< integerToString(ibody+1)<<".spi"<< " ibody= " << ibody << std::endl;
				}
#endif
			}
		}

	} // end loop ipart
	transformer.clear();

}

void MlOptimiser::precalculateShiftedImagesCtfsAndInvSigma2s(bool do_also_unmasked,
		long int my_ori_particle, int exp_current_image_size, int exp_current_oversampling, int metadata_offset,
		int exp_itrans_min, int exp_itrans_max,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
		std::vector<MultidimArray<RFLOAT> > &exp_Fctfs,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
		std::vector<MultidimArray<RFLOAT> > &exp_local_Fctfs,
		std::vector<RFLOAT> &exp_local_sqrtXi2,
		std::vector<MultidimArray<RFLOAT> > &exp_local_Minvsigma2s)
{

#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
	{
		if (do_also_unmasked)
			timer.tic(TIMING_ESP_PRECW);
		else if (exp_current_oversampling == 0) timer.tic(TIMING_ESP_PREC1);
		else timer.tic(TIMING_ESP_PREC2);
	}
#endif

	int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	int nr_shifts = (do_shifts_onthefly || do_skip_align) ? exp_nr_particles : exp_nr_particles * sampling.NrTranslationalSamplings(exp_current_oversampling);
	// Don't re-do if nothing has changed....
	bool do_ctf_invsig = (exp_local_Fctfs.size() > 0) ? YSIZE(exp_local_Fctfs[0])  != exp_current_image_size : true; // size has changed
	bool do_masked_shifts = (do_ctf_invsig || nr_shifts != exp_local_Fimgs_shifted.size()); // size or nr_shifts has changed

	// Use pre-sized vectors instead of push_backs!!
	exp_local_Fimgs_shifted.resize(nr_shifts);
	if (do_also_unmasked)
		exp_local_Fimgs_shifted_nomask.resize(nr_shifts);
	exp_local_Minvsigma2s.resize(exp_nr_particles);
	exp_local_Fctfs.resize(exp_nr_particles);
	exp_local_sqrtXi2.resize(exp_nr_particles);

	MultidimArray<Complex > Fimg, Fimg_nomask;
	for (int ipart = 0, my_trans_image = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		int group_id = mydata.getGroupId(part_id);

		if (do_masked_shifts)
			windowFourierTransform(exp_Fimgs[ipart], Fimg, exp_current_image_size);
		if (do_also_unmasked)
			windowFourierTransform(exp_Fimgs_nomask[ipart], Fimg_nomask, exp_current_image_size);

		if (do_ctf_invsig)
		{
			// Also precalculate the sqrt of the sum of all Xi2
			// Could exp_current_image_size ever be different from mymodel.current_size?
			// Probably therefore do it here rather than in getFourierTransforms
			if ((iter == 1 && do_firstiter_cc) || do_always_cc)
			{
				RFLOAT sumxi2 = 0.;
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
				{
					sumxi2 += norm(DIRECT_MULTIDIM_ELEM(Fimg, n));
				}
				// Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
				exp_local_sqrtXi2[ipart] = sqrt(sumxi2);
			}

			// Also store downsized Fctfs
			// In the second pass of the adaptive approach this will have no effect,
			// since then exp_current_image_size will be the same as the size of exp_Fctfs
#ifdef DEBUG_CHECKSIZES
			if (ipart >= exp_local_Fctfs.size())
			{
				std::cerr<< "ipart= "<<ipart<<" exp_local_Fctfs.size()= "<< exp_local_Fctfs.size() <<std::endl;
				REPORT_ERROR("ipart >= exp_local_Fctfs.size()");
			}
#endif
			windowFourierTransform(exp_Fctfs[ipart], exp_local_Fctfs[ipart], exp_current_image_size);

			// Also prepare Minvsigma2
#ifdef DEBUG_CHECKSIZES
			if (ipart >= exp_local_Minvsigma2s.size())
			{
				std::cerr<< "ipart= "<<ipart<<" exp_local_Minvsigma2s.size()= "<< exp_local_Minvsigma2s.size() <<std::endl;
				REPORT_ERROR("ipart >= exp_local_Minvsigma2s.size()");
			}
#endif
			if (mymodel.data_dim == 3)
				exp_local_Minvsigma2s[ipart].initZeros(ZSIZE(Fimg), YSIZE(Fimg), XSIZE(Fimg));
			else
				exp_local_Minvsigma2s[ipart].initZeros(YSIZE(Fimg), XSIZE(Fimg));

			int *myMresol = (YSIZE(Fimg) == coarse_size) ? Mresol_coarse.data : Mresol_fine.data;
			// With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Minvsigma2s[ipart])
			{
				int ires = *(myMresol + n);
				// Exclude origin (ires==0) from the Probability-calculation
				// This way we are invariant to additive factors
				if (ires > 0)
					DIRECT_MULTIDIM_ELEM(exp_local_Minvsigma2s[ipart], n) = 1. / (sigma2_fudge * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], ires));
			}
		}

		//Shifts are done on the fly on the gpu
		if (do_shifts_onthefly && !do_gpu)
		{
			// Store a single, down-sized version of exp_Fimgs[ipart] in exp_local_Fimgs_shifted
#ifdef DEBUG_CHECKSIZES
			if (ipart >= exp_local_Fimgs_shifted.size())
			{
				std::cerr<< "ipart= "<<ipart<<" exp_local_Fimgs_shifted.size()= "<< exp_local_Fimgs_shifted.size() <<std::endl;
				REPORT_ERROR("ipart >= exp_local_Fimgs_shifted.size()");
			}
#endif
			if (do_masked_shifts)
				exp_local_Fimgs_shifted[ipart] = Fimg;
			if (do_also_unmasked)
				exp_local_Fimgs_shifted_nomask[ipart] = Fimg_nomask;
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
			std::cerr << " MlOptimiser::precalculateShiftedImagesCtfsAndInvSigma2s(): do_shifts_onthefly && !do_gpu" << std::endl;
#endif
		}
		else if(!do_gpu)
		{
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
			Image<RFLOAT> img_save_ori, img_save_mask, img_save_nomask;
			img_save_ori.clear();
			img_save_ori().resize((mymodel.data_dim == 3) ? (mymodel.ori_size) : (1), mymodel.ori_size, mymodel.ori_size);
			img_save_ori().initZeros();
			img_save_mask.clear();
			img_save_mask().resize((mymodel.data_dim == 3) ? (mymodel.ori_size) : (1), mymodel.ori_size, mymodel.ori_size);
			img_save_mask().initZeros();
			img_save_nomask.clear();
			img_save_nomask().resize((mymodel.data_dim == 3) ? (mymodel.ori_size) : (1), mymodel.ori_size, mymodel.ori_size);
			img_save_nomask().initZeros();
#endif
			// Store all translated variants of Fimg
			for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++)
			{

				// First get the non-oversampled translations as defined by the sampling object
				std::vector<RFLOAT > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
				// Jun01,2014 - Shaoda & Sjors, Helical refinement
				sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
						(do_helical_refine) && (!ignore_helical_symmetry), helical_rise_initial / mymodel.pixel_size, helical_twist_initial);
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
				std::cerr << "MlOptimiser::precalculateShiftedImagesCtfsAndInvSigma2s(): Store all translated variants of Fimg" << std::endl;
#endif
				// Then loop over all its oversampled relatives
				for (long int iover_trans = 0; iover_trans < oversampled_translations_x.size(); iover_trans++, my_trans_image++)
				{
					// Helical reconstruction: rotate oversampled_translations_x[iover_trans] and oversampled_translations_y[iover_trans] according to rlnAnglePsi of this particle!
					RFLOAT xshift = 0., yshift = 0., zshift = 0.;

					xshift = oversampled_translations_x[iover_trans];
					yshift = oversampled_translations_y[iover_trans];
					if (mymodel.data_dim == 3)
						zshift = oversampled_translations_z[iover_trans];

					if ( (do_helical_refine) && (!ignore_helical_symmetry) )
					{
						RFLOAT rot_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT);
						RFLOAT tilt_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT);
						RFLOAT psi_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI);
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
						std::cerr << "Helical xyz shift = (" << xshift << ", " << yshift << ", " << zshift << ")" << std::endl;
#endif
						transformCartesianAndHelicalCoords(xshift, yshift, zshift, xshift, yshift, zshift, rot_deg, tilt_deg, psi_deg, mymodel.data_dim, HELICAL_TO_CART_COORDS);
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
						std::cerr << "Cartesian xyz shift = (" << xshift << ", " << yshift << ", " << zshift << ")" << std::endl;
#endif
					}

					// Shift through phase-shifts in the Fourier transform
					// Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
					if (do_masked_shifts)
					{
						exp_local_Fimgs_shifted[my_trans_image].resize(Fimg);
						shiftImageInFourierTransform(Fimg, exp_local_Fimgs_shifted[my_trans_image], (RFLOAT)mymodel.ori_size, xshift, yshift, zshift);
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
						if ( (do_helical_refine) && (!ignore_helical_symmetry) )  // Shall we let 2D classification do this as well?
						{
							std::cerr << " Size of Fourier map (Z, Y, X) = "
									<< ZSIZE(exp_local_Fimgs_shifted[my_trans_image]) << ", "
									<< YSIZE(exp_local_Fimgs_shifted[my_trans_image]) << ", "
									<< XSIZE(exp_local_Fimgs_shifted[my_trans_image]) << std::endl;
							std::cerr << " mymodel.ori_size = " << mymodel.ori_size << std::endl;
							MultidimArray<Complex> Faux, Fo;
							Image<RFLOAT> tt;
							FourierTransformer transformer;
							tt().resize((mymodel.data_dim == 3) ? (mymodel.ori_size) : (1), mymodel.ori_size, mymodel.ori_size);
							Faux = exp_local_Fimgs_shifted[my_trans_image];
							windowFourierTransform(Faux, Fo, mymodel.ori_size);
							transformer.inverseFourierTransform(Fo, tt());
							CenterFFT(tt(), false);
							img_save_mask() += tt();
							img_save_mask.write("translational_searches_mask_helix.spi");
							std::cerr << " written translational_searches_mask_helix.spi; press any key to continue..." << std::endl;
							std::string str;
							std::cin >> str;
						}
#endif
					}
					if (do_also_unmasked)
					{
						exp_local_Fimgs_shifted_nomask[my_trans_image].resize(Fimg_nomask);
						shiftImageInFourierTransform(Fimg_nomask, exp_local_Fimgs_shifted_nomask[my_trans_image], (RFLOAT)mymodel.ori_size, xshift, yshift, zshift);
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
						if ( (do_helical_refine) && (!ignore_helical_symmetry) )
						{
							std::cerr << " Size of Fourier map (Z, Y, X) = "
									<< ZSIZE(exp_local_Fimgs_shifted_nomask[my_trans_image]) << ", "
									<< YSIZE(exp_local_Fimgs_shifted_nomask[my_trans_image]) << ", "
									<< XSIZE(exp_local_Fimgs_shifted_nomask[my_trans_image]) << std::endl;
							std::cerr << " mymodel.ori_size = " << mymodel.ori_size << std::endl;
							MultidimArray<Complex> Faux, Fo;
							Image<RFLOAT> tt;
							FourierTransformer transformer;
							tt().resize((mymodel.data_dim == 3) ? (mymodel.ori_size) : (1), mymodel.ori_size, mymodel.ori_size);
							Faux = exp_local_Fimgs_shifted_nomask[my_trans_image];
							windowFourierTransform(Faux, Fo, mymodel.ori_size);
							transformer.inverseFourierTransform(Fo, tt());
							CenterFFT(tt(), false);
							img_save_nomask() += tt();
							img_save_nomask.write("translational_searches_nomask_helix.spi");
							std::cerr << " written translational_searches_nomask_helix.spi; press any key to continue..." << std::endl;
							std::string str;
							std::cin >> str;
						}
#endif
					}
				}
			}
		}
	}
#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
	{
		if (do_also_unmasked)
			timer.toc(TIMING_ESP_PRECW);
		else if (exp_current_oversampling == 0) timer.toc(TIMING_ESP_PREC1);
		else timer.toc(TIMING_ESP_PREC2);
	}
#endif

}

bool MlOptimiser::isSignificantAnyParticleAnyTranslation(long int iorient, int exp_itrans_min, int exp_itrans_max, MultidimArray<bool> &exp_Mcoarse_significant)
{

	long int exp_nr_trans = exp_itrans_max - exp_itrans_min + 1;
	for (long int ipart = 0; ipart < YSIZE(exp_Mcoarse_significant); ipart++)
	{
		long int ihidden = iorient * exp_nr_trans;
		for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++)
		{
#ifdef DEBUG_CHECKSIZES
			if (ihidden >= XSIZE(exp_Mcoarse_significant))
			{
				std::cerr << " ihidden= " << ihidden << " XSIZE(exp_Mcoarse_significant)= " << XSIZE(exp_Mcoarse_significant) << std::endl;
				std::cerr << " iorient= " << iorient << " itrans= " << itrans << " exp_nr_trans= " << exp_nr_trans << std::endl;
				REPORT_ERROR("ihidden > XSIZE: ");
			}
#endif
			if (DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden))
				return true;
		}
	}
	return false;

}


void MlOptimiser::getAllSquaredDifferences(long int my_ori_particle, int ibody,  int exp_current_image_size,
		int exp_ipass, int exp_current_oversampling, int metadata_offset,
		int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
		int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
		std::vector<RFLOAT> &exp_min_diff2,
		std::vector<RFLOAT> &exp_highres_Xi2_imgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<RFLOAT> > &exp_Fctfs,
		MultidimArray<RFLOAT> &exp_Mweight,
		MultidimArray<bool> &exp_Mcoarse_significant,
		std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
		std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
		std::vector<MultidimArray<RFLOAT> > &exp_local_Minvsigma2s,
		std::vector<MultidimArray<RFLOAT> > &exp_local_Fctfs,
		std::vector<RFLOAT> &exp_local_sqrtXi2)
{

#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
	{
		if (exp_ipass == 0) timer.tic(TIMING_ESP_DIFF1);
		else timer.tic(TIMING_ESP_DIFF2);
	}
#endif

//#define DEBUG_GETALLDIFF2
#ifdef DEBUG_GETALLDIFF2
	std::cerr << " ipass= " << exp_ipass << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
	std::cerr << " sampling.NrPsiSamplings(exp_current_oversampling)= " << sampling.NrPsiSamplings(exp_current_oversampling) << std::endl;
	std::cerr << " sampling.NrTranslationalSamplings(exp_current_oversampling)= " << sampling.NrTranslationalSamplings(exp_current_oversampling) << std::endl;
	std::cerr << " sampling.NrSamplingPoints(exp_current_oversampling)= " << sampling.NrSamplingPoints(exp_current_oversampling) << std::endl;
	std::cerr << " sampling.oversamplingFactorOrientations(exp_current_oversampling)= "<<sampling.oversamplingFactorOrientations(exp_current_oversampling) << std::endl;
	std::cerr << " sampling.oversamplingFactorTranslations(exp_current_oversampling)= "<<sampling.oversamplingFactorTranslations(exp_current_oversampling) << std::endl;
#endif

	// Initialise min_diff and exp_Mweight for this pass

	int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	long int exp_nr_dir = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
	long int exp_nr_psi = (do_skip_align || do_skip_rotate || do_only_sample_tilt) ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
	long int exp_nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings();
	long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
	long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

	exp_Mweight.resize(exp_nr_particles, mymodel.nr_classes * exp_nr_dir * exp_nr_psi * exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans);
	exp_Mweight.initConstant(-999.);
	if (exp_ipass==0)
		exp_Mcoarse_significant.clear();

	exp_min_diff2.clear();
	exp_min_diff2.resize(exp_nr_particles, LARGE_NUMBER);

	std::vector<MultidimArray<Complex > > dummy;
	precalculateShiftedImagesCtfsAndInvSigma2s(false, my_ori_particle, exp_current_image_size, exp_current_oversampling, metadata_offset,
			exp_itrans_min, exp_itrans_max, exp_Fimgs, dummy, exp_Fctfs, exp_local_Fimgs_shifted, dummy,
			exp_local_Fctfs, exp_local_sqrtXi2, exp_local_Minvsigma2s);

	// Loop only from exp_iclass_min to exp_iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{
		if (mymodel.pdf_class[exp_iclass] > 0.)
		{
			// Local variables
			std::vector< RFLOAT > oversampled_rot, oversampled_tilt, oversampled_psi;
			std::vector< RFLOAT > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
			MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_otfshift;
			RFLOAT *Minvsigma2;
			Matrix2D<RFLOAT> A, Abody, Aori;

			if (mymodel.nr_bodies > 1)
			{
				// ipart=0 because in multi-body refinement we do not do movie frames!
				RFLOAT rot_ori = DIRECT_A2D_ELEM(exp_metadata, metadata_offset, METADATA_ROT);
				RFLOAT tilt_ori = DIRECT_A2D_ELEM(exp_metadata, metadata_offset, METADATA_TILT);
				RFLOAT psi_ori = DIRECT_A2D_ELEM(exp_metadata, metadata_offset, METADATA_PSI);
				Euler_angles2matrix(rot_ori, tilt_ori, psi_ori, Aori, false);
			}

			Fref.resize(exp_local_Minvsigma2s[0]);
			Frefctf.resize(exp_local_Minvsigma2s[0]);
			if (do_shifts_onthefly)
				Fimg_otfshift.resize(Frefctf);

            for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
			{
				for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++)
				{
					long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;

					// Get prior for this direction and skip calculation if prior==0
					RFLOAT pdf_orientation;
					if (do_skip_align || do_skip_rotate)
					{
						#ifdef DEBUG_CHECKSIZES
						if (exp_iclass >= mymodel.pdf_class.size())
						{
							std::cerr<< "exp_iclass= "<<exp_iclass<<" mymodel.pdf_class.size()= "<< mymodel.pdf_class.size() <<std::endl;
							REPORT_ERROR("exp_iclass >= mymodel.pdf_class.size()");
						}
						#endif
						pdf_orientation = mymodel.pdf_class[exp_iclass];
					}
					else if (mymodel.orientational_prior_mode == NOPRIOR)
					{
#ifdef DEBUG_CHECKSIZES
						if (idir >= XSIZE(mymodel.pdf_direction[exp_iclass]))
						{
							std::cerr<< "idir= "<<idir<<" XSIZE(mymodel.pdf_direction[exp_iclass])= "<< XSIZE(mymodel.pdf_direction[exp_iclass]) <<std::endl;
							REPORT_ERROR("idir >= mymodel.pdf_direction[exp_iclass].size()");
						}
#endif
						pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
					}
					else
					{
						pdf_orientation = exp_directions_prior[idir] * exp_psi_prior[ipsi];
					}
					// In the first pass, always proceed
					// In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
					// if so, proceed with projecting the reference in that direction
					bool do_proceed = (exp_ipass==0) ? true :
						isSignificantAnyParticleAnyTranslation(iorientclass, exp_itrans_min, exp_itrans_max, exp_Mcoarse_significant);
					if (do_proceed && pdf_orientation > 0.)
					{
						// Now get the oversampled (rot, tilt, psi) triplets
						// This will be only the original (rot,tilt,psi) triplet in the first pass (exp_current_oversampling==0)
						sampling.getOrientations(idir, ipsi, exp_current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
								exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior);
						// Loop over all oversampled orientations (only a single one in the first pass)
						for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
						{
							// Get the Euler matrix
							Euler_angles2matrix(oversampled_rot[iover_rot],
												oversampled_tilt[iover_rot],
												oversampled_psi[iover_rot], A, false);

							// For multi-body refinements, A are only 'residual' orientations, Abody is the complete Euler matrix
							if (mymodel.nr_bodies > 1)
								Abody = A * Aori;

							// Project the reference map (into Fref)
#ifdef TIMING
							// Only time one thread, as I also only time one MPI process
							if (my_ori_particle == exp_my_first_ori_particle)
								timer.tic(TIMING_DIFF_PROJ);
#endif
							if (mymodel.nr_bodies > 1)
							{
								(mymodel.PPref[ibody]).get2DFourierTransform(Fref, Abody, IS_NOT_INV);
								// 19may2015: inside PPref, each body is centered at its own COM
								// Put it back to where it was!
								Matrix1D<RFLOAT> comp(3);
								comp = Abody * (mymodel.com_bodies[ibody]);
								shiftImageInFourierTransform(Fref, Fref, (RFLOAT)mymodel.ori_size, XX(comp), YY(comp), ZZ(comp));
							}
							else
								(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);

#ifdef TIMING
							// Only time one thread, as I also only time one MPI process
							if (my_ori_particle == exp_my_first_ori_particle)
								timer.toc(TIMING_DIFF_PROJ);
#endif
							/// Now that reference projection has been made loop over someParticles!
							// loop over all particles inside this ori_particle
							for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
							{
#ifdef DEBUG_CHECKSIZES
								if (my_ori_particle >= mydata.ori_particles.size())
								{
									std::cerr<< "my_ori_particle= "<<my_ori_particle<<" mydata.ori_particles.size()= "<< mydata.ori_particles.size() <<std::endl;
									REPORT_ERROR("my_ori_particle >= mydata.ori_particles.size()");
								}
								if (ipart >= mydata.ori_particles[my_ori_particle].particles_id.size())
								{
									std::cerr<< "ipart= "<<ipart<<" mydata.ori_particles[my_ori_particle].particles_id.size()= "<< mydata.ori_particles[my_ori_particle].particles_id.size() <<std::endl;
									REPORT_ERROR("ipart >= mydata.ori_particles[my_ori_particle].particles_id.size()");
								}
#endif
								long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
#ifdef DEBUG_CHECKSIZES
								if (ipart >= exp_local_Minvsigma2s.size())
								{
									std::cerr<< "ipart= "<<ipart<<" exp_local_Minvsigma2s.size()= "<< exp_local_Minvsigma2s.size() <<std::endl;
									REPORT_ERROR("ipart >= exp_local_Minvsigma2s.size()");
								}
#endif
								Minvsigma2 = exp_local_Minvsigma2s[ipart].data;

								// Apply CTF to reference projection
								if (do_ctf_correction && refs_are_ctf_corrected)
								{
									if (ctf_premultiplied)
									{
										// TODO: ignore CTF until first peak of premultiplied CTF?
										FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
										{
											DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n) * DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n);
										}
									}
									else
									{
										FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
										{
											DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], n);
										}
									}
								}
								else
									Frefctf = Fref;

								if (do_scale_correction)
								{
									int group_id = mydata.getGroupId(part_id);
#ifdef DEBUG_CHECKSIZES
									if (group_id >= mymodel.scale_correction.size())
									{
										std::cerr<< "group_id= "<<group_id<<" mymodel.scale_correction.size()= "<< mymodel.scale_correction.size() <<std::endl;
										REPORT_ERROR("group_id >= mymodel.scale_correction.size()");
									}
#endif
									RFLOAT myscale = mymodel.scale_correction[group_id];
									FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
									{
										DIRECT_MULTIDIM_ELEM(Frefctf, n) *= myscale;
									}
								}
								long int ihidden = iorientclass * exp_nr_trans;
								for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++)
								{
#ifdef DEBUG_CHECKSIZES
									if (exp_ipass > 0 && ihidden >= XSIZE(exp_Mcoarse_significant))
									{
										std::cerr<< "ihidden= "<<ihidden<<" XSIZE(exp_Mcoarse_significant)= "<< XSIZE(exp_Mcoarse_significant) <<std::endl;
										REPORT_ERROR("ihidden >= XSIZE(exp_Mcoarse_significant)");
									}
#endif
									// In the first pass, always proceed
									// In the second pass, check whether this translations (&orientation) had a significant weight in the first pass
									bool do_proceed = (exp_ipass == 0) ? true : DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden);
									if (do_proceed)
									{
										// Jun01,2015 - Shaoda & Sjors, Helical refinement
										sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
												(do_helical_refine) && (!ignore_helical_symmetry), helical_rise_initial / mymodel.pixel_size, helical_twist_initial);
										for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++)
										{
#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
												timer.tic(TIMING_DIFF2_GETSHIFT);
#endif
											/// Now get the shifted image
											// Use a pointer to avoid copying the entire array again in this highly expensive loop
											Complex *Fimg_shift;
											if (!do_shifts_onthefly)
											{
												long int ishift = ipart * exp_nr_oversampled_trans * exp_nr_trans +
														(itrans - exp_itrans_min) * exp_nr_oversampled_trans + iover_trans;
												if (do_skip_align)
													ishift = ipart;
#ifdef DEBUG_CHECKSIZES
												if (ishift >= exp_local_Fimgs_shifted.size())
												{
													std::cerr<< "ishift= "<<ishift<<" exp_local_Fimgs_shifted.size()= "<< exp_local_Fimgs_shifted.size() <<std::endl;
													std::cerr << " itrans= " << itrans << std::endl;
													std::cerr << " ipart= " << ipart << std::endl;
													std::cerr << " exp_nr_oversampled_trans= " << exp_nr_oversampled_trans << " exp_nr_trans= " << exp_nr_trans << " iover_trans= " << iover_trans << std::endl;
													REPORT_ERROR("ishift >= exp_local_Fimgs_shifted.size()");
												}
#endif
												Fimg_shift = exp_local_Fimgs_shifted[ishift].data;
											}
											else
											{
												// Calculate shifted image on-the-fly to save replicating memory in multi-threaded jobs.
												// Feb01,2017 - Shaoda, on-the-fly shifts in helical reconstuctions (2D and 3D)
												if ( (do_helical_refine) && (!ignore_helical_symmetry) )
												{
													bool use_coarse_size = false;
													RFLOAT xshift = 0., yshift = 0., zshift = 0.;

													xshift = (exp_current_oversampling == 0) ? (oversampled_translations_x[0]) : (oversampled_translations_x[iover_trans]);
													yshift = (exp_current_oversampling == 0) ? (oversampled_translations_y[0]) : (oversampled_translations_y[iover_trans]);
													if (mymodel.data_dim == 3)
														zshift = (exp_current_oversampling == 0) ? (oversampled_translations_z[0]) : (oversampled_translations_z[iover_trans]);

													RFLOAT rot_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT);
													RFLOAT tilt_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT);
													RFLOAT psi_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI);
													transformCartesianAndHelicalCoords(
																xshift, yshift, zshift,
																xshift, yshift, zshift,
																rot_deg, tilt_deg, psi_deg,
																mymodel.data_dim,
																HELICAL_TO_CART_COORDS);

													use_coarse_size = ((exp_current_oversampling == 0) && (YSIZE(Frefctf) == coarse_size)) || ((exp_current_oversampling > 0) && (strict_highres_exp > 0.));
													shiftImageInFourierTransformWithTabSincos(
															exp_local_Fimgs_shifted[ipart],
															Fimg_otfshift,
															(RFLOAT)mymodel.ori_size,
															(use_coarse_size) ? (coarse_size) : (mymodel.current_size),
															tab_sin, tab_cos,
															xshift, yshift, zshift);
												}
												else
												{
													Complex *myAB;
													if (exp_current_oversampling == 0)
													{
														#ifdef DEBUG_CHECKSIZES
														if (YSIZE(Frefctf) == coarse_size && itrans >= global_fftshifts_ab_coarse.size())
														{
															std::cerr<< "itrans= "<<itrans<<" global_fftshifts_ab_coarse.size()= "<< global_fftshifts_ab_coarse.size() <<std::endl;
															REPORT_ERROR("itrans >= global_fftshifts_ab_coarse.size()");
														}
														if (YSIZE(Frefctf) != coarse_size && itrans >= global_fftshifts_ab_current.size())
														{
															std::cerr<< "itrans= "<<itrans<<" global_fftshifts_ab_current.size()= "<< global_fftshifts_ab_current.size() <<std::endl;
															REPORT_ERROR("itrans >= global_fftshifts_ab_current.size()");
														}
														#endif
														myAB = (YSIZE(Frefctf) == coarse_size) ? global_fftshifts_ab_coarse[itrans].data
														        : global_fftshifts_ab_current[itrans].data;
													}
													else
													{
														int iitrans = itrans * exp_nr_oversampled_trans +  iover_trans;
														myAB = (strict_highres_exp > 0.) ? global_fftshifts_ab2_coarse[iitrans].data
																: global_fftshifts_ab2_current[iitrans].data;
													}
													FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
													{
														RFLOAT real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real
																- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag;
														RFLOAT imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag
																+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real;
														DIRECT_MULTIDIM_ELEM(Fimg_otfshift, n) = Complex(real, imag);
													}
												}
												Fimg_shift = Fimg_otfshift.data;
											}
#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
												timer.toc(TIMING_DIFF2_GETSHIFT);
#endif

//#define DEBUG_GETALLDIFF2
#ifdef DEBUG_GETALLDIFF2
											pthread_mutex_lock(&global_mutex);
											//if (verb> 0)
											{
												std::cerr << " A= " << A << std::endl;

												FourierTransformer transformer;
												MultidimArray<Complex> Fish;
												Fish.resize(exp_local_Minvsigma2s[0]);
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fish)
												{
													DIRECT_MULTIDIM_ELEM(Fish, n) = *(Fimg_shift + n);
												}
												Image<RFLOAT> tt;
												if (mymodel.data_dim == 3)
													tt().resize(exp_current_image_size, exp_current_image_size, exp_current_image_size);
												else
													tt().resize(exp_current_image_size, exp_current_image_size);
												transformer.inverseFourierTransform(Fish, tt());
												CenterFFT(tt(),false);
												tt.write("Fimg_shift.spi");

												transformer.inverseFourierTransform(Frefctf, tt());
												CenterFFT(tt(),false);
												tt.write("Fref.spi");
												char c;
												std::cerr << " ipart " << ipart << " DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], 12)= " << DIRECT_MULTIDIM_ELEM(exp_local_Fctfs[ipart], 12) << std::endl;
												std::cerr << " ipart " << ipart << " DIRECT_MULTIDIM_ELEM(exp_Fctfs[ipart], 12)= " << DIRECT_MULTIDIM_ELEM(exp_Fctfs[ipart], 12) << std::endl;

												int group_id = mydata.getGroupId(part_id);
												for (int i = 0; i< mymodel.scale_correction.size(); i++)
													std::cerr << i << " scale="<<mymodel.scale_correction[i]<<std::endl;
												RFLOAT myscale = mymodel.scale_correction[group_id];
												std::cerr << " oversampled_rot[iover_rot]= " << oversampled_rot[iover_rot] << " oversampled_tilt[iover_rot]= " << oversampled_tilt[iover_rot] << " oversampled_psi[iover_rot]= " << oversampled_psi[iover_rot] << std::endl;
												std::cerr << " group_id= " << group_id << " myscale= " << myscale <<std::endl;
												std::cerr << " itrans= " << itrans << " itrans * exp_nr_oversampled_trans +  iover_trans= " << itrans * exp_nr_oversampled_trans +  iover_trans << " ihidden= " << ihidden << std::endl;
												std::cerr <<" my_ori_particle= "<<my_ori_particle<<" name= "<< mydata.ori_particles[my_ori_particle].name << std::endl;
												//std::cerr << " myrank= "<< myrank<<std::endl;
												std::cerr << "Written Fimg_shift.spi and Fref.spi. Press any key to continue... my_ori_particle= " << my_ori_particle<< std::endl;
												exit(0);
											}
											pthread_mutex_unlock(&global_mutex);

#endif
#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
												timer.tic(TIMING_DIFF_DIFF2);
#endif
											RFLOAT diff2;
											if ((iter == 1 && do_firstiter_cc) || do_always_cc)
											{
												// Do not calculate squared-differences, but signal product
												// Negative values because smaller is worse in this case
												diff2 = 0.;
												RFLOAT suma2 = 0.;
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
												{
													diff2 -= (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (*(Fimg_shift + n)).real;
												    diff2 -= (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (*(Fimg_shift + n)).imag;
													suma2 += norm(DIRECT_MULTIDIM_ELEM(Frefctf, n));
												}
												// Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
												diff2 /= sqrt(suma2) * exp_local_sqrtXi2[ipart];
											}
											else
											{
#ifdef DEBUG_CHECKSIZES
												if (ipart >= exp_highres_Xi2_imgs.size())
												{
													std::cerr<< "ipart= "<<ipart<<" exp_highres_Xi2_imgs.size()= "<< exp_highres_Xi2_imgs.size() <<std::endl;
													REPORT_ERROR("ipart >= exp_highres_Xi2_imgs.size()");
												}
#endif
												// Calculate the actual squared difference term of the Gaussian probability function
												// If current_size < mymodel.ori_size diff2 is initialised to the sum of
												// all |Xij|2 terms that lie between current_size and ori_size
												// Factor two because of factor 2 in division below, NOT because of 2-dimensionality of the complex plane!
												diff2 = exp_highres_Xi2_imgs[ipart] / 2.;
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
												{
													RFLOAT diff_real = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real - (*(Fimg_shift + n)).real;
													RFLOAT diff_imag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag - (*(Fimg_shift + n)).imag;
													diff2 += (diff_real * diff_real + diff_imag * diff_imag) * 0.5 * (*(Minvsigma2 + n));
												}
											}
#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
												timer.toc(TIMING_DIFF_DIFF2);
#endif

											// Store all diff2 in exp_Mweight
											long int ihidden_over = sampling.getPositionOversampledSamplingPoint(ihidden, exp_current_oversampling,
																											iover_rot, iover_trans);
//#define DEBUG_DIFF2_ISNAN
#ifdef DEBUG_DIFF2_ISNAN
											if (std::isnan(diff2))
											{
												pthread_mutex_lock(&global_mutex);
												std::cerr << " ipart= " << ipart << std::endl;
												std::cerr << " exp_iclass= " << exp_iclass << std::endl;
												std::cerr << " diff2= " << diff2 << " ipart= " << ipart << std::endl;
												std::cerr << " exp_highres_Xi2_imgs[ipart]= " << exp_highres_Xi2_imgs[ipart] << std::endl;
												std::cerr<< " exp_nr_oversampled_trans="<<exp_nr_oversampled_trans<<std::endl;
												std::cerr<< " exp_nr_oversampled_rot="<<exp_nr_oversampled_rot<<std::endl;
												std::cerr << " iover_rot= " << iover_rot << " iover_trans= " << iover_trans << " ihidden= " << ihidden << std::endl;
												std::cerr << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
												std::cerr << " ihidden_over= " << ihidden_over << " XSIZE(Mweight)= " << XSIZE(exp_Mweight) << std::endl;
												int group_id = mydata.getGroupId(part_id);
												std::cerr << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
												if (std::isnan(mymodel.scale_correction[group_id]))
												{
													for (int i=0; i < mymodel.scale_correction.size(); i++)
														std::cerr << " i= " << i << " mymodel.scale_correction[i]= " << mymodel.scale_correction[i] << std::endl;
												}
												std::cerr << " group_id= " << group_id << std::endl;
												Image<RFLOAT> It;
												std::cerr << "Frefctf shape= "; Frefctf.printShape(std::cerr);
												MultidimArray<Complex> Fish;
												Fish.resize(exp_local_Minvsigma2s[0]);
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fish)
												{
													DIRECT_MULTIDIM_ELEM(Fish, n) = *(Fimg_shift + n);
												}
												std::cerr << "Fimg_shift shape= "; (Fish).printShape(std::cerr);
												It()=exp_local_Fctfs[ipart];
												It.write("exp_local_Fctf.spi");
												std::cerr << "written exp_local_Fctf.spi" << std::endl;
												FourierTransformer transformer;
												Image<RFLOAT> tt;
												tt().resize(exp_current_image_size, exp_current_image_size);
												transformer.inverseFourierTransform(Fish, tt());
												CenterFFT(tt(),false);
												tt.write("Fimg_shift.spi");
												std::cerr << "written Fimg_shift.spi" << std::endl;
												FourierTransformer transformer2;
												tt().initZeros();
												transformer2.inverseFourierTransform(Frefctf, tt());
												CenterFFT(tt(),false);
												tt.write("Frefctf.spi");
												std::cerr << "written Frefctf.spi" << std::endl;
												FourierTransformer transformer3;
												tt().initZeros();
												transformer3.inverseFourierTransform(Fref, tt());
												CenterFFT(tt(),false);
												tt.write("Fref.spi");
												std::cerr << "written Fref.spi" << std::endl;
												std::cerr << " A= " << A << std::endl;
												std::cerr << "written Frefctf.spi" << std::endl;

												std::cerr << " exp_iclass= " << exp_iclass << std::endl;
												Fref.resize(exp_local_Minvsigma2s[0]);
												(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);
												transformer3.inverseFourierTransform(Fref, tt());
												CenterFFT(tt(),false);
												tt.write("Fref2.spi");
												std::cerr << "written Fref2.spi" << std::endl;
												Image<RFLOAT> Itt;
												Itt().resize(ZSIZE(mymodel.PPref[exp_iclass].data), YSIZE(mymodel.PPref[exp_iclass].data), XSIZE(mymodel.PPref[exp_iclass].data));
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Itt())
												{
													DIRECT_MULTIDIM_ELEM(Itt(), n) = abs(DIRECT_MULTIDIM_ELEM(mymodel.PPref[exp_iclass].data, n));
												}
												Itt.write("PPref_data.spi");
												REPORT_ERROR("diff2 is not a number");
												pthread_mutex_unlock(&global_mutex);
											}
#endif
//#define DEBUG_VERBOSE
#ifdef DEBUG_VERBOSE
											pthread_mutex_lock(&global_mutex);
											std::cout <<" name= "<< mydata.ori_particles[my_ori_particle].name << " rot= " << oversampled_rot[iover_rot] << " tilt= "<< oversampled_tilt[iover_rot] << " psi= " << oversampled_psi[iover_rot] << std::endl;
											std::cout <<" name= "<< mydata.ori_particles[my_ori_particle].name << " ihidden_over= " << ihidden_over << " diff2= " << diff2 << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
											pthread_mutex_unlock(&global_mutex);
#endif
#ifdef DEBUG_CHECKSIZES
											if (ihidden_over >= XSIZE(exp_Mweight) )
											{
												std::cerr<< " exp_nr_oversampled_trans="<<exp_nr_oversampled_trans<<std::endl;
												std::cerr<< " exp_nr_oversampled_rot="<<exp_nr_oversampled_rot<<std::endl;
												std::cerr << " iover_rot= " << iover_rot << " iover_trans= " << iover_trans << " ihidden= " << ihidden << std::endl;
												std::cerr << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
												std::cerr << " exp_itrans_min= " << exp_itrans_min <<" exp_nr_trans= " << exp_nr_trans << std::endl;
												std::cerr << " exp_itrans_max= " << exp_itrans_max << " iorientclass= " << iorientclass << " itrans= " << itrans << std::endl;
												std::cerr << " exp_nr_dir= " << exp_nr_dir << " exp_idir_min= " << exp_idir_min << " exp_idir_max= " << exp_idir_max << std::endl;
												std::cerr << " exp_nr_psi= " << exp_nr_psi << " exp_ipsi_min= " << exp_ipsi_min << " exp_ipsi_max= " << exp_ipsi_max << std::endl;
												std::cerr << " exp_iclass= " << exp_iclass << std::endl;
												std::cerr << " iorient= " << iorient << std::endl;
												std::cerr << " ihidden_over= " << ihidden_over << " XSIZE(Mweight)= " << XSIZE(exp_Mweight) << std::endl;
												REPORT_ERROR("ihidden_over >= XSIZE(Mweight)");
											}
#endif
											//std::cerr << " my_ori_particle= " << my_ori_particle<< " ipart= " << ipart << " ihidden_over= " << ihidden_over << " diff2= " << diff2 << std::endl;
											DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = diff2;
#ifdef DEBUG_CHECKSIZES
											if (ipart >= exp_min_diff2.size())
											{
												std::cerr<< "ipart= "<<ipart<<" exp_min_diff2.size()= "<< exp_min_diff2.size() <<std::endl;
												REPORT_ERROR("ipart >= exp_min_diff2.size() ");
											}
#endif
											// Keep track of minimum of all diff2, only for the last image in this series
											if (diff2 < exp_min_diff2[ipart])
												exp_min_diff2[ipart] = diff2;

										} // end loop iover_trans
									} // end if do_proceed translations
								} // end loop itrans
							} // end loop part_id
						}// end loop iover_rot
					} // end if do_proceed orientations
				} // end loop ipsi
			} // end loop idir
		} // end if mymodel.pdf_class[iclass] > 0.
	} // end loop iclass

#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
	{
		if (exp_ipass == 0) timer.toc(TIMING_ESP_DIFF1);
		else timer.toc(TIMING_ESP_DIFF2);
	}
#endif

}


void MlOptimiser::convertAllSquaredDifferencesToWeights(long int my_ori_particle, int exp_ipass,
		int exp_current_oversampling, int metadata_offset,
		int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
		int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
		MultidimArray<RFLOAT> &exp_Mweight, MultidimArray<bool> &exp_Mcoarse_significant,
		std::vector<RFLOAT> &exp_significant_weight, std::vector<RFLOAT> &exp_sum_weight,
		std::vector<Matrix1D<RFLOAT> > &exp_old_offset, std::vector<Matrix1D<RFLOAT> > &exp_prior,
		std::vector<RFLOAT> &exp_min_diff2,
		std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
		std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior)
{

#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
	{
		if (exp_ipass == 0) timer.tic(TIMING_ESP_WEIGHT1);
		else timer.tic(TIMING_ESP_WEIGHT2);
	}
#endif

	// Convert the squared differences into weights
	// Note there is only one weight for each part_id, because a whole series of images is treated as one particle

	long int exp_nr_dir = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
	long int exp_nr_psi = (do_skip_align || do_skip_rotate || do_only_sample_tilt) ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
	long int exp_nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings();
	long int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
	long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

	// Initialising...
	exp_sum_weight.clear();
	exp_sum_weight.resize(exp_nr_particles, 0.);

//#define DEBUG_CONVERTDIFF2W
#ifdef DEBUG_CONVERTDIFF2W
	RFLOAT max_weight = -1.;
	RFLOAT opt_psi, opt_xoff, opt_yoff;
	int opt_iover_rot, opt_iover_trans, opt_ipsi, opt_itrans;
	long int opt_ihidden, opt_ihidden_over;
#endif

	// loop over all particles inside this ori_particle
	for (long int ipart = 0; ipart < exp_nr_particles; ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		RFLOAT exp_thisparticle_sumweight = 0.;

		RFLOAT old_offset_z;
		RFLOAT old_offset_x = XX(exp_old_offset[ipart]);
		RFLOAT old_offset_y = YY(exp_old_offset[ipart]);
		if (mymodel.data_dim == 3)
			old_offset_z = ZZ(exp_old_offset[ipart]);

		if ((iter == 1 && do_firstiter_cc) || do_always_cc)
		{
			// Binarize the squared differences array to skip marginalisation
			RFLOAT mymindiff2 = 99.e10;
			long int myminidx = -1;
			// Find the smallest element in this row of exp_Mweight
			for (long int i = 0; i < XSIZE(exp_Mweight); i++)
			{

				RFLOAT cc = DIRECT_A2D_ELEM(exp_Mweight, ipart, i);
				// ignore non-determined cc
				if (cc == -999.)
					continue;

				// just search for the maximum
				if (cc < mymindiff2)
				{
					mymindiff2 = cc;
					myminidx = i;
				}
			}
			// Set all except for the best hidden variable to zero and the smallest element to 1
			for (long int i = 0; i < XSIZE(exp_Mweight); i++)
				DIRECT_A2D_ELEM(exp_Mweight, ipart, i)= 0.;

			DIRECT_A2D_ELEM(exp_Mweight, ipart, myminidx)= 1.;
			exp_thisparticle_sumweight += 1.;

		}
		else
		{
			// Extra normalization
			RFLOAT pdf_orientation_mean(0),pdf_offset_mean(0);
			unsigned long pdf_orientation_count(0), pdf_offset_count(0);
			for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
			{
				for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
					for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++)
					{
						if (do_skip_align || do_skip_rotate)
							pdf_orientation_mean += mymodel.pdf_class[exp_iclass];
						else if (mymodel.orientational_prior_mode == NOPRIOR)
							pdf_orientation_mean += DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
						else
							pdf_orientation_mean += exp_directions_prior[idir] * exp_psi_prior[ipsi];
						pdf_orientation_count++;
					}

				RFLOAT myprior_x, myprior_y, myprior_z;
				if (mymodel.ref_dim == 2)
				{
					myprior_x = XX(mymodel.prior_offset_class[exp_iclass]);
					myprior_y = YY(mymodel.prior_offset_class[exp_iclass]);
				}
				else
				{
					myprior_x = XX(exp_prior[ipart]);
					myprior_y = YY(exp_prior[ipart]);
					if (mymodel.data_dim == 3)
						myprior_z = ZZ(exp_prior[ipart]);
				}
				for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++)
				{
					RFLOAT offset_x = old_offset_x + sampling.translations_x[itrans];
					RFLOAT offset_y = old_offset_y + sampling.translations_y[itrans];
					RFLOAT tdiff2 = 0.;
					if ( (!do_helical_refine) || (ignore_helical_symmetry) || (mymodel.data_dim == 3) )
						tdiff2 += (offset_x - myprior_x) * (offset_x - myprior_x);
					tdiff2 += (offset_y - myprior_y) * (offset_y - myprior_y);
					if (mymodel.data_dim == 3)
					{
						RFLOAT offset_z = old_offset_z + sampling.translations_z[itrans];
						if ( (!do_helical_refine) || (ignore_helical_symmetry) )
							tdiff2 += (offset_z - myprior_z) * (offset_z - myprior_z);
					}
					// P(offset|sigma2_offset)
					// This is the probability of the offset, given the model offset and variance.
					RFLOAT pdf_offset;
					if (mymodel.sigma2_offset < 0.0001)
						pdf_offset_mean += ( tdiff2 > 0.) ? 0. : 1.;
					else
						pdf_offset_mean += exp ( tdiff2 / (-2. * mymodel.sigma2_offset) ) / ( 2. * PI * mymodel.sigma2_offset );
					pdf_offset_count ++;
				}
			}
			pdf_orientation_mean /= (RFLOAT) pdf_orientation_count;
			pdf_offset_mean /= (RFLOAT) pdf_offset_count;

			// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
			for (int exp_iclass = 0; exp_iclass < mymodel.nr_classes; exp_iclass++)
			{

				// Make PdfOffset calculation much faster...
				RFLOAT myprior_x, myprior_y, myprior_z;
				if (mymodel.ref_dim == 2)
				{
					myprior_x = XX(mymodel.prior_offset_class[exp_iclass]);
					myprior_y = YY(mymodel.prior_offset_class[exp_iclass]);
				}
				else
				{
					myprior_x = XX(exp_prior[ipart]);
					myprior_y = YY(exp_prior[ipart]);
					if (mymodel.data_dim == 3)
						myprior_z = ZZ(exp_prior[ipart]);
				}
				for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
				{
					for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++)
					{
						long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;
						RFLOAT pdf_orientation;

						// Get prior for this direction
						if (do_skip_align || do_skip_rotate)
						{
							pdf_orientation = mymodel.pdf_class[exp_iclass];
						}
						else if (mymodel.orientational_prior_mode == NOPRIOR)
						{
#ifdef DEBUG_CHECKSIZES
							if (idir >= XSIZE(mymodel.pdf_direction[exp_iclass]))
							{
								std::cerr<< "idir= "<<idir<<" XSIZE(mymodel.pdf_direction[exp_iclass])= "<< XSIZE(mymodel.pdf_direction[exp_iclass]) <<std::endl;
								REPORT_ERROR("idir >= mymodel.pdf_direction[exp_iclass].size()");
							}
#endif
							pdf_orientation = DIRECT_MULTIDIM_ELEM(mymodel.pdf_direction[exp_iclass], idir);
						}
						else
						{
							// P(orientation) = P(idir|dir_prior) * P(ipsi|psi_prior)
							// This is the probability of the orientation, given the gathered
							// statistics of all assigned orientations of the dataset, since we
							// are assigning a gaussian prior to all parameters.
							pdf_orientation = exp_directions_prior[idir] * exp_psi_prior[ipsi];
						}

						if (pdf_orientation_mean != 0.)
							pdf_orientation /= pdf_orientation_mean;

						// Loop over all translations
						long int ihidden = iorientclass * exp_nr_trans;
						for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++)
						{
							// May18,2015 - Shaoda & Sjors - Helical refinement (translational searches)
							// Calculate the vector length of myprior
							RFLOAT mypriors_len2 = myprior_x * myprior_x + myprior_y * myprior_y;
							if (mymodel.data_dim == 3)
								mypriors_len2 += myprior_z * myprior_z;
							// If it is doing helical refinement AND Cartesian vector myprior has a length > 0, transform the vector to its helical coordinates
							if ( (do_helical_refine) && (!ignore_helical_symmetry) && (mypriors_len2 > 0.00001) )
							{
								RFLOAT rot_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT);
								RFLOAT tilt_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT);
								RFLOAT psi_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI);
								transformCartesianAndHelicalCoords(myprior_x, myprior_y, myprior_z, myprior_x, myprior_y, myprior_z, rot_deg, tilt_deg, psi_deg, mymodel.data_dim, CART_TO_HELICAL_COORDS);
							}
							// (For helical refinement) Now offset, old_offset, sampling.translations and myprior are all in helical coordinates

							// To speed things up, only calculate pdf_offset at the coarse sampling.
							// That should not matter much, and that way one does not need to calculate all the OversampledTranslations
							RFLOAT offset_x = old_offset_x + sampling.translations_x[itrans];
							RFLOAT offset_y = old_offset_y + sampling.translations_y[itrans];
							RFLOAT tdiff2 = 0.;
							if ( (!do_helical_refine) || (ignore_helical_symmetry) || (mymodel.data_dim == 3) )
								tdiff2 += (offset_x - myprior_x) * (offset_x - myprior_x);
							tdiff2 += (offset_y - myprior_y) * (offset_y - myprior_y);
							if (mymodel.data_dim == 3)
							{
								RFLOAT offset_z = old_offset_z + sampling.translations_z[itrans];
								if ( (!do_helical_refine) || (ignore_helical_symmetry) )
									tdiff2 += (offset_z - myprior_z) * (offset_z - myprior_z);
							}
							// P(offset|sigma2_offset)
							// This is the probability of the offset, given the model offset and variance.
							RFLOAT pdf_offset;
							if (mymodel.sigma2_offset < 0.0001)
								pdf_offset = ( tdiff2 > 0.) ? 0. : 1.;
							else
								pdf_offset = exp ( tdiff2 / (-2. * mymodel.sigma2_offset) ) / ( 2. * PI * mymodel.sigma2_offset );

							if (pdf_offset_mean != 0.)
								pdf_offset /= pdf_offset_mean;

							/*
							// TMP DEBUGGING
							if (mymodel.orientational_prior_mode != NOPRIOR && (pdf_offset==0. || pdf_orientation==0.))
							{
								pthread_mutex_lock(&global_mutex);
								std::cerr << " pdf_offset= " << pdf_offset << " pdf_orientation= " << pdf_orientation << std::endl;
								std::cerr << " ipart= " << ipart << " part_id= " << part_id << std::endl;
								std::cerr << " iorient= " << iorient << " idir= " << idir << " ipsi= " << ipsi << std::endl;
								//std::cerr << " exp_nr_psi= " << exp_nr_psi << " exp_nr_dir= " << exp_nr_dir << " exp_nr_trans= " << exp_nr_trans << std::endl;
								for (long int i = 0; i < exp_directions_prior.size(); i++)
									std::cerr << " exp_directions_prior["<<i<<"]= " << exp_directions_prior[i] << std::endl;
								for (long int i = 0; i < exp_psi_prior.size(); i++)
									std::cerr << " exp_psi_prior["<<i<<"]= " << exp_psi_prior[i] << std::endl;
								REPORT_ERROR("ERROR! pdf_offset==0.|| pdf_orientation==0.");
								pthread_mutex_unlock(&global_mutex);
							}
							if (exp_nr_oversampled_rot == 0)
								REPORT_ERROR("exp_nr_oversampled_rot == 0");
							if (exp_nr_oversampled_trans == 0)
								REPORT_ERROR("exp_nr_oversampled_trans == 0");
							*/
#ifdef TIMING
							// Only time one thread, as I also only time one MPI process
							if (my_ori_particle == exp_my_first_ori_particle)
								timer.tic(TIMING_WEIGHT_EXP);
#endif
							// Now first loop over iover_rot, because that is the order in exp_Mweight as well
							long int ihidden_over = ihidden * exp_nr_oversampled_rot * exp_nr_oversampled_trans;
							for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
							{
								// Then loop over iover_trans
								for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++, ihidden_over++)
								{
#ifdef DEBUG_CHECKSIZES
									if (ihidden_over >= XSIZE(exp_Mweight))
									{
										std::cerr<< "ihidden_over= "<<ihidden_over<<" XSIZE(Mweight)= "<< XSIZE(exp_Mweight) <<std::endl;
										REPORT_ERROR("ihidden_over >= XSIZE(exp_Mweight)");
									}
#endif
									// Only exponentiate for determined values of exp_Mweight
									// (this is always true in the first pass, but not so in the second pass)
									// Only deal with this sampling point if its weight was significant
#ifdef DEBUG_CHECKSIZES
									if (ipart >= YSIZE(exp_Mweight))
									{
										std::cerr << " YSIZE(exp_Mweight)= "<< YSIZE(exp_Mweight) <<std::endl;
										std::cerr << " ipart= " << ipart << std::endl;
										REPORT_ERROR("ipart >= YSIZE(exp_Mweight)");
									}
#endif
									if (DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) < 0.)
									{
										DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = 0.;
									}
									else
									{
										// Set the weight base to the probability of the parameters given the prior
										RFLOAT weight = pdf_orientation * pdf_offset;
										RFLOAT diff2 = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) - exp_min_diff2[ipart];
										// next line because of numerical precision of exp-function
#ifdef RELION_SINGLE_PRECISION
										if (diff2 > 88.)
											weight = 0.;
#else
										if (diff2 > 700.)
											weight = 0.;
#endif
										// TODO: use tabulated exp function?
										else weight *= exp(-diff2);
//#define DEBUG_PSIANGLE_PDISTRIBUTION
#ifdef DEBUG_PSIANGLE_PDISTRIBUTION
										std::cout << ipsi*360./sampling.NrPsiSamplings() << " "<< weight << std::endl;
#endif
										// Store the weight
										DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) = weight;
#ifdef DEBUG_CHECKSIZES
										if (std::isnan(weight))
										{
											pthread_mutex_lock(&global_mutex);
											std::cerr<< "weight= "<<weight<<" is not a number! " <<std::endl;
											std::cerr << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
											std::cerr << " ipart= " << ipart << std::endl;
											std::cerr << " part_id= " << part_id << std::endl;
											std::cerr << " DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over)= " << DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over) << std::endl;
											REPORT_ERROR("weight is not a number");
											pthread_mutex_unlock(&global_mutex);
										}
#endif
										// Keep track of sum and maximum of all weights for this particle
										// Later add all to exp_thisparticle_sumweight, but inside this loop sum to local thisthread_sumweight first
										exp_thisparticle_sumweight += weight;
									} // end if/else exp_Mweight < 0.
								} // end loop iover_trans
							}// end loop iover_rot
#ifdef TIMING
							// Only time one thread, as I also only time one MPI process
							if (my_ori_particle == exp_my_first_ori_particle)
								timer.toc(TIMING_WEIGHT_EXP);
#endif
						} // end loop itrans
					} // end loop ipsi
				} // end loop idir
			} // end loop exp_iclass
		} // end if iter==1

		//Store parameters for this particle
		exp_sum_weight[ipart] = exp_thisparticle_sumweight;

		// Check the sum of weights is not zero
// On a Mac, the isnan function does not compile. Just uncomment the define statement, as this is merely a debugging statement
//#define MAC_OSX
#ifndef MAC_OSX
		if (exp_thisparticle_sumweight == 0. || std::isnan(exp_thisparticle_sumweight))
		{
			std::cerr << " exp_thisparticle_sumweight= " << exp_thisparticle_sumweight << std::endl;
			Image<RFLOAT> It;
			It() = exp_Mweight;
			It.write("Mweight.spi");
			//It() = DEBUGGING_COPY_exp_Mweight;
			//It.write("Mweight_copy.spi");
			It().resize(exp_Mcoarse_significant);
			if (MULTIDIM_SIZE(It()) > 0)
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(It())
				{
					if (DIRECT_MULTIDIM_ELEM(exp_Mcoarse_significant, n))
						DIRECT_MULTIDIM_ELEM(It(), n) = 1.;
					else
						DIRECT_MULTIDIM_ELEM(It(), n) = 0.;
				}
				It.write("Mcoarse_significant.spi");
			}
			std::cerr << " part_id= " << part_id << std::endl;
			std::cerr << " ipart= " << ipart << std::endl;
			/*
			MultidimArray<Complex> Faux;
			FourierTransformer transformer;
			windowFourierTransform(exp_Fimgs[ipart], Faux, mymodel.ori_size);
			It().resize(mymodel.ori_size, mymodel.ori_size);
			transformer.inverseFourierTransform(Faux, It());
			CenterFFT(It(), false);
			It.write("exp_Fimgs.spi");
			std::cerr << "written exp_Fimgs.spi " << std::endl;
		    */
			int group_id = mydata.getGroupId(part_id);
			std::cerr << " group_id= " << group_id << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
			std::cerr << " exp_ipass= " << exp_ipass << std::endl;
			std::cerr << " sampling.NrDirections(0, true)= " << sampling.NrDirections()
					<< " sampling.NrDirections(0, false)= " << sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior) << std::endl;
			std::cerr << " sampling.NrPsiSamplings(0, true)= " << sampling.NrPsiSamplings()
					<< " sampling.NrPsiSamplings(0, false)= " << sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior) << std::endl;
			std::cerr << " mymodel.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
			//std::cerr << " wsum_model.sigma2_noise[group_id]= " << wsum_model.sigma2_noise[group_id] << std::endl;
			if (do_norm_correction)
			{
				std::cerr << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
				std::cerr << " wsum_model.avg_norm_correction= " << wsum_model.avg_norm_correction << std::endl;
			}

			std::cerr << "written out Mweight.spi" << std::endl;
			std::cerr << " exp_thisparticle_sumweight= " << exp_thisparticle_sumweight << std::endl;
			std::cerr << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
			REPORT_ERROR("ERROR!!! zero sum of weights....");
		}
#endif

	} // end loop ipart

	// Initialise exp_Mcoarse_significant
	if (exp_ipass==0)
		exp_Mcoarse_significant.resize(exp_nr_particles, XSIZE(exp_Mweight));

	// Now, for each particle,  find the exp_significant_weight that encompasses adaptive_fraction of exp_sum_weight
	exp_significant_weight.clear();
	exp_significant_weight.resize(exp_nr_particles, 0.);
	for (long int ipart = 0; ipart < exp_nr_particles; ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];

#ifdef TIMING
		if (my_ori_particle == exp_my_first_ori_particle)
			timer.tic(TIMING_WEIGHT_SORT);
#endif
		MultidimArray<RFLOAT> sorted_weight;
		// Get the relevant row for this particle
		exp_Mweight.getRow(ipart, sorted_weight);

		// Only select non-zero probabilities to speed up sorting
		long int np = 0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
		{
			if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) > 0.)
			{
				DIRECT_MULTIDIM_ELEM(sorted_weight, np) = DIRECT_MULTIDIM_ELEM(sorted_weight, n);
				np++;
			}
		}
		sorted_weight.resize(np);

		// Sort from low to high values
		sorted_weight.sort();

#ifdef TIMING
		if (my_ori_particle == exp_my_first_ori_particle)
			timer.toc(TIMING_WEIGHT_SORT);
#endif
		RFLOAT frac_weight = 0.;
		RFLOAT my_significant_weight;
		long int my_nr_significant_coarse_samples = 0;
		for (long int i = XSIZE(sorted_weight) - 1; i >= 0; i--)
		{
			if (maximum_significants > 0 )
			{
				if(my_nr_significant_coarse_samples < maximum_significants)
				{
					if (exp_ipass==0)
						my_nr_significant_coarse_samples++;
					my_significant_weight = DIRECT_A1D_ELEM(sorted_weight, i);
				}
			}
			else
			{
				if (exp_ipass==0)
					my_nr_significant_coarse_samples++;
				my_significant_weight = DIRECT_A1D_ELEM(sorted_weight, i);
			}
			frac_weight += DIRECT_A1D_ELEM(sorted_weight, i);
			if (frac_weight > adaptive_fraction * exp_sum_weight[ipart])
				break;
		}


#ifdef DEBUG_SORT
		// Check sorted array is really sorted
		RFLOAT prev = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
		{
			if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) < prev)
			{
				Image<RFLOAT> It;
				It()=sorted_weight;
				It() *= 10000;
				It.write("sorted_weight.spi");
				std::cerr << "written sorted_weight.spi" << std::endl;
				REPORT_ERROR("Error in sorting!");
			}
			prev=DIRECT_MULTIDIM_ELEM(sorted_weight, n);
		}
#endif

		if (exp_ipass==0 && my_nr_significant_coarse_samples == 0)
		{
			std::cerr << " ipart= " << ipart << " adaptive_fraction= " << adaptive_fraction << std::endl;
			std::cerr << " frac-weight= " << frac_weight << std::endl;
			std::cerr << " exp_sum_weight[ipart]= " << exp_sum_weight[ipart] << std::endl;
			Image<RFLOAT> It;
			std::cerr << " XSIZE(exp_Mweight)= " << XSIZE(exp_Mweight) << std::endl;
			It()=exp_Mweight;
			It() *= 10000;
			It.write("Mweight2.spi");
			std::cerr << "written Mweight2.spi" << std::endl;
			std::cerr << " np= " << np << std::endl;
			It()=sorted_weight;
			It() *= 10000;
			std::cerr << " XSIZE(sorted_weight)= " << XSIZE(sorted_weight) << std::endl;
			if (XSIZE(sorted_weight) > 0)
			{
				It.write("sorted_weight.spi");
				std::cerr << "written sorted_weight.spi" << std::endl;
			}
			REPORT_ERROR("my_nr_significant_coarse_samples == 0");
		}

		if (exp_ipass==0)
		{
			// Store nr_significant_coarse_samples for this particle
			DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NR_SIGN) = (RFLOAT)my_nr_significant_coarse_samples;

			// Keep track of which coarse samplings were significant were significant for this particle
			for (int ihidden = 0; ihidden < XSIZE(exp_Mcoarse_significant); ihidden++)
			{
				if (DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden) >= my_significant_weight)
					DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden) = true;
				else
					DIRECT_A2D_ELEM(exp_Mcoarse_significant, ipart, ihidden) = false;
			}

		}
		exp_significant_weight[ipart] = my_significant_weight;
#ifdef DEBUG_OVERSAMPLING
		std::cerr << " sum_weight[ipart]= " << exp_sum_weight[ipart] << " my_significant_weight= " << my_significant_weight << std::endl;
		std::cerr << " my_nr_significant_coarse_samples= " << my_nr_significant_coarse_samples << std::endl;
		std::cerr << " ipass= " << exp_ipass << " Pmax="<<DIRECT_A1D_ELEM(sorted_weight,XSIZE(sorted_weight) - 1)/frac_weight
				<<" nr_sign_sam= "<<nr_significant_samples<<" sign w= "<<exp_significant_weight<< "sum_weight= "<<exp_sum_weight<<std::endl;
#endif

	} // end loop ipart


#ifdef DEBUG_CONVERTDIFF2W
	//Image<RFLOAT> tt;
	//tt()=sorted_weight;
	//tt.write("sorted_weight.spi");
	//std::cerr << "written sorted_weight.spi" << std::endl;
	std::cerr << " ipass= " << exp_ipass << " part_id= " << part_id << std::endl;
	std::cerr << " diff2w: opt_xoff= " << opt_xoff << " opt_yoff= " << opt_yoff << " opt_psi= " << opt_psi << std::endl;
	std::cerr << " diff2w: opt_iover_rot= " << opt_iover_rot << " opt_iover_trans= " << opt_iover_trans << " opt_ipsi= " << opt_ipsi << std::endl;
	std::cerr << " diff2w: opt_itrans= " << opt_itrans << " opt_ihidden= " << opt_ihidden << " opt_ihidden_over= " << opt_ihidden_over << std::endl;
	std::cerr << "significant_weight= " << exp_significant_weight << " max_weight= " << max_weight << std::endl;
	std::cerr << "nr_significant_coarse_samples= " << nr_significant_coarse_samples <<std::endl;
	debug2 = (RFLOAT)opt_ihidden_over;
#endif

#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
	{
		if (exp_ipass == 0) timer.toc(TIMING_ESP_WEIGHT1);
		else timer.toc(TIMING_ESP_WEIGHT2);
	}
#endif

}

void MlOptimiser::storeWeightedSums(long int my_ori_particle, int ibody, int exp_current_image_size,
		int exp_current_oversampling, int metadata_offset,
		int exp_idir_min, int exp_idir_max, int exp_ipsi_min, int exp_ipsi_max,
		int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
		std::vector<RFLOAT> &exp_min_diff2,
		std::vector<RFLOAT> &exp_highres_Xi2_imgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs,
		std::vector<MultidimArray<Complex > > &exp_Fimgs_nomask,
		std::vector<MultidimArray<RFLOAT> > &exp_Fctfs,
		std::vector<MultidimArray<RFLOAT> > &exp_power_imgs,
		std::vector<Matrix1D<RFLOAT> > &exp_old_offset,
		std::vector<Matrix1D<RFLOAT> > &exp_prior,
		MultidimArray<RFLOAT> &exp_Mweight,
		MultidimArray<bool> &exp_Mcoarse_significant,
		std::vector<RFLOAT> &exp_significant_weight,
		std::vector<RFLOAT> &exp_sum_weight,
		std::vector<RFLOAT> &exp_max_weight,
		std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
		std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted,
		std::vector<MultidimArray<Complex > > &exp_local_Fimgs_shifted_nomask,
		std::vector<MultidimArray<RFLOAT> > &exp_local_Minvsigma2s,
		std::vector<MultidimArray<RFLOAT> > &exp_local_Fctfs,
		std::vector<RFLOAT> &exp_local_sqrtXi2)
{

#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
		timer.tic(TIMING_ESP_WSUM);
#endif

	int exp_nr_particles = mydata.ori_particles[my_ori_particle].particles_id.size();
	long int exp_nr_dir = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
	long int exp_nr_psi = (do_skip_align || do_skip_rotate || do_only_sample_tilt) ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
	long int exp_nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings();
	long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
	long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

	// Re-do below because now also want unmasked images AND if (stricht_highres_exp >0.) then may need to resize
	precalculateShiftedImagesCtfsAndInvSigma2s(true, my_ori_particle, exp_current_image_size, exp_current_oversampling, metadata_offset,
			exp_itrans_min, exp_itrans_max, exp_Fimgs, exp_Fimgs_nomask, exp_Fctfs, exp_local_Fimgs_shifted, exp_local_Fimgs_shifted_nomask,
			exp_local_Fctfs, exp_local_sqrtXi2, exp_local_Minvsigma2s);

	// In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the exp_local_Minvsigma2s was omitted.
	// Set those back here
	for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		int group_id = mydata.getGroupId(part_id);
		DIRECT_MULTIDIM_ELEM(exp_local_Minvsigma2s[ipart], 0) = 1. / (sigma2_fudge * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], 0));
	}

	// Initialise the maximum of all weights to a negative value
	exp_max_weight.clear();
	exp_max_weight.resize(exp_nr_particles, -1.);

	// For norm_correction and scale_correction of all particles of this ori_particle
	std::vector<RFLOAT> exp_wsum_norm_correction;
	std::vector<MultidimArray<RFLOAT> > exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA;
	std::vector<MultidimArray<RFLOAT> > thr_wsum_signal_product_spectra, thr_wsum_reference_power_spectra;
	exp_wsum_norm_correction.resize(exp_nr_particles, 0.);

	// For scale_correction
	if (do_scale_correction)
	{
		MultidimArray<RFLOAT> aux;
		aux.initZeros(mymodel.ori_size/2 + 1);
		exp_wsum_scale_correction_XA.resize(exp_nr_particles, aux);
		exp_wsum_scale_correction_AA.resize(exp_nr_particles, aux);
		thr_wsum_signal_product_spectra.resize(mymodel.nr_groups, aux);
		thr_wsum_reference_power_spectra.resize(mymodel.nr_groups, aux);
	}

	std::vector< RFLOAT> oversampled_rot, oversampled_tilt, oversampled_psi;
	std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
	Matrix2D<RFLOAT> A, Abody, Aori;
	MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_otfshift, Fimg_otfshift_nomask, Fimg_store_sgd;
	MultidimArray<RFLOAT> Minvsigma2, Mctf, Fweight;
	RFLOAT rot, tilt, psi;
	bool have_warned_small_scale = false;
	// Initialising... exp_Fimgs[0] has mymodel.current_size (not coarse_size!)
	Fref.resize(exp_Fimgs[0]);
	Frefctf.resize(exp_Fimgs[0]);
	Fweight.resize(exp_Fimgs[0]);
	Fimg.resize(exp_Fimgs[0]);
	// Initialise Mctf to all-1 for if !do_ctf_corection
	Mctf.resize(exp_Fimgs[0]);
	Mctf.initConstant(1.);
	// Initialise Minvsigma2 to all-1 for if !do_map
	Minvsigma2.resize(exp_Fimgs[0]);
	Minvsigma2.initConstant(1.);
	if (do_shifts_onthefly)
	{
		Fimg_otfshift.resize(Frefctf);
		Fimg_otfshift_nomask.resize(Frefctf);
	}
	if (do_sgd)
	{
		Fimg_store_sgd.resize(Frefctf);
	}


	if (mymodel.nr_bodies > 1)
	{
		// ipart=0 because in multi-body refinement we do not do movie frames!
		RFLOAT rot_ori = DIRECT_A2D_ELEM(exp_metadata, metadata_offset, METADATA_ROT);
		RFLOAT tilt_ori = DIRECT_A2D_ELEM(exp_metadata, metadata_offset, METADATA_TILT);
		RFLOAT psi_ori = DIRECT_A2D_ELEM(exp_metadata, metadata_offset, METADATA_PSI);
		Euler_angles2matrix(rot_ori, tilt_ori, psi_ori, Aori, false);
	}

	// Make local copies of weighted sums (except BPrefs, which are too big)
	// so that there are not too many mutex locks below
	std::vector<MultidimArray<RFLOAT> > thr_wsum_sigma2_noise, thr_wsum_pdf_direction;
	std::vector<RFLOAT> thr_wsum_norm_correction, thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class;
	RFLOAT thr_wsum_sigma2_offset;
	MultidimArray<RFLOAT> thr_metadata, zeroArray;
	// Wsum_sigma_noise2 is a 1D-spectrum for each group
	zeroArray.initZeros(mymodel.ori_size/2 + 1);
	thr_wsum_sigma2_noise.resize(mymodel.nr_groups, zeroArray);
	// wsum_pdf_direction is a 1D-array (of length sampling.NrDirections()) for each class
	zeroArray.initZeros(sampling.NrDirections());
	thr_wsum_pdf_direction.resize(mymodel.nr_classes, zeroArray);
	// sumw_group is a RFLOAT for each group
	thr_sumw_group.resize(mymodel.nr_groups, 0.);
	// wsum_pdf_class is a RFLOAT for each class
	thr_wsum_pdf_class.resize(mymodel.nr_classes, 0.);
	if (mymodel.ref_dim == 2)
	{
		thr_wsum_prior_offsetx_class.resize(mymodel.nr_classes, 0.);
		thr_wsum_prior_offsety_class.resize(mymodel.nr_classes, 0.);
	}
	// wsum_sigma2_offset is just a RFLOAT
	thr_wsum_sigma2_offset = 0.;

	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++)
	{
		for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
		{
			for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++)
			{
				long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;

				// Only proceed if any of the particles had any significant coarsely sampled translation
				if (isSignificantAnyParticleAnyTranslation(iorientclass, exp_itrans_min, exp_itrans_max, exp_Mcoarse_significant))
				{
					// Now get the oversampled (rot, tilt, psi) triplets
					// This will be only the original (rot,tilt,psi) triplet if (adaptive_oversampling==0)
					sampling.getOrientations(idir, ipsi, adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
							exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior);
					// Loop over all oversampled orientations (only a single one in the first pass)
					for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++)
					{
						rot = oversampled_rot[iover_rot];
						tilt = oversampled_tilt[iover_rot];
						psi = oversampled_psi[iover_rot];
						// Get the Euler matrix
						Euler_angles2matrix(rot, tilt, psi, A, false);

						// For multi-body refinements, A are only 'residual' orientations, Abody is the complete Euler matrix
						if (mymodel.nr_bodies > 1)
							Abody = A * Aori;

#ifdef TIMING
						// Only time one thread, as I also only time one MPI process
						if (my_ori_particle == exp_my_first_ori_particle)
							timer.tic(TIMING_WSUM_PROJ);
#endif
						// Project the reference map (into Fref)
						if (!do_skip_maximization)
						{
							if (mymodel.nr_bodies > 1)
							{
								(mymodel.PPref[ibody]).get2DFourierTransform(Fref, Abody, IS_NOT_INV);
								// 19may2015: inside PPref, each body is centered at its own COM
								// Put it back to where it was!
								Matrix1D<RFLOAT> comp(3);
								comp = Abody * (mymodel.com_bodies[ibody]);
								shiftImageInFourierTransform(Fref, Fref, (RFLOAT)mymodel.ori_size, XX(comp), YY(comp), ZZ(comp));
							}
							else
								(mymodel.PPref[exp_iclass]).get2DFourierTransform(Fref, A, IS_NOT_INV);
						}
#ifdef TIMING
						// Only time one thread, as I also only time one MPI process
						if (my_ori_particle == exp_my_first_ori_particle)
							timer.toc(TIMING_WSUM_PROJ);
#endif
						// Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
						// Then outside this loop do the actual backprojection
						Fimg.initZeros();
						Fweight.initZeros();
						/// Now that reference projection has been made loop over all particles inside this ori_particle
						for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
						{
							// This is an attempt to speed up illogically slow updates of wsum_sigma2_offset....
							// It seems to make a big difference!
							RFLOAT myprior_x, myprior_y, myprior_z, old_offset_z;
							RFLOAT old_offset_x = XX(exp_old_offset[ipart]);
							RFLOAT old_offset_y = YY(exp_old_offset[ipart]);
							if (mymodel.ref_dim == 2)
							{
								myprior_x = XX(mymodel.prior_offset_class[exp_iclass]);
								myprior_y = YY(mymodel.prior_offset_class[exp_iclass]);
							}
							else
							{
								myprior_x = XX(exp_prior[ipart]);
								myprior_y = YY(exp_prior[ipart]);
								if (mymodel.data_dim == 3)
								{
									myprior_z = ZZ(exp_prior[ipart]);
									old_offset_z = ZZ(exp_old_offset[ipart]);
								}
							}

							long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
							int group_id = mydata.getGroupId(part_id);
#ifdef DEBUG_CHECKSIZES
							if (group_id >= mymodel.nr_groups)
							{
								std::cerr<< "group_id= "<<group_id<<" ml_model.nr_groups= "<< mymodel.nr_groups <<std::endl;
								REPORT_ERROR("group_id >= ml_model.nr_groups");
							}
#endif
							if (!do_skip_maximization)
							{
								if (do_map)
									Minvsigma2 = exp_local_Minvsigma2s[ipart];
								// else Minvsigma2 was initialised to ones
								// Apply CTF to reference projection
								if (do_ctf_correction)
								{
									Mctf = exp_local_Fctfs[ipart];
									if (refs_are_ctf_corrected)
									{
										if (ctf_premultiplied)
										{
											FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
											{
												DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(Mctf, n) * DIRECT_MULTIDIM_ELEM(Mctf, n);
											}
										}
										else
										{
											FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
											{
												DIRECT_MULTIDIM_ELEM(Frefctf, n) = DIRECT_MULTIDIM_ELEM(Fref, n) * DIRECT_MULTIDIM_ELEM(Mctf, n);
											}
										}
									}
									else
									{
										Frefctf = Fref;
									}
								}
								else
								{
									// initialise because there are multiple particles and Mctf gets selfMultiplied for scale_correction
									Mctf.initConstant(1.);
									Frefctf = Fref;
								}
								if (do_scale_correction)
								{
									RFLOAT myscale = mymodel.scale_correction[group_id];
									if (myscale > 10000.)
									{
										std::cerr << " rlnMicrographScaleCorrection= " << myscale << " group= " << group_id + 1 << std::endl;
										REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
									}
									else if (myscale < 0.001)
									{
										if (!have_warned_small_scale)
										{
											std::cout << " WARNING: ignoring group " << group_id + 1 << " with very small or negative scale (" << myscale <<
													"); Use larger groups for more stable scale estimates." << std::endl;
											have_warned_small_scale = true;
										}
										myscale = 0.001;
									}
									FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
									{
										DIRECT_MULTIDIM_ELEM(Frefctf, n) *= myscale;
									}
									// For CTF-terms in BP
									Mctf *= myscale;
								}
							} // end if !do_skip_maximization

							long int ihidden = iorientclass * exp_nr_trans;
							for (long int itrans = exp_itrans_min, iitrans = 0; itrans <= exp_itrans_max; itrans++, ihidden++)
							{
								// Jun01,2015 - Shaoda & Sjors, Helical refinement
								sampling.getTranslations(itrans, exp_current_oversampling, oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
										(do_helical_refine) && (!ignore_helical_symmetry), helical_rise_initial / mymodel.pixel_size, helical_twist_initial);
								for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++, iitrans++)
								{
#ifdef DEBUG_CHECKSIZES
									if (iover_trans >= oversampled_translations_x.size())
									{
										std::cerr<< "iover_trans= "<<iover_trans<<" oversampled_translations_x.size()= "<< oversampled_translations_x.size() <<std::endl;
										REPORT_ERROR("iover_trans >= oversampled_translations_x.size()");
									}
#endif
									// Only deal with this sampling point if its weight was significant
									long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
											iover_rot * exp_nr_oversampled_trans + iover_trans;
#ifdef DEBUG_CHECKSIZES
									if (ihidden_over >= XSIZE(exp_Mweight))
									{
										std::cerr<< "ihidden_over= "<<ihidden_over<<" XSIZE(exp_Mweight)= "<< XSIZE(exp_Mweight) <<std::endl;
										REPORT_ERROR("ihidden_over >= XSIZE(exp_Mweight)");
									}
									if (ipart >= exp_significant_weight.size())
									{
										std::cerr<< "ipart= "<<ipart<<" exp_significant_weight.size()= "<< exp_significant_weight.size() <<std::endl;
										REPORT_ERROR("ipart >= significant_weight.size()");
									}
									if (ipart >= exp_max_weight.size())
									{
										std::cerr<< "ipart= "<<ipart<<" exp_max_weight.size()= "<< exp_max_weight.size() <<std::endl;
										REPORT_ERROR("ipart >= exp_max_weight.size()");
									}
									if (ipart >= exp_sum_weight.size())
									{
										std::cerr<< "ipart= "<<ipart<<" exp_max_weight.size()= "<< exp_sum_weight.size() <<std::endl;
										REPORT_ERROR("ipart >= exp_sum_weight.size()");
									}
#endif
									RFLOAT weight = DIRECT_A2D_ELEM(exp_Mweight, ipart, ihidden_over);
									// Only sum weights for non-zero weights
									if (weight >= exp_significant_weight[ipart])
									{
										// Normalise the weight (do this after the comparison with exp_significant_weight!)
										weight /= exp_sum_weight[ipart];

										if (!do_skip_maximization)
										{

#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
												timer.tic(TIMING_WSUM_GETSHIFT);
#endif

											/// Now get the shifted image
											// Use a pointer to avoid copying the entire array again in this highly expensive loop
											Complex *Fimg_shift, *Fimg_shift_nomask;
											if (!do_shifts_onthefly)
											{
												long int ishift = ipart * exp_nr_oversampled_trans * exp_nr_trans + iitrans;
												Fimg_shift = exp_local_Fimgs_shifted[ishift].data;
												Fimg_shift_nomask = exp_local_Fimgs_shifted_nomask[ishift].data;
											}
											else
											{
												// Feb01,2017 - Shaoda, on-the-fly shifts in helical reconstuctions (2D and 3D)
												if ( (do_helical_refine) && (!ignore_helical_symmetry) )
												{
													RFLOAT xshift = 0., yshift = 0., zshift = 0.;

													xshift = oversampled_translations_x[iover_trans];
													yshift = oversampled_translations_y[iover_trans];
													if (mymodel.data_dim == 3)
														zshift = oversampled_translations_z[iover_trans];

													RFLOAT rot_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT);
													RFLOAT tilt_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT);
													RFLOAT psi_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI);
													transformCartesianAndHelicalCoords(
																xshift, yshift, zshift,
																xshift, yshift, zshift,
																rot_deg, tilt_deg, psi_deg,
																mymodel.data_dim,
																HELICAL_TO_CART_COORDS);

													// Fimg_shift
													shiftImageInFourierTransformWithTabSincos(
															exp_local_Fimgs_shifted[ipart],
															Fimg_otfshift,
															(RFLOAT)mymodel.ori_size,
															mymodel.current_size,
															tab_sin, tab_cos,
															xshift, yshift, zshift);
													// Fimg_shift_nomask
													shiftImageInFourierTransformWithTabSincos(
															exp_local_Fimgs_shifted_nomask[ipart],
															Fimg_otfshift_nomask,
															(RFLOAT)mymodel.ori_size,
															mymodel.current_size,
															tab_sin, tab_cos,
															xshift, yshift, zshift);
												}
												else
												{
													Complex* myAB;
													myAB = (adaptive_oversampling == 0 ) ? global_fftshifts_ab_current[iitrans].data : global_fftshifts_ab2_current[iitrans].data;
													FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(exp_local_Fimgs_shifted[ipart])
													{
														RFLOAT a = (*(myAB + n)).real;
														RFLOAT b = (*(myAB + n)).imag;
														// Fimg_shift
														RFLOAT real = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real
																- b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag;
														RFLOAT imag = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).imag
																+ b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted[ipart], n)).real;
														DIRECT_MULTIDIM_ELEM(Fimg_otfshift, n) = Complex(real, imag);
														// Fimg_shift_nomask
														real = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).real
																- b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).imag;
														imag = a * (DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).imag
																+ b *(DIRECT_MULTIDIM_ELEM(exp_local_Fimgs_shifted_nomask[ipart], n)).real;
														DIRECT_MULTIDIM_ELEM(Fimg_otfshift_nomask, n) = Complex(real, imag);
													}
												}
												Fimg_shift = Fimg_otfshift.data;
												Fimg_shift_nomask = Fimg_otfshift_nomask.data;
											}
#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
											{
												timer.toc(TIMING_WSUM_GETSHIFT);
												timer.tic(TIMING_WSUM_DIFF2);
											}
#endif

											// Store weighted sum of squared differences for sigma2_noise estimation
											// Suggestion Robert Sinkovitz: merge difference and scale steps to make better use of cache
											FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
											{
												int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
												if (ires > -1)
												{
													// Use FT of masked image for noise estimation!
													RFLOAT diff_real = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real - (*(Fimg_shift + n)).real;
													RFLOAT diff_imag = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag - (*(Fimg_shift + n)).imag;
													RFLOAT wdiff2 = weight * (diff_real*diff_real + diff_imag*diff_imag);
													// group-wise sigma2_noise
													DIRECT_MULTIDIM_ELEM(thr_wsum_sigma2_noise[group_id], ires) += wdiff2;
													// For norm_correction
													exp_wsum_norm_correction[ipart] += wdiff2;
													if (do_scale_correction  && DIRECT_A1D_ELEM(mymodel.data_vs_prior_class[exp_iclass], ires) > 3.)
													{
														RFLOAT sumXA, sumA2;
														sumXA = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (*(Fimg_shift + n)).real;
														sumXA += (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (*(Fimg_shift + n)).imag;
														DIRECT_A1D_ELEM(exp_wsum_scale_correction_XA[ipart], ires) += weight * sumXA;
														sumA2 = (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real * (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
														sumA2 += (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag * (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
														DIRECT_A1D_ELEM(exp_wsum_scale_correction_AA[ipart], ires) += weight * sumA2;
													}
												}
											}
#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
											{
												timer.toc(TIMING_WSUM_DIFF2);
												timer.tic(TIMING_WSUM_LOCALSUMS);
											}
#endif

											// Store sum of weights for this group
											thr_sumw_group[group_id] += weight;
											// Store weights for this class and orientation
											thr_wsum_pdf_class[exp_iclass] += weight;

											// The following goes MUCH faster than the original lines below....
											if (mymodel.ref_dim == 2)
											{
												thr_wsum_prior_offsetx_class[exp_iclass] += weight * (old_offset_x + oversampled_translations_x[iover_trans]);
												thr_wsum_prior_offsety_class[exp_iclass] += weight * (old_offset_y + oversampled_translations_y[iover_trans]);
											}
											// May18,2015 - Shaoda & Sjors, Helical refinement (translational searches)
											// Calculate the vector length of myprior
											RFLOAT mypriors_len2 = myprior_x * myprior_x + myprior_y * myprior_y;
											if (mymodel.data_dim == 3)
												mypriors_len2 += myprior_z * myprior_z;
											// If it is doing helical refinement AND Cartesian vector myprior has a length > 0, transform the vector to its helical coordinates
											if ( (do_helical_refine) && (!ignore_helical_symmetry) && (mypriors_len2 > 0.00001) )
											{
												RFLOAT rot_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_ROT);
												RFLOAT tilt_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_TILT);
												RFLOAT psi_deg = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PSI);
												transformCartesianAndHelicalCoords(myprior_x, myprior_y, myprior_z, myprior_x, myprior_y, myprior_z, rot_deg, tilt_deg, psi_deg, mymodel.data_dim, CART_TO_HELICAL_COORDS);
											}

											if ( (!do_helical_refine) || (ignore_helical_symmetry) || (mymodel.data_dim == 3) )
											{
												RFLOAT diffx = myprior_x - old_offset_x - oversampled_translations_x[iover_trans];
												thr_wsum_sigma2_offset += weight * diffx * diffx;
											}
											RFLOAT diffy = myprior_y - old_offset_y - oversampled_translations_y[iover_trans];
											thr_wsum_sigma2_offset += weight * diffy * diffy;
											if (mymodel.data_dim == 3)
											{
												RFLOAT diffz  = myprior_z - old_offset_z - oversampled_translations_z[iover_trans];
												if ( (!do_helical_refine) || (ignore_helical_symmetry) )
													thr_wsum_sigma2_offset += weight * diffz * diffz;
											}

											// Store weight for this direction of this class
											if (do_skip_align || do_skip_rotate )
											{
												//ignore pdf_direction
											}
											else if (mymodel.orientational_prior_mode == NOPRIOR)
											{
#ifdef DEBUG_CHECKSIZES
												if (idir >= XSIZE(thr_wsum_pdf_direction[exp_iclass]))
												{
													std::cerr<< "idir= "<<idir<<" XSIZE(thr_wsum_pdf_direction[exp_iclass])= "<< XSIZE(thr_wsum_pdf_direction[exp_iclass]) <<std::endl;
													REPORT_ERROR("idir >= XSIZE(thr_wsum_pdf_direction[iclass])");
												}
#endif
												DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], idir) += weight;
											}
											else
											{
												// In the case of orientational priors, get the original number of the direction back
												long int mydir = exp_pointer_dir_nonzeroprior[idir];
												DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += weight;
											}

#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
											{
												timer.toc(TIMING_WSUM_LOCALSUMS);
												timer.tic(TIMING_WSUM_SUMSHIFT);
											}
#endif

											Complex *Fimg_store;
											if (do_sgd)
											{
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Frefctf)
												{
													(DIRECT_MULTIDIM_ELEM(Fimg_store_sgd, n)).real = (*(Fimg_shift_nomask + n)).real - (DIRECT_MULTIDIM_ELEM(Frefctf, n)).real;
													(DIRECT_MULTIDIM_ELEM(Fimg_store_sgd, n)).imag = (*(Fimg_shift_nomask + n)).imag - (DIRECT_MULTIDIM_ELEM(Frefctf, n)).imag;
												}
												Fimg_store = Fimg_store_sgd.data;
											}
											else
											{
												Fimg_store = Fimg_shift_nomask;
											}
//#define DEBUG_SGD3
#ifdef DEBUG_SGD3
											FourierTransformer transformer;
											MultidimArray<Complex> Ftt(Frefctf);
											FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Ftt)
												DIRECT_MULTIDIM_ELEM(Ftt, n) = *(Fimg_store + n);

											Image<RFLOAT> tt;
											tt().resize(exp_current_image_size, exp_current_image_size);
											transformer.inverseFourierTransform(Ftt, tt());
											CenterFFT(tt(),false);
											tt.write("BPimg.spi");
											std::cerr << " rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;
											std::cerr << " itrans= " << itrans << " iover_trans= " << iover_trans << std::endl;
											std::cerr << "written BPimg.spi, press any key to continue" << std::endl;
											char c;
											std::cin >> c;
#endif

											// Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
											// Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
											if (ctf_premultiplied)
											{
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
												{
													RFLOAT myctf = DIRECT_MULTIDIM_ELEM(Mctf, n);
													RFLOAT weightxinvsigma2 = weight * DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
													// now Fimg stores sum of all shifted w*Fimg
													(DIRECT_MULTIDIM_ELEM(Fimg, n)).real += (*(Fimg_store + n)).real * weightxinvsigma2;
													(DIRECT_MULTIDIM_ELEM(Fimg, n)).imag += (*(Fimg_store + n)).imag * weightxinvsigma2;
													// now Fweight stores sum of all w and multiply by CTF^2
													DIRECT_MULTIDIM_ELEM(Fweight, n) += weightxinvsigma2 * myctf * myctf;
												}
											}
											else
											{
												FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fimg)
												{
													RFLOAT myctf = DIRECT_MULTIDIM_ELEM(Mctf, n);
													RFLOAT weightxinvsigma2 = weight * myctf * DIRECT_MULTIDIM_ELEM(Minvsigma2, n);
													// now Fimg stores sum of all shifted w*Fimg
													(DIRECT_MULTIDIM_ELEM(Fimg, n)).real += (*(Fimg_store + n)).real * weightxinvsigma2;
													(DIRECT_MULTIDIM_ELEM(Fimg, n)).imag += (*(Fimg_store + n)).imag * weightxinvsigma2;
													// now Fweight stores sum of all w
													// Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
													DIRECT_MULTIDIM_ELEM(Fweight, n) += weightxinvsigma2 * myctf;

												}
											}

#ifdef TIMING
											// Only time one thread, as I also only time one MPI process
											if (my_ori_particle == exp_my_first_ori_particle)
												timer.toc(TIMING_WSUM_SUMSHIFT);
#endif
										} // end if !do_skip_maximization

										// Keep track of max_weight and the corresponding optimal hidden variables
										if (weight > exp_max_weight[ipart])
										{
											// Store optimal image parameters
											exp_max_weight[ipart] = weight;

											Euler_matrix2angles(A, rot, tilt, psi);

											int icol_rot  = (mymodel.nr_bodies == 1) ? METADATA_ROT  : 0 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
											int icol_tilt = (mymodel.nr_bodies == 1) ? METADATA_TILT : 1 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
											int icol_psi  = (mymodel.nr_bodies == 1) ? METADATA_PSI  : 2 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
											int icol_xoff = (mymodel.nr_bodies == 1) ? METADATA_XOFF : 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
											int icol_yoff = (mymodel.nr_bodies == 1) ? METADATA_YOFF : 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
											int icol_zoff = (mymodel.nr_bodies == 1) ? METADATA_ZOFF : 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;

											RFLOAT old_rot = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_rot);
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_rot) = rot;
											RFLOAT old_tilt = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_tilt);
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_tilt) = tilt;
											RFLOAT old_psi = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_psi);
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_psi) = psi;

											Matrix1D<RFLOAT> shifts(2);
											//21may2015
											if (mymodel.nr_bodies == 1)
											{
												// include old_offsets for normal refinement (i.e. non multi-body)
												XX(shifts) = XX(exp_old_offset[ipart]) + oversampled_translations_x[iover_trans];
												YY(shifts) = YY(exp_old_offset[ipart]) + oversampled_translations_y[iover_trans];
											}
											else
											{
												// For multi-body refinements, only store 'residual' translations
												XX(shifts) = oversampled_translations_x[iover_trans];
												YY(shifts) = oversampled_translations_y[iover_trans];
											}
											if (mymodel.data_dim == 3)
											{
												shifts.resize(3);
												if (mymodel.nr_bodies == 1)
													ZZ(shifts) = ZZ(exp_old_offset[ipart]) + oversampled_translations_z[iover_trans];
												else
													ZZ(shifts) = oversampled_translations_z[iover_trans];
											}
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
											std::cerr << "MlOptimiser::storeWeightedSums()" << std::endl;
											if (mymodel.data_dim == 2)
											{
												std::cerr << " exp_old_offset = (" << XX(exp_old_offset[ipart]) << ", " << YY(exp_old_offset[ipart]) << ")" << std::endl;
												std::cerr << " Oversampled trans = (" << oversampled_translations_x[iover_trans] << ", " << oversampled_translations_y[iover_trans] << ")" << std::endl;
												std::cerr << " shifts = (" << XX(shifts) << ", " << YY(shifts) << ")" << std::endl;
											}
											else
											{
												std::cerr << " exp_old_offset = (" << XX(exp_old_offset[ipart]) << ", " << YY(exp_old_offset[ipart]) << ", " << ZZ(exp_old_offset[ipart]) << ")" << std::endl;
												std::cerr << " Oversampled trans = (" << oversampled_translations_x[iover_trans] << ", " << oversampled_translations_y[iover_trans] << ", " << oversampled_translations_z[iover_trans] << ")" << std::endl;
												std::cerr << " shifts = (" << XX(shifts) << ", " << YY(shifts) << ", " << ZZ(shifts) << ")" << std::endl;
											}
#endif

											// Helical reconstruction: use oldpsi-angle to rotate back the XX(exp_old_offset[ipart]) + oversampled_translations_x[iover_trans] and
											if ( (do_helical_refine) && (!ignore_helical_symmetry) )
											{
												// Bring xshift, yshift and zshift back to cartesian coords for outputting in the STAR file
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
												std::cerr << "MlOptimiser::storeWeightedSums()" << std::endl;
												std::cerr << "Bring xy(z) shifts back to Cartesian coordinates for output in the STAR file" << std::endl;
												std::cerr << " itrans = " << itrans << ", iover_trans = " << iover_trans << std::endl;
												if(shifts.size() == 2)
												{
													std::cerr << "  old_psi = " << old_psi << " degrees" << std::endl;
													std::cerr << "  Helical offsets (r, p) = (" << XX(shifts) << ", " << YY(shifts) << ")" << std::endl;
												}
												else
												{
													std::cerr << "  old_psi = " << old_psi << " degrees, old_tilt = " << old_tilt << " degrees" << std::endl;
													std::cerr << "  Helical offsets (p1, p2, r) = (" << XX(shifts) << ", " << YY(shifts) << ", " << ZZ(shifts) << ")" << std::endl;
												}
#endif
												transformCartesianAndHelicalCoords(shifts, shifts, old_rot, old_tilt, old_psi, HELICAL_TO_CART_COORDS);
#ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
												if(shifts.size() == 2)
													std::cerr << "  Cartesian offsets (x, y) = (" << XX(shifts) << ", " << YY(shifts) << ")" << std::endl;
												else
													std::cerr << "  Cartesian offsets (x, y, z) = (" << XX(shifts) << ", " << YY(shifts) << ", " << ZZ(shifts) << ")" << std::endl;
#endif
											}

											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_xoff) = XX(shifts);
											DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_yoff) = YY(shifts);
											if (mymodel.data_dim == 3)
												DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, icol_zoff) = ZZ(shifts);

											if (ibody == 0)
											{
												DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_CLASS) = (RFLOAT)exp_iclass + 1;
												DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PMAX) = exp_max_weight[ipart];
											}
										}
									} // end if weight >= exp_significant_weight
								} // end loop iover_trans
							} // end loop itrans
						} // end loop ipart
#ifdef RELION_TESTING
		std::string fnm = std::string("cpu_out_exp_wsum_norm_correction.txt");
		char *text = &fnm[0];
		freopen(text,"w",stdout);
		for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
		{
			printf("%4.8f \n",exp_wsum_norm_correction[ipart]);
		}
		fclose(stdout);
		//----------
		fnm = std::string("cpu_out_thr_wsum_sigma2_noise.txt");
		text = &fnm[0];
		freopen(text,"w",stdout);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
		{
			printf("%4.8f \n",thr_wsum_sigma2_noise[0].data[n]);
		}
		fclose(stdout);
		//----------
		fnm = std::string("cpu_out_Fweights.txt");
		text = &fnm[0];
		freopen(text,"w",stdout);
		for(int n = 0; n < 1000; n++)
		{
			printf("%4.8f \n",*(Fweight.data + n*60+50));
		}
		fclose(stdout);
#endif

#ifdef DEBUG_SGD3
					FourierTransformer transformer;
					Image<RFLOAT> tt;
					tt().resize(exp_current_image_size, exp_current_image_size);
					transformer.inverseFourierTransform(Fimg, tt());
					CenterFFT(tt(),false);
					tt.write("BPimg_sig2.spi");
					std::cerr << "written BPimg_sig2.spi, press any key to continue" << std::endl;
					char c;
					std::cin >> c;
#endif
						if (!do_skip_maximization)
						{
#ifdef TIMING
                                                    // Only time one thread, as I also only time one MPI process
                                                    if (my_ori_particle == exp_my_first_ori_particle)
								timer.tic(TIMING_WSUM_BACKPROJ);
#endif
							// Perform the actual back-projection.
							// This is done with the sum of all (in-plane) shifted Fimg's
							// Perform this inside a mutex
							int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
							pthread_mutex_lock(&global_mutex2[my_mutex]);
							if (mymodel.nr_bodies > 1)
							{
								//19may2015: place Fimg so that it is centered at the COM of the body
								Matrix1D<RFLOAT> comp(3);
								comp = Abody * (-mymodel.com_bodies[ibody]);
								shiftImageInFourierTransform(Fimg, Fimg, (RFLOAT)mymodel.ori_size, XX(comp), YY(comp), ZZ(comp));
								(wsum_model.BPref[ibody]).set2DFourierTransform(Fimg, Abody, IS_NOT_INV, &Fweight);
							}
							else
								(wsum_model.BPref[exp_iclass]).set2DFourierTransform(Fimg, A, IS_NOT_INV, &Fweight);
							pthread_mutex_unlock(&global_mutex2[my_mutex]);
#ifdef TIMING
							// Only time one thread, as I also only time one MPI process
                                                        if (my_ori_particle == exp_my_first_ori_particle)
								timer.toc(TIMING_WSUM_BACKPROJ);
#endif
						} // end if !do_skip_maximization
					} // end loop iover_rot
				}// end loop do_proceed
			} // end loop ipsi
		} // end loop idir
	} // end loop iclass


	// Extend norm_correction and sigma2_noise estimation to higher resolutions for all particles
	// Also calculate dLL for each particle and store in metadata
	// loop over all particles inside this ori_particle
	RFLOAT thr_avg_norm_correction = 0.;
	RFLOAT thr_sum_dLL = 0., thr_sum_Pmax = 0.;
	for (long int ipart = 0; ipart < mydata.ori_particles[my_ori_particle].particles_id.size(); ipart++)
	{
		long int part_id = mydata.ori_particles[my_ori_particle].particles_id[ipart];
		int group_id = mydata.getGroupId(part_id);

		// If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
		for (int ires = mymodel.current_size/2 + 1; ires < mymodel.ori_size/2 + 1; ires++)
		{
			DIRECT_A1D_ELEM(thr_wsum_sigma2_noise[group_id], ires) += DIRECT_A1D_ELEM(exp_power_imgs[ipart], ires);
			// Also extend the weighted sum of the norm_correction
			exp_wsum_norm_correction[ipart] += DIRECT_A1D_ELEM(exp_power_imgs[ipart], ires);
		}

		// Store norm_correction
		// Multiply by old value because the old norm_correction term was already applied to the image
		if (do_norm_correction)
		{
			RFLOAT old_norm_correction = DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NORM);
			old_norm_correction /= mymodel.avg_norm_correction;
			// The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
			// The variance of the total image (on which one normalizes) is twice this value!
			RFLOAT normcorr = old_norm_correction * sqrt(exp_wsum_norm_correction[ipart] * 2.);
			thr_avg_norm_correction += normcorr;
			// Now set the new norm_correction in the relevant position of exp_metadata
			DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NORM) = normcorr;

			// Print warning for strange norm-correction values
			if (!((iter == 1 && do_firstiter_cc) || do_always_cc) && DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NORM) > 10.)
			{
				std::cout << " WARNING: norm_correction= "<< DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_NORM)
						<< " for particle " << part_id << " in group " << group_id + 1
						<< "; Are your groups large enough?  Or is the reference on the correct greyscale?" << std::endl;
			}

			//TMP DEBUGGING
			/*
			if (!(iter == 1 && do_firstiter_cc) && DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) > 10.)
			{
				std::cerr << " mymodel.current_size= " << mymodel.current_size << " mymodel.ori_size= " << mymodel.ori_size << " part_id= " << part_id << std::endl;
				std::cerr << " DIRECT_A2D_ELEM(exp_metadata, ipart, METADATA_NORM)= " << DIRECT_A2D_ELEM(exp_metadata, ipart, METADATA_NORM) << std::endl;
				std::cerr << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
				std::cerr << " exp_wsum_norm_correction[ipart]= " << exp_wsum_norm_correction[ipart] << std::endl;
				std::cerr << " old_norm_correction= " << old_norm_correction << std::endl;
				std::cerr << " wsum_model.avg_norm_correction= " << wsum_model.avg_norm_correction << std::endl;
				std::cerr << " group_id= " << group_id << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
				std::cerr << " mymodel.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
				std::cerr << " wsum_model.sigma2_noise[group_id]= " << wsum_model.sigma2_noise[group_id] << std::endl;
				std::cerr << " exp_power_imgs[ipart]= " << exp_power_imgs[ipart] << std::endl;
				std::cerr << " exp_wsum_scale_correction_XA[ipart]= " << exp_wsum_scale_correction_XA[ipart] << " exp_wsum_scale_correction_AA[ipart]= " << exp_wsum_scale_correction_AA[ipart] << std::endl;
				std::cerr << " wsum_model.wsum_signal_product_spectra[group_id]= " << wsum_model.wsum_signal_product_spectra[group_id] << " wsum_model.wsum_reference_power_spectra[group_id]= " << wsum_model.wsum_reference_power_spectra[group_id] << std::endl;
				std::cerr << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
				std::cerr << " ml_model.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
				std::cerr << " exp_significant_weight[ipart]= " << exp_significant_weight[ipart] << std::endl;
				std::cerr << " exp_max_weight[ipart]= " << exp_max_weight[ipart] << std::endl;
				mymodel.write("debug");
				std::cerr << "written debug_model.star" << std::endl;
				REPORT_ERROR("MlOptimiser::storeWeightedSums ERROR: normalization is larger than 10");
			}
			*/

		}

		// Store weighted sums for scale_correction
		if (do_scale_correction)
		{
			// Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
			exp_wsum_scale_correction_XA[ipart] /= mymodel.scale_correction[group_id];
			exp_wsum_scale_correction_AA[ipart] /= mymodel.scale_correction[group_id] * mymodel.scale_correction[group_id];

			thr_wsum_signal_product_spectra[group_id] += exp_wsum_scale_correction_XA[ipart];
			thr_wsum_reference_power_spectra[group_id] += exp_wsum_scale_correction_AA[ipart];
		}

		// Calculate DLL for each particle
		RFLOAT logsigma2 = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mresol_fine)
		{
			int ires = DIRECT_MULTIDIM_ELEM(Mresol_fine, n);
			// Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
			// Also exclude origin from logsigma2, as this will not be considered in the P-calculations
			if (ires > 0)
				logsigma2 += log( 2. * PI * DIRECT_A1D_ELEM(mymodel.sigma2_noise[group_id], ires));
		}
		if (exp_sum_weight[ipart]==0)
		{
			std::cerr << " part_id= " << part_id << std::endl;
			std::cerr << " ipart= " << ipart << std::endl;
			std::cerr << " exp_min_diff2[ipart]= " << exp_min_diff2[ipart] << std::endl;
			std::cerr << " logsigma2= " << logsigma2 << std::endl;
			int group_id = mydata.getGroupId(part_id);
			std::cerr << " group_id= " << group_id << std::endl;
			std::cerr << " ml_model.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
			std::cerr << " exp_significant_weight[ipart]= " << exp_significant_weight[ipart] << std::endl;
			std::cerr << " exp_max_weight[ipart]= " << exp_max_weight[ipart] << std::endl;
			std::cerr << " ml_model.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
			REPORT_ERROR("ERROR: exp_sum_weight[ipart]==0");
		}
		RFLOAT dLL;
		if ((iter==1 && do_firstiter_cc) || do_always_cc)
			dLL = -exp_min_diff2[ipart];
		else
			dLL = log(exp_sum_weight[ipart]) - exp_min_diff2[ipart] - logsigma2;

		// Store dLL of each image in the output array, and keep track of total sum
		DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_DLL) = dLL;
		thr_sum_dLL += dLL;

		// Also store sum of Pmax
		thr_sum_Pmax += DIRECT_A2D_ELEM(exp_metadata, metadata_offset + ipart, METADATA_PMAX);

	}

	// Now, inside a global_mutex, update the other weighted sums among all threads
	if (!do_skip_maximization)
	{
		pthread_mutex_lock(&global_mutex);
		for (int n = 0; n < mymodel.nr_groups; n++)
		{
			wsum_model.sigma2_noise[n] += thr_wsum_sigma2_noise[n];
			wsum_model.sumw_group[n] += thr_sumw_group[n];
			if (do_scale_correction)
			{
				wsum_model.wsum_signal_product_spectra[n] += thr_wsum_signal_product_spectra[n];
				wsum_model.wsum_reference_power_spectra[n] += thr_wsum_reference_power_spectra[n];
			}
		}
		for (int n = 0; n < mymodel.nr_classes; n++)
		{
			wsum_model.pdf_class[n] += thr_wsum_pdf_class[n];
			if (mymodel.ref_dim == 2)
			{
				XX(wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsetx_class[n];
				YY(wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsety_class[n];
			}
#ifdef CHECKSIZES
			if (XSIZE(wsum_model.pdf_direction[n]) != XSIZE(thr_wsum_pdf_direction[n]))
			{
				std::cerr << " XSIZE(wsum_model.pdf_direction[n])= " << XSIZE(wsum_model.pdf_direction[n]) << " XSIZE(thr_wsum_pdf_direction[n])= " << XSIZE(thr_wsum_pdf_direction[n]) << std::endl;
				REPORT_ERROR("XSIZE(wsum_model.pdf_direction[n]) != XSIZE(thr_wsum_pdf_direction[n])");
			}
#endif
			if (!(do_skip_align || do_skip_rotate) )
				wsum_model.pdf_direction[n] += thr_wsum_pdf_direction[n];
		}
		wsum_model.sigma2_offset += thr_wsum_sigma2_offset;
		if (do_norm_correction)
			wsum_model.avg_norm_correction += thr_avg_norm_correction;
		wsum_model.LL += thr_sum_dLL;
		wsum_model.ave_Pmax += thr_sum_Pmax;
		pthread_mutex_unlock(&global_mutex);
	} // end if !do_skip_maximization


#ifdef DEBUG_OVERSAMPLING
	std::cerr << " max_weight= " << max_weight << " nr_sign_sam= "<<nr_significant_samples<<" sign w= "<<exp_significant_weight<<std::endl;
#endif

#ifdef TIMING
	if (my_ori_particle == exp_my_first_ori_particle)
		timer.toc(TIMING_ESP_WSUM);
#endif

}

/** Monitor the changes in the optimal translations, orientations and class assignments for some particles */
void MlOptimiser::monitorHiddenVariableChanges(long int my_first_ori_particle, long int my_last_ori_particle)
{

	for (long int ori_part_id = my_first_ori_particle, my_image_no = 0; ori_part_id <= my_last_ori_particle; ori_part_id++)
	{

#ifdef DEBUG_CHECKSIZES
		if (ori_part_id >= mydata.ori_particles.size())
		{
			std::cerr<< "ori_part_id= "<<ori_part_id<<" mydata.ori_particles.size()= "<< mydata.ori_particles.size() <<std::endl;
			REPORT_ERROR("ori_part_id >= mydata.ori_particles.size()");
		}
#endif

		// loop over all particles inside this ori_particle
		for (long int ipart = 0; ipart < mydata.ori_particles[ori_part_id].particles_id.size(); ipart++, my_image_no++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[ipart];

#ifdef DEBUG_CHECKSIZES
			if (part_id >= mydata.MDimg.numberOfObjects())
			{
				std::cerr<< "part_id= "<<part_id<<" mydata.MDimg.numberOfObjects()= "<< mydata.MDimg.numberOfObjects() <<std::endl;
				REPORT_ERROR("part_id >= mydata.MDimg.numberOfObjects()");
			}
			if (my_image_no >= YSIZE(exp_metadata))
			{
				std::cerr<< "my_image_no= "<<my_image_no<<" YSIZE(exp_metadata)= "<< YSIZE(exp_metadata) <<std::endl;
				REPORT_ERROR("my_image_no >= YSIZE(exp_metadata)");
			}
#endif

			// Old optimal parameters
			RFLOAT old_rot, old_tilt, old_psi, old_xoff, old_yoff, old_zoff = 0.;
			int old_iclass;
			mydata.MDimg.getValue(EMDL_ORIENT_ROT,  old_rot, part_id);
			mydata.MDimg.getValue(EMDL_ORIENT_TILT, old_tilt, part_id);
			mydata.MDimg.getValue(EMDL_ORIENT_PSI,  old_psi, part_id);
			mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, old_xoff, part_id);
			mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, old_yoff, part_id);
			if (mymodel.data_dim == 3)
				mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Z, old_zoff, part_id);
			mydata.MDimg.getValue(EMDL_PARTICLE_CLASS, old_iclass, part_id);

			// New optimal parameters
			RFLOAT rot = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT);
			RFLOAT tilt = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT);
			RFLOAT psi = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI);
			RFLOAT xoff = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF);
			RFLOAT yoff = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF);
			RFLOAT zoff = 0.;
			if (mymodel.data_dim == 3)
				zoff = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ZOFF);
			int iclass = (int)DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CLASS);

			// Some orientational distance....
			sum_changes_optimal_orientations += sampling.calculateAngularDistance(rot, tilt, psi, old_rot, old_tilt, old_psi);
			sum_changes_optimal_offsets += (xoff-old_xoff)*(xoff-old_xoff) + (yoff-old_yoff)*(yoff-old_yoff) + (zoff-old_zoff)*(zoff-old_zoff);
			if (iclass != old_iclass)
				sum_changes_optimal_classes += 1.;
			sum_changes_count += 1.;
		} // end loop part_id (i)
	} //end loop ori_part_id


}

void MlOptimiser::updateOverallChangesInHiddenVariables()
{

	// Calculate hidden variable changes
	current_changes_optimal_classes = sum_changes_optimal_classes / sum_changes_count;
	current_changes_optimal_orientations = sum_changes_optimal_orientations / sum_changes_count;
	current_changes_optimal_offsets = sqrt(sum_changes_optimal_offsets / (2. * sum_changes_count));

	// Reset the sums
	sum_changes_optimal_classes = 0.;
	sum_changes_optimal_orientations = 0.;
	sum_changes_optimal_offsets = 0.;
	sum_changes_count = 0.;

	RFLOAT ratio_orient_changes = current_changes_optimal_orientations /  sampling.getAngularSampling(adaptive_oversampling);
	RFLOAT ratio_trans_changes = current_changes_optimal_offsets /  sampling.getTranslationalSampling(adaptive_oversampling);

	// Update nr_iter_wo_large_hidden_variable_changes if all three assignment types are within 3% of the smallest thus far
	// Or if changes in offsets or orientations are smaller than 40% of the current sampling
	if (1.03 * current_changes_optimal_classes >= smallest_changes_optimal_classes &&
		(ratio_trans_changes < 0.40 || 1.03 * current_changes_optimal_offsets >= smallest_changes_optimal_offsets) &&
		(ratio_orient_changes < 0.40 || 1.03 * current_changes_optimal_orientations >= smallest_changes_optimal_orientations) )
		nr_iter_wo_large_hidden_variable_changes++;
	else
		nr_iter_wo_large_hidden_variable_changes = 0;

	// Update smallest changes in hidden variables thus far
	if (current_changes_optimal_classes < smallest_changes_optimal_classes)
		smallest_changes_optimal_classes = ROUND(current_changes_optimal_classes);
	if (current_changes_optimal_offsets < smallest_changes_optimal_offsets)
		smallest_changes_optimal_offsets = current_changes_optimal_offsets;
	if (current_changes_optimal_orientations < smallest_changes_optimal_orientations)
		smallest_changes_optimal_orientations = current_changes_optimal_orientations;


}


void MlOptimiser::calculateExpectedAngularErrors(long int my_first_ori_particle, long int my_last_ori_particle)
{

	long int n_trials = 0;
	for (long int ori_part_id = my_first_ori_particle; ori_part_id <= my_last_ori_particle; ori_part_id++)
    {
		n_trials +=  mydata.ori_particles[ori_part_id].particles_id.size();
    }

	//int current_image_size = (strict_highres_exp > 0. && !do_acc_currentsize_despite_highres_exp) ? coarse_size : mymodel.current_size;
	// Set current_image_size to the coarse_size to calculate exepcted angular errors
	int current_image_size;
	if (strict_highres_exp > 0. && !do_acc_currentsize_despite_highres_exp)
	{
		// Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
		current_image_size = coarse_size;
	}
	else
	{
		// Use smaller images in the first pass, but larger ones in the second pass
		current_image_size = mymodel.current_size;
	}

	// Separate angular error estimate for each of the classes
	acc_rot = acc_trans = 999.; // later XMIPP_MIN will be taken to find the best class...

	// P(X | X_1) / P(X | X_2) = exp ( |F_1 - F_2|^2 / (-2 sigma2) )
	// exp(-4.60517) = 0.01
	RFLOAT pvalue = 4.60517;
	//if (mymodel.data_dim == 3)
	//	pvalue *= 2.;

	std::cout << " Estimating accuracies in the orientational assignment ... " << std::endl;
	init_progress_bar(n_trials * mymodel.nr_classes);
	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
	{

		// Don't do this for (almost) empty classes
		if (mymodel.pdf_class[iclass] < 0.01)
		{
			mymodel.acc_rot[iclass]   = 999.;
			mymodel.acc_trans[iclass] = 999.;
			continue;
		}

		// Initialise the orientability arrays that will be written out in the model.star file
		// These are for the user's information only: nothing will be actually done with them
#ifdef DEBUG_CHECKSIZES
		if (iclass >= (mymodel.orientability_contrib).size())
		{
			std::cerr<< "iclass= "<<iclass<<" (mymodel.orientability_contrib).size()= "<< (mymodel.orientability_contrib).size() <<std::endl;
			REPORT_ERROR("iclass >= (mymodel.orientability_contrib).size()");
		}
#endif
		(mymodel.orientability_contrib)[iclass].initZeros(mymodel.ori_size/2 + 1);

		RFLOAT acc_rot_class = 0.;
		RFLOAT acc_trans_class = 0.;
		// Particles are already in random order, so just move from 0 to n_trials
		for (long int ori_part_id = my_first_ori_particle, my_metadata_entry = 0, ipart = 0; ori_part_id <= my_last_ori_particle; ori_part_id++)
	    {
			for (long int ipart = 0; ipart < mydata.ori_particles[ori_part_id].particles_id.size(); ipart++, my_metadata_entry++)
			{
				long int part_id = mydata.ori_particles[ori_part_id].particles_id[ipart];


				MultidimArray<RFLOAT> Fctf;
				// Get CTF for this particle
				if (do_ctf_correction)
				{
					if (mymodel.data_dim == 3)
					{
						Image<RFLOAT> Ictf;
						// Read CTF-image from disc
						FileName fn_ctf;
						if (!mydata.getImageNameOnScratch(part_id, fn_ctf, true))
						{
							std::istringstream split(exp_fn_ctf);
							// Get the right line in the exp_fn_img string
							for (int i = 0; i <= my_metadata_entry; i++)
								getline(split, fn_ctf);
						}
						Ictf.read(fn_ctf);

						// Set the CTF-image in Fctf
						Ictf().setXmippOrigin();
						Fctf.resize(current_image_size, current_image_size, current_image_size/ 2 + 1);
						FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf)
						{
							// Use negative kp, ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
							DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
						}
					}
					else
					{
						CTF ctf;
						ctf.setValues(DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_CTF_DEFOCUS_U),
									  DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_CTF_DEFOCUS_V),
									  DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_CTF_DEFOCUS_ANGLE),
									  DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_CTF_VOLTAGE),
									  DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_CTF_CS),
									  DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_CTF_Q0),
									  DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_CTF_BFAC),
									  1.,
									  DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_CTF_PHASE_SHIFT));

						Fctf.resize(current_image_size, current_image_size/ 2 + 1);
						ctf.getFftwImage(Fctf, mymodel.ori_size, mymodel.ori_size, mymodel.pixel_size, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);
					}
				}

				// Search 2 times: ang and off
				// Don't estimate rotational accuracies if we're doing do_skip_rotate (for faster movie-frame alignment)
				int imode_start = (do_skip_rotate) ? 1 : 0;
				for (int imode = imode_start; imode < 2; imode++)
				{
					RFLOAT ang_error = 0.;
					RFLOAT sh_error = 0.;
					RFLOAT ang_step;
					RFLOAT sh_step;
					RFLOAT my_snr = 0.;

					// Search for ang_error and sh_error where there are at least 3-sigma differences!
					// 13feb12: change for explicit probability at P=0.01
					while (my_snr <= pvalue)
					{
						// Graduallly increase the step size
						if (ang_error < 0.2)
							ang_step = 0.05;
						else if (ang_error < 1.)
							ang_step = 0.1;
						else if (ang_error < 2.)
							ang_step = 0.2;
						else if (ang_error < 5.)
							ang_step = 0.5;
						else if (ang_error < 10.)
							ang_step = 1.0;
						else if (ang_error < 20.)
							ang_step = 2;
						else
							ang_step = 5.0;

						if (sh_error < 0.2)
							sh_step = 0.05;
						else if (sh_error < 1.)
							sh_step = 0.1;
						else if (sh_error < 2.)
							sh_step = 0.2;
						else if (sh_error < 5.)
							sh_step = 0.5;
						else if (sh_error < 10.)
							sh_step = 1.0;
						else
							sh_step = 2.0;

						ang_error += ang_step;
						sh_error += sh_step;

						// Prevent an endless while by putting boundaries on ang_error and sh_error
						if ( (imode == 0 && ang_error > 30.) || (imode == 1 && sh_error > 10.) )
							break;

						init_random_generator(random_seed + part_id);

						int group_id = mydata.getGroupId(part_id);
#ifdef DEBUG_CHECKSIZES
						if (group_id  >= mymodel.sigma2_noise.size())
						{
							std::cerr<< "group_id = "<<group_id <<" mymodel.sigma2_noise.size()= "<< mymodel.sigma2_noise.size() <<std::endl;
							REPORT_ERROR("group_id  >= mymodel.sigma2_noise.size()");
						}
#endif
						MultidimArray<Complex > F1, F2;
						Matrix2D<RFLOAT> A1, A2;

						RFLOAT rot1 = DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_ROT);
						RFLOAT tilt1 = DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_TILT);
						RFLOAT psi1 = DIRECT_A2D_ELEM(exp_metadata, my_metadata_entry, METADATA_PSI);
						RFLOAT xoff1 = 0.;
						RFLOAT yoff1 = 0.;
						RFLOAT zoff1 = 0.;

						if (mymodel.data_dim == 2)
							F1.initZeros(current_image_size, current_image_size/ 2 + 1);
						else
							F1.initZeros(current_image_size, current_image_size, current_image_size/ 2 + 1);

						// Get the FT of the first image
						Euler_angles2matrix(rot1, tilt1, psi1, A1, false);
						(mymodel.PPref[iclass]).get2DFourierTransform(F1, A1, IS_NOT_INV);

						// Apply the angular or shift error
						RFLOAT rot2 = rot1;
						RFLOAT tilt2 = tilt1;
						RFLOAT psi2 = psi1;
						RFLOAT xshift = xoff1;
						RFLOAT yshift = yoff1;
						RFLOAT zshift = zoff1;

						// Perturb psi or xoff , depending on the mode
						if (imode == 0)
						{
							if (mymodel.ref_dim == 3)
							{
								// Randomly change rot, tilt or psi
								RFLOAT ran = rnd_unif();
								if (ran < 0.3333)
									rot2 = rot1 + ang_error;
								else if (ran < 0.6667)
									tilt2 = tilt1 + ang_error;
								else
									psi2  = psi1 + ang_error;
							}
							else
							{
								psi2  = psi1 + ang_error;
							}
						}
						else
						{
							// Randomly change xoff or yoff
							RFLOAT ran = rnd_unif();
							if (mymodel.data_dim == 3)
							{
								if (ran < 0.3333)
									xshift = xoff1 + sh_error;
								else if (ran < 0.6667)
									yshift = yoff1 + sh_error;
								else
									zshift = zoff1 + sh_error;
							}
							else
							{
								if (ran < 0.5)
									xshift = xoff1 + sh_error;
								else
									yshift = yoff1 + sh_error;
							}
						}
						// Get the FT of the second image
						if (mymodel.data_dim == 2)
							F2.initZeros(current_image_size, current_image_size/ 2 + 1);
						else
							F2.initZeros(current_image_size, current_image_size, current_image_size/ 2 + 1);

						if (imode == 0)
						{
							// Get new rotated version of reference
							Euler_angles2matrix(rot2, tilt2, psi2, A2, false);
							(mymodel.PPref[iclass]).get2DFourierTransform(F2, A2, IS_NOT_INV);
						}
						else
						{
							// Get shifted version
							shiftImageInFourierTransform(F1, F2, (RFLOAT) mymodel.ori_size, -xshift, -yshift, -zshift);
						}

						// Apply CTF to F1 and F2 if necessary
						if (do_ctf_correction)
						{
#ifdef DEBUG_CHECKSIZES
							if (!Fctf.sameShape(F1) || !Fctf.sameShape(F2))
							{
								std::cerr<<" Fctf: "; Fctf.printShape(std::cerr);
								std::cerr<<" F1:   "; F1.printShape(std::cerr);
								std::cerr<<" F2:   "; F2.printShape(std::cerr);
								REPORT_ERROR("ERROR: Fctf has a different shape from F1 and F2");
							}
#endif
							FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
							{
								DIRECT_MULTIDIM_ELEM(F1, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
								DIRECT_MULTIDIM_ELEM(F2, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
							}
							if (ctf_premultiplied)
							{
								FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
								{
									DIRECT_MULTIDIM_ELEM(F1, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
									DIRECT_MULTIDIM_ELEM(F2, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
								}
							}
						}

						MultidimArray<int> * myMresol = (YSIZE(F1) == coarse_size) ? &Mresol_coarse : &Mresol_fine;
						my_snr = 0.;
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
						{
							int ires = DIRECT_MULTIDIM_ELEM(*myMresol, n);
							if (ires > 0)
							{
								my_snr += norm(DIRECT_MULTIDIM_ELEM(F1, n) - DIRECT_MULTIDIM_ELEM(F2, n)) / (2 * sigma2_fudge * mymodel.sigma2_noise[group_id](ires) );
							}
						}
//#define DEBUG_ANGACC
#ifdef DEBUG_ANGACC
							if (imode==0)
							{
								std::cerr << " ang_error= " << ang_error << std::endl;
								std::cerr << " rot1= " << rot1 << " tilt1= " << tilt1 << " psi1= " << psi1 << std::endl;
								std::cerr << " rot2= " << rot2 << " tilt2= " << tilt2 << " psi2= " << psi2 << std::endl;

							}
							else
							{
								std::cerr << " xshift= " << xshift << " yshift= " << yshift << " zshift= " << zshift << std::endl;
								std::cerr << " sh_error= " << sh_error << std::endl;
							}
							std::cerr << " my_snr= " << my_snr << std::endl;
							FourierTransformer transformer;
							MultidimArray<RFLOAT> spec_img(mymodel.ori_size), spec_diff(mymodel.ori_size), count(mymodel.ori_size);
							FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F1)
							{
								long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
								spec_img(idx) += norm(dAkij(F1, k, i, j));
								spec_diff(idx) += norm(dAkij(F1, k, i, j) - dAkij(F2, k, i, j));
								count(idx) += 1.;
							}
							spec_img /= count;
							spec_diff /= count;
							for (int i=0; i < XSIZE(F1); i++)
								std::cerr << " i= " << i << " spec_img(i)= "<< spec_img(i) << " spec_diff(i)= "<< spec_diff(i) << " sigma2_noise(i)= "<<  mymodel.sigma2_noise[group_id](i)
								<< " count(i)= " << count(i) << " sum-diff-norm=" << count(i)*spec_diff(i)<< " sum-diff-norm/sigma2= " << count(i)*(spec_diff(i)/mymodel.sigma2_noise[group_id](i))<< std::endl;
							Image<RFLOAT> tt;
							if (mymodel.data_dim == 3)
								tt().resize(YSIZE(F1), YSIZE(F1), YSIZE(F1));
							else
								tt().resize(YSIZE(F1), YSIZE(F1));
							transformer.inverseFourierTransform(F1, tt());
							CenterFFT(tt(),false);
							tt.write("F1.spi");
							transformer.inverseFourierTransform(F2, tt());
							CenterFFT(tt(),false);
							tt.write("F2.spi");
							std::cerr << "Written F1.spi and F2.spi. Press any key to continue... "<<std::endl;
							char c;
							std::cin >> c;
#endif

						// Only for the psi-angle and the translations, and only when my_prob < 0.01 calculate a histogram of the contributions at each resolution shell
						if (my_snr > pvalue && imode == 0)
						{
							FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F1)
							{
								int ires = DIRECT_MULTIDIM_ELEM(*myMresol, n);
								if (ires > 0)
									mymodel.orientability_contrib[iclass](ires) +=
											norm(DIRECT_MULTIDIM_ELEM(F1, n) - DIRECT_MULTIDIM_ELEM(F2, n)) / ( (2 * sigma2_fudge * mymodel.sigma2_noise[group_id](ires) ) );
							}
						}

					} // end while my_snr >= pvalue
					if (imode == 0)
						acc_rot_class += ang_error;
					else if (imode == 1)
						acc_trans_class += sh_error;
				} // end for imode

			}// end for part_id

			progress_bar(n_trials*iclass + my_metadata_entry);
		} // end for ori_part_id

		mymodel.acc_rot[iclass]   = acc_rot_class / (RFLOAT)n_trials;
		mymodel.acc_trans[iclass] = acc_trans_class / (RFLOAT)n_trials;

		// Store normalised spectral contributions to orientability
		if (mymodel.orientability_contrib[iclass].sum() > 0.)
			mymodel.orientability_contrib[iclass]   /= mymodel.orientability_contrib[iclass].sum();

		// Keep the orientational accuracy of the best class for the auto-sampling approach
		acc_rot     = XMIPP_MIN(mymodel.acc_rot[iclass], acc_rot);
		acc_trans   = XMIPP_MIN(mymodel.acc_trans[iclass], acc_trans);


		// Richard's formula with Greg's constant
		//RFLOAT b_orient = (acc_rot_class*acc_rot_class* particle_diameter*particle_diameter) / 3000.;
		//std::cout << " + expected B-factor from the orientational errors = "
		//		<< b_orient<<std::endl;
		// B=8 PI^2 U^2
		//std::cout << " + expected B-factor from the translational errors = "
		//		<< 8 * PI * PI * mymodel.pixel_size * mymodel.pixel_size * acc_trans_class * acc_trans_class << std::endl;

	} // end loop iclass
	progress_bar(n_trials * mymodel.nr_classes);


	std::cout << " Auto-refine: Estimated accuracy angles= " << acc_rot<< " degrees; offsets= " << acc_trans << " pixels" << std::endl;
	// Warn for inflated resolution estimates
	if (acc_rot > 10. && do_auto_refine)
	{
		std::cout << " Auto-refine: WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles (yet)!" << std::endl;
		std::cout << " Auto-refine: WARNING: You probably need not worry if the accuracy improves during the next few iterations." << std::endl;
		std::cout << " Auto-refine: WARNING: However, if the problem persists it may lead to spurious FSC curves, so be wary of inflated resolution estimates..." << std::endl;
		std::cout << " Auto-refine: WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
	}

}

void MlOptimiser::updateAngularSampling(bool myverb)
{

	if (!do_split_random_halves)
		REPORT_ERROR("MlOptimiser::updateAngularSampling: BUG! updating of angular sampling should only happen for gold-standard (auto-) refinements.");

	if (do_realign_movies)
	{

		// A. Adjust translational sampling to 75% of estimated accuracy
		RFLOAT new_step = XMIPP_MIN(1.5, 0.75 * acc_trans) * std::pow(2., adaptive_oversampling);

		// Search ranges are three times the estimates std.dev. in the offsets
		RFLOAT new_range = 3. * sqrt(mymodel.sigma2_offset);

		// Prevent too narrow searches: always at least 3x3 pixels in the coarse search
		if (new_range < 1.5 * new_step)
			new_range = 1.5 * new_step;

		// Also prevent too wide searches: that will lead to memory problems:
		// Just use coarser step size and hope things will settle down later...
		if (new_range > 4. * new_step)
			new_step = new_range / 4.;

		// Jun08,2015 Shaoda & Sjors, Helical refinement
		RFLOAT new_helical_offset_step = sampling.helical_offset_step;
		bool do_local_searches = true;
		if ( (do_helical_refine) && (!ignore_helical_symmetry) )
		{
			if (new_step < new_helical_offset_step)
				new_helical_offset_step = new_step;
		}
		sampling.setTranslations(new_step, new_range, do_local_searches, (do_helical_refine) && (!ignore_helical_symmetry), new_helical_offset_step, helical_rise_initial / mymodel.pixel_size, helical_twist_initial);

		if (!do_skip_rotate || do_only_sample_tilt)
		{
			// B. Find the healpix order that corresponds to at least 50% of the estimated rotational accuracy
			RFLOAT angle_range = sqrt(mymodel.sigma2_tilt) * 3.;
			RFLOAT new_ang_step, new_ang_step_wo_over;
			int new_hp_order;
			for (new_hp_order = 0; new_hp_order < 8; new_hp_order++)
			{
				new_ang_step = 360. / (6 * ROUND(std::pow(2., new_hp_order + adaptive_oversampling)));
				new_ang_step_wo_over = 2. * new_ang_step;
				// Only consider healpix orders that gives at least more than one (non-oversampled) samplings within the local angular searches
				if (new_ang_step_wo_over > angle_range)
					continue;
				// If sampling is at least twice as fine as the estimated rotational accuracy, then use this sampling
				if (new_ang_step < 0.50 * acc_rot)
					break;
			}

			if (new_hp_order != sampling.healpix_order)
			{
				// Set the new sampling in the sampling-object
				sampling.setOrientations(new_hp_order, new_ang_step * std::pow(2., adaptive_oversampling));
				// Resize the pdf_direction arrays to the correct size and fill with an even distribution
				mymodel.initialisePdfDirection(sampling.NrDirections());
				// Also reset the nr_directions in wsum_model
				wsum_model.nr_directions = mymodel.nr_directions;
				// Also resize and initialise wsum_model.pdf_direction for each class!
				for (int iclass=0; iclass < mymodel.nr_classes; iclass++)
					wsum_model.pdf_direction[iclass].initZeros(mymodel.nr_directions);
			}
		}
	}
	else
	{

		if (do_skip_rotate)
			REPORT_ERROR("ERROR: --skip_rotate can only be used in movie-frame refinement ...");

		// Only change the sampling if the resolution has not improved during the last 2 iterations
		// AND the hidden variables have not changed during the last 2 iterations
		RFLOAT old_rottilt_step = sampling.getAngularSampling(adaptive_oversampling);

		// Only use a finer angular sampling is the angular accuracy is still above 75% of the estimated accuracy
		// If it is already below, nothing will change and eventually nr_iter_wo_resol_gain or nr_iter_wo_large_hidden_variable_changes will go above MAX_NR_ITER_WO_RESOL_GAIN
		//std::cerr << "iter= " << iter << " nr_iter_wo_resol_gain= " << nr_iter_wo_resol_gain << " nr_iter_wo_large_hidden_variable_changes= " << nr_iter_wo_large_hidden_variable_changes << " acc_rot= " << acc_rot << std::endl;
		//std::cerr  << "iter= " << iter << " acc_trans= " << acc_trans << " minimum_angular_sampling= " << minimum_angular_sampling << " current_changes_optimal_offsets= " << current_changes_optimal_offsets << std::endl;
		if (nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN && nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES)
		{
			// Old rottilt step is already below 75% of estimated accuracy: have to stop refinement?
			// If a minimum_angular_sampling is given and we're not there yet, also just continue
			if (old_rottilt_step < 0.75 * acc_rot && !(minimum_angular_sampling > 0. && old_rottilt_step > minimum_angular_sampling) )
			{
				// don't change angular sampling, as it is already fine enough
				has_fine_enough_angular_sampling = true;
			}
			else
			{
				has_fine_enough_angular_sampling = false;

				// A. Use translational sampling as suggested by acc_trans

				// Prevent very coarse translational samplings: max 1.5
				// Also stay a bit on the safe side with the translational sampling: 75% of estimated accuracy
				RFLOAT new_step = XMIPP_MIN(1.5, 0.75 * acc_trans) * std::pow(2., adaptive_oversampling);
				// For subtomogram averaging: use at least half times previous step size
				if (mymodel.data_dim == 3) // TODO: check: this might just as well work for 2D data...
					new_step = XMIPP_MAX(sampling.offset_step / 2., new_step);
				// Search ranges are five times the last observed changes in offsets
				// Only 3x for subtomogram averaging....
				RFLOAT new_range = (mymodel.data_dim == 2) ? 5. * current_changes_optimal_offsets : 3 * current_changes_optimal_offsets;
				// New range can only become 30% bigger than the previous range (to prevent very slow iterations in the beginning)
				new_range = XMIPP_MIN(1.3*sampling.offset_range, new_range);
				// Prevent too narrow searches: always at least 3x3 pixels in the coarse search
				if (new_range < 1.5 * new_step)
					new_range = 1.5 * new_step;
				// Also prevent too wide searches: that will lead to memory problems:
				// If steps size < 1/4th of search range, then decrease search range by 50%
				if (new_range > 4. * new_step)
					new_range /= 2.;
				//If even that was not enough: use coarser step size and hope things will settle down later...
				if (new_range > 4. * new_step)
					new_step = new_range / 4.;
				// Jun08,2015 Shaoda & Sjors, Helical refinement
				RFLOAT new_helical_offset_step = sampling.helical_offset_step;
				if (mymodel.ref_dim == 3)
				{
					// Jun08,2015 Shaoda & Sjors, Helical refinement
					//new_helical_offset_step /= 2.;

					// AFTER AUG17,2015
					// ??? new_step ~= 1/4 * helical_offset_step ??? still divide helical_offset_step by 2 ! That is reasonable... Because it cannot happen (see above)
					if (new_step < new_helical_offset_step)
						new_helical_offset_step /= 2.;
				}

				// B. Use twice as fine angular sampling
				int new_hp_order;
				RFLOAT new_rottilt_step, new_psi_step;
				if (mymodel.ref_dim == 3)
				{
					new_hp_order = sampling.healpix_order + 1;
					new_rottilt_step = new_psi_step = 360. / (6 * ROUND(std::pow(2., new_hp_order + adaptive_oversampling)));

					// Set the new sampling in the sampling-object
					sampling.setOrientations(new_hp_order, new_psi_step * std::pow(2., adaptive_oversampling));

					// Resize the pdf_direction arrays to the correct size and fill with an even distribution
					mymodel.initialisePdfDirection(sampling.NrDirections());

					// Also reset the nr_directions in wsum_model
					wsum_model.nr_directions = mymodel.nr_directions;

					// Also resize and initialise wsum_model.pdf_direction for each class!
					for (int iclass=0; iclass < mymodel.nr_classes; iclass++)
						wsum_model.pdf_direction[iclass].initZeros(mymodel.nr_directions);

				}
				else if (mymodel.ref_dim == 2)
				{
					sampling.psi_step /= 2.;
				}
				else
					REPORT_ERROR("MlOptimiser::autoAdjustAngularSampling BUG: ref_dim should be two or three");

				// Jun08,2015 Shaoda & Sjors, Helical refinement
				bool do_local_searches = ((do_auto_refine) && (sampling.healpix_order >= autosampling_hporder_local_searches));

				// Don't go to coarse angular samplings. Then just keep doing as it was
				if (new_step > sampling.offset_step)
				{
					new_step = sampling.offset_step;
					new_range = sampling.offset_range;
				}

				sampling.setTranslations(new_step, new_range, do_local_searches, (do_helical_refine) && (!ignore_helical_symmetry), new_helical_offset_step, helical_rise_initial / mymodel.pixel_size, helical_twist_initial);

				// Reset iteration counters
				nr_iter_wo_resol_gain = 0;
				nr_iter_wo_large_hidden_variable_changes = 0;

				// Reset smallest changes hidden variables
				smallest_changes_optimal_classes = 9999999;
				smallest_changes_optimal_offsets = 999.;
				smallest_changes_optimal_orientations = 999.;

				// If the angular sampling is smaller than autosampling_hporder_local_searches, then use local searches of +/- 6 times the angular sampling
				if (mymodel.ref_dim == 3 && new_hp_order >= autosampling_hporder_local_searches)
				{
					// Switch ON local angular searches
					mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
					sampling.orientational_prior_mode = PRIOR_ROTTILT_PSI;
					mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 2. * 2. * new_rottilt_step * new_rottilt_step;

					// Aug20,2015 - Shaoda, Helical refinement
					if ( (do_helical_refine) && (!ignore_helical_symmetry) )
						mymodel.sigma2_rot = getHelicalSigma2Rot((helical_rise_initial / mymodel.pixel_size), helical_twist_initial, sampling.helical_offset_step, new_rottilt_step, mymodel.sigma2_rot);
				}

			}
		}
	}

	// Print to screen
	if (myverb)
	{
		std::cout << " Auto-refine: Angular step= " << sampling.getAngularSampling(adaptive_oversampling) << " degrees; local searches= ";
		if (sampling.orientational_prior_mode == NOPRIOR)
			std:: cout << "false" << std::endl;
		else
			std:: cout << "true" << std::endl;
		// Jun08,2015 Shaoda & Sjors, Helical refine
		if ( (do_helical_refine) && (!ignore_helical_symmetry) )
		{
			std::cout << " Auto-refine: Helical refinement... Local translational searches along helical axis= ";
			if ( (mymodel.ref_dim == 3) && (do_auto_refine) && (sampling.healpix_order >= autosampling_hporder_local_searches) )
				std:: cout << "true" << std::endl;
			else
				std:: cout << "false" << std::endl;
		}
		std::cout << " Auto-refine: Offset search range= " << sampling.offset_range << " pixels; offset step= " << sampling.getTranslationalSampling(adaptive_oversampling) << " pixels";
		if ( (do_helical_refine) && (!ignore_helical_symmetry) )
			std::cout << "; offset step along helical axis= " << sampling.getHelicalTranslationalSampling(adaptive_oversampling) << " pixels";
		std::cout << std::endl;
	}

}

void MlOptimiser::checkConvergence(bool myverb)
{

	if (do_realign_movies)
		return;

	if ( has_fine_enough_angular_sampling && nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN && nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES )
	{
		has_converged = true;
		do_join_random_halves = true;
		// In the last iteration, include all data until Nyquist
		do_use_all_data = true;
	}

        if (myverb)
        {
            std::cout << " Auto-refine: Iteration= "<< iter<< std::endl;
            std::cout << " Auto-refine: Resolution= "<< 1./mymodel.current_resolution<< " (no gain for " << nr_iter_wo_resol_gain << " iter) "<< std::endl;
            std::cout << " Auto-refine: Changes in angles= " << current_changes_optimal_orientations << " degrees; and in offsets= " << current_changes_optimal_offsets
			<< " pixels (no gain for " << nr_iter_wo_large_hidden_variable_changes << " iter) "<< std::endl;

            if (has_converged)
            {
		std::cout << " Auto-refine: Refinement has converged, entering last iteration where two halves will be combined..."<<std::endl;
		std::cout << " Auto-refine: The last iteration will use data to Nyquist frequency, which may take more CPU and RAM."<<std::endl;
            }
        }

}

void MlOptimiser::setMetaDataSubset(int first_ori_particle_id, int last_ori_particle_id)
{

	for (long int ori_part_id = first_ori_particle_id, my_image_no = 0; ori_part_id <= last_ori_particle_id; ori_part_id++)
    {

#ifdef DEBUG_CHECKSIZES
		if (ori_part_id >= mydata.ori_particles.size())
		{
			std::cerr<< "ori_part_id= "<<ori_part_id<<" mydata.ori_particles.size()= "<< mydata.ori_particles.size() <<std::endl;
			REPORT_ERROR("ori_part_id >= mydata.ori_particles.size()");
		}
#endif

		for (long int ipart = 0; ipart < mydata.ori_particles[ori_part_id].particles_id.size(); ipart++, my_image_no++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[ipart];

#ifdef DEBUG_CHECKSIZES
			if (part_id >= mydata.MDimg.numberOfObjects())
			{
				std::cerr<< "part_id= "<<part_id<<" mydata.MDimg.numberOfObjects()= "<< mydata.MDimg.numberOfObjects() <<std::endl;
				REPORT_ERROR("part_id >= mydata.MDimg.numberOfObjects()");
			}
			if (my_image_no >= YSIZE(exp_metadata))
			{
				std::cerr<< "my_image_no= "<<my_image_no<<" YSIZE(exp_metadata)= "<< YSIZE(exp_metadata) <<std::endl;
				REPORT_ERROR("my_image_no >= YSIZE(exp_metadata)");
			}
#endif
			// SHWS: Upon request of Juha Huiskonen, 5apr2016
			if (mymodel.ref_dim > 2)
			{
				mydata.MDimg.setValue(EMDL_ORIENT_ROT,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT), part_id);
				mydata.MDimg.setValue(EMDL_ORIENT_TILT, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT), part_id);
			}
			mydata.MDimg.setValue(EMDL_ORIENT_PSI,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI), part_id);
			mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_X, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF), part_id);
			mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF), part_id);
			if (mymodel.data_dim == 3)
				mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_Z, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ZOFF), part_id);
			mydata.MDimg.setValue(EMDL_PARTICLE_CLASS, (int)DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CLASS) , part_id);
			mydata.MDimg.setValue(EMDL_PARTICLE_DLL,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_DLL), part_id);
			mydata.MDimg.setValue(EMDL_PARTICLE_PMAX, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PMAX), part_id);
			mydata.MDimg.setValue(EMDL_PARTICLE_NR_SIGNIFICANT_SAMPLES,(int)DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NR_SIGN), part_id);
			mydata.MDimg.setValue(EMDL_IMAGE_NORM_CORRECTION, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM), part_id);

			// For the moment, CTF, prior and transformation matrix info is NOT updated...
			RFLOAT prior_x = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF_PRIOR);
			RFLOAT prior_y = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF_PRIOR);
			if (prior_x < 999.)
				mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_X_PRIOR, prior_x, part_id);
			if (prior_y < 999.)
				mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, prior_y, part_id);
			if (mymodel.data_dim == 3)
			{
				RFLOAT prior_z = DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ZOFF_PRIOR);
				if (prior_z < 999.)
					mydata.MDimg.setValue(EMDL_ORIENT_ORIGIN_Z_PRIOR, prior_z, part_id);
			}

			// For multi-body refinement
			if (mymodel.nr_bodies > 1)
			{
				for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++)
				{
					int icol_rot  = 0 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_tilt = 1 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_psi  = 2 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_xoff = 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_yoff = 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_zoff = 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_norm = 6 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					RFLOAT rot =  DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_rot);
					RFLOAT tilt = DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_tilt);
					RFLOAT psi =  DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_psi);
					RFLOAT xoff = DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_xoff);
					RFLOAT yoff = DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_yoff);
					RFLOAT zoff = DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_zoff);
					RFLOAT norm = DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_norm);
					mydata.MDbodies[ibody].setValue(EMDL_ORIENT_ROT, rot, part_id);
					mydata.MDbodies[ibody].setValue(EMDL_ORIENT_TILT, tilt, part_id);
					mydata.MDbodies[ibody].setValue(EMDL_ORIENT_PSI,  psi, part_id);
					mydata.MDbodies[ibody].setValue(EMDL_ORIENT_ORIGIN_X, xoff, part_id);
					mydata.MDbodies[ibody].setValue(EMDL_ORIENT_ORIGIN_Y, yoff, part_id);
					if (mymodel.data_dim == 3)
						mydata.MDbodies[ibody].setValue(EMDL_ORIENT_ORIGIN_Z, zoff, part_id);
					mydata.MDbodies[ibody].setValue(EMDL_IMAGE_NORM_CORRECTION, norm, part_id);
				}
			}

		}
	}

}

void MlOptimiser::getMetaAndImageDataSubset(int first_ori_particle_id, int last_ori_particle_id, bool do_also_imagedata)
{

    // In case we're reading images here, only open stacks once and then read multiple images
	fImageHandler hFile;
	FileName fn_img, fn_stack, fn_open_stack="";;
	long int dump;

    // Initialise filename strings if not reading imagedata here
	if (!do_also_imagedata)
	{
		exp_fn_img = "";
		exp_fn_ctf = "";
		exp_fn_recimg = "";
	}

	int nr_images = 0;
	for (long int ori_part_id = first_ori_particle_id; ori_part_id <= last_ori_particle_id; ori_part_id++)
	{
		nr_images += mydata.ori_particles[ori_part_id].particles_id.size();
	}

	exp_metadata.initZeros(nr_images, METADATA_LINE_LENGTH_BEFORE_BODIES + (mymodel.nr_bodies) * METADATA_NR_BODY_PARAMS);
	if (do_also_imagedata)
	{
		if (mymodel.data_dim == 3)
		{
			if (nr_images > 1)
				REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: cannot get multiple images for 3D data!");

			if (do_ctf_correction)
			{
				if (has_converged && do_use_reconstruct_images)
					exp_imagedata.resize(3*mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
				else
					exp_imagedata.resize(2*mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
			}
			else
			{
				if (has_converged && do_use_reconstruct_images)
					exp_imagedata.resize(2*mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
				else
					exp_imagedata.resize(mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
			}
		}
		else
		{
			if (has_converged && do_use_reconstruct_images)
				exp_imagedata.resize(2*nr_images, mymodel.ori_size, mymodel.ori_size);
			else
				exp_imagedata.resize(nr_images, mymodel.ori_size, mymodel.ori_size);
		}
	}

	for (long int ori_part_id = first_ori_particle_id, my_image_no = 0; ori_part_id <= last_ori_particle_id; ori_part_id++)
    {
		for (long int ipart = 0; ipart < mydata.ori_particles[ori_part_id].particles_id.size(); ipart++, my_image_no++)
		{
			long int part_id = mydata.ori_particles[ori_part_id].particles_id[ipart];

#ifdef DEBUG_CHECKSIZES
			if (part_id >= mydata.MDimg.numberOfObjects())
			{
				std::cerr<< "part_id= "<<part_id<<" mydata.MDimg.numberOfObjects()= "<< mydata.MDimg.numberOfObjects() <<std::endl;
				REPORT_ERROR("part_id >= mydata.MDimg.numberOfObjects()");
			}
			if (my_image_no >= YSIZE(exp_metadata))
			{
				std::cerr<< "my_image_no= "<<my_image_no<<" YSIZE(exp_metadata)= "<< YSIZE(exp_metadata) <<std::endl;
				REPORT_ERROR("my_image_no >= YSIZE(exp_metadata)");
			}
			if (my_image_no >= nr_images)
			{
				std::cerr<< "my_image_no= "<<my_image_no<<" nr_images= "<< nr_images <<std::endl;
				REPORT_ERROR("my_image_no >= nr_images");
			}
#endif
			// Get the image names from the MDimg table
			FileName fn_img="", fn_rec_img="", fn_ctf="";
			if (!mydata.getImageNameOnScratch(part_id, fn_img))
				mydata.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
			if (mymodel.data_dim == 3 && do_ctf_correction)
			{
				// Also read the CTF image from disc
				if (!mydata.getImageNameOnScratch(part_id, fn_ctf, true))
				{
					if (!mydata.MDimg.getValue(EMDL_CTF_IMAGE, fn_ctf, part_id))
						REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: cannot find rlnCtfImage for 3D CTF correction!");
				}
			}
			if (has_converged && do_use_reconstruct_images)
			{
				mydata.MDimg.getValue(EMDL_IMAGE_RECONSTRUCT_NAME, fn_rec_img, part_id);
			}

			if (do_also_imagedata)
			{
				// First read the image from disc or get it from the preread images in the mydata structure
				Image<RFLOAT> img, rec_img;
				if (do_preread_images)
				{
					img().reshape(mydata.particles[part_id].img);
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(mydata.particles[part_id].img)
					{
						DIRECT_MULTIDIM_ELEM(img(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(mydata.particles[part_id].img, n);
					}
				}
				else
				{
					// only open new stacks
					fn_img.decompose(dump, fn_stack);
					if (fn_stack != fn_open_stack)
					{
						hFile.openFile(fn_stack, WRITE_READONLY);
						fn_open_stack = fn_stack;
					}
					img.readFromOpenFile(fn_img, hFile, -1, false);
					img().setXmippOrigin();
				}
				if (XSIZE(img()) != XSIZE(exp_imagedata) || YSIZE(img()) != YSIZE(exp_imagedata) )
				{
					std::cerr << " fn_img= " << fn_img << " XSIZE(img())= " << XSIZE(img()) << " YSIZE(img())= " << YSIZE(img()) << std::endl;
					std::cerr << " while XSIZE(exp_imagedata)= " << XSIZE(exp_imagedata) << " and YSIZE(exp_imagedata)= " << YSIZE(exp_imagedata) << std::endl;
					REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: incorrect image size");
				}
				if (has_converged && do_use_reconstruct_images)
				{
					rec_img.read(fn_rec_img);
					if (XSIZE(rec_img()) != XSIZE(exp_imagedata) || YSIZE(rec_img()) != YSIZE(exp_imagedata) )
					{
						std::cerr << " fn_rec_img= " << fn_rec_img << " XSIZE(rec_img())= " << XSIZE(rec_img()) << " YSIZE(rec_img())= " << YSIZE(rec_img()) << std::endl;
						REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: incorrect reconstruct_image size");
					}
				}
				if (mymodel.data_dim == 3)
				{

					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
					{
						DIRECT_A3D_ELEM(exp_imagedata, k, i, j) = DIRECT_A3D_ELEM(img(), k, i, j);
					}

					if (do_ctf_correction)
					{
						img.read(fn_ctf);
						FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
						{
							DIRECT_A3D_ELEM(exp_imagedata, mymodel.ori_size + k, i, j) = DIRECT_A3D_ELEM(img(), k, i, j);
						}
					}

					if (has_converged && do_use_reconstruct_images)
					{
						int offset = (do_ctf_correction) ? 2 * mymodel.ori_size : mymodel.ori_size;
						FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
						{
							DIRECT_A3D_ELEM(exp_imagedata, offset + k, i, j) = DIRECT_A3D_ELEM(rec_img(), k, i, j);
						}
					}

				}
				else
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img())
					{
						DIRECT_A3D_ELEM(exp_imagedata, my_image_no, i, j) = DIRECT_A2D_ELEM(img(), i, j);
					}

					if (has_converged && do_use_reconstruct_images)
					{
						FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(rec_img())
						{
							DIRECT_A3D_ELEM(exp_imagedata, nr_images + my_image_no, i, j) = DIRECT_A2D_ELEM(rec_img(), i, j);
						}
					}

				}
			}
			else
			{
				exp_fn_img += fn_img + "\n";
				if (fn_ctf != "")
					exp_fn_ctf += fn_ctf + "\n";
				if (fn_rec_img != "")
					exp_fn_recimg += fn_rec_img + "\n";
			}

			// Now get the metadata
			int iaux;
			mydata.MDimg.getValue(EMDL_ORIENT_ROT,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT), part_id);
			mydata.MDimg.getValue(EMDL_ORIENT_TILT, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT), part_id);
			mydata.MDimg.getValue(EMDL_ORIENT_PSI,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI), part_id);
			mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF), part_id);
			mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF), part_id);
			if (mymodel.data_dim == 3)
				mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Z, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ZOFF), part_id);
			mydata.MDimg.getValue(EMDL_PARTICLE_CLASS, iaux, part_id);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CLASS) = (RFLOAT)iaux;
			mydata.MDimg.getValue(EMDL_PARTICLE_DLL,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_DLL), part_id);
			mydata.MDimg.getValue(EMDL_PARTICLE_PMAX, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PMAX), part_id);
			mydata.MDimg.getValue(EMDL_PARTICLE_NR_SIGNIFICANT_SAMPLES, iaux, part_id);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NR_SIGN) = (RFLOAT)iaux;
			if (!mydata.MDimg.getValue(EMDL_IMAGE_NORM_CORRECTION, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM), part_id))
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_NORM) = 1.;
			if (do_ctf_correction)
			{
				long int mic_id = mydata.getMicrographId(part_id);
				RFLOAT kV, DeltafU, DeltafV, azimuthal_angle, Cs, Bfac, Q0, phase_shift;
				if (!mydata.MDimg.getValue(EMDL_CTF_VOLTAGE, kV, part_id))
					if (!mydata.MDmic.getValue(EMDL_CTF_VOLTAGE, kV, mic_id))
						kV=200;

				if (!mydata.MDimg.getValue(EMDL_CTF_DEFOCUSU, DeltafU, part_id))
					if (!mydata.MDmic.getValue(EMDL_CTF_DEFOCUSU, DeltafU, mic_id))
						DeltafU=0;

				if (!mydata.MDimg.getValue(EMDL_CTF_DEFOCUSV, DeltafV, part_id))
					if (!mydata.MDmic.getValue(EMDL_CTF_DEFOCUSV, DeltafV, mic_id))
						DeltafV=DeltafU;

				if (!mydata.MDimg.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, part_id))
					if (!mydata.MDmic.getValue(EMDL_CTF_DEFOCUS_ANGLE, azimuthal_angle, mic_id))
						azimuthal_angle=0;

				if (!mydata.MDimg.getValue(EMDL_CTF_CS, Cs, part_id))
					if (!mydata.MDmic.getValue(EMDL_CTF_CS, Cs, mic_id))
						Cs=0;

				if (!mydata.MDimg.getValue(EMDL_CTF_BFACTOR, Bfac, part_id))
					if (!mydata.MDmic.getValue(EMDL_CTF_BFACTOR, Bfac, mic_id))
						Bfac=0;

				if (!mydata.MDimg.getValue(EMDL_CTF_Q0, Q0, part_id))
					if (!mydata.MDmic.getValue(EMDL_CTF_Q0, Q0, mic_id))
						Q0=0;

				if (!mydata.MDimg.getValue(EMDL_CTF_PHASESHIFT, phase_shift, part_id))
					if (!mydata.MDmic.getValue(EMDL_CTF_PHASESHIFT, phase_shift, mic_id))
						phase_shift=0;

				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_VOLTAGE) = kV;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_U) = DeltafU;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_V) = DeltafV;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_DEFOCUS_ANGLE) = azimuthal_angle;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_CS) = Cs;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_BFAC) = Bfac;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_Q0) = Q0;
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_CTF_PHASE_SHIFT) = phase_shift;

			}

			// beamtilt
			RFLOAT beamtilt_x = 0., beamtilt_y = 0.;
			if (mydata.MDimg.containsLabel(EMDL_IMAGE_BEAMTILT_X))
				mydata.MDimg.getValue(EMDL_IMAGE_BEAMTILT_X, beamtilt_x, part_id);
			if (mydata.MDimg.containsLabel(EMDL_IMAGE_BEAMTILT_Y))
				mydata.MDimg.getValue(EMDL_IMAGE_BEAMTILT_Y, beamtilt_y, part_id);
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_BEAMTILT_X) = beamtilt_x;
			DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_BEAMTILT_Y) = beamtilt_y;

			// If the priors are NOT set, then set their values to 999.
			if (!mydata.MDimg.getValue(EMDL_ORIENT_ROT_PRIOR,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT_PRIOR), part_id))
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ROT_PRIOR) = 999.;
			if (!mydata.MDimg.getValue(EMDL_ORIENT_TILT_PRIOR, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT_PRIOR), part_id))
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_TILT_PRIOR) = 999.;
			if (!mydata.MDimg.getValue(EMDL_ORIENT_PSI_PRIOR,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI_PRIOR), part_id))
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI_PRIOR) = 999.;
			if (!mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_X_PRIOR, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF_PRIOR), part_id))
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_XOFF_PRIOR) = 999.;
			if (!mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Y_PRIOR, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF_PRIOR), part_id))
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_YOFF_PRIOR) = 999.;
			if (mymodel.data_dim == 3)
				if (!mydata.MDimg.getValue(EMDL_ORIENT_ORIGIN_Z_PRIOR, DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ZOFF_PRIOR), part_id))
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_ZOFF_PRIOR) = 999.;
			if (!mydata.MDimg.getValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO,  DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI_PRIOR_FLIP_RATIO), part_id))
				DIRECT_A2D_ELEM(exp_metadata, my_image_no, METADATA_PSI_PRIOR_FLIP_RATIO) = 999.;

			// For multi-body refinement
			if (mymodel.nr_bodies > 1)
			{
				for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++)
				{
					int icol_rot  = 0 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_tilt = 1 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_psi  = 2 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_xoff = 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_yoff = 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_zoff = 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					int icol_norm = 6 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
					RFLOAT rot, tilt, psi, xoff, yoff, zoff, norm;
					mydata.MDbodies[ibody].getValue(EMDL_ORIENT_ROT, rot, part_id);
					mydata.MDbodies[ibody].getValue(EMDL_ORIENT_TILT, tilt, part_id);
					mydata.MDbodies[ibody].getValue(EMDL_ORIENT_PSI,  psi, part_id);
					mydata.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_X, xoff, part_id);
					mydata.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_Y, yoff, part_id);
					if (mymodel.data_dim == 3)
						mydata.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_Z, zoff, part_id);
					mydata.MDbodies[ibody].getValue(EMDL_IMAGE_NORM_CORRECTION, norm, part_id);
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_rot)  = rot;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_tilt) = tilt;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_psi)  = psi;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_xoff) = xoff;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_yoff) = yoff;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_zoff) = zoff;
					DIRECT_A2D_ELEM(exp_metadata, my_image_no, icol_norm) = norm;
				}
			}
		}
    }

}


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
#include "src/ml_optimiser_mpi.h"
#include "src/ml_optimiser.h"
#include "src/gpu_utils/cuda_ml_optimiser.h"

//#define DEBUG
//#define DEBUG_MPIEXP2

#ifdef TIMING
        int TIMING_MPIPACK, TIMING_MPIWAIT, TIMING_MPICOMBINEDISC, TIMING_MPICOMBINENETW, TIMING_MPISLAVEWORK;
        int TIMING_MPISLAVEWAIT1, TIMING_MPISLAVEWAIT2, TIMING_MPISLAVEWAIT3;
#endif

void MlOptimiserMpi::read(int argc, char **argv)
{
#ifdef DEBUG
    std::cerr<<"MlOptimiserMpi::read Entering "<<std::endl;
#endif

    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    MlOptimiser::read(argc, argv, node->rank);
    fn_scratch = parser.getOption("--scratchdir", "Directory (with absolute path, and visible to all nodes) for temporary files", "");

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? 1 : 0;

//#define DEBUG_BODIES
#ifdef DEBUG_BODIES
    verb = (node->rank==1) ? 1 : 0;
#endif

    // TMP for debugging only
    //if (node->rank==1)
    //	verb = 1;
    // Possibly also read parallelisation-dependent variables here

#ifdef DEBUG
    std::cerr<<"MlOptimiserMpi::read done"<<std::endl;
#endif

}
void MlOptimiserMpi::finalise()
{
	delete node;
}

void MlOptimiserMpi::initialise()
{

#ifdef DEBUG
    std::cerr<<"MlOptimiserMpi::initialise Entering"<<std::endl;
#endif

    // Print information about MPI nodes:
    printMpiNodesMachineNames(*node, nr_threads);

    MlOptimiser::initialiseGeneral(node->rank);

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
		// This is done in parallel
		//std::cout << " Hello world1! I am node " << node->rank << " out of " << node->size <<" and my hostname= "<< getenv("HOSTNAME")<< std::endl;
		calculateSumOfPowerSpectraAndAverageImage(Mavg);

		// Set sigma2_noise and Iref from averaged poser spectra and Mavg
		if (!node->isMaster())
			MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage(Mavg);
		//std::cout << " Hello world3! I am node " << node->rank << " out of " << node->size <<" and my hostname= "<< getenv("HOSTNAME")<< std::endl;
	}

   MlOptimiser::initialLowPassFilterReferences();

	// Initialise the data_versus_prior ratio to get the initial current_size right
	if (iter == 0)
		mymodel.initialiseDataVersusPrior(fix_tau); // fix_tau was set in initialiseGeneral

	//std::cout << " Hello world! I am node " << node->rank << " out of " << node->size <<" and my hostname= "<< getenv("HOSTNAME")<< std::endl;

	// Only master writes out initial mymodel (do not gather metadata yet)
	int nr_subsets = (do_split_random_halves) ? 2 : 1;
	if (node->isMaster())
		MlOptimiser::write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, node->rank);
	else if (node->rank <= nr_subsets)
	{
		//Only the first_slave of each subset writes model to disc
		MlOptimiser::write(DO_WRITE_SAMPLING, DONT_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, node->rank);

		bool do_warn = false;
		for (int igroup = 0; igroup< mymodel.nr_groups; igroup++)
		{
			if (mymodel.nr_particles_group[igroup] < 5 && node->rank == 1) // only warn for half1 to avoid messy output
			{
				if (nr_subsets == 1)
					std:: cout << "WARNING: There are only " << mymodel.nr_particles_group[igroup] << " particles in group " << igroup + 1 << std::endl;
				else
					std:: cout << "WARNING: There are only " << mymodel.nr_particles_group[igroup] << " particles in group " << igroup + 1 << " of half-set " << node->rank << std::endl;
				do_warn = true;
			}
		}
		if (do_warn)
		{
			std:: cout << "WARNING: You may want to consider joining some micrographs into larger groups to obtain more robust noise estimates. " << std::endl;
			std:: cout << "         You can do so by using the same rlnMicrographName for particles from multiple different micrographs in the input STAR file. " << std::endl;
            std:: cout << "         It is then best to join micrographs with similar defocus values and similar apparent signal-to-noise ratios. " << std::endl;
		}
	}

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


	if (do_gpu)
	{
		do_shifts_onthefly = true;
		// Sequential initialisation of GPUs on all ranks
		for (int rank = 0; rank < node->size; rank++)
		{
			if (rank == node->rank && !node->isMaster()) //Device not initialized on master rank
			{
				if (!std::isdigit(*gpu_ids.begin()))
					std::cout << " No gpu-ids specified, threads will automatically be mapped to devices (incrementally)."<< std::endl;
				else if(gpu_ids.length()<nr_threads)
					REPORT_ERROR("You did not supply enough gpu ids to supply all the threads you wanted");

				int devCount;
				HANDLE_ERROR(cudaGetDeviceCount(&devCount));

				for (int i = 0; i < nr_threads; i ++)
				{
					int dev_id;
					if (!std::isdigit(*gpu_ids.begin()))
					{
						// Sjors: hack to make use of several cards; will only work if all MPI slaves are on the same node!
						if (node->size > 1)
							dev_id = ( (i * (node->size - 1)) + (node->rank - 1) )%devCount;
						else
							dev_id = i%devCount;
					}
					else
						dev_id = (int)(gpu_ids[i]-'0');

					std::cout << " Thread " << i << " on slave " << node->rank << " mapped to device " << dev_id << std::endl;
					cudaMlOptimisers.push_back((void *) new MlOptimiserCuda(this, dev_id));
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}


#ifdef DEBUG
    std::cerr<<"MlOptimiserMpi::initialise Done"<<std::endl;
#endif
}

void MlOptimiserMpi::initialiseWorkLoad()
{

    if (do_split_random_halves && node->size <= 2)
    	REPORT_ERROR("MlOptimiserMpi::initialiseWorkLoad: at least 3 MPI processes are required when splitting data into random halves");
    else if(node->size <= 1)
    	REPORT_ERROR("MlOptimiserMpi::initialiseWorkLoad: at least 2 MPI processes are required, otherwise use the sequential program");

	// Get the same random number generator seed for all mpi nodes
	if (random_seed == -1)
	{
		if (node->isMaster())
		{
			random_seed = time(NULL);
	        for (int slave = 1; slave < node->size; slave++)
	        	node->relion_MPI_Send(&random_seed, 1, MPI_INT, slave, MPITAG_RANDOMSEED, MPI_COMM_WORLD);
		}
		else
		{
			MPI_Status status;
			node->relion_MPI_Recv(&random_seed, 1, MPI_INT, 0, MPITAG_RANDOMSEED, MPI_COMM_WORLD, status);
		}
	}

    // First split the data into two random halves and then randomise the particle order
	if (do_split_random_halves)
		mydata.divideOriginalParticlesInRandomHalves(random_seed);

	// Randomise the order of the particles
    mydata.randomiseOriginalParticlesOrder(random_seed, do_split_random_halves);

    // Also randomize random-number-generator for perturbations on the angles
    init_random_generator(random_seed);


	if (node->isMaster())
	{
		// The master never participates in any actual work
		my_first_ori_particle_id = 0;
		my_last_ori_particle_id = -1;
	}
	else
	{
		if (do_split_random_halves)
		{
	    	int nr_slaves_subset1 = (node->size - 1) / 2;
	    	int nr_slaves_subset2 = nr_slaves_subset1;
	    	if ( (node->size - 1) % 2 != 0)
	    		nr_slaves_subset1 += 1;
	    	if (node->myRandomSubset() == 1)
	    	{
	    		// Divide first half of the images
	    		divide_equally(mydata.numberOfOriginalParticles(1), nr_slaves_subset1, node->rank / 2, my_first_ori_particle_id, my_last_ori_particle_id);
	    	}
	    	else
	    	{
	    		// Divide second half of the images
	    		divide_equally(mydata.numberOfOriginalParticles(2), nr_slaves_subset2, node->rank / 2 - 1, my_first_ori_particle_id, my_last_ori_particle_id);
	    		my_first_ori_particle_id += mydata.numberOfOriginalParticles(1);
	    		my_last_ori_particle_id += mydata.numberOfOriginalParticles(1);
	    	}
		}
		else
		{
			int nr_slaves = (node->size - 1);
			divide_equally(mydata.numberOfOriginalParticles(), nr_slaves, node->rank - 1, my_first_ori_particle_id, my_last_ori_particle_id);
		}

	}
//#define DEBUG_WORKLOAD
#ifdef DEBUG_WORKLOAD
	std::cerr << " node->rank= " << node->rank << " my_first_ori_particle_id= " << my_first_ori_particle_id << " my_last_ori_particle_id= " << my_last_ori_particle_id << std::endl;
#endif
}

void MlOptimiserMpi::calculateSumOfPowerSpectraAndAverageImage(MultidimArray<RFLOAT> &Mavg)
{

	// First calculate the sum of all individual power spectra on each subset
	MlOptimiser::calculateSumOfPowerSpectraAndAverageImage(Mavg, node->rank == 1);
	//std::cout << " Hello world22! I am node " << node->rank << " out of " << node->size <<" and my hostname= "<< getenv("HOSTNAME")<< std::endl;

	// When splitting the data into two random halves, perform two passes: one for each subset
	int nr_subsets = (do_split_random_halves) ? 2 : 1;

	MultidimArray<RFLOAT> Msum, MsumI;
	MultidimArray<int> Mnr, Msumnr;
	RFLOAT dsum;
	int isum;
	MPI_Status status;


	for (int isubset = 1; isubset <= nr_subsets; isubset++)
	{
		int my_first_slave = isubset; // first pass subset1: my_first_rank = 1, second pass subset2: my_first_rank = 2
		int my_first_other_slave = XMIPP_MIN(my_first_slave + nr_subsets, node->size);
		//std::cerr << " my_first_slave= " << my_first_slave << " my_first_other_slave= " << my_first_other_slave << std::endl;
		// Initialise on the first slave
		if (node->rank == my_first_slave)
		{
			Msum = wsum_model.sigma2_noise[0];
			MsumI = Mavg;
			dsum = wsum_model.sumw_group[0];
			Mnr.resize(mymodel.nr_groups);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Mnr)
			{
				DIRECT_A1D_ELEM(Mnr, i) = mymodel.nr_particles_group[i];
			}
			Msumnr = Mnr;

    	}
		for (int other_slave = my_first_other_slave; other_slave < node->size; other_slave += nr_subsets)
		{
			if (node->rank == other_slave)
			{
				//std::cerr << "Sending from "<<other_slave<< " to "<<my_first_slave << std::endl;
				node->relion_MPI_Send(MULTIDIM_ARRAY(wsum_model.sigma2_noise[0]), MULTIDIM_SIZE(wsum_model.sigma2_noise[0]), MY_MPI_DOUBLE, my_first_slave, MPITAG_PACK, MPI_COMM_WORLD);
				node->relion_MPI_Send(MULTIDIM_ARRAY(Mavg), MULTIDIM_SIZE(Mavg), MY_MPI_DOUBLE, my_first_slave, MPITAG_IMAGE, MPI_COMM_WORLD);
				Mnr.resize(mymodel.nr_groups);
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Mnr)
				{
					DIRECT_A1D_ELEM(Mnr, i) = mymodel.nr_particles_group[i];
        		}
				node->relion_MPI_Send(MULTIDIM_ARRAY(Mnr), MULTIDIM_SIZE(Mnr), MPI_INT, my_first_slave, MPITAG_METADATA, MPI_COMM_WORLD);
				node->relion_MPI_Send(&wsum_model.sumw_group[0], 1, MY_MPI_DOUBLE, my_first_slave, MPITAG_RFLOAT, MPI_COMM_WORLD);

    		}
			else if (node->rank == my_first_slave)
			{
				//std::cerr << "Receiving at "<<my_first_slave<< " from "<<other_slave<<std::endl;
				node->relion_MPI_Recv(MULTIDIM_ARRAY(wsum_model.sigma2_noise[0]), MULTIDIM_SIZE(wsum_model.sigma2_noise[0]), MY_MPI_DOUBLE, other_slave, MPITAG_PACK, MPI_COMM_WORLD, status);
				node->relion_MPI_Recv(MULTIDIM_ARRAY(Mavg), MULTIDIM_SIZE(Mavg), MY_MPI_DOUBLE, other_slave, MPITAG_IMAGE, MPI_COMM_WORLD, status);
				node->relion_MPI_Recv(MULTIDIM_ARRAY(Mnr), MULTIDIM_SIZE(Mnr), MPI_INT, other_slave, MPITAG_METADATA, MPI_COMM_WORLD, status);
				node->relion_MPI_Recv(&wsum_model.sumw_group[0], 1, MY_MPI_DOUBLE, other_slave, MPITAG_RFLOAT, MPI_COMM_WORLD, status);
				// Add the contribution of this other slave
				Msum  += wsum_model.sigma2_noise[0];
				MsumI += Mavg;
				Msumnr += Mnr;
				dsum += wsum_model.sumw_group[0];
    		}
		}

		// Now the first slave has the sum of all other slaves in its Msum
		// Send this sum to all relevant other slaves
		for (int other_slave = my_first_other_slave; other_slave < node->size; other_slave += nr_subsets)
		{
			if (node->rank == my_first_slave)
			{
				//std::cerr << "Sending from my_first_slave"<<my_first_slave<< " to other_slave "<<other_slave<< std::endl;
				node->relion_MPI_Send(MULTIDIM_ARRAY(Msum), MULTIDIM_SIZE(Msum), MY_MPI_DOUBLE, other_slave, MPITAG_PACK, MPI_COMM_WORLD);
				node->relion_MPI_Send(MULTIDIM_ARRAY(Msumnr), MULTIDIM_SIZE(Msumnr), MPI_INT, other_slave, MPITAG_METADATA, MPI_COMM_WORLD);
				node->relion_MPI_Send(MULTIDIM_ARRAY(MsumI), MULTIDIM_SIZE(MsumI), MY_MPI_DOUBLE, other_slave, MPITAG_IMAGE, MPI_COMM_WORLD);
				node->relion_MPI_Send(&dsum, 1, MY_MPI_DOUBLE, other_slave, MPITAG_RFLOAT, MPI_COMM_WORLD);
			}
			else if (node->rank == other_slave)
			{
				//std::cerr << "Receiving at other_slave "<<other_slave<<" from my_first_slave "<<my_first_slave<< std::endl;
				node->relion_MPI_Recv(MULTIDIM_ARRAY(wsum_model.sigma2_noise[0]), MULTIDIM_SIZE(wsum_model.sigma2_noise[0]), MY_MPI_DOUBLE, my_first_slave, MPITAG_PACK, MPI_COMM_WORLD, status);
				node->relion_MPI_Recv(MULTIDIM_ARRAY(Mnr), MULTIDIM_SIZE(Mnr), MPI_INT, my_first_slave, MPITAG_METADATA, MPI_COMM_WORLD, status);
				node->relion_MPI_Recv(MULTIDIM_ARRAY(Mavg), MULTIDIM_SIZE(Mavg), MY_MPI_DOUBLE, my_first_slave, MPITAG_IMAGE, MPI_COMM_WORLD, status);
				node->relion_MPI_Recv(&wsum_model.sumw_group[0], 1, MY_MPI_DOUBLE, my_first_slave, MPITAG_RFLOAT, MPI_COMM_WORLD, status);
				// Unpack Mnr on the other slave
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Mnr)
				{
					mymodel.nr_particles_group[i] = DIRECT_A1D_ELEM(Mnr, i);
				}
			}
		}

		// Also set sums on the first slave
		if (node->rank == my_first_slave)
		{
			wsum_model.sigma2_noise[0] = Msum;
			Mavg = MsumI;
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(Msumnr)
			{
				mymodel.nr_particles_group[i] = DIRECT_A1D_ELEM(Msumnr, i);
			}
			wsum_model.sumw_group[0] = dsum;
		}

	}
	//std::cout << " Hello world23! I am node " << node->rank << " out of " << node->size <<" and my hostname= "<< getenv("HOSTNAME")<< std::endl;

}

void MlOptimiserMpi::expectation()
{

#ifdef DEBUG
	std::cerr << "MlOptimiserMpi::expectation: Entering " << std::endl;
#endif

	MultidimArray<long int> first_last_nr_images(6);
	MultidimArray<RFLOAT> metadata;
	int first_slave = 1;
	// Use maximum of 100 particles for 3D and 10 particles for 2D estimations
	int n_trials_acc = (mymodel.ref_dim==3 && mymodel.data_dim != 3) ? 100 : 10;
	n_trials_acc = XMIPP_MIN(n_trials_acc, mydata.numberOfOriginalParticles());
	MPI_Status status;

	// Initialise some stuff
	// A. Update current size (may have been changed to ori_size in autoAdjustAngularSampling) and resolution pointers
	updateImageSizeAndResolutionPointers();

	// B. Set the PPref Fourier transforms, initialise wsum_model, etc.
	// The master only holds metadata, it does not set up the wsum_model (to save memory)
	if (!node->isMaster())
	{
		MlOptimiser::expectationSetup();

		// All slaves no longer need mydata.MD tables
		mydata.MDimg.clear();
		mydata.MDmic.clear();

		// Many small new's are not returned to the OS upon free-ing them. To force this, use the following call
		// from http://stackoverflow.com/questions/10943907/linux-allocator-does-not-release-small-chunks-of-memory
#if !defined(__APPLE__)
		malloc_trim(0);
#endif

	}

	// C. Calculate expected angular errors
	// Do not do this for maxCC
	// Only the first (reconstructing) slave (i.e. from half1) calculates expected angular errors
	////TMP if (!(iter==1 && do_firstiter_cc) && (mymodel.data_dim == 2) && !do_skip_align )
	if (!(iter==1 && do_firstiter_cc) && !(do_skip_align || do_skip_rotate) )
	{
		int my_nr_images, length_fn_ctf;
		if (node->isMaster())
		{
			// Master sends metadata (but not imagedata) for first 100 particles to first_slave (for calculateExpectedAngularErrors)
			MlOptimiser::getMetaAndImageDataSubset(0, n_trials_acc-1, false);
			my_nr_images = YSIZE(exp_metadata);
			node->relion_MPI_Send(&my_nr_images, 1, MPI_INT, first_slave, MPITAG_JOB_REQUEST, MPI_COMM_WORLD);
			node->relion_MPI_Send(MULTIDIM_ARRAY(exp_metadata), MULTIDIM_SIZE(exp_metadata), MY_MPI_DOUBLE, first_slave, MPITAG_METADATA, MPI_COMM_WORLD);
			// Also send exp_fn_ctfs if necessary
			length_fn_ctf = exp_fn_ctf.length() + 1; // +1 to include \0 at the end of the string
			node->relion_MPI_Send(&length_fn_ctf, 1, MPI_INT, first_slave, MPITAG_JOB_REQUEST, MPI_COMM_WORLD);
			if (length_fn_ctf > 1)
				node->relion_MPI_Send((void*)exp_fn_ctf.c_str(), length_fn_ctf, MPI_CHAR, first_slave, MPITAG_METADATA, MPI_COMM_WORLD);
		}
		else if (node->rank == first_slave)
		{
			// Slave has to receive all metadata from the master!
			node->relion_MPI_Recv(&my_nr_images, 1, MPI_INT, 0, MPITAG_JOB_REQUEST, MPI_COMM_WORLD, status);
			exp_metadata.resize(my_nr_images, METADATA_LINE_LENGTH_BEFORE_BODIES + (mymodel.nr_bodies) * METADATA_NR_BODY_PARAMS);
			node->relion_MPI_Recv(MULTIDIM_ARRAY(exp_metadata), MULTIDIM_SIZE(exp_metadata), MY_MPI_DOUBLE, 0, MPITAG_METADATA, MPI_COMM_WORLD, status);
			node->relion_MPI_Recv(&length_fn_ctf, 1, MPI_INT, 0, MPITAG_JOB_REQUEST, MPI_COMM_WORLD, status);
			if (length_fn_ctf > 1)
			{
				char* rec_buf2;
				rec_buf2 = (char *) malloc(length_fn_ctf);
				node->relion_MPI_Recv(rec_buf2, length_fn_ctf, MPI_CHAR, 0, MPITAG_METADATA, MPI_COMM_WORLD, status);
				exp_fn_ctf = rec_buf2;
				free(rec_buf2);
			}
			calculateExpectedAngularErrors(0, n_trials_acc-1);
		}

		// The reconstructing slave Bcast acc_rottilt, acc_psi, acc_trans to all other nodes!
		node->relion_MPI_Bcast(&acc_rot, 1, MY_MPI_DOUBLE, first_slave, MPI_COMM_WORLD);
		node->relion_MPI_Bcast(&acc_trans, 1, MY_MPI_DOUBLE, first_slave, MPI_COMM_WORLD);
	}

	// D. Update the angular sampling (all nodes except master)
	if (!node->isMaster() && do_auto_refine && iter > 1 )
	{
		updateAngularSampling(node->rank == 1);
	}
	node->relion_MPI_Bcast(&has_fine_enough_angular_sampling, 1, MPI_INT, first_slave, MPI_COMM_WORLD);

	// E. All nodes, except the master, check memory and precalculate AB-matrices for on-the-fly shifts
	if (!node->isMaster())
	{
		// Check whether everything fits into memory
		MlOptimiser::expectationSetupCheckMemory(node->rank == first_slave);

		// F. Precalculate AB-matrices for on-the-fly shifts
		if (do_shifts_onthefly)
			precalculateABMatrices();
	}
	// Slave 1 sends has_converged to everyone else (in particular the master needs it!)
	node->relion_MPI_Bcast(&has_converged, 1, MPI_INT, first_slave, MPI_COMM_WORLD);
	node->relion_MPI_Bcast(&do_join_random_halves, 1, MPI_INT, first_slave, MPI_COMM_WORLD);


	// Wait until expected angular errors have been calculated
	MPI_Barrier(MPI_COMM_WORLD);
	sleep(1);

	// Now perform real expectation step in parallel, use an on-demand master-slave system
#define JOB_FIRST (first_last_nr_images(0))
#define JOB_LAST  (first_last_nr_images(1))
#define JOB_NIMG  (first_last_nr_images(2))
#define JOB_LEN_FN_IMG  (first_last_nr_images(3))
#define JOB_LEN_FN_CTF  (first_last_nr_images(4))
#define JOB_LEN_FN_RECIMG  (first_last_nr_images(5))
#define JOB_NPAR  (JOB_LAST - JOB_FIRST + 1)

    if (node->isMaster())
    {
        try
        {
    		if (verb > 0)
    		{
				std::cout << " Expectation iteration " << iter;
				if (!do_auto_refine)
					std::cout << " of " << nr_iter;
				std::cout << std::endl;
				init_progress_bar(mydata.numberOfOriginalParticles());
    		}
			// Master distributes all packages of SomeParticles
			int nr_slaves_done = 0;
			int random_subset = 0;
			long int nr_ori_particles_done = 0;
			long int prev_step_done = nr_ori_particles_done;
			long int progress_bar_step_size = ROUND(mydata.numberOfOriginalParticles() / 80);
			long int nr_ori_particles_done_subset1 = 0;
			long int nr_ori_particles_done_subset2 = 0;
			long int my_nr_ori_particles_done = 0;

			while (nr_slaves_done < node->size - 1)
			{
				// Receive a job request from a slave
				node->relion_MPI_Recv(MULTIDIM_ARRAY(first_last_nr_images), MULTIDIM_SIZE(first_last_nr_images), MPI_LONG, MPI_ANY_SOURCE, MPITAG_JOB_REQUEST, MPI_COMM_WORLD, status);
				// Which slave sent this request?
				int this_slave = status.MPI_SOURCE;

				//#define DEBUG_MPIEXP2
#ifdef DEBUG_MPIEXP2
				std::cerr << " MASTER RECEIVING from slave= " << this_slave<< " JOB_FIRST= " << JOB_FIRST << " JOB_LAST= " << JOB_LAST
						<< " JOB_NIMG= "<<JOB_NIMG<< " JOB_NPAR= "<<JOB_NPAR<< std::endl;
#endif
				// The first time a slave reports it only asks for input, but does not send output of a previous processing task. In that case JOB_NIMG==0
				// Otherwise, the master needs to receive and handle the updated metadata from the slaves
				if (JOB_NIMG > 0)
				{
					exp_metadata.resize(JOB_NIMG, METADATA_LINE_LENGTH_BEFORE_BODIES + (mymodel.nr_bodies) * METADATA_NR_BODY_PARAMS);
					node->relion_MPI_Recv(MULTIDIM_ARRAY(exp_metadata), MULTIDIM_SIZE(exp_metadata), MY_MPI_DOUBLE, this_slave, MPITAG_METADATA, MPI_COMM_WORLD, status);

					// The master monitors the changes in the optimal orientations and classes
					monitorHiddenVariableChanges(JOB_FIRST, JOB_LAST);

					// The master then updates the mydata.MDimg table
					MlOptimiser::setMetaDataSubset(JOB_FIRST, JOB_LAST);
					if (verb > 0 && nr_ori_particles_done - prev_step_done > progress_bar_step_size)
					{
						prev_step_done = nr_ori_particles_done;
						progress_bar(nr_ori_particles_done + JOB_NPAR);
					}
				}

				// See which random_subset this slave belongs to, and keep track of the number of ori_particles that have been processed already
				if (do_split_random_halves)
				{
					random_subset = (this_slave % 2 == 1) ? 1 : 2;
					if (random_subset == 1)
					{
						my_nr_ori_particles_done = nr_ori_particles_done_subset1;
						// random_subset1 is stored in second half of OriginalParticles
						JOB_FIRST = nr_ori_particles_done_subset1;
						JOB_LAST  = XMIPP_MIN(mydata.numberOfOriginalParticles(1) - 1, JOB_FIRST + nr_pool - 1);
					}
					else
					{
						my_nr_ori_particles_done = nr_ori_particles_done_subset2;
						// random_subset2 is stored in second half of OriginalParticles
						JOB_FIRST = mydata.numberOfOriginalParticles(1) + nr_ori_particles_done_subset2;
						JOB_LAST  = XMIPP_MIN(mydata.numberOfOriginalParticles() - 1, JOB_FIRST + nr_pool - 1);
					}
				}
				else
				{
					random_subset = 0;
					my_nr_ori_particles_done = nr_ori_particles_done;
					JOB_FIRST = nr_ori_particles_done;
					JOB_LAST  = XMIPP_MIN(mydata.numberOfOriginalParticles() - 1, JOB_FIRST + nr_pool - 1);
				}

				// Now send out a new job
				if (my_nr_ori_particles_done < mydata.numberOfOriginalParticles(random_subset))
				{

					MlOptimiser::getMetaAndImageDataSubset(JOB_FIRST, JOB_LAST, !do_parallel_disc_io);
					JOB_NIMG = YSIZE(exp_metadata);
					JOB_LEN_FN_IMG = exp_fn_img.length() + 1; // +1 to include \0 at the end of the string
					JOB_LEN_FN_CTF = exp_fn_ctf.length() + 1;
					JOB_LEN_FN_RECIMG = exp_fn_recimg.length() + 1;
				}
				else
				{
					// There are no more particles in the list
					JOB_FIRST = -1;
					JOB_LAST = -1;
					JOB_NIMG = 0;
					JOB_LEN_FN_IMG = 0;
					JOB_LEN_FN_CTF = 0;
					JOB_LEN_FN_RECIMG = 0;
					exp_metadata.clear();
					exp_imagedata.clear();

					// No more particles, this slave is done now
					nr_slaves_done++;
				}
#ifdef DEBUG_MPIEXP2
				std::cerr << " MASTER SENDING to slave= " << this_slave<< " JOB_FIRST= " << JOB_FIRST << " JOB_LAST= " << JOB_LAST
								<< " JOB_NIMG= "<<JOB_NIMG<< " JOB_NPAR= "<<JOB_NPAR<< std::endl;
#endif
				node->relion_MPI_Send(MULTIDIM_ARRAY(first_last_nr_images), MULTIDIM_SIZE(first_last_nr_images), MPI_LONG, this_slave, MPITAG_JOB_REPLY, MPI_COMM_WORLD);

				// Master also sends the required metadata and imagedata for this job
				if (JOB_NIMG > 0)
				{
					node->relion_MPI_Send(MULTIDIM_ARRAY(exp_metadata), MULTIDIM_SIZE(exp_metadata), MY_MPI_DOUBLE, this_slave, MPITAG_METADATA, MPI_COMM_WORLD);
					if (do_parallel_disc_io)
					{
						node->relion_MPI_Send((void*)exp_fn_img.c_str(), JOB_LEN_FN_IMG, MPI_CHAR, this_slave, MPITAG_METADATA, MPI_COMM_WORLD);
						// Send filenames of images to the slaves
						if (JOB_LEN_FN_CTF > 1)
							node->relion_MPI_Send((void*)exp_fn_ctf.c_str(), JOB_LEN_FN_CTF, MPI_CHAR, this_slave, MPITAG_METADATA, MPI_COMM_WORLD);
						if (JOB_LEN_FN_RECIMG > 1)
							node->relion_MPI_Send((void*)exp_fn_recimg.c_str(), JOB_LEN_FN_RECIMG, MPI_CHAR, this_slave, MPITAG_METADATA, MPI_COMM_WORLD);
					}
					else
					{
						// Send imagedata to the slaves
						node->relion_MPI_Send(MULTIDIM_ARRAY(exp_imagedata), MULTIDIM_SIZE(exp_imagedata), MY_MPI_DOUBLE, this_slave, MPITAG_IMAGE, MPI_COMM_WORLD);
					}
				}

				// Update the total number of particles that has been done already
				nr_ori_particles_done += JOB_NPAR;
				if (do_split_random_halves)
				{
					// Also update the number of particles that has been done for each subset
					if (random_subset == 1)
						nr_ori_particles_done_subset1 += JOB_NPAR;
					else
						nr_ori_particles_done_subset2 += JOB_NPAR;
				}
			}
        }
        catch (RelionError XE)
        {
            std::cerr << "master encountered error: " << XE;
            MlOptimiser::usage();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

    }
    else
    {

    	try
    	{
			// Slaves do the real work (The slave does not need to know to which random_subset he belongs)

    		if (do_gpu)
    			for (int i = 0; i < cudaMlOptimisers.size(); i ++)
    				((MlOptimiserCuda *) cudaMlOptimisers[i])->resetData();

    		// Start off with an empty job request
			JOB_FIRST = 0;
			JOB_LAST = -1; // So that initial nr_particles (=JOB_LAST-JOB_FIRST+1) is zero!
			JOB_NIMG = 0;
			JOB_LEN_FN_IMG = 0;
			JOB_LEN_FN_CTF = 0;
			JOB_LEN_FN_RECIMG = 0;
			node->relion_MPI_Send(MULTIDIM_ARRAY(first_last_nr_images), MULTIDIM_SIZE(first_last_nr_images), MPI_LONG, 0, MPITAG_JOB_REQUEST, MPI_COMM_WORLD);

			while (true)
			{
#ifdef TIMING
				timer.tic(TIMING_MPISLAVEWAIT1);
#endif

				//Receive a new bunch of particles
				node->relion_MPI_Recv(MULTIDIM_ARRAY(first_last_nr_images), MULTIDIM_SIZE(first_last_nr_images), MPI_LONG, 0, MPITAG_JOB_REPLY, MPI_COMM_WORLD, status);
#ifdef TIMING
				timer.toc(TIMING_MPISLAVEWAIT1);
#endif

				//Check whether I am done
				if (JOB_NIMG <= 0)
				{
#ifdef DEBUG
					std::cerr <<" slave "<< node->rank << " has finished expectation.."<<std::endl;
#endif
					exp_imagedata.clear();
					exp_metadata.clear();
					break;
				}
				else
				{
#ifdef TIMING
					timer.tic(TIMING_MPISLAVEWAIT2);
#endif
					// Also receive the imagedata and the metadata for these images from the master
					exp_metadata.resize(JOB_NIMG, METADATA_LINE_LENGTH_BEFORE_BODIES + (mymodel.nr_bodies) * METADATA_NR_BODY_PARAMS);
					node->relion_MPI_Recv(MULTIDIM_ARRAY(exp_metadata), MULTIDIM_SIZE(exp_metadata), MY_MPI_DOUBLE, 0, MPITAG_METADATA, MPI_COMM_WORLD, status);

					// Receive the image filenames or the exp_imagedata
					if (do_parallel_disc_io)
					{
						// Resize the exp_fn_img strings
				        char* rec_buf;
				        rec_buf = (char *) malloc(JOB_LEN_FN_IMG);
				        node->relion_MPI_Recv(rec_buf, JOB_LEN_FN_IMG, MPI_CHAR, 0, MPITAG_METADATA, MPI_COMM_WORLD, status);
				        exp_fn_img = rec_buf;
				        free(rec_buf);
						if (JOB_LEN_FN_CTF > 1)
						{
					        char* rec_buf2;
					        rec_buf2 = (char *) malloc(JOB_LEN_FN_CTF);
					        node->relion_MPI_Recv(rec_buf2, JOB_LEN_FN_CTF, MPI_CHAR, 0, MPITAG_METADATA, MPI_COMM_WORLD, status);
					        exp_fn_ctf = rec_buf2;
					        free(rec_buf2);

						}
						if (JOB_LEN_FN_RECIMG > 1)
						{
					        char* rec_buf3;
					        rec_buf3 = (char *) malloc(JOB_LEN_FN_RECIMG);
					        node->relion_MPI_Recv(rec_buf3, JOB_LEN_FN_RECIMG, MPI_CHAR, 0, MPITAG_METADATA, MPI_COMM_WORLD, status);
					        exp_fn_recimg = rec_buf3;
					        free(rec_buf3);
						}
					}
					else
					{
						// resize the exp_imagedata array
						if (mymodel.data_dim == 3)
						{

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
								exp_imagedata.resize(2*JOB_NIMG, mymodel.ori_size, mymodel.ori_size);
							else
								exp_imagedata.resize(JOB_NIMG, mymodel.ori_size, mymodel.ori_size);
						}
						node->relion_MPI_Recv(MULTIDIM_ARRAY(exp_imagedata), MULTIDIM_SIZE(exp_imagedata), MY_MPI_DOUBLE, 0, MPITAG_IMAGE, MPI_COMM_WORLD, status);
					}

					// Now process these images
#ifdef DEBUG_MPIEXP
					std::cerr << " SLAVE EXECUTING node->rank= " << node->rank << " JOB_FIRST= " << JOB_FIRST << " JOB_LAST= " << JOB_LAST << std::endl;
#endif
#ifdef TIMING
					timer.toc(TIMING_MPISLAVEWAIT2);
					timer.tic(TIMING_MPISLAVEWORK);
#endif
					expectationSomeParticles(JOB_FIRST, JOB_LAST);
#ifdef TIMING
					timer.toc(TIMING_MPISLAVEWORK);
					timer.tic(TIMING_MPISLAVEWAIT3);
#endif

					// Report to the master how many particles I have processed
					node->relion_MPI_Send(MULTIDIM_ARRAY(first_last_nr_images), MULTIDIM_SIZE(first_last_nr_images), MPI_LONG, 0, MPITAG_JOB_REQUEST, MPI_COMM_WORLD);
					// Also send the metadata belonging to those
					node->relion_MPI_Send(MULTIDIM_ARRAY(exp_metadata), MULTIDIM_SIZE(exp_metadata), MY_MPI_DOUBLE, 0, MPITAG_METADATA, MPI_COMM_WORLD);

#ifdef TIMING
					timer.toc(TIMING_MPISLAVEWAIT3);
#endif

				}

			}

			if (do_gpu)
			{
				for (int i = 0; i < cudaMlOptimisers.size(); i ++)
				{
					for (int iclass = 0; iclass < wsum_model.nr_classes; iclass++)
					{
						unsigned long s = wsum_model.BPref[iclass].data.nzyxdim;
						XFLOAT *reals = new XFLOAT[s];
						XFLOAT *imags = new XFLOAT[s];
						XFLOAT *weights = new XFLOAT[s];

						( (MlOptimiserCuda*) cudaMlOptimisers[i])->cudaBackprojectors[iclass].getMdlData(reals, imags, weights);

						for (unsigned long n = 0; n < s; n++)
						{
							wsum_model.BPref[iclass].data.data[n].real += (RFLOAT) reals[n];
							wsum_model.BPref[iclass].data.data[n].imag += (RFLOAT) imags[n];
							wsum_model.BPref[iclass].weight.data[n] += (RFLOAT) weights[n];
						}

						delete [] reals;
						delete [] imags;
						delete [] weights;
					}
				}
			}


    	}
        catch (RelionError XE)
        {
            std::cerr << "slave "<< node->rank << " encountered error: " << XE;
            MlOptimiser::usage();
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }

    }

    // Just make sure the temporary arrays are empty...
	exp_imagedata.clear();
	exp_metadata.clear();

	if (verb > 0)
		progress_bar(mydata.numberOfOriginalParticles());

#ifdef TIMING
    // Measure how long I have to wait for the rest
    timer.tic(TIMING_MPIWAIT);
    node->barrierWait();
    timer.toc(TIMING_MPIWAIT);
#endif

	// Wait until expected angular errors have been calculated
	MPI_Barrier(MPI_COMM_WORLD);

	// All slaves reset the size of their projector to zero tosave memory
	if (!node->isMaster())
	{
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
			mymodel.PPref[iclass].initialiseData(0);
	}


#ifdef DEBUG
	std::cerr << "MlOptimiserMpi::expectation: done" << std::endl;
#endif

}


void MlOptimiserMpi::combineAllWeightedSumsViaFile()
{

#ifdef TIMING
    timer.tic(TIMING_MPICOMBINEDISC);
#endif
	MultidimArray<RFLOAT> Mpack;
	FileName fn_pack;

	int nr_subsets = (do_split_random_halves) ? 2 : 1;

	// Only need to combine if there are more than one slaves per subset!
	if ((node->size - 1)/nr_subsets > 1)
	{
		// A. First all slaves pack up their wsum_model (this is done simultaneously)
		if (!node->isMaster())
		{
			wsum_model.pack(Mpack); // use negative piece and nr_pieces to only make a single Mpack, i.e. do not split into multiple pieces
		}

		// B. All slaves write their Mpack to disc. Do this SEQUENTIALLY to prevent heavy load on disc I/O
		for (int this_slave = 1; this_slave < node->size; this_slave++ )
		{
			if (this_slave == node->rank)
			{
				fn_pack.compose(fn_out+"_rank", node->rank, "tmp");
				if (fn_scratch != "")
					fn_pack = fn_scratch + "/" + fn_pack;

				Mpack.writeBinary(fn_pack);
				//std::cerr << "Rank "<< node->rank <<" has written: "<<fn_pack << " sum= "<<Mpack.sum()<< std::endl;
			}
			if (!do_parallel_disc_io)
				MPI_Barrier(MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// C. First slave of each subset reads all other slaves' Mpack; sum; and write sum to disc
		// Again, do this SEQUENTIALLY to prevent heavy load on disc I/O
		for (int first_slave = 1; first_slave <= nr_subsets; first_slave++)
		{
			if (node->rank == first_slave)
			{
				for (int other_slave = first_slave + nr_subsets; other_slave < node->size; other_slave+= nr_subsets )
				{
					fn_pack.compose(fn_out+"_rank", other_slave, "tmp");
					if (fn_scratch != "")
						fn_pack = fn_scratch + "/" + fn_pack;

					Mpack.readBinaryAndSum(fn_pack);
					//std::cerr << "Slave "<<node->rank<<" has read "<<fn_pack<< " sum= "<<Mpack.sum() << std::endl;
				}
				// After adding all Mpacks together: write the sum to disc
				fn_pack.compose(fn_out+"_rank", node->rank, "tmp");
				if (fn_scratch != "")
					fn_pack = fn_scratch + "/" + fn_pack;

				Mpack.writeBinary(fn_pack);
				//std::cerr << "Slave "<<node->rank<<" is writing total SUM in "<<fn_pack << std::endl;
			}
			if (!do_parallel_disc_io)
				MPI_Barrier(MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);


		// D. All other slaves read the summed Mpack from their first_slave
		// Again, do this SEQUENTIALLY to prevent heavy load on disc I/O
		for (int this_slave = nr_subsets + 1; this_slave < node->size; this_slave++ )
		{
			if (this_slave == node->rank)
			{
				// Find out who is the first slave in this subset
				int first_slave;
				if (!do_split_random_halves)
					first_slave = 1;
				else
					first_slave = (this_slave % 2 == 1) ? 1 : 2;

				// Read the corresponding Mpack (which now contains the sum of all Mpacks)
				fn_pack.compose(fn_out+"_rank", first_slave, "tmp");
				if (fn_scratch != "")
					fn_pack = fn_scratch + "/" + fn_pack;

				Mpack.readBinary(fn_pack);
				//std::cerr << "Rank "<< node->rank <<" has read: "<<fn_pack << " sum= "<<Mpack.sum()<< std::endl;
			}
			if (!do_parallel_disc_io)
				MPI_Barrier(MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// E. All slaves delete their own temporary file
		// Again, do this SEQUENTIALLY to prevent heavy load on disc I/O
		for (int this_slave = 1; this_slave < node->size; this_slave++ )
		{
			if (this_slave == node->rank)
			{
				fn_pack.compose(fn_out+"_rank", node->rank, "tmp");
				if (fn_scratch != "")
					fn_pack = fn_scratch + "/" + fn_pack;

				remove((fn_pack).c_str());
				//std::cerr << "Rank "<< node->rank <<" has deleted: "<<fn_pack << std::endl;
			}
			if (!do_parallel_disc_io)
				MPI_Barrier(MPI_COMM_WORLD);
		}

		// F. Finally all slaves unpack Msum into their wsum_model (do this simultaneously)
		if (!node->isMaster())
			wsum_model.unpack(Mpack);

	} // end if ((node->size - 1)/nr_subsets > 1)
#ifdef TIMING
    timer.toc(TIMING_MPICOMBINEDISC);
#endif

}

void MlOptimiserMpi::combineAllWeightedSums()
{

#ifdef TIMING
    timer.tic(TIMING_MPICOMBINENETW);
#endif

    // Pack all weighted sums in Mpack
	MultidimArray<RFLOAT> Mpack, Msum;
	MPI_Status status;

	// First slave manually sums over all other slaves of it's subset
	// And then sends results back to all those slaves
	// When splitting the data into two random halves, perform two passes: one for each subset
	int nr_subsets = (do_split_random_halves) ? 2 : 1;

	// Only combine weighted sums if there are more than one slaves per subset!
	if ((node->size - 1)/nr_subsets > 1)
	{
		// Loop over possibly multiple instances of Mpack of maximum size
		int piece = 0;
		int nr_pieces = 1;
		long int pack_size;
		while (piece < nr_pieces)
		{
			// All nodes except those who will reset nr_pieces piece will pass while loop in next pass
			nr_pieces = 0;

			// First all slaves pack up their wsum_model
			if (!node->isMaster())
			{
				wsum_model.pack(Mpack, piece, nr_pieces);
				// The first slave(s) set Msum equal to Mpack, the others initialise to zero
				if (node->rank <= nr_subsets)
					Msum = Mpack;
				else
					Msum.initZeros(Mpack);
			}


			// Loop through all slaves: each slave sends its Msum to the next slave for its subset.
			// Each next slave sums its own Mpack to the received Msum and sends it on to the next slave

			for (int this_slave = 1; this_slave < node->size; this_slave++ )
			{
				// Find out who is the first slave in this subset
				int first_slave;
				if (!do_split_random_halves)
					first_slave = 1;
				else
					first_slave = (this_slave % 2 == 1) ? 1 : 2;

				// Find out who is the next slave in this subset
				int other_slave = this_slave + nr_subsets;

				if (other_slave < node->size)
				{
					if (node->rank == this_slave)
					{
#ifdef DEBUG
						std::cerr << " AA SEND node->rank= " << node->rank << " MULTIDIM_SIZE(Msum)= "<< MULTIDIM_SIZE(Msum)
								<< " this_slave= " << this_slave << " other_slave= "<<other_slave << std::endl;
#endif
						node->relion_MPI_Send(MULTIDIM_ARRAY(Msum), MULTIDIM_SIZE(Msum), MY_MPI_DOUBLE, other_slave, MPITAG_PACK, MPI_COMM_WORLD);
					}
					else if (node->rank == other_slave)
					{
						MPI_Recv(MULTIDIM_ARRAY(Msum), MULTIDIM_SIZE(Msum), MY_MPI_DOUBLE, this_slave, MPITAG_PACK, MPI_COMM_WORLD, &status);
#ifdef DEBUG
						std::cerr << " AA RECV node->rank= " << node->rank  << " MULTIDIM_SIZE(Msum)= "<< MULTIDIM_SIZE(Msum)
								<< " this_slave= " << this_slave << " other_slave= "<<other_slave << std::endl;
#endif
						// Add my own Mpack to send onwards in the next step
						Msum += Mpack;
					}
				}
				else
				{
					// Now this_slave has reached the last slave, which passes the final Msum to the first one (i.e. first_slave)
					if (node->rank == this_slave)
					{
#ifdef DEBUG
						std::cerr << " BB SEND node->rank= " << node->rank  << " MULTIDIM_SIZE(Msum)= "<< MULTIDIM_SIZE(Msum)
								<< " this_slave= " << this_slave << " first_slave= "<<first_slave << std::endl;
#endif
						node->relion_MPI_Send(MULTIDIM_ARRAY(Msum), MULTIDIM_SIZE(Msum), MY_MPI_DOUBLE, first_slave, MPITAG_PACK, MPI_COMM_WORLD);
					}
					else if (node->rank == first_slave)
					{
						node->relion_MPI_Recv(MULTIDIM_ARRAY(Msum), MULTIDIM_SIZE(Msum), MY_MPI_DOUBLE, this_slave, MPITAG_PACK, MPI_COMM_WORLD, status);
#ifdef DEBUG
						std::cerr << " BB RECV node->rank= " << node->rank  << " MULTIDIM_SIZE(Msum)= "<< MULTIDIM_SIZE(Msum)
								<< " this_slave= " << this_slave << " first_slave= "<<first_slave << std::endl;
#endif
					}
				}
			} // end for this_slave

			// Now loop through all slaves again to pass around the Msum
			for (int this_slave = 1; this_slave < node->size; this_slave++ )
			{
				// Find out who is the next slave in this subset
				int other_slave = this_slave + nr_subsets;

				// Do not send to the last slave, because it already had its Msum from the cycle above, therefore subtract nr_subsets from node->size
				if (other_slave < node->size - nr_subsets)
				{
					if (node->rank == this_slave)
					{
#ifdef DEBUG
						std::cerr << " CC SEND node->rank= " << node->rank << " MULTIDIM_SIZE(Msum)= "<< MULTIDIM_SIZE(Msum)
								<< " this_slave= " << this_slave << " other_slave= "<<other_slave << std::endl;
#endif
						node->relion_MPI_Send(MULTIDIM_ARRAY(Msum), MULTIDIM_SIZE(Msum), MY_MPI_DOUBLE, other_slave, MPITAG_PACK, MPI_COMM_WORLD);
					}
					else if (node->rank == other_slave)
					{
						node->relion_MPI_Recv(MULTIDIM_ARRAY(Msum), MULTIDIM_SIZE(Msum), MY_MPI_DOUBLE, this_slave, MPITAG_PACK, MPI_COMM_WORLD, status);
#ifdef DEBUG
						std::cerr << " CC RECV node->rank= " << node->rank << " MULTIDIM_SIZE(Msum)= "<< MULTIDIM_SIZE(Msum)
								<< " this_slave= " << this_slave << " other_slave= "<<other_slave << std::endl;
#endif
					}
				}
			} // end for this_slave


			// Finally all slaves unpack Msum into their wsum_model
			if (!node->isMaster())
			{
				// Subtract 1 from piece because it was incremented already...
				wsum_model.unpack(Msum, piece - 1);
			}


		} // end for piece

		MPI_Barrier(MPI_COMM_WORLD);
	}

#ifdef TIMING
    timer.toc(TIMING_MPICOMBINENETW);
#endif

}

void MlOptimiserMpi::combineWeightedSumsTwoRandomHalvesViaFile()
{
	// Just sum the weighted halves from slave 1 and slave 2 and Bcast to everyone else
	if (!do_split_random_halves)
		REPORT_ERROR("MlOptimiserMpi::combineWeightedSumsTwoRandomHalvesViaFile BUG: you cannot combineWeightedSumsTwoRandomHalves if you have not split random halves");

	MultidimArray<RFLOAT> Mpack;
	FileName fn_pack = fn_out + ".tmp";
        if (fn_scratch != "")
            fn_pack = fn_scratch + "/" + fn_pack;

	// Everyone packs up his wsum_model (simultaneously)
	// The slaves from 3 and onwards also need this in order to have the correct Mpack size to be able to read in the summed Mpack
	if (!node->isMaster())
		wsum_model.pack(Mpack);

	// Rank 2 writes it Mpack to file
	if (node->rank == 2)
	{
		Mpack.writeBinary(fn_pack);
		//std::cerr << "Rank "<< node->rank <<" has written: "<<fn_pack << " sum= "<<Mpack.sum()<< std::endl;
	}

	// Wait until rank2 is ready
	MPI_Barrier(MPI_COMM_WORLD);

	// Now rank1 reads the file of rank2 and adds it to its own Mpack and (over)write the sum to disc
	if (node->rank == 1)
	{
		Mpack.readBinaryAndSum(fn_pack);
		Mpack.writeBinary(fn_pack);
	}

	// Now all slaves read the sum and unpack
	// Do this sequentially in order not to have very heavy disc I/O
	for (int this_slave = 2; this_slave < node->size; this_slave++)
	{
		if (node->rank == this_slave)
			Mpack.readBinary(fn_pack);
		if (!do_parallel_disc_io)
			MPI_Barrier(MPI_COMM_WORLD);
	}

	// in case we're doing do_parallel_disc_io
	MPI_Barrier(MPI_COMM_WORLD);

	// Remove temporary file
	if (node->rank == 1)
		remove(fn_pack.c_str());

	// Then everyone except the master unpacks
	if (!node->isMaster())
		wsum_model.unpack(Mpack);

}

void MlOptimiserMpi::combineWeightedSumsTwoRandomHalves()
{
	// Just sum the weighted halves from slave 1 and slave 2 and Bcast to everyone else
	if (!do_split_random_halves)
		REPORT_ERROR("MlOptimiserMpi::combineWeightedSumsTwoRandomHalves BUG: you cannot combineWeightedSumsTwoRandomHalves if you have not split random halves");

	MultidimArray<RFLOAT> Mpack, Msum;
	MPI_Status status;

	int piece = 0;
	int nr_pieces = 1;
	long int pack_size;
	while (piece < nr_pieces)
	{
		// All nodes except those who will reset nr_pieces will pass while next time
		nr_pieces = 0;

		if (node->rank == 2)
		{
			wsum_model.pack(Mpack, piece, nr_pieces, false); // do not clear the model!
			node->relion_MPI_Send(MULTIDIM_ARRAY(Mpack), MULTIDIM_SIZE(Mpack), MY_MPI_DOUBLE, 1, MPITAG_PACK, MPI_COMM_WORLD);
			Mpack.clear();
		}
		else if (node->rank == 1)
		{
			if (verb > 0) std::cout << " Combining two random halves ..."<< std::endl;
			wsum_model.pack(Mpack, piece, nr_pieces);
			Msum.initZeros(Mpack);
			node->relion_MPI_Recv(MULTIDIM_ARRAY(Msum), MULTIDIM_SIZE(Msum), MY_MPI_DOUBLE, 2, MPITAG_PACK, MPI_COMM_WORLD, status);
			Msum += Mpack;
			// Unpack the sum (subtract 1 from piece because it was incremented already...)
			wsum_model.unpack(Msum, piece - 1);
			Msum.clear();
			Mpack.clear();
		}
	}

	// Now slave 1 has the sum of the two halves
	MPI_Barrier(MPI_COMM_WORLD);

	// Bcast to everyone
	piece = 0;
	nr_pieces = 1;
	while (piece < nr_pieces)
	{

		// All nodes except those who will reset nr_pieces will pass while next time
		nr_pieces = 0;

		// The master does not have a wsum_model!
		if (!node->isMaster())
		{
			// Let's have everyone repack their Mpack, so we know the size etc
			wsum_model.pack(Mpack, piece, nr_pieces);

			// rank one sends Mpack to everyone else
			if (node->rank == 1)
			{
				for (int other_slave = 2; other_slave < node->size; other_slave++)
					node->relion_MPI_Send(MULTIDIM_ARRAY(Mpack), MULTIDIM_SIZE(Mpack), MY_MPI_DOUBLE, other_slave, MPITAG_PACK, MPI_COMM_WORLD);
			}
			else
			{
				node->relion_MPI_Recv(MULTIDIM_ARRAY(Mpack), MULTIDIM_SIZE(Mpack), MY_MPI_DOUBLE, 1, MPITAG_PACK, MPI_COMM_WORLD, status);
			}

			// Everyone unpacks the new Mpack
			wsum_model.unpack(Mpack, piece - 1);
			Mpack.clear();

		}

	}

}

void MlOptimiserMpi::maximization()
{
#ifdef DEBUG
	std::cerr << "MlOptimiserMpi::maximization: Entering " << std::endl;
#endif

#ifdef TIMING
		timer.tic(TIMING_RECONS);
#endif

	if (verb > 0)
	{
		std::cout << " Maximization ..."<< std::endl;
		init_progress_bar(mymodel.nr_classes);
	}

	// Aug05,2015 - Shaoda, helical refinement
	helical_twist_half1 = helical_twist_half2 = helical_twist;
	helical_rise_half1 = helical_rise_half2 = helical_rise;

	// First reconstruct all classes in parallel
	for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++)
	{
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		{
			// either ibody or iclass can be larger than 0, never 2 at the same time!
			int ith_recons = (mymodel.nr_bodies > 1) ? ibody : iclass;

			if (mymodel.pdf_class[iclass] > 0.)
			{
				// Parallelise: each MPI-node has a different reference
				int reconstruct_rank1;
				if (do_split_random_halves)
					reconstruct_rank1 = 2 * (ith_recons % ( (node->size - 1)/2 ) ) + 1;
				else
					reconstruct_rank1 = ith_recons % (node->size - 1) + 1;

				if (node->rank == reconstruct_rank1)
				{
					wsum_model.BPref[ith_recons].reconstruct(mymodel.Iref[ith_recons], gridding_nr_iter, do_map,
							mymodel.tau2_fudge_factor, mymodel.tau2_class[ith_recons], mymodel.sigma2_class[ith_recons],
							mymodel.data_vs_prior_class[ith_recons], mymodel.fsc_halves_class, wsum_model.pdf_class[iclass],
							do_split_random_halves, (do_join_random_halves || do_always_join_random_halves), nr_threads, minres_map);

					// Also perform the unregularized reconstruction
					if (do_auto_refine && has_converged)
						readTemporaryDataAndWeightArraysAndReconstruct(ith_recons, 1);

					// Apply the body mask
					if (mymodel.nr_bodies > 1)
					{
						// 19may2015 translate the reconstruction back to its C.O.M.
						selfTranslate(mymodel.Iref[ith_recons], mymodel.com_bodies[ibody], DONT_WRAP);

						// Also write out unmasked body recontruction
						FileName fn_tmp;
						fn_tmp.compose(fn_out + "_unmasked_half1_body", ibody+1,"mrc");
						Image<RFLOAT> Itmp;
						Itmp()=mymodel.Iref[ith_recons];
						Itmp.write(fn_tmp);
						mymodel.Iref[ith_recons].setXmippOrigin();
						mymodel.Iref[ith_recons] *= mymodel.masks_bodies[ith_recons];
					}

					// Shaoda Jul26,2015 - Helical symmetry local refinement
					if ( (iter > 1) && (do_auto_refine) && (do_helical_refine) && (do_helical_symmetry_local_refinement) )
					{
						localSearchHelicalSymmetry(
								mymodel.Iref[ith_recons],
								mymodel.pixel_size,
								(particle_diameter / 2.),
								(helical_mask_tube_inner_diameter / 2.),
								(helical_mask_tube_outer_diameter / 2.),
								helical_central_proportion,
								helical_rise,
								helical_twist,
								helical_symmetry_local_refinement_max_dev,
								helical_symmetry_local_refinement_max_dev,
								helical_rise_half1,
								helical_twist_half1);
					}
					// Sjors & Shaoda Apr 2015 - Apply real space helical symmetry and real space Z axis expansion.
					if( (do_helical_refine) && (!has_converged) )
					{
						makeHelicalReferenceInRealSpace(
								mymodel.Iref[ith_recons],
								mymodel.pixel_size,
								helical_twist_half1,
								helical_rise_half1,
								helical_central_proportion,
								(particle_diameter / 2.),
								(helical_mask_tube_inner_diameter / 2.),
								(helical_mask_tube_outer_diameter / 2.),
								width_mask_edge);
					}
				}

				// In some cases there is not enough memory to reconstruct two random halves in parallel
				// Therefore the following option exists to perform them sequentially
				if (do_sequential_halves_recons)
					MPI_Barrier(MPI_COMM_WORLD);

				// When splitting the data into two random halves, perform two reconstructions in parallel: one for each subset
				if (do_split_random_halves)
				{
					int reconstruct_rank2 = 2 * (ith_recons % ( (node->size - 1)/2 ) ) + 2;

					if (node->rank == reconstruct_rank2)
					{
						// Rank 2 does not need to do the joined reconstruction
						if (!do_join_random_halves)
							wsum_model.BPref[ith_recons].reconstruct(mymodel.Iref[ith_recons], gridding_nr_iter, do_map,
									mymodel.tau2_fudge_factor, mymodel.tau2_class[ith_recons], mymodel.sigma2_class[ith_recons],
									mymodel.data_vs_prior_class[ith_recons], mymodel.fsc_halves_class, wsum_model.pdf_class[iclass],
									do_split_random_halves, do_join_random_halves, nr_threads, minres_map);

						// But rank 2 always does the unfiltered reconstruction
						if (do_auto_refine && has_converged)
							readTemporaryDataAndWeightArraysAndReconstruct(ith_recons, 2);

						// Apply the body mask
						if (mymodel.nr_bodies > 1)
						{
							// 19may2015 translate the reconstruction back to its C.O.M.
							selfTranslate(mymodel.Iref[ith_recons], mymodel.com_bodies[ibody], DONT_WRAP);

							FileName fn_tmp;
							fn_tmp.compose(fn_out + "_unmasked_half2_body", ibody+1,"mrc");
							Image<RFLOAT> Itmp;
							Itmp()=mymodel.Iref[ith_recons];
							Itmp.write(fn_tmp);
							mymodel.Iref[ith_recons].setXmippOrigin();
							mymodel.Iref[ith_recons] *= mymodel.masks_bodies[ith_recons];
						}

						// Shaoda Jul26,2015 - Helical symmetry local refinement
						if ( (iter > 1) && (do_auto_refine) && (do_helical_refine) && (do_helical_symmetry_local_refinement) )
						{
							localSearchHelicalSymmetry(
									mymodel.Iref[ith_recons],
									mymodel.pixel_size,
									(particle_diameter / 2.),
									(helical_mask_tube_inner_diameter / 2.),
									(helical_mask_tube_outer_diameter / 2.),
									helical_central_proportion,
									helical_rise,
									helical_twist,
									helical_symmetry_local_refinement_max_dev,
									helical_symmetry_local_refinement_max_dev,
									helical_rise_half2,
									helical_twist_half2);
						}
						// Sjors & Shaoda Apr 2015 - Apply real space helical symmetry and real space Z axis expansion.
						if( (do_helical_refine) && (!has_converged) )
						{
							makeHelicalReferenceInRealSpace(
									mymodel.Iref[ith_recons],
									mymodel.pixel_size,
									helical_twist_half2,
									helical_rise_half2,
									helical_central_proportion,
									(particle_diameter / 2.),
									(helical_mask_tube_inner_diameter / 2.),
									(helical_mask_tube_outer_diameter / 2.),
									width_mask_edge);
						}

					}
				}

			} // endif pdf_class[iclass] > 0.
			else
			{
				mymodel.Iref[ith_recons].initZeros();
			}

//#define DEBUG_RECONSTRUCT
#ifdef DEBUG_RECONSTRUCT
			MPI_Barrier( MPI_COMM_WORLD);
#endif
		} // end for iclass
	} // end for ibody

#ifdef DEBUG
	std::cerr << "rank= "<<node->rank<<" has reached barrier of reconstruction" << std::endl;
#endif
	MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
	std::cerr << "All classes have been reconstructed" << std::endl;
#endif

	// Once reconstructed, broadcast new models to all other nodes
	// This cannot be done in the reconstruction loop itself because then it will be executed sequentially
	for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++)
	{
		for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
		{
			// either ibody or iclass can be larger than 0, never 2 at the same time!
			int ith_recons = (mymodel.nr_bodies > 1) ? ibody : iclass;

			if (do_split_random_halves && !do_join_random_halves)
			{
				MPI_Status status;
				// Make sure I am sending from the rank where the reconstruction was done (see above) to all other slaves of this subset
				// Loop twice through this, as each class was reconstructed by two different slaves!!
				int nr_subsets = 2;
				for (int isubset = 1; isubset <= nr_subsets; isubset++)
				{
					if (node->myRandomSubset() == isubset)
					{
						int reconstruct_rank = 2 * (ith_recons % ( (node->size - 1)/2 ) ) + isubset; // first pass subset1, second pass subset2
						int my_first_recv = node->myRandomSubset();

						for (int recv_node = my_first_recv; recv_node < node->size; recv_node += nr_subsets)
						{
							if (node->rank == reconstruct_rank && recv_node != node->rank)
							{
#ifdef DEBUG
								std::cerr << "isubset= "<<isubset<<" Sending iclass="<<iclass<<" from node "<<reconstruct_rank<<" to node "<<recv_node << std::endl;
#endif
								node->relion_MPI_Send(MULTIDIM_ARRAY(mymodel.Iref[ith_recons]), MULTIDIM_SIZE(mymodel.Iref[ith_recons]), MY_MPI_DOUBLE, recv_node, MPITAG_IMAGE, MPI_COMM_WORLD);
								node->relion_MPI_Send(MULTIDIM_ARRAY(mymodel.data_vs_prior_class[iclass]), MULTIDIM_SIZE(mymodel.data_vs_prior_class[iclass]), MY_MPI_DOUBLE, recv_node, MPITAG_METADATA, MPI_COMM_WORLD);
								node->relion_MPI_Send(MULTIDIM_ARRAY(mymodel.sigma2_class[iclass]), MULTIDIM_SIZE(mymodel.sigma2_class[iclass]), MY_MPI_DOUBLE, recv_node, MPITAG_RFLOAT, MPI_COMM_WORLD);
								node->relion_MPI_Send(MULTIDIM_ARRAY(mymodel.fsc_halves_class), MULTIDIM_SIZE(mymodel.fsc_halves_class), MY_MPI_DOUBLE, recv_node, MPITAG_RANDOMSEED, MPI_COMM_WORLD);
							}
							else if (node->rank != reconstruct_rank && node->rank == recv_node)
							{
								//std::cerr << "isubset= "<<isubset<< " Receiving iclass="<<iclass<<" from node "<<reconstruct_rank<<" at node "<<node->rank<< std::endl;
								node->relion_MPI_Recv(MULTIDIM_ARRAY(mymodel.Iref[ith_recons]), MULTIDIM_SIZE(mymodel.Iref[ith_recons]), MY_MPI_DOUBLE, reconstruct_rank, MPITAG_IMAGE, MPI_COMM_WORLD, status);
								node->relion_MPI_Recv(MULTIDIM_ARRAY(mymodel.data_vs_prior_class[iclass]), MULTIDIM_SIZE(mymodel.data_vs_prior_class[iclass]), MY_MPI_DOUBLE, reconstruct_rank, MPITAG_METADATA, MPI_COMM_WORLD, status);
								node->relion_MPI_Recv(MULTIDIM_ARRAY(mymodel.sigma2_class[iclass]), MULTIDIM_SIZE(mymodel.sigma2_class[iclass]), MY_MPI_DOUBLE, reconstruct_rank, MPITAG_RFLOAT, MPI_COMM_WORLD, status);
								node->relion_MPI_Recv(MULTIDIM_ARRAY(mymodel.fsc_halves_class), MULTIDIM_SIZE(mymodel.fsc_halves_class), MY_MPI_DOUBLE, reconstruct_rank, MPITAG_RANDOMSEED, MPI_COMM_WORLD, status);
#ifdef DEBUG
								std::cerr << "isubset= "<<isubset<< " Received!!!="<<iclass<<" from node "<<reconstruct_rank<<" at node "<<node->rank<< std::endl;
#endif
							}
						}

					}
				}
				// No one should continue until we're all here
				MPI_Barrier(MPI_COMM_WORLD);

				// Now all slaves have all relevant reconstructions
				// TODO: someone should also send reconstructions to the master (for comparison with other subset?)
			}
			else
			{
				int reconstruct_rank = (ith_recons % (node->size - 1) ) + 1;
				// Broadcast the reconstructed references to all other MPI nodes
				node->relion_MPI_Bcast(MULTIDIM_ARRAY(mymodel.Iref[ith_recons]),
						MULTIDIM_SIZE(mymodel.Iref[ith_recons]), MY_MPI_DOUBLE, reconstruct_rank, MPI_COMM_WORLD);
				// Broadcast the data_vs_prior spectra to all other MPI nodes
				node->relion_MPI_Bcast(MULTIDIM_ARRAY(mymodel.data_vs_prior_class[iclass]),
						MULTIDIM_SIZE(mymodel.data_vs_prior_class[iclass]), MY_MPI_DOUBLE, reconstruct_rank, MPI_COMM_WORLD);
				// Broadcast the sigma2_class spectra to all other MPI nodes
				node->relion_MPI_Bcast(MULTIDIM_ARRAY(mymodel.sigma2_class[iclass]),
						MULTIDIM_SIZE(mymodel.sigma2_class[iclass]), MY_MPI_DOUBLE, reconstruct_rank, MPI_COMM_WORLD);

			}

			// Re-set the origin (this may be lost in some cases??)
			mymodel.Iref[ith_recons].setXmippOrigin();

			// Aug05,2015 - Shaoda, helical symmetry refinement, broadcast refined helical parameters
			if ( (mymodel.pdf_class[iclass] > 0.)
					&& (iter > 1) && (do_auto_refine) && (do_helical_refine) && (do_helical_symmetry_local_refinement) )
			{
				int reconstruct_rank1;
				if (do_split_random_halves)
					reconstruct_rank1 = 2 * (ith_recons % ( (node->size - 1)/2 ) ) + 1;
				else
					reconstruct_rank1 = ith_recons % (node->size - 1) + 1;
				node->relion_MPI_Bcast(&helical_twist_half1, 1, MY_MPI_DOUBLE, reconstruct_rank1, MPI_COMM_WORLD);
				node->relion_MPI_Bcast(&helical_rise_half1, 1, MY_MPI_DOUBLE, reconstruct_rank1, MPI_COMM_WORLD);

				// When splitting the data into two random halves, perform two reconstructions in parallel: one for each subset
				if (do_split_random_halves)
				{
					int reconstruct_rank2 = 2 * (ith_recons % ( (node->size - 1)/2 ) ) + 2;
					node->relion_MPI_Bcast(&helical_twist_half2, 1, MY_MPI_DOUBLE, reconstruct_rank2, MPI_COMM_WORLD);
					node->relion_MPI_Bcast(&helical_rise_half2, 1, MY_MPI_DOUBLE, reconstruct_rank2, MPI_COMM_WORLD);
				}
			}
		}
	}
#ifdef TIMING
		timer.toc(TIMING_RECONS);
#endif

	if (node->isMaster())
	{
		// If we're realinging movies then the master needs to update the priors
		// These will then be re-distributed to the slaves in the expectation step (through get/setMetadataSubset)
		// Switch this off, as we're only doing a single round of movie processing anyway!
		//if (do_realign_movies)
		//	updatePriorsForMovieFrames();

		// The master also updates the changes in hidden variables
		updateOverallChangesInHiddenVariables();
	}
	else
	{
		// Now do the maximisation of all other parameters (and calculate the tau2_class-spectra of the reconstructions
		// The lazy master never does this: it only handles metadata and does not have the weighted sums
		maximizationOtherParameters();
	}

	// The master broadcasts the changes in hidden variables to all other nodes
	node->relion_MPI_Bcast(&current_changes_optimal_classes, 1, MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);
	node->relion_MPI_Bcast(&current_changes_optimal_orientations, 1, MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);
	node->relion_MPI_Bcast(&current_changes_optimal_offsets, 1, MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);
	node->relion_MPI_Bcast(&nr_iter_wo_large_hidden_variable_changes, 1, MPI_INT, 0, MPI_COMM_WORLD);
	node->relion_MPI_Bcast(&smallest_changes_optimal_classes, 1, MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);
	node->relion_MPI_Bcast(&smallest_changes_optimal_offsets, 1, MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);
	node->relion_MPI_Bcast(&smallest_changes_optimal_orientations, 1, MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if (verb > 0)
		progress_bar(mymodel.nr_classes);

	// Aug05,2015 - Shaoda, Output for helical refinement
	if ( (verb > 0) && (do_helical_refine) )
	{
		std::cout << " Helical parameters (before refinement): twist = " << helical_twist
				<< " degrees, rise = " << helical_rise << " Angstroms ("
				<< (helical_rise / mymodel.pixel_size) << " pixels), asymmetrical units = "
				<< nr_helical_asu << ", central Z proportion = "
				<< helical_central_proportion << ", spherical diameter = "
				<< particle_diameter << " Angstroms ("
				<< (particle_diameter / mymodel.pixel_size) << " pixels), tube inner diameter = "
				<< helical_mask_tube_inner_diameter << " Angstroms ("
				<< (helical_mask_tube_inner_diameter / mymodel.pixel_size) << " pixels), tube outer diameter = "
				<< helical_mask_tube_outer_diameter << " Angstroms ("
				<< (helical_mask_tube_outer_diameter / mymodel.pixel_size) << " pixels), edge width (mask) = "
				<< width_mask_edge << " pixels." << std::endl;
	}
	if ( (iter > 1) && (do_auto_refine) && (do_helical_refine) && (do_helical_symmetry_local_refinement) )
	{
		if (verb > 0)
		{
			RFLOAT nr_units_min = 2.;
			RFLOAT twist_deg_min, twist_deg_max, rise_A_min, rise_A_max, rise_A_threshold;

			twist_deg_min = helical_twist * (1. - helical_symmetry_local_refinement_max_dev);
			twist_deg_max = helical_twist * (1. + helical_symmetry_local_refinement_max_dev);
			if (fabs(twist_deg_min) < 0.01)
				twist_deg_min = (twist_deg_min > 0.) ? (0.01) : (-0.01);
			if (fabs(twist_deg_max) > 179.99)
				twist_deg_max = (twist_deg_max > 0.) ? (179.99) : (-179.99);

			rise_A_min = helical_rise * (1. - helical_symmetry_local_refinement_max_dev);
			rise_A_max = helical_rise * (1. + helical_symmetry_local_refinement_max_dev);
			rise_A_threshold = get_rise_A_max(mymodel.ori_size, mymodel.pixel_size, helical_central_proportion, nr_units_min);
			if ( (rise_A_min / mymodel.pixel_size) < 0.001)
				rise_A_min = 0.001 * mymodel.pixel_size;
			if (rise_A_max > rise_A_threshold)
			{
				rise_A_max = rise_A_threshold;
				std::cout << " WARNING: Increase of lenZ_percentage will enable larger range for helical rise searches." << std::endl;
			}

			std::cout << " Helical twist searched from "
					<< twist_deg_min << " to " << twist_deg_max << " degrees, rise from "
					<< rise_A_min << " to " << rise_A_max << " Angstroms ("
					<< (rise_A_min / mymodel.pixel_size) << " to " << (rise_A_max / mymodel.pixel_size)
					<< " pixels)." << std::endl;

			if (do_split_random_halves)
			{
				std::cout << " Refined helical parameters (half1): twist = "
						<< helical_twist_half1 << " degrees, rise = "
						<< helical_rise_half1 << " Angstroms ("
						<< (helical_rise_half1 / mymodel.pixel_size) << " pixels)" << std::endl;
				std::cout << " Refined helical parameters (half2): twist = "
						<< helical_twist_half2 << " degrees, rise = "
						<< helical_rise_half2 << " Angstroms ("
						<< (helical_rise_half2 / mymodel.pixel_size) << " pixels)" << std::endl;
			}
			else
			{
				std::cout << " Refined helical parameters: twist = "
						<< helical_twist_half1 << " degrees, rise = "
						<< helical_rise_half1 << " Angstroms ("
						<< (helical_rise_half1 / mymodel.pixel_size) << " pixels)" << std::endl;
			}
		}

		// No need to average/modify helical parameters if nr_asu==1 (when no helical symmetry is applied in Fourier space)
		if (nr_helical_asu > 1)
		{
			if (do_split_random_halves)
			{
				if (verb > 0)
				{
					if ( (fabs((helical_twist_half1 / helical_twist_half2) - 1.) > 0.01)
							|| (fabs((helical_rise_half1 / helical_rise_half2) - 1.) > 0.01) )
						std::cout << " WARNING: Refined helical parameters from the two halves are different ( >1% apart)! Averaged parameters might be wrong! Please make sure that the data set is large enough and does not contain multiple symmetries! You probably need not worry if they converge during the next few iterations." << std::endl;
				}
				helical_twist = (helical_twist_half1 + helical_twist_half2) / 2.;
				helical_rise = (helical_rise_half1 + helical_rise_half2) / 2.;
				if (verb > 0)
				{
					std::cout << " Averaged helical parameters will be used for the next iteration: twist = "
							<< helical_twist << " degrees, rise = "
							<< helical_rise << " Angstroms ("
							<< (helical_rise / mymodel.pixel_size) << " pixels)" << std::endl;
				}
			}
			else
			{
				helical_twist = helical_twist_half1;
				helical_rise = helical_rise_half1;
			}
		}
	}

#ifdef DEBUG
	std::cerr << "MlOptimiserMpi::maximization: done" << std::endl;
#endif

}

void MlOptimiserMpi::joinTwoHalvesAtLowResolution()
{
#ifdef DEBUG
	std::cerr << "MlOptimiserMpi::joinTwoHalvesAtLowResolution: Entering " << std::endl;
#endif

	if (!do_split_random_halves)
		REPORT_ERROR("BUG: you should not be in MlOptimiserMpi::joinTwoHalvesAtLowResolution!");

	// Loop over all classes (this will be just one class for now...)
	RFLOAT myres = XMIPP_MAX(low_resol_join_halves, 1./mymodel.current_resolution);
	int lowres_r_max = CEIL(mymodel.ori_size * mymodel.pixel_size / myres);

	for (int iclass = 0; iclass < mymodel.nr_classes; iclass++ )
	{
		if (node->rank == 1 || node->rank == 2)
		{
			MultidimArray<Complex > lowres_data;
			MultidimArray<RFLOAT > lowres_weight;
			wsum_model.BPref[iclass].getLowResDataAndWeight(lowres_data, lowres_weight, lowres_r_max);

			if (node->rank == 2)
			{
				MPI_Status status;

				// The second slave sends its lowres_data and lowres_weight to the first slave
				node->relion_MPI_Send(MULTIDIM_ARRAY(lowres_data), 2*MULTIDIM_SIZE(lowres_data), MY_MPI_DOUBLE, 1, MPITAG_IMAGE, MPI_COMM_WORLD);
				node->relion_MPI_Send(MULTIDIM_ARRAY(lowres_weight), MULTIDIM_SIZE(lowres_weight), MY_MPI_DOUBLE, 1, MPITAG_RFLOAT, MPI_COMM_WORLD);

				// Now the first slave is calculating the average....

				// Then the second slave receives the average back from the first slave
				node->relion_MPI_Recv(MULTIDIM_ARRAY(lowres_data), 2*MULTIDIM_SIZE(lowres_data), MY_MPI_DOUBLE, 1, MPITAG_IMAGE, MPI_COMM_WORLD, status);
				node->relion_MPI_Recv(MULTIDIM_ARRAY(lowres_weight), MULTIDIM_SIZE(lowres_weight), MY_MPI_DOUBLE, 1, MPITAG_RFLOAT, MPI_COMM_WORLD, status);


			}
			else if (node->rank == 1)
			{

				std::cout << " Averaging half-reconstructions up to " << myres << " Angstrom resolution to prevent diverging orientations ..." << std::endl;
				std::cout << " Note that only for higher resolutions the FSC-values are according to the gold-standard!" << std::endl;
				MPI_Status status;
				MultidimArray<Complex > lowres_data_half2;
				MultidimArray<RFLOAT > lowres_weight_half2;
				lowres_data_half2.resize(lowres_data);
				lowres_weight_half2.resize(lowres_weight);
#ifdef DEBUG
				std::cerr << "AAArank=1 lowresdata: "; lowres_data.printShape();
				std::cerr << "AAArank=1 lowresdata_half2: "; lowres_data_half2.printShape();
#endif
				// The first slave receives the average from the second slave
				node->relion_MPI_Recv(MULTIDIM_ARRAY(lowres_data_half2), 2*MULTIDIM_SIZE(lowres_data_half2), MY_MPI_DOUBLE, 2, MPITAG_IMAGE, MPI_COMM_WORLD, status);
				node->relion_MPI_Recv(MULTIDIM_ARRAY(lowres_weight_half2), MULTIDIM_SIZE(lowres_weight_half2), MY_MPI_DOUBLE, 2, MPITAG_RFLOAT, MPI_COMM_WORLD, status);

				// The first slave calculates the average of the two lowres_data and lowres_weight arrays
#ifdef DEBUG
				std::cerr << "BBBrank=1 lowresdata: "; lowres_data.printShape();
				std::cerr << "BBBrank=1 lowresdata_half2: "; lowres_data_half2.printShape();
#endif
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(lowres_data)
				{
					DIRECT_MULTIDIM_ELEM(lowres_data, n)   += DIRECT_MULTIDIM_ELEM(lowres_data_half2, n) ;
					DIRECT_MULTIDIM_ELEM(lowres_data, n)   /= 2.;
					DIRECT_MULTIDIM_ELEM(lowres_weight, n) += DIRECT_MULTIDIM_ELEM(lowres_weight_half2, n) ;
					DIRECT_MULTIDIM_ELEM(lowres_weight, n) /= 2.;
				}

				// The first slave sends the average lowres_data and lowres_weight also back to the second slave
				node->relion_MPI_Send(MULTIDIM_ARRAY(lowres_data), 2*MULTIDIM_SIZE(lowres_data), MY_MPI_DOUBLE, 2, MPITAG_IMAGE, MPI_COMM_WORLD);
				node->relion_MPI_Send(MULTIDIM_ARRAY(lowres_weight), MULTIDIM_SIZE(lowres_weight), MY_MPI_DOUBLE, 2, MPITAG_RFLOAT, MPI_COMM_WORLD);

			}

			// Now that both slaves have the average lowres arrays, set them back into the backprojector
			wsum_model.BPref[iclass].setLowResDataAndWeight(lowres_data, lowres_weight, lowres_r_max);
		}
	}


#ifdef DEBUG
	std::cerr << "MlOptimiserMpi::joinTwoHalvesAtLowResolution: done" << std::endl;
#endif

}

void MlOptimiserMpi::reconstructUnregularisedMapAndCalculateSolventCorrectedFSC()
{

	if (fn_mask == "")
		return;

	if (mymodel.ref_dim == 3 && (node->rank == 1 || (do_split_random_halves && node->rank == 2) ) )
	{
		Image<RFLOAT> Iunreg;
		MultidimArray<RFLOAT> dummy;
		FileName fn_root;
		if (iter > -1)
			fn_root.compose(fn_out+"_it", iter, "", 3);
		else
			fn_root = fn_out;
		fn_root += "_half" + integerToString(node->rank);;
		fn_root.compose(fn_root+"_class", 1, "", 3);

		// This only works for do_auto_refine, so iclass=0
		BackProjector BPextra(wsum_model.BPref[0]);

		BPextra.reconstruct(Iunreg(), gridding_nr_iter, false, 1., dummy, dummy, dummy, dummy, 1., false, true, nr_threads, -1);

		// Update header information
		RFLOAT avg, stddev, minval, maxval;
	    Iunreg().setXmippOrigin();
		Iunreg().computeStats(avg, stddev, minval, maxval);
	    Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
	    Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
	    Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
	    Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
		Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, mymodel.pixel_size);
		Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, mymodel.pixel_size);
		Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, mymodel.pixel_size);
		// And write the resulting model to disc
		Iunreg.write(fn_root+"_unfil.mrc");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	if (node->rank == 0) // Let's do this on the master (hopefully it has more memory)
	{
		std::cout << " Calculating solvent-corrected gold-standard FSC ..."<< std::endl;

		// Read in the half-reconstruction from rank2 and perform the postprocessing-like FSC correction
		Image<RFLOAT> Iunreg1, Iunreg2;
		FileName fn_root1, fn_root2;
		if (iter > -1)
			fn_root1.compose(fn_out+"_it", iter, "", 3);
		else
			fn_root1 = fn_out;
		fn_root2.compose(fn_root1+"_half2_class", 1, "", 3);
		fn_root1.compose(fn_root1+"_half1_class", 1, "", 3);
		fn_root1 += "_unfil.mrc";
		fn_root2 += "_unfil.mrc";
		Iunreg1.read(fn_root1);
		Iunreg2.read(fn_root2);
		Iunreg1().setXmippOrigin();
		Iunreg2().setXmippOrigin();

		// Now do phase-randomisation FSC-correction for the solvent mask
		MultidimArray<RFLOAT> fsc_unmasked, fsc_masked, fsc_random_masked, fsc_true;

		// Calculate FSC of the unmasked maps
		getFSC(Iunreg1(), Iunreg2(), fsc_unmasked);

		Image<RFLOAT> Imask;
		Imask.read(fn_mask);
		Imask().setXmippOrigin();
		Iunreg1() *= Imask();
		Iunreg2() *= Imask();
		getFSC(Iunreg1(), Iunreg2(), fsc_masked);

		// To save memory re-read the same input maps again and randomize phases before masking
		Iunreg1.read(fn_root1);
		Iunreg2.read(fn_root2);
		Iunreg1().setXmippOrigin();
		Iunreg2().setXmippOrigin();

		// Check at which resolution shell the FSC drops below 0.8
		int randomize_at = -1;
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_unmasked)
		{
			if (i > 0 && DIRECT_A1D_ELEM(fsc_unmasked, i) < 0.8)
			{
				randomize_at = i;
				break;
			}
		}
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left << "  + randomize phases beyond: "; std::cout << XSIZE(Iunreg1())* mymodel.pixel_size / randomize_at << " Angstroms" << std::endl;
		}
		if (randomize_at > 0)
		{
			randomizePhasesBeyond(Iunreg1(), randomize_at);
			randomizePhasesBeyond(Iunreg2(), randomize_at);
			// Mask randomized phases maps and calculated fsc_random_masked
			Iunreg1() *= Imask();
			Iunreg2() *= Imask();
			getFSC(Iunreg1(), Iunreg2(), fsc_random_masked);

		}
		else
			REPORT_ERROR("Postprocessing::run ERROR: FSC curve never drops below randomize_fsc_at.  You may want to check your mask.");

		// Now that we have fsc_masked and fsc_random_masked, calculate fsc_true according to Richard's formula
		// FSC_true = FSC_t - FSC_n / ( )
		fsc_true.resize(fsc_masked);
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
		{
					// 29jan2015: let's move this 2 shells upwards, because of small artefacts near the resolution of randomisation!
			if (i < randomize_at + 2)
			{
				DIRECT_A1D_ELEM(fsc_true, i) = DIRECT_A1D_ELEM(fsc_masked, i);
			}
			else
			{
				RFLOAT fsct = DIRECT_A1D_ELEM(fsc_masked, i);
				RFLOAT fscn = DIRECT_A1D_ELEM(fsc_random_masked, i);
				if (fscn > fsct)
					DIRECT_A1D_ELEM(fsc_true, i) = 0.;
				else
					DIRECT_A1D_ELEM(fsc_true, i) = (fsct - fscn) / (1. - fscn);
			}
		}
		mymodel.fsc_halves_class = fsc_true;

	}

	// Now the master sends the fsc curve to everyone else
	node->relion_MPI_Bcast(MULTIDIM_ARRAY(mymodel.fsc_halves_class), MULTIDIM_SIZE(mymodel.fsc_halves_class), MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);


}

void MlOptimiserMpi::writeTemporaryDataAndWeightArrays()
{

	if ( (node->rank == 1 || (do_split_random_halves && node->rank == 2) ) )
	{
		Image<RFLOAT> It;
		FileName fn_root = fn_out + "_half" + integerToString(node->rank);;

		// Write out temporary arrays for all classes
		for (int iclass = 0; iclass < mymodel.nr_bodies * mymodel.nr_classes; iclass++)
		{
			FileName fn_tmp;
			fn_tmp.compose(fn_root+"_class", iclass+1, "", 3);
			if (mymodel.pdf_class[iclass] > 0.)
			{
				It().resize(wsum_model.BPref[iclass].data);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(It())
				{
					DIRECT_MULTIDIM_ELEM(It(), n) = (DIRECT_MULTIDIM_ELEM(wsum_model.BPref[iclass].data, n)).real;
				}
				It.write(fn_tmp+"_data_real.mrc");
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(It())
				{
					DIRECT_MULTIDIM_ELEM(It(), n) = (DIRECT_MULTIDIM_ELEM(wsum_model.BPref[iclass].data, n)).imag;
				}
				It.write(fn_tmp+"_data_imag.mrc");
				It()=wsum_model.BPref[iclass].weight;
				It.write(fn_tmp+"_weight.mrc");
			}
		}
    }

}

void MlOptimiserMpi::readTemporaryDataAndWeightArraysAndReconstruct(int iclass, int ihalf)
{
	MultidimArray<RFLOAT> dummy;
	Image<RFLOAT> Iunreg, Itmp;
	FileName fn_root = fn_out + "_half" + integerToString(ihalf);;
	fn_root.compose(fn_root+"_class", iclass+1, "", 3);

	// Read temporary arrays back in
	Itmp.read(fn_root+"_data_real.mrc");
	Itmp().setXmippOrigin();
	Itmp().xinit=0;
	if (!Itmp().sameShape(wsum_model.BPref[iclass].data))
	{
		wsum_model.BPref[iclass].data.printShape(std::cerr);
		Itmp().printShape(std::cerr);
		REPORT_ERROR("Incompatible size of "+fn_root+"_data_real.mrc");
	}
	FOR_ALL_ELEMENTS_IN_ARRAY3D(Itmp())
	{
		A3D_ELEM(wsum_model.BPref[iclass].data, k, i, j).real = A3D_ELEM(Itmp(), k, i, j);
	}

	Itmp.read(fn_root+"_data_imag.mrc");
	Itmp().setXmippOrigin();
	Itmp().xinit=0;
	if (!Itmp().sameShape(wsum_model.BPref[iclass].data))
	{
		wsum_model.BPref[iclass].data.printShape(std::cerr);
		Itmp().printShape(std::cerr);
		REPORT_ERROR("Incompatible size of "+fn_root+"_data_imag.mrc");
	}
	FOR_ALL_ELEMENTS_IN_ARRAY3D(Itmp())
	{
		A3D_ELEM(wsum_model.BPref[iclass].data, k, i, j).imag = A3D_ELEM(Itmp(), k, i, j);
	}

	Itmp.read(fn_root+"_weight.mrc");
	Itmp().setXmippOrigin();
	Itmp().xinit=0;
	if (!Itmp().sameShape(wsum_model.BPref[iclass].weight))
	{
		wsum_model.BPref[iclass].weight.printShape(std::cerr);
		Itmp().printShape(std::cerr);
		REPORT_ERROR("Incompatible size of "+fn_root+"_weight.mrc");
	}
	FOR_ALL_ELEMENTS_IN_ARRAY3D(Itmp())
	{
		A3D_ELEM(wsum_model.BPref[iclass].weight, k, i, j) = A3D_ELEM(Itmp(), k, i, j);
	}

	// Now perform the unregularized reconstruction
	wsum_model.BPref[iclass].reconstruct(Iunreg(), gridding_nr_iter, false, 1., dummy, dummy, dummy, dummy, 1., false, true, nr_threads, -1);

	// Update header information
	RFLOAT avg, stddev, minval, maxval;
	Iunreg().computeStats(avg, stddev, minval, maxval);
	Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
	Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
	Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
	Iunreg.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
	Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, mymodel.pixel_size);
	Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, mymodel.pixel_size);
	Iunreg.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, mymodel.pixel_size);
	// And write the resulting model to disc
	Iunreg.write(fn_root+"_unfil.mrc");


	// remove temporary arrays from the disc
	remove((fn_root+"_data_real.mrc").c_str());
	remove((fn_root+"_data_imag.mrc").c_str());
	remove((fn_root+"_weight.mrc").c_str());

}

void MlOptimiserMpi::compareTwoHalves()
{
#ifdef DEBUG
	std::cerr << "MlOptimiserMpi::compareTwoHalves: Entering " << std::endl;
#endif

	if (!do_split_random_halves)
		REPORT_ERROR("ERROR: you should not be in MlOptimiserMpi::compareTwoHalves if !do_split_random_halves");

	if (mymodel.nr_classes > 1)
		REPORT_ERROR("ERROR: you should not be in MlOptimiserMpi::compareTwoHalves if mymodel.nr_classes > 1");

	// Only do gold-standard FSC comparisons for single-class refinements
	int iclass = 0;
	// The first two slaves calculate the sum of the downsampled average of all bodies
	if (node->rank == 1 || node->rank == 2)
	{
		// Sum over all bodies
		MultidimArray<Complex > avg1;
		for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++ )
		{
			MultidimArray<Complex > tmp1;
			wsum_model.BPref[ibody].getDownsampledAverage(tmp1);

			if (mymodel.nr_bodies > 1)
			{
				// 19may2015 Shift each downsampled avg in the FourierTransform to its original COM!!
				// shiftImageInFourierTransform, but tmp1 is already centered (ie not in FFTW-arrangement)
				RFLOAT dotp, a, b, c, d, ac, bd, ab_cd, x, y, z;
				RFLOAT xshift = -XX(mymodel.com_bodies[ibody])/(RFLOAT)mymodel.ori_size;
				RFLOAT yshift = -YY(mymodel.com_bodies[ibody])/(RFLOAT)mymodel.ori_size;
				RFLOAT zshift = -ZZ(mymodel.com_bodies[ibody])/(RFLOAT)mymodel.ori_size;
				FOR_ALL_ELEMENTS_IN_ARRAY3D(tmp1)
				{
					dotp = 2 * PI * (j * xshift + i * yshift + k * zshift);
					a = cos(dotp);
					b = sin(dotp);
					c = A3D_ELEM(tmp1, k, i, j).real;
					d = A3D_ELEM(tmp1, k, i, j).imag;
					ac = a * c;
					bd = b * d;
					ab_cd = (a + b) * (c + d);
					A3D_ELEM(tmp1, k, i, j) = Complex(ac - bd, ab_cd - ac - bd);
				}
			}

			if (ibody == 0)
				avg1 = tmp1;
			else
				avg1 += tmp1;
		}

//#define DEBUG_FSC
#ifdef DEBUG_FSC
		MultidimArray<Complex > avg;
		MultidimArray<RFLOAT> Mavg;
		if (mymodel.ref_dim == 2)
			Mavg.resize(mymodel.ori_size, mymodel.ori_size);
		else
			Mavg.resize(mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);

		FourierTransformer transformer_debug;
		transformer_debug.setReal(Mavg);
		transformer_debug.getFourierAlias(avg);
		wsum_model.BPref[0].decenter(avg1, avg, wsum_model.BPref[0].r_max * wsum_model.BPref[0].r_max);
		transformer_debug.inverseFourierTransform();
		FileName fnt;
		fnt.compose("downsampled_avg_half",node->rank,"spi");
		Image<RFLOAT> It;
		CenterFFT(Mavg, true);
		It()=Mavg;
		It.write(fnt);
#endif

		if (node->rank == 2)
		{
			// The second slave sends its average to the first slave
			node->relion_MPI_Send(MULTIDIM_ARRAY(avg1), 2*MULTIDIM_SIZE(avg1), MY_MPI_DOUBLE, 1, MPITAG_IMAGE, MPI_COMM_WORLD);
		}
		else if (node->rank == 1)
		{

			std::cout << " Calculating gold-standard FSC ..."<< std::endl;
			// The first slave receives the average from the second slave and calculates the FSC between them
			MPI_Status status;
			MultidimArray<Complex > avg2;
			avg2.resize(avg1);
			node->relion_MPI_Recv(MULTIDIM_ARRAY(avg2), 2*MULTIDIM_SIZE(avg2), MY_MPI_DOUBLE, 2, MPITAG_IMAGE, MPI_COMM_WORLD, status);
			wsum_model.BPref[iclass].calculateDownSampledFourierShellCorrelation(avg1, avg2, mymodel.fsc_halves_class);
		}

	}

	// Now slave 1 sends the fsc curve to everyone else
	node->relion_MPI_Bcast(MULTIDIM_ARRAY(mymodel.fsc_halves_class), MULTIDIM_SIZE(mymodel.fsc_halves_class), MY_MPI_DOUBLE, 1, MPI_COMM_WORLD);

#ifdef DEBUG
	std::cerr << "MlOptimiserMpi::compareTwoHalves: done" << std::endl;
#endif
}

void MlOptimiserMpi::iterate()
{
#ifdef TIMING
	// MPI-specific timing stuff goes here...
	TIMING_MPIWAIT= timer.setNew("mpiWaitEndOfExpectation");
	TIMING_MPICOMBINEDISC= timer.setNew("mpiCombineThroughDisc");
	TIMING_MPICOMBINENETW= timer.setNew("mpiCombineThroughNetwork");
	TIMING_MPISLAVEWORK= timer.setNew("mpiSlaveWorking");
	TIMING_MPISLAVEWAIT1= timer.setNew("mpiSlaveWaiting1");
	TIMING_MPISLAVEWAIT2= timer.setNew("mpiSlaveWaiting2");
	TIMING_MPISLAVEWAIT3= timer.setNew("mpiSlaveWaiting3");
#endif


	// Launch threads etc.
	MlOptimiser::iterateSetup();

	// Initialize the current resolution
	updateCurrentResolution();

	for (iter = iter + 1; iter <= nr_iter; iter++)
    {
#ifdef TIMING
		timer.tic(TIMING_EXP);
#endif

		// Nobody can start the next iteration until everyone has finished
		MPI_Barrier(MPI_COMM_WORLD);

		if (do_auto_refine && node->rank == 1)
			printConvergenceStats();

		expectation();

		MPI_Barrier(MPI_COMM_WORLD);

		if (do_skip_maximization)
		{
			// Only write data.star file and break from the iteration loop
			if (node->isMaster())
			{
				// The master only writes the data file (he's the only one who has and manages these data!)
				iter = -1; // write output file without iteration number
				MlOptimiser::write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, node->rank);
				std::cout << " Auto-refine: Skipping maximization step, so stopping now... " << std::endl;
			}
			break;
		}

		// Now combine all weighted sums
		// Leave the option ot both for a while. Then, if there are no problems with the system via files keep that one and remove the MPI version from the code
		if (combine_weights_thru_disc)
			combineAllWeightedSumsViaFile();
		else
			combineAllWeightedSums();

		MPI_Barrier(MPI_COMM_WORLD);

		// Sjors & Shaoda Apr 2015
		// This function does enforceHermitianSymmetry, applyHelicalSymmetry and applyPointGroupSymmetry sequentially.
		// First it enforces Hermitian symmetry to the back-projected Fourier 3D matrix.
		// Then helical symmetry is applied in Fourier space. It does rise and twist for all asymmetrical units in Fourier space.
		// Finally it applies point group symmetry (such as Cn, ...).
		// DEBUG
		if ( (verb > 0) && (node->isMaster()) )
		{
			if (do_helical_refine)
			{
				if (nr_helical_asu > 1)
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

		// Write out data and weight arrays to disc in order to also do an unregularized reconstruction
		if (do_auto_refine && has_converged)
			writeTemporaryDataAndWeightArrays();

		// Inside iterative refinement: do FSC-calculation BEFORE the solvent flattening, otherwise over-estimation of resolution
		// anyway, now that this is done inside BPref, there would be no other way...
		if (do_split_random_halves)
		{

			// For asymmetric molecules, join 2 half-reconstructions at the lowest resolutions to prevent them from diverging orientations
			if (low_resol_join_halves > 0.)
				joinTwoHalvesAtLowResolution();

			// Sjors 27-oct-2015
			// Calculate gold-standard FSC curve
			if (do_phase_random_fsc && fn_mask != "None")
				reconstructUnregularisedMapAndCalculateSolventCorrectedFSC();
			else
				compareTwoHalves();

			// For automated sampling procedure
			if (!node->isMaster()) // the master does not have the correct mymodel.current_size, it only handles metadata!
			{
				// Check that incr_size is at least the number of shells as between FSC=0.5 and FSC=0.143
				int fsc05   = -1;
				int fsc0143 = -1;
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(mymodel.fsc_halves_class)
				{
					if (DIRECT_A1D_ELEM(mymodel.fsc_halves_class, i) < 0.5 && fsc05 < 0)
						fsc05 = i;
					if (DIRECT_A1D_ELEM(mymodel.fsc_halves_class, i) < 0.143 && fsc0143 < 0)
						fsc0143 = i;
				}
				// At least fsc05 - fsc0143 + 5 shells as incr_size
				incr_size = XMIPP_MAX(incr_size, fsc0143 - fsc05 + 5);
				has_high_fsc_at_limit = (DIRECT_A1D_ELEM(mymodel.fsc_halves_class, mymodel.current_size/2 - 1) > 0.2);
			}

			// Upon convergence join the two random halves
			if (do_join_random_halves || do_always_join_random_halves)
			{
				if (combine_weights_thru_disc)
					combineWeightedSumsTwoRandomHalvesViaFile();
				else
					combineWeightedSumsTwoRandomHalves();

			}
		}

#ifdef TIMING
		timer.toc(TIMING_EXP);
		timer.tic(TIMING_MAX);
#endif

		maximization();

		// Make sure all nodes have the same resolution, set the data_vs_prior array from half1 also for half2
		if (do_split_random_halves)
			for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
				node->relion_MPI_Bcast(MULTIDIM_ARRAY(mymodel.data_vs_prior_class[iclass]), MULTIDIM_SIZE(mymodel.data_vs_prior_class[iclass]), MY_MPI_DOUBLE, 1, MPI_COMM_WORLD);

#ifdef TIMING
		timer.toc(TIMING_MAX);
#endif

		MPI_Barrier(MPI_COMM_WORLD);

		// Mask the reconstructions to get rid of noisy solvent areas
		// Skip masking upon convergence (for validation purposes)
		if (do_solvent && !has_converged)
			solventFlatten();

		// Re-calculate the current resolution, do this before writing to get the correct values in the output files
		updateCurrentResolution();

		// If we are joining random halves, then do not write an optimiser file so that it cannot be restarted!
		bool do_write_optimiser = !do_join_random_halves;
		// Write out final map without iteration number in the filename
		if (do_join_random_halves)
			iter = -1;

		if (node->rank == 1 || (do_split_random_halves && !do_join_random_halves && node->rank == 2) )
			//Only the first_slave of each subset writes model to disc (do not write the data.star file, only master will do this)
			MlOptimiser::write(DO_WRITE_SAMPLING, DONT_WRITE_DATA, do_write_optimiser, DO_WRITE_MODEL, node->rank);
		else if (node->isMaster())
			// The master only writes the data file (he's the only one who has and manages these data!)
			MlOptimiser::write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, node->rank);

		if (do_auto_refine && has_converged)
		{
			if (verb > 0)
			{
				std::cout << " Auto-refine: Refinement has converged, stopping now... " << std::endl;
				std::cout << " Auto-refine: + Final reconstruction from all particles is saved as: " <<  fn_out << "_class001.mrc" << std::endl;
				std::cout << " Auto-refine: + Final model parameters are stored in: " << fn_out << "_model.star" << std::endl;
				std::cout << " Auto-refine: + Final data parameters are stored in: " << fn_out << "_data.star" << std::endl;
				std::cout << " Auto-refine: + Final resolution (without masking) is: " << 1./mymodel.current_resolution << std::endl;
				if (acc_rot < 10.)
					std::cout << " Auto-refine: + But you may want to run relion_postprocess to mask the unfil.mrc maps and calculate a higher resolution FSC" << std::endl;
				else
				{
					std::cout << " Auto-refine: + WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles!" << std::endl;
					std::cout << " Auto-refine: + WARNING: This has been observed to lead to spurious FSC curves, so be VERY wary of inflated resolution estimates..." << std::endl;
					std::cout << " Auto-refine: + WARNING: You most probably do NOT want to publish these results!" << std::endl;
					std::cout << " Auto-refine: + WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
				}
				if (do_use_reconstruct_images)
					std::cout << " Auto-refine: + Used rlnReconstructImageName images for final reconstruction. Ignore filtered map, and only assess the unfiltered half-reconstructions!" << std::endl;
			}
			break;
		}

		// Check whether we have converged by now
		// Master does not have all info, only slaves do this and then broadcast has_converged so also master has it
		// If we have, set do_join_random_halves and do_use_all_data for the next iteration
		if (!node->isMaster() && do_auto_refine)
			checkConvergence();
		node->relion_MPI_Bcast(&has_converged, 1, MPI_INT, 1, MPI_COMM_WORLD);

#ifdef TIMING
		// Only first slave prints it timing information
		if (node->rank == 1)
			timer.printTimes(false);
#endif

    }

	// delete threads etc.
	MlOptimiser::iterateWrapUp();
	MPI_Barrier(MPI_COMM_WORLD);

}


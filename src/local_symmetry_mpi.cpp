#include "src/local_symmetry_mpi.h"

void local_symmetry_parameters_mpi::read(int argc, char **argv)
{
    // Define a new MpiNode
    node = new MpiNode(argc, argv);

    // First read in non-parallelisation-dependent variables
    local_symmetry_parameters::read(argc, argv);

    // Don't put any output to screen for mpi slaves
    verb = (node->isMaster()) ? 1 : 0;

    // Possibly also read parallelisation-dependent variables here

    // Print out MPI info
	printMpiNodesMachineNames(*node);
}

void local_symmetry_parameters_mpi::run()
{
	Image<RFLOAT> map, mask;
	Matrix1D<RFLOAT> op_search_ranges, op;
	std::vector<FileName> fn_mask_list;
	std::vector<std::vector<Matrix1D<RFLOAT> > > op_list;
	std::vector<Matrix1D<RFLOAT> > op_samplings, op_samplings_batch;
	MultidimArray<RFLOAT> op_samplings_batch_packed;
	int nr_masks = 0, nr_ops = 0, nr_total_samplings = 0;
	long int olddim = 0, newdim = 0;
	RFLOAT mask_val = 0.;
	bool use_healpix = true;

	map.clear(); mask.clear();
	op_search_ranges.clear(); op.clear();
	fn_mask_list.clear();
	op_list.clear();
	op_samplings.clear(); op_samplings_batch.clear();
	op_samplings_batch_packed.clear();

	if ( (do_apply_local_symmetry) || (do_duplicate_local_symmetry) || (do_txt2rln) || (do_debug) || (!do_local_search_local_symmetry_ops) )
		REPORT_ERROR("ERROR: Please specify '--search' as the only option! For other options use non-parallel version (without '_mpi') instead!");

	if (node->isMaster())
	{
		bool do_verb = true;
		FileName fn_parsed;

#ifdef DEBUG
		std::cout << " relion_localsym_mpi is running ..." << std::endl;
#endif

		if (angpix_image < 0.001)
			REPORT_ERROR("Invalid pixel size!");

		// Master parses and reads mask info file
		if (fn_info_in.getExtension() == "star")
		{
			readRelionFormatMasksAndOperators(fn_info_in, fn_mask_list, op_list, angpix_image, do_verb);
		}
		else
		{
			fn_parsed = fn_info_in + std::string(".") + fn_info_in_parsed_ext;
			parseDMFormatMasksAndOperators(fn_info_in, fn_parsed);
			readDMFormatMasksAndOperators(fn_parsed, fn_mask_list, op_list, angpix_image, do_verb);
		}

		// Master set total number of masks
		nr_masks = fn_mask_list.size();

		// Master reads input map
		map.read(fn_unsym);
		if ((NSIZE(map()) != 1) || (ZSIZE(map()) <= 1) || (YSIZE(map()) <= 1) || (XSIZE(map()) <= 1)
				|| (ZSIZE(map()) != YSIZE(map())) || (ZSIZE(map()) != XSIZE(map()))
				|| (ZSIZE(map()) % 2) )
			REPORT_ERROR("ERROR: Input map is not 3D cube!");
		olddim = newdim = ZSIZE(map());
		if ((binning_factor - 1.) > XMIPP_EQUAL_ACCURACY)
		{
			newdim = (long int)(ceil(RFLOAT(olddim) / binning_factor));
			if (newdim < 2)
				REPORT_ERROR("ERROR: Binning factor is too large!");
			if ((newdim + 1) < olddim) // Need rescaling
			{
				// Dimension should always be even
				if (newdim % 2)
					newdim++;
				resizeMap(map(), newdim);
			}
			else
				newdim = olddim;
			binning_factor = RFLOAT(olddim) / RFLOAT(newdim);
			if (newdim != olddim)
				std::cout << " + Rescale box size from " << olddim << " to " << newdim << ". Binning factor = " << binning_factor << std::endl;
		}
		else
			binning_factor = 1.;

#ifdef DEBUG
		std::cout << " I am master. The nxyzdim of map() is " << map().nzyxdim << std::endl;
#endif
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// Master broadcasts size of the mask
	node->relion_MPI_Bcast(&olddim, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	node->relion_MPI_Bcast(&newdim, 1, MPI_LONG, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	// Slave allocate space for MultidimArray
	// TODO: check whether the space is allocated and the map is read and successfully broadcast!!!
	if (! node->isMaster())
	{
		map.clear();
		map().initZeros(1, newdim, newdim, newdim);
#ifdef DEBUG
		if (node->rank == 1)
			std::cout << " I am node rank 1. The nxyzdim of map() is " << map().nzyxdim << std::endl;
#endif
	}
	MPI_Barrier(MPI_COMM_WORLD);

	// Master broadcasts input map to all nodes
#ifdef DEBUG
	if (node->isMaster())
		std::cout << " Master is broadcasting input map ..." << std::endl;
#endif
	node->relion_MPI_Bcast(MULTIDIM_ARRAY(map()), MULTIDIM_SIZE(map()), MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef DEBUG
	if (node->isMaster())
		std::cout << " Master has completed broadcasting." << std::endl;
#endif

	// Master gets search ranges
	if (node->isMaster())
	{
		Localsym_composeOperator(
				op_search_ranges,
				(ang_range > (XMIPP_EQUAL_ACCURACY)) ? (ang_range) : (ang_rot_range),
				(ang_range > (XMIPP_EQUAL_ACCURACY)) ? (ang_range) : (ang_tilt_range),
				(ang_range > (XMIPP_EQUAL_ACCURACY)) ? (ang_range) : (ang_psi_range),
				(offset_range > (XMIPP_EQUAL_ACCURACY)) ? (offset_range) : (offset_x_range),
				(offset_range > (XMIPP_EQUAL_ACCURACY)) ? (offset_range) : (offset_y_range),
				(offset_range > (XMIPP_EQUAL_ACCURACY)) ? (offset_range) : (offset_z_range) );
		Localsym_scaleTranslations(op_search_ranges, 1. / angpix_image);
		offset_step /= angpix_image;

		// Rescale translational search ranges and steps
		if (newdim != olddim)
		{
			Localsym_scaleTranslations(op_search_ranges, 1. / binning_factor);
			offset_step /= binning_factor;

			for (int imask = 0; imask < op_list.size(); imask++)
			{
				for (int iop = 0; iop < op_list[imask].size(); iop++)
					Localsym_scaleTranslations(op_list[imask][iop], 1. / binning_factor);
			}
		}
	}

	// Master broadcasts total number of masks
	node->relion_MPI_Bcast(&nr_masks, 1, MPI_INT, 0, MPI_COMM_WORLD);

	// All nodes loop over all masks
	for (int imask = 0; imask < nr_masks; imask++)
	{
		MPI_Barrier(MPI_COMM_WORLD);

		// Master reads the mask
		if (node->isMaster())
		{
			mask.read(fn_mask_list[imask]);
			if ((NSIZE(mask()) != 1) || (ZSIZE(mask()) <= 1) || (YSIZE(mask()) <= 1) || (XSIZE(mask()) <= 1)
					|| (ZSIZE(mask()) != YSIZE(mask())) || (ZSIZE(mask()) != XSIZE(mask()))
					|| (ZSIZE(mask()) % 2) )
				REPORT_ERROR("ERROR: Input mask is not 3D cube!");
			if (olddim != ZSIZE(mask()))
				REPORT_ERROR("ERROR: Input map and masks should have the same size!");
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mask())
			{
				mask_val = DIRECT_A3D_ELEM(mask(), k, i, j);
				if ((mask_val < -(XMIPP_EQUAL_ACCURACY)) || ((mask_val - 1.) > (XMIPP_EQUAL_ACCURACY)))
					REPORT_ERROR("ERROR: mask " + std::string(fn_mask_list[imask]) + " - values are not in range [0,1]!");
			}
			if (newdim != olddim)
			{
				resizeMap(mask(), newdim);

				// Mask values might go out of range after rescaling. Fix it if it happens
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mask())
				{
					mask_val = DIRECT_A3D_ELEM(mask(), k, i, j);
					if (mask_val < 0.)
						DIRECT_A3D_ELEM(mask(), k, i, j) = 0.;
					if (mask_val > 1.)
						DIRECT_A3D_ELEM(mask(), k, i, j) = 1.;
				}
			}
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// Slave allocate space for MultidimArray
		// TODO: check whether the space is allocated and the map is read and successfully broadcast!!!
		if (! node->isMaster())
		{
			mask.clear();
			mask().initZeros(1, newdim, newdim, newdim);
#ifdef DEBUG
			if (node->rank == 1)
				std::cout << " I am node rank 1. The nxyzdim of mask() is " << mask().nzyxdim << std::endl;
#endif
		}
		MPI_Barrier(MPI_COMM_WORLD);

		// Master broadcasts the mask to all nodes
#ifdef DEBUG
		if (node->isMaster())
			std::cout << " Master is broadcasting mask #" << (imask + 1) << " ..." << std::endl;
#endif
		node->relion_MPI_Bcast(MULTIDIM_ARRAY(mask()), MULTIDIM_SIZE(mask()), MY_MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef DEBUG
		if (node->isMaster())
			std::cout << " Master has completed broadcasting mask #" << (imask + 1) << "." << std::endl;
#endif

		// Master checks sampling rates for this mask
		if (node->isMaster())
			checkSamplingRatesForMask(mask(), ang_step, offset_step);

		// Master reads total number of operators for this mask
		if (node->isMaster())
			nr_ops = op_list[imask].size();
		MPI_Barrier(MPI_COMM_WORLD);

		// Master broadcasts total number of operators for this mask to all slaves
		node->relion_MPI_Bcast(&nr_ops, 1, MPI_INT, 0, MPI_COMM_WORLD);

		// All nodes loop over all operators of this mask
		for (int iop = 0; iop < nr_ops; iop++)
		{
			MPI_Barrier(MPI_COMM_WORLD);

			// Master gets sampling points
			if (node->isMaster())
			{
				// Check whether use Healpix
				use_healpix = useHealpixAngularSamplings(op_list[imask][iop], op_search_ranges, use_healpix_tilt_min);

				// Get sampling points
				getLocalSearchOperatorSamplings(
						op_list[imask][iop],
						op_search_ranges,
						op_samplings,
						ang_step,
						offset_step,
						use_healpix,
						true);

				if (op_samplings.size() <= (node->size))
					REPORT_ERROR("ERROR: Too few sampling points! Use non-parallel version (without '_mpi') instead!");
			}
			MPI_Barrier(MPI_COMM_WORLD);

			// Master sends the number of total samplings to all slaves
			if (node->isMaster())
				nr_total_samplings = op_samplings.size();
			node->relion_MPI_Bcast(&nr_total_samplings, 1, MPI_INT, 0, MPI_COMM_WORLD);

			// All nodes allocate space for op_samplings_batch_packed
			// Allow 10 more empty units to prevent segmentation fault
			op_samplings_batch_packed.initZeros((nr_total_samplings / (node->size)) + 10, NR_LOCALSYM_PARAMETERS);

			MPI_Barrier(MPI_COMM_WORLD);

			// Master distributes sampling points to all slaves
			long int first = 0, last = 0;
			if (node->isMaster())
			{
				for (int id_rank = (node->size) - 1; id_rank >= 0; id_rank--)
				{
					divide_equally(nr_total_samplings, node->size, id_rank, first, last);

					// Beware: YSIZE(op_samplings_batch_packed) is larger than (last - first + 1)
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(op_samplings_batch_packed)
					{
						if ( (i >= 0) && (i <= (last - first)) )
							DIRECT_A2D_ELEM(op_samplings_batch_packed, i, j) = VEC_ELEM(op_samplings[i + first], j);
					}

					// Master distributes sampling points to all slaves
					if (id_rank > 0)
						node->relion_MPI_Send(MULTIDIM_ARRAY(op_samplings_batch_packed), (last - first + 1) * NR_LOCALSYM_PARAMETERS, MY_MPI_DOUBLE, id_rank, MPITAG_LOCALSYM_SAMPLINGS_PACK, MPI_COMM_WORLD);
					// If id_rank == 0 (master), just keep op_samplings_batch_packed to the master itself
				}
			}
			else
			{
				MPI_Status status;
				// Slaves receive sampling points from master
				// Important: Slaves calculate first and last subscripts!
				divide_equally(nr_total_samplings, node->size, node->rank, first, last);
				node->relion_MPI_Recv(MULTIDIM_ARRAY(op_samplings_batch_packed), (last - first + 1) * NR_LOCALSYM_PARAMETERS, MY_MPI_DOUBLE, 0, MPITAG_LOCALSYM_SAMPLINGS_PACK, MPI_COMM_WORLD, status);
			}
			MPI_Barrier(MPI_COMM_WORLD);

			// All nodes unpack sampling points
			op_samplings_batch.clear();
			for (long int i=0; i<YSIZE(op_samplings_batch_packed); i++)
			{
				op.initZeros(NR_LOCALSYM_PARAMETERS);
				for (long int j=0; j<XSIZE(op_samplings_batch_packed); j++)
				{
					VEC_ELEM(op, j) = DIRECT_A2D_ELEM(op_samplings_batch_packed, i, j);
				}
				op_samplings_batch.push_back(op);
			}
			MPI_Barrier(MPI_COMM_WORLD);

			// All nodes calculate CC, with master profiling (DONT SORT!)
			calculateOperatorCC(map(), mask(), op_samplings_batch, false, node->isMaster());
			for (int op_id = 0; op_id < op_samplings_batch.size(); op_id++)
			{
				DIRECT_A2D_ELEM(op_samplings_batch_packed, op_id, CC_POS) = VEC_ELEM(op_samplings_batch[op_id], CC_POS);
			}
			MPI_Barrier(MPI_COMM_WORLD);

			// Slaves send their results back to master
			if (! node->isMaster())
			{
				node->relion_MPI_Send(MULTIDIM_ARRAY(op_samplings_batch_packed), (last - first + 1) * NR_LOCALSYM_PARAMETERS, MY_MPI_DOUBLE, 0, MPITAG_LOCALSYM_SAMPLINGS_PACK, MPI_COMM_WORLD);
			}
			else
			{
				MPI_Status status;

				for (int id_rank = 0; id_rank < (node->size); id_rank++)
				{
					divide_equally(op_samplings.size(), node->size, id_rank, first, last);

					// Master receives op_samplings_batch_packed from slaves
					if (id_rank > 0)
						node->relion_MPI_Recv(MULTIDIM_ARRAY(op_samplings_batch_packed), (last - first + 1) * NR_LOCALSYM_PARAMETERS, MY_MPI_DOUBLE, id_rank, MPITAG_LOCALSYM_SAMPLINGS_PACK, MPI_COMM_WORLD, status);

					// Master does something for itself if id_rank == 0
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(op_samplings_batch_packed)
					{
						// Beware: YSIZE(op_samplings_batch_packed) is larger than (last - first + 1)
						if ( (i >= 0) && (i <= (last - first)) )
							VEC_ELEM(op_samplings[i + first], CC_POS) = DIRECT_A2D_ELEM(op_samplings_batch_packed, i, CC_POS);
					}
				}
			}
			MPI_Barrier(MPI_COMM_WORLD);

			if (node->isMaster())
			{
				FileName fn_tmp, fn_searched_op_samplings;

				// For rescaled maps
				if (newdim != olddim)
				{
					for (int isamp = 0; isamp < op_samplings.size(); isamp++)
						Localsym_scaleTranslations(op_samplings[isamp], binning_factor);
				}

				// Master sorts the results
				std::stable_sort(op_samplings.begin(), op_samplings.end(), compareOperatorsByCC);

				// Master outputs the local searches results
				fn_tmp.compose("cc_mask", imask + 1, "tmp", 3);  // "cc_mask001.tmp"
				fn_tmp = fn_tmp.withoutExtension(); // "cc_mask001"
				fn_searched_op_samplings.compose(fn_tmp + "_op", iop + 1, "star", 3); // "cc_mask001_op001.star"
				writeRelionFormatLocalSearchOperatorResults(fn_searched_op_samplings, op_samplings, angpix_image);

				// Master updates this operator and do screen output
				op_list[imask][iop] = op_samplings[0];
				std::cout << " + Done! Local refined local symmetry operator: Angles = ("
						<< VEC_ELEM(op_samplings[0], AA_POS) << ", " << VEC_ELEM(op_samplings[0], BB_POS) << ", " << VEC_ELEM(op_samplings[0], GG_POS) << "), Translations = ("
						<< angpix_image * VEC_ELEM(op_samplings[0], DX_POS) << ", " << angpix_image * VEC_ELEM(op_samplings[0], DY_POS) << ", " << angpix_image * VEC_ELEM(op_samplings[0], DZ_POS) << ")." << std::endl;
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (node->isMaster())
	{
		// Master writes out new mask info file
		if (fn_info_out.getExtension() == "star")
			writeRelionFormatMasksAndOperators(fn_info_out, fn_mask_list, op_list, angpix_image);
		else
			writeDMFormatMasksAndOperators(fn_info_out, fn_mask_list, op_list, angpix_image);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

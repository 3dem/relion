static omp_lock_t global_mutex;

#include "src/ml_optimiser_mpi.h"

// ----------------------------------------------------------------------------
// -------------------- getFourierTransformsAndCtfs ---------------------------
// ----------------------------------------------------------------------------
template <class MlClass>
void getFourierTransformsAndCtfs(long int part_id,
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		MlClass *accMLO,
		AccPtrFactory ptrFactory,
		int ibody = 0
		)
{
		GTIC(accMLO->timer,"getFourierTransformsAndCtfs");
#ifdef TIMING
	if (part_id == baseMLO->exp_my_first_part_id)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_FT);
#endif

	CUSTOM_ALLOCATOR_REGION_NAME("GFTCTF");

	for (int img_id = 0; img_id < sp.nr_images; img_id++)
	{
		CTIC(accMLO->timer,"init");
		FileName fn_img;
		Image<RFLOAT> img, rec_img;
		MultidimArray<Complex > Fimg;
		MultidimArray<Complex > Faux;
		MultidimArray<RFLOAT> Fctf, FstMulti;
		Matrix2D<RFLOAT> Aori;
		Matrix1D<RFLOAT> my_projected_com(baseMLO->mymodel.data_dim), my_refined_ibody_offset(baseMLO->mymodel.data_dim);

		// Which group do I belong?
		int group_id =baseMLO->mydata.getGroupId(part_id, img_id);
		RFLOAT my_pixel_size = baseMLO->mydata.getImagePixelSize(part_id, img_id);
		// What is my optics group?
		int optics_group = baseMLO->mydata.getOpticsGroup(part_id, img_id);
		bool ctf_premultiplied = baseMLO->mydata.obsModel.getCtfPremultiplied(optics_group);

		// metadata offset for this image in the particle
		int my_metadata_offset = op.metadata_offset + img_id;

		// Get the right line in the exp_fn_img strings (also exp_fn_recimg and exp_fn_ctfs)
		int istop = 0;
		for (long int ii = baseMLO->exp_my_first_part_id; ii < part_id; ii++)
                    istop += baseMLO->mydata.numberOfImagesInParticle(part_id);
		istop += img_id;

		if (!baseMLO->mydata.getImageNameOnScratch(part_id, img_id, fn_img))
		{
			std::istringstream split(baseMLO->exp_fn_img);
			for (int i = 0; i <= my_metadata_offset; i++)
				getline(split, fn_img);
		}
		sp.current_img = fn_img;

		// Get the norm_correction
		RFLOAT normcorr = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_NORM);

		// Safeguard against gold-standard separation
		if (baseMLO->do_split_random_halves)
		{
			int halfset = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_NR_SIGN);
			if (halfset != baseMLO->my_halfset)
			{
				std::cerr << "BUG!!! halfset= " << halfset << " my_halfset= " << baseMLO->my_halfset << " part_id= " << part_id << std::endl;
				REPORT_ERROR("BUG! Mixing gold-standard separation!!!!");
			}

		}

		// Get the optimal origin offsets from the previous iteration
		// Sjors 5mar18: it is very important that my_old_offset has baseMLO->mymodel.data_dim and not just (3), as transformCartesianAndHelicalCoords will give different results!!!
		Matrix1D<RFLOAT> my_old_offset(baseMLO->mymodel.data_dim), my_prior(baseMLO->mymodel.data_dim), my_old_offset_ori;
		int icol_rot, icol_tilt, icol_psi, icol_xoff, icol_yoff, icol_zoff;
		XX(my_old_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_XOFF);
		YY(my_old_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_YOFF);
		XX(my_prior)      = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_XOFF_PRIOR);
		YY(my_prior)      = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_YOFF_PRIOR);
		// Uninitialised priors were set to 999.
		if (XX(my_prior) > 998.99 && XX(my_prior) < 999.01)
			XX(my_prior) = 0.;
		if (YY(my_prior) > 998.99 && YY(my_prior) < 999.01)
			YY(my_prior) = 0.;

		if (accMLO->dataIs3D)
		{
			ZZ(my_old_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ZOFF);
			ZZ(my_prior)      = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ZOFF_PRIOR);
			// Unitialised priors were set to 999.
			if (ZZ(my_prior) > 998.99 && ZZ(my_prior) < 999.01)
				ZZ(my_prior) = 0.;
		}

		if (baseMLO->mymodel.nr_bodies > 1)
		{

			// 17May2017: Shift image to the projected COM for this body!
			// Aori is the original transformation matrix of the consensus refinement
			Euler_angles2matrix(DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT),
					            DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT),
								DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI), Aori, false);
			my_projected_com = Aori * baseMLO->mymodel.com_bodies[ibody];
			// This will have made my_projected_com of size 3 again! resize to mymodel.data_dim
			my_projected_com.resize(baseMLO->mymodel.data_dim);

			// Subtract the projected COM offset, to position this body in the center
			// Also keep the my_old_offset in my_old_offset_ori
			my_old_offset_ori = my_old_offset;
			my_old_offset -= my_projected_com;

			// Also get refined offset for this body
			icol_xoff = 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
			icol_yoff = 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
			icol_zoff = 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
			XX(my_refined_ibody_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_xoff);
			YY(my_refined_ibody_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_yoff);
			if (baseMLO->mymodel.data_dim == 3)
				ZZ(my_refined_ibody_offset) = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_zoff);

			// For multi-body refinement: set the priors of the translations to zero (i.e. everything centred around consensus offset)
			my_prior.initZeros();
		}

		CTOC(accMLO->timer,"init");

		CTIC(accMLO->timer,"nonZeroProb");
		// Orientational priors
		if (baseMLO->mymodel.nr_bodies > 1 )
		{

			// Centre local searches around the orientation from the previous iteration, this one goes with overall sigma2_ang
			// On top of that, apply prior on the deviation from (0,0,0) with mymodel.sigma_tilt_bodies[ibody] and mymodel.sigma_psi_bodies[ibody]
			icol_rot  = 0 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
			icol_tilt = 1 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
			icol_psi  = 2 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
			RFLOAT prior_rot  = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_rot);
			RFLOAT prior_tilt = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_tilt);
			RFLOAT prior_psi =  DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_psi);
			baseMLO->sampling.selectOrientationsWithNonZeroPriorProbability(
					prior_rot, prior_tilt, prior_psi,
					sqrt(baseMLO->mymodel.sigma2_rot),
					sqrt(baseMLO->mymodel.sigma2_tilt),
					sqrt(baseMLO->mymodel.sigma2_psi),
					op.pointer_dir_nonzeroprior, op.directions_prior,
					op.pointer_psi_nonzeroprior, op.psi_prior, false, 3.,
					baseMLO->mymodel.sigma_tilt_bodies[ibody],
					baseMLO->mymodel.sigma_psi_bodies[ibody]);

		}
		else if (baseMLO->mymodel.orientational_prior_mode != NOPRIOR && !(baseMLO->do_skip_align ||baseMLO-> do_skip_rotate))
		{
			// First try if there are some fixed prior angles
			// For multi-body refinements, ignore the original priors and get the refined residual angles from the previous iteration
			RFLOAT prior_rot = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT_PRIOR);
			RFLOAT prior_tilt = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT_PRIOR);
			RFLOAT prior_psi = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI_PRIOR);
			RFLOAT prior_psi_flip_ratio = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI_PRIOR_FLIP_RATIO);

			bool do_auto_refine_local_searches = (baseMLO->do_auto_refine) && (baseMLO->sampling.healpix_order >= baseMLO->autosampling_hporder_local_searches);
			bool do_classification_local_searches = (! baseMLO->do_auto_refine) && (baseMLO->mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI)
					&& (baseMLO->mymodel.sigma2_rot > 0.) && (baseMLO->mymodel.sigma2_tilt > 0.) && (baseMLO->mymodel.sigma2_psi > 0.);
			bool do_local_angular_searches = (do_auto_refine_local_searches) || (do_classification_local_searches);

			// If there were no defined priors (i.e. their values were 999.), then use the "normal" angles
			if (prior_rot > 998.99 && prior_rot < 999.01)
				prior_rot = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT);
			if (prior_tilt > 998.99 && prior_tilt < 999.01)
				prior_tilt = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT);
			if (prior_psi > 998.99 && prior_psi < 999.01)
				prior_psi = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI);
			if (prior_psi_flip_ratio > 998.99 && prior_psi_flip_ratio < 999.01)
				prior_psi_flip_ratio = 0.5;

			////////// How does this work now: each particle has a different sampling object?!!!
			// Select only those orientations that have non-zero prior probability

			if (baseMLO->do_helical_refine && baseMLO->mymodel.ref_dim == 3)
			{
				baseMLO->sampling.selectOrientationsWithNonZeroPriorProbabilityFor3DHelicalReconstruction(prior_rot, prior_tilt, prior_psi,
										sqrt(baseMLO->mymodel.sigma2_rot), sqrt(baseMLO->mymodel.sigma2_tilt), sqrt(baseMLO->mymodel.sigma2_psi),
										op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior,
										do_local_angular_searches, prior_psi_flip_ratio);
			}
			else
			{
				baseMLO->sampling.selectOrientationsWithNonZeroPriorProbability(prior_rot, prior_tilt, prior_psi,
						sqrt(baseMLO->mymodel.sigma2_rot), sqrt(baseMLO->mymodel.sigma2_tilt), sqrt(baseMLO->mymodel.sigma2_psi),
						op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);
			}

			long int nr_orients = baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior) * baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior);
			if (nr_orients == 0)
			{
				std::cerr << " sampling.NrDirections()= " << baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior)
						<< " sampling.NrPsiSamplings()= " << baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior) << std::endl;
				REPORT_ERROR("Zero orientations fall within the local angular search. Increase the sigma-value(s) on the orientations!");
			}

		}
		CTOC(accMLO->timer,"nonZeroProb");

		// ------------------------------------------------------------------------------------------

		CTIC(accMLO->timer,"readData");
		// Get the image and recimg data
		if (baseMLO->do_parallel_disc_io)
		{

			// If all followers had preread images into RAM: get those now
			if (baseMLO->do_preread_images)
			{

                img().reshape(baseMLO->mydata.particles[part_id].images[img_id].img);
                CTIC(accMLO->timer,"ParaReadPrereadImages");
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(baseMLO->mydata.particles[part_id].images[img_id].img)
				{
                	DIRECT_MULTIDIM_ELEM(img(), n) = (RFLOAT)DIRECT_MULTIDIM_ELEM(baseMLO->mydata.particles[part_id].images[img_id].img, n);
				}
				CTOC(accMLO->timer,"ParaReadPrereadImages");
			}
			else
			{
				if (accMLO->dataIs3D)
				{
					CTIC(accMLO->timer,"ParaRead3DImages");
					img.read(fn_img);
					img().setXmippOrigin();
					CTOC(accMLO->timer,"ParaRead3DImages");
				}
				else
				{
					CTIC(accMLO->timer,"ParaRead2DImages");
					img() = baseMLO->exp_imgs[my_metadata_offset];
					CTOC(accMLO->timer,"ParaRead2DImages");
				}
			}
			if (baseMLO->has_converged && baseMLO->do_use_reconstruct_images)
			{
				FileName fn_recimg;
				std::istringstream split2(baseMLO->exp_fn_recimg);
				// Get the right line in the exp_fn_img string
				for (int i = 0; i <= my_metadata_offset; i++)
					getline(split2, fn_recimg);
				rec_img.read(fn_recimg);
				rec_img().setXmippOrigin();
			}
		}
		else
		{
			// Unpack the image from the imagedata
			if (accMLO->dataIs3D)
			{
				CTIC(accMLO->timer,"Read3DImages");
				CTIC(accMLO->timer,"resize");
				img().resize(baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group]);
				CTOC(accMLO->timer,"resize");
				// Only allow a single image per call of this function!!! nr_pool needs to be set to 1!!!!
				// This will save memory, as we'll need to store all translated images in memory....
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
				{
					DIRECT_A3D_ELEM(img(), k, i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, k, i, j);
				}
				img().setXmippOrigin();

				if (baseMLO->has_converged && baseMLO->do_use_reconstruct_images)
				{
					rec_img().resize(baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group]);
					int offset = (baseMLO->do_ctf_correction) ? 2 * baseMLO->image_full_size[optics_group] : baseMLO->image_full_size[optics_group];
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(rec_img())
					{
						DIRECT_A3D_ELEM(rec_img(), k, i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, offset + k, i, j);
					}
					rec_img().setXmippOrigin();

				}
				CTOC(accMLO->timer,"Read3DImages");

			}
			else
			{
				CTIC(accMLO->timer,"Read2DImages");
				img().resize(baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group]);
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img())
				{
					DIRECT_A2D_ELEM(img(), i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, my_metadata_offset, i, j);
				}
				img().setXmippOrigin();
				if (baseMLO->has_converged && baseMLO->do_use_reconstruct_images)
				{

					/// TODO: this will be WRONG for multi-image particles, but I guess that's not going to happen anyway...
					int my_nr_particles = baseMLO->exp_my_last_part_id - baseMLO->exp_my_first_part_id + 1;
					rec_img().resize(baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group]);
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(rec_img())
					{
						DIRECT_A2D_ELEM(rec_img(), i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, my_nr_particles + my_metadata_offset, i, j);
					}
					rec_img().setXmippOrigin();
				}
				CTOC(accMLO->timer,"Read2DImages");
			}
		}
		CTOC(accMLO->timer,"readData");

		// ------------------------------------------------------------------------------------------

		size_t current_size_x = baseMLO->image_current_size[optics_group] / 2 + 1;
		size_t current_size_y = baseMLO->image_current_size[optics_group];
		size_t current_size_z = (accMLO->dataIs3D) ? baseMLO->image_current_size[optics_group] : 1;
		accMLO->transformer1.setSize(img().xdim,img().ydim,img().zdim);
		Fimg.initZeros(current_size_z, current_size_y, current_size_x);

		// ------------------------------------------------------------------------------------------

		CTIC(cudaMLO->timer,"makeNoiseMask");
        // Either mask with zeros or noise. Here, make a noise-image that will be optional in the softMask-kernel.
		AccDataTypes::Image<XFLOAT> RandomImage(img(),ptrFactory);

        if (!baseMLO->do_zero_mask) // prepare a acc-side Random image
        {
        		if(RandomImage.is3D())
        				CRITICAL("Noise-masking not supported with acceleration and 3D input: Noise-kernel(s) is hard-coded 2D");

                // Make a F-space image to hold generate and modulate noise
                RandomImage.accAlloc();

                // Set up scalar adjustment factor and random seed
                XFLOAT temp_sigmaFudgeFactor = baseMLO->sigma2_fudge;
                int seed(baseMLO->random_seed + part_id);


    			// Remap mymodel.sigma2_noise[optics_group] onto remapped_sigma2_noise for this images's size and angpix
    			MultidimArray<RFLOAT > remapped_sigma2_noise;
    			remapped_sigma2_noise.initZeros(XSIZE(img())/2+1);
    			RFLOAT remap_image_sizes = (baseMLO->image_full_size[optics_group] * my_pixel_size) / (baseMLO->mymodel.ori_size * baseMLO->mymodel.pixel_size);
    			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(baseMLO->mymodel.sigma2_noise[optics_group])
    			{
    				int i_remap = ROUND(remap_image_sizes * i);
    				if (i_remap < XSIZE(remapped_sigma2_noise))
    					DIRECT_A1D_ELEM(remapped_sigma2_noise, i_remap) = DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[optics_group], i);
    			}


                LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
                // construct the noise-image
                AccUtilities::makeNoiseImage<MlClass>(	temp_sigmaFudgeFactor,
                								remapped_sigma2_noise,
												seed,
												accMLO,
												RandomImage,
												RandomImage.is3D());
                LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
        }
        CTOC(cudaMLO->timer,"makeNoiseMask");

		// ------------------------------------------------------------------------------------------

		CTIC(accMLO->timer,"HelicalPrep");

		/* FIXME :  For some reason the device-allocation inside "selfTranslate" takes a much longer time than expected.
		 * 			I tried moving it up and placing the size under a bunch of if()-cases, but this simply transferred the
		 * 			allocation-cost to that region. /BjoernF,160129
		 */

		// Apply (rounded) old offsets first
		my_old_offset.selfROUND();

		// Helical reconstruction: calculate old_offset in the system of coordinates of the helix, i.e. parallel & perpendicular, depending on psi-angle!
		// For helices do NOT apply old_offset along the direction of the helix!!
		Matrix1D<RFLOAT> my_old_offset_helix_coords;
		RFLOAT rot_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT);
		RFLOAT tilt_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT);
		RFLOAT psi_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI);
		if ( (baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry) )
		{
			// Calculate my_old_offset_helix_coords from my_old_offset and psi angle
			transformCartesianAndHelicalCoords(my_old_offset, my_old_offset_helix_coords, rot_deg, tilt_deg, psi_deg, CART_TO_HELICAL_COORDS);
			// We do NOT want to accumulate the offsets in the direction along the helix (which is X in the helical coordinate system!)
			// However, when doing helical local searches, we accumulate offsets
			// Do NOT accumulate offsets in 3D classification of helices
			if ( (! baseMLO->do_skip_align) && (! baseMLO->do_skip_rotate) )
			{
				// TODO: check whether the following lines make sense
				bool do_auto_refine_local_searches = (baseMLO->do_auto_refine) && (baseMLO->sampling.healpix_order >= baseMLO->autosampling_hporder_local_searches);
				bool do_classification_local_searches = (! baseMLO->do_auto_refine) && (baseMLO->mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI)
						&& (baseMLO->mymodel.sigma2_rot > 0.) && (baseMLO->mymodel.sigma2_tilt > 0.) && (baseMLO->mymodel.sigma2_psi > 0.);
				bool do_local_angular_searches = (do_auto_refine_local_searches) || (do_classification_local_searches);
				if (!do_local_angular_searches)
				{
					if (! accMLO->dataIs3D)
						XX(my_old_offset_helix_coords) = 0.;
					else
						ZZ(my_old_offset_helix_coords) = 0.;
				}
			}
			// TODO: Now re-calculate the my_old_offset in the real (or image) system of coordinate (rotate -psi angle)
			transformCartesianAndHelicalCoords(my_old_offset_helix_coords, my_old_offset, rot_deg, tilt_deg, psi_deg, HELICAL_TO_CART_COORDS);
		}
		CTOC(accMLO->timer,"HelicalPrep");

		// ------------------------------------------------------------------------------------------

		my_old_offset.selfROUND();

		// ------------------------------------------------------------------------------------------

		CTIC(accMLO->timer,"TranslateAndNormCorrect");

		AccDataTypes::Image<XFLOAT> d_img(img.data, ptrFactory);
		AccDataTypes::Image<XFLOAT> d_rec_img(img.data, ptrFactory);

		d_img.allAlloc();
		d_img.allInit(0);

		XFLOAT normcorr_val = baseMLO->do_norm_correction ? (XFLOAT)(baseMLO->mymodel.avg_norm_correction / normcorr) : 1;
		AccUtilities::TranslateAndNormCorrect(	img.data,	// input   	host-side 	MultidimArray
												d_img,		// output  	acc-side  	Array
												normcorr_val,
												XX(my_old_offset),
												YY(my_old_offset),
												(accMLO->dataIs3D) ? ZZ(my_old_offset) : 0.,
												accMLO->dataIs3D);
		LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

		CTOC(accMLO->timer,"TranslateAndNormCorrect");

		// Set up the UNMASKED image to use for reconstruction, which may be a separate image altogether (rec_img)
		//
		//			d_img has the image information which will be masked
		//
		if(baseMLO->has_converged && baseMLO->do_use_reconstruct_images)
		{
			CTIC(accMLO->timer,"TranslateAndNormCorrect_recImg");
			d_rec_img.allAlloc();
			d_rec_img.allInit(0);
			AccUtilities::TranslateAndNormCorrect(	rec_img.data,	// input   	host-side 	MultidimArray
													d_rec_img,		// output  	acc-side  	Array
													normcorr_val,
													XX(my_old_offset),
													YY(my_old_offset),
													(accMLO->dataIs3D) ? ZZ(my_old_offset) : 0.,
													accMLO->dataIs3D);
			LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
			CTOC(accMLO->timer,"TranslateAndNormCorrect_recImg");

			CTIC(cudaMLO->timer,"normalizeAndTransform_recImg");
			// The image used to reconstruct is not masked, so we transform and beam-tilt it
			AccUtilities::normalizeAndTransformImage<MlClass>(d_rec_img,		// input  acc-side  Array
															  Fimg,			// output host-side MultidimArray
															  accMLO,
															  current_size_x,
															  current_size_y,
															  current_size_z);
			LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
			CTOC(cudaMLO->timer,"normalizeAndTransform_recImg");
		}
		else // if we don't have special images, just use the same as for alignment. But do it here, *before masking*
		{
			CTIC(cudaMLO->timer,"normalizeAndTransform_recImg");
			// The image used to reconstruct is not masked, so we transform and beam-tilt it
			AccUtilities::normalizeAndTransformImage<MlClass>(	 d_img,		// input  acc-side  Array
																 Fimg,		// output host-side MultidimArray
																 accMLO,
																 current_size_x,
																 current_size_y,
																 current_size_z);
			LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
			CTOC(cudaMLO->timer,"normalizeAndTransform_recImg");
		}

		// ------------------------------------------------------------------------------------------

		if ( (baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry) )
		{
			// Transform rounded Cartesian offsets to corresponding helical ones
			transformCartesianAndHelicalCoords(my_old_offset, my_old_offset_helix_coords, rot_deg, tilt_deg, psi_deg, CART_TO_HELICAL_COORDS);
			op.old_offset[img_id] = my_old_offset_helix_coords;
		}
		else
		{
			// For multi-bodies: store only the old refined offset, not the constant consensus offset or the projected COM of this body
			if (baseMLO->mymodel.nr_bodies > 1)
				op.old_offset[img_id] = my_refined_ibody_offset;
			else
				op.old_offset[img_id] = my_old_offset;  // Not doing helical refinement. Rounded Cartesian offsets are stored.
		}
		// Also store priors on translations
		op.prior[img_id] = my_prior;

		// ------------------------------------------------------------------------------------------

		CTIC(accMLO->timer,"selfApplyBeamTilt");
		baseMLO->mydata.obsModel.demodulatePhase(optics_group, Fimg);
		baseMLO->mydata.obsModel.divideByMtf(optics_group, Fimg);
		CTOC(accMLO->timer,"selfApplyBeamTilt");

		op.Fimg_nomask.at(img_id) = Fimg;

		// ------------------------------------------------------------------------------------------

		MultidimArray<RFLOAT> Mnoise;
		bool is_helical_segment = (baseMLO->do_helical_refine) || ((baseMLO->mymodel.ref_dim == 2) && (baseMLO->helical_tube_outer_diameter > 0.));

		// For multibodies: have the mask radius equal to maximum radius within body mask plus the translational offset search range
		RFLOAT my_mask_radius = (baseMLO->mymodel.nr_bodies > 1 ) ?
						(baseMLO->mymodel.max_radius_mask_bodies[ibody] + baseMLO->sampling.offset_range) / my_pixel_size :
						baseMLO->particle_diameter / (2. * my_pixel_size);

		// ------------------------------------------------------------------------------------------

		// We are now done with the unmasked image used for reconstruction.
		// Now make the masked image used for alignment and classification.

		if (is_helical_segment)
		{
			CTIC(accMLO->timer,"applyHelicalMask");

			// download img...
			d_img.cpToHost();
			d_img.streamSync();
			d_img.getHost(img());

			// ...modify img...
			if(baseMLO->do_zero_mask)
			{
				softMaskOutsideMapForHelix(img(), psi_deg, tilt_deg, my_mask_radius,
						(baseMLO->helical_tube_outer_diameter / (2. * my_pixel_size)),
						baseMLO->width_mask_edge);
			}
			else
			{
				MultidimArray<RFLOAT> Mnoise;
				RandomImage.hostAlloc();
				RandomImage.cpToHost();
				Mnoise.resize(img());
				RandomImage.getHost(Mnoise);
				softMaskOutsideMapForHelix(img(), psi_deg, tilt_deg, my_mask_radius,
						(baseMLO->helical_tube_outer_diameter / (2. * my_pixel_size)),
						baseMLO->width_mask_edge,
						&Mnoise);
			}

			// ... and re-upload img
			d_img.setHost(img());
			d_img.cpToDevice();
			CTOC(accMLO->timer,"applyHelicalMask");
		}
		else // this is not a helical segment
		{
			CTIC(accMLO->timer,"applyMask");

			// Shared parameters for noise/zero masking
			XFLOAT cosine_width = baseMLO->width_mask_edge;
			XFLOAT radius = (XFLOAT) my_mask_radius;
			if (radius < 0)
				radius = ((RFLOAT)img.data.xdim)/2.;
			XFLOAT radius_p = radius + cosine_width;

			// For zero-masking, we need the background-value
			XFLOAT bg_val(0.);
			if(baseMLO->do_zero_mask)
			{
				AccPtr<XFLOAT> softMaskSum    = ptrFactory.make<XFLOAT>((size_t)SOFTMASK_BLOCK_SIZE, 0);
				AccPtr<XFLOAT> softMaskSum_bg = ptrFactory.make<XFLOAT>((size_t)SOFTMASK_BLOCK_SIZE, 0);
				softMaskSum.accAlloc();
				softMaskSum_bg.accAlloc();
				softMaskSum.accInit(0);
				softMaskSum_bg.accInit(0);

				// Calculate the background value
				AccUtilities::softMaskBackgroundValue(
						d_img,
						radius,
						radius_p,
						cosine_width,
						softMaskSum,
						softMaskSum_bg);

				LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
				softMaskSum.streamSync();

				// Finalize the background value
				bg_val = (RFLOAT) AccUtilities::getSumOnDevice<XFLOAT>(softMaskSum_bg) /
						 (RFLOAT) AccUtilities::getSumOnDevice<XFLOAT>(softMaskSum);
				softMaskSum.streamSync();
			}

			//avoid kernel-calls warning about null-pointer for RandomImage
			if (baseMLO->do_zero_mask)
				RandomImage.setAccPtr(d_img);

			// Apply a cosine-softened mask, using either the background value or the noise-image outside of the radius
			AccUtilities::cosineFilter(
					d_img,
					baseMLO->do_zero_mask,
					RandomImage,
					radius,
					radius_p,
					cosine_width,
					bg_val);

			LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

			CTOC(accMLO->timer,"applyMask");
		}

		// ------------------------------------------------------------------------------------------

		CTIC(cudaMLO->timer,"normalizeAndTransform");
		AccUtilities::normalizeAndTransformImage<MlClass>(	 d_img,		// input
															 Fimg,		// output
															 accMLO,
															 current_size_x,
															 current_size_y,
															 current_size_z);
		LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
		CTOC(cudaMLO->timer,"normalizeAndTransform");

		// ------------------------------------------------------------------------------------------

		CTIC(accMLO->timer,"powerClass");
		// Store the power_class spectrum of the whole image (to fill sigma2_noise between current_size and full_size
		if (baseMLO->image_current_size[optics_group] < baseMLO->image_full_size[optics_group])
		{
			AccPtr<XFLOAT> spectrumAndXi2 = ptrFactory.make<XFLOAT>((size_t)((baseMLO->image_full_size[optics_group]/2+1)+1), 0); // last +1 is the Xi2, to remove an expensive memcpy
			spectrumAndXi2.allAlloc();
			spectrumAndXi2.accInit(0);
			spectrumAndXi2.streamSync();

			int gridSize = CEIL((float)(accMLO->transformer1.fouriers.getSize()) / (float)POWERCLASS_BLOCK_SIZE);
			if(accMLO->dataIs3D)
				AccUtilities::powerClass<true>(gridSize,POWERCLASS_BLOCK_SIZE,
					~accMLO->transformer1.fouriers,
					~spectrumAndXi2,
					accMLO->transformer1.fouriers.getSize(),
					spectrumAndXi2.getSize()-1,
					accMLO->transformer1.xFSize,
					accMLO->transformer1.yFSize,
					accMLO->transformer1.zFSize,
					(baseMLO->image_current_size[optics_group]/2)+1, // note: NOT baseMLO->image_full_size[optics_group]/2+1
					&(~spectrumAndXi2)[spectrumAndXi2.getSize()-1]); // last element is the hihgres_Xi2
			else
				AccUtilities::powerClass<false>(gridSize,POWERCLASS_BLOCK_SIZE,
					~accMLO->transformer1.fouriers,
					~spectrumAndXi2,
					accMLO->transformer1.fouriers.getSize(),
					spectrumAndXi2.getSize()-1,
					accMLO->transformer1.xFSize,
					accMLO->transformer1.yFSize,
					accMLO->transformer1.zFSize,
					(baseMLO->image_current_size[optics_group]/2)+1, // note: NOT baseMLO->image_full_size[optics_group]/2+1
					&(~spectrumAndXi2)[spectrumAndXi2.getSize()-1]); // last element is the hihgres_Xi2

			LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

			spectrumAndXi2.streamSync();
			spectrumAndXi2.cpToHost();
			spectrumAndXi2.streamSync();

			op.power_img.at(img_id).resize(baseMLO->image_full_size[optics_group]/2 + 1);

			for (int i = 0; i<(spectrumAndXi2.getSize()-1); i ++)
				op.power_img.at(img_id).data[i] = spectrumAndXi2[i];
			op.highres_Xi2_img.at(img_id) = spectrumAndXi2[spectrumAndXi2.getSize()-1];
		}
		else
		{
			op.highres_Xi2_img.at(img_id) = 0.;
		}
		CTOC(accMLO->timer,"powerClass");

		Fctf.resize(Fimg);
		// Now calculate the actual CTF
		if (baseMLO->do_ctf_correction)
		{
			if (accMLO->dataIs3D)
			{
				Image<RFLOAT> Ictf;
				if (baseMLO->do_parallel_disc_io)
				{
					CTIC(accMLO->timer,"CTFRead3D_disk");
					// Read CTF-image from disc
					FileName fn_ctf;
					if (!baseMLO->mydata.getImageNameOnScratch(part_id, img_id, fn_ctf, true))
					{
						std::istringstream split(baseMLO->exp_fn_ctf);
						// Get the right line in the exp_fn_img string
						for (int i = 0; i <= my_metadata_offset; i++)
							getline(split, fn_ctf);
					}
					Ictf.read(fn_ctf);
					CTOC(accMLO->timer,"CTFRead3D_disk");
				}
				else
				{
					CTIC(accMLO->timer,"CTFRead3D_array");
					// Unpack the CTF-image from the exp_imagedata array
					Ictf().resize(baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group]);
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Ictf())
					{
						DIRECT_A3D_ELEM(Ictf(), k, i, j) = DIRECT_A3D_ELEM(baseMLO->exp_imagedata, baseMLO->image_full_size[optics_group] + k, i, j);
					}
					CTOC(accMLO->timer,"CTFRead3D_array");
				}
				// Set the CTF-image in Fctf
				CTIC(accMLO->timer,"CTFSet3D_array");
				baseMLO->get3DCTFAndMulti(Ictf(), Fctf, FstMulti, ctf_premultiplied);
				CTOC(accMLO->timer,"CTFSet3D_array");
			}
			else
			{
				CTIC(accMLO->timer,"CTFRead2D");
				CTF ctf;
				ctf.setValuesByGroup(
                                        &(baseMLO->mydata).obsModel, optics_group,
					DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_CTF_DEFOCUS_U),
					DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_CTF_DEFOCUS_V),
					DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_CTF_DEFOCUS_ANGLE),
					DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_CTF_BFACTOR),
					DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_CTF_KFACTOR),
					DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_CTF_PHASE_SHIFT));

				ctf.getFftwImage(Fctf, baseMLO->image_full_size[optics_group], baseMLO->image_full_size[optics_group], my_pixel_size,
						baseMLO->ctf_phase_flipped, baseMLO->only_flip_phases, baseMLO->intact_ctf_first_peak, true, baseMLO->do_ctf_padding);

				// SHWS 13feb2020: when using CTF-premultiplied, from now on use the normal kernels, but replace ctf by ctf^2
				if (ctf_premultiplied)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fctf)
					{
						DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}

				CTOC(accMLO->timer,"CTFRead2D");
			}
		}
		else
		{
			Fctf.initConstant(1.);
		}
		CTOC(accMLO->timer,"ctfCorr");

		CTIC(accMLO->timer,"selfApplyBeamTilt");
		baseMLO->mydata.obsModel.demodulatePhase(optics_group, Fimg);
		baseMLO->mydata.obsModel.divideByMtf(optics_group, Fimg);
		CTOC(accMLO->timer,"selfApplyBeamTilt");

		// Store Fimg and Fctf
		op.Fimg.at(img_id) = Fimg;
		op.Fctf.at(img_id) = Fctf;

		// Correct images and CTFs by Multiplicity, if required, and store it
		if ( NZYXSIZE(FstMulti) > 0 )
		{
			baseMLO->applySubtomoCorrection(op.Fimg.at(img_id), op.Fimg_nomask.at(img_id), op.Fctf.at(img_id), FstMulti);
			op.FstMulti.at(img_id) = FstMulti;
		}

		// If we're doing multibody refinement, now subtract projections of the other bodies from both the masked and the unmasked particle
		if (baseMLO->mymodel.nr_bodies > 1)
		{
			MultidimArray<Complex> Fsum_obody;
			Fsum_obody.initZeros(Fimg);

			for (int obody = 0; obody < baseMLO->mymodel.nr_bodies; obody++)
			{
				if (obody != ibody) // Only subtract if other body is not this body....
				{
					// Get the right metadata
					int ocol_rot  = 0 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_tilt = 1 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_psi  = 2 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_xoff = 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_yoff = 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					int ocol_zoff = 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;
					//int ocol_norm = 6 + METADATA_LINE_LENGTH_BEFORE_BODIES + (obody) * METADATA_NR_BODY_PARAMS;

					Matrix2D<RFLOAT> Aresi,  Abody;
					// Aresi is the residual orientation for this obody
					Euler_angles2matrix(DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, ocol_rot),
										DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, ocol_tilt),
										DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, ocol_psi), Aresi, false);
					// The real orientation to be applied is the obody transformation applied and the original one
					Abody = Aori * (baseMLO->mymodel.orient_bodies[obody]).transpose() * baseMLO->A_rot90 * Aresi * baseMLO->mymodel.orient_bodies[obody];

					// Apply anisotropic mag and scaling
					Abody = baseMLO->mydata.obsModel.applyAnisoMag(Abody, optics_group);
					Abody = baseMLO->mydata.obsModel.applyScaleDifference(Abody, optics_group, baseMLO->mymodel.ori_size, baseMLO->mymodel.pixel_size);

					// Get the FT of the projection in the right direction
					MultidimArray<Complex> FTo;
					FTo.initZeros(Fimg);
					// The following line gets the correct pointer to account for overlap in the bodies
					int oobody = DIRECT_A2D_ELEM(baseMLO->mymodel.pointer_body_overlap, ibody, obody);
					baseMLO->mymodel.PPref[oobody].get2DFourierTransform(FTo, Abody);

					/********************************************************************************
					 * Currently CPU-memory for projectors is not deallocated when doing multibody
					 * due to the previous line. See cpu_ml_optimiser.cpp and cuda_ml_optimiser.cu
					 ********************************************************************************/

					// 17May2017: Body is centered at its own COM
					// move it back to its place in the original particle image
					Matrix1D<RFLOAT> other_projected_com(baseMLO->mymodel.data_dim);

					// Projected COM for this body (using Aori, just like above for ibody and my_projected_com!!!)
					other_projected_com = Aori * (baseMLO->mymodel.com_bodies[obody]);
					// This will have made other_projected_com of size 3 again! resize to mymodel.data_dim
					other_projected_com.resize(baseMLO->mymodel.data_dim);

					// Do the exact same as was done for the ibody, but DONT selfROUND here, as later phaseShift applied to ibody below!!!
					other_projected_com -= my_old_offset_ori;

					// Subtract refined obody-displacement
					XX(other_projected_com) -= DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, ocol_xoff);
					YY(other_projected_com) -= DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, ocol_yoff);
					if (baseMLO->mymodel.data_dim == 3)
						ZZ(other_projected_com) -= DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, ocol_zoff);

					// Add the my_old_offset=selfRound(my_old_offset_ori - my_projected_com) already applied to this image for ibody
					other_projected_com += my_old_offset;

					shiftImageInFourierTransform(FTo, Faux, (RFLOAT)baseMLO->image_full_size[optics_group],
							XX(other_projected_com), YY(other_projected_com), (accMLO->dataIs3D) ? ZZ(other_projected_com) : 0.);

					// Sum the Fourier transforms of all the obodies
					Fsum_obody += Faux;

				} // end if obody != ibody
			} // end for obody

			// Now that we have all the summed projections of the obodies, apply CTF, masks etc
			// Apply the CTF to this reference projection
			if (baseMLO->do_ctf_correction)
			{
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsum_obody)
				{
					DIRECT_MULTIDIM_ELEM(Fsum_obody, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}

				// Also do phase modulation, for beam tilt correction and other asymmetric aberrations
				baseMLO->mydata.obsModel.demodulatePhase(optics_group, Fsum_obody, true); // true means do_modulate_instead
				baseMLO->mydata.obsModel.divideByMtf(optics_group, Fsum_obody, true); // true means do_multiply_instead

			}

			// Subtract the other-body FT from the current image FT
			// First the unmasked one, which will be used for reconstruction
			// Only do this if the flag below is true. Otherwise, use the original particles for reconstruction
			if (baseMLO->do_reconstruct_subtracted_bodies)
				op.Fimg_nomask.at(img_id) -= Fsum_obody;

			// For the masked one, have to mask outside the circular mask to prevent negative values outside the mask in the subtracted image!
			CenterFFTbySign(Fsum_obody);
			windowFourierTransform(Fsum_obody, Faux, baseMLO->image_full_size[optics_group]);
			accMLO->transformer.inverseFourierTransform(Faux, img());

			softMaskOutsideMap(img(), my_mask_radius, (RFLOAT)baseMLO->width_mask_edge);

			// And back to Fourier space now
			accMLO->transformer.FourierTransform(img(), Faux);
			windowFourierTransform(Faux, Fsum_obody, baseMLO->image_current_size[optics_group]);
			CenterFFTbySign(Fsum_obody);

			// Subtract the other-body FT from the masked exp_Fimgs
			op.Fimg.at(img_id) -= Fsum_obody;

			// 23jul17: NEW: as we haven't applied the (nonROUNDED!!)  my_refined_ibody_offset yet, do this now in the FourierTransform
			Faux = op.Fimg.at(img_id);
			shiftImageInFourierTransform(Faux, op.Fimg.at(img_id), (RFLOAT)baseMLO->image_full_size[optics_group],
					XX(my_refined_ibody_offset), YY(my_refined_ibody_offset), (accMLO->dataIs3D) ? ZZ(my_refined_ibody_offset) : 0);
			Faux = op.Fimg_nomask.at(img_id);
			shiftImageInFourierTransform(Faux, op.Fimg_nomask.at(img_id), (RFLOAT)baseMLO->image_full_size[optics_group],
					XX(my_refined_ibody_offset), YY(my_refined_ibody_offset), (accMLO->dataIs3D) ? ZZ(my_refined_ibody_offset) : 0);
		} // end if mymodel.nr_bodies > 1


	} // end loop img_id
	//accMLO->transformer.clear();
#ifdef TIMING
	if (part_id == baseMLO->exp_my_first_part_id)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_FT);
#endif
	GTOC(accMLO->timer,"getFourierTransformsAndCtfs");
	GATHERGPUTIMINGS(accMLO->timer);
}

// ----------------------------------------------------------------------------
// ------------------ getAllSquaredDifferencesCoarse --------------------------
// ----------------------------------------------------------------------------
template <class MlClass>
void getAllSquaredDifferencesCoarse(
		unsigned exp_ipass,
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		MlClass *accMLO,
	 	AccPtr<XFLOAT> &Mweight,
	 	AccPtrFactory ptrFactory,
	 	int ibody = 0)
{

#ifdef TIMING
	if (op.part_id == baseMLO->exp_my_first_part_id)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF1);
#endif

	CUSTOM_ALLOCATOR_REGION_NAME("DIFF_COARSE");

	CTIC(accMLO->timer,"diff_pre_gpu");
	unsigned long weightsPerPart(baseMLO->mymodel.nr_classes * sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.nr_oversampled_rot * sp.nr_oversampled_trans);

	std::vector<MultidimArray<Complex > > dummy;
	std::vector<std::vector<MultidimArray<Complex > > > dummy2;
	std::vector<MultidimArray<RFLOAT> > dummyRF;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, false, op.part_id, sp.current_oversampling, op.metadata_offset, // inserted SHWS 12112015
			sp.itrans_min, sp.itrans_max, op.Fimg, dummy, op.Fctf, dummy2, dummy2,
			op.local_Fctf, op.local_sqrtXi2, op.local_Minvsigma2, op.FstMulti, dummyRF);

	CTOC(accMLO->timer,"diff_pre_gpu");

	std::vector< AccProjectorPlan > projectorPlans(0, (CudaCustomAllocator *)accMLO->getAllocator());

	//If particle specific sampling plan required
	if (accMLO->generateProjectionPlanOnTheFly)
	{
		CTIC(accMLO->timer,"generateProjectionSetupCoarse");

		projectorPlans.resize(baseMLO->mymodel.nr_classes, (CudaCustomAllocator *)accMLO->getAllocator());


		for (unsigned long iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
		{
			if (baseMLO->mymodel.pdf_class[iclass] > 0.)
			{
				Matrix2D<RFLOAT> MBL, MBR, Aori;

				if (baseMLO->mymodel.nr_bodies > 1)
				{
					// img_id=0 because in multi-body refinement we do not do movie frames!
					RFLOAT rot_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_ROT);
					RFLOAT tilt_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_TILT);
					RFLOAT psi_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_PSI);
					Euler_angles2matrix(rot_ori, tilt_ori, psi_ori, Aori, false);

					MBL = Aori * (baseMLO->mymodel.orient_bodies[ibody]).transpose() * baseMLO->A_rot90;
					MBR = baseMLO->mymodel.orient_bodies[ibody];
				}

				int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, 0); // get optics group of first image for this particle...
				Matrix2D<RFLOAT> mag;
				mag.initIdentity(3);
				mag = baseMLO->mydata.obsModel.applyAnisoMag(mag, optics_group);
				mag = baseMLO->mydata.obsModel.applyScaleDifference(mag, optics_group, baseMLO->mymodel.ori_size, baseMLO->mymodel.pixel_size);
				if (!mag.isIdentity())
				{
					if (MBL.mdimx == 3 && MBL.mdimx ==3) MBL = mag * MBL;
					else MBL = mag;
				}

				projectorPlans[iclass].setup(
						baseMLO->sampling,
						op.directions_prior,
						op.psi_prior,
						op.pointer_dir_nonzeroprior,
						op.pointer_psi_nonzeroprior,
						NULL, //Mcoarse_significant
						baseMLO->mymodel.pdf_class,
						baseMLO->mymodel.pdf_direction,
						sp.nr_dir,
						sp.nr_psi,
						sp.idir_min,
						sp.idir_max,
						sp.ipsi_min,
						sp.ipsi_max,
						sp.itrans_min,
						sp.itrans_max,
						0, //current_oversampling
						1, //nr_oversampled_rot
						iclass,
						true, //coarse
						!IS_NOT_INV,
						baseMLO->do_skip_align,
						baseMLO->do_skip_rotate,
						baseMLO->mymodel.orientational_prior_mode,
						MBL,
						MBR
						);
			}
		}
		CTOC(accMLO->timer,"generateProjectionSetupCoarse");
	}
	else
		projectorPlans = accMLO->bundle->coarseProjectionPlans;

	// Loop only from sp.iclass_min to sp.iclass_max to deal with seed generation in first iteration
	size_t allWeights_size(0);
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		allWeights_size += projectorPlans[exp_iclass].orientation_num * sp.nr_trans*sp.nr_oversampled_trans;

	AccPtr<XFLOAT> allWeights = ptrFactory.make<XFLOAT>(allWeights_size);

	allWeights.accAlloc();
	deviceInitValue<XFLOAT>(allWeights, 0);  // Make sure entire array initialized

	long int allWeights_pos=0;	bool do_CC = (baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc;

	for (int img_id = 0; img_id < sp.nr_images; img_id++)
	{
		int my_metadata_offset = op.metadata_offset + img_id;
		long int group_id = baseMLO->mydata.getGroupId(op.part_id, img_id);
		RFLOAT my_pixel_size = baseMLO->mydata.getImagePixelSize(op.part_id, img_id);
		int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);
		unsigned long image_size = op.local_Minvsigma2[img_id].nzyxdim;

		/*====================================
				Generate Translations
		======================================*/

		CTIC(accMLO->timer,"translation_1");

		long unsigned translation_num((sp.itrans_max - sp.itrans_min + 1) * sp.nr_oversampled_trans);
		// here we introduce offsets for the trans_ and img_ in an array as it is more efficient to
		// copy one big array to/from GPU rather than four small arrays
		size_t trans_x_offset = 0*(size_t)translation_num;
		size_t trans_y_offset = 1*(size_t)translation_num;
		size_t trans_z_offset = 2*(size_t)translation_num;
		size_t img_re_offset = 0*(size_t)image_size;
		size_t img_im_offset = 1*(size_t)image_size;

		AccPtr<XFLOAT> Fimg_ = ptrFactory.make<XFLOAT>((size_t)image_size*2);
		AccPtr<XFLOAT> trans_xyz = ptrFactory.make<XFLOAT>((size_t)translation_num*3);

		Fimg_.allAlloc();
		trans_xyz.allAlloc();

		std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;

		for (long int itrans = 0; itrans < translation_num; itrans++)
		{
			baseMLO->sampling.getTranslationsInPixel(itrans, 0, my_pixel_size, oversampled_translations_x,
					oversampled_translations_y, oversampled_translations_z,
					(baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry));

			RFLOAT xshift = 0., yshift = 0., zshift = 0.;

			xshift = oversampled_translations_x[0];
			yshift = oversampled_translations_y[0];
			if (accMLO->dataIs3D)
				zshift = oversampled_translations_z[0];

			if ( (baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry) )
			{
				RFLOAT rot_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT);
				RFLOAT tilt_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT);
				RFLOAT psi_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata,my_metadata_offset, METADATA_PSI);
				transformCartesianAndHelicalCoords(xshift, yshift, zshift, xshift, yshift, zshift, rot_deg, tilt_deg, psi_deg, (accMLO->dataIs3D) ? (3) : (2), HELICAL_TO_CART_COORDS);
			}

			trans_xyz[trans_x_offset+itrans] = -2 * PI * xshift / (double)baseMLO->image_full_size[optics_group];
			trans_xyz[trans_y_offset+itrans] = -2 * PI * yshift / (double)baseMLO->image_full_size[optics_group];
			trans_xyz[trans_z_offset+itrans] = -2 * PI * zshift / (double)baseMLO->image_full_size[optics_group];
		}

		XFLOAT scale_correction = baseMLO->do_scale_correction ? baseMLO->mymodel.scale_correction[group_id] : 1;

		int exp_current_image_size = (baseMLO->strict_highres_exp > 0.|| baseMLO->adaptive_oversampling > 0) ?
				baseMLO->image_coarse_size[optics_group] : baseMLO->image_current_size[optics_group];
		MultidimArray<Complex > Fimg;
		windowFourierTransform(op.Fimg[img_id], Fimg, exp_current_image_size);

		for (unsigned long i = 0; i < image_size; i ++)
		{
			XFLOAT pixel_correction = 1.0/scale_correction;
			if (baseMLO->do_ctf_correction && fabs(op.local_Fctf[img_id].data[i]) > 1e-8)
			{
				// if ctf[i]==0, pix_corr[i] becomes NaN.
				// However, corr_img[i]==0, so pix-diff in kernel==0.
				// This is ok since originally, pix-diff==Img.real^2 + Img.imag^2,
				// which is ori-indep, and we subtract min_diff form ALL orients.

				if (baseMLO->refs_are_ctf_corrected)
				{
					pixel_correction /= op.local_Fctf[img_id].data[i];
				}
			}
			Fimg_[img_re_offset+i] = Fimg.data[i].real * pixel_correction;
			Fimg_[img_im_offset+i] = Fimg.data[i].imag * pixel_correction;
		}

		trans_xyz.cpToDevice();

		Fimg_.cpToDevice();

		CTOC(accMLO->timer,"translation_1");

		// To speed up calculation, several image-corrections are grouped into a single pixel-wise "filter", or image-correciton

		AccPtr<XFLOAT> corr_img = ptrFactory.make<XFLOAT>((size_t)image_size);

		corr_img.allAlloc();

		buildCorrImage(baseMLO,op,corr_img,img_id,group_id);
		corr_img.cpToDevice();

		deviceInitValue<XFLOAT>(allWeights, (XFLOAT) (op.highres_Xi2_img[img_id] / 2.));
		allWeights_pos = 0;

		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[exp_iclass]));
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

		for (unsigned long iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
		{
			int iproj;
			if (baseMLO->mymodel.nr_bodies > 1) iproj = ibody;
			else                                iproj = iclass;

			if ( projectorPlans[iclass].orientation_num > 0 )
			{
				AccProjectorKernel projKernel = AccProjectorKernel::makeKernel(
						accMLO->bundle->projectors[iproj],
						op.local_Minvsigma2[img_id].xdim,
						op.local_Minvsigma2[img_id].ydim,
						op.local_Minvsigma2[img_id].zdim,
						op.local_Minvsigma2[img_id].xdim-1);

				runDiff2KernelCoarse(
						projKernel,
						&(~trans_xyz)[trans_x_offset], //~trans_x,
						&(~trans_xyz)[trans_y_offset], //~trans_y,
						&(~trans_xyz)[trans_z_offset], //~trans_z,
						~corr_img,
						&(~Fimg_)[img_re_offset], //~Fimg_real,
						&(~Fimg_)[img_im_offset], //~Fimg_imag,
						~projectorPlans[iclass].eulers,
						&(~allWeights)[allWeights_pos],
						(XFLOAT) op.local_sqrtXi2[img_id],
						projectorPlans[iclass].orientation_num,
						translation_num,
						image_size,
						accMLO->classStreams[iclass],
						do_CC,
						accMLO->dataIs3D);

				mapAllWeightsToMweights(
						~projectorPlans[iclass].iorientclasses,
						&(~allWeights)[allWeights_pos],
						&(~Mweight)[img_id*weightsPerPart],
						projectorPlans[iclass].orientation_num,
						translation_num,
						accMLO->classStreams[iclass]
						);

				/*====================================
				    	   Retrieve Results
				======================================*/
				allWeights_pos += projectorPlans[iclass].orientation_num*translation_num;

			}
		}

		for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[exp_iclass]));
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread)); // does not appear to be NEEDED FOR NON-BLOCKING CLASS STREAMS in tests, but should be to sync against classStreams

		op.min_diff2[img_id] = AccUtilities::getMinOnDevice<XFLOAT>(allWeights);

	} // end loop img_id

#ifdef TIMING
	if (op.part_id == baseMLO->exp_my_first_part_id)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF1);
#endif
}

// ----------------------------------------------------------------------------
// -------------------- getAllSquaredDifferencesFine --------------------------
// ----------------------------------------------------------------------------
template <class MlClass>
void getAllSquaredDifferencesFine(
		unsigned exp_ipass,
	OptimisationParamters &op,
	SamplingParameters &sp,
	MlOptimiser *baseMLO,
	MlClass *accMLO,
	std::vector<IndexedDataArray > &FinePassWeights,
	std::vector<std::vector< IndexedDataArrayMask > > &FPCMasks,
	std::vector<ProjectionParams> &FineProjectionData,
	AccPtrFactory ptrFactory,
	int ibody,
	std::vector<AccPtrBundle > &bundleD2)
{
#ifdef TIMING
	if (op.part_id == baseMLO->exp_my_first_part_id)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2);
#endif

	CUSTOM_ALLOCATOR_REGION_NAME("DIFF_FINE");
	CTIC(accMLO->timer,"diff_pre_gpu");

	CTIC(accMLO->timer,"precalculateShiftedImagesCtfsAndInvSigma2s");
	std::vector<MultidimArray<Complex > > dummy;
	std::vector<std::vector<MultidimArray<Complex > > > dummy2;
	std::vector<MultidimArray<RFLOAT> > dummyRF;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, false, op.part_id, sp.current_oversampling, op.metadata_offset, // inserted SHWS 12112015
			sp.itrans_min, sp.itrans_max, op.Fimg, dummy, op.Fctf, dummy2, dummy2,
			op.local_Fctf, op.local_sqrtXi2, op.local_Minvsigma2, op.FstMulti, dummyRF);
	CTOC(accMLO->timer,"precalculateShiftedImagesCtfsAndInvSigma2s");


	CTOC(accMLO->timer,"diff_pre_gpu");

	/*=======================================================================================
										  Particle Iteration
	=========================================================================================*/
	for (int img_id = 0; img_id < sp.nr_images; img_id++)
	{
		// Reset size without de-allocating: we will append everything significant within
		// the current allocation and then re-allocate the then determined (smaller) volume

		int my_metadata_offset = op.metadata_offset + img_id;
		long int group_id = baseMLO->mydata.getGroupId(op.part_id, img_id);
		RFLOAT my_pixel_size = baseMLO->mydata.getImagePixelSize(op.part_id, img_id);
		int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);
		unsigned long image_size = op.local_Minvsigma2[img_id].nzyxdim;

		MultidimArray<Complex > Fref;
		Fref.resize(op.local_Minvsigma2[img_id]);

		/*====================================
				Generate Translations
		======================================*/

		CTIC(accMLO->timer,"translation_2");

		long unsigned translation_num((sp.itrans_max - sp.itrans_min + 1) * sp.nr_oversampled_trans);
		// here we introduce offsets for the trans_ and img_ in an array as it is more efficient to
		// copy one big array to/from GPU rather than four small arrays
		size_t trans_x_offset = 0*(size_t)translation_num;
		size_t trans_y_offset = 1*(size_t)translation_num;
		size_t trans_z_offset = 2*(size_t)translation_num;
		size_t img_re_offset = 0*(size_t)image_size;
		size_t img_im_offset = 1*(size_t)image_size;

		AccPtr<XFLOAT> Fimg_     = ptrFactory.make<XFLOAT>((size_t)image_size*2);
		AccPtr<XFLOAT> trans_xyz = ptrFactory.make<XFLOAT>((size_t)translation_num*3);

		Fimg_.allAlloc();
		trans_xyz.allAlloc();

		std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;

		int j = 0;
		for (long int itrans = 0; itrans < (sp.itrans_max - sp.itrans_min + 1); itrans++)
		{
			baseMLO->sampling.getTranslationsInPixel(itrans, baseMLO->adaptive_oversampling, my_pixel_size, oversampled_translations_x,
					oversampled_translations_y, oversampled_translations_z,
					(baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry));

			for (long int iover_trans = 0; iover_trans < oversampled_translations_x.size(); iover_trans++)
			{
				RFLOAT xshift = 0., yshift = 0., zshift = 0.;

				xshift = oversampled_translations_x[iover_trans];
				yshift = oversampled_translations_y[iover_trans];
				if (accMLO->dataIs3D)
					zshift = oversampled_translations_z[iover_trans];

				if ( (baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry) )
				{
					RFLOAT rot_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT);
					RFLOAT tilt_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT);
					RFLOAT psi_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI);
					transformCartesianAndHelicalCoords(xshift, yshift, zshift, xshift, yshift, zshift, rot_deg, tilt_deg, psi_deg, (accMLO->dataIs3D) ? (3) : (2), HELICAL_TO_CART_COORDS);
				}

				trans_xyz[trans_x_offset+j] = -2 * PI * xshift / (double)baseMLO->image_full_size[optics_group];
				trans_xyz[trans_y_offset+j] = -2 * PI * yshift / (double)baseMLO->image_full_size[optics_group];
				trans_xyz[trans_z_offset+j] = -2 * PI * zshift / (double)baseMLO->image_full_size[optics_group];
				j ++;
			}
		}

		XFLOAT scale_correction = baseMLO->do_scale_correction ? baseMLO->mymodel.scale_correction[group_id] : 1;

		int exp_current_image_size = (baseMLO->strict_highres_exp > 0.) ? baseMLO->image_coarse_size[optics_group] : baseMLO->image_current_size[optics_group];
		MultidimArray<Complex > Fimg, Fimg_nomask;
		windowFourierTransform(op.Fimg[img_id], Fimg, exp_current_image_size);

		for (unsigned long i = 0; i < image_size; i ++)
		{
			XFLOAT pixel_correction = 1.0/scale_correction;
			if (baseMLO->do_ctf_correction && fabs(op.local_Fctf[img_id].data[i]) > 1e-8)
			{
				// if ctf[i]==0, pix_corr[i] becomes NaN.
				// However, corr_img[i]==0, so pix-diff in kernel==0.
				// This is ok since originally, pix-diff==Img.real^2 + Img.imag^2,
				// which is ori-indep, and we subtract min_diff form ALL orients.

				if (baseMLO->refs_are_ctf_corrected)
				{
					pixel_correction /= op.local_Fctf[img_id].data[i];
				}
			}

			Fimg_[img_re_offset+i] = Fimg.data[i].real * pixel_correction;
			Fimg_[img_im_offset+i] = Fimg.data[i].imag * pixel_correction;
		}

		CTOC(accMLO->timer,"translation_2");


		CTIC(accMLO->timer,"kernel_init_1");

		AccPtr<XFLOAT> corr_img = ptrFactory.make<XFLOAT>((size_t)image_size);

		corr_img.allAlloc();
		buildCorrImage(baseMLO,op,corr_img,img_id,group_id);

		trans_xyz.cpToDevice();


		Fimg_.cpToDevice();
		corr_img.cpToDevice();

		CTOC(accMLO->timer,"kernel_init_1");

		std::vector< AccPtr<XFLOAT> > eulers((size_t)(sp.iclass_max-sp.iclass_min+1), ptrFactory.make<XFLOAT>());

		AccPtrBundle AllEulers = ptrFactory.makeBundle();
		AllEulers.setSize(9*FineProjectionData[img_id].orientationNumAllClasses*sizeof(XFLOAT));
		AllEulers.allAlloc();

		unsigned long newDataSize(0);

		for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			FPCMasks[img_id][exp_iclass].weightNum=0;

			if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) && (FineProjectionData[img_id].class_entries[exp_iclass] > 0) )
			{
				// use "slice" constructor with class-specific parameters to retrieve a temporary ProjectionParams with data for this class
				ProjectionParams thisClassProjectionData(	FineProjectionData[img_id],
															FineProjectionData[img_id].class_idx[exp_iclass],
															FineProjectionData[img_id].class_idx[exp_iclass]+FineProjectionData[img_id].class_entries[exp_iclass]);
				// since we retrieved the ProjectionParams for *the whole* class the orientation_num is also equal.

				thisClassProjectionData.orientation_num[0] = FineProjectionData[img_id].class_entries[exp_iclass];
				long unsigned orientation_num  = thisClassProjectionData.orientation_num[0];

				if(orientation_num==0)
					continue;

				CTIC(accMLO->timer,"pair_list_1");
				long unsigned significant_num(0);
				long int nr_over_orient = baseMLO->sampling.oversamplingFactorOrientations(sp.current_oversampling);
				long int nr_over_trans = baseMLO->sampling.oversamplingFactorTranslations(sp.current_oversampling);
				// Prepare the mask of the weight-array for this class
				if (FPCMasks[img_id][exp_iclass].weightNum==0)
					FPCMasks[img_id][exp_iclass].firstPos = newDataSize;

				long unsigned ihidden(0);
				std::vector< long unsigned > iover_transes, ihiddens;

				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++, ihidden++)
				{
					for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++)
					{
						ihiddens.push_back(ihidden);
						iover_transes.push_back(iover_trans);
					}
				}

				int chunkSize(0);
				if(accMLO->dataIs3D)
					chunkSize = D2F_CHUNK_DATA3D;
				else if(accMLO->refIs3D)
					chunkSize = D2F_CHUNK_DATA3D;
				else
					chunkSize = D2F_CHUNK_2D;

				// Do more significance checks on translations and create jobDivision
				significant_num = makeJobsForDiff2Fine(	op,	sp,												// alot of different type inputs...
														orientation_num, translation_num,
														thisClassProjectionData,
														iover_transes, ihiddens,
														nr_over_orient, nr_over_trans, img_id,
														FinePassWeights[img_id],
														FPCMasks[img_id][exp_iclass],   // ..and output into index-arrays mask...
														chunkSize);                    // ..based on a given maximum chunk-size

				// extend size by number of significants found this class
				newDataSize += significant_num;
				FPCMasks[img_id][exp_iclass].weightNum = significant_num;
				FPCMasks[img_id][exp_iclass].lastPos = FPCMasks[img_id][exp_iclass].firstPos + significant_num;
				CTOC(accMLO->timer,"pair_list_1");

				CTIC(accMLO->timer,"IndexedArrayMemCp2");
				bundleD2[img_id].pack(FPCMasks[img_id][exp_iclass].jobOrigin);
				bundleD2[img_id].pack(FPCMasks[img_id][exp_iclass].jobExtent);
				CTOC(accMLO->timer,"IndexedArrayMemCp2");

				Matrix2D<RFLOAT> MBL, MBR;

				if (baseMLO->mymodel.nr_bodies > 1)
				{
					Matrix2D<RFLOAT> Aori;
					RFLOAT rot_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_ROT);
					RFLOAT tilt_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_TILT);
					RFLOAT psi_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_PSI);
					Euler_angles2matrix(rot_ori, tilt_ori, psi_ori, Aori, false);

					MBL = Aori * (baseMLO->mymodel.orient_bodies[ibody]).transpose() * baseMLO->A_rot90;
					MBR = baseMLO->mymodel.orient_bodies[ibody];
				}

				CTIC(accMLO->timer,"generateEulerMatrices");
				eulers[exp_iclass-sp.iclass_min].setSize(9*FineProjectionData[img_id].class_entries[exp_iclass]);
				eulers[exp_iclass-sp.iclass_min].hostAlloc();

				Matrix2D<RFLOAT> mag;
				mag.initIdentity(3);
				mag = baseMLO->mydata.obsModel.applyAnisoMag(mag, optics_group);
				mag = baseMLO->mydata.obsModel.applyScaleDifference(mag, optics_group, baseMLO->mymodel.ori_size, baseMLO->mymodel.pixel_size);
				if (!mag.isIdentity())
				{
					if (MBL.mdimx == 3 && MBL.mdimx ==3) MBL = mag * MBL;
					else MBL = mag;
				}

				generateEulerMatrices(
						thisClassProjectionData,
						&(eulers[exp_iclass-sp.iclass_min])[0],
						true,
						MBL,
						MBR);

				AllEulers.pack(eulers[exp_iclass-sp.iclass_min]);

				CTOC(accMLO->timer,"generateEulerMatrices");
			}
		}

		bundleD2[img_id].cpToDevice();
		AllEulers.cpToDevice();

		FinePassWeights[img_id].rot_id.cpToDevice(); //FIXME this is not used
		FinePassWeights[img_id].rot_idx.cpToDevice();
		FinePassWeights[img_id].trans_idx.cpToDevice();

		for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[exp_iclass]));
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

		for (unsigned long iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
		{
			int iproj;
			if (baseMLO->mymodel.nr_bodies > 1) iproj = ibody;
			else                                iproj = iclass;

			if ((baseMLO->mymodel.pdf_class[iclass] > 0.) && (FineProjectionData[img_id].class_entries[iclass] > 0) )
			{
				long unsigned orientation_num  = FineProjectionData[img_id].class_entries[iclass];
				if(orientation_num==0)
					continue;

				long unsigned significant_num(FPCMasks[img_id][iclass].weightNum);
				if(significant_num==0)
					continue;

				CTIC(accMLO->timer,"Diff2MakeKernel");
				AccProjectorKernel projKernel = AccProjectorKernel::makeKernel(
						accMLO->bundle->projectors[iproj],
						op.local_Minvsigma2[img_id].xdim,
						op.local_Minvsigma2[img_id].ydim,
						op.local_Minvsigma2[img_id].zdim,
						op.local_Minvsigma2[img_id].xdim-1);
				CTOC(accMLO->timer,"Diff2MakeKernel");

				// Use the constructed mask to construct a partial class-specific input
				IndexedDataArray thisClassFinePassWeights(FinePassWeights[img_id],FPCMasks[img_id][iclass]);

				CTIC(accMLO->timer,"Diff2CALL");

				runDiff2KernelFine(
						projKernel,
						~corr_img,
						&(~Fimg_)[img_re_offset], //~Fimg_real,
						&(~Fimg_)[img_im_offset], //~Fimg_imag,
						&(~trans_xyz)[trans_x_offset], //~trans_x,
						&(~trans_xyz)[trans_y_offset], //~trans_y,
						&(~trans_xyz)[trans_z_offset], //~trans_z,
						~eulers[iclass-sp.iclass_min],
						~thisClassFinePassWeights.rot_id,
						~thisClassFinePassWeights.rot_idx,
						~thisClassFinePassWeights.trans_idx,
						~FPCMasks[img_id][iclass].jobOrigin,
						~FPCMasks[img_id][iclass].jobExtent,
						~thisClassFinePassWeights.weights,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						significant_num,
						image_size,
						img_id,
						iclass,
						accMLO->classStreams[iclass],
						FPCMasks[img_id][iclass].jobOrigin.getSize(),
						((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc),
						accMLO->dataIs3D
						);

//				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));
				CTOC(accMLO->timer,"Diff2CALL");

			} // end if class significant
		} // end loop iclass

		for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[exp_iclass]));
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

		FinePassWeights[img_id].setDataSize( newDataSize );

		CTIC(accMLO->timer,"collect_data_1");
		if(baseMLO->adaptive_oversampling!=0)
		{
			op.min_diff2[img_id] = (RFLOAT) AccUtilities::getMinOnDevice<XFLOAT>(FinePassWeights[img_id].weights);
		}
		CTOC(accMLO->timer,"collect_data_1");
//		std::cerr << "  fine pass minweight  =  " << op.min_diff2[img_id] << std::endl;

	}// end loop img_id
#ifdef TIMING
	if (op.part_id == baseMLO->exp_my_first_part_id)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2);
#endif
}

// ----------------------------------------------------------------------------
// -------------- convertAllSquaredDifferencesToWeights -----------------------
// ----------------------------------------------------------------------------
template<class MlClass>
void convertAllSquaredDifferencesToWeights(unsigned exp_ipass,
											OptimisationParamters &op,
											SamplingParameters &sp,
											MlOptimiser *baseMLO,
											MlClass *accMLO,
											std::vector< IndexedDataArray > &PassWeights,
											std::vector< std::vector< IndexedDataArrayMask > > &FPCMasks,
											AccPtr<XFLOAT> &Mweight, // FPCMasks = Fine-Pass Class-Masks
											AccPtrFactory ptrFactory,
											int ibody)
{
#ifdef TIMING
	if (op.part_id == baseMLO->exp_my_first_part_id)
	{
		if (exp_ipass == 0) baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT1);
		else baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT2);
	}
#endif

	RFLOAT my_sigma2_offset = (baseMLO->mymodel.nr_bodies > 1) ?
			baseMLO->mymodel.sigma_offset_bodies[ibody]*baseMLO->mymodel.sigma_offset_bodies[ibody] : baseMLO->mymodel.sigma2_offset;

	// Ready the "prior-containers" for all classes (remake every img_id)
	AccPtr<XFLOAT>  pdf_orientation       = ptrFactory.make<XFLOAT>((size_t)((sp.iclass_max-sp.iclass_min+1) * sp.nr_dir * sp.nr_psi));
	AccPtr<bool>    pdf_orientation_zeros = ptrFactory.make<bool>(pdf_orientation.getSize());
	AccPtr<XFLOAT>  pdf_offset            = ptrFactory.make<XFLOAT>((size_t)((sp.iclass_max-sp.iclass_min+1)*sp.nr_trans));
	AccPtr<bool>    pdf_offset_zeros      = ptrFactory.make<bool>(pdf_offset.getSize());

	pdf_orientation.accAlloc();
	pdf_orientation_zeros.accAlloc();
	pdf_offset.allAlloc();
	pdf_offset_zeros.allAlloc();

	CUSTOM_ALLOCATOR_REGION_NAME("CASDTW_PDF");

	// pdf_orientation is img_id-independent, so we keep it above img_id scope
	CTIC(accMLO->timer,"get_orient_priors");
	AccPtr<RFLOAT> pdfs				= ptrFactory.make<RFLOAT>((size_t)((sp.iclass_max-sp.iclass_min+1) * sp.nr_dir * sp.nr_psi));
	pdfs.allAlloc();

	for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		for (unsigned long idir = sp.idir_min, iorientclass = (exp_iclass-sp.iclass_min) * sp.nr_dir * sp.nr_psi; idir <=sp.idir_max; idir++)
			for (unsigned long ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorientclass++)
			{
				RFLOAT pdf(0);

				if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
					pdf = baseMLO->mymodel.pdf_class[exp_iclass];
				else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
					pdf = DIRECT_MULTIDIM_ELEM(baseMLO->mymodel.pdf_direction[exp_iclass], idir);
				else
					pdf = op.directions_prior[idir] * op.psi_prior[ipsi];

				pdfs[iorientclass] = pdf;
			}

	pdfs.cpToDevice();
	AccUtilities::initOrientations(pdfs, pdf_orientation, pdf_orientation_zeros);
	CTOC(accMLO->timer,"get_orient_priors");

	if(exp_ipass==0 || baseMLO->adaptive_oversampling!=0)
	{
		op.sum_weight.clear();
		op.sum_weight.resize(sp.nr_images, (RFLOAT)(sp.nr_images));
		op.max_weight.clear();
		op.max_weight.resize(sp.nr_images, (RFLOAT)-1);
	}

	if (exp_ipass==0)
		op.Mcoarse_significant.resizeNoCp(1,1,sp.nr_images, XSIZE(op.Mweight));

	XFLOAT my_significant_weight;
	op.significant_weight.clear();
	op.significant_weight.resize(sp.nr_images, 0.);

	// loop over all images inside this particle
	for (int img_id = 0; img_id < sp.nr_images; img_id++)
	{
		int my_metadata_offset = op.metadata_offset + img_id;
		RFLOAT my_pixel_size = baseMLO->mydata.getImagePixelSize(op.part_id, img_id);

		RFLOAT old_offset_x, old_offset_y, old_offset_z;

		if (baseMLO->mymodel.nr_bodies > 1)
		{
			old_offset_x = old_offset_y = old_offset_z = 0.;
		}
		else
		{
			old_offset_x = XX(op.old_offset[img_id]);
			old_offset_y = YY(op.old_offset[img_id]);
			if (accMLO->dataIs3D)
				old_offset_z = ZZ(op.old_offset[img_id]);
		}

		if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
		{
			if(exp_ipass==0)
			{
				int nr_coarse_weights = (sp.iclass_max-sp.iclass_min+1)*sp.nr_images * sp.nr_dir * sp.nr_psi * sp.nr_trans;
				PassWeights[img_id].weights.setAccPtr(&(~Mweight)[img_id*nr_coarse_weights]);
				PassWeights[img_id].weights.setHostPtr(&Mweight[img_id*nr_coarse_weights]);
				PassWeights[img_id].weights.setSize(nr_coarse_weights);
			}
			PassWeights[img_id].weights.doFreeHost=false;

			std::pair<size_t, XFLOAT> min_pair=AccUtilities::getArgMinOnDevice<XFLOAT>(PassWeights[img_id].weights);
			PassWeights[img_id].weights.cpToHost();
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

			//Set all device-located weights to zero, and only the smallest one to 1.
#ifdef _CUDA_ENABLED
			DEBUG_HANDLE_ERROR(cudaMemsetAsync(~(PassWeights[img_id].weights), 0.f, PassWeights[img_id].weights.getSize()*sizeof(XFLOAT),0));

			XFLOAT unity=1;
			DEBUG_HANDLE_ERROR(cudaMemcpyAsync( &(PassWeights[img_id].weights(min_pair.first) ), &unity, sizeof(XFLOAT), cudaMemcpyHostToDevice, 0));

			PassWeights[img_id].weights.cpToHost();
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));
#else
			deviceInitValue<XFLOAT>(PassWeights[img_id].weights, (XFLOAT)0.0);
			PassWeights[img_id].weights[min_pair.first] = (XFLOAT)1.0;
#endif

			my_significant_weight = 0.999;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_NR_SIGN) = (RFLOAT) 1.;
			if (exp_ipass==0) // TODO better memset, 0 => false , 1 => true
				for (int ihidden = 0; ihidden < XSIZE(op.Mcoarse_significant); ihidden++)
					if (DIRECT_A2D_ELEM(op.Mweight, img_id, ihidden) >= my_significant_weight)
						DIRECT_A2D_ELEM(op.Mcoarse_significant, img_id, ihidden) = true;
					else
						DIRECT_A2D_ELEM(op.Mcoarse_significant, img_id, ihidden) = false;
			else
			{
				std::pair<size_t, XFLOAT> max_pair = AccUtilities::getArgMaxOnDevice<XFLOAT>(PassWeights[img_id].weights);
				op.max_index[img_id].fineIdx = PassWeights[img_id].ihidden_overs[max_pair.first];
				op.max_weight[img_id] = max_pair.second;
			}

		}
		else
		{


			long int sumRedSize=0;
			for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
				sumRedSize+= (exp_ipass==0) ? ceilf((float)(sp.nr_dir*sp.nr_psi)/(float)SUMW_BLOCK_SIZE) : ceil((float)FPCMasks[img_id][exp_iclass].jobNum / (float)SUMW_BLOCK_SIZE);

			// loop through making translational priors for all classes this img_id - then copy all at once - then loop through kernel calls ( TODO: group kernel calls into one big kernel)
			CTIC(accMLO->timer,"get_offset_priors");

			for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{
				RFLOAT myprior_x, myprior_y, myprior_z;
				if (baseMLO->mymodel.nr_bodies > 1)
				{
					myprior_x = myprior_y = myprior_z = 0.;
				}
				else if (baseMLO->mymodel.ref_dim == 2 && !baseMLO->do_helical_refine)
				{
					myprior_x = XX(baseMLO->mymodel.prior_offset_class[exp_iclass]);
					myprior_y = YY(baseMLO->mymodel.prior_offset_class[exp_iclass]);
				}
				else
				{
					myprior_x = XX(op.prior[img_id]);
					myprior_y = YY(op.prior[img_id]);
					if (accMLO->dataIs3D)
						myprior_z = ZZ(op.prior[img_id]);
				}

				for (unsigned long itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++)
				{

					// If it is doing helical refinement AND Cartesian vector myprior has a length > 0, transform the vector to its helical coordinates
					if ( (baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry))
					{
						RFLOAT mypriors_len2 = myprior_x * myprior_x + myprior_y * myprior_y;
						if (accMLO->dataIs3D)
							mypriors_len2 += myprior_z * myprior_z;

						if (mypriors_len2 > 0.00001)
						{
							RFLOAT rot_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT);
							RFLOAT tilt_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT);
							RFLOAT psi_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI);
							transformCartesianAndHelicalCoords(myprior_x, myprior_y, myprior_z, myprior_x, myprior_y, myprior_z, rot_deg, tilt_deg, psi_deg, (accMLO->dataIs3D) ? (3) : (2), CART_TO_HELICAL_COORDS);
						}
					}
					// (For helical refinement) Now offset, old_offset, sampling.translations and myprior are all in helical coordinates

					// To speed things up, only calculate pdf_offset at the coarse sampling.
					// That should not matter much, and that way one does not need to calculate all the OversampledTranslations
					double pdf(0), pdf_zeros(0);
					RFLOAT offset_x = old_offset_x + baseMLO->sampling.translations_x[itrans];
					RFLOAT offset_y = old_offset_y + baseMLO->sampling.translations_y[itrans];
					double tdiff2 = 0.;

					if ( (! baseMLO->do_helical_refine) || (baseMLO->ignore_helical_symmetry) || (accMLO->dataIs3D) )
						tdiff2 += (offset_x - myprior_x) * (offset_x - myprior_x);
					tdiff2 += (offset_y - myprior_y) * (offset_y - myprior_y);
					if (accMLO->dataIs3D)
					{
						RFLOAT offset_z = old_offset_z + baseMLO->sampling.translations_z[itrans];
						if ( (! baseMLO->do_helical_refine) || (baseMLO->ignore_helical_symmetry) )
							tdiff2 += (offset_z - myprior_z) * (offset_z - myprior_z);
					}

					// As of version 3.1, sigma_offsets are in Angstroms!
					tdiff2 *= my_pixel_size * my_pixel_size;

					// P(offset|sigma2_offset)
					// This is the probability of the offset, given the model offset and variance.
					if (my_sigma2_offset < 0.0001)
					{
						pdf_zeros = tdiff2 > 0.;
						pdf = pdf_zeros ? 0. : 1.;

					}
					else
					{
						pdf_zeros = false;
						pdf = tdiff2 / (-2. * my_sigma2_offset);
					}

					pdf_offset_zeros[(exp_iclass-sp.iclass_min)*sp.nr_trans + itrans] = pdf_zeros;
					pdf_offset     [(exp_iclass-sp.iclass_min)*sp.nr_trans + itrans] = pdf;
				}
			}

			pdf_offset_zeros.cpToDevice();
			pdf_offset.cpToDevice();

			CTOC(accMLO->timer,"get_offset_priors");
			CTIC(accMLO->timer,"sumweight1");

			if(exp_ipass==0)
			{
				AccPtr<XFLOAT>  ipartMweight(
						Mweight,
						img_id * op.Mweight.xdim + sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.iclass_min,
						(sp.iclass_max-sp.iclass_min+1) * sp.nr_dir * sp.nr_psi * sp.nr_trans);

				pdf_offset.streamSync();

				AccUtilities::kernel_weights_exponent_coarse(
						sp.iclass_max-sp.iclass_min+1,
						pdf_orientation,
						pdf_orientation_zeros,
						pdf_offset,
						pdf_offset_zeros,
						ipartMweight,
						(XFLOAT)op.min_diff2[img_id],
						sp.nr_dir*sp.nr_psi,
						sp.nr_trans);


				XFLOAT weights_max = AccUtilities::getMaxOnDevice<XFLOAT>(ipartMweight);

				/*
				 * Add 50 since we want to stay away from e^88, which approaches the single precision limit.
				 * We still want as high numbers as possible to utilize most of the single precision span.
				 * Dari - 201710
				*/
				AccUtilities::kernel_exponentiate( ipartMweight, 50 - weights_max);

				CTIC(accMLO->timer,"sort");
				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

				unsigned long ipart_length = (sp.iclass_max-sp.iclass_min+1) * sp.nr_dir * sp.nr_psi * sp.nr_trans;
				size_t offset = img_id * op.Mweight.xdim + sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.iclass_min;

				if (ipart_length > 1)
				{
					//Wrap the current ipart data in a new pointer
					AccPtr<XFLOAT> unsorted_ipart(
							Mweight,
							offset,
							ipart_length);

					AccPtr<XFLOAT> filtered = ptrFactory.make<XFLOAT>((size_t)unsorted_ipart.getSize());

					CUSTOM_ALLOCATOR_REGION_NAME("CASDTW_SORTSUM");

					filtered.deviceAlloc();

#ifdef DEBUG_CUDA
					if (unsorted_ipart.getSize()==0)
						ACC_PTR_DEBUG_FATAL("Unsorted array size zero.\n");  // Hopefully Impossible
#endif
					size_t filteredSize = AccUtilities::filterGreaterZeroOnDevice<XFLOAT>(unsorted_ipart, filtered);

					if (filteredSize == 0)
					{
						std::cerr << std::endl;
						std::cerr << " fn_img= " << sp.current_img << std::endl;
						std::cerr << " img_id= " << img_id << " adaptive_fraction= " << baseMLO->adaptive_fraction << std::endl;
						std::cerr << " min_diff2= " << op.min_diff2[img_id] << std::endl;

						pdf_orientation.dumpAccToFile("error_dump_pdf_orientation");
						pdf_offset.dumpAccToFile("error_dump_pdf_offset");
						unsorted_ipart.dumpAccToFile("error_dump_filtered");

						std::cerr << "Dumped data: error_dump_pdf_orientation, error_dump_pdf_orientation and error_dump_unsorted." << std::endl;

						CRITICAL(ERRFILTEREDZERO); // "filteredSize == 0"
					}
					filtered.setSize(filteredSize);

					AccPtr<XFLOAT> sorted =         ptrFactory.make<XFLOAT>((size_t)filteredSize);
					AccPtr<XFLOAT> cumulative_sum = ptrFactory.make<XFLOAT>((size_t)filteredSize);

					sorted.accAlloc();
					cumulative_sum.accAlloc();

					AccUtilities::sortOnDevice<XFLOAT>(filtered, sorted);
					AccUtilities::scanOnDevice<XFLOAT>(sorted, cumulative_sum);

					CTOC(accMLO->timer,"sort");

					op.sum_weight[img_id] = cumulative_sum.getAccValueAt(cumulative_sum.getSize() - 1);

					long int my_nr_significant_coarse_samples;
					size_t thresholdIdx = findThresholdIdxInCumulativeSum<XFLOAT>(cumulative_sum,
							(1 - baseMLO->adaptive_fraction) * op.sum_weight[img_id]);

					my_nr_significant_coarse_samples = filteredSize - thresholdIdx;

					if (my_nr_significant_coarse_samples == 0)
					{
						std::cerr << std::endl;
						std::cerr << " fn_img= " << sp.current_img << std::endl;
						std::cerr << " img_id= " << img_id << " adaptive_fraction= " << baseMLO->adaptive_fraction << std::endl;
						std::cerr << " threshold= " << (1 - baseMLO->adaptive_fraction) * op.sum_weight[img_id] << " thresholdIdx= " << thresholdIdx << std::endl;
						std::cerr << " op.sum_weight[img_id]= " << op.sum_weight[img_id] << std::endl;
						std::cerr << " min_diff2= " << op.min_diff2[img_id] << std::endl;

						unsorted_ipart.dumpAccToFile("error_dump_unsorted");
						filtered.dumpAccToFile("error_dump_filtered");
						sorted.dumpAccToFile("error_dump_sorted");
						cumulative_sum.dumpAccToFile("error_dump_cumulative_sum");

						std::cerr << "Written error_dump_unsorted, error_dump_filtered, error_dump_sorted, and error_dump_cumulative_sum." << std::endl;

						CRITICAL(ERRNOSIGNIFS); // "my_nr_significant_coarse_samples == 0"
					}

					if (baseMLO->maximum_significants > 0 &&
							my_nr_significant_coarse_samples > baseMLO->maximum_significants)
					{
						my_nr_significant_coarse_samples = baseMLO->maximum_significants;
						thresholdIdx = filteredSize - my_nr_significant_coarse_samples;
					}

					XFLOAT significant_weight = sorted.getAccValueAt(thresholdIdx);

					CTIC(accMLO->timer,"getArgMaxOnDevice");
					std::pair<size_t, XFLOAT> max_pair = AccUtilities::getArgMaxOnDevice<XFLOAT>(unsorted_ipart);
					CTOC(accMLO->timer,"getArgMaxOnDevice");
					op.max_index[img_id].coarseIdx = max_pair.first;
					op.max_weight[img_id] = max_pair.second;

					// Store nr_significant_coarse_samples for this particle
					// Don't do this for multibody, as it would be overwritten for each body,
					// and we also use METADATA_NR_SIGN in the new safeguard for the gold-standard separation
					if (baseMLO->mymodel.nr_bodies == 1)
						DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_NR_SIGN) = (RFLOAT) my_nr_significant_coarse_samples;

					AccPtr<bool> Mcoarse_significant = ptrFactory.make<bool>(ipart_length);
					Mcoarse_significant.setHostPtr(&op.Mcoarse_significant.data[offset]);

					CUSTOM_ALLOCATOR_REGION_NAME("CASDTW_SIG");
					Mcoarse_significant.deviceAlloc();

					DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));
					arrayOverThreshold<XFLOAT>(unsorted_ipart, Mcoarse_significant, significant_weight);
					Mcoarse_significant.cpToHost();
					DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));
				}
				else if (ipart_length == 1)
				{
					op.Mcoarse_significant.data[img_id * op.Mweight.xdim + sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.iclass_min] = 1;
				}
				else
					CRITICAL(ERRNEGLENGTH);
			}
			else
			{

				for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
					DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[exp_iclass]));
				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

				XFLOAT weights_max = -std::numeric_limits<XFLOAT>::max();

				pdf_offset.streamSync();

				for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++) // TODO could use classStreams
				{
					if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) && (FPCMasks[img_id][exp_iclass].weightNum > 0) )
					{
						// Use the constructed mask to build a partial (class-specific) input
						// (until now, PassWeights has been an empty placeholder. We now create class-partials pointing at it, and start to fill it with stuff)

						IndexedDataArray thisClassPassWeights(PassWeights[img_id],FPCMasks[img_id][exp_iclass]);

						AccPtr<XFLOAT> pdf_orientation_class =       ptrFactory.make<XFLOAT>(sp.nr_dir*sp.nr_psi),
						               pdf_offset_class =            ptrFactory.make<XFLOAT>(sp.nr_trans);
						AccPtr<bool>   pdf_orientation_zeros_class = ptrFactory.make<bool>(sp.nr_dir*sp.nr_psi),
						               pdf_offset_zeros_class =      ptrFactory.make<bool>(sp.nr_trans);

						pdf_orientation_class      .setAccPtr(&((~pdf_orientation)      [(exp_iclass-sp.iclass_min)*sp.nr_dir*sp.nr_psi]));
						pdf_orientation_zeros_class.setAccPtr(&((~pdf_orientation_zeros)[(exp_iclass-sp.iclass_min)*sp.nr_dir*sp.nr_psi]));

						pdf_offset_class           .setAccPtr(&((~pdf_offset)           [(exp_iclass-sp.iclass_min)*sp.nr_trans]));
						pdf_offset_zeros_class     .setAccPtr(&((~pdf_offset_zeros)     [(exp_iclass-sp.iclass_min)*sp.nr_trans]));

						thisClassPassWeights.weights.setStream(accMLO->classStreams[exp_iclass]);

						AccUtilities::kernel_exponentiate_weights_fine(
								~pdf_orientation_class,
								~pdf_orientation_zeros_class,
								~pdf_offset_class,
								~pdf_offset_zeros_class,
								~thisClassPassWeights.weights,
								(XFLOAT)op.min_diff2[img_id],
								sp.nr_oversampled_rot,
								sp.nr_oversampled_trans,
								~thisClassPassWeights.rot_id,
								~thisClassPassWeights.trans_idx,
								~FPCMasks[img_id][exp_iclass].jobOrigin,
								~FPCMasks[img_id][exp_iclass].jobExtent,
								FPCMasks[img_id][exp_iclass].jobNum,
								accMLO->classStreams[exp_iclass]);

						XFLOAT m = AccUtilities::getMaxOnDevice<XFLOAT>(thisClassPassWeights.weights);

						if (m > weights_max)
							weights_max = m;
					}
				}

				for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++) // TODO could use classStreams
				{
					if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) && (FPCMasks[img_id][exp_iclass].weightNum > 0) )
					{
						IndexedDataArray thisClassPassWeights(PassWeights[img_id],FPCMasks[img_id][exp_iclass]);

						thisClassPassWeights.weights.setStream(accMLO->classStreams[exp_iclass]);
						/*
						 * Add 50 since we want to stay away from e^88, which approaches the single precision limit.
						 * We still want as high numbers as possible to utilize most of the single precision span.
						 * Dari - 201710
						*/
						AccUtilities::kernel_exponentiate( thisClassPassWeights.weights, 50 - weights_max );
					}
				}

				op.min_diff2[img_id] += 50 - weights_max;

				for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
					DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[exp_iclass]));
				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

				if(baseMLO->is_som_iter) {
					op.sum_weight_class[img_id].resize(baseMLO->mymodel.nr_classes, 0);

					for (unsigned long exp_iclass = sp.iclass_min;
					     exp_iclass <= sp.iclass_max; exp_iclass++) // TODO could use classStreams
					{
						if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) &&
						    (FPCMasks[img_id][exp_iclass].weightNum > 0)) {
							IndexedDataArray thisClassPassWeights(PassWeights[img_id], FPCMasks[img_id][exp_iclass]);
							op.sum_weight_class[img_id][exp_iclass] = AccUtilities::getSumOnDevice(thisClassPassWeights.weights);
						}
					}
				}

				PassWeights[img_id].weights.cpToHost(); // note that the host-pointer is shared: we're copying to Mweight.


				CTIC(accMLO->timer,"sort");
				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));
				size_t weightSize = PassWeights[img_id].weights.getSize();

				AccPtr<XFLOAT> sorted =         ptrFactory.make<XFLOAT>((size_t)weightSize);
				AccPtr<XFLOAT> cumulative_sum = ptrFactory.make<XFLOAT>((size_t)weightSize);

				CUSTOM_ALLOCATOR_REGION_NAME("CASDTW_FINE");

				sorted.accAlloc();
				cumulative_sum.accAlloc();

				AccUtilities::sortOnDevice<XFLOAT>(PassWeights[img_id].weights, sorted);
				AccUtilities::scanOnDevice<XFLOAT>(sorted, cumulative_sum);
				CTOC(accMLO->timer,"sort");

				if(baseMLO->adaptive_oversampling!=0)
				{
					op.sum_weight[img_id] = cumulative_sum.getAccValueAt(cumulative_sum.getSize() - 1);

					if (op.sum_weight[img_id]==0)
					{
						std::cerr << std::endl;
						std::cerr << " fn_img= " << sp.current_img << std::endl;
						std::cerr << " op.part_id= " << op.part_id << std::endl;
						std::cerr << " img_id= " << img_id << std::endl;
						std::cerr << " op.min_diff2[img_id]= " << op.min_diff2[img_id] << std::endl;
						int group_id = baseMLO->mydata.getGroupId(op.part_id, img_id);
						std::cerr << " group_id= " << group_id << std::endl;
						int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);
						std::cerr << " optics_group= " << optics_group << std::endl;
						std::cerr << " ml_model.scale_correction[group_id]= " << baseMLO->mymodel.scale_correction[group_id] << std::endl;
						std::cerr << " exp_significant_weight[img_id]= " << op.significant_weight[img_id] << std::endl;
						std::cerr << " exp_max_weight[img_id]= " << op.max_weight[img_id] << std::endl;
						std::cerr << " ml_model.sigma2_noise[optics_group]= " << baseMLO->mymodel.sigma2_noise[optics_group] << std::endl;
						CRITICAL(ERRSUMWEIGHTZERO); //"op.sum_weight[img_id]==0"
					}

					size_t thresholdIdx = findThresholdIdxInCumulativeSum<XFLOAT>(cumulative_sum, (1 - baseMLO->adaptive_fraction) * op.sum_weight[img_id]);
					my_significant_weight = sorted.getAccValueAt(thresholdIdx);

					CTIC(accMLO->timer,"getArgMaxOnDevice");
					std::pair<size_t, XFLOAT> max_pair = AccUtilities::getArgMaxOnDevice<XFLOAT>(PassWeights[img_id].weights);
					CTOC(accMLO->timer,"getArgMaxOnDevice");
					op.max_index[img_id].fineIdx = PassWeights[img_id].ihidden_overs[max_pair.first];
					op.max_weight[img_id] = max_pair.second;
				}
				else
				{
					my_significant_weight = sorted.getAccValueAt(0);
				}
			}
			CTOC(accMLO->timer,"sumweight1");
		}

		op.significant_weight[img_id] = (RFLOAT) my_significant_weight;
		} // end loop img_id


#ifdef TIMING
	if (op.part_id == baseMLO->exp_my_first_part_id)
	{
		if (exp_ipass == 0) baseMLO->timer.toc(baseMLO->TIMING_ESP_WEIGHT1);
		else baseMLO->timer.toc(baseMLO->TIMING_ESP_WEIGHT2);
	}
#endif
}

// ----------------------------------------------------------------------------
// -------------------------- storeWeightedSums -------------------------------
// ----------------------------------------------------------------------------
template<class MlClass>
void storeWeightedSums(OptimisationParamters &op, SamplingParameters &sp,
						MlOptimiser *baseMLO,
						MlClass *accMLO,
						std::vector<IndexedDataArray > &FinePassWeights,
						std::vector<ProjectionParams> &ProjectionData,
						std::vector<std::vector<IndexedDataArrayMask > > &FPCMasks,
						AccPtrFactory ptrFactory,
						int ibody,
						std::vector< AccPtrBundle > &bundleSWS)
{
#ifdef TIMING
	if (op.part_id == baseMLO->exp_my_first_part_id)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_WSUM);
#endif
	CTIC(accMLO->timer,"store_init");

	// Re-do below because now also want unmasked images AND if (stricht_highres_exp >0.) then may need to resize
	std::vector<MultidimArray<Complex > > dummy;
	std::vector<std::vector<MultidimArray<Complex > > > dummy2;
	std::vector<MultidimArray<RFLOAT> > exp_local_STMulti;
	bool do_subtomo_correction = op.FstMulti.size() > 0 && NZYXSIZE(op.FstMulti[0]) > 0;
	if (do_subtomo_correction)
		exp_local_STMulti.resize(sp.nr_images);

	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, true, op.part_id, sp.current_oversampling, op.metadata_offset, // inserted SHWS 12112015
			sp.itrans_min, sp.itrans_max, op.Fimg, op.Fimg_nomask, op.Fctf, dummy2, dummy2,
			op.local_Fctf, op.local_sqrtXi2, op.local_Minvsigma2, op.FstMulti, exp_local_STMulti);

	// In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the op.local_Minvsigma2s was omitted.
	// Set those back here
	for (int img_id = 0; img_id < sp.nr_images; img_id++)
	{
		int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);
		DIRECT_MULTIDIM_ELEM(op.local_Minvsigma2[img_id], 0) = 1. / (baseMLO->sigma2_fudge * DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[optics_group], 0));
	}

	// For norm_correction and scale_correction of all images of this particle
	std::vector<RFLOAT> exp_wsum_norm_correction;
	std::vector<RFLOAT> exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA;
	std::vector<RFLOAT> thr_wsum_signal_product_spectra, thr_wsum_reference_power_spectra;
	exp_wsum_norm_correction.resize(sp.nr_images, 0.);
	std::vector<MultidimArray<RFLOAT> > thr_wsum_sigma2_noise, thr_wsum_ctf2, thr_wsum_stMulti;

	// for noise estimation (per image)
	thr_wsum_sigma2_noise.resize(sp.nr_images);
    thr_wsum_ctf2.resize(sp.nr_images);
    thr_wsum_stMulti.resize(sp.nr_images);

	// For scale_correction
	if (baseMLO->do_scale_correction)
	{
		exp_wsum_scale_correction_XA.resize(sp.nr_images);
		exp_wsum_scale_correction_AA.resize(sp.nr_images);
		thr_wsum_signal_product_spectra.resize(sp.nr_images);
		thr_wsum_reference_power_spectra.resize(sp.nr_images);
	}

	// Possibly different array sizes in different optics groups!
	for (int img_id = 0; img_id < sp.nr_images; img_id++)
	{
		int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);
		thr_wsum_sigma2_noise[img_id].initZeros(baseMLO->image_full_size[optics_group]/2 + 1);
        thr_wsum_stMulti[img_id].initZeros(baseMLO->image_full_size[optics_group]/2 + 1);
        thr_wsum_ctf2[img_id].initZeros(baseMLO->image_full_size[optics_group]/2 + 1);
		if (baseMLO->do_scale_correction)
		{
			exp_wsum_scale_correction_AA[img_id] = 0.;
			exp_wsum_scale_correction_XA[img_id] = 0.;
			thr_wsum_signal_product_spectra[img_id] = 0.;
			thr_wsum_reference_power_spectra[img_id] = 0.;
		}
	}

	std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
	bool have_warned_small_scale = false;

	// Make local copies of weighted sums (except BPrefs, which are too big)
	// so that there are not too many mutex locks below
	std::vector<MultidimArray<RFLOAT> > thr_wsum_pdf_direction;
	std::vector<RFLOAT> thr_wsum_norm_correction, thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class;
	RFLOAT thr_wsum_sigma2_offset;
	MultidimArray<RFLOAT> thr_metadata, zeroArray;
	// wsum_pdf_direction is a 1D-array (of length sampling.NrDirections()) for each class
	zeroArray.initZeros(baseMLO->sampling.NrDirections());
	thr_wsum_pdf_direction.resize(baseMLO->mymodel.nr_classes * baseMLO->mymodel.nr_bodies, zeroArray);
	// sumw_group is a RFLOAT for each group
	thr_sumw_group.resize(sp.nr_images, 0.);
	// wsum_pdf_class is a RFLOAT for each class
	thr_wsum_pdf_class.resize(baseMLO->mymodel.nr_classes, 0.);
	if (baseMLO->mymodel.ref_dim == 2)
	{
		thr_wsum_prior_offsetx_class.resize(baseMLO->mymodel.nr_classes, 0.);
		thr_wsum_prior_offsety_class.resize(baseMLO->mymodel.nr_classes, 0.);
	}
	// wsum_sigma2_offset is just a RFLOAT
	thr_wsum_sigma2_offset = 0.;
	CTOC(accMLO->timer,"store_init");

	/*=======================================================================================
	                           COLLECT 2 AND SET METADATA
	=======================================================================================*/

	CTIC(accMLO->timer,"collect_data_2");
	unsigned long nr_transes = sp.nr_trans*sp.nr_oversampled_trans;
	unsigned long nr_fake_classes = (sp.iclass_max-sp.iclass_min+1);
	unsigned long oversamples = sp.nr_oversampled_trans * sp.nr_oversampled_rot;
	std::vector<long int> block_nums(sp.nr_images*nr_fake_classes);

	for (int img_id = 0; img_id < sp.nr_images; img_id++)
	{
		// here we introduce offsets for the oo_transes in an array as it is more efficient to
		// copy one big array to/from GPU rather than four small arrays
		size_t otrans_x      = 0*(size_t)nr_fake_classes*nr_transes;
		size_t otrans_y      = 1*(size_t)nr_fake_classes*nr_transes;
		size_t otrans_z      = 2*(size_t)nr_fake_classes*nr_transes;
		size_t otrans_x2y2z2 = 3*(size_t)nr_fake_classes*nr_transes;

		// Allocate space for all classes, so that we can pre-calculate data for all classes, copy in one operation, call kenrels on all classes, and copy back in one operation
		AccPtr<XFLOAT>          oo_otrans = ptrFactory.make<XFLOAT>((size_t)nr_fake_classes*nr_transes*4);

		oo_otrans.allAlloc();

		int sumBlockNum =0;
		int my_metadata_offset = op.metadata_offset + img_id;
		int group_id = baseMLO->mydata.getGroupId(op.part_id, img_id);
		const int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);
		RFLOAT my_pixel_size = baseMLO->mydata.getImagePixelSize(op.part_id, img_id);

		CTIC(accMLO->timer,"collect_data_2_pre_kernel");
		for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			unsigned long fake_class = exp_iclass-sp.iclass_min; // if we only have the third class to do, the third class will be the "first" we do, i.e. the "fake" first.
			if ((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[img_id].class_entries[exp_iclass] == 0) )
				continue;

			// Use the constructed mask to construct a partial class-specific input
			IndexedDataArray thisClassFinePassWeights(FinePassWeights[img_id],FPCMasks[img_id][exp_iclass]);

			// Re-define the job-partition of the indexedArray of weights so that the collect-kernel can work with it.
			block_nums[nr_fake_classes*img_id + fake_class] = makeJobsForCollect(thisClassFinePassWeights, FPCMasks[img_id][exp_iclass], ProjectionData[img_id].orientation_num[exp_iclass]);

			bundleSWS[img_id].pack(FPCMasks[img_id][exp_iclass].jobOrigin);
			bundleSWS[img_id].pack(FPCMasks[img_id][exp_iclass].jobExtent);

			sumBlockNum+=block_nums[nr_fake_classes*img_id + fake_class];

			RFLOAT myprior_x, myprior_y, myprior_z, old_offset_z;
			RFLOAT old_offset_x = XX(op.old_offset[img_id]);
			RFLOAT old_offset_y = YY(op.old_offset[img_id]);

			if (baseMLO->mymodel.ref_dim == 2 && baseMLO->mymodel.nr_bodies == 1)
			{
				myprior_x = XX(baseMLO->mymodel.prior_offset_class[exp_iclass]);
				myprior_y = YY(baseMLO->mymodel.prior_offset_class[exp_iclass]);
			}
			else
			{
				myprior_x = XX(op.prior[img_id]);
				myprior_y = YY(op.prior[img_id]);
				if (baseMLO->mymodel.data_dim == 3)
				{
					myprior_z = ZZ(op.prior[img_id]);
					old_offset_z = ZZ(op.old_offset[img_id]);
				}
			}

			/*======================================================
								COLLECT 2
			======================================================*/

			//Pregenerate oversampled translation objects for kernel-call
			for (long int itrans = 0, iitrans = 0; itrans < sp.nr_trans; itrans++)
			{
				baseMLO->sampling.getTranslationsInPixel(itrans, baseMLO->adaptive_oversampling, my_pixel_size,
						oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
						(baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry));
				for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++, iitrans++)
				{
					oo_otrans[otrans_x+fake_class*nr_transes+iitrans] = old_offset_x + oversampled_translations_x[iover_trans];
					oo_otrans[otrans_y+fake_class*nr_transes+iitrans] = old_offset_y + oversampled_translations_y[iover_trans];
					if (accMLO->dataIs3D)
						oo_otrans[otrans_z+fake_class*nr_transes+iitrans] = old_offset_z + oversampled_translations_z[iover_trans];

					// Calculate the vector length of myprior
					RFLOAT mypriors_len2 = myprior_x * myprior_x + myprior_y * myprior_y;
					if (accMLO->dataIs3D)
						mypriors_len2 += myprior_z * myprior_z;

					// If it is doing helical refinement AND Cartesian vector myprior has a length > 0, transform the vector to its helical coordinates
					if ( (baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry) && (mypriors_len2 > 0.00001) )
					{
						RFLOAT rot_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT);
						RFLOAT tilt_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT);
						RFLOAT psi_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI);
						transformCartesianAndHelicalCoords(myprior_x, myprior_y, myprior_z, myprior_x, myprior_y, myprior_z, rot_deg, tilt_deg, psi_deg, (accMLO->dataIs3D) ? (3) : (2), CART_TO_HELICAL_COORDS);
					}

					if ( (! baseMLO->do_helical_refine) || (baseMLO->ignore_helical_symmetry) )
						RFLOAT diffx = myprior_x - oo_otrans[otrans_x+fake_class*nr_transes+iitrans];
					RFLOAT diffx = myprior_x - oo_otrans[otrans_x+fake_class*nr_transes+iitrans];
					RFLOAT diffy = myprior_y - oo_otrans[otrans_y+fake_class*nr_transes+iitrans];
					RFLOAT diffz = 0;
					if (accMLO->dataIs3D)
						diffz = myprior_z - (old_offset_z + oversampled_translations_z[iover_trans]);

					oo_otrans[otrans_x2y2z2+fake_class*nr_transes+iitrans] = diffx*diffx + diffy*diffy + diffz*diffz;
				}
			}
		}

		bundleSWS[img_id].cpToDevice();
		oo_otrans.cpToDevice();

		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

		// here we introduce offsets for the clases in an array as it is more efficient to
		// copy one big array to/from GPU rather than four small arrays
		size_t offsetx_class = 0*(size_t)sumBlockNum;
		size_t offsety_class = 1*(size_t)sumBlockNum;
		size_t offsetz_class = 2*(size_t)sumBlockNum;
		size_t sigma2_offset = 3*(size_t)sumBlockNum;

		AccPtr<XFLOAT>                      p_weights = ptrFactory.make<XFLOAT>((size_t)sumBlockNum);
		AccPtr<XFLOAT> p_thr_wsum_prior_offsetxyz_class = ptrFactory.make<XFLOAT>((size_t)sumBlockNum*4);

		p_weights.allAlloc();
		p_thr_wsum_prior_offsetxyz_class.allAlloc();
		CTOC(accMLO->timer,"collect_data_2_pre_kernel");
		int partial_pos=0;

		for (long int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			long int fake_class = exp_iclass-sp.iclass_min; // if we only have the third class to do, the third class will be the "first" we do, i.e. the "fake" first.
			if ((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[img_id].class_entries[exp_iclass] == 0) )
				continue;

			// Use the constructed mask to construct a partial class-specific input
			IndexedDataArray thisClassFinePassWeights(FinePassWeights[img_id],FPCMasks[img_id][exp_iclass]);

			long int cpos=fake_class*nr_transes;
			int block_num = block_nums[nr_fake_classes*img_id + fake_class];

			runCollect2jobs(block_num,
						&(~oo_otrans)[otrans_x+cpos],          // otrans-size -> make const
						&(~oo_otrans)[otrans_y+cpos],          // otrans-size -> make const
						&(~oo_otrans)[otrans_z+cpos],          // otrans-size -> make const
						&(~oo_otrans)[otrans_x2y2z2+cpos], // otrans-size -> make const
						~thisClassFinePassWeights.weights,
						(XFLOAT)op.significant_weight[img_id],
						(XFLOAT)op.sum_weight[img_id],
						sp.nr_trans,
						sp.nr_oversampled_trans,
						sp.nr_oversampled_rot,
						oversamples,
						(baseMLO->do_skip_align || baseMLO->do_skip_rotate ),
						&(~p_weights)[partial_pos],
						&(~p_thr_wsum_prior_offsetxyz_class)[offsetx_class+partial_pos],
						&(~p_thr_wsum_prior_offsetxyz_class)[offsety_class+partial_pos],
						&(~p_thr_wsum_prior_offsetxyz_class)[offsetz_class+partial_pos],
						&(~p_thr_wsum_prior_offsetxyz_class)[sigma2_offset+partial_pos],
						~thisClassFinePassWeights.rot_idx,
						~thisClassFinePassWeights.trans_idx,
						~FPCMasks[img_id][exp_iclass].jobOrigin,
						~FPCMasks[img_id][exp_iclass].jobExtent,
						accMLO->dataIs3D);
			LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

			partial_pos+=block_num;
		}

		CTIC(accMLO->timer,"collect_data_2_post_kernel");
		p_weights.cpToHost();
		p_thr_wsum_prior_offsetxyz_class.cpToHost();

		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));
		int iorient = 0;
		partial_pos=0;
		for (long int iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
		{
			long int fake_class = iclass-sp.iclass_min; // if we only have the third class to do, the third class will be the "first" we do, i.e. the "fake" first.
			if ((baseMLO->mymodel.pdf_class[iclass] == 0.) || (ProjectionData[img_id].class_entries[iclass] == 0) )
				continue;
			int block_num = block_nums[nr_fake_classes*img_id + fake_class];

			for (long int n = partial_pos; n < partial_pos+block_num; n++)
			{
				iorient= FinePassWeights[img_id].rot_id[FPCMasks[img_id][iclass].jobOrigin[n-partial_pos]+FPCMasks[img_id][iclass].firstPos];

				long int mydir, idir=floor(iorient/sp.nr_psi);
				if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
					mydir = idir;
				else
					mydir = op.pointer_dir_nonzeroprior[idir];

				// store partials according to indices of the relevant dimension
				unsigned ithr_wsum_pdf_direction = baseMLO->mymodel.nr_bodies > 1 ? ibody : iclass;
				DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[ithr_wsum_pdf_direction], mydir) += p_weights[n];
				thr_sumw_group[img_id]                                                       += p_weights[n];
				thr_wsum_pdf_class[iclass]                                                   += p_weights[n];

				thr_wsum_sigma2_offset                                                       += my_pixel_size * my_pixel_size * p_thr_wsum_prior_offsetxyz_class[sigma2_offset+n];

				if (baseMLO->mymodel.ref_dim == 2)
				{
					thr_wsum_prior_offsetx_class[iclass] += my_pixel_size * p_thr_wsum_prior_offsetxyz_class[offsetx_class+n];
					thr_wsum_prior_offsety_class[iclass] += my_pixel_size * p_thr_wsum_prior_offsetxyz_class[offsety_class+n];
				}
			}
			partial_pos+=block_num;
		} // end loop iclass
		CTOC(accMLO->timer,"collect_data_2_post_kernel");
	} // end loop img_id

	/*======================================================
	                     SET METADATA
	======================================================*/

	std::vector< RFLOAT> oversampled_rot, oversampled_tilt, oversampled_psi;
	for (long int img_id = 0; img_id < sp.nr_images; img_id++)
	{
		int my_metadata_offset = op.metadata_offset + img_id;
		RFLOAT my_pixel_size = baseMLO->mydata.getImagePixelSize(op.part_id, img_id);

		CTIC(accMLO->timer,"setMetadata");

		if(baseMLO->adaptive_oversampling!=0)
			op.max_index[img_id].fineIndexToFineIndices(sp); // set partial indices corresponding to the found max_index, to be used below
		else
			op.max_index[img_id].coarseIndexToCoarseIndices(sp);

		baseMLO->sampling.getTranslationsInPixel(op.max_index[img_id].itrans, baseMLO->adaptive_oversampling, my_pixel_size,
				oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
				(baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry));

		//TODO We already have rot, tilt and psi don't calculated them again
		if(baseMLO->do_skip_align || baseMLO->do_skip_rotate)
			   baseMLO->sampling.getOrientations(sp.idir_min, sp.ipsi_min, baseMLO->adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
					   op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);
		else
			   baseMLO->sampling.getOrientations(op.max_index[img_id].idir, op.max_index[img_id].ipsi, baseMLO->adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
					op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

		baseMLO->sampling.getOrientations(op.max_index[img_id].idir, op.max_index[img_id].ipsi, baseMLO->adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
				op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);

		RFLOAT rot = oversampled_rot[op.max_index[img_id].ioverrot];
		RFLOAT tilt = oversampled_tilt[op.max_index[img_id].ioverrot];
		RFLOAT psi = oversampled_psi[op.max_index[img_id].ioverrot];

		int icol_rot  = (baseMLO->mymodel.nr_bodies == 1) ? METADATA_ROT  : 0 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
		int icol_tilt = (baseMLO->mymodel.nr_bodies == 1) ? METADATA_TILT : 1 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
		int icol_psi  = (baseMLO->mymodel.nr_bodies == 1) ? METADATA_PSI  : 2 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
		int icol_xoff = (baseMLO->mymodel.nr_bodies == 1) ? METADATA_XOFF : 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
		int icol_yoff = (baseMLO->mymodel.nr_bodies == 1) ? METADATA_YOFF : 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;
		int icol_zoff = (baseMLO->mymodel.nr_bodies == 1) ? METADATA_ZOFF : 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + (ibody) * METADATA_NR_BODY_PARAMS;

		RFLOAT old_rot = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_rot);
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_rot) = rot;
		RFLOAT old_tilt = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_tilt);
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_tilt) = tilt;
		RFLOAT old_psi = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_psi);
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_psi) = psi;

		Matrix1D<RFLOAT> shifts(baseMLO->mymodel.data_dim);

		XX(shifts) = XX(op.old_offset[img_id]) + oversampled_translations_x[op.max_index[img_id].iovertrans];
		YY(shifts) = YY(op.old_offset[img_id]) + oversampled_translations_y[op.max_index[img_id].iovertrans];
		if (accMLO->dataIs3D)
		{
			ZZ(shifts) = ZZ(op.old_offset[img_id]) + oversampled_translations_z[op.max_index[img_id].iovertrans];
		}

		// Use oldpsi-angle to rotate back the XX(exp_old_offset[img_id]) + oversampled_translations_x[iover_trans] and
		if ( (baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry) )
			transformCartesianAndHelicalCoords(shifts, shifts, old_rot, old_tilt, old_psi, HELICAL_TO_CART_COORDS);

		DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_xoff) = XX(shifts);
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_yoff) = YY(shifts);
		if (accMLO->dataIs3D)
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, icol_zoff) = ZZ(shifts);

		if (ibody == 0)
		{
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_CLASS) = (RFLOAT)op.max_index[img_id].iclass + 1;
			RFLOAT pmax = op.max_weight[img_id]/op.sum_weight[img_id];
			if(pmax>1) //maximum normalised probability weight is (unreasonably) larger than unity
				CRITICAL("Relion is finding a normalised probability greater than 1");
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PMAX) = pmax;
		}
		CTOC(accMLO->timer,"setMetadata");
	}
	CTOC(accMLO->timer,"collect_data_2");


	/*=======================================================================================
	                                   MAXIMIZATION
	=======================================================================================*/

	if (!baseMLO->do_skip_maximization)
	{
		CTIC(accMLO->timer,"maximization");

		for (int img_id = 0; img_id < sp.nr_images; img_id++)
		{
			int my_metadata_offset = op.metadata_offset + img_id;
			int group_id = baseMLO->mydata.getGroupId(op.part_id, img_id);
			const int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);
			RFLOAT my_pixel_size = baseMLO->mydata.getImagePixelSize(op.part_id, img_id);
			bool ctf_premultiplied = baseMLO->mydata.obsModel.getCtfPremultiplied(optics_group);

			/*======================================================
			                     TRANSLATIONS
			======================================================*/

			long unsigned translation_num((sp.itrans_max - sp.itrans_min + 1) * sp.nr_oversampled_trans);

			size_t trans_x_offset = 0*(size_t)translation_num;
			size_t trans_y_offset = 1*(size_t)translation_num;
			size_t trans_z_offset = 2*(size_t)translation_num;

			AccPtr<XFLOAT> trans_xyz = ptrFactory.make<XFLOAT>((size_t)translation_num*3);

			trans_xyz.allAlloc();

			int j = 0;
			for (long int itrans = 0; itrans < (sp.itrans_max - sp.itrans_min + 1); itrans++)
			{
				//TODO Called multiple time to generate same list, reuse the same list
				baseMLO->sampling.getTranslationsInPixel(itrans, baseMLO->adaptive_oversampling, my_pixel_size,
						oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
						(baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry));

				for (long int iover_trans = 0; iover_trans < oversampled_translations_x.size(); iover_trans++)
				{
					RFLOAT xshift = 0., yshift = 0., zshift = 0.;

					xshift = oversampled_translations_x[iover_trans];
					yshift = oversampled_translations_y[iover_trans];
					if (accMLO->dataIs3D)
						zshift = oversampled_translations_z[iover_trans];

					if ( (baseMLO->do_helical_refine) && (! baseMLO->ignore_helical_symmetry) )
					{
						RFLOAT rot_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_ROT);
						RFLOAT tilt_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_TILT);
						RFLOAT psi_deg = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PSI);
						transformCartesianAndHelicalCoords(xshift, yshift, zshift, xshift, yshift, zshift, rot_deg, tilt_deg, psi_deg, (accMLO->dataIs3D) ? (3) : (2), HELICAL_TO_CART_COORDS);
					}

					trans_xyz[trans_x_offset+j] = -2 * PI * xshift / (double)baseMLO->image_full_size[optics_group];
					trans_xyz[trans_y_offset+j] = -2 * PI * yshift / (double)baseMLO->image_full_size[optics_group];
					trans_xyz[trans_z_offset+j] = -2 * PI * zshift / (double)baseMLO->image_full_size[optics_group];
					j ++;
				}
			}

			trans_xyz.cpToDevice();


			/*======================================================
			                     IMAGES
			======================================================*/

			CUSTOM_ALLOCATOR_REGION_NAME("TRANS_3");

			CTIC(accMLO->timer,"translation_3");

			MultidimArray<Complex > Fimg, Fimg_nonmask;
			windowFourierTransform(op.Fimg[img_id], Fimg, baseMLO->image_current_size[optics_group]); //TODO PO isen't this already done in getFourierTransformsAndCtfs?
			windowFourierTransform(op.Fimg_nomask[img_id], Fimg_nonmask, baseMLO->image_current_size[optics_group]);
			unsigned long image_size = Fimg.nzyxdim;

			size_t re_offset = 0*(size_t)image_size;
			size_t im_offset = 1*(size_t)image_size;
			size_t re_nomask_offset = 2*(size_t)image_size;
			size_t im_nomask_offset = 3*(size_t)image_size;

			AccPtr<XFLOAT> Fimgs = ptrFactory.make<XFLOAT>(4*(size_t)image_size);

			Fimgs.allAlloc();

			for (unsigned long i = 0; i < image_size; i ++)
			{
				Fimgs[re_offset+i] = Fimg.data[i].real;
				Fimgs[im_offset+i] = Fimg.data[i].imag;
				Fimgs[re_nomask_offset+i] = Fimg_nonmask.data[i].real;
				Fimgs[im_nomask_offset+i] = Fimg_nonmask.data[i].imag;
			}

			Fimgs.cpToDevice();

			CTOC(accMLO->timer,"translation_3");


			/*======================================================
			                       SCALE
			======================================================*/

			XFLOAT part_scale(1.);

			if (baseMLO->do_scale_correction)
			{
				part_scale = baseMLO->mymodel.scale_correction[group_id];
				if (part_scale > 10000.)
				{
					std::cerr << " rlnMicrographScaleCorrection= " << part_scale << " group= " << group_id + 1 << std::endl;
					CRITICAL(ERRHIGHSCALE);
				}
				else if (part_scale < 0.001)
				{
					if (!have_warned_small_scale)
					{
						std::cout << " WARNING: ignoring group " << group_id + 1 << " with very small or negative scale (" << part_scale <<
								"); Use larger groups for more stable scale estimates." << std::endl;
						have_warned_small_scale = true;
					}
					part_scale = 0.001;
				}
			}

			AccPtr<XFLOAT> ctfs = ptrFactory.make<XFLOAT>((size_t)image_size);
			ctfs.allAlloc();

			if (baseMLO->do_ctf_correction)
			{
				for (unsigned long i = 0; i < image_size; i++)
					ctfs[i] = (XFLOAT) op.local_Fctf[img_id].data[i] * part_scale;
			}
			else //TODO should be handled by memset
			{
				for (unsigned long i = 0; i < image_size; i++)
					ctfs[i] = part_scale;
			}

			ctfs.cpToDevice();

			/*======================================================
			                       MINVSIGMA
			======================================================*/

			AccPtr<XFLOAT> Minvsigma2s = ptrFactory.make<XFLOAT>((size_t)image_size);
			Minvsigma2s.allAlloc();

			if (baseMLO->do_map)
				for (unsigned long i = 0; i < image_size; i++)
					Minvsigma2s[i] = op.local_Minvsigma2[img_id].data[i];
			else
				for (unsigned long i = 0; i < image_size; i++)
					Minvsigma2s[i] = 1;

			Minvsigma2s.cpToDevice();

			/*======================================================
			                      CLASS LOOP
			======================================================*/

			CUSTOM_ALLOCATOR_REGION_NAME("wdiff2s");

			size_t wdiff2s_buf = (size_t)(baseMLO->mymodel.nr_classes*image_size)*2+(size_t)image_size;
			size_t AA_offset =  0*(size_t)(baseMLO->mymodel.nr_classes*image_size);
			size_t XA_offset =  1*(size_t)(baseMLO->mymodel.nr_classes*image_size);
			size_t sum_offset = 2*(size_t)(baseMLO->mymodel.nr_classes*image_size);

			AccPtr<XFLOAT> wdiff2s    = ptrFactory.make<XFLOAT>(wdiff2s_buf);

			wdiff2s.allAlloc();
			wdiff2s.accInit(0);

			unsigned long AAXA_pos=0;

			CUSTOM_ALLOCATOR_REGION_NAME("BP_data");

			// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
			AccPtr<XFLOAT> sorted_weights = ptrFactory.make<XFLOAT>((size_t)(ProjectionData[img_id].orientationNumAllClasses * translation_num));
			sorted_weights.allAlloc();
			std::vector<AccPtr<XFLOAT> > eulers(baseMLO->mymodel.nr_classes, ptrFactory.make<XFLOAT>());

			unsigned long classPos = 0;

			for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[exp_iclass]));
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

			for (unsigned long iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
			{
				if((baseMLO->mymodel.pdf_class[iclass] == 0.) || (ProjectionData[img_id].class_entries[iclass] == 0))
					continue;

				// Use the constructed mask to construct a partial class-specific input
				IndexedDataArray thisClassFinePassWeights(FinePassWeights[img_id],FPCMasks[img_id][iclass]);

				CTIC(accMLO->timer,"thisClassProjectionSetupCoarse");
				// use "slice" constructor with class-specific parameters to retrieve a temporary ProjectionParams with data for this class
				ProjectionParams thisClassProjectionData(	ProjectionData[img_id],
															ProjectionData[img_id].class_idx[iclass],
															ProjectionData[img_id].class_idx[iclass]+ProjectionData[img_id].class_entries[iclass]);

				thisClassProjectionData.orientation_num[0] = ProjectionData[img_id].orientation_num[iclass];
				CTOC(accMLO->timer,"thisClassProjectionSetupCoarse");

				long unsigned orientation_num(thisClassProjectionData.orientation_num[0]);

				/*======================================================
									PROJECTIONS
				======================================================*/


				Matrix2D<RFLOAT> MBL, MBR;

				if (baseMLO->mymodel.nr_bodies > 1)
				{
					Matrix2D<RFLOAT> Aori;
					RFLOAT rot_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_ROT);
					RFLOAT tilt_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_TILT);
					RFLOAT psi_ori = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset, METADATA_PSI);
					Euler_angles2matrix(rot_ori, tilt_ori, psi_ori, Aori, false);

					MBL = Aori * (baseMLO->mymodel.orient_bodies[ibody]).transpose() * baseMLO->A_rot90;
					MBR = baseMLO->mymodel.orient_bodies[ibody];
				}

				eulers[iclass].setSize(orientation_num * 9);
				eulers[iclass].setStream(accMLO->classStreams[iclass]);
				eulers[iclass].hostAlloc();

				CTIC(accMLO->timer,"generateEulerMatricesProjector");

				Matrix2D<RFLOAT> mag;
				mag.initIdentity(3);
				mag = baseMLO->mydata.obsModel.applyAnisoMag(mag, optics_group);
				mag = baseMLO->mydata.obsModel.applyScaleDifference(mag, optics_group, baseMLO->mymodel.ori_size, baseMLO->mymodel.pixel_size);
				if (!mag.isIdentity())
				{
					if (MBL.mdimx == 3 && MBL.mdimx ==3) MBL = mag * MBL;
					else MBL = mag;
				}

				generateEulerMatrices(
						thisClassProjectionData,
						&eulers[iclass][0],
						true,
						MBL,
						MBR);

				eulers[iclass].deviceAlloc();
				eulers[iclass].cpToDevice();

				CTOC(accMLO->timer,"generateEulerMatricesProjector");

				/*======================================================
									 MAP WEIGHTS
				======================================================*/

				CTIC(accMLO->timer,"pre_wavg_map");

				for (long unsigned i = 0; i < orientation_num*translation_num; i++)
					sorted_weights[classPos+i] = -std::numeric_limits<XFLOAT>::max();

				for (long unsigned i = 0; i < thisClassFinePassWeights.weights.getSize(); i++)
					sorted_weights[classPos+(thisClassFinePassWeights.rot_idx[i]) * translation_num + thisClassFinePassWeights.trans_idx[i] ]
									= thisClassFinePassWeights.weights[i];

				classPos+=orientation_num*translation_num;
				CTOC(accMLO->timer,"pre_wavg_map");
			}
			sorted_weights.cpToDevice();

			// These syncs are necessary (for multiple ranks on the same GPU), and (assumed) low-cost.
			for (unsigned long iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[iclass]));

			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

			classPos = 0;
			for (unsigned long iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
			{
				int iproj;
				if (baseMLO->mymodel.nr_bodies > 1) iproj = ibody;
				else                                iproj = iclass;

				if((baseMLO->mymodel.pdf_class[iclass] == 0.) || (ProjectionData[img_id].class_entries[iclass] == 0))
					continue;
				/*======================================================
									 KERNEL CALL
				======================================================*/

				long unsigned orientation_num(ProjectionData[img_id].orientation_num[iclass]);

				AccProjectorKernel projKernel = AccProjectorKernel::makeKernel(
						accMLO->bundle->projectors[iproj],
						op.local_Minvsigma2[img_id].xdim,
						op.local_Minvsigma2[img_id].ydim,
						op.local_Minvsigma2[img_id].zdim,
						op.local_Minvsigma2[img_id].xdim-1);

				runWavgKernel(
						projKernel,
						~eulers[iclass],
						&(~Fimgs)[re_offset], //~Fimgs_real,
						&(~Fimgs)[im_offset], //~Fimgs_imag,
						&(~trans_xyz)[trans_x_offset], //~trans_x,
						&(~trans_xyz)[trans_y_offset], //~trans_y,
						&(~trans_xyz)[trans_z_offset], //~trans_z,
						&(~sorted_weights)[classPos],
						~ctfs,
						&(~wdiff2s)[sum_offset],
						&(~wdiff2s)[AA_offset+AAXA_pos],
						&(~wdiff2s)[XA_offset+AAXA_pos],
						op,
						orientation_num,
						translation_num,
						image_size,
						img_id,
						group_id,
						iclass,
						part_scale,
						baseMLO->refs_are_ctf_corrected,
						ctf_premultiplied,
						accMLO->dataIs3D,
						accMLO->classStreams[iclass]);

				AAXA_pos += image_size;
				classPos += orientation_num*translation_num;
			}

			/*======================================================
								      SOM
			======================================================*/

			int nr_classes = baseMLO->mymodel.nr_classes;
			std::vector<RFLOAT> class_sum_weight(nr_classes, baseMLO->is_som_iter ? 0 : op.sum_weight[img_id]);

			if (baseMLO->is_som_iter)
			{
				std::vector<unsigned> s = SomGraph::arg_sort(op.sum_weight_class[img_id], false);
				unsigned bpu = s[0];
				unsigned sbpu = s[1];

				baseMLO->wsum_model.som.add_edge_activity(bpu, sbpu);

				class_sum_weight[bpu] = op.sum_weight_class[img_id][bpu];
				thr_wsum_pdf_class[bpu] += 1;
				baseMLO->wsum_model.som.add_node_activity(bpu);
				baseMLO->mymodel.som.add_node_age(bpu);

				std::vector<std::pair<unsigned, float> > weights = baseMLO->mymodel.som.get_neighbours(bpu);

				for (int i = 0; i < weights.size(); i++) {
					unsigned idx = weights[i].first;
					float w = weights[i].second * baseMLO->som_neighbour_pull;
					class_sum_weight[idx] = op.sum_weight_class[img_id][idx] / w;
					thr_wsum_pdf_class[idx] += w;
					baseMLO->wsum_model.som.add_node_activity(idx, w);
					baseMLO->mymodel.som.add_node_age(idx, w);
				}
			}

			/*======================================================
								BACKPROJECTION
			======================================================*/

			classPos = 0;
			for (unsigned long iclass = sp.iclass_min; iclass <= sp.iclass_max; iclass++)
			{

				int iproj;
				if (baseMLO->mymodel.nr_bodies > 1) iproj = ibody;
				else                                iproj = iclass;

				if((baseMLO->mymodel.pdf_class[iclass] == 0.) || (ProjectionData[img_id].class_entries[iclass] == 0))
					continue;

				if ( baseMLO->is_som_iter && class_sum_weight[iclass] == 0)
					continue;

				long unsigned orientation_num(ProjectionData[img_id].orientation_num[iclass]);

				AccProjectorKernel projKernel = AccProjectorKernel::makeKernel(
						accMLO->bundle->projectors[iproj],
						op.local_Minvsigma2[img_id].xdim,
						op.local_Minvsigma2[img_id].ydim,
						op.local_Minvsigma2[img_id].zdim,
						op.local_Minvsigma2[img_id].xdim - 1);

	#ifdef TIMING
				if (op.part_id == baseMLO->exp_my_first_part_id)
					baseMLO->timer.tic(baseMLO->TIMING_WSUM_BACKPROJ);
	#endif

				// If doing pseudo gold standard select random half-model
				int iproj_offset = 0;
				if (baseMLO->grad_pseudo_halfsets)
					// Backproject every other particle into separate volumes
					iproj_offset = (op.part_id % 2) * baseMLO->mymodel.nr_classes

				CTIC(accMLO->timer,"backproject");
				runBackProjectKernel(
					accMLO->bundle->backprojectors[iproj + iproj_offset],
					projKernel,
					&(~Fimgs)[re_nomask_offset], //~Fimgs_nomask_real,
					&(~Fimgs)[im_nomask_offset], //~Fimgs_nomask_imag,
					&(~trans_xyz)[trans_x_offset], //~trans_x,
					&(~trans_xyz)[trans_y_offset], //~trans_y,
					&(~trans_xyz)[trans_z_offset], //~trans_z,
					&(~sorted_weights)[classPos],
					~Minvsigma2s,
					~ctfs,
					translation_num,
					(XFLOAT) op.significant_weight[img_id],
					(XFLOAT) (baseMLO->is_som_iter ? class_sum_weight[iclass] : op.sum_weight[img_id]),
					~eulers[iclass],
					op.local_Minvsigma2[img_id].xdim,
					op.local_Minvsigma2[img_id].ydim,
					op.local_Minvsigma2[img_id].zdim,
					orientation_num,
					accMLO->dataIs3D,
					(baseMLO->do_grad),
					ctf_premultiplied,
					accMLO->classStreams[iclass]);

				CTOC(accMLO->timer,"backproject");

	#ifdef TIMING
				if (op.part_id == baseMLO->exp_my_first_part_id)
					baseMLO->timer.toc(baseMLO->TIMING_WSUM_BACKPROJ);
	#endif

				//Update indices
				AAXA_pos += image_size;
				classPos += orientation_num*translation_num;

			} // end loop iclass

			CUSTOM_ALLOCATOR_REGION_NAME("UNSET");

			// NOTE: We've never seen that this sync is necessary, but it is needed in principle, and
			// its absence in other parts of the code has caused issues. It is also very low-cost.
			for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
				DEBUG_HANDLE_ERROR(cudaStreamSynchronize(accMLO->classStreams[exp_iclass]));
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

			wdiff2s.cpToHost();
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(cudaStreamPerThread));

			AAXA_pos=0;

			for (unsigned long exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{
				if((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[img_id].class_entries[exp_iclass] == 0))
					continue;
				for (long int j = 0; j < image_size; j++)
				{
					int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine[optics_group], j);
					if (ires > -1 && baseMLO->do_scale_correction &&
							DIRECT_A1D_ELEM(baseMLO->mymodel.data_vs_prior_class[exp_iclass], ires) > 3.)
					{
						exp_wsum_scale_correction_AA[img_id] += wdiff2s[AA_offset+AAXA_pos+j];
						exp_wsum_scale_correction_XA[img_id] += wdiff2s[XA_offset+AAXA_pos+j];
					}
				}
				AAXA_pos += image_size;
			} // end loop iclass

			for (unsigned long j = 0; j < image_size; j++)
			{
				int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine[optics_group], j);
				if (ires > -1)
				{
					thr_wsum_sigma2_noise[img_id].data[ires] += (RFLOAT) wdiff2s[sum_offset+j];
					exp_wsum_norm_correction[img_id] += (RFLOAT) wdiff2s[sum_offset+j]; //TODO could be gpu-reduced
				}
			}

		} // end loop img_id
		CTOC(accMLO->timer,"maximization");

		CTIC(accMLO->timer,"store_post_gpu");

		// Extend norm_correction and sigma2_noise estimation to higher resolutions for all particles
		// Also calculate dLL for each particle and store in metadata
		// loop over all images inside this particle
		RFLOAT thr_avg_norm_correction = 0.;
		RFLOAT thr_sum_dLL = 0., thr_sum_Pmax = 0.;
		for (int img_id = 0; img_id < sp.nr_images; img_id++)
		{

			int my_metadata_offset = op.metadata_offset + img_id;
			int group_id = baseMLO->mydata.getGroupId(op.part_id, img_id);
			const int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);
			RFLOAT my_pixel_size = baseMLO->mydata.getOpticsPixelSize(optics_group);
			int my_image_size = baseMLO->mydata.getOpticsImageSize(optics_group);

			// If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
			for (unsigned long ires = baseMLO->image_current_size[optics_group]/2 + 1; ires < baseMLO->image_full_size[optics_group]/2 + 1; ires++)
			{
				DIRECT_A1D_ELEM(thr_wsum_sigma2_noise[img_id], ires) += DIRECT_A1D_ELEM(op.power_img[img_id], ires);
				// Also extend the weighted sum of the norm_correction
				exp_wsum_norm_correction[img_id] += DIRECT_A1D_ELEM(op.power_img[img_id], ires);
			}

			// Store norm_correction
			// Multiply by old value because the old norm_correction term was already applied to the image
			if (baseMLO->do_norm_correction && baseMLO->mymodel.nr_bodies == 1)
			{
				RFLOAT old_norm_correction = DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_NORM);
				old_norm_correction /= baseMLO->mymodel.avg_norm_correction;
				// The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
				// The variance of the total image (on which one normalizes) is twice this value!
				RFLOAT normcorr = old_norm_correction * sqrt(exp_wsum_norm_correction[img_id] * 2.);
				thr_avg_norm_correction += normcorr;

				// Now set the new norm_correction in the relevant position of exp_metadata
				DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_NORM) = normcorr;


				// Print warning for strange norm-correction values
				if (!((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc) && DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_NORM) > 10.)
				{
					std::cout << " WARNING: norm_correction= "<< DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_NORM)
							<< " for particle " << op.part_id << " in group " << group_id + 1
							<< "; Are your groups large enough? Or is the reference on the correct greyscale?" << std::endl;
				}

			}

			// Store weighted sums for scale_correction
			if (baseMLO->do_scale_correction)
			{
				// Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
				exp_wsum_scale_correction_XA[img_id] /= baseMLO->mymodel.scale_correction[group_id];
				exp_wsum_scale_correction_AA[img_id] /= baseMLO->mymodel.scale_correction[group_id] * baseMLO->mymodel.scale_correction[group_id];

				thr_wsum_signal_product_spectra[img_id] += exp_wsum_scale_correction_XA[img_id];
				thr_wsum_reference_power_spectra[img_id] += exp_wsum_scale_correction_AA[img_id];
			}

			// Calculate DLL for each particle
			RFLOAT logsigma2 = 0.;
			RFLOAT remap_image_sizes = (baseMLO->mymodel.ori_size * baseMLO->mymodel.pixel_size) / (my_image_size * my_pixel_size);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(baseMLO->Mresol_fine[optics_group])
			{
				int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine[optics_group], n);
				int ires_remapped = ROUND(remap_image_sizes * ires);
				// Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
				// Also exclude origin from logsigma2, as this will not be considered in the P-calculations
				if (ires > 0 && ires_remapped < XSIZE(baseMLO->mymodel.sigma2_noise[optics_group]))
					logsigma2 += log( 2. * PI * DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[optics_group], ires_remapped));
			}
			RFLOAT dLL;

			if ((baseMLO->iter==1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
				dLL = -op.min_diff2[img_id];
			else
				dLL = log(op.sum_weight[img_id]) - op.min_diff2[img_id] - logsigma2;

			// Store dLL of each image in the output array, and keep track of total sum
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_DLL) = dLL;
			thr_sum_dLL += dLL;

			// Also store sum of Pmax
			thr_sum_Pmax += DIRECT_A2D_ELEM(baseMLO->exp_metadata, my_metadata_offset, METADATA_PMAX);

		}

		// Now, inside a global_mutex, update the other weighted sums among all threads
		#pragma omp critical(AccMLO_global)
		{
			for (int img_id = 0; img_id < sp.nr_images; img_id++)
			{
				long int igroup = baseMLO->mydata.getGroupId(op.part_id, img_id);
				int optics_group = baseMLO->mydata.getOpticsGroup(op.part_id, img_id);


                if (baseMLO->mydata.obsModel.getCtfPremultiplied(optics_group))
                {
                    RFLOAT myscale = XMIPP_MAX(0.001, baseMLO->mymodel.scale_correction[igroup]);
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(baseMLO->Mresol_fine[optics_group])
                    {
                        int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine[optics_group], n);
                        if (ires > -1)
                            DIRECT_MULTIDIM_ELEM(thr_wsum_ctf2[img_id], ires) += myscale * DIRECT_MULTIDIM_ELEM(op.local_Fctf[img_id], n);
                    }
                }

                if (do_subtomo_correction)
                {
                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(baseMLO->Mresol_fine[optics_group])
                    {
                        int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine[optics_group], n);
                        if (ires > -1)
                            DIRECT_MULTIDIM_ELEM(thr_wsum_stMulti[img_id], ires) += DIRECT_MULTIDIM_ELEM(exp_local_STMulti[img_id], n);
                    }
                }

                int my_image_size = baseMLO->mydata.getOpticsImageSize(optics_group);
				RFLOAT my_pixel_size = baseMLO->mydata.getOpticsPixelSize(optics_group);
				RFLOAT remap_image_sizes = (baseMLO->mymodel.ori_size * baseMLO->mymodel.pixel_size) / (my_image_size * my_pixel_size);
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(thr_wsum_sigma2_noise[img_id])
				{
					int i_resam = ROUND(i * remap_image_sizes);
					if (i_resam < XSIZE(baseMLO->wsum_model.sigma2_noise[optics_group]))
					{
						DIRECT_A1D_ELEM(baseMLO->wsum_model.sigma2_noise[optics_group], i_resam) += DIRECT_A1D_ELEM(thr_wsum_sigma2_noise[img_id], i);
                        DIRECT_A1D_ELEM(baseMLO->wsum_model.sumw_ctf2[optics_group], i_resam) += DIRECT_A1D_ELEM(thr_wsum_ctf2[img_id], i);

                        if (do_subtomo_correction)
                           DIRECT_A1D_ELEM(baseMLO->wsum_model.sumw_stMulti[optics_group], i_resam) += DIRECT_A1D_ELEM(thr_wsum_stMulti[img_id], i);
                    }
				}
				baseMLO->wsum_model.sumw_group[optics_group] += thr_sumw_group[img_id];
				if (baseMLO->do_scale_correction)
				{
					baseMLO->wsum_model.wsum_signal_product[igroup] += thr_wsum_signal_product_spectra[img_id];
					baseMLO->wsum_model.wsum_reference_power[igroup] += thr_wsum_reference_power_spectra[img_id];
				}
			}
			for (int n = 0; n < baseMLO->mymodel.nr_classes; n++)
			{
				baseMLO->wsum_model.pdf_class[n] += thr_wsum_pdf_class[n];
				if (baseMLO->mymodel.ref_dim == 2)
				{
					XX(baseMLO->wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsetx_class[n];
					YY(baseMLO->wsum_model.prior_offset_class[n]) += thr_wsum_prior_offsety_class[n];
				}
			}

			for (int n = 0; n < baseMLO->mymodel.nr_classes * baseMLO->mymodel.nr_bodies; n++)
			{
				if (!(baseMLO->do_skip_align || baseMLO->do_skip_rotate) )
					baseMLO->wsum_model.pdf_direction[n] += thr_wsum_pdf_direction[n];
			}

			baseMLO->wsum_model.sigma2_offset += thr_wsum_sigma2_offset;

			if (baseMLO->do_norm_correction && baseMLO->mymodel.nr_bodies == 1)
				baseMLO->wsum_model.avg_norm_correction += thr_avg_norm_correction;

			baseMLO->wsum_model.LL += thr_sum_dLL;
			baseMLO->wsum_model.ave_Pmax += thr_sum_Pmax;
		}	
	} // end if !do_skip_maximization

	CTOC(accMLO->timer,"store_post_gpu");
#ifdef TIMING
	if (op.part_id == baseMLO->exp_my_first_part_id)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_WSUM);
#endif
}

// ----------------------------------------------------------------------------
// -------------------- accDoExpectationOneParticle ---------------------------
// ----------------------------------------------------------------------------
template <class MlClass>
void accDoExpectationOneParticle(MlClass *myInstance, unsigned long part_id_sorted, int thread_id, AccPtrFactory ptrFactory)
{
	SamplingParameters sp;
	MlOptimiser *baseMLO = myInstance->baseMLO;

	CTIC(timer,"oneParticle");
#ifdef TIMING
	// Only time one thread
	if (thread_id == 0)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_A);
#endif

	long int part_id = baseMLO->mydata.sorted_idx[part_id_sorted];
	sp.nr_images = baseMLO->mydata.numberOfImagesInParticle(part_id);

	OptimisationParamters op(sp.nr_images, part_id);
	if (baseMLO->mydata.is_3D)
		op.FstMulti.resize(sp.nr_images);

	// In the first iteration, multiple seeds will be generated
	// A single random class is selected for each pool of images, and one does not marginalise over the orientations
	// The optimal orientation is based on signal-product (rather than the signal-intensity sensitive Gaussian)
	// If do_firstiter_cc, then first perform a single iteration with K=1 and cross-correlation criteria, afterwards

	// Decide which classes to integrate over (for random class assignment in 1st iteration)
	sp.iclass_min = 0;
	sp.iclass_max = baseMLO->mymodel.nr_classes - 1;
	// low-pass filter again and generate the seeds
	if (baseMLO->do_generate_seeds)
	{
		if (baseMLO->do_firstiter_cc && baseMLO->iter == 1)
		{
			// In first (CC) iter, use a single reference (and CC)
			sp.iclass_min = sp.iclass_max = 0;
		}
		else if ( (baseMLO->do_firstiter_cc && baseMLO->iter == 2) ||

				(!baseMLO->do_firstiter_cc && baseMLO->iter == 1))
		{
			// In second CC iter, or first iter without CC: generate the seeds
			// Now select a single random class
			// exp_part_id is already in randomized order (controlled by -seed)
			// WARNING: USING SAME iclass_min AND iclass_max FOR SomeParticles!!
			// Make sure random division is always the same with the same seed
			long int idx = part_id_sorted - baseMLO->exp_my_first_part_id;
			if (idx >= baseMLO->exp_random_class_some_particles.size())
				REPORT_ERROR("BUG: expectationOneParticle idx>random_class_some_particles.size()");
			sp.iclass_min = sp.iclass_max = baseMLO->exp_random_class_some_particles[idx];
		}
	}
    // Loop over all bodies of the multi-body refinement
    // Basically, subsequently align and store weighted sums for each body
    for (int ibody = 0; ibody < baseMLO->mymodel.nr_bodies; ibody++)
    {

    	OptimisationParamters op(sp.nr_images, part_id);
		if (baseMLO->mydata.is_3D)
			op.FstMulti.resize(sp.nr_images);

		// Skip this body if keep_fixed_bodies[ibody] or if it's angular accuracy is worse than 1.5x the sampling rate
    	if ( baseMLO->mymodel.nr_bodies > 1 && baseMLO->mymodel.keep_fixed_bodies[ibody] > 0)
			continue;


		// Global exp_metadata array has metadata of all particles. Where does part_id start?
		for (long int iori = baseMLO->exp_my_first_part_id; iori <= baseMLO->exp_my_last_part_id; iori++)
		{
			if (iori == part_id_sorted) break;
			op.metadata_offset += baseMLO->mydata.numberOfImagesInParticle(iori);
		}
#ifdef TIMING
// Only time one thread
if (thread_id == 0)
baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2_A);
#endif
		CTIC(timer,"getFourierTransformsAndCtfs");
		getFourierTransformsAndCtfs<MlClass>(part_id, op, sp, baseMLO, myInstance, ptrFactory, ibody);
		CTOC(timer,"getFourierTransformsAndCtfs");

		// To deal with skipped alignments/rotations
		if (baseMLO->do_skip_align)
		{
			sp.itrans_min = sp.itrans_max = sp.idir_min = sp.idir_max = sp.ipsi_min = sp.ipsi_max =
					part_id_sorted - baseMLO->exp_my_first_part_id;
		}
		else
		{
			sp.itrans_min = 0;
			sp.itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;
		}
		if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
		{
			sp.idir_min = sp.idir_max = sp.ipsi_min = sp.ipsi_max =
					part_id_sorted - baseMLO->exp_my_first_part_id;
		}
		else if (baseMLO->do_only_sample_tilt)
		{
			sp.idir_min = 0;
			sp.idir_max = baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior) - 1;
			sp.ipsi_min = sp.ipsi_max = part_id_sorted - baseMLO->exp_my_first_part_id;

		}
		else
		{
			sp.idir_min = sp.ipsi_min = 0;
			sp.idir_max = baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior) - 1;
			sp.ipsi_max = baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior ) - 1;
		}

		// Initialise significant weight to minus one, so that all coarse sampling points will be handled in the first pass
		op.significant_weight.resize(sp.nr_images, -1.);

		// Only perform a second pass when using adaptive oversampling
		//int nr_sampling_passes = (baseMLO->adaptive_oversampling > 0) ? 2 : 1;
		// But on the gpu the data-structures are different between passes, so we need to make a symbolic pass to set the weights up for storeWS
		int nr_sampling_passes = 2;

		/// -- This is a iframe-indexed vector, each entry of which is a dense data-array. These are replacements to using
		//    Mweight in the sparse (Fine-sampled) pass, coarse is unused but created empty input for convert ( FIXME )
		std::vector <IndexedDataArray > CoarsePassWeights(1, ptrFactory);
		std::vector <IndexedDataArray > FinePassWeights(sp.nr_images, ptrFactory);

		// -- This is a iframe-indexed vector, each entry of which is a class-indexed vector of masks, one for each
		//    class in FinePassWeights
		std::vector < std::vector <IndexedDataArrayMask > > FinePassClassMasks(sp.nr_images, std::vector <IndexedDataArrayMask >(baseMLO->mymodel.nr_classes, ptrFactory));

		// -- This is a iframe-indexed vector, each entry of which is parameters used in the projection-operations *after* the
		//    coarse pass, declared here to keep scope to storeWS
		std::vector < ProjectionParams > FineProjectionData(sp.nr_images, baseMLO->mymodel.nr_classes);

		std::vector < AccPtrBundle > bundleD2(sp.nr_images, ptrFactory.makeBundle());
		std::vector < AccPtrBundle > bundleSWS(sp.nr_images, ptrFactory.makeBundle());

		for (int ipass = 0; ipass < nr_sampling_passes; ipass++)
		{
			CTIC(timer,"weightPass");
#ifdef TIMING
// Only time one thread
if (thread_id == 0)
baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_B);
#endif

			// Use coarse sampling in the first pass, oversampled one the second pass
			sp.current_oversampling = (ipass == 0) ? 0 : baseMLO->adaptive_oversampling;

			sp.nr_dir = (baseMLO->do_skip_align || baseMLO->do_skip_rotate) ? 1 : baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior);
			sp.nr_psi = (baseMLO->do_skip_align || baseMLO->do_skip_rotate) ? 1 : baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior);
			sp.nr_trans = (baseMLO->do_skip_align) ? 1 : baseMLO->sampling.NrTranslationalSamplings();
			sp.nr_oversampled_rot = baseMLO->sampling.oversamplingFactorOrientations(sp.current_oversampling);
			sp.nr_oversampled_trans = baseMLO->sampling.oversamplingFactorTranslations(sp.current_oversampling);
#ifdef TIMING
// Only time one thread
if (thread_id == 0)
baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2_B);
#endif

			op.min_diff2.resize(sp.nr_images, 0);

			if (ipass == 0)
			{
				unsigned long weightsPerPart(baseMLO->mymodel.nr_classes * sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.nr_oversampled_rot * sp.nr_oversampled_trans);

				op.Mweight.resizeNoCp(1,1,sp.nr_images, weightsPerPart);

				AccPtr<XFLOAT> Mweight = ptrFactory.make<XFLOAT>();

				Mweight.setSize(sp.nr_images * weightsPerPart);
				Mweight.setHostPtr(op.Mweight.data);
				Mweight.deviceAlloc();
				deviceInitValue<XFLOAT>(Mweight, -std::numeric_limits<XFLOAT>::max());
				Mweight.streamSync();

				CTIC(timer,"getAllSquaredDifferencesCoarse");
				getAllSquaredDifferencesCoarse<MlClass>(ipass, op, sp, baseMLO, myInstance, Mweight, ptrFactory, ibody);
				CTOC(timer,"getAllSquaredDifferencesCoarse");

				CTIC(timer,"convertAllSquaredDifferencesToWeightsCoarse");
				convertAllSquaredDifferencesToWeights<MlClass>(ipass, op, sp, baseMLO, myInstance, CoarsePassWeights, FinePassClassMasks, Mweight, ptrFactory, ibody);
				CTOC(timer,"convertAllSquaredDifferencesToWeightsCoarse");
			}
			else
			{
#ifdef TIMING
// Only time one thread
if (thread_id == 0)
baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_D);
#endif
//					// -- go through all classes and generate projectionsetups for all classes - to be used in getASDF and storeWS below --
//					// the reason to do this globally is subtle - we want the orientation_num of all classes to estimate a largest possible
//					// weight-array, which would be insanely much larger than necessary if we had to assume the worst.
				for (int img_id = 0; img_id < sp.nr_images; img_id++)
				{
					FineProjectionData[img_id].orientationNumAllClasses = 0;
					for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
					{
						if(exp_iclass>0)
							FineProjectionData[img_id].class_idx[exp_iclass] = FineProjectionData[img_id].rots.size();
						FineProjectionData[img_id].class_entries[exp_iclass] = 0;

						CTIC(timer,"generateProjectionSetup");
						FineProjectionData[img_id].orientationNumAllClasses += generateProjectionSetupFine(
								op,
								sp,
								baseMLO,
								exp_iclass,
								FineProjectionData[img_id]);
						CTOC(timer,"generateProjectionSetup");

					}
					//set a maximum possible size for all weights (to be reduced by significance-checks)
					size_t dataSize = FineProjectionData[img_id].orientationNumAllClasses*sp.nr_trans*sp.nr_oversampled_trans;
					FinePassWeights[img_id].setDataSize(dataSize);
					FinePassWeights[img_id].dual_alloc_all();

					bundleD2[img_id].setSize(2*(FineProjectionData[img_id].orientationNumAllClasses*sp.nr_trans*sp.nr_oversampled_trans)*sizeof(unsigned long));
					bundleD2[img_id].allAlloc();
				}
#ifdef TIMING
// Only time one thread
if (thread_id == 0)
baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2_D);
#endif

				CTIC(timer,"getAllSquaredDifferencesFine");
				getAllSquaredDifferencesFine<MlClass>(ipass, op, sp, baseMLO, myInstance, FinePassWeights, FinePassClassMasks, FineProjectionData, ptrFactory, ibody, bundleD2);
				CTOC(timer,"getAllSquaredDifferencesFine");
				FinePassWeights[0].weights.cpToHost();

				AccPtr<XFLOAT> Mweight = ptrFactory.make<XFLOAT>(); //DUMMY

				CTIC(timer,"convertAllSquaredDifferencesToWeightsFine");
				convertAllSquaredDifferencesToWeights<MlClass>(ipass, op, sp, baseMLO, myInstance, FinePassWeights, FinePassClassMasks, Mweight, ptrFactory, ibody);
				CTOC(timer,"convertAllSquaredDifferencesToWeightsFine");

			}

			CTOC(timer,"weightPass");
		}
#ifdef TIMING
// Only time one thread
if (thread_id == 0)
baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_E);
#endif

		// For the reconstruction step use mymodel.current_size!
		// as of 3.1, no longer necessary?
		sp.current_image_size = baseMLO->mymodel.current_size;

	for (unsigned long img_id = 0; img_id < sp.nr_images; img_id++)
	{
		bundleSWS[img_id].setSize(2*(FineProjectionData[img_id].orientationNumAllClasses)*sizeof(unsigned long));
		bundleSWS[img_id].allAlloc();
	}

#ifdef TIMING
// Only time one thread
if (thread_id == 0)
baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2_E);
#endif
		CTIC(timer,"storeWeightedSums");
		storeWeightedSums<MlClass>(op, sp, baseMLO, myInstance, FinePassWeights, FineProjectionData, FinePassClassMasks, ptrFactory, ibody, bundleSWS);
		CTOC(timer,"storeWeightedSums");

		for (long int img_id = 0; img_id < sp.nr_images; img_id++)
		{
			FinePassWeights[img_id].dual_free_all();
		}
    }

	CTOC(timer,"oneParticle");
}

#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_benchmark_utils.cuh"
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include "src/gpu_utils/cuda_kernels/diff2.cuh"
#include "src/gpu_utils/cuda_kernels/wavg.cuh"
#include "src/gpu_utils/cuda_utils.cuh"
#include "src/gpu_utils/cuda_helper_functions.cu"
#include "src/complex.h"
#include <fstream>
#include <cuda_runtime.h>
#include "src/parallel.h"
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <signal.h>

void getAllSquaredDifferencesCoarse(unsigned exp_ipass, OptimisationParamters &op, SamplingParameters &sp, MlOptimiser *baseMLO)
{

	CUDA_CPU_TIC("diff_pre_gpu");

	op.Mweight.resize(sp.nr_particles, baseMLO->mymodel.nr_classes * sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.nr_oversampled_rot * sp.nr_oversampled_trans);
	op.Mweight.initConstant(-999.);

	op.Mcoarse_significant.resize(sp.nr_dir*sp.nr_psi*sp.nr_trans, 1);
	op.Mcoarse_significant.clear();

	op.min_diff2.clear();
	op.min_diff2.resize(sp.nr_particles, 99.e99);

	std::vector<MultidimArray<Complex > > dummy;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, dummy, op.Fctfs, op.local_Fimgs_shifted, dummy,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);

	unsigned image_size = op.local_Minvsigma2s[0].nzyxdim;

	CUDA_CPU_TOC("diff_pre_gpu");

	// Loop only from sp.iclass_min to sp.iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
	{
		if (baseMLO->mymodel.pdf_class[exp_iclass] > 0.)
		{
			// Local variables
			CudaGlobalPtr<FLOAT> gpuMinvsigma2(image_size);
			gpuMinvsigma2.device_alloc();

			// Mapping index look-up table
			std::vector< long unsigned > iorientclasses, iover_rots;
			std::vector< double > rots, tilts, psis;

			CUDA_CPU_TIC("generateProjectionSetupCoarse");
			long unsigned orientation_num = generateProjectionSetup(
					op,
					sp,
					baseMLO,
					exp_ipass == 0, //coarse
					exp_iclass,
					rots, tilts, psis,
					iorientclasses,
					iover_rots);

			CUDA_CPU_TOC("generateProjectionSetupCoarse");


			CUDA_CPU_TIC("generateEulerMatricesCoarse");
			CudaGlobalPtr<FLOAT> eulers(9 * orientation_num);

			generateEulerMatrices(
					baseMLO->mymodel.PPref[exp_iclass].padding_factor,
					rots,
					tilts,
					psis,
					eulers,
					!IS_NOT_INV);

			CUDA_CPU_TOC("generateEulerMatricesCoarse");

			CUDA_GPU_TIC("eulersMemCpCoarse");
			eulers.device_alloc();
			eulers.cp_to_device();
			CUDA_GPU_TAC("eulersMemCpCoarse");

			/*=======================================================================================
			                                  	  Particle Iteration
			=========================================================================================*/

			for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
			{
				long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
				long int group_id = baseMLO->mydata.getGroupId(part_id);

				/*====================================
				        Generate Translations
				======================================*/

				CUDA_CPU_TIC("translation_1");

				CudaGlobalPtr<FLOAT> Fimgs_real(image_size * sp.nr_trans * sp.nr_oversampled_trans);
				CudaGlobalPtr<FLOAT> Fimgs_imag(image_size * sp.nr_trans * sp.nr_oversampled_trans);

				std::vector< double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
				long unsigned translation_num(0), ihidden(0);
				std::vector< long unsigned > iover_transes, itranses, ihiddens;

				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++, ihidden++)
				{
					baseMLO->sampling.getTranslations(itrans, sp.current_oversampling,
							oversampled_translations_x, oversampled_translations_y, oversampled_translations_z );

					for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++)
					{
						/// Now get the shifted image
						// Use a pointer to avoid copying the entire array again in this highly expensive loop
						Complex *myAB;
						if (sp.current_oversampling == 0)
						{
							myAB = (op.local_Minvsigma2s[0].ydim == baseMLO->coarse_size) ? baseMLO->global_fftshifts_ab_coarse[itrans].data
									: baseMLO->global_fftshifts_ab_current[itrans].data;
						}
						else
						{
							int iitrans = itrans * sp.nr_oversampled_trans +  iover_trans;
							myAB = (baseMLO->strict_highres_exp > 0.) ? baseMLO->global_fftshifts_ab2_coarse[iitrans].data
									: baseMLO->global_fftshifts_ab2_current[iitrans].data;
						}


						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
						{
							FLOAT real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real
									- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag;
							FLOAT imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag
									+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real;

							//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
							if (baseMLO->do_scale_correction)
							{
								FLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
								real /= myscale;
								imag /= myscale;
							}
							if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
							{
								real /= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
								imag /= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
							}
							Fimgs_real[translation_num * image_size + n] = real;
							Fimgs_imag[translation_num * image_size + n] = imag;
						}

						translation_num ++;

						ihiddens.push_back(ihidden);
						itranses.push_back(itrans);
						iover_transes.push_back(iover_trans);
					}
				}

				Fimgs_real.size = translation_num * image_size;
				Fimgs_imag.size = translation_num * image_size;

				Fimgs_real.device_alloc();
				Fimgs_imag.device_alloc();

				Fimgs_real.cp_to_device();
				Fimgs_imag.cp_to_device();

				CUDA_CPU_TOC("translation_1");

				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
				{
					gpuMinvsigma2[n] = op.local_Minvsigma2s[ipart].data[n];
				}

				if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
					{
						gpuMinvsigma2[n] *= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n)*DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
					}
				}
				if (baseMLO->do_scale_correction)
				{
					FLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
					{
						gpuMinvsigma2[n] *= myscale * myscale;
					}
				}

				gpuMinvsigma2.cp_to_device();

				CudaGlobalPtr<FLOAT> diff2s(orientation_num*translation_num);
				diff2s.device_alloc();

				/*====================================
				    	   Kernel Call
				======================================*/

				Cuda3DProjectorKernel projKernel = Cuda3DProjectorKernel::makeKernel(
						baseMLO->cudaProjectors[exp_iclass],
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						op.local_Minvsigma2s[0].xdim-1);

				runDiff2KernelCoarse(
						projKernel,
						gpuMinvsigma2,
						Fimgs_real,
						Fimgs_imag,
						eulers,
						diff2s,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						image_size,
						ipart,
						group_id,
						exp_iclass);

				/*====================================
				    	   Retrieve Results
				======================================*/

				op.min_diff2[ipart] = thrustGetMinVal(diff2s);

				CUDA_GPU_TIC("diff2sMemCpCoarse");
				diff2s.cp_to_host();
				CUDA_GPU_TAC("diff2sMemCpCoarse");

				HANDLE_ERROR(cudaDeviceSynchronize());
				CUDA_GPU_TOC();

				for (unsigned i = 0; i < orientation_num; i ++)
				{
					unsigned iorientclass = iorientclasses[i];
					for (unsigned j = 0; j < translation_num; j ++)
						DIRECT_A2D_ELEM(op.Mweight, ipart, iorientclass * translation_num + j) = diff2s[i * translation_num + j];
				}
			} // end loop ipart
		} // end if class significant
	} // end loop iclass
}

void getAllSquaredDifferencesFine(unsigned exp_ipass,
		 	 	 	 	 	 	  OptimisationParamters &op,
		 	 	 	 	 	 	  SamplingParameters &sp,
		 	 	 	 	 	 	  MlOptimiser *baseMLO,
		 	 	 	 	 	 	  CudaGlobalPtr <long unsigned> &rot_id,
		 	 	 	 	 	 	  CudaGlobalPtr <long unsigned> &rot_idx,
		 	 	 	 	 	 	  CudaGlobalPtr <long unsigned> &trans_idx,
		 	 	 	 	 	 	  CudaGlobalPtr <long unsigned> &ihidden_overs,
		 	 	 	 	 	 	  CudaGlobalPtr <long unsigned> &job_idx,
		 	 	 	 	 	 	  CudaGlobalPtr <long unsigned> &job_num,
		 	 	 	 	 	 	  CudaGlobalPtr <FLOAT> &weights)
{

	CUDA_CPU_TIC("diff_pre_gpu");

	//for scale_correction
	int group_id;

	//printf("sp.nr_oversampled_rot=%d\n", (unsigned)sp.nr_oversampled_rot);

	op.Mweight.resize(sp.nr_particles, baseMLO->mymodel.nr_classes * sp.nr_dir * sp.nr_psi * sp.nr_trans * sp.nr_oversampled_rot * sp.nr_oversampled_trans);
	op.Mweight.initConstant(-999.);

	op.min_diff2.clear();
	op.min_diff2.resize(sp.nr_particles, 99.e99);

	std::vector<MultidimArray<Complex > > dummy;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, dummy, op.Fctfs, op.local_Fimgs_shifted, dummy,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);

	MultidimArray<Complex > Fref;
	Fref.resize(op.local_Minvsigma2s[0]);

	unsigned image_size = op.local_Minvsigma2s[0].nzyxdim;

	CUDA_CPU_TOC("diff_pre_gpu");

	// Loop only from sp.iclass_min to sp.iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
	{
		if (baseMLO->mymodel.pdf_class[exp_iclass] > 0.)
		{
			// Local variables
			std::vector< double > oversampled_rot, oversampled_tilt, oversampled_psi;
			std::vector< double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
			CudaGlobalPtr<FLOAT> gpuMinvsigma2(image_size);
			gpuMinvsigma2.device_alloc();

			// Mapping index look-up table
			std::vector< long unsigned > iorientclasses, iover_rots;
			std::vector< double > rots, tilts, psis;

			CUDA_CPU_TIC("generateProjectionSetup");
			long unsigned orientation_num = generateProjectionSetup(
					op,
					sp,
					baseMLO,
					exp_ipass == 0, //coarse
					exp_iclass,
					rots, tilts, psis,
					iorientclasses,
					iover_rots);

			CUDA_CPU_TOC("generateProjectionSetup");


			CUDA_CPU_TIC("generateEulerMatrices");
			CudaGlobalPtr<FLOAT> eulers(9 * orientation_num);

			generateEulerMatrices(
					baseMLO->mymodel.PPref[exp_iclass].padding_factor,
					rots,
					tilts,
					psis,
					eulers,
					!IS_NOT_INV);

			CUDA_CPU_TOC("generateEulerMatrices");

			CUDA_GPU_TIC("eulersMemCp_1");
			eulers.device_alloc();
			eulers.cp_to_device();
			CUDA_GPU_TAC("eulersMemCp_1");

			/*=======================================================================================
			                                  	  Particle Iteration
			=========================================================================================*/

			for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
			{
				/*====================================
				        Generate Translations
				======================================*/

				CUDA_CPU_TIC("translation_1");

				CudaGlobalPtr<FLOAT> Fimgs_real(image_size * sp.nr_trans * sp.nr_oversampled_trans);
				CudaGlobalPtr<FLOAT> Fimgs_imag(image_size * sp.nr_trans * sp.nr_oversampled_trans);

				long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
				long unsigned translation_num(0), ihidden(0);
				std::vector< long unsigned > iover_transes, itranses, ihiddens;

				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++, ihidden++)
				{
					baseMLO->sampling.getTranslations(itrans, sp.current_oversampling,
							oversampled_translations_x, oversampled_translations_y, oversampled_translations_z );

					for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++)
					{
						/// Now get the shifted image
						// Use a pointer to avoid copying the entire array again in this highly expensive loop
						Complex *myAB;
						if (sp.current_oversampling == 0)
						{
							myAB = (Fref.ydim == baseMLO->coarse_size) ? baseMLO->global_fftshifts_ab_coarse[itrans].data
									: baseMLO->global_fftshifts_ab_current[itrans].data;
						}
						else
						{
							int iitrans = itrans * sp.nr_oversampled_trans +  iover_trans;
							myAB = (baseMLO->strict_highres_exp > 0.) ? baseMLO->global_fftshifts_ab2_coarse[iitrans].data
									: baseMLO->global_fftshifts_ab2_current[iitrans].data;
						}


						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
						{
							FLOAT real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real
									- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag;
							FLOAT imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag
									+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real;

							//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
							if (baseMLO->do_scale_correction)
							{
								//group_id = mydata.getGroupId(part_id);
								FLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
								real /= myscale;
								imag /= myscale;
							}
							if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
							{
								real /= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
								imag /= DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n);
							}
							Fimgs_real[translation_num * image_size + n] = real;
							Fimgs_imag[translation_num * image_size + n] = imag;
						}
						translation_num ++;

						ihiddens.push_back(ihidden);
						itranses.push_back(itrans);
						iover_transes.push_back(iover_trans);
					}
				}
				CUDA_CPU_TOC("translation_1");

				/*===========================================
				   Determine significant comparison indices
				=============================================*/
				//      This section is annoying to test because
				//		it can't complete on first pass, since
				//		the significance has never been set

				CUDA_CPU_TIC("pair_list_1");
				long unsigned coarse_num = sp.nr_dir*sp.nr_psi*sp.nr_trans, significant_num(0), k=0;
				long int nr_over_orient = baseMLO->sampling.oversamplingFactorOrientations(sp.current_oversampling);
				long int nr_over_trans = baseMLO->sampling.oversamplingFactorTranslations(sp.current_oversampling);

				significant_num = divideOrientationsIntoBlockjobs(	op,	sp,															   // alot of different type inputs...
																	orientation_num, translation_num,
																	iorientclasses,	iover_rots,	iover_transes, ihiddens,
																	nr_over_orient, nr_over_trans, ipart,
																	rot_id, rot_idx, trans_idx, ihidden_overs, job_idx,job_num);               // ..and output into index-arrays
				CUDA_CPU_TOC("pair_list_1");

				weights.size=significant_num; // use the weights array from the above level
				weights.host_alloc();
				weights.device_alloc();

				Cuda3DProjectorKernel projKernel = Cuda3DProjectorKernel::makeKernel(
						baseMLO->cudaProjectors[exp_iclass],
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						op.local_Minvsigma2s[0].xdim-1);

				runDiff2KernelFine(
						projKernel,
						gpuMinvsigma2,
						Fimgs_real,
						Fimgs_imag,
						eulers,
						rot_id,
						rot_idx,
						trans_idx,
						job_idx,
						job_num,
						weights,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						significant_num,
						image_size,
						ipart,
						group_id,
						exp_iclass
						);
				eulers.free_device();

				HANDLE_ERROR(cudaDeviceSynchronize());
				CUDA_GPU_TOC();

				/*====================================
				    	Write To Destination
				======================================*/

				CUDA_CPU_TIC("collect_data_1");
				op.min_diff2[ipart] = thrustGetMinVal(weights);
				if (exp_ipass == 0)
				{
					CUDA_GPU_TIC("diff2sMemCp");
					weights.cp_to_host();
					CUDA_GPU_TAC("diff2sMemCp");
					for (long unsigned k = 0; k<significant_num; k++)
					DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_overs[k]) = weights[k];
				}
				CUDA_CPU_TOC("collect_data_1");
			} // end loop ipart
		} // end if class significant
	} // end loop iclass
}


void convertAllSquaredDifferencesToWeights(unsigned exp_ipass,
											OptimisationParamters &op,
											SamplingParameters &sp,
											MlOptimiser *baseMLO,
											CudaGlobalPtr <long unsigned> &rot_id,
											CudaGlobalPtr <long unsigned> &rot_idx,
											CudaGlobalPtr <long unsigned> &trans_idx,
											CudaGlobalPtr <long unsigned> &ihidden_overs,
											CudaGlobalPtr <long unsigned> &job_idx,
											CudaGlobalPtr <long unsigned> &job_num,
											CudaGlobalPtr <FLOAT> 		  &weights)
{
	op.sum_weight.clear();
	op.sum_weight.resize(sp.nr_particles, 0.);

	// loop over all particles inside this ori_particle
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		double exp_thisparticle_sumweight = 0.;

		double old_offset_z;
		double old_offset_x = XX(op.old_offset[ipart]);
		double old_offset_y = YY(op.old_offset[ipart]);
		if (baseMLO->mymodel.data_dim == 3)
			old_offset_z = ZZ(op.old_offset[ipart]);

		if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
		{
			// Binarize the squared differences array to skip marginalisation
			double mymindiff2 = 99.e10;
			long int myminidx = -1;
			// Find the smallest element in this row of op.Mweight
			for (long int i = 0; i < XSIZE(op.Mweight); i++)
			{

				double cc = DIRECT_A2D_ELEM(op.Mweight, ipart, i);
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
			for (long int i = 0; i < XSIZE(op.Mweight); i++)
				DIRECT_A2D_ELEM(op.Mweight, ipart, i)= 0.;

			DIRECT_A2D_ELEM(op.Mweight, ipart, myminidx)= 1.;
			exp_thisparticle_sumweight += 1.;

		}
		else
		{
			// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{
				// TODO Move to RANK LEVEL
				/*=========================================
						Fetch+generate Orientation data
				===========================================*/
				CudaGlobalPtr<FLOAT >  pdf_orientation(sp.nr_dir * sp.nr_psi);
				pdf_orientation.size = sp.nr_dir * sp.nr_psi;
				for (long int idir = sp.idir_min, iorient = 0; idir <= sp.idir_max; idir++)
				{
					for (long int ipsi = sp.ipsi_min; ipsi <= sp.ipsi_max; ipsi++, iorient++)
					{
						//std::cerr << "orient "  << idir << "," << iorient <<  std::endl;
						// Get prior for this direction
						if (baseMLO->do_skip_align || baseMLO->do_skip_rotate)
						{
							pdf_orientation[iorient] = baseMLO->mymodel.pdf_class[exp_iclass];
						}
						else if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
						{
							pdf_orientation[iorient] = DIRECT_MULTIDIM_ELEM(baseMLO->mymodel.pdf_direction[exp_iclass], idir);
						}
						else
						{
							// P(orientation) = P(idir|dir_prior) * P(ipsi|psi_prior)
							// This is the probability of the orientation, given the gathered
							// statistics of all assigned orientations of the dataset, since we
							// are assigning a gaussian prior to all parameters.
							pdf_orientation[iorient] = op.directions_prior[idir] * op.psi_prior[ipsi];
						}
					}
				}

				// TODO Move to EXPECTATION-LEVEL ( potentially depenends on priors from getFourierTransformsAndCtfs() )
				// TODO Also hide under GetAllSq... at later stage (if neccessary)
				/*=========================================
						Fetch+generate Translation data
				===========================================*/
				CudaGlobalPtr<FLOAT >  pdf_offset(sp.nr_trans);
				double myprior_x, myprior_y, myprior_z;
				if (baseMLO->mymodel.ref_dim == 2)
				{
					myprior_x = XX(baseMLO->mymodel.prior_offset_class[exp_iclass]);
					myprior_y = YY(baseMLO->mymodel.prior_offset_class[exp_iclass]);
				}
				else
				{
					myprior_x = XX(op.prior[ipart]);
					myprior_y = YY(op.prior[ipart]);
					if (baseMLO->mymodel.data_dim == 3)
						myprior_z = ZZ(op.prior[ipart]);
				}

				for (long int itrans = sp.itrans_min; itrans <= sp.itrans_max; itrans++)
				{
					//std::cerr << "trans " << itrans << "," << jtrans <<  std::endl;
			        // To speed things up, only calculate pdf_offset at the coarse sampling.
					// That should not matter much, and that way one does not need to calculate all the OversampledTranslations
					double offset_x = old_offset_x + baseMLO->sampling.translations_x[itrans];
					double offset_y = old_offset_y + baseMLO->sampling.translations_y[itrans];
					double tdiff2 = (offset_x - myprior_x) * (offset_x - myprior_x) + (offset_y - myprior_y) * (offset_y - myprior_y);
					if (baseMLO->mymodel.data_dim == 3)
					{
						double offset_z = old_offset_z + baseMLO->sampling.translations_z[itrans];
						tdiff2 += (offset_z - myprior_z) * (offset_z - myprior_z);
					}
					// P(offset|sigma2_offset)
					// This is the probability of the offset, given the model offset and variance.
					if (baseMLO->mymodel.sigma2_offset < 0.0001)
						pdf_offset[itrans] = ( tdiff2 > 0.) ? 0. : 1.;
					else
						pdf_offset[itrans] = exp ( tdiff2 / (-2. * baseMLO->mymodel.sigma2_offset) ) / ( 2. * PI * baseMLO->mymodel.sigma2_offset );
				}

// TODO : Put back when  convertAllSquaredDifferencesToWeights is GPU-parallel.
//				// TMP DEBUGGING
//				if (baseMLO->mymodel.orientational_prior_mode != NOPRIOR && (pdf_offset==0. || pdf_orientation==0.))
//				{
//					pthread_mutex_lock(&global_mutex);
//					std::cerr << " pdf_offset= " << pdf_offset << " pdf_orientation= " << pdf_orientation << std::endl;
//					std::cerr << " ipart= " << ipart << " part_id= " << part_id << std::endl;
//					std::cerr << " iorient= " << iorient << " idir= " << idir << " ipsi= " << ipsi << std::endl;
//					//std::cerr << " sp.nr_psi= " << sp.nr_psi << " exp_nr_dir= " << exp_nr_dir << " sp.nr_trans= " << sp.nr_trans << std::endl;
//					for (long int i = 0; i < op.directions_prior.size(); i++)
//						std::cerr << " op.directions_prior["<<i<<"]= " << op.directions_prior[i] << std::endl;
//					for (long int i = 0; i < op.psi_prior.size(); i++)
//						std::cerr << " op.psi_prior["<<i<<"]= " << op.psi_prior[i] << std::endl;
//					REPORT_ERROR("ERROR! pdf_offset==0.|| pdf_orientation==0.");
//					//pthread_mutex_unlock(&global_mutex);
//				}
//				if (sp.nr_oversampled_rot == 0)
//					REPORT_ERROR("sp.nr_oversampled_rot == 0");
//				if (sp.nr_oversampled_trans == 0)
//					REPORT_ERROR("sp.nr_oversampled_trans == 0");

				// Now first loop over iover_rot, because that is the order in op.Mweight as well
//				long int ihidden_over = ihidden * sp.nr_oversampled_rot * sp.nr_oversampled_trans;

				/*================================================
					 Sumweights - exponentiation and reduction
				==================================================*/
				int oversamples = sp.nr_oversampled_trans * sp.nr_oversampled_rot;

				CUDA_CPU_TIC("sumweight1");
				CUDA_GPU_TIC("sumweightMemCp1");

				std::cerr << "summing weights on GPU... " << std::endl;
				pdf_orientation.put_on_device();
				pdf_offset.put_on_device();

				CUDA_GPU_TAC("sumweightMemCp1");

				CudaGlobalPtr<FLOAT > thisparticle_sumweight;
				long int block_num;

				if(exp_ipass==0)  //use Mweight for now
				{

					CudaGlobalPtr<FLOAT >  Mweight( &(op.Mweight.data[(ipart)*(op.Mweight).xdim]),
												sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples);
					Mweight.device_alloc();
					Mweight.cp_to_device();
					block_num = sp.nr_dir*sp.nr_psi*sp.nr_oversampled_rot/SUM_BLOCK_SIZE;

					dim3 block_dim(block_num);
					thisparticle_sumweight.size=block_num;
					thisparticle_sumweight.host_alloc();
					thisparticle_sumweight.device_alloc();
					CUDA_GPU_TIC("cuda_kernel_sumweight");
					cuda_kernel_sumweightCoarse<<<block_dim,SUM_BLOCK_SIZE>>>(	~pdf_orientation,
																			    ~pdf_offset,
																			    ~Mweight,
																			    ~thisparticle_sumweight,
																			    (FLOAT)op.min_diff2[ipart],
																			    sp.nr_oversampled_rot,
																			    sp.nr_oversampled_trans,
																			    sp.nr_trans );
					CUDA_GPU_TAC("cuda_kernel_sumweight");
					CUDA_GPU_TIC("sumweightMemCp2");
					Mweight.cp_to_host();  //FIXME remove when mapping is eliminated
					Mweight.free_device();
					CUDA_GPU_TAC("sumweightMemCp2");
				}
				else
				{
					block_num = rot_idx.size / SUM_BLOCK_SIZE;
					dim3 block_dim(block_num);
					thisparticle_sumweight.size=block_num;
					thisparticle_sumweight.host_alloc();
					thisparticle_sumweight.device_alloc();

					CUDA_GPU_TIC("cuda_kernel_sumweight");
					cuda_kernel_sumweightFine<<<block_dim,SUM_BLOCK_SIZE>>>(	~pdf_orientation,
																			    ~pdf_offset,
																			    ~weights,
																			    ~thisparticle_sumweight,
																			    (FLOAT)op.min_diff2[ipart],
																			    sp.nr_oversampled_rot,
																			    sp.nr_oversampled_trans,
																			    ~rot_id,
																			    ~trans_idx,
																				~job_idx,
																			 	~job_num,
																			 	job_num.size);
					CUDA_GPU_TAC("cuda_kernel_sumweight");
					CUDA_GPU_TIC("sumweightMemCp2");
					weights.cp_to_host();  //FIXME remove when mapping is eliminated - NOTE ALOT OF MWEIGHT-DEPS  BELOW
					//weights.free_device();
					CUDA_GPU_TAC("sumweightMemCp2");

					for (long unsigned k = 0; k< weights.size; k++)
						DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden_overs[k]) = weights[k];
				}



				thrust::device_ptr<FLOAT> dp = thrust::device_pointer_cast(~thisparticle_sumweight);
				exp_thisparticle_sumweight += thrust::reduce(dp, dp + block_num);
				thisparticle_sumweight.free_device();

				CUDA_CPU_TOC("sumweight1");
			} // end loop exp_iclass
		} // end if iter==1

		//Store parameters for this particle
		op.sum_weight[ipart] = exp_thisparticle_sumweight;
		std::cerr << "  sumweight =  " << exp_thisparticle_sumweight << std::endl;

#if defined(DEBUG_CUDA) && defined(__linux__)
		if (exp_thisparticle_sumweight == 0. || std::isnan(exp_thisparticle_sumweight))
		{
			printf("DEBUG_ERROR: zero sum of weights.\n");
			raise(SIGSEGV);
		}
#endif

	} // end loop ipart
	if (exp_ipass==0)
	{
		op.Mcoarse_significant.resize(sp.nr_particles, XSIZE(op.Mweight));
	}

	CUDA_CPU_TIC("convert_post_kernel");
	// Now, for each particle,  find the exp_significant_weight that encompasses adaptive_fraction of op.sum_weight
	op.significant_weight.clear();
	op.significant_weight.resize(sp.nr_particles, 0.);
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		MultidimArray<FLOAT> sorted_weight;
		// Get the relevant row for this particle
		op.Mweight.getRow(ipart, sorted_weight);

		// Only select non-zero probabilities to speed up sorting // TODO Remove when mapping is eliminated
		long int np = 0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
		{
			if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) > 0.)
			{
//				std::cerr << "weights = " << DIRECT_A2D_ELEM(op.Mweight, ipart, n) << std::endl;
				DIRECT_MULTIDIM_ELEM(sorted_weight, np) = DIRECT_MULTIDIM_ELEM(sorted_weight, n);
				np++;
			}
		}
		sorted_weight.resize(np);

		// Sort from low to high values
		CUDA_CPU_TIC("sort");
//		std::cerr << "sort on " << sorted_weight.xdim << " which should have np = " << np << std::endl;
#if  defined(USE_THRUST) && !defined(CUDA_DOUBLE_PRECISION) // Thrust seems incredibly slow in debug build this is clearly a FIXME
		thrust::sort(sorted_weight.data, sorted_weight.data + np );
#else
		sorted_weight.sort();
#endif
		CUDA_CPU_TOC("sort");

		double frac_weight = 0.;
		double my_significant_weight;
		long int my_nr_significant_coarse_samples = 0;
		for (long int i = XSIZE(sorted_weight) - 1; i >= 0; i--)
		{
			if (exp_ipass==0) my_nr_significant_coarse_samples++;
			my_significant_weight = DIRECT_A1D_ELEM(sorted_weight, i);
			//std::cerr << "thisweight = " << my_significant_weight << std::endl;
			frac_weight += my_significant_weight;
			if (frac_weight > baseMLO->adaptive_fraction * op.sum_weight[ipart])
				break;
		}

		if (exp_ipass==0 && my_nr_significant_coarse_samples == 0)
		{
			std::cerr << " ipart= " << ipart << " adaptive_fraction= " << baseMLO->adaptive_fraction << std::endl;
			std::cerr << " frac-weight= " << frac_weight << std::endl;
			std::cerr << " op.sum_weight[ipart]= " << op.sum_weight[ipart] << std::endl;
			Image<FLOAT> It;
			std::cerr << " XSIZE(op.Mweight)= " << XSIZE(op.Mweight) << std::endl;
			It()=op.Mweight;
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
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NR_SIGN) = (double)my_nr_significant_coarse_samples;

			// Keep track of which coarse samplings were significant were significant for this particle
			for (int ihidden = 0; ihidden < XSIZE(op.Mcoarse_significant); ihidden++)
			{
				if (DIRECT_A2D_ELEM(op.Mweight, ipart, ihidden) >= my_significant_weight)
					DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden) = true;
				else
					DIRECT_A2D_ELEM(op.Mcoarse_significant, ipart, ihidden) = false;
			}

		}
		op.significant_weight[ipart] = my_significant_weight;
		//std::cerr << "@sort op.significant_weight[ipart]= " << (FLOAT)op.significant_weight[ipart] << std::endl;

	} // end loop ipart
	CUDA_CPU_TOC("convert_post_kernel");

}

void storeWeightedSums(OptimisationParamters &op, SamplingParameters &sp,
		MlOptimiser *baseMLO,
		CudaGlobalPtr <long unsigned> &rot_id,
		CudaGlobalPtr <long unsigned> &rot_idx,
		CudaGlobalPtr <long unsigned> &trans_idx,
		CudaGlobalPtr <long unsigned> &ihidden_overs,
		CudaGlobalPtr <long unsigned> &job_idx,
		CudaGlobalPtr <long unsigned> &job_num,
		CudaGlobalPtr <FLOAT> 		  &weights)
{
	CUDA_CPU_TIC("store_init");

	// Re-do below because now also want unmasked images AND if (stricht_highres_exp >0.) then may need to resize
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(true, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, op.Fimgs_nomask, op.Fctfs, op.local_Fimgs_shifted, op.local_Fimgs_shifted_nomask,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);

	// In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the op.local_Minvsigma2s was omitted.
	// Set those back here
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);
		DIRECT_MULTIDIM_ELEM(op.local_Minvsigma2s[ipart], 0) = 1. / (baseMLO->sigma2_fudge * DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[group_id], 0));
	}

	// Initialise the maximum of all weights to a negative value
	op.max_weight.clear();
	op.max_weight.resize(sp.nr_particles, -1.);

	// For norm_correction and scale_correction of all particles of this ori_particle
	std::vector<double> exp_wsum_norm_correction;
	std::vector<MultidimArray<double> > exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA;
	std::vector<MultidimArray<double> > thr_wsum_signal_product_spectra, thr_wsum_reference_power_spectra;
	exp_wsum_norm_correction.resize(sp.nr_particles, 0.);

	// For scale_correction
	if (baseMLO->do_scale_correction)
	{
		MultidimArray<double> aux;
		aux.initZeros(baseMLO->mymodel.ori_size/2 + 1);
		exp_wsum_scale_correction_XA.resize(sp.nr_particles, aux);
		exp_wsum_scale_correction_AA.resize(sp.nr_particles, aux);
		thr_wsum_signal_product_spectra.resize(baseMLO->mymodel.nr_groups, aux);
		thr_wsum_reference_power_spectra.resize(baseMLO->mymodel.nr_groups, aux);
	}


	std::vector<double> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
	Matrix2D<double> A;
	MultidimArray<FLOAT > Fimg_real, Fimg_imag;
	MultidimArray<Complex > Fimg, Fimg_otfshift_nomask;  //TODO remove, currently needed for Fourier stuff, which is based on the complex class
	MultidimArray<double> Fweight, Minvsigma2, Mctf;
	bool have_warned_small_scale = false;

	Fimg_real.resize(op.Fimgs[0]);
	Fimg_imag.resize(op.Fimgs[0]);
	Fimg.resize(op.Fimgs[0]);
	Fweight.resize(op.Fimgs[0]);

	// Initialise Mctf to all-1 for if !do_ctf_corection
	Mctf.resize(op.Fimgs[0]);
	Mctf.initConstant(1.);
	// Initialise Minvsigma2 to all-1 for if !do_map
	Minvsigma2.resize(op.Fimgs[0]);
	Minvsigma2.initConstant(1.);

	// Make local copies of weighted sums (except BPrefs, which are too big)
	// so that there are not too many mutex locks below
	std::vector<MultidimArray<double> > thr_wsum_sigma2_noise, thr_wsum_pdf_direction;
	std::vector<double> thr_wsum_norm_correction, thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class;
	double thr_wsum_sigma2_offset;
	MultidimArray<double> thr_metadata, zeroArray;
	// Wsum_sigma_noise2 is a 1D-spectrum for each group
	zeroArray.initZeros(baseMLO->mymodel.ori_size/2 + 1);
	thr_wsum_sigma2_noise.resize(baseMLO->mymodel.nr_groups, zeroArray);
	// wsum_pdf_direction is a 1D-array (of length sampling.NrDirections()) for each class
	zeroArray.initZeros(baseMLO->sampling.NrDirections());
	thr_wsum_pdf_direction.resize(baseMLO->mymodel.nr_classes, zeroArray);
	// sumw_group is a double for each group
	thr_sumw_group.resize(baseMLO->mymodel.nr_groups, 0.);
	// wsum_pdf_class is a double for each class
	thr_wsum_pdf_class.resize(baseMLO->mymodel.nr_classes, 0.);
	if (baseMLO->mymodel.ref_dim == 2)
	{
		thr_wsum_prior_offsetx_class.resize(baseMLO->mymodel.nr_classes, 0.);
		thr_wsum_prior_offsety_class.resize(baseMLO->mymodel.nr_classes, 0.);
	}
	// wsum_sigma2_offset is just a double
	thr_wsum_sigma2_offset = 0.;

	unsigned image_size = op.Fimgs[0].nzyxdim;

	unsigned proj_div_max_count(4096*2);

	CUDA_CPU_TOC("store_init");

	// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
	{

		for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
		{
			long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
			int group_id = baseMLO->mydata.getGroupId(part_id);

			double myprior_x, myprior_y, myprior_z;
			double old_offset_x = XX(op.old_offset[ipart]);
			double old_offset_y = YY(op.old_offset[ipart]);
			double old_offset_z;

			if (baseMLO->mymodel.ref_dim == 2)
			{
				myprior_x = XX(baseMLO->mymodel.prior_offset_class[exp_iclass]);
				myprior_y = YY(baseMLO->mymodel.prior_offset_class[exp_iclass]);
			}
			else
			{
				myprior_x = XX(op.prior[ipart]);
				myprior_y = YY(op.prior[ipart]);
				if (baseMLO->mymodel.data_dim == 3)
				{
					myprior_z = ZZ(op.prior[ipart]);
					old_offset_z = ZZ(op.old_offset[ipart]);
				}
			}

			/*=======================================================================================
												   COLLECT 2
			=======================================================================================*/


			CUDA_CPU_TIC("collect_data_2");

			CUDA_CPU_TIC("collect_data_2_pre_kernel");
			//TODO should be replaced with loop over pairs of projections and translations (like in the getAllSquaredDifferences-function)

			std::vector< double> oversampled_rot, oversampled_tilt, oversampled_psi;

			int oversamples = sp.nr_oversampled_trans * sp.nr_oversampled_rot;
			CudaGlobalPtr<FLOAT >  Mweight( &(op.Mweight.data[(ipart)*(op.Mweight).xdim]),
											sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples);
			int nr_transes = sp.nr_trans*sp.nr_oversampled_trans;
			CudaGlobalPtr<FLOAT>     oo_otrans_x(nr_transes); // old_offset_oversampled_trans_x
			CudaGlobalPtr<FLOAT>     oo_otrans_y(nr_transes);
			CudaGlobalPtr<FLOAT> myp_oo_otrans_x2y2z2(nr_transes); // my_prior_old_offs....x^2*y^2*z^2

			//Pregenerate oversampled translation objects for kernel-call
			for (long int itrans = 0, iitrans = 0; itrans < sp.nr_trans; itrans++)
			{
				baseMLO->sampling.getTranslations(itrans, baseMLO->adaptive_oversampling,
						oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
				for (long int iover_trans = 0; iover_trans < sp.nr_oversampled_trans; iover_trans++, iitrans++)
				{
					oo_otrans_x[iitrans] = old_offset_x + oversampled_translations_x[iover_trans];
					oo_otrans_y[iitrans] = old_offset_y + oversampled_translations_y[iover_trans];
					double diffx = myprior_x - oo_otrans_x[iitrans];
					double diffy = myprior_y - oo_otrans_y[iitrans];
					if (baseMLO->mymodel.data_dim == 3)
					{
						double diffz = myprior_z - (old_offset_z + oversampled_translations_z[iover_trans]);
						myp_oo_otrans_x2y2z2[iitrans] = diffx*diffx + diffy*diffy + diffz*diffz ;
					}
					else
					{
						myp_oo_otrans_x2y2z2[iitrans] = diffx*diffx + diffy*diffy ;
					}
				}
			}
			long int jobid=0;
			job_idx[jobid]=0;
			job_num[jobid]=1;
			long int crot = rot_idx[jobid]; // set current rot
			for(long int n=1; n<rot_idx.size; n++)
			{
				if(rot_idx[n]==crot)
				{
					job_num[jobid]++;
				}
				else
				{
					jobid++;
					job_num[jobid]=1;
					job_idx[jobid]=n;
					crot=rot_idx[n];
				}
			}
			job_idx.size=jobid;
			job_num.size=jobid;
			job_idx.cp_to_device();
			job_num.cp_to_device();

			oo_otrans_x.put_on_device();
			oo_otrans_y.put_on_device();
			myp_oo_otrans_x2y2z2.put_on_device();

			int block_num = job_idx.size;

			std::cerr << "block_num = " << block_num << std::endl;
			CudaGlobalPtr<FLOAT>                      p_weights(block_num);
			CudaGlobalPtr<FLOAT> p_thr_wsum_prior_offsetx_class(block_num);
			CudaGlobalPtr<FLOAT> p_thr_wsum_prior_offsety_class(block_num);
			CudaGlobalPtr<FLOAT>       p_thr_wsum_sigma2_offset(block_num);
			unsigned long coarse_nr = sp.nr_dir*sp.nr_psi;
			p_weights.device_alloc();
			p_thr_wsum_prior_offsetx_class.device_alloc();
			p_thr_wsum_prior_offsety_class.device_alloc();
			p_thr_wsum_sigma2_offset.device_alloc();
//			weights.device_alloc();
//			weights.cp_to_device();

			CUDA_CPU_TOC("collect_data_2_pre_kernel");
			CUDA_GPU_TIC("collect2-kernel");

			dim3 grid_dim_collect2 = splitCudaBlocks(block_num,false);
			cuda_kernel_collect2jobs<<<grid_dim_collect2,SUM_BLOCK_SIZE>>>(
					~oo_otrans_x,          // otrans-size -> make const
					~oo_otrans_y,          // otrans-size -> make const
					~myp_oo_otrans_x2y2z2, // otrans-size -> make const
					~weights,
					(FLOAT)op.significant_weight[ipart],
					(FLOAT)op.sum_weight[ipart],
					sp.nr_trans,
					sp.nr_oversampled_trans,
					sp.nr_oversampled_rot,
					oversamples,
					(baseMLO->do_skip_align || baseMLO->do_skip_rotate ),
					~p_weights,
					~p_thr_wsum_prior_offsetx_class,
					~p_thr_wsum_prior_offsety_class,
					~p_thr_wsum_sigma2_offset,
					~rot_idx,
					~trans_idx,
					~job_idx,
					~job_num
						);
			CUDA_GPU_TAC("collect2-kernel");

			CUDA_GPU_TIC("cuda_kernel_collect2_memcpy2");
			// TODO further reduce the below 4 arrays while data is still on gpu
			p_weights.cp_to_host();
			p_thr_wsum_prior_offsetx_class.cp_to_host();
			p_thr_wsum_prior_offsety_class.cp_to_host();
			p_thr_wsum_sigma2_offset.cp_to_host();
			CUDA_GPU_TAC("cuda_kernel_collect2_memcpy2");

			CUDA_CPU_TIC("collect_data_2_post_kernel");
			HANDLE_ERROR(cudaDeviceSynchronize());
			CUDA_GPU_TOC();
			thr_wsum_sigma2_offset = 0.0;
			int iorient = 0;
			for (long int n = 0; n < block_num; n++)
			{
				iorient= rot_id[job_idx[n]];
				long int iorientclass = exp_iclass * sp.nr_dir * sp.nr_psi + iorient;
				// Only proceed if any of the particles had any significant coarsely sampled translation

				if (baseMLO->isSignificantAnyParticleAnyTranslation(iorientclass, sp.itrans_min, sp.itrans_max, op.Mcoarse_significant))
				{
					long int mydir, idir=floor(iorient/sp.nr_psi);
					if (baseMLO->mymodel.orientational_prior_mode == NOPRIOR)
						mydir = idir;
					else
						mydir = op.pointer_dir_nonzeroprior[idir];

					// store partials according to indices of the relevant dimension
					DIRECT_MULTIDIM_ELEM(thr_wsum_pdf_direction[exp_iclass], mydir) += p_weights[n];
					thr_sumw_group[group_id]                 						+= p_weights[n];
					thr_wsum_pdf_class[exp_iclass]           						+= p_weights[n];
					thr_wsum_sigma2_offset                   						+= p_thr_wsum_sigma2_offset[n];

					if (baseMLO->mymodel.ref_dim == 2)
					{
						thr_wsum_prior_offsetx_class[exp_iclass] 	+= p_thr_wsum_prior_offsetx_class[iorient];
						thr_wsum_prior_offsety_class[exp_iclass] 	+= p_thr_wsum_prior_offsety_class[iorient];
					}
				}
			}
			Mweight.free_device();
			p_weights.free();
			p_thr_wsum_sigma2_offset.free();
			p_thr_wsum_prior_offsetx_class.free();
			p_thr_wsum_prior_offsety_class.free();
			oo_otrans_y.free();
			oo_otrans_x.free();
			myp_oo_otrans_x2y2z2.free();
			CUDA_CPU_TOC("collect_data_2_post_kernel");

			/*=======================================================================================
												  SET META DATA
			=======================================================================================*/

			CUDA_CPU_TIC("setMetadata");

			//Get index of max element using GPU-tool thrust
			Indices max_index;
#if !defined(CUDA_DOUBLE_PRECISION) && defined(USE_THRUST)
			max_index.fineIdx = thrust::max_element(&DIRECT_A2D_ELEM(op.Mweight, ipart, 0),&DIRECT_A2D_ELEM(op.Mweight, ipart+1, 0)) - &DIRECT_A2D_ELEM(op.Mweight, ipart, 0);
#else
			max_index.fineIdx = std::max_element(&DIRECT_A2D_ELEM(op.Mweight, ipart, 0),&DIRECT_A2D_ELEM(op.Mweight, ipart+1, 0)) - &DIRECT_A2D_ELEM(op.Mweight, ipart, 0);
#endif
			//std::cerr << "max index = " << max_index.fineIdx << std::endl;

			op.max_weight[ipart] = DIRECT_A2D_ELEM(op.Mweight, ipart, max_index.fineIdx);
			max_index.fineIndexToFineIndices(sp); // set partial indices corresponding to the found max_index, to be used below

			baseMLO->sampling.getTranslations(max_index.itrans, baseMLO->adaptive_oversampling,
					oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);
			baseMLO->sampling.getOrientations(max_index.idir, max_index.ipsi, baseMLO->adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
					op.pointer_dir_nonzeroprior, op.directions_prior, op.pointer_psi_nonzeroprior, op.psi_prior);
			double rot = oversampled_rot[max_index.ioverrot];
			double tilt = oversampled_tilt[max_index.ioverrot];
			double psi = oversampled_psi[max_index.ioverrot];
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ROT) = rot;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_TILT) = tilt;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PSI) = psi;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_XOFF) = XX(op.old_offset[ipart]) + oversampled_translations_x[max_index.iovertrans];
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_YOFF) = YY(op.old_offset[ipart]) + oversampled_translations_y[max_index.iovertrans];
			if (baseMLO->mymodel.data_dim == 3)
				DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_ZOFF) = ZZ(op.old_offset[ipart]) + oversampled_translations_z[max_index.iovertrans];
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_CLASS) = (double)max_index.iclass + 1;
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PMAX) = op.max_weight[ipart]/op.sum_weight[ipart];

			CUDA_CPU_TOC("setMetadata");


			CUDA_CPU_TOC("collect_data_2");

		}

		/*=======================================================================================
									          MAXIMIZATION
		=======================================================================================*/


		std::vector< long unsigned > iorientclasses_all, iover_rots_all;
		std::vector< double > rots_all, tilts_all, psis_all;

		CUDA_CPU_TIC("generateProjectionSetup_2");

		long unsigned orientation_num_all = generateProjectionSetup(
					op,
					sp,
					baseMLO,
					false,  //coarse
					exp_iclass,
					rots_all, tilts_all, psis_all,
					iorientclasses_all,
					iover_rots_all);

		CUDA_CPU_TOC("generateProjectionSetup_2");

		unsigned proj_div_nr = ceil((float)orientation_num_all / (float)proj_div_max_count);

		for (int iproj_div = 0; iproj_div < proj_div_nr; iproj_div++)
		{
			unsigned long proj_div_start(proj_div_max_count * iproj_div),
					proj_div_end;

			if (iproj_div < proj_div_nr - 1)
				proj_div_end = proj_div_start + proj_div_max_count;
			else
				proj_div_end = orientation_num_all;

			long unsigned orientation_num(proj_div_end - proj_div_start);

			std::vector< long unsigned >
				iorientclasses(&iorientclasses_all[proj_div_start],&iorientclasses_all[proj_div_end]),
				iover_rots(&iover_rots_all[proj_div_start],&iover_rots_all[proj_div_end]);

			std::vector< double >
				rots(&rots_all[proj_div_start],&rots_all[proj_div_end]),
				tilts(&tilts_all[proj_div_start],&tilts_all[proj_div_end]),
				psis(&psis_all[proj_div_start],&psis_all[proj_div_end]);


			/*======================================================
								 PROJECTIONS
			======================================================*/

			// Since we will need the euler matrices for all projections in the data_collect stage,
			// we might as well make it wider in scope and retain it on the GPU until then. When we
			// switch from pair to bool, there won't be any need to remake it every class, but for
			// now we create only those matrices corresponding to significant orientations, which IS  * class-specific *


			CUDA_CPU_TIC("generateEulerMatrices");
			CudaGlobalPtr<FLOAT> eulers(9 * orientation_num);

			generateEulerMatrices(
					baseMLO->mymodel.PPref[exp_iclass].padding_factor,
					rots,
					tilts,
					psis,
					eulers,
					!IS_NOT_INV);



			eulers.device_alloc();
			eulers.cp_to_device();

			CUDA_CPU_TOC("generateEulerMatrices");

			CudaGlobalPtr<FLOAT> wavgs_real(orientation_num * image_size);
			wavgs_real.device_alloc();
			wavgs_real.device_init(0);
			CudaGlobalPtr<FLOAT> wavgs_imag(orientation_num * image_size);
			wavgs_imag.device_alloc();
			wavgs_imag.device_init(0);
			CudaGlobalPtr<FLOAT> Fweights(orientation_num * image_size);
			Fweights.device_alloc();
			Fweights.device_init(0);

			/// Now that reference projection has been made loop over all particles inside this ori_particle
			for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
			{
				long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
				int group_id = baseMLO->mydata.getGroupId(part_id);

				/*======================================================
									 TRANSLATIONS
				======================================================*/

				CUDA_CPU_TIC("translation_2");

				CudaGlobalPtr<FLOAT> Fimgs_real(image_size * sp.nr_trans * sp.nr_oversampled_trans);
				CudaGlobalPtr<FLOAT> Fimgs_imag(Fimgs_real.size);
				CudaGlobalPtr<FLOAT> Fimgs_nomask_real(Fimgs_real.size);
				CudaGlobalPtr<FLOAT> Fimgs_nomask_imag(Fimgs_real.size);

				std::vector< long unsigned > iover_transes, itranses, ihiddens;

				long unsigned translation_num = imageTranslation(
						Fimgs_real,
						Fimgs_imag,
						Fimgs_nomask_real,
						Fimgs_nomask_imag,
						sp.itrans_min,
						sp.itrans_max,
						baseMLO->adaptive_oversampling,
						baseMLO->sampling,
						oversampled_translations_x,
						oversampled_translations_y,
						oversampled_translations_z,
						sp.nr_oversampled_trans,
						baseMLO->global_fftshifts_ab_current,
						baseMLO->global_fftshifts_ab2_current,
						op.local_Fimgs_shifted[ipart],
						op.local_Fimgs_shifted_nomask[ipart],
						iover_transes,
						itranses,
						ihiddens,
						image_size);

				Fimgs_real.device_alloc();
				Fimgs_real.cp_to_device();
				Fimgs_imag.device_alloc();
				Fimgs_imag.cp_to_device();

				Fimgs_nomask_real.device_alloc();
				Fimgs_nomask_real.cp_to_device();
				Fimgs_nomask_imag.device_alloc();
				Fimgs_nomask_imag.cp_to_device();

				CUDA_CPU_TOC("translation_2");


				/*======================================================
										SCALE
				======================================================*/

				CUDA_CPU_TIC("scale_ctf");
				FLOAT part_scale(1.);

				if (baseMLO->do_scale_correction)
				{
					part_scale = baseMLO->mymodel.scale_correction[group_id];
					if (part_scale > 10000.)
					{
						std::cerr << " rlnMicrographScaleCorrection= " << part_scale << " group= " << group_id + 1 << std::endl;
						REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
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

				CudaGlobalPtr<FLOAT> ctfs(image_size); //TODO Same size for all iparts, should be allocated once
				ctfs.device_alloc();

				if (baseMLO->do_ctf_correction)
				{
					for (unsigned i = 0; i < image_size; i++)
						ctfs[i] = (FLOAT) op.local_Fctfs[ipart].data[i] * part_scale;
				}
				else //TODO should be handled by memset
					for (unsigned i = 0; i < image_size; i++)
						ctfs[i] = part_scale;

				ctfs.cp_to_device();
				CUDA_CPU_TOC("scale_ctf");


				/*======================================================
									MAP WEIGHTS
				======================================================*/

				CUDA_CPU_TIC("map");
				CudaGlobalPtr<FLOAT> sorted_weights(orientation_num * translation_num);

				mapWeights(
						sorted_weights,
						orientation_num,
						translation_num,
						baseMLO->sampling,
						ipart,
						iover_transes,
						ihiddens,
						iorientclasses,
						iover_rots,
						op.Mweight,
						sp.current_oversampling,
						sp.nr_trans);

				sorted_weights.device_alloc();
				sorted_weights.cp_to_device();

				CUDA_CPU_TOC("map");

				/*======================================================
									KERNEL CALL
				======================================================*/

				CudaGlobalPtr<FLOAT> Minvsigma2s(image_size); //TODO Same size for all iparts, should be allocated once
				Minvsigma2s.device_alloc();
				CudaGlobalPtr<FLOAT> wdiff2s_parts(orientation_num * image_size);
				wdiff2s_parts.device_alloc();

				if (baseMLO->do_map)
					for (unsigned i = 0; i < image_size; i++)
						Minvsigma2s[i] = op.local_Minvsigma2s[ipart].data[i];
				else //TODO should be handled by memset
					for (unsigned i = 0; i < image_size; i++)
						Minvsigma2s[i] = 1;

				Minvsigma2s.cp_to_device();

				Cuda3DProjectorKernel projKernel = Cuda3DProjectorKernel::makeKernel(
						baseMLO->cudaProjectors[exp_iclass],
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						op.local_Minvsigma2s[0].xdim-1);

				runWavgKernel(
						projKernel,
						eulers,
						Fimgs_real,
						Fimgs_imag,
						Fimgs_nomask_real,
						Fimgs_nomask_imag,
						sorted_weights,
						ctfs,
						Minvsigma2s,
						wdiff2s_parts,
						wavgs_real,
						wavgs_imag,
						Fweights,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						image_size,
						ipart,
						group_id,
						exp_iclass
						);

				Fimgs_real.free_device();
				Fimgs_imag.free_device();
				Fimgs_nomask_real.free_device();
				Fimgs_nomask_imag.free_device();

				sorted_weights.free_device();
				ctfs.free_device();
				Minvsigma2s.free_device();

				/*======================================================
									REDUCE WDIFF2S
				======================================================*/

				CUDA_CPU_TIC("reduce_wdiff2s");
				// reduction_block_num = the highest possible power of two that covers more than or exactly half of all images to be reduced
				int num_reductions = (int)floor(log2((float)orientation_num));
				int reduction_block_num = pow(2,num_reductions);
				if(reduction_block_num==orientation_num) // (possibly) very special case where orientation_num is a power of 2
					reduction_block_num /= 2;

				CUDA_GPU_TIC("cuda_kernels_reduce_wdiff2s");

				for(int k=reduction_block_num; k>=1; k/=2) //invoke kernel repeatedly until all images have been stacked into the first image position
				{

					dim3 block_dim_wd = splitCudaBlocks(k,true);

					// TODO **OF VERY LITTLE IMPORTANCE**  One block treating just 2 images is a very inefficient amount of loads per store
					cuda_kernel_reduce_wdiff2s<<<block_dim_wd,BLOCK_SIZE>>>(~wdiff2s_parts,orientation_num,image_size,k);
				}

				CUDA_GPU_TAC("cuda_kernels_reduce_wdiff2s");

				wdiff2s_parts.size = image_size; //temporarily set the size to the single image we have now reduced, to not copy more than necessary
				wdiff2s_parts.cp_to_host();
				wdiff2s_parts.size = orientation_num * image_size;
				wdiff2s_parts.free_device();

				HANDLE_ERROR(cudaDeviceSynchronize());

				CUDA_GPU_TOC();

				for (long int j = 0; j < image_size; j++)
				{
					int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine, j);
					if (ires > -1)
					{
						thr_wsum_sigma2_noise[group_id].data[ires] += (double) wdiff2s_parts[j];
						exp_wsum_norm_correction[ipart] += (double) wdiff2s_parts[j]; //TODO could be thrust-reduced
					}
				}
				wdiff2s_parts.free_host();

				CUDA_CPU_TOC("reduce_wdiff2s");

			} // end loop ipart

			/*======================================================
								BACKPROJECTION
			======================================================*/

			CUDA_CPU_TIC("backprojection");

			CudaGlobalPtr<FLOAT> bp_eulers(9 * orientation_num);

			FLOAT padding_factor = baseMLO->wsum_model.BPref[exp_iclass].padding_factor;

			generateEulerMatrices(
					1/padding_factor, //Why squared scale factor is given in backprojection
					rots,
					tilts,
					psis,
					bp_eulers,
					IS_NOT_INV);

			bp_eulers.device_alloc();
			bp_eulers.cp_to_device();
			bp_eulers.free_host();

			int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
			pthread_mutex_lock(&global_mutex2[my_mutex]);

			baseMLO->cudaBackprojectors[exp_iclass].backproject(
					~wavgs_real,
					~wavgs_imag,
					~Fweights,
					~bp_eulers,
					op.local_Minvsigma2s[0].xdim,
					op.local_Minvsigma2s[0].ydim,
					orientation_num);

			pthread_mutex_unlock(&global_mutex2[my_mutex]);



			CUDA_CPU_TOC("backprojection");
		}

	} // end loop iclass

	CUDA_CPU_TIC("store_post_gpu");

	// Extend norm_correction and sigma2_noise estimation to higher resolutions for all particles
	// Also calculate dLL for each particle and store in metadata
	// loop over all particles inside this ori_particle
	double thr_avg_norm_correction = 0.;
	double thr_sum_dLL = 0., thr_sum_Pmax = 0.;
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);

		// If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
		for (int ires = baseMLO->mymodel.current_size/2 + 1; ires < baseMLO->mymodel.ori_size/2 + 1; ires++)
		{
			DIRECT_A1D_ELEM(thr_wsum_sigma2_noise[group_id], ires) += DIRECT_A1D_ELEM(op.power_imgs[ipart], ires);
			// Also extend the weighted sum of the norm_correction
			exp_wsum_norm_correction[ipart] += DIRECT_A1D_ELEM(op.power_imgs[ipart], ires);
		}

		// Store norm_correction
		// Multiply by old value because the old norm_correction term was already applied to the image
		if (baseMLO->do_norm_correction)
		{
			double old_norm_correction = DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM);
			old_norm_correction /= baseMLO->mymodel.avg_norm_correction;
			// The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
			// The variance of the total image (on which one normalizes) is twice this value!
			double normcorr = old_norm_correction * sqrt(exp_wsum_norm_correction[ipart] * 2.);
			thr_avg_norm_correction += normcorr;

			// Now set the new norm_correction in the relevant position of exp_metadata
			DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) = normcorr;


			// Print warning for strange norm-correction values
			if (!(baseMLO->iter == 1 && baseMLO->do_firstiter_cc) && DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) > 10.)
			{
				std::cout << " WARNING: norm_correction= "<< DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_NORM) << " for particle " << part_id << " in group " << group_id + 1 << "; Are your groups large enough?" << std::endl;
			}

		}

		// Store weighted sums for scale_correction
		if (baseMLO->do_scale_correction)
		{
			// Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
			exp_wsum_scale_correction_XA[ipart] /= baseMLO->mymodel.scale_correction[group_id];
			exp_wsum_scale_correction_AA[ipart] /= baseMLO->mymodel.scale_correction[group_id] * baseMLO->mymodel.scale_correction[group_id];

			thr_wsum_signal_product_spectra[group_id] += exp_wsum_scale_correction_XA[ipart];
			thr_wsum_reference_power_spectra[group_id] += exp_wsum_scale_correction_AA[ipart];
		}

		// Calculate DLL for each particle
		double logsigma2 = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(baseMLO->Mresol_fine)
		{
			int ires = DIRECT_MULTIDIM_ELEM(baseMLO->Mresol_fine, n);
			// Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
			// Also exclude origin from logsigma2, as this will not be considered in the P-calculations
			if (ires > 0)
				logsigma2 += log( 2. * PI * DIRECT_A1D_ELEM(baseMLO->mymodel.sigma2_noise[group_id], ires));
		}
		if (op.sum_weight[ipart]==0)
		{
			std::cerr << " part_id= " << part_id << std::endl;
			std::cerr << " ipart= " << ipart << std::endl;
			std::cerr << " op.min_diff2[ipart]= " << op.min_diff2[ipart] << std::endl;
			std::cerr << " logsigma2= " << logsigma2 << std::endl;
			int group_id = baseMLO->mydata.getGroupId(part_id);
			std::cerr << " group_id= " << group_id << std::endl;
			std::cerr << " ml_model.scale_correction[group_id]= " << baseMLO->mymodel.scale_correction[group_id] << std::endl;
			std::cerr << " exp_significant_weight[ipart]= " << op.significant_weight[ipart] << std::endl;
			std::cerr << " exp_max_weight[ipart]= " << op.max_weight[ipart] << std::endl;
			std::cerr << " ml_model.sigma2_noise[group_id]= " << baseMLO->mymodel.sigma2_noise[group_id] << std::endl;
			REPORT_ERROR("ERROR: op.sum_weight[ipart]==0");
		}
		double dLL;
		if ((baseMLO->iter==1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc)
			dLL = -op.min_diff2[ipart];
		else
			dLL = log(op.sum_weight[ipart]) - op.min_diff2[ipart] - logsigma2;

		// Store dLL of each image in the output array, and keep track of total sum
		DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_DLL) = dLL;
		thr_sum_dLL += dLL;

		// Also store sum of Pmax
		thr_sum_Pmax += DIRECT_A2D_ELEM(baseMLO->exp_metadata, op.metadata_offset + ipart, METADATA_PMAX);

	}

	// Now, inside a global_mutex, update the other weighted sums among all threads
	if (!baseMLO->do_skip_maximization)
	{
		pthread_mutex_lock(&global_mutex);
		for (int n = 0; n < baseMLO->mymodel.nr_groups; n++)
		{
			baseMLO->wsum_model.sigma2_noise[n] += thr_wsum_sigma2_noise[n];
			baseMLO->wsum_model.sumw_group[n] += thr_sumw_group[n];
			if (baseMLO->do_scale_correction)
			{
				baseMLO->wsum_model.wsum_signal_product_spectra[n] += thr_wsum_signal_product_spectra[n];
				baseMLO->wsum_model.wsum_reference_power_spectra[n] += thr_wsum_reference_power_spectra[n];
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

			if (!(baseMLO->do_skip_align || baseMLO->do_skip_rotate) )
				baseMLO->wsum_model.pdf_direction[n] += thr_wsum_pdf_direction[n];
		}
		baseMLO->wsum_model.sigma2_offset += thr_wsum_sigma2_offset;
		if (baseMLO->do_norm_correction)
			baseMLO->wsum_model.avg_norm_correction += thr_avg_norm_correction;
		baseMLO->wsum_model.LL += thr_sum_dLL;
		baseMLO->wsum_model.ave_Pmax += thr_sum_Pmax;
		pthread_mutex_unlock(&global_mutex);
	} // end if !do_skip_maximization

	CUDA_CPU_TOC("store_post_gpu");
}


void MlOptimiserCuda::doThreadExpectationSomeParticles(unsigned thread_id)
{

	CUDA_CPU_TOC("interParticle");

	//put mweight allocation here
	size_t first_ipart = 0, last_ipart = 0;
	while (baseMLO->exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
		for (long unsigned ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
			CUDA_CPU_TIC("oneParticle");

			unsigned my_ori_particle = baseMLO->exp_my_first_ori_particle + ipart;
			SamplingParameters sp;
			sp.nr_particles = baseMLO->mydata.ori_particles[my_ori_particle].particles_id.size();

			OptimisationParamters op(sp.nr_particles, my_ori_particle);

			// In the first iteration, multiple seeds will be generated
			// A single random class is selected for each pool of images, and one does not marginalise over the orientations
			// The optimal orientation is based on signal-product (rather than the signal-intensity sensitive Gaussian)
			// If do_firstiter_cc, then first perform a single iteration with K=1 and cross-correlation criteria, afterwards

			// Decide which classes to integrate over (for random class assignment in 1st iteration)
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
					sp.iclass_min = sp.iclass_max = divide_equally_which_group(baseMLO->mydata.numberOfOriginalParticles(), baseMLO->mymodel.nr_classes, op.my_ori_particle);
				}
			}
			// Global exp_metadata array has metadata of all ori_particles. Where does my_ori_particle start?
			for (long int iori = baseMLO->exp_my_first_ori_particle; iori <= baseMLO->exp_my_last_ori_particle; iori++)
			{
				if (iori == my_ori_particle) break;
				op.metadata_offset += baseMLO->mydata.ori_particles[iori].particles_id.size();
			}

			baseMLO->getFourierTransformsAndCtfs(my_ori_particle, op.metadata_offset, op.Fimgs, op.Fimgs_nomask, op.Fctfs,
					op.old_offset, op.prior, op.power_imgs, op.highres_Xi2_imgs,
					op.pointer_dir_nonzeroprior, op.pointer_psi_nonzeroprior, op.directions_prior, op.psi_prior);

			if (baseMLO->do_realign_movies && baseMLO->movie_frame_running_avg_side > 0)
			{
				baseMLO->calculateRunningAveragesOfMovieFrames(my_ori_particle, op.Fimgs, op.power_imgs, op.highres_Xi2_imgs);
			}

			// To deal with skipped alignments/rotations
			if (baseMLO->do_skip_align)
			{
				sp.itrans_min = sp.itrans_max = sp.idir_min = sp.idir_max = sp.ipsi_min = sp.ipsi_max =
						my_ori_particle - baseMLO->exp_my_first_ori_particle;
			}
			else
			{
				sp.itrans_min = 0;
				sp.itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;

				if (baseMLO->do_skip_rotate)
				{
					sp.idir_min = sp.idir_max = sp.ipsi_min = sp.ipsi_max =
							my_ori_particle - baseMLO->exp_my_first_ori_particle;
				}
				else
				{
					sp.idir_min = sp.ipsi_min = 0;
					sp.idir_max = baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior) - 1;
					sp.ipsi_max = baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior ) - 1;
				}
			}

			// Initialise significant weight to minus one, so that all coarse sampling points will be handled in the first pass
			op.significant_weight.resize(sp.nr_particles, -1.);

			// Only perform a second pass when using adaptive oversampling
			int nr_sampling_passes = (baseMLO->adaptive_oversampling > 0) ? 2 : 1;

			CudaGlobalPtr<long unsigned> rot_id_C, 			rot_id_F,
			 	 	 	 	 	 	 	 rot_idx_C, 		rot_idx_F,
										 trans_idx_C,		trans_idx_F,
										 ihidden_overs_C,  	ihidden_overs_F,      // TODO to be removed WHEN mapping is eliminated
										 job_idx_C, 		job_idx_F,
										 job_num_C, 		job_num_F;

			CudaGlobalPtr<FLOAT> 		 weights_C,			weights_F;

			for (int ipass = 0; ipass < nr_sampling_passes; ipass++)
			{
				CUDA_CPU_TIC("weightPass");

				if (baseMLO->strict_highres_exp > 0.)
					// Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
					sp.current_image_size = baseMLO->coarse_size;
				else if (baseMLO->adaptive_oversampling > 0)
					// Use smaller images in the first pass, larger ones in the second pass
					sp.current_image_size = (ipass == 0) ? baseMLO->coarse_size : baseMLO->mymodel.current_size;
				else
					sp.current_image_size = baseMLO->mymodel.current_size;

				// Use coarse sampling in the first pass, oversampled one the second pass
				sp.current_oversampling = (ipass == 0) ? 0 : baseMLO->adaptive_oversampling;

				sp.nr_dir = (baseMLO->do_skip_align || baseMLO->do_skip_rotate) ? 1 : baseMLO->sampling.NrDirections(0, &op.pointer_dir_nonzeroprior);
				sp.nr_psi = (baseMLO->do_skip_align || baseMLO->do_skip_rotate) ? 1 : baseMLO->sampling.NrPsiSamplings(0, &op.pointer_psi_nonzeroprior);
				sp.nr_trans = (baseMLO->do_skip_align) ? 1 : baseMLO->sampling.NrTranslationalSamplings();
				sp.nr_oversampled_rot = baseMLO->sampling.oversamplingFactorOrientations(sp.current_oversampling);
				sp.nr_oversampled_trans = baseMLO->sampling.oversamplingFactorTranslations(sp.current_oversampling);

				CUDA_CPU_TIC("getAllSquaredDifferences");
				if (ipass == 0)
					getAllSquaredDifferencesCoarse(ipass, op, sp, baseMLO);
				else
					getAllSquaredDifferencesFine(ipass, op, sp, baseMLO, rot_id_F,
																		 rot_idx_F,
																		 trans_idx_F,
																		 ihidden_overs_F,
																		 job_idx_F,
																		 job_num_F,
																		 weights_F);

				CUDA_CPU_TOC("getAllSquaredDifferences");

				CUDA_CPU_TIC("convertAllSquaredDifferencesToWeights");
				if (ipass == 0)
					convertAllSquaredDifferencesToWeights(ipass, op, sp, baseMLO, rot_id_C,
																				  rot_idx_C,
																				  trans_idx_C,
																				  ihidden_overs_C,
																				  job_idx_C,
																				  job_num_C,
																				  weights_C);
				else
					convertAllSquaredDifferencesToWeights(ipass, op, sp, baseMLO, rot_id_F,
																				  rot_idx_F,
																				  trans_idx_F,
																				  ihidden_overs_F,
																				  job_idx_F,
																				  job_num_F,
																				  weights_F);

				CUDA_CPU_TOC("convertAllSquaredDifferencesToWeights");

				CUDA_CPU_TOC("weightPass");
			}

			// For the reconstruction step use mymodel.current_size!
			sp.current_image_size = baseMLO->mymodel.current_size;

			CUDA_CPU_TIC("storeWeightedSums");
			storeWeightedSums(op, sp, baseMLO, rot_id_F, rot_idx_F,trans_idx_F,ihidden_overs_F,job_idx_F,job_num_F,weights_F);
			CUDA_CPU_TOC("storeWeightedSums");

			CUDA_CPU_TOC("oneParticle");
		}
	}

	CUDA_CPU_TIC("interParticle");
//	exit(0);
}


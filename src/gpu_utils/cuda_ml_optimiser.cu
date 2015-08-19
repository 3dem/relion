#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_projector_plan.h"
#include "src/gpu_utils/cuda_benchmark_utils.cuh"
#include "src/gpu_utils/cuda_ml_optimiser.h"
#include "src/gpu_utils/cuda_kernels/helper.cuh"
#include "src/gpu_utils/cuda_kernels/diff2.cuh"
#include "src/gpu_utils/cuda_kernels/wavg.cuh"
#include "src/gpu_utils/cuda_utils.cuh"
#include "src/gpu_utils/cuda_helper_functions.cu"
#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/complex.h"
#include <fstream>
#include <cuda_runtime.h>
#include "src/parallel.h"
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <signal.h>

static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;

void getAllSquaredDifferencesCoarse(
		unsigned exp_ipass,
		OptimisationParamters &op,
		SamplingParameters &sp,
		MlOptimiser *baseMLO,
		MlOptimiserCuda *cudaMLO,
		std::vector<CudaProjectorPlan> &projectorPlans)
{

#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF1);
#endif

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
	CudaGlobalPtr<XFLOAT> allWeights(cudaMLO->allocator);
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		allWeights.size+=projectorPlans[exp_iclass].orientation_num * sp.nr_trans*sp.nr_oversampled_trans * sp.nr_particles;
	allWeights.device_alloc();

	long int allWeights_pos=0;
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		long int group_id = baseMLO->mydata.getGroupId(part_id);

		/*====================================
				Generate Translations
		======================================*/

		CUDA_CPU_TIC("translation_1");

		CudaGlobalPtr<XFLOAT> Fimgs_real(image_size * sp.nr_trans * sp.nr_oversampled_trans, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_imag(image_size * sp.nr_trans * sp.nr_oversampled_trans, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> gpuMinvsigma2(image_size, cudaMLO->allocator);
		gpuMinvsigma2.device_alloc();

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
					XFLOAT real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real
							- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag;
					XFLOAT imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag
							+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real;

					//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
					if (baseMLO->do_scale_correction)
					{
						XFLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
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
			XFLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
			{
				gpuMinvsigma2[n] *= myscale * myscale;
			}
		}
		gpuMinvsigma2.cp_to_device();
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			CudaProjectorPlan projectorPlan = projectorPlans[exp_iclass];

			if ( projectorPlan.orientation_num > 0 )
			{
				CudaGlobalPtr<XFLOAT> diff2s(projectorPlan.orientation_num*translation_num,cudaMLO->allocator);
//				diff2s.h_ptr = &allWeights.h_ptr[allWeights_pos];
				diff2s.d_ptr = &allWeights.d_ptr[allWeights_pos];
				diff2s.h_do_free=false;
				diff2s.d_do_free=false;


				/*====================================
				    	   Kernel Call
				======================================*/

				CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
						cudaMLO->cudaProjectors[exp_iclass],
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						op.local_Minvsigma2s[0].xdim-1);

				runDiff2KernelCoarse(
						projKernel,
						~gpuMinvsigma2,
						~Fimgs_real,
						~Fimgs_imag,
						~(*projectorPlan.eulers),
						~diff2s,
						op,
						baseMLO,
						projectorPlan.orientation_num,
						translation_num,
						image_size,
						ipart,
						group_id,
						exp_iclass);

				/*====================================
				    	   Retrieve Results
				======================================*/
				allWeights_pos+=projectorPlan.orientation_num*translation_num;
				HANDLE_ERROR(cudaStreamSynchronize(0));

				CUDA_GPU_TIC("diff2sMemCpCoarse");
				diff2s.cp_to_host();
				CUDA_GPU_TAC("diff2sMemCpCoarse");

				HANDLE_ERROR(cudaStreamSynchronize(0));
				CUDA_GPU_TOC();

				for (unsigned i = 0; i < projectorPlan.orientation_num; i ++)
				{
					unsigned iorientclass = projectorPlan.iorientclasses[i];
					for (unsigned j = 0; j < translation_num; j ++)
						DIRECT_A2D_ELEM(op.Mweight, ipart, iorientclass * translation_num + j) = diff2s[i * translation_num + j];
				}
			} // end if class significant
		} // end loop iclass
		op.min_diff2[ipart] = thrustGetMinVal(~allWeights, allWeights.size); // class
	} // end loop ipart
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF1);
#endif
}

void getAllSquaredDifferencesFine(unsigned exp_ipass,
		 	 	 	 	 	 	  OptimisationParamters &op,
		 	 	 	 	 	 	  SamplingParameters &sp,
		 	 	 	 	 	 	  MlOptimiser *baseMLO,
		 	 	 	 	 	 	  MlOptimiserCuda *cudaMLO,
		 	 	 	 	 	 	  std::vector<IndexedDataArray> &FinePassWeights,
		 	 	 	 	 	 	  std::vector<std::vector< IndexedDataArrayMask > > &FPCMasks,
		 	 	 	 	 	 	  std::vector<ProjectionParams> &FineProjectionData)
{
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2);
#endif
	CUDA_CPU_TIC("diff_pre_gpu");

	op.min_diff2.clear();
	op.min_diff2.resize(sp.nr_particles, 99.e99);
	CUDA_CPU_TIC("precalculateShiftedImagesCtfsAndInvSigma2s");
	std::vector<MultidimArray<Complex > > dummy;
	baseMLO->precalculateShiftedImagesCtfsAndInvSigma2s(false, op.my_ori_particle, sp.current_image_size, sp.current_oversampling,
			sp.itrans_min, sp.itrans_max, op.Fimgs, dummy, op.Fctfs, op.local_Fimgs_shifted, dummy,
			op.local_Fctfs, op.local_sqrtXi2, op.local_Minvsigma2s);
	CUDA_CPU_TOC("precalculateShiftedImagesCtfsAndInvSigma2s");
	MultidimArray<Complex > Fref;
	Fref.resize(op.local_Minvsigma2s[0]);

	unsigned image_size = op.local_Minvsigma2s[0].nzyxdim;

	CUDA_CPU_TOC("diff_pre_gpu");

	/*=======================================================================================
										  Particle Iteration
	=========================================================================================*/
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		// Reset size without de-allocating: we will append everything significant within
		// the current allocation and then re-allocate the then determined (smaller) volume
		unsigned long newDataSize(0);

		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		long int group_id = baseMLO->mydata.getGroupId(part_id);
		std::vector< double > oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
		/*====================================
				Generate Translations
		======================================*/

		CUDA_CPU_TIC("translation_1");

		CudaGlobalPtr<XFLOAT> Fimgs_real(image_size * sp.nr_trans * sp.nr_oversampled_trans, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_imag(image_size * sp.nr_trans * sp.nr_oversampled_trans, cudaMLO->allocator);

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
					XFLOAT real = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real
							- (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag;
					XFLOAT imag = (*(myAB + n)).real * (DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).imag
							+ (*(myAB + n)).imag *(DIRECT_MULTIDIM_ELEM(op.local_Fimgs_shifted[ipart], n)).real;

					//When on gpu, it makes more sense to ctf-correct translated images, rather than anti-ctf-correct ref-projections
					if (baseMLO->do_scale_correction)
					{
						//group_id = baseMLO->mydata.getGroupId(part_id);
						XFLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
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
		CUDA_CPU_TIC("kernel_init_1");

		CudaGlobalPtr<XFLOAT> gpuMinvsigma2(image_size, cudaMLO->allocator);
		gpuMinvsigma2.device_alloc();
		// Since we hijack Minvsigma to carry a bit more info into the GPU-kernel
		// we need to make a modified copy, since the global object shouldn't be
		// changed
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
		{
			gpuMinvsigma2[n] = *(op.local_Minvsigma2s[ipart].data + n );
		}

		if (baseMLO->do_ctf_correction && baseMLO->refs_are_ctf_corrected)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
			{
				gpuMinvsigma2[n] *= (DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n)*DIRECT_MULTIDIM_ELEM(op.local_Fctfs[ipart], n));
			}
		}

		// TODO :    + Assure accuracy with the implemented GPU-based ctf-scaling
		//           + Make setting of myscale robust between here and above.
		//  (scale_correction turns off by default with only one group: ml_optimiser-line 1067,
		//   meaning small-scale test will probably not catch this malfunctioning when/if it breaks.)
		if (baseMLO->do_scale_correction)
		{
			XFLOAT myscale = baseMLO->mymodel.scale_correction[group_id];
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(op.local_Fimgs_shifted[ipart])
			{
				gpuMinvsigma2[n] *= (myscale*myscale);
			}
		}

		CUDA_GPU_TIC("imagMemCp");
		gpuMinvsigma2.cp_to_device();
		Fimgs_real.put_on_device(translation_num * image_size);
		Fimgs_imag.put_on_device(translation_num * image_size);
		CUDA_GPU_TAC("imagMemCp");
		CUDA_CPU_TOC("kernel_init_1");
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			FPCMasks[ipart][exp_iclass].weightNum=0;

			if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) && (FineProjectionData[ipart].class_entries[exp_iclass] > 0) )
			{
				// use "slice" constructor with class-specific parameters to retrieve a temporary ProjectionParams with data for this class
				ProjectionParams thisClassProjectionData(	FineProjectionData[ipart],
															FineProjectionData[ipart].class_idx[exp_iclass],
															FineProjectionData[ipart].class_idx[exp_iclass]+FineProjectionData[ipart].class_entries[exp_iclass]);
				// since we retrieved the ProjectionParams for *the whole* class the orientation_num is also equal.

				thisClassProjectionData.orientation_num[0] = FineProjectionData[ipart].class_entries[exp_iclass];
				long unsigned orientation_num  = thisClassProjectionData.orientation_num[0];

				if(orientation_num==0)
					continue;

				// Local variables
				std::vector< double > oversampled_rot, oversampled_tilt, oversampled_psi;

				CUDA_CPU_TIC("generateEulerMatrices");
				CudaGlobalPtr<XFLOAT> eulers(9 * orientation_num, cudaMLO->allocator);

				generateEulerMatrices(
						baseMLO->mymodel.PPref[exp_iclass].padding_factor,
						thisClassProjectionData,
						&eulers[0],
						!IS_NOT_INV);

				CUDA_CPU_TOC("generateEulerMatrices");

				CUDA_GPU_TIC("eulersMemCp_1");
				eulers.device_alloc();
				eulers.cp_to_device();
				CUDA_GPU_TAC("eulersMemCp_1");

				/*===========================================
				   Determine significant comparison indices
				=============================================*/
				//      This section is annoying to test because
				//		it can't complete on first pass, since
				//		the significance has never been set

				CUDA_CPU_TIC("pair_list_1");
				long unsigned significant_num(0);
				long int nr_over_orient = baseMLO->sampling.oversamplingFactorOrientations(sp.current_oversampling);
				long int nr_over_trans = baseMLO->sampling.oversamplingFactorTranslations(sp.current_oversampling);

				// Prepare the mask of the weight-array for this class
				if (FPCMasks[ipart][exp_iclass].weightNum==0)
					FPCMasks[ipart][exp_iclass].firstPos = newDataSize;

				// Do more significance checks on translations and create jobDivision
				significant_num = makeJobsForDiff2Fine(	op,	sp,												// alot of different type inputs...
														orientation_num, translation_num,
														thisClassProjectionData,
														iover_transes, ihiddens,
														nr_over_orient, nr_over_trans, ipart,
														FinePassWeights[ipart],
														FPCMasks[ipart][exp_iclass]);                // ..and output into index-arrays mask

				// extend size by number of significants found this class
				newDataSize += significant_num;
				FPCMasks[ipart][exp_iclass].weightNum = significant_num;
				FPCMasks[ipart][exp_iclass].lastPos = FPCMasks[ipart][exp_iclass].firstPos + significant_num;

				if(significant_num==0)
					continue;

				// Use the constructed mask to construct a partial class-specific input
				IndexedDataArray thisClassFinePassWeights(FinePassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);

				CUDA_CPU_TOC("pair_list_1");
				CUDA_CPU_TIC("Diff2MakeKernel");
				CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
						cudaMLO->cudaProjectors[exp_iclass],
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						op.local_Minvsigma2s[0].xdim-1);
				CUDA_CPU_TOC("Diff2MakeKernel");
				CUDA_CPU_TIC("Diff2CALL");

				CUDA_GPU_TIC("IndexedArrayMemCp");
				thisClassFinePassWeights.weights.cp_to_device();
				thisClassFinePassWeights.rot_id.cp_to_device(); //FIXME this is not used
				thisClassFinePassWeights.rot_idx.cp_to_device();
				thisClassFinePassWeights.trans_idx.cp_to_device();
				FPCMasks[ipart][exp_iclass].jobOrigin.cp_to_device();
				FPCMasks[ipart][exp_iclass].jobExtent.cp_to_device();
				CUDA_GPU_TAC("IndexedArrayMemCp");


				CUDA_GPU_TIC("kernel_diff_proj");


				// Could be used to automate __ldg() fallback runtime within cuda_kernel_diff2.
				//				cudaDeviceProp dP;
				//				cudaGetDeviceProperties(&dP, 0);
				//				printf("-arch=sm_%d%d\n", dP.major, dP.minor);

				if ((baseMLO->iter == 1 && baseMLO->do_firstiter_cc) || baseMLO->do_always_cc) // do cross-correlation instead of diff
				{
					// FIXME  make _CC
					printf("Cross correlation is not supported yet.");
					exit(0);
				}

				runDiff2KernelFine(
						projKernel,
						~gpuMinvsigma2,
						~Fimgs_real,
						~Fimgs_imag,
						~eulers,
						~thisClassFinePassWeights.rot_id,
						~thisClassFinePassWeights.rot_idx,
						~thisClassFinePassWeights.trans_idx,
						~FPCMasks[ipart][exp_iclass].jobOrigin,
						~FPCMasks[ipart][exp_iclass].jobExtent,
						~thisClassFinePassWeights.weights,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						significant_num,
						image_size,
						ipart,
						FPCMasks[ipart][exp_iclass].jobOrigin.size
						);

				HANDLE_ERROR(cudaStreamSynchronize(0));
				CUDA_GPU_TOC();
				CUDA_CPU_TOC("Diff2CALL");

			} // end if class significant
		} // end loop iclass
		FinePassWeights[ipart].setDataSize( newDataSize );
		CUDA_CPU_TIC("collect_data_1");
		op.min_diff2[ipart] = std::min(op.min_diff2[ipart],(double)thrustGetMinVal(~FinePassWeights[ipart].weights, newDataSize));
		CUDA_CPU_TOC("collect_data_1");
	}// end loop ipart
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_DIFF2);
#endif
}


void convertAllSquaredDifferencesToWeights(unsigned exp_ipass,
											OptimisationParamters &op,
											SamplingParameters &sp,
											MlOptimiser *baseMLO,
											MlOptimiserCuda *cudaMLO,
											std::vector< IndexedDataArray> &PassWeights,
											std::vector< std::vector< IndexedDataArrayMask > > &FPCMasks) // FPCMasks = Fine-Pass Class-Masks
{
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
	{
		if (exp_ipass == 0) baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT1);
		else baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT2);
	}
#endif
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
			std::cerr << "the gpu inplementation cannot handle the new dense arrays but instead needs the old Mweight. Go maketh if thous needeth it." << std::endl;
			exit(1);
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
			long int sumRedSize=0;
			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
				sumRedSize+= (exp_ipass==0) ? sp.nr_dir*sp.nr_psi*sp.nr_oversampled_rot/SUM_BLOCK_SIZE : ceil((float)FPCMasks[ipart][exp_iclass].jobNum / (float)SUM_BLOCK_SIZE);

			CudaGlobalPtr<XFLOAT> thisparticle_sumweight(sumRedSize, cudaMLO->allocator);
			thisparticle_sumweight.host_alloc();
			thisparticle_sumweight.device_alloc();

			long int sumweight_pos=0;

			// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
			for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
			{
				// TODO Move to RANK LEVEL
				/*=========================================
						Fetch+generate Orientation data
				===========================================*/
				CudaGlobalPtr<XFLOAT>  pdf_orientation(sp.nr_dir * sp.nr_psi, cudaMLO->allocator);
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
				CudaGlobalPtr<XFLOAT>  pdf_offset(sp.nr_trans, cudaMLO->allocator);
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

				CUDA_CPU_TIC("sumweight1");
				CUDA_GPU_TIC("sumweightMemCp1");

				//std::cerr << "summing weights on GPU... " << std::endl;
				pdf_orientation.put_on_device();
				pdf_offset.put_on_device();

				CUDA_GPU_TAC("sumweightMemCp1");
				CUDA_GPU_TOC();

				long int block_num;

				if(exp_ipass==0)  //use Mweight for now - FIXME use PassWeights.weights (ignore indexArrays)
				{
					CudaGlobalPtr<XFLOAT>  Mweight( &(op.Mweight.data[(ipart)*(op.Mweight).xdim+
					                                  exp_iclass * sp.nr_dir * sp.nr_psi * sp.nr_trans]),
													  sp.nr_dir * sp.nr_psi * sp.nr_trans, cudaMLO->allocator);
					Mweight.put_on_device();
					block_num = sp.nr_dir*sp.nr_psi*sp.nr_oversampled_rot/SUM_BLOCK_SIZE;

					dim3 block_dim(block_num);
					CUDA_GPU_TIC("cuda_kernel_sumweight");
					cuda_kernel_sumweightCoarse<<<block_dim,SUM_BLOCK_SIZE>>>(	~pdf_orientation,
																			    ~pdf_offset,
																			    ~Mweight,
																			    ~thisparticle_sumweight,
																			    (XFLOAT)op.min_diff2[ipart],
																			    sp.nr_oversampled_rot,
																			    sp.nr_oversampled_trans,
																			    sp.nr_trans,
																			    sumweight_pos);
					CUDA_GPU_TAC("cuda_kernel_sumweight");
					CUDA_GPU_TIC("sumweightMemCp2");
					Mweight.cp_to_host();  //FIXME remove when mapping is eliminated
					HANDLE_ERROR(cudaStreamSynchronize(0));
					CUDA_GPU_TAC("sumweightMemCp2");
					sumweight_pos+=block_num;
				}
				else if ((baseMLO->mymodel.pdf_class[exp_iclass] > 0.) && (FPCMasks[ipart][exp_iclass].weightNum > 0) )
				{
					// Use the constructed mask to build a partial (class-specific) input
					// (until now, PassWeights has been an empty placeholder. We now create class-paritals pointing at it, and start to fill it with stuff)
					IndexedDataArray thisClassPassWeights(PassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);

					block_num = ceil((float)FPCMasks[ipart][exp_iclass].jobNum / (float)SUM_BLOCK_SIZE); //thisClassPassWeights.rot_idx.size / SUM_BLOCK_SIZE;
					dim3 block_dim(block_num);
					thisClassPassWeights.weights.cp_to_host();
					CUDA_GPU_TIC("cuda_kernel_sumweight");
					cuda_kernel_sumweightFine<<<block_dim,SUM_BLOCK_SIZE>>>(	~pdf_orientation,
																			    ~pdf_offset,
																			    ~thisClassPassWeights.weights,
																			    ~thisparticle_sumweight,
																			    (XFLOAT)op.min_diff2[ipart],
																			    sp.nr_oversampled_rot,
																			    sp.nr_oversampled_trans,
																			    ~thisClassPassWeights.rot_id,
																			    ~thisClassPassWeights.trans_idx,
																				~FPCMasks[ipart][exp_iclass].jobOrigin,
																			 	~FPCMasks[ipart][exp_iclass].jobExtent,
																			 	FPCMasks[ipart][exp_iclass].jobNum,
																			 	sumweight_pos);
					CUDA_GPU_TAC("cuda_kernel_sumweight");

					CUDA_GPU_TIC("sumweightMemCp2");
					thisparticle_sumweight.cp_to_host();
					thisClassPassWeights.weights.cp_to_host();  //FIXME remove when mapping is eliminated - NOTE ALOT OF MWEIGHT-DEPS  BELOW
					HANDLE_ERROR(cudaStreamSynchronize(0));
					CUDA_GPU_TAC("sumweightMemCp2");
					sumweight_pos+=block_num;
				}
				CUDA_CPU_TOC("sumweight1");
			} // end loop exp_iclass
			thrust::device_ptr<XFLOAT> dp = thrust::device_pointer_cast(~thisparticle_sumweight);
			exp_thisparticle_sumweight += thrust::reduce(dp, dp + sumweight_pos);
		} // end if iter==1

		//Store parameters for this particle
		op.sum_weight[ipart] = exp_thisparticle_sumweight;
		//std::cerr << "  sumweight =  " << exp_thisparticle_sumweight << std::endl;

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

	CUDA_CPU_TIC("convertPostKernel");
	// Now, for each particle,  find the exp_significant_weight that encompasses adaptive_fraction of op.sum_weight
	op.significant_weight.clear();
	op.significant_weight.resize(sp.nr_particles, 0.);
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];

		MultidimArray<XFLOAT> sorted_weight;
		CudaGlobalPtr<XFLOAT> sorted_weight_new(cudaMLO->allocator);  // make new sorted weights
		double frac_weight = 0.;
		double my_significant_weight;
		long int my_nr_significant_coarse_samples = 0;
		long int np = 0;
		if (exp_ipass!=0)
		{
			CUDA_CPU_TIC("sort");
			sorted_weight_new.size = PassWeights[ipart].weights.size;
			sorted_weight_new.host_alloc();
			sorted_weight_new.d_ptr = PassWeights[ipart].weights.d_ptr;			    // set pointer to weights
			sorted_weight_new.cp_to_host();							// make host-copy
			sorted_weight_new.d_do_free = false;

			thrust::sort(sorted_weight_new.h_ptr, sorted_weight_new.h_ptr + sorted_weight_new.size );
			CUDA_CPU_TOC("sort");
			for (long int i=sorted_weight_new.size-1; i>=0; i--)
			{
//				if (exp_ipass==0) my_nr_significant_coarse_samples++;
					my_significant_weight = sorted_weight_new[i];
				//std::cerr << "thisweight = " << my_significant_weight << std::endl;
				frac_weight += my_significant_weight;
				if (frac_weight > baseMLO->adaptive_fraction * op.sum_weight[ipart])
					break;
			}
		}
		else
		{
			// Get the relevant row for this particle
			CUDA_CPU_TIC("getRow");
			op.Mweight.getRow(ipart, sorted_weight);
			CUDA_CPU_TOC("getRow");

			// Only select non-zero probabilities to speed up sorting // TODO Remove when mapping is eliminated
			CUDA_CPU_TIC("nonZero");
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sorted_weight)
			{
				if (DIRECT_MULTIDIM_ELEM(sorted_weight, n) > 0.)
				{
//					std::cerr << "weights = " << DIRECT_A2D_ELEM(op.Mweight, ipart, n) << std::endl;
					DIRECT_MULTIDIM_ELEM(sorted_weight, np) = DIRECT_MULTIDIM_ELEM(sorted_weight, n);
					np++;
				}
			}
			sorted_weight.resize(np);
			CUDA_CPU_TOC("nonZero");

			// Sort from low to high values
			CUDA_CPU_TIC("sort");
//			std::cerr << "sort on " << sorted_weight.xdim << " which should have np = " << np << std::endl;
#if  defined(USE_THRUST) && !defined(CUDA_DOUBLE_PRECISION) // Thrust seems incredibly slow in debug build this is clearly a FIXME
			thrust::sort(sorted_weight.data, sorted_weight.data + np );
#else
			sorted_weight.sort();
#endif
			CUDA_CPU_TOC("sort");

			for (long int i = XSIZE(sorted_weight) - 1; i >= 0; i--)
			{
				if (exp_ipass==0) my_nr_significant_coarse_samples++;
				my_significant_weight = DIRECT_A1D_ELEM(sorted_weight, i);
				//std::cerr << "thisweight = " << my_significant_weight << std::endl;
				frac_weight += my_significant_weight;
				if (frac_weight > baseMLO->adaptive_fraction * op.sum_weight[ipart])
					break;
			}
		}



		if (exp_ipass==0 && my_nr_significant_coarse_samples == 0)
		{
			std::cerr << " ipart= " << ipart << " adaptive_fraction= " << baseMLO->adaptive_fraction << std::endl;
			std::cerr << " frac-weight= " << frac_weight << std::endl;
			std::cerr << " op.sum_weight[ipart]= " << op.sum_weight[ipart] << std::endl;
			Image<XFLOAT> It;
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
		//std::cerr << "@sort op.significant_weight[ipart]= " << (XFLOAT)op.significant_weight[ipart] << std::endl;

	} // end loop ipart
	CUDA_CPU_TOC("convertPostKernel");
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
	{
		if (exp_ipass == 0) baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT1);
		else baseMLO->timer.tic(baseMLO->TIMING_ESP_WEIGHT2);
	}
#endif
}

void storeWeightedSums(OptimisationParamters &op, SamplingParameters &sp,
						MlOptimiser *baseMLO,
						MlOptimiserCuda *cudaMLO,
						std::vector<IndexedDataArray> &FinePassWeights,
						std::vector<ProjectionParams> &ProjectionData,
						std::vector<std::vector<IndexedDataArrayMask> > FPCMasks)  // FIXME? by ref?
{
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.tic(baseMLO->TIMING_ESP_WSUM);
#endif
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
	bool have_warned_small_scale = false;

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
	CUDA_CPU_TIC("maximization");

	/*=======================================================================================
	                                   MAXIMIZATION
	=======================================================================================*/

	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		long int part_id = baseMLO->mydata.ori_particles[op.my_ori_particle].particles_id[ipart];
		int group_id = baseMLO->mydata.getGroupId(part_id);

		/*======================================================
		                     TRANSLATIONS
		======================================================*/

		CUDA_CPU_TIC("translation_2");

		CudaGlobalPtr<XFLOAT> Fimgs_real(image_size * sp.nr_trans * sp.nr_oversampled_trans, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_imag(Fimgs_real.size, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_nomask_real(Fimgs_real.size, cudaMLO->allocator);
		CudaGlobalPtr<XFLOAT> Fimgs_nomask_imag(Fimgs_real.size, cudaMLO->allocator);

		std::vector< long unsigned > iover_transes, itranses, ihiddens;

		long unsigned translation_num = imageTranslation(
				Fimgs_real.h_ptr,
				Fimgs_imag.h_ptr,
				Fimgs_nomask_real.h_ptr,
				Fimgs_nomask_imag.h_ptr,
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


		Fimgs_real.put_on_device(translation_num * image_size);
		Fimgs_imag.put_on_device(Fimgs_real.size);
		Fimgs_nomask_real.put_on_device(Fimgs_real.size);
		Fimgs_nomask_imag.put_on_device(Fimgs_real.size);

		CUDA_CPU_TOC("translation_2");


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

		CudaGlobalPtr<XFLOAT> ctfs(image_size, cudaMLO->allocator); //TODO Same size for all iparts, should be allocated once

		if (baseMLO->do_ctf_correction)
		{
			for (unsigned i = 0; i < image_size; i++)
				ctfs[i] = (XFLOAT) op.local_Fctfs[ipart].data[i] * part_scale;
		}
		else //TODO should be handled by memset
			for (unsigned i = 0; i < image_size; i++)
				ctfs[i] = part_scale;

		ctfs.put_on_device();

		/*======================================================
		                       MINVSIGMA
		======================================================*/

		CudaGlobalPtr<XFLOAT> Minvsigma2s(image_size, cudaMLO->allocator); //TODO Same size for all iparts, should be allocated once

		if (baseMLO->do_map)
			for (unsigned i = 0; i < image_size; i++)
				Minvsigma2s[i] = op.local_Minvsigma2s[ipart].data[i];
		else
			for (unsigned i = 0; i < image_size; i++)
				Minvsigma2s[i] = 1;

		Minvsigma2s.put_on_device();

		/*======================================================
		                      CLASS LOOP
		======================================================*/


		// Loop from iclass_min to iclass_max to deal with seed generation in first iteration
		for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
		{
			if((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[ipart].class_entries[exp_iclass] == 0))
				continue;

			// Use the constructed mask to construct a partial class-specific input
			IndexedDataArray thisClassFinePassWeights(FinePassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);

			CUDA_CPU_TIC("thisClassProjectionSetupCoarse");
			// use "slice" constructor with class-specific parameters to retrieve a temporary ProjectionParams with data for this class
			ProjectionParams thisClassProjectionData(	ProjectionData[ipart],
														ProjectionData[ipart].class_idx[exp_iclass],
														ProjectionData[ipart].class_idx[exp_iclass]+ProjectionData[ipart].class_entries[exp_iclass]);

			thisClassProjectionData.orientation_num[0] = ProjectionData[ipart].orientation_num[exp_iclass];
			CUDA_CPU_TOC("thisClassProjectionSetupCoarse");

			unsigned proj_div_nr = ceil((float)thisClassProjectionData.orientation_num[0] / (float)proj_div_max_count);
			unsigned long idxArrPos_start(0),idxArrPos_end(0);

			/// Now that reference projection has been made loop over all particles inside this ori_particle
			for (int iproj_div = 0; iproj_div < proj_div_nr; iproj_div++)
			{
				CUDA_CPU_TIC("BP-ProjectionDivision");
				unsigned long proj_div_start(proj_div_max_count * iproj_div), proj_div_end;

				if (iproj_div < proj_div_nr - 1)
					proj_div_end = proj_div_start + proj_div_max_count;
				else
					proj_div_end = thisClassProjectionData.orientation_num[0];

				long unsigned orientation_num(proj_div_end - proj_div_start);

				// use "slice" constructor to slice out betwen start and stop in the specified region of thisClassProjectionData
				ProjectionParams ProjectionData_projdiv(thisClassProjectionData,proj_div_start,proj_div_end);
				idxArrPos_start=idxArrPos_end;
				while((thisClassFinePassWeights.rot_idx[idxArrPos_end]-thisClassFinePassWeights.rot_idx[idxArrPos_start]<orientation_num) && (idxArrPos_end < thisClassFinePassWeights.rot_idx.size))
					idxArrPos_end++;

				CUDA_CPU_TOC("BP-ProjectionDivision");

				/*======================================================
				                    PROJECTIONS
				======================================================*/

				cudaStream_t currentBPStream = cudaMLO->cudaBackprojectors[exp_iclass].getStream();

				CUDA_CPU_TIC("generateEulerMatrices");
				CudaGlobalPtr<XFLOAT> eulers(9 * orientation_num, currentBPStream, cudaMLO->allocator);

				generateEulerMatrices(
						baseMLO->mymodel.PPref[exp_iclass].padding_factor,
						ProjectionData_projdiv,
						&eulers[0],
						!IS_NOT_INV);

				eulers.put_on_device();

				CUDA_CPU_TOC("generateEulerMatrices");

				CudaGlobalPtr<XFLOAT> wavgs_real(orientation_num * image_size, currentBPStream, cudaMLO->allocator);
				wavgs_real.device_alloc();
				wavgs_real.device_init(0);
				CudaGlobalPtr<XFLOAT> wavgs_imag(orientation_num * image_size, currentBPStream, cudaMLO->allocator);
				wavgs_imag.device_alloc();
				wavgs_imag.device_init(0);
				CudaGlobalPtr<XFLOAT> Fweights(orientation_num * image_size, currentBPStream, cudaMLO->allocator);
				Fweights.device_alloc();
				Fweights.device_init(0);


				/*======================================================
				                     MAP WEIGHTS
				======================================================*/

				CudaGlobalPtr<XFLOAT> sorted_weights(orientation_num * translation_num, currentBPStream, cudaMLO->allocator);

				mapWeights(
						proj_div_start,
						&sorted_weights[0],
						orientation_num,
						idxArrPos_start,
						idxArrPos_end,
						translation_num,
						&thisClassFinePassWeights.weights[0],
						&thisClassFinePassWeights.rot_idx[0],
						&thisClassFinePassWeights.trans_idx[0],
						baseMLO->sampling,
						ipart,
						iover_transes,
						ihiddens,
						ProjectionData_projdiv.iorientclasses,
						ProjectionData_projdiv.iover_rots,
						op.Mweight,
						sp.current_oversampling,
						sp.nr_trans);

				sorted_weights.put_on_device();

				/*======================================================
				                     KERNEL CALL
				======================================================*/

				CudaGlobalPtr<XFLOAT> wdiff2s_parts(orientation_num * image_size, currentBPStream, cudaMLO->allocator);
				wdiff2s_parts.device_alloc();

				CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
						cudaMLO->cudaProjectors[exp_iclass],
						op.local_Minvsigma2s[0].xdim,
						op.local_Minvsigma2s[0].ydim,
						op.local_Minvsigma2s[0].xdim-1);

				runWavgKernel(
						projKernel,
						~eulers,
						~Fimgs_real,
						~Fimgs_imag,
						~Fimgs_nomask_real,
						~Fimgs_nomask_imag,
						~sorted_weights,
						~ctfs,
						~Minvsigma2s,
						~wdiff2s_parts,
						~wavgs_real,
						~wavgs_imag,
						~Fweights,
						op,
						baseMLO,
						orientation_num,
						translation_num,
						image_size,
						ipart,
						group_id,
						exp_iclass,
						currentBPStream);

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
					cuda_kernel_reduce_wdiff2s<<<block_dim_wd,BLOCK_SIZE,0,currentBPStream>>>(
							~wdiff2s_parts,
							orientation_num,
							image_size,
							k);
				}

				CUDA_GPU_TAC("cuda_kernels_reduce_wdiff2s");

				wdiff2s_parts.size = image_size; //temporarily set the size to the single image we have now reduced, to not copy more than necessary
				wdiff2s_parts.cp_to_host();
				wdiff2s_parts.size = orientation_num * image_size;

				HANDLE_ERROR(cudaStreamSynchronize(currentBPStream));

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

				/*======================================================
				                    BACKPROJECTION
				======================================================*/

				CudaGlobalPtr<XFLOAT> bp_eulers(9*orientation_num, currentBPStream, cudaMLO->allocator);

				XFLOAT padding_factor = baseMLO->wsum_model.BPref[exp_iclass].padding_factor;

				generateEulerMatrices(
						1/padding_factor, //Why squared scale factor is given in backprojection
						ProjectionData_projdiv,
						&bp_eulers[0],
						IS_NOT_INV);

#ifdef TIMING
				if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
					baseMLO->timer.tic(baseMLO->TIMING_WSUM_BACKPROJ);
#endif

				if (cudaMLO->refIs3D)
				{
					bp_eulers.device_alloc();
					bp_eulers.cp_to_device();

					CUDA_GPU_TIC("cuda_kernels_backproject");
					cudaMLO->cudaBackprojectors[exp_iclass].backproject(
						~wavgs_real,
						~wavgs_imag,
						~Fweights,
						~bp_eulers,
							op.local_Minvsigma2s[0].xdim,
							op.local_Minvsigma2s[0].ydim,
							orientation_num);

					HANDLE_ERROR(cudaStreamSynchronize(currentBPStream));
					CUDA_GPU_TAC("cuda_kernels_backproject");
				}
				else
				{
					CUDA_CPU_TIC("cpu_backproject");

					wavgs_real.cp_to_host();
					wavgs_imag.cp_to_host();
					Fweights.cp_to_host();

					MultidimArray<Complex > Fimg;
					MultidimArray<double > Fweight;
					Fimg.resize(op.Fimgs[0]);
					Fweight.resize(op.Fimgs[0]);
					Matrix2D<double> A(3,3);

					int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
					pthread_mutex_lock(&global_mutex2[my_mutex]);
					for (int i = 0; i < orientation_num; i ++)
					{
						for (int j = 0; j < image_size; j ++)
						{
							Fimg.data[j].real = (double) wavgs_real[i * image_size + j];
							Fimg.data[j].imag = (double) wavgs_imag[i * image_size + j];
							Fweight.data[j]   = (double) Fweights[i * image_size + j];
						}

						for (int j = 0; j < 9; j ++)
							A.mdata[j] = bp_eulers[i * 9 + j];

						baseMLO->wsum_model.BPref[exp_iclass].backrotate2D(Fimg, A, IS_NOT_INV, &Fweight);
					}
					pthread_mutex_unlock(&global_mutex2[my_mutex]);
					CUDA_CPU_TOC("cpu_backproject");
				}

#ifdef TIMING
				if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
					baseMLO->timer.toc(baseMLO->TIMING_WSUM_BACKPROJ);
#endif
			} // end loop proj_div
		} // end loop iclass
	} // end loop ipart
	CUDA_CPU_TOC("maximization");
	CUDA_CPU_TIC("collect_data_2");
	for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
	{
		/*=======================================================================================
		                           COLLECT 2 AND SET METADATA
		=======================================================================================*/
		for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
		{
			if ((baseMLO->mymodel.pdf_class[exp_iclass] == 0.) || (ProjectionData[ipart].class_entries[exp_iclass] == 0) )
				continue;
			// Use the constructed mask to construct a partial class-specific input
			IndexedDataArray thisClassFinePassWeights(FinePassWeights[ipart],FPCMasks[ipart][exp_iclass], cudaMLO->allocator);

			CUDA_CPU_TIC("thisClassProjectionSetupCoarse");
			// use "slice" constructor with class-specific parameters to retrieve a temporary ProjectionParams with data for this class
			ProjectionParams thisClassProjectionData(	ProjectionData[ipart],
														ProjectionData[ipart].class_idx[exp_iclass],
														ProjectionData[ipart].class_idx[exp_iclass]+ProjectionData[ipart].class_entries[exp_iclass]);

			thisClassProjectionData.orientation_num[0] = ProjectionData[ipart].orientation_num[exp_iclass];
			CUDA_CPU_TOC("thisClassProjectionSetupCoarse");
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

			/*======================================================
								COLLECT 2
			======================================================*/

			CUDA_CPU_TIC("collect_data_2_pre_kernel");
			//TODO should be replaced with loop over pairs of projections and translations (like in the getAllSquaredDifferences-function)

			int oversamples = sp.nr_oversampled_trans * sp.nr_oversampled_rot;
//				CudaGlobalPtr<XFLOAT >  Mweight( &(op.Mweight.data[(ipart)*(op.Mweight).xdim]),
//												sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples);
			int nr_transes = sp.nr_trans*sp.nr_oversampled_trans;
				CudaGlobalPtr<XFLOAT>     oo_otrans_x(nr_transes, cudaMLO->allocator); // old_offset_oversampled_trans_x
				CudaGlobalPtr<XFLOAT>     oo_otrans_y(nr_transes, cudaMLO->allocator);
				CudaGlobalPtr<XFLOAT> myp_oo_otrans_x2y2z2(nr_transes, cudaMLO->allocator); // my_prior_old_offs....x^2*y^2*z^2

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
			// Re-define the job-partition of the indexedArray of weights so that the collect-kernel can work with it.
			int block_num = makeJobsForCollect(thisClassFinePassWeights, FPCMasks[ipart][exp_iclass]);

				oo_otrans_x.put_on_device();
				oo_otrans_y.put_on_device();
				myp_oo_otrans_x2y2z2.put_on_device();

//				std::cerr << "block_num = " << block_num << std::endl;
				CudaGlobalPtr<XFLOAT>                      p_weights(block_num, cudaMLO->allocator);
				CudaGlobalPtr<XFLOAT> p_thr_wsum_prior_offsetx_class(block_num, cudaMLO->allocator);
				CudaGlobalPtr<XFLOAT> p_thr_wsum_prior_offsety_class(block_num, cudaMLO->allocator);
				CudaGlobalPtr<XFLOAT>       p_thr_wsum_sigma2_offset(block_num, cudaMLO->allocator);

				p_weights.device_alloc();
				p_thr_wsum_prior_offsetx_class.device_alloc();
				p_thr_wsum_prior_offsety_class.device_alloc();
				p_thr_wsum_sigma2_offset.device_alloc();
			CUDA_CPU_TOC("collect_data_2_pre_kernel");

			CUDA_GPU_TIC("collect2-kernel");
			dim3 grid_dim_collect2 = splitCudaBlocks(block_num,false);
			cuda_kernel_collect2jobs<<<grid_dim_collect2,SUM_BLOCK_SIZE>>>(
						~oo_otrans_x,          // otrans-size -> make const
						~oo_otrans_y,          // otrans-size -> make const
						~myp_oo_otrans_x2y2z2, // otrans-size -> make const
						~thisClassFinePassWeights.weights,
					(XFLOAT)op.significant_weight[ipart],
					(XFLOAT)op.sum_weight[ipart],
					sp.nr_trans,
					sp.nr_oversampled_trans,
					sp.nr_oversampled_rot,
					oversamples,
					(baseMLO->do_skip_align || baseMLO->do_skip_rotate ),
						~p_weights,
						~p_thr_wsum_prior_offsetx_class,
						~p_thr_wsum_prior_offsety_class,
						~p_thr_wsum_sigma2_offset,
					~thisClassFinePassWeights.rot_idx,
					~thisClassFinePassWeights.trans_idx,
					~FPCMasks[ipart][exp_iclass].jobOrigin,
					~FPCMasks[ipart][exp_iclass].jobExtent
						);
			CUDA_GPU_TAC("collect2-kernel");

			CUDA_CPU_TIC("collect_data_2_post_kernel");
			// TODO further reduce the below 4 arrays while data is still on gpu
				p_weights.cp_to_host();
				p_thr_wsum_prior_offsetx_class.cp_to_host();
				p_thr_wsum_prior_offsety_class.cp_to_host();
				p_thr_wsum_sigma2_offset.cp_to_host();

			HANDLE_ERROR(cudaStreamSynchronize(0));
			CUDA_GPU_TOC();

			thr_wsum_sigma2_offset = 0.0;
			int iorient = 0;
			for (long int n = 0; n < block_num; n++)
			{
				iorient= thisClassFinePassWeights.rot_id[FPCMasks[ipart][exp_iclass].jobOrigin[n]];
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
						thr_wsum_prior_offsetx_class[exp_iclass] += p_thr_wsum_prior_offsetx_class[n];
						thr_wsum_prior_offsety_class[exp_iclass] += p_thr_wsum_prior_offsety_class[n];
					}
				}
			}
				p_weights.free();
				p_thr_wsum_sigma2_offset.free();
				p_thr_wsum_prior_offsetx_class.free();
				p_thr_wsum_prior_offsety_class.free();
				oo_otrans_y.free();
				oo_otrans_x.free();
				myp_oo_otrans_x2y2z2.free();
			CUDA_CPU_TOC("collect_data_2_post_kernel");
		}
	} // end loop iclass

	std::vector< double> oversampled_rot, oversampled_tilt, oversampled_psi;
	for (long int ipart = 0; ipart < sp.nr_particles; ipart++)
	{
		CUDA_CPU_TIC("setMetadata");
		/*======================================================
							SET METADATA
		======================================================*/

		//Get index of max element using GPU-tool thrust
		Indices max_index;
		thrust::device_ptr<XFLOAT> dp = thrust::device_pointer_cast(~FinePassWeights[ipart].weights);
		thrust::device_ptr<XFLOAT> pos = thrust::max_element(dp, dp + FinePassWeights[ipart].weights.size);
		unsigned int pos_idx = thrust::distance(dp, pos);

		XFLOAT max_val;
		HANDLE_ERROR(cudaMemcpyAsync(&max_val, &FinePassWeights[ipart].weights.d_ptr[pos_idx], sizeof(XFLOAT), cudaMemcpyDeviceToHost, 0));

		if(max_val>op.max_weight[ipart])
		{
			op.max_weight[ipart] = max_val;
			max_index.fineIdx =FinePassWeights[ipart].ihidden_overs[pos_idx];
			//std::cerr << "max val = " << op.max_weight[ipart] << std::endl;
			//std::cerr << "max index = " << max_index.fineIdx << std::endl;
			max_index.fineIndexToFineIndices(sp); // set partial indices corresponding to the found max_index, to be used below

			baseMLO->sampling.getTranslations(max_index.itrans, baseMLO->adaptive_oversampling,
					oversampled_translations_x, oversampled_translations_y, oversampled_translations_z);

			//TODO We already have rot, tilt and psi don't calculated them again
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
		}
		CUDA_CPU_TOC("setMetadata");
	}
	CUDA_CPU_TOC("collect_data_2");
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
#ifdef TIMING
	if (op.my_ori_particle == baseMLO->exp_my_first_ori_particle)
		baseMLO->timer.toc(baseMLO->TIMING_ESP_WSUM);
#endif
}

MlOptimiserCuda::MlOptimiserCuda(MlOptimiser *baseMLOptimiser, int dev_id) : baseMLO(baseMLOptimiser)
{
	unsigned nr_classes = baseMLOptimiser->mymodel.nr_classes;

	/*======================================================
					DEVICE MEM OBJ SETUP
	======================================================*/
	device_id=dev_id;
	int devCount;
	cudaGetDeviceCount(&devCount);
	if(dev_id>=devCount)
	{
		std::cerr << " using device_id=" << dev_id << " (device no. " << dev_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
		REPORT_ERROR("ERROR: Too many MPI threads using GPUs");
	}
	else
	{
		cudaSetDevice(dev_id);
	}

	allocator = new CudaCustomAllocator(1024*1024*512);
	
	/*======================================================
	   PROJECTOR, PROJECTOR PLAN AND BACKPROJECTOR SETUP
	======================================================*/

	refIs3D = baseMLO->mymodel.ref_dim == 3;

	cudaProjectors.resize(nr_classes);
	cudaBackprojectors.resize(nr_classes);

	//Can we pre-generate projector plan and corresponding euler matrices for all particles
	if (baseMLO->do_skip_align || baseMLO->do_skip_rotate || baseMLO->do_auto_refine || baseMLO->mymodel.orientational_prior_mode != NOPRIOR)
		generateProjectionPlanOnTheFly = true;
	else
	{
		generateProjectionPlanOnTheFly = false;
		cudaCoarseProjectionPlans.resize(nr_classes);
	}

	//Loop over classes
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		cudaProjectors[iclass].setMdlDim(
				baseMLO->mymodel.PPref[iclass].data.xdim,
				baseMLO->mymodel.PPref[iclass].data.ydim,
				baseMLO->mymodel.PPref[iclass].data.zdim,
				baseMLO->mymodel.PPref[iclass].data.yinit,
				baseMLO->mymodel.PPref[iclass].data.zinit,
				baseMLO->mymodel.PPref[iclass].r_max,
				baseMLO->mymodel.PPref[iclass].padding_factor);

		cudaProjectors[iclass].setMdlData(baseMLO->mymodel.PPref[iclass].data.data);

		cudaBackprojectors[iclass].setMdlDim(
				baseMLO->wsum_model.BPref[iclass].data.xdim,
				baseMLO->wsum_model.BPref[iclass].data.ydim,
				baseMLO->wsum_model.BPref[iclass].data.zdim,
				baseMLO->wsum_model.BPref[iclass].data.yinit,
				baseMLO->wsum_model.BPref[iclass].data.zinit,
				baseMLO->wsum_model.BPref[iclass].r_max,
				baseMLO->wsum_model.BPref[iclass].padding_factor);

		cudaBackprojectors[iclass].initMdl();

		//If doing predefined projector plan at all and is this class significant
		if (!generateProjectionPlanOnTheFly && baseMLO->mymodel.pdf_class[iclass] > 0.)
		{
			std::vector<int> exp_pointer_dir_nonzeroprior;
			std::vector<int> exp_pointer_psi_nonzeroprior;
			std::vector<double> exp_directions_prior;
			std::vector<double> exp_psi_prior;

			long unsigned itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;
			long unsigned nr_idir = baseMLO->sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
			long unsigned nr_ipsi = baseMLO->sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior );

			cudaCoarseProjectionPlans[iclass].setup(
					baseMLO->sampling,
					exp_directions_prior,
					exp_psi_prior,
					exp_pointer_dir_nonzeroprior,
					exp_pointer_psi_nonzeroprior,
					NULL, //Mcoarse_significant
					baseMLO->mymodel.pdf_class,
					baseMLO->mymodel.pdf_direction,
					nr_idir,
					nr_ipsi,
					0, //idir_min
					nr_idir - 1, //idir_max
					0, //ipsi_min
					nr_ipsi - 1, //ipsi_max
					0, //itrans_min
					itrans_max,
					0, //current_oversampling
					1, //nr_oversampled_rot
					iclass,
					true, //coarse
					!IS_NOT_INV,
					baseMLO->do_skip_align,
					baseMLO->do_skip_rotate,
					baseMLO->mymodel.orientational_prior_mode
					);
		}
	}
};

void MlOptimiserCuda::doThreadExpectationSomeParticles(unsigned thread_id)
{

//	CUDA_CPU_TOC("interParticle");
	CUDA_CPU_TIC("oneTask");
	cudaSetDevice(device_id);
	//std::cerr << " calling on device " << device_id << std::endl;
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
			CUDA_CPU_TIC("getFourierTransformsAndCtfs");
			baseMLO->getFourierTransformsAndCtfs(my_ori_particle, op.metadata_offset, op.Fimgs, op.Fimgs_nomask, op.Fctfs,
					op.old_offset, op.prior, op.power_imgs, op.highres_Xi2_imgs,
					op.pointer_dir_nonzeroprior, op.pointer_psi_nonzeroprior, op.directions_prior, op.psi_prior);
			CUDA_CPU_TOC("getFourierTransformsAndCtfs");
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

			/// -- This is a iframe-indexed vector, each entry of which is a dense data-array. These are replacements to using
			//    Mweight in the sparse (Fine-sampled) pass, coarse is unused but created empty input for convert ( FIXME )
			std::vector <IndexedDataArray> CoarsePassWeights(1, allocator) ,FinePassWeights(sp.nr_particles, allocator);
			// -- This is a iframe-indexed vector, each entry of which is a class-indexed vector of masks, one for each
			//    class in FinePassWeights
			std::vector < std::vector <IndexedDataArrayMask> > FinePassClassMasks(sp.nr_particles, std::vector <IndexedDataArrayMask>(baseMLO->mymodel.nr_classes, allocator));
			// -- This is a iframe-indexed vector, each entry of which is parameters used in the projection-operations *after* the
			//    coarse pass, declared here to keep scope to storeWS
			std::vector < ProjectionParams > FineProjectionData(sp.nr_particles, baseMLO->mymodel.nr_classes);

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

				if (ipass == 0)
				{
					std::vector< CudaProjectorPlan > coarseProjectionPlans;

					 //If particle specific sampling plan required
					if (generateProjectionPlanOnTheFly)
					{
						CUDA_CPU_TIC("generateProjectionSetupCoarse");

						coarseProjectionPlans.clear();
						coarseProjectionPlans.resize(baseMLO->mymodel.nr_classes);

						for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
						{
							if (baseMLO->mymodel.pdf_class[exp_iclass] > 0.)
							{
								coarseProjectionPlans[exp_iclass].setup(
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
										exp_iclass,
										true, //coarse
										!IS_NOT_INV,
										baseMLO->do_skip_align,
										baseMLO->do_skip_rotate,
										baseMLO->mymodel.orientational_prior_mode
										);
							}
						}
						CUDA_CPU_TOC("generateProjectionSetupCoarse");
					}
					else //Otherwise use precalculated plan
						coarseProjectionPlans = cudaCoarseProjectionPlans;

					CUDA_CPU_TIC("getAllSquaredDifferencesCoarse");
					getAllSquaredDifferencesCoarse(ipass, op, sp, baseMLO, this, coarseProjectionPlans);
					CUDA_CPU_TOC("getAllSquaredDifferencesCoarse");
					CUDA_CPU_TIC("convertAllSquaredDifferencesToWeightsCoarse");
					convertAllSquaredDifferencesToWeights(ipass, op, sp, baseMLO, this, CoarsePassWeights, FinePassClassMasks);
					CUDA_CPU_TOC("convertAllSquaredDifferencesToWeightsCoarse");
				}
				else
				{
//					// -- go through all classes and generate projectionsetups for all classes - to be used in getASDF and storeWS below --
//					// the reason to do this globally is subtle - we want the orientation_num of all classes to estimate a largest possible
//					// weight-array, which would be insanely much larger than necessary if we had to assume the worst.
					for (long int iframe = 0; iframe < sp.nr_particles; iframe++)
					{
						FineProjectionData[iframe].orientationNumAllClasses = 0;
						for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
						{
							if(exp_iclass>0)
								FineProjectionData[iframe].class_idx[exp_iclass]=FineProjectionData[iframe].rots.size();
							FineProjectionData[iframe].class_entries[exp_iclass]=0;

							CUDA_CPU_TIC("generateProjectionSetup");
							FineProjectionData[iframe].orientationNumAllClasses += generateProjectionSetup(
									op,
									sp,
									baseMLO,
									false, //not coarse
									exp_iclass,
									FineProjectionData[iframe]);
							CUDA_CPU_TOC("generateProjectionSetup");

						}
						//set a maximum possible size for all weights (to be reduced by significance-checks)
						FinePassWeights[iframe].setDataSize(FineProjectionData[iframe].orientationNumAllClasses*sp.nr_trans*sp.nr_oversampled_trans);
						FinePassWeights[iframe].dual_alloc_all();
					}
					CUDA_CPU_TIC("getAllSquaredDifferencesFine");
					getAllSquaredDifferencesFine(ipass, op, sp, baseMLO, this, FinePassWeights, FinePassClassMasks, FineProjectionData);
					CUDA_CPU_TOC("getAllSquaredDifferencesFine");
					CUDA_CPU_TIC("convertAllSquaredDifferencesToWeightsFine");
					convertAllSquaredDifferencesToWeights(ipass, op, sp, baseMLO, this, FinePassWeights, FinePassClassMasks);
					CUDA_CPU_TOC("convertAllSquaredDifferencesToWeightsFine");
				}

				CUDA_CPU_TOC("weightPass");
			}

			// For the reconstruction step use mymodel.current_size!
			sp.current_image_size = baseMLO->mymodel.current_size;

			CUDA_CPU_TIC("storeWeightedSums");
			storeWeightedSums(op, sp, baseMLO, this, FinePassWeights, FineProjectionData, FinePassClassMasks);
			CUDA_CPU_TOC("storeWeightedSums");

			CUDA_CPU_TIC("Freefalse");
			for (long int iframe = 0; iframe < sp.nr_particles; iframe++)
			{
				for (int exp_iclass = sp.iclass_min; exp_iclass <= sp.iclass_max; exp_iclass++)
				{
					FinePassClassMasks[iframe][exp_iclass].jobOrigin.d_do_free=false;
					FinePassClassMasks[iframe][exp_iclass].jobOrigin.h_do_free=false;
					FinePassClassMasks[iframe][exp_iclass].jobExtent.d_do_free=false;
					FinePassClassMasks[iframe][exp_iclass].jobExtent.h_do_free=false;
				}
			}
			CUDA_CPU_TOC("Freefalse");
			CUDA_CPU_TOC("oneParticle");
		}
	}
	CUDA_CPU_TOC("oneTask");
//	CUDA_CPU_TIC("interParticle");
//	exit(0);
}


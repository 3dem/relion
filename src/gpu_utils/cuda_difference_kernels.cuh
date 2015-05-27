/*
 * cuda_difference_kernels.cuh
 *
 *  Created on: May 26, 2015
 *      Author: bjornf
 */

#ifndef CUDA_DIFFERENCE_KERNELS_CUH_
#define CUDA_DIFFERENCE_KERNELS_CUH_

#include <cuda.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_utils.cuh"

__global__ void cuda_kernel_D2(	    FLOAT *g_refs_real,
										FLOAT *g_refs_imag,
										FLOAT *g_imgs_real,
										FLOAT *g_imgs_imag,
										FLOAT *g_Minvsigma2, FLOAT *g_diff2s,
										unsigned img_size, FLOAT sum_init,
										unsigned long significant_num,
										unsigned long translation_num,
										unsigned long *d_rotidx,
										unsigned long *d_transidx,
									    unsigned long *d_ihidden_overs // TODO use it to map in here, get rid of collect_data_1
									);

__global__ void cuda_kernel_D2_CC(	FLOAT *g_refs_real,
										FLOAT *g_refs_imag,
										FLOAT *g_imgs_real,
										FLOAT *g_imgs_imag,
										FLOAT *g_Minvsigma2, FLOAT *g_diff2s,
										unsigned img_size, FLOAT exp_local_sqrtXi2,
										unsigned long significant_num,
										unsigned long translation_num,
										unsigned long *d_rotidx,
										unsigned long *d_transidx);

__global__ void cuda_kernel_wavg(       FLOAT *g_refs_real,
			                            FLOAT *g_refs_imag,
			                            FLOAT *g_imgs_real,
			                            FLOAT *g_imgs_imag,
			                            FLOAT *g_imgs_nomask_real,
			                            FLOAT *g_imgs_nomask_imag,
			                            FLOAT* g_weights,
			                            FLOAT* g_ctfs,
			                            FLOAT* g_Minvsigma2s,
			                            FLOAT *g_wdiff2s_parts,
			                            FLOAT *g_wavgs_real,
			                            FLOAT *g_wavgs_imag,
			                            FLOAT* g_Fweights,
			                            unsigned long orientation_num,
			                            unsigned long translation_num,
			                            FLOAT weight_norm,
			                            FLOAT significant_weight,
			                            unsigned image_size,
			                            bool refs_are_ctf_corrected);

#endif /* CUDA_DIFFERENCE_KERNELS_CUH_ */

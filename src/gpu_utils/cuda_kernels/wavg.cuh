#ifndef CUDA_WAVG_KERNEL_CUH_
#define CUDA_WAVG_KERNEL_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_utils.cuh"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_settings.h"

__global__ void cuda_kernel_wavg(
		FLOAT *g_eulers,
		Cuda3DProjectorKernel projector,
		unsigned image_size,
		unsigned long orientation_num,
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
		unsigned long translation_num,
		FLOAT weight_norm,
		FLOAT significant_weight,
		bool refs_are_ctf_corrected);

#endif /* CUDA_WAVG_KERNEL_CUH_ */

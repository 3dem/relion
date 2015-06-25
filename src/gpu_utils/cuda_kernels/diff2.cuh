#ifndef CUDA_DIFF2_KERNELS_CUH_
#define CUDA_DIFF2_KERNELS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_utils.cuh"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_settings.h"

// ===================================================
//     Combined Projection+Difference kernels are
// ===================================================
//	-PAV_TTI_D2      Texture Implicit   - single prec. only
//	-PAV_TTE_D2      Texture Explicit   - single prec. only (?)
//	-PAV_TGE_D2      Global  Explicit
//
//   PAV  =   Project All Views
//
//   FIXME: All should be available with suffix _CC  (cross-correlation algorithm)

__global__ void cuda_kernel_diff2_course(
		FLOAT *g_eulers,
		FLOAT *g_imgs_real,
		FLOAT *g_imgs_imag,
		Cuda3DProjectorKernel projector,
		FLOAT *g_Minvsigma2,
		FLOAT *g_diff2s,
		unsigned translation_num,
		int image_size,
		FLOAT sum_init
		);

__global__ void cuda_kernel_diff2_fine(
		FLOAT *g_eulers,
		FLOAT *g_imgs_real,
		FLOAT *g_imgs_imag,
		Cuda3DProjectorKernel projector,
		FLOAT *g_Minvsigma2,
		FLOAT *g_diff2s,
		unsigned image_size,
		FLOAT sum_init,
		unsigned long orientation_num,
		unsigned long translation_num,
		unsigned long todo_blocks,
		unsigned long *d_rot_idx,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num
		);

#endif /* CUDA_DIFF2_KERNELS_CUH_ */

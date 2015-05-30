#ifndef CUDA_PROJDIFF_KERNELS_CUH_
#define CUDA_PROJDIFF_KERNELS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_utils.cuh"

#if !defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)
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
__global__ void cuda_kernel_PAV_TTI_D2( FLOAT *g_eulers,
		                                FLOAT *g_imgs_real,
		                                FLOAT *g_imgs_imag,
										cudaTextureObject_t texModel_real,
										cudaTextureObject_t texModel_imag,
										FLOAT *g_Minvsigma2,
										FLOAT *g_diff2s,
										unsigned image_size,
										FLOAT sum_init,
										unsigned long orientation_num,
										unsigned long translation_num,
										unsigned long todo_blocks,
										unsigned long *d_rotidx,
										unsigned long *d_transidx,
										unsigned long *d_trans_num,
										unsigned long *d_ihidden_overs,
										unsigned my_r_max,
										int max_r2,
										int min_r2_nn,
										long int img_x,
										long int img_y,
										long int mdl_init_y,
										long int mdl_init_z
										);
#elif !defined(CUDA_DOUBLE_PRECISION)
// __global__ void cuda_kernel_PAV_TTE_D2
#else
// __global__ void cuda_kernel_PAV_TGE_D2
#endif // !defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)


#endif /* CUDA_PROJDIFF_KERNELS_CUH_ */

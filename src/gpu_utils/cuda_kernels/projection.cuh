#ifndef CUDA_PROJECTION_KERNELS_CUH_
#define CUDA_PROJECTION_KERNELS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_utils.cuh"

// ===================================================
//     Projection kernels are
// ===================================================
//	-PAV_TTI      Texture Implicit   - single prec. only
//	-PAV_TTE      Texture Explicit   - single prec. only (?)
//	-PAV_TGE      Global  Explicit
//
//   PAV  =   Project All Views


// uses global memory and explicit interpolation = can do double precision.
__global__ void cuda_kernel_PAV_TGE( FLOAT *g_model_real,
									 FLOAT *g_model_imag,
									 FLOAT *g_eulers,
									 FLOAT *g_Frefs_real,
									 FLOAT *g_Frefs_imag,
									 unsigned  my_r_max,
									 int max_r2,
									 int min_r2_nn,
									 unsigned  image_size,
									 unsigned long orientation_num,
									 long int XSIZE_img,
									 long int YSIZE_img,
									 long int XSIZE_mdl,
									 long int YSIZE_mdl,
									 long int STARTINGY_mdl,
									 long int STARTINGZ_mdl,
									 float padding_factor);


#if !defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)
// uses texture memory and implicit (texture) interpolation = requires float precision.
__global__ void cuda_kernel_PAV_TTI( FLOAT *g_eulers,
									 FLOAT *g_Frefs_real,
									 FLOAT *g_Frefs_imag,
									 cudaTextureObject_t texModel_real,
									 cudaTextureObject_t texModel_imag,
									 unsigned  my_r_max,
									 int max_r2,
									 int min_r2_nn,
									 unsigned  image_size,
									 unsigned long orientation_num,
									 long int XSIZE_img,
									 long int YSIZE_img,
									 long int STARTINGY_mdl,
									 long int STARTINGZ_mdl,
									 float padding_factor);

#elif !defined(CUDA_DOUBLE_PRECISION)
// uses texture memory and explicit interpolation = requires float precision.
__global__ void cuda_kernel_PAV_TTE( FLOAT *g_eulers,
									 FLOAT *g_Frefs_real,
									 FLOAT *g_Frefs_imag,
									 cudaTextureObject_t texModel_real,
									 cudaTextureObject_t texModel_imag,
									 unsigned  my_r_max,
									 int max_r2,
									 int min_r2_nn,
									 unsigned  image_size,
									 unsigned long orientation_num,
									 long int XSIZE_img,
									 long int YSIZE_img,
									 long int STARTINGY_mdl,
									 long int STARTINGZ_mdl,
									 float padding_factor);

#endif //!defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)

#endif /* CUDA_PROJECTION_KERNELS_CUH_ */

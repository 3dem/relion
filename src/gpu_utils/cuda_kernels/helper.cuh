#ifndef CUDA_HELPER_KERNELS_CUH_
#define CUDA_HELPER_KERNELS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_device_utils.cuh"
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_projector.h"

__global__ void cuda_kernel_exponentiate_weights_coarse(  XFLOAT *g_pdf_orientation,
									    	  XFLOAT *g_pdf_offset,
									    	  XFLOAT *g_Mweight,
									    	  XFLOAT min_diff2,
									     	  int nr_coarse_orient,
									     	  int nr_coarse_trans);

__global__ void cuda_kernel_exponentiate_weights_fine(    XFLOAT *g_pdf_orientation,
											  XFLOAT *g_pdf_offset,
											  XFLOAT *g_weights,
											  XFLOAT min_diff2,
											  int oversamples_orient,
											  int oversamples_trans,
									     	  unsigned long *d_rot_id,
											  unsigned long *d_trans_idx,
											  unsigned long *d_job_idx,
											  unsigned long *d_job_num,
									     	  long int job_num);

__global__ void cuda_kernel_collect2jobs(	XFLOAT *g_oo_otrans_x,          // otrans-size -> make const
										XFLOAT *g_oo_otrans_y,          // otrans-size -> make const
										XFLOAT *g_myp_oo_otrans_x2y2z2, // otrans-size -> make const
										XFLOAT *g_i_weights,
										XFLOAT op_significant_weight,    // TODO Put in const
										XFLOAT op_sum_weight,            // TODO Put in const
										int   coarse_trans,
										int   oversamples_trans,
										int   oversamples_orient,
										int   oversamples,
										bool  do_ignore_pdf_direction,
										XFLOAT *g_o_weights,
										XFLOAT *g_thr_wsum_prior_offsetx_class,
										XFLOAT *g_thr_wsum_prior_offsety_class,
										XFLOAT *g_thr_wsum_sigma2_offset,
								     	unsigned long *d_rot_idx,
								     	unsigned long *d_trans_idx,
								     	unsigned long *d_job_idx,
								     	unsigned long *d_job_num
								     	);

__global__ void cuda_kernel_softMaskOutsideMap(	XFLOAT *vol,
												long int vol_size,
												long int xdim,
												long int ydim,
												long int zdim,
												long int xinit,
												long int yinit,
												long int zinit,
												bool do_Mnoise,
												XFLOAT radius,
												XFLOAT radius_p,
												XFLOAT cosine_width	);

__global__ void cuda_kernel_centerFFT_2D(XFLOAT *img_in,
										 int image_size,
										 int xdim,
										 int ydim,
										 int xshift,
										 int yshift);

__global__ void cuda_kernel_probRatio(  XFLOAT *d_Mccf,
										XFLOAT *d_Mpsi,
										XFLOAT *d_Maux,
										XFLOAT *d_Mmean,
										XFLOAT *d_Mstddev,
										int image_size,
										XFLOAT normfft,
										XFLOAT sum_ref_under_circ_mask,
										XFLOAT sum_ref2_under_circ_mask,
										XFLOAT expected_Pratio,
										XFLOAT psi);

__global__ void cuda_kernel_rotateAndCtf( CUDACOMPLEX *d_Faux,
						  	  	  	  	  XFLOAT *d_ctf,
						  	  	  	  	  XFLOAT psi,
						  	  			  CudaProjectorKernel projector);

__global__ void cuda_kernel_convol(	 CUDACOMPLEX *d_A,
									 CUDACOMPLEX *d_B,
									 int image_size);

#endif /* CUDA_HELPER_KERNELS_CUH_ */

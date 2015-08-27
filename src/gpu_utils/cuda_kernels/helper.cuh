#ifndef CUDA_HELPER_KERNELS_CUH_
#define CUDA_HELPER_KERNELS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_utils_stl.cuh"

__global__ void cuda_kernel_sumweightCoarse(  XFLOAT *g_pdf_orientation,
									    	  XFLOAT *g_pdf_offset,
									    	  XFLOAT *g_Mweight,
									    	  XFLOAT *g_thisparticle_sumweight,
									    	  XFLOAT min_diff2,
									    	  int oversamples_orient,
									    	  int oversamples_trans,
									    	  int coarse_trans,
									     	  long int sumweight_pos);

__global__ void cuda_kernel_sumweightFine(    XFLOAT *g_pdf_orientation,
											  XFLOAT *g_pdf_offset,
											  XFLOAT *g_weights,
											  XFLOAT *g_thisparticle_sumweight,
											  XFLOAT min_diff2,
											  int oversamples_orient,
											  int oversamples_trans,
									     	  unsigned long *d_rot_id,
											  unsigned long *d_trans_idx,
											  unsigned long *d_job_idx,
											  unsigned long *d_job_num,
									     	  long int job_num,
									     	  long int sumweight_pos);

__global__ void cuda_kernel_collect2(	XFLOAT *g_oo_otrans_x,
										XFLOAT *g_oo_otrans_y,
										XFLOAT *g_myp_oo_otrans_x2y2z2,
										XFLOAT *g_Mweight,
										XFLOAT op_significant_weight,
										XFLOAT op_sum_weight,
										int   coarse_trans,
										int   oversamples_trans,
										int   oversamples_orient,
										int   oversamples,
										bool  do_ignore_pdf_direction,
										XFLOAT *g_weights,
										XFLOAT *g_thr_wsum_prior_offsetx_class,
										XFLOAT *g_thr_wsum_prior_offsety_class,
										XFLOAT *g_thr_wsum_sigma2_offset
										);

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

// Stacks images in place, reducing at most 2*gridDim.x images down to gridDim.x images.
// Ex; 19 -> 16 or 32 -> 16,
__global__ void cuda_kernel_reduce_wdiff2s(XFLOAT *g_wdiff2s_parts,
										   long int orientation_num,
										   int image_size,
										   int current_block_num);

__global__ void cuda_kernel_wdparts_to_wdsum(XFLOAT *g_wdiff2s_parts, XFLOAT *g_wdiff2s_sum, unsigned long image_size);

#endif /* CUDA_HELPER_KERNELS_CUH_ */

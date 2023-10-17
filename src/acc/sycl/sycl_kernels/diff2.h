#ifndef DIFF2_KERNELS_H_
#define DIFF2_KERNELS_H_

#include "src/acc/sycl/sycl_virtual_dev.h"

class AccProjectorKernel;

namespace syclKernels
{
	template<bool REF3D, bool DATA3D, int block_sz, int eulers_per_block, int prefetch_fraction>
#ifdef __INTEL_COMPILER
	inline
#else
	__attribute__((always_inline)) inline
#endif
	void diff2_coarse(
		unsigned long grid_size,
		XFLOAT *g_eulers, XFLOAT *trans_x, XFLOAT *trans_y, XFLOAT *trans_z,
		XFLOAT *g_real,	XFLOAT *g_imag,
		AccProjectorKernel &projector,
		XFLOAT *g_corr, XFLOAT *g_diff2s,
		unsigned long translation_num, unsigned long image_size,
		unsigned long orientation_num, virtualSYCL *devAcc
		);

	template<bool REF3D, bool DATA3D, int block_sz, int chunk_sz>
#ifdef __INTEL_COMPILER
	inline
#else
	__attribute__((always_inline)) inline
#endif
	void diff2_fine(
		unsigned long grid_size,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,		
		AccProjectorKernel &projector,
		XFLOAT *g_corr_img, XFLOAT *g_diff2s,
		unsigned long image_size, XFLOAT sum_init,
		unsigned long orientation_num, unsigned long translation_num, unsigned long num_jobs,
		unsigned long *d_rot_idx, unsigned long *d_trans_idx,
		unsigned long *d_job_idx, unsigned long *d_job_num,
		size_t d_size, virtualSYCL *devAcc
		);

/*
 *   	CROSS-CORRELATION-BASED KERNELS
 */
	template<bool REF3D, bool DATA3D, int block_sz>
#ifdef __INTEL_COMPILER
	inline
#else
	__attribute__((always_inline)) inline
#endif
	void diff2_CC_coarse(
		unsigned long grid_size,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y,
		AccProjectorKernel &projector,
		XFLOAT *g_corr_img, XFLOAT *g_diff2s,
		unsigned long trans_num, unsigned long image_size, XFLOAT exp_local_sqrtXi2,
		unsigned long orientation_num, virtualSYCL *devAcc
		);

	template<bool REF3D, bool DATA3D, int block_sz,int chunk_sz>
#ifdef __INTEL_COMPILER
	inline
#else
	__attribute__((always_inline)) inline
#endif
	void diff2_CC_fine(
		unsigned long grid_size,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,		
		AccProjectorKernel &projector,
		XFLOAT *g_corr_img, XFLOAT *g_diff2s,
		unsigned long image_size, XFLOAT sum_init, XFLOAT exp_local_sqrtXi2,
		unsigned long orientation_num, unsigned long translation_num, unsigned long num_jobs,
		unsigned long *d_rot_idx, unsigned long *d_trans_idx,
		unsigned long *d_job_idx, unsigned long *d_job_num,
		size_t d_size, virtualSYCL *devAcc
		);
} // end of namespace syclKernels

#endif /* DIFF2_KERNELS_H_ */

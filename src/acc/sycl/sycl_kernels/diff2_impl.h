#ifndef DIFF2_IMPL_KERNELS_H_
#define DIFF2_IMPL_KERNELS_H_

#include <vector>
#include <iostream>
#include <limits>
#include <chrono>
#include <algorithm>
#include <cassert>
#include <sycl/sycl.hpp>

#include "src/acc/acc_projectorkernel_impl.h"
#include "src/acc/sycl/sycl_dev.h"
#include "src/acc/sycl/sycl_kernels/sycl_utils.h"
#include "src/acc/sycl/sycl_kernels/diff2_gpu.h"

namespace syclKernels
{
	template<bool REF3D, bool DATA3D, int block_sz, int eulers_per_block, int prefetch_fraction>
	__attribute__((always_inline))
	inline
	void diff2_coarse(
		unsigned long grid_size,
		XFLOAT *g_eulers, XFLOAT *trans_x, XFLOAT *trans_y, XFLOAT *trans_z,
		XFLOAT *g_real, XFLOAT *g_imag,
		AccProjectorKernel &projector,
		XFLOAT *g_corr, XFLOAT *g_diff2s,
		unsigned long trans_num, unsigned long image_size,
		virtualSYCL *devAcc
		)
	{
		devSYCL *dGPU = dynamic_cast<devSYCL*>(devAcc);
		assert( trans_num <= std::numeric_limits<int>::max());
		assert(image_size <= std::numeric_limits<int>::max());
		assert( grid_size <= std::numeric_limits<int>::max());
		assert( grid_size <= dGPU->maxWorkGroup[1]);
		assert(  block_sz <= dGPU->maxItem[2]);
		assert(eulers_per_block*9 + 2*(block_sz/prefetch_fraction*eulers_per_block) + 3*block_sz <= dGPU->localMem);

		sycl::range<3> wi (1,1,block_sz);
		sycl::range<3> wg (1,grid_size,block_sz);
		dGPU->reCalculateRange(wg, wi);
		auto event = dGPU->syclSubmit
		(
			[&](sycl::handler &cgh) //
			{
				using namespace sycl;
				// accessors to device memory
				local_accessor<XFLOAT, 1> s_eulers_acc(range<1>(eulers_per_block*9), cgh);
				local_accessor<XFLOAT, 1> s_ref_real_acc(range<1>(block_sz/prefetch_fraction*eulers_per_block), cgh);
				local_accessor<XFLOAT, 1> s_ref_imag_acc(range<1>(block_sz/prefetch_fraction*eulers_per_block), cgh);
				local_accessor<XFLOAT, 1> s_real_acc(range<1>(block_sz), cgh);
				local_accessor<XFLOAT, 1> s_imag_acc(range<1>(block_sz), cgh);
				local_accessor<XFLOAT, 1> s_corr_acc(range<1>(block_sz), cgh);

				cgh.parallel_for
				(
					nd_range<3>(wg, wi), [=](nd_item<3> nit)
 #if defined(__INTEL_LLVM_COMPILER) && defined(INTEL_SG_SIZE)
					[[intel::reqd_sub_group_size(INTEL_SG_SIZE)]]
 #endif
					{
						syclGpuKernels::sycl_kernel_diff2_coarse<REF3D, DATA3D, block_sz, eulers_per_block, prefetch_fraction>
						( //
							nit,
							const_cast<AccProjectorKernel&>(projector), // Why const_cast is needed for compilation?
							g_eulers, trans_x, trans_y, trans_z, g_real, g_imag,
							g_corr, g_diff2s,
							static_cast<int>(trans_num), static_cast<int>(image_size),
							s_eulers_acc.get_multi_ptr<access::decorated::yes>().get_raw(),
							s_ref_real_acc.get_multi_ptr<access::decorated::yes>().get_raw(),
							s_ref_imag_acc.get_multi_ptr<access::decorated::yes>().get_raw(),
							s_real_acc.get_multi_ptr<access::decorated::yes>().get_raw(),
							s_imag_acc.get_multi_ptr<access::decorated::yes>().get_raw(),
							s_corr_acc.get_multi_ptr<access::decorated::yes>().get_raw()
						);	// End of sycl_kernel_diff2_coarse
					}	// End of cgh.parallel_for Lamda function
				);	// End of cgh.parallel_for
			}	// End of dGPU->syclSubmit Lamda function
		);	// End of dGPU->syclSubmit
	}

	template<bool REF3D, bool DATA3D, int block_sz,int chunk_sz>
	__attribute__((always_inline))
	inline
	void diff2_fine(
		unsigned long grid_size,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,		
		AccProjectorKernel &projector,
		XFLOAT *g_corr_img, XFLOAT *g_diff2s,
		unsigned long image_size, XFLOAT sum_init,
		unsigned long orient_num, unsigned long trans_num, unsigned long num_jobs,
		unsigned long *d_rot_idx, unsigned long *d_trans_idx,
		unsigned long *d_job_idx, unsigned long *d_job_num,
		virtualSYCL *devAcc
		)
	{
		devSYCL *dGPU = dynamic_cast<devSYCL*>(devAcc);
		assert(orient_num <= std::numeric_limits<int>::max());
		assert( trans_num <= std::numeric_limits<int>::max());
		assert(image_size <= std::numeric_limits<int>::max());
		assert(  num_jobs <= std::numeric_limits<int>::max());
		assert( grid_size <= std::numeric_limits<int>::max());
		assert( grid_size <= dGPU->maxWorkGroup[1]);
		assert(  block_sz <= dGPU->maxItem[2]);
		assert(block_sz * chunk_sz <= dGPU->localMem);

		sycl::range<3> wi (1,1,block_sz);
		sycl::range<3> wg (1,grid_size,block_sz);
		dGPU->reCalculateRange(wg, wi);
		auto event = dGPU->syclSubmit
		(
			[&](sycl::handler &cgh) //
			{
				using namespace sycl;
				// accessors to device memory
				local_accessor<XFLOAT, 1> s_acc(range<1>(block_sz * chunk_sz), cgh);

				cgh.parallel_for
				(
					nd_range<3>(wg, wi), [=](nd_item<3> nit)
 #if defined(__INTEL_LLVM_COMPILER) && defined(INTEL_SG_SIZE)
					[[intel::reqd_sub_group_size(INTEL_SG_SIZE)]]
 #endif
					{
						syclGpuKernels::sycl_kernel_diff2_fine<REF3D, DATA3D, block_sz>
						( //
							nit,
							const_cast<AccProjectorKernel&>(projector),
							g_eulers, g_imgs_real, g_imgs_imag,
							g_trans_x, g_trans_y, g_trans_z,
							g_corr_img, g_diff2s,
							static_cast<int>(image_size), sum_init,
							static_cast<int>(num_jobs),
							d_rot_idx, d_trans_idx, d_job_idx, d_job_num,
							s_acc.get_multi_ptr<access::decorated::yes>().get_raw()
						);	// End of sycl_kernel_diff2_fine
					}	// End of cgh.parallel_for Lamda function
				);	// End of cgh.parallel_for
			}	// End of dGPU->syclSubmit Lamda function
		);	// End of dGPU->syclSubmit
	}


/*
 *   	CROSS-CORRELATION-BASED KERNELS
 */
	template<bool REF3D, bool DATA3D, int block_sz>
	__attribute__((always_inline))
	inline
	void diff2_CC_coarse(
		unsigned long grid_size,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
		AccProjectorKernel &projector,
		XFLOAT *g_corr_img, XFLOAT *g_diff2s,
		unsigned long trans_num, unsigned long image_size, XFLOAT exp_local_sqrtXi2,
		virtualSYCL *devAcc
		)
	{
		devSYCL *dGPU = dynamic_cast<devSYCL*>(devAcc);
		assert(grid_size*trans_num <= std::numeric_limits<int>::max());
		assert(image_size <= std::numeric_limits<int>::max());
		assert( grid_size <= dGPU->maxWorkGroup[0]);
		assert( trans_num <= dGPU->maxWorkGroup[1]);
		assert(  block_sz <= dGPU->maxItem[2]);
		assert(2*block_sz <= dGPU->localMem);

		sycl::range<3> wi {1,1,block_sz};
		sycl::range<3> wg {grid_size,trans_num,block_sz};
		auto event = dGPU->syclSubmit
		(
			[&](sycl::handler &cgh) //
			{
				using namespace sycl;
				// accessors to device memory
				local_accessor<XFLOAT, 1> s_weight_acc(range<1>(block_sz), cgh);
				local_accessor<XFLOAT, 1> s_norm_acc(range<1>(block_sz), cgh);

				cgh.parallel_for
				(
					nd_range<3>(wg, wi), [=](nd_item<3> nit)
 #if defined(__INTEL_LLVM_COMPILER) && defined(INTEL_SG_SIZE)
					[[intel::reqd_sub_group_size(INTEL_SG_SIZE)]]
 #endif
					{
						syclGpuKernels::sycl_kernel_diff2_CC_coarse<REF3D, DATA3D, block_sz>
						( //
							nit,
							const_cast<AccProjectorKernel&>(projector), // Why const_cast is needed for compilation?
							g_eulers, g_imgs_real, g_imgs_imag, g_trans_x, g_trans_y, g_trans_z,
							g_corr_img, g_diff2s,
							static_cast<int>(trans_num), static_cast<int>(image_size),
							s_weight_acc.get_multi_ptr<access::decorated::yes>().get_raw(),
							s_norm_acc.get_multi_ptr<access::decorated::yes>().get_raw()
						);	// End of sycl_kernel_diff2_CC_coarse
					}	// End of cgh.parallel_for Lamda function
				);	// End of cgh.parallel_for
			}	// End of dGPU->syclSubmit Lamda function
		);	// End of dGPU->syclSubmit
	}

	template<bool REF3D, bool DATA3D, int block_sz,int chunk_sz>
	__attribute__((always_inline))
	inline
	void diff2_CC_fine(
		unsigned long grid_size,
		XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,		
		AccProjectorKernel &projector,
		XFLOAT *g_corr_img, XFLOAT *g_diff2s,
		unsigned long image_size, XFLOAT sum_init, XFLOAT exp_local_sqrtXi2,
		unsigned long orient_num, unsigned long trans_num, unsigned long num_jobs,
		unsigned long *d_rot_idx, unsigned long *d_trans_idx,
		unsigned long *d_job_idx, unsigned long *d_job_num,
		virtualSYCL *devAcc
		)
	{
		devSYCL *dGPU = dynamic_cast<devSYCL*>(devAcc);
		assert(orient_num <= std::numeric_limits<int>::max());
		assert( trans_num <= std::numeric_limits<int>::max());
		assert(image_size <= std::numeric_limits<int>::max());
		assert(  num_jobs <= std::numeric_limits<int>::max());
		assert( grid_size <= std::numeric_limits<int>::max());
		assert( grid_size <= dGPU->maxWorkGroup[1]);
		assert(  block_sz <= dGPU->maxItem[2]);
		assert(2 * block_sz * chunk_sz <= dGPU->localMem);

		sycl::range<3> wi (1,1,block_sz);
		sycl::range<3> wg (1,grid_size,block_sz);
		dGPU->reCalculateRange(wg, wi);
		auto event = dGPU->syclSubmit
		(
			[&](sycl::handler &cgh) //
			{
				using namespace sycl;
				// accessors to device memory
				local_accessor<XFLOAT, 1> s_acc(range<1>(block_sz * chunk_sz), cgh);
				local_accessor<XFLOAT, 1> s_cc_acc(range<1>(block_sz * chunk_sz), cgh);

				cgh.parallel_for
				(
					nd_range<3>(wg, wi), [=](nd_item<3> nit)
 #if defined(__INTEL_LLVM_COMPILER) && defined(INTEL_SG_SIZE)
					[[intel::reqd_sub_group_size(INTEL_SG_SIZE)]]
 #endif
					{
						syclGpuKernels::sycl_kernel_diff2_CC_fine<REF3D, DATA3D, block_sz>
						( //
							nit,
							const_cast<AccProjectorKernel&>(projector),
							g_eulers, g_imgs_real, g_imgs_imag,
							g_trans_x, g_trans_y, g_trans_z,
							g_corr_img, g_diff2s,
							static_cast<int>(image_size), sum_init, static_cast<int>(num_jobs),
							d_rot_idx, d_trans_idx, d_job_idx, d_job_num,
							s_acc.get_multi_ptr<access::decorated::yes>().get_raw(), s_cc_acc.get_multi_ptr<access::decorated::yes>().get_raw()
						);	// End of sycl_kernel_diff2_CC_fine
					}	// End of cgh.parallel_for Lamda function
				);	// End of cgh.parallel_for
			}	// End of dGPU->syclSubmit Lamda function
		);	// End of dGPU->syclSubmit
	}

} // end of namespace syclKernels

#endif /* DIFF2_IMPL_KERNELS_H_ */

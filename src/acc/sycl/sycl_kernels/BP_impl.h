#ifndef BP_IMPL_KERNELS_H_
#define BP_IMPL_KERNELS_H_


#include <vector>
#include <iostream>
#include <limits>
#include <chrono>
#include <cassert>
#include <sycl/sycl.hpp>

#include "src/acc/acc_backprojector.h"
#include "src/acc/sycl/sycl_dev.h"
#include "src/acc/sycl/sycl_kernels/sycl_utils.h"
#include "src/acc/sycl/sycl_kernels/BP_gpu.h"

namespace syclKernels
{

template <bool CTF_PREMULTIPLIED, bool SGD>
__attribute__((always_inline))
inline
void backproject2D(
		unsigned long imageCount, int block_size,
		AccProjectorKernel &projector,
		XFLOAT *g_img_real, XFLOAT *g_img_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y,
		XFLOAT *g_weights, XFLOAT *g_Minvsigma2s, XFLOAT *g_ctfs,
		unsigned long trans_num, XFLOAT significant_weight, XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real, XFLOAT *g_model_imag, XFLOAT *g_model_weight,
		int max_r, int max_r2, XFLOAT padding_factor,
		int img_x, int img_y, unsigned img_xy, int mdl_x,
		int mdl_inity, int mdl_y,
		virtualSYCL *devAcc)
{
	devSYCL *dGPU = dynamic_cast<devSYCL*>(devAcc);
	assert(  trans_num <= std::numeric_limits<int>::max());
	assert(     img_xy <= std::numeric_limits<int>::max());
	assert(mdl_x*mdl_y <= std::numeric_limits<int>::max());
	assert( imageCount <= std::numeric_limits<int>::max());
	assert( imageCount <= dGPU->maxWorkGroup[1]);
	assert( block_size <= dGPU->maxItem[2]);

	sycl::range<3> wi (1,1,block_size);
	sycl::range<3> wg (1,imageCount,block_size);
	dGPU->reCalculateRange(wg, wi);
	auto event = dGPU->syclSubmit
	(
		[&](sycl::handler &cgh) //
		{
			using namespace sycl;
			// accessors to device memory
			local_accessor<XFLOAT, 1> s_eulers_acc(range<1>(4), cgh);

			cgh.parallel_for
			(
				nd_range<3>(wg, wi), [=](nd_item<3> nit)
 #if defined(__INTEL_LLVM_COMPILER) && defined(INTEL_SG_SIZE)
				[[intel::reqd_sub_group_size(INTEL_SG_SIZE)]]
 #endif
				{
					syclGpuKernels::sycl_kernel_backproject2D<CTF_PREMULTIPLIED, SGD>
					( //
						nit, const_cast<AccProjectorKernel&>(projector),
						g_img_real, g_img_imag, g_trans_x, g_trans_y,
						g_weights, g_Minvsigma2s, g_ctfs,
						static_cast<int>(trans_num), significant_weight, weight_norm,
						g_eulers, g_model_real, g_model_imag, g_model_weight,
						max_r, max_r2, padding_factor,
						img_x, img_y, static_cast<int>(img_xy),
						mdl_x, mdl_inity,
						s_eulers_acc.get_multi_ptr<access::decorated::yes>().get_raw()
					);  // End of sycl_kernel_backproject2D
				}   // End of cgh.parallel_for Lamda function
			);  // End of cgh.parallel_for
		}   // End of dGPU->syclSubmit Lamda function
	);  // End of dGPU->syclSubmit
}

template <bool DATA3D, bool CTF_PREMULTIPLIED, bool SGD>
__attribute__((always_inline))
inline
void backproject3D(
		unsigned long imageCount, int block_size,
		AccProjectorKernel &projector,
		XFLOAT *g_img_real, XFLOAT *g_img_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
		XFLOAT *g_weights, XFLOAT *g_Minvsigma2s, XFLOAT *g_ctfs,
		unsigned long trans_num, XFLOAT significant_weight, XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real, XFLOAT *g_model_imag, XFLOAT *g_model_weight,
		int max_r, int max_r2, XFLOAT padding_factor,
		int img_x, int img_y, int img_z, unsigned img_xyz,
		int mdl_x, int mdl_y, int mdl_inity,
		int mdl_initz, size_t mdl_xyz,
		virtualSYCL *devAcc)
{
	devSYCL *dGPU = dynamic_cast<devSYCL*>(devAcc);
	assert( trans_num <= std::numeric_limits<int>::max());
	assert(   img_xyz <= std::numeric_limits<int>::max());
	assert(   mdl_xyz <= std::numeric_limits<int>::max());
	assert(imageCount <= std::numeric_limits<int>::max());
	assert(imageCount <= dGPU->maxWorkGroup[1]);
	assert(block_size <= dGPU->maxItem[2]);

	sycl::range<3> wi (1,1,block_size);
	sycl::range<3> wg (1,imageCount,block_size);
	dGPU->reCalculateRange(wg, wi);
	auto event = dGPU->syclSubmit
	(
		[&](sycl::handler &cgh) //
		{
			using namespace sycl;
			// accessors to device memory
			local_accessor<XFLOAT, 1> s_eulers_acc(range<1>(9), cgh);

			cgh.parallel_for
			(
				nd_range<3>(wg, wi), [=](nd_item<3> nit)
 #if defined(__INTEL_LLVM_COMPILER) && defined(INTEL_SG_SIZE)
				[[intel::reqd_sub_group_size(INTEL_SG_SIZE)]]
 #endif
				{
					syclGpuKernels::sycl_kernel_backproject3D<DATA3D, CTF_PREMULTIPLIED, SGD>
					( //
						nit, const_cast<AccProjectorKernel&>(projector),
						g_img_real, g_img_imag, g_trans_x, g_trans_y, g_trans_z,
						g_weights, g_Minvsigma2s, g_ctfs,
						static_cast<int>(trans_num), significant_weight, weight_norm,
						g_eulers, g_model_real, g_model_imag, g_model_weight,
						max_r, max_r2, padding_factor,
						img_x, img_y, img_z, static_cast<int>(img_xyz),
						mdl_x, mdl_y, mdl_inity, mdl_initz,
						s_eulers_acc.get_multi_ptr<access::decorated::yes>().get_raw()
					);  // End of sycl_kernel_backproject3D
				}   // End of cgh.parallel_for Lamda function
			);  // End of cgh.parallel_for
		}   // End of dGPU->syclSubmit Lamda function
	);  // End of dGPU->syclSubmit
}

} // namespace

#endif

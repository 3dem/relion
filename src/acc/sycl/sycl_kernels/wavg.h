#ifndef WAVG_KERNEL_H_
#define WAVG_KERNEL_H_

#include "src/acc/sycl/sycl_virtual_dev.h"

class AccProjectorKernel;

namespace syclKernels
{

template<bool REFCTF, bool REF3D, bool DATA3D, int block_sz>
#ifdef __INTEL_COMPILER
inline
#else
__attribute__((always_inline)) inline
#endif
void wavg(
		XFLOAT *g_eulers,
		AccProjectorKernel &projector,
		unsigned long image_size,
		unsigned long orientation_num,
		XFLOAT *g_img_real,
		XFLOAT *g_img_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT *g_trans_z,
		XFLOAT *g_weights,
		XFLOAT *g_ctfs,
		XFLOAT *g_wdiff2s_parts,
		XFLOAT *g_wdiff2s_AA,
		XFLOAT *g_wdiff2s_XA,
		unsigned long trans_num,
		XFLOAT  weight_norm,
		XFLOAT  significant_weight,
		XFLOAT  part_scale,
		virtualSYCL *devACC);

} // end of namespace syclKernels

#endif /* WAVG_KERNEL_H_ */

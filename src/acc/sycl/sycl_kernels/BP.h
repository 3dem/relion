#ifndef BP_KERNELS_H_
#define BP_KERNELS_H_

#include "src/acc/acc_backprojector.h"

namespace syclKernels
{

template <bool CTF_PREMULTIPLIED, bool SGD>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline
void backproject2D(
		unsigned long imageCount, int block_size,
		AccProjectorKernel &projector,
		XFLOAT *g_img_real, XFLOAT *g_img_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y,
		XFLOAT *g_weights, XFLOAT *g_Minvsigma2s, XFLOAT *g_ctfs,
		unsigned long translation_num, XFLOAT significant_weight, XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real, XFLOAT *g_model_imag, XFLOAT *g_model_weight,
		int max_r, int max_r2, XFLOAT padding_factor,
		unsigned img_x, unsigned img_y, unsigned img_xy, unsigned mdl_x,
		int mdl_inity, size_t mdl_xyz,
		virtualSYCL *devAcc);

template <bool DATA3D, bool CTF_PREMULTIPLIED, bool SGD>
#ifndef __INTEL_COMPILER
__attribute__((always_inline))
#endif
inline
void backproject3D(
		unsigned long imageCount, int block_size,
		AccProjectorKernel &projector,
		XFLOAT *g_img_real, XFLOAT *g_img_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
		XFLOAT *g_weights, XFLOAT *g_Minvsigma2s, XFLOAT *g_ctfs,
		unsigned long translation_num, XFLOAT significant_weight, XFLOAT weight_norm,
		XFLOAT *g_eulers,
		XFLOAT *g_model_real, XFLOAT *g_model_imag, XFLOAT *g_model_weight,
		int max_r, int max_r2, XFLOAT padding_factor,
		unsigned img_x, unsigned img_y, unsigned img_z, unsigned img_xyz,
		unsigned mdl_x, unsigned mdl_y, int mdl_inity,
		int mdl_initz, size_t mdl_xyz,
		virtualSYCL *devAcc);

} // namespace

#endif

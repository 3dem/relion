#ifndef WAVG_GPU_KERNEL_H_
#define WAVG_GPU_KERNEL_H_

#include <vector>
#include <iostream>
#include <fstream>
#include <sycl/sycl.hpp>

#include "src/acc/acc_projector.h"
#include "src/acc/sycl/sycl_settings.h"
#include "src/acc/sycl/sycl_kernels/sycl_utils.h"
#include "src/acc/sycl/sycl_kernels/helper.h"
#include "src/acc/sycl/sycl_kernels/helper_gpu.h"

namespace syclGpuKernels
{

template <bool REFCTF, bool REF3D, bool DATA3D, int block_sz>
void sycl_kernel_wavg(
		sycl::nd_item<3> nit, AccProjectorKernel &projector,
		XFLOAT *g_eulers, int image_size, int orientation_num,
		XFLOAT *g_img_real, XFLOAT *g_img_imag,
		XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
		XFLOAT *g_weights, XFLOAT *g_ctfs,
		XFLOAT *g_wdiff2s_parts, XFLOAT *g_wdiff2s_AA, XFLOAT *g_wdiff2s_XA,
		int translation_num, XFLOAT weight_norm, XFLOAT significant_weight, XFLOAT part_scale,
		XFLOAT *s_parts, XFLOAT *s_sumAA, XFLOAT *s_sumXA, XFLOAT *s_eulers
		)
{
	const int bid = nit.get_group_linear_id();
	const int tid = nit.get_local_id(2);

	const int xSize = projector.imgX;
	const int ySize = projector.imgY;
	const int zSize = projector.imgZ;
	const int maxR = projector.maxR;

	const XFLOAT inv_weight_norm = 1.0f / weight_norm;

	if (tid < 9)
		s_eulers[tid] = g_eulers[bid*9+tid];

	__group_barrier(nit);

	int pass_num {image_size/block_sz + 1};
	for (int pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
	{
		s_parts[tid] = 0.0f;
		s_sumXA[tid] = 0.0f;
		s_sumAA[tid] = 0.0f;

		int pixel = pass*block_sz + tid;
		if (pixel < image_size)
		{
			int x, y, z, xy;
			if (DATA3D)
			{
				z =  pixel / (xSize*ySize);
				xy = pixel % (xSize*ySize);
				x =     xy % xSize;
				y =     xy / xSize;
				if (z > maxR)
				{
					if (z >= zSize - maxR)
						z = z - zSize;
					else
						x = maxR;
				}
			}
			else
			{
				x = pixel % xSize;
				y = pixel / xSize;
			}
			if (y > maxR)
			{
				if (y >= ySize - maxR)
					y = y - ySize;
				else
					x = maxR;
			}

			XFLOAT ref_real, ref_imag;
			if (DATA3D)
				projector.project3Dmodel(
					x,y,z,
					s_eulers[0], s_eulers[1], s_eulers[2],
					s_eulers[3], s_eulers[4], s_eulers[5],
					s_eulers[6], s_eulers[7], s_eulers[8],
					ref_real, ref_imag);
			else if (REF3D)
				projector.project3Dmodel(
					x,y,
					s_eulers[0], s_eulers[1],
					s_eulers[3], s_eulers[4],
					s_eulers[6], s_eulers[7],
					ref_real, ref_imag);
			else
				projector.project2Dmodel(
					x,y,
					s_eulers[0], s_eulers[1],
					s_eulers[3], s_eulers[4],
					ref_real, ref_imag);

			if (REFCTF)
			{
				ref_real *= g_ctfs[pixel];
				ref_imag *= g_ctfs[pixel];
			}
			else
			{
				ref_real *= part_scale;
				ref_imag *= part_scale;
			}

			XFLOAT img_real = g_img_real[pixel];
			XFLOAT img_imag = g_img_imag[pixel];
			for (int itrans = 0; itrans < translation_num; itrans++)
			{
				XFLOAT weight = g_weights[bid*translation_num + itrans];

				if (weight >= significant_weight)
				{
					weight *= inv_weight_norm;

					XFLOAT trans_real, trans_imag;
					if (DATA3D)
						translatePixel(x, y, z, g_trans_x[itrans], g_trans_y[itrans], g_trans_z[itrans], img_real, img_imag, trans_real, trans_imag);
					else
						translatePixel(x, y,    g_trans_x[itrans], g_trans_y[itrans],                    img_real, img_imag, trans_real, trans_imag);

					XFLOAT diff_real = ref_real - trans_real;
					XFLOAT diff_imag = ref_imag - trans_imag;

					s_parts[tid] += weight * (diff_real*diff_real  + diff_imag*diff_imag);
					s_sumXA[tid] += weight * ( ref_real*trans_real + ref_imag*trans_imag);
					s_sumAA[tid] += weight * ( ref_real*ref_real   +  ref_imag*ref_imag );
				}
			}

			g_wdiff2s_XA[pixel] += s_sumXA[tid];
			g_wdiff2s_AA[pixel] += s_sumAA[tid];
//			atomic_ref_dev(    g_wdiff2s_XA[pixel] ).fetch_add(s_sumXA[tid]);
//			atomic_ref_dev(    g_wdiff2s_AA[pixel] ).fetch_add(s_sumAA[tid]);
			atomic_ref_dev( g_wdiff2s_parts[pixel] ).fetch_add(s_parts[tid]);
		}
	}
}

} // end of namespace syclGpuKernels

#endif /* WAVG_GPU_KERNEL_H_ */

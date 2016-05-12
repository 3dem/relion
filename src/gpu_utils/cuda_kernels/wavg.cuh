#ifndef CUDA_WAVG_KERNEL_CUH_
#define CUDA_WAVG_KERNEL_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_device_utils.cuh"

template<bool do_3DProjection>
__global__ void cuda_kernel_wavg(
		XFLOAT *g_eulers,
		CudaProjectorKernel projector,
		unsigned image_size,
		unsigned long orientation_num,
		XFLOAT *g_imgs_real,
		XFLOAT *g_imgs_imag,
		XFLOAT *g_trans_x,
		XFLOAT *g_trans_y,
		XFLOAT* g_weights,
		XFLOAT* g_ctfs,
		XFLOAT *g_wdiff2s_parts,
		XFLOAT *g_wdiff2s_AA,
		XFLOAT *g_wdiff2s_XA,
		unsigned long translation_num,
		XFLOAT weight_norm,
		XFLOAT significant_weight,
		bool refs_are_ctf_corrected,
		XFLOAT part_scale)
{
	XFLOAT ref_real, ref_imag, real, imag;

	int bid = blockIdx.x; //block ID
	int tid = threadIdx.x;

	__shared__ XFLOAT s_eulers[9];

	if (tid < 9)
		s_eulers[tid] = g_eulers[bid*9+tid];
	__syncthreads();

	unsigned pass_num(ceilfracf(image_size,WAVG_BLOCK_SIZE)),pixel;
	__shared__ XFLOAT s_wdiff2s_parts[WAVG_BLOCK_SIZE];
	__shared__ XFLOAT s_sumXA[WAVG_BLOCK_SIZE];
	__shared__ XFLOAT s_sumA2[WAVG_BLOCK_SIZE];

	for (unsigned pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
	{
		s_wdiff2s_parts[tid] = 0.0f;
		s_sumXA[tid] = 0.0f;
		s_sumA2[tid] = 0.0f;

		pixel = pass * WAVG_BLOCK_SIZE + tid;

		if(pixel<image_size)
		{
			int x = pixel % projector.imgX;
			int y = floorfracf(pixel,projector.imgX);

			if (y > projector.maxR)
			{
				if (y >= projector.imgY - projector.maxR)
					y = y - projector.imgY;
				else
					x = projector.maxR;
			}

			if(do_3DProjection)
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

			if (refs_are_ctf_corrected) //FIXME Create two kernels for the different cases
			{
				ref_real *= g_ctfs[pixel];
				ref_imag *= g_ctfs[pixel];
			}
			else
			{
				ref_real *= part_scale;
				ref_imag *= part_scale;
			}

			for (unsigned long itrans = 0; itrans < translation_num; itrans++)
			{
				XFLOAT weight = __ldg(&g_weights[bid * translation_num + itrans]);

				if (weight >= significant_weight)
				{
					weight /= weight_norm;

					unsigned long img_pixel_idx = itrans * image_size + pixel;

					translatePixel(x, y, g_trans_x[itrans], g_trans_y[itrans], g_imgs_real[pixel], g_imgs_imag[pixel], real, imag);

					XFLOAT diff_real = ref_real - real;
					XFLOAT diff_imag = ref_imag - imag;

					s_wdiff2s_parts[tid] += weight * (diff_real*diff_real + diff_imag*diff_imag);

					s_sumXA[tid] +=  weight * ( ref_real * real + ref_imag * imag);
					s_sumA2[tid] +=  weight * ( ref_real*ref_real  +  ref_imag*ref_imag );
				}
			}

			cuda_atomic_add(&g_wdiff2s_XA[pixel], s_sumXA[tid]);
			cuda_atomic_add(&g_wdiff2s_AA[pixel], s_sumA2[tid]);
			cuda_atomic_add(&g_wdiff2s_parts[pixel], s_wdiff2s_parts[tid]);
		}
	}
}

#endif /* CUDA_WAVG_KERNEL_CUH_ */

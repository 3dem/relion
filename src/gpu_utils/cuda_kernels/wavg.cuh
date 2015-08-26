#ifndef CUDA_WAVG_KERNEL_CUH_
#define CUDA_WAVG_KERNEL_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_projector.cuh"
#include "src/gpu_utils/cuda_settings.h"

template<bool do_3DProjection>
__global__ void cuda_kernel_wavg(
		XFLOAT *g_eulers,
		CudaProjectorKernel projector,
		unsigned image_size,
		unsigned long orientation_num,
		XFLOAT *g_imgs_real,
		XFLOAT *g_imgs_imag,
		XFLOAT *g_imgs_nomask_real,
		XFLOAT *g_imgs_nomask_imag,
		XFLOAT* g_weights,
		XFLOAT* g_ctfs,
		XFLOAT* g_Minvsigma2s,
		XFLOAT *g_wdiff2s_parts,
		XFLOAT *g_wavgs_real,
		XFLOAT *g_wavgs_imag,
		XFLOAT* g_Fweights,
		unsigned long translation_num,
		XFLOAT weight_norm,
		XFLOAT significant_weight,
		bool refs_are_ctf_corrected)
{
	XFLOAT ref_real, ref_imag;
	int bid = blockIdx.y * gridDim.x + blockIdx.x; //block ID
	int tid = threadIdx.x;
	// inside the padded 2D orientation grid
//	if( bid < orientation_num )
//	{
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  )),pixel;
		XFLOAT Fweight;
		__shared__ XFLOAT s_wavgs_real[BLOCK_SIZE];
		__shared__ XFLOAT s_wavgs_imag[BLOCK_SIZE];
		__shared__ XFLOAT s_wdiff2s_parts[BLOCK_SIZE];
		__shared__ XFLOAT s_Minvsigma2s[BLOCK_SIZE];
		for (unsigned pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
		{
			s_wavgs_real[tid]  = 0.0f;
			s_wavgs_imag[tid]  = 0.0f;
			s_wdiff2s_parts[tid] = 0.0f;
			Fweight = 0.0f;

			pixel = pass * BLOCK_SIZE + tid;
			s_Minvsigma2s[tid]=g_Minvsigma2s[pixel];

			if(pixel<image_size)
			{
				unsigned long ref_pixel = bid * image_size + pixel;
				if(do_3DProjection)
					projector.project3Dmodel(
						pixel,
						__ldg(&g_eulers[bid*9  ]), __ldg(&g_eulers[bid*9+1]),
						__ldg(&g_eulers[bid*9+3]), __ldg(&g_eulers[bid*9+4]),
						__ldg(&g_eulers[bid*9+6]), __ldg(&g_eulers[bid*9+7]),
						ref_real, ref_imag);
				else
					projector.project2Dmodel(
						pixel,
						__ldg(&g_eulers[bid*9  ]), __ldg(&g_eulers[bid*9+1]),
						__ldg(&g_eulers[bid*9+3]), __ldg(&g_eulers[bid*9+4]),
						__ldg(&g_eulers[bid*9+6]), __ldg(&g_eulers[bid*9+7]),
						ref_real, ref_imag);

				if (refs_are_ctf_corrected) //FIXME Create two kernels for the different cases
				{
					ref_real *= __ldg(&g_ctfs[pixel]);
					ref_imag *= __ldg(&g_ctfs[pixel]);
				}

				for (unsigned long itrans = 0; itrans < translation_num; itrans++)
				{
					XFLOAT weight = __ldg(&g_weights[bid * translation_num + itrans]);

					if (weight >= significant_weight)
					{
						weight /= weight_norm;

						unsigned long img_pixel_idx = itrans * image_size + pixel;

						XFLOAT diff_real = ref_real - g_imgs_real[img_pixel_idx];    // TODO  Put in texture (in such a way that fetching of next image might hit in cache)
						XFLOAT diff_imag = ref_imag - g_imgs_imag[img_pixel_idx];

						s_wdiff2s_parts[tid] += weight * (diff_real*diff_real + diff_imag*diff_imag);

						XFLOAT weightxinvsigma2 = weight * __ldg(&g_ctfs[pixel]) * s_Minvsigma2s[tid];

						s_wavgs_real[tid] += g_imgs_nomask_real[img_pixel_idx] * weightxinvsigma2;    // TODO  Put in texture (in such a way that fetching of next image might hit in cache)
						s_wavgs_imag[tid] += g_imgs_nomask_imag[img_pixel_idx] * weightxinvsigma2;

						Fweight += weightxinvsigma2 * __ldg(&g_ctfs[pixel]);
					}
				}
				g_wavgs_real[ref_pixel] += s_wavgs_real[tid];
				g_wavgs_imag[ref_pixel] += s_wavgs_imag[tid];
				g_wdiff2s_parts[ref_pixel] = s_wdiff2s_parts[tid]; //TODO this could be further reduced in here
				g_Fweights[ref_pixel] += Fweight; //TODO should be buffered into shared
			}
		}
//	}
}

#endif /* CUDA_WAVG_KERNEL_CUH_ */

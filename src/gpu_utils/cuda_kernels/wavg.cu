#include "src/gpu_utils/cuda_kernels/wavg.cuh"
#include <vector>
#include <iostream>



__global__ void cuda_kernel_wavg(
		FLOAT *g_eulers,
		Cuda3DProjectorKernel projector,
		unsigned image_size,
		unsigned long orientation_num,
		FLOAT *g_imgs_real,
		FLOAT *g_imgs_imag,
		FLOAT *g_imgs_nomask_real,
		FLOAT *g_imgs_nomask_imag,
		FLOAT* g_weights,
		FLOAT* g_ctfs,
		FLOAT* g_Minvsigma2s,
		FLOAT *g_wdiff2s_parts,
		FLOAT *g_wavgs_real,
		FLOAT *g_wavgs_imag,
		FLOAT* g_Fweights,
		unsigned long translation_num,
		FLOAT weight_norm,
		FLOAT significant_weight,
		bool refs_are_ctf_corrected)
{
	FLOAT ref_real, ref_imag;
	int bid = blockIdx.y * gridDim.x + blockIdx.x; //block ID
	int tid = threadIdx.x;
	// inside the padded 2D orientation grid
//	if( bid < orientation_num )
//	{
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  )),pixel;
		FLOAT Fweight;
		__shared__ FLOAT s_wavgs_real[BLOCK_SIZE];
		__shared__ FLOAT s_wavgs_imag[BLOCK_SIZE];
		__shared__ FLOAT s_wdiff2s_parts[BLOCK_SIZE];
		__shared__ FLOAT s_Minvsigma2s[BLOCK_SIZE];
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

				projector.project3Dmodel(
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
					FLOAT weight = __ldg(&g_weights[bid * translation_num + itrans]);

					if (weight >= significant_weight)
					{
						weight /= weight_norm;

						unsigned long img_pixel_idx = itrans * image_size + pixel;

						FLOAT diff_real = ref_real - g_imgs_real[img_pixel_idx];    // TODO  Put in texture (in such a way that fetching of next image might hit in cache)
						FLOAT diff_imag = ref_imag - g_imgs_imag[img_pixel_idx];

						s_wdiff2s_parts[tid] += weight * (diff_real*diff_real + diff_imag*diff_imag);

						FLOAT weightxinvsigma2 = weight * __ldg(&g_ctfs[pixel]) * s_Minvsigma2s[tid];

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

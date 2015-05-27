/*
 * cuda_difference_kernels.cu

 *
 *  Created on: May 26, 2015
 *      Author: bjornf
 */
#include "src/gpu_utils/cuda_difference_kernels.cuh"
#include <vector>
#include <iostream>

__global__ void cuda_kernel_D2(	FLOAT *g_refs_real,
									FLOAT *g_refs_imag,
									FLOAT *g_imgs_real,
									FLOAT *g_imgs_imag,
									FLOAT *g_Minvsigma2, FLOAT *g_diff2s,
									unsigned img_size, FLOAT sum_init,
									unsigned long significant_num,
									unsigned long translation_num,
									unsigned long *d_rotidx,
									unsigned long *d_transidx,
									unsigned long *d_ihidden_overs // TODO use it to map in here, get rid of collect_data_1
									)
{
	// blockid
	int ex = blockIdx.y * gridDim.x + blockIdx.x;
    int tid = threadIdx.x;

	// inside the padded 2D orientation grid
	if( ex < significant_num )
	{
		// index of comparison
		unsigned long int ix=d_rotidx[ex];
		unsigned long int iy=d_transidx[ex];

		__shared__ FLOAT s[BLOCK_SIZE];
		s[tid] = 0.0f;

		unsigned pass_num(ceilf((float)img_size/(float)BLOCK_SIZE)), pixel;

		unsigned long ref_start(ix * img_size);
		unsigned long img_start(iy * img_size);
		unsigned long ref_pixel_idx;
		unsigned long img_pixel_idx;

		for (unsigned pass = 0; pass < pass_num; pass ++)
		{
			pixel = pass * BLOCK_SIZE + tid;

			if (pixel < img_size) //Is inside image
			{
				ref_pixel_idx = ref_start + pixel;
				img_pixel_idx = img_start + pixel;

				FLOAT diff_real = __ldg(&g_refs_real[ref_pixel_idx]) - __ldg(&g_imgs_real[img_pixel_idx]); // TODO  Put g_img_* in texture (in such a way that fetching of next image might hit in cache)
				FLOAT diff_imag = __ldg(&g_refs_imag[ref_pixel_idx]) - __ldg(&g_imgs_imag[img_pixel_idx]);

				s[tid] += (diff_real * diff_real + diff_imag * diff_imag) * 0.5f * __ldg(&g_Minvsigma2[pixel]);
			}
		}
		__syncthreads();

		for(int j=(BLOCK_SIZE/2); j>0; j>>=1)
		{
			if(tid<j)
			{
				s[tid] += s[tid+j];
			}
			__syncthreads();
		}
//		if (threadIdx.x*ex == 0)
		{
			g_diff2s[ix * translation_num + iy] = s[0]+sum_init;
		}
		// -------------------------------------------------------------------------
	}
}

__global__ void cuda_kernel_D2_CC(	FLOAT *g_refs_real,
										FLOAT *g_refs_imag,
										FLOAT *g_imgs_real,
										FLOAT *g_imgs_imag,
										FLOAT *g_Minvsigma2, FLOAT *g_diff2s,
										unsigned img_size, FLOAT exp_local_sqrtXi2,
										unsigned long significant_num,
										unsigned long translation_num,
										unsigned long *d_rotidx,
										unsigned long *d_transidx)
{
	// blockid
	int ex = blockIdx.y * gridDim.x + blockIdx.x;
	// inside the padded 2D orientation grid
	if( ex < significant_num )
	{
		// index of comparison
		unsigned long int ix=d_rotidx[ex];
		unsigned long int iy=d_transidx[ex];
		__shared__ double    s[BLOCK_SIZE];
		__shared__ double norm[BLOCK_SIZE];
		s[threadIdx.x] = 0;
		unsigned pass_num(ceilf((float)img_size/(float)BLOCK_SIZE));
		unsigned long pixel,
		ref_start(ix * img_size),
		img_start(iy * img_size);
		unsigned long ref_pixel_idx;
		unsigned long img_pixel_idx;
		for (unsigned pass = 0; pass < pass_num; pass ++)
		{
			pixel = pass * BLOCK_SIZE + threadIdx.x;

			if (pixel < img_size) //Is inside image
			{
				ref_pixel_idx = ref_start + pixel;
				img_pixel_idx = img_start + pixel;

				double diff_real = g_refs_real[ref_pixel_idx] * g_imgs_real[img_pixel_idx];
				double diff_imag = g_refs_imag[ref_pixel_idx] * g_imgs_imag[img_pixel_idx];

				double nR = g_refs_real[ref_pixel_idx]*g_refs_real[ref_pixel_idx];
				double nI = g_refs_imag[ref_pixel_idx]*g_refs_imag[ref_pixel_idx];

				s[threadIdx.x] -= (diff_real + diff_imag);
				norm[threadIdx.x] += nR+nI;
			}
		}
		// -------------------------------------------------------------------------
		__syncthreads();
		int trads = 32;
		int itr = BLOCK_SIZE/trads;
		if(threadIdx.x<trads)
		{
			for(int i=1; i<itr; i++)
			{
				s[threadIdx.x] += s[i*trads + threadIdx.x];
				norm[threadIdx.x] += norm[i*trads + threadIdx.x];
			}
		}
		for(int j=(trads/2); j>0; j/=2)
		{
			if(threadIdx.x<j)
			{
				s[threadIdx.x] += s[threadIdx.x+j];
				norm[threadIdx.x] += norm[threadIdx.x+j];
			}
		}
		__syncthreads();
		// -------------------------------------------------------------------------
		g_diff2s[ix * translation_num + iy] = s[0]/(sqrt(norm[0])*exp_local_sqrtXi2);
	}
}

__global__ void cuda_kernel_wavg(
		FLOAT *g_refs_real,
		FLOAT *g_refs_imag,
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
		unsigned long orientation_num,
		unsigned long translation_num,
		FLOAT weight_norm,
		FLOAT significant_weight,
		unsigned image_size,
		bool refs_are_ctf_corrected)
{
	unsigned long bid = blockIdx.y*gridDim.x + blockIdx.x;
	unsigned tid = threadIdx.x;

	unsigned pass_num(ceilf((float)image_size/(float)BLOCK_SIZE)),pixel;
	FLOAT Fweight;
	__shared__ FLOAT s_wavgs_real[BLOCK_SIZE];
	__shared__ FLOAT s_wavgs_imag[BLOCK_SIZE];
	__shared__ FLOAT s_wdiff2s_parts[BLOCK_SIZE];
	__shared__ FLOAT s_Minvsigma2s[BLOCK_SIZE];
//	if( bid < orientation_num )
//	{
		for (unsigned pass = 0; pass < pass_num; pass ++)
		{
			s_wavgs_real[tid]  = 0.0f;
			s_wavgs_imag[tid]  = 0.0f;
			s_wdiff2s_parts[tid] = 0.0f;
			Fweight = 0;

			pixel = pass * BLOCK_SIZE + tid;
			s_Minvsigma2s[tid]=g_Minvsigma2s[pixel];

			if (pixel < image_size)
			{
				unsigned long ref_pixel = bid * image_size + pixel;
				FLOAT ref_real = g_refs_real[ref_pixel];
				FLOAT ref_imag = g_refs_imag[ref_pixel];
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


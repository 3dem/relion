/*
 * cuda_ProjDiff_kenrels.cu

 *
 *  Created on: May 27, 2015
 *      Author: bjornf
 */
#include "src/gpu_utils/cuda_ProjDiff_kernels.cuh"
#include <vector>
#include <iostream>

#if !defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)

__global__ void cuda_kernel_PAV_TTI_D2( FLOAT *g_eulers,
		                                FLOAT *g_imgs_real,
		                                FLOAT *g_imgs_imag,
										cudaTextureObject_t texModel_real,
										cudaTextureObject_t texModel_imag,
										FLOAT *g_Minvsigma2,
										FLOAT *g_diff2s,
										int image_size,
										FLOAT sum_init,
										int orientation_num,
										int translation_num,
										int significant_num,
										unsigned long *d_rotidx,
										unsigned long *d_transidx,
										unsigned long *d_trans_num,
										unsigned long *d_ihidden_overs,
										int my_r_max,
										int max_r2,
										int min_r2_nn,
										long int img_x,
										long int img_y,
										long int mdl_init_y,
										long int mdl_init_z
										)
{
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;

	FLOAT xp, yp, zp;
	long int r2;
	int pixel;
	bool is_neg_x;
	FLOAT ref_real;
	FLOAT ref_imag;

	// inside the padded 2D orientation grid
	if( bid < significant_num ) // we only need to make
	{
		__shared__ FLOAT s[BLOCK_SIZE];
		s[tid] = 0.0f;

		// index of comparison
		unsigned long int ix=d_rotidx[bid];
		unsigned long int iy=d_transidx[bid];

		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  ));
		unsigned long img_start(iy * image_size);
		unsigned long img_pixel_idx;

		for (unsigned pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
		{
			pixel = (pass * BLOCK_SIZE) + tid;
			if(pixel<image_size)
			{
				int x = pixel % img_x;
				int y = (int)floorf( (float)pixel / (float)img_x);
				img_pixel_idx = img_start + pixel;

				// Dont search beyond square with side max_r
				if (y > my_r_max)
				{
					if (y >= img_y - my_r_max)
						y = y - img_y ;
					else
						x=r2;
				}

				r2 = x*x + y*y;
				if (r2 <= max_r2)
				{
					xp = __ldg(&g_eulers[ix*9])   * x + __ldg(&g_eulers[ix*9+1]) * y;  // FIXME: xp,yp,zp has has accuracy loss
					yp = __ldg(&g_eulers[ix*9+3]) * x + __ldg(&g_eulers[ix*9+4]) * y;  // compared to CPU-based projection. This
					zp = __ldg(&g_eulers[ix*9+6]) * x + __ldg(&g_eulers[ix*9+7]) * y;  // propagates to dx00, dx10, and so on.
					// Only asymmetric half is stored
					if (xp < 0)
					{
						// Get complex conjugated hermitian symmetry pair
						xp = -xp;
						yp = -yp;
						zp = -zp;
						is_neg_x = true;
					}
					else
					{
						is_neg_x = false;
					}
					yp -= mdl_init_y;
					zp -= mdl_init_z;

					ref_real=tex3D<FLOAT>(texModel_real,xp+0.5f,yp+0.5f,zp+0.5f);
					ref_imag=tex3D<FLOAT>(texModel_imag,xp+0.5f,yp+0.5f,zp+0.5f);

//					printf("%i, %i", x,y);
//					printf("%f, %f,%f", xp,yp,zp);
					if (is_neg_x)
					{
						ref_imag = -ref_imag;
					}
				}
				else
				{
					ref_real=0.0f;
					ref_imag=0.0f;
				}
				FLOAT diff_real =  ref_real - __ldg(&g_imgs_real[img_pixel_idx]); // TODO  Put g_img_* in texture (in such a way that fetching of next image might hit in cache)
				FLOAT diff_imag =  ref_imag - __ldg(&g_imgs_imag[img_pixel_idx]);

				s[tid] += (diff_real * diff_real + diff_imag * diff_imag) * 0.5f * __ldg(&g_Minvsigma2[pixel]);
//				printf(" diffs = %f, %f \n",ref_real,img_pixel_idx);
//				printf(" diffs = %i, %i ,%i \n",x,y);
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
		if (tid == 0)
		{
			g_diff2s[ix * translation_num + iy] = s[0]+sum_init;
		}
	}
}
#elif !defined(CUDA_DOUBLE_PRECISION)
// __global__ void cuda_kernel_PAV_TTE_D2
#else
// __global__ void cuda_kernel_PAV_TGE_D2
#endif // !defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)

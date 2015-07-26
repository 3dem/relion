#include "src/gpu_utils/cuda_kernels/diff2.cuh"
#include <vector>
#include <iostream>

__global__ void cuda_kernel_diff2_coarse(
		FLOAT *g_eulers,
		FLOAT *g_imgs_real,
		FLOAT *g_imgs_imag,
		Cuda3DProjectorKernel projector,
		FLOAT *g_Minvsigma2,
		FLOAT *g_diff2s,
		unsigned translation_num,
		int image_size,
		FLOAT sum_init
		)
{
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;

	FLOAT ref_real;
	FLOAT ref_imag;

	FLOAT e0,e1,e3,e4,e6,e7;
	e0 = __ldg(&g_eulers[bid*9  ]);
	e1 = __ldg(&g_eulers[bid*9+1]);
	e3 = __ldg(&g_eulers[bid*9+3]);
	e4 = __ldg(&g_eulers[bid*9+4]);
	e6 = __ldg(&g_eulers[bid*9+6]);
	e7 = __ldg(&g_eulers[bid*9+7]);

	extern __shared__ FLOAT s_cuda_kernel_diff2s[];

	for (unsigned i = 0; i < translation_num; i++)
		s_cuda_kernel_diff2s[translation_num * tid + i] = 0.0f;

	unsigned pixel_pass_num( ceilf( (float)image_size / (float)BLOCK_SIZE ) );
	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
	{
		unsigned pixel = (pass * BLOCK_SIZE) + tid;

		if(pixel < image_size)
		{
			projector.project3Dmodel(
					pixel,
					e0,e1,e3,e4,e6,e7,
					ref_real, ref_imag);

			for (int itrans = 0; itrans < translation_num; itrans ++)
			{
				unsigned long img_pixel_idx = itrans * image_size + pixel;

				FLOAT diff_real =  ref_real - __ldg(&g_imgs_real[img_pixel_idx]);
				FLOAT diff_imag =  ref_imag - __ldg(&g_imgs_imag[img_pixel_idx]);
				FLOAT diff2 = (diff_real * diff_real + diff_imag * diff_imag) * 0.5f * __ldg(&g_Minvsigma2[pixel]);

				s_cuda_kernel_diff2s[translation_num * tid + itrans] += diff2;
			}
		}
	}

	__syncthreads();

	unsigned trans_pass_num( ceilf( (float)translation_num / (float)BLOCK_SIZE ) );
	for (unsigned pass = 0; pass < trans_pass_num; pass++)
	{
		unsigned itrans = (pass * BLOCK_SIZE) + tid;
		if (itrans < translation_num)
		{
			FLOAT sum(sum_init);
			for (unsigned i = 0; i < BLOCK_SIZE; i++)
				sum += s_cuda_kernel_diff2s[i * translation_num + itrans];

			g_diff2s[bid * translation_num + itrans] = sum;
		}
	}
}

__global__ void cuda_kernel_diff2_fine(
		FLOAT *g_eulers,
		FLOAT *g_imgs_real,
		FLOAT *g_imgs_imag,
		Cuda3DProjectorKernel projector,
		FLOAT *g_Minvsigma2,
		FLOAT *g_diff2s,
		unsigned image_size,
		FLOAT sum_init,
		unsigned long orientation_num,
		unsigned long translation_num,
		unsigned long todo_blocks,
		unsigned long *d_rot_idx,
		unsigned long *d_trans_idx,
		unsigned long *d_job_idx,
		unsigned long *d_job_num
		)
{
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;

//    // Specialize BlockReduce for a 1D block of 128 threads on type FLOAT
//    typedef cub::BlockReduce<FLOAT, 128> BlockReduce;
//    // Allocate shared memory for BlockReduce
//    __shared__ typename BlockReduce::TempStorage temp_storage;

	int pixel;
	FLOAT ref_real;
	FLOAT ref_imag;

	__shared__ FLOAT s[BLOCK_SIZE*PROJDIFF_CHUNK_SIZE]; //We MAY have to do up to PROJDIFF_CHUNK_SIZE translations in each block
	__shared__ FLOAT s_outs[PROJDIFF_CHUNK_SIZE];
	// inside the padded 2D orientation gri
	if( bid < todo_blocks ) // we only need to make
	{
		unsigned trans_num   = d_job_num[bid]; //how many transes we have for this rot
		for (int itrans=0; itrans<trans_num; itrans++)
		{
			s[itrans*BLOCK_SIZE+tid] = 0.0f;
		}
		// index of comparison
		unsigned long int ix = d_rot_idx[d_job_idx[bid]];
		unsigned long int iy;
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  ));

		for (unsigned pass = 0; pass < pass_num; pass++) // finish an entire ref image each block
		{
			pixel = (pass * BLOCK_SIZE) + tid;

			if(pixel < image_size)
			{
				projector.project3Dmodel(
						pixel,
						__ldg(&g_eulers[ix*9  ]), __ldg(&g_eulers[ix*9+1]),
						__ldg(&g_eulers[ix*9+3]), __ldg(&g_eulers[ix*9+4]),
						__ldg(&g_eulers[ix*9+6]), __ldg(&g_eulers[ix*9+7]),
						ref_real, ref_imag);

				FLOAT diff_real;
				FLOAT diff_imag;
				for (int itrans=0; itrans<trans_num; itrans++) // finish all translations in each partial pass
				{
					iy=d_trans_idx[d_job_idx[bid]]+itrans;
					unsigned long img_start(iy * image_size);
					unsigned long img_pixel_idx = img_start + pixel;
					diff_real =  ref_real - __ldg(&g_imgs_real[img_pixel_idx]); // TODO  Put g_img_* in texture (in such a way that fetching of next image might hit in cache)
					diff_imag =  ref_imag - __ldg(&g_imgs_imag[img_pixel_idx]);
					s[itrans*BLOCK_SIZE + tid] += (diff_real * diff_real + diff_imag * diff_imag) * 0.5f * __ldg(&g_Minvsigma2[pixel]);
				}
				__syncthreads();
			}
		}
		for(int j=(BLOCK_SIZE/2); j>0; j/=2)
		{
			if(tid<j)
			{
				for (int itrans=0; itrans<trans_num; itrans++) // finish all translations in each partial pass
				{
					s[itrans*BLOCK_SIZE+tid] += s[itrans*BLOCK_SIZE+tid+j];
				}
			}
			__syncthreads();
		}
		if (tid < trans_num)
		{
			s_outs[tid]=s[tid*BLOCK_SIZE]+sum_init;
		}
		if (tid < trans_num)
		{
			iy=d_job_idx[bid]+tid;
			g_diff2s[iy] = s_outs[tid];
		}
	}
}

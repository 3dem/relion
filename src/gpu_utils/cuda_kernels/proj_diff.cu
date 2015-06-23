#include "src/gpu_utils/cuda_kernels/proj_diff.cuh"
//#include <external/cub-1.4.1/cub/block/block_reduce.cuh>
#include <vector>
#include <iostream>

#ifndef CUDA_DOUBLE_PRECISION

__global__ void cuda_kernel_PAV_TTI_D2( FLOAT *g_eulers,
		                                FLOAT *g_imgs_real,
		                                FLOAT *g_imgs_imag,
										cudaTextureObject_t texModel_real,
										cudaTextureObject_t texModel_imag,
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
										unsigned long *d_job_num,
										int my_r_max,
										int max_r2,
										int img_x,
										int img_y,
										int mdl_init_y,
										int mdl_init_z,
										float padding_factor
										)
{
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;

//    // Specialize BlockReduce for a 1D block of 128 threads on type FLOAT
//    typedef cub::BlockReduce<FLOAT, 128> BlockReduce;
//    // Allocate shared memory for BlockReduce
//    __shared__ typename BlockReduce::TempStorage temp_storage;

	FLOAT xp, yp, zp;
	long int r2;
	int pixel;
	bool is_neg_x;
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
			if(pixel<image_size)
			{
				int x = pixel % img_x;
				int y = (int)floorf( (float)pixel / (float)img_x);

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
					xp = (__ldg(&g_eulers[ix*9])   * x + __ldg(&g_eulers[ix*9+1]) * y ) * padding_factor;  // FIXME: xp,yp,zp has has accuracy loss
					yp = (__ldg(&g_eulers[ix*9+3]) * x + __ldg(&g_eulers[ix*9+4]) * y ) * padding_factor;  // compared to CPU-based projection. This
					zp = (__ldg(&g_eulers[ix*9+6]) * x + __ldg(&g_eulers[ix*9+7]) * y ) * padding_factor;  // propagates to dx00, dx10, and so on.
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
//				printf(" diffs = %f, %f \n",ref_real,img_pixel_idx);
//				printf(" diffs = %i, %i ,%i \n",x,y);
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
#elif !defined(CUDA_DOUBLE_PRECISION)
__global__ void cuda_kernel_PAV_TTE_D2( FLOAT *g_eulers,
		                                FLOAT *g_imgs_real,
		                                FLOAT *g_imgs_imag,
										cudaTextureObject_t texModel_real,
										cudaTextureObject_t texModel_imag,
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
										unsigned long *d_job_num,
										unsigned my_r_max,
										int max_r2,
										int min_r2_nn,
										long int img_x,
										long int img_y,
										long int mdl_init_y,
										long int mdl_init_z,
										float padding_factor
										)
{
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;

//    // Specialize BlockReduce for a 1D block of 128 threads on type FLOAT
//    typedef cub::BlockReduce<FLOAT, 128> BlockReduce;
//    // Allocate shared memory for BlockReduce
//    __shared__ typename BlockReduce::TempStorage temp_storage;
	CudaComplex d000, d001, d010, d011, d100, d101, d110, d111;
	CudaComplex dx00, dx01, dx10, dx11, dxy0, dxy1, val;
	FLOAT fx, fy, fz, xp, yp, zp;
	long int r2;
	int x0, x1, y0, y1, z0, z1, pixel;
	bool is_neg_x;
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
			if(pixel<image_size)
			{
				int x = pixel % img_x;
				int y = (int)floorf( (float)pixel / (float)img_x);

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
					xp = (__ldg(&g_eulers[ix*9])   * x + __ldg(&g_eulers[ix*9+1]) * y ) * padding_factor;  // FIXME: xp,yp,zp has has accuracy loss
					yp = (__ldg(&g_eulers[ix*9+3]) * x + __ldg(&g_eulers[ix*9+4]) * y ) * padding_factor;  // compared to CPU-based projection. This
					zp = (__ldg(&g_eulers[ix*9+6]) * x + __ldg(&g_eulers[ix*9+7]) * y ) * padding_factor;  // propagates to dx00, dx10, and so on.
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
					// Trilinear interpolation (with physical coords)
					// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
					// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
					x0 = floorf(xp);
					fx = xp - x0;
					x1 = x0 + 1;
					xp = fx + x0;

					y0 = floorf(yp);
					fy = yp - y0;
					y0 -= mdl_init_y;
					y1 = y0 + 1;
					yp -= mdl_init_y;

					z0 = floorf(zp);
					fz = zp - z0;
					z0 -= mdl_init_z;
					z1 = z0 + 1;
					zp -= mdl_init_z;

					d000.real = tex3D<FLOAT>(texModel_real,x0+0.5f,y0+0.5f,z0+0.5f);
					d001.real = tex3D<FLOAT>(texModel_real,x1+0.5f,y0+0.5f,z0+0.5f);
					d010.real = tex3D<FLOAT>(texModel_real,x0+0.5f,y1+0.5f,z0+0.5f);
					d011.real = tex3D<FLOAT>(texModel_real,x1+0.5f,y1+0.5f,z0+0.5f);
					d100.real = tex3D<FLOAT>(texModel_real,x0+0.5f,y0+0.5f,z1+0.5f);
					d101.real = tex3D<FLOAT>(texModel_real,x1+0.5f,y0+0.5f,z1+0.5f);
					d110.real = tex3D<FLOAT>(texModel_real,x0+0.5f,y1+0.5f,z1+0.5f);
					d111.real = tex3D<FLOAT>(texModel_real,x1+0.5f,y1+0.5f,z1+0.5f);

					d000.imag = tex3D<FLOAT>(texModel_imag,x0+0.5f,y0+0.5f,z0+0.5f);
					d001.imag = tex3D<FLOAT>(texModel_imag,x1+0.5f,y0+0.5f,z0+0.5f);
					d010.imag = tex3D<FLOAT>(texModel_imag,x0+0.5f,y1+0.5f,z0+0.5f);
					d011.imag = tex3D<FLOAT>(texModel_imag,x1+0.5f,y1+0.5f,z0+0.5f);
					d100.imag = tex3D<FLOAT>(texModel_imag,x0+0.5f,y0+0.5f,z1+0.5f);
					d101.imag = tex3D<FLOAT>(texModel_imag,x1+0.5f,y0+0.5f,z1+0.5f);
					d110.imag = tex3D<FLOAT>(texModel_imag,x0+0.5f,y1+0.5f,z1+0.5f);
					d111.imag = tex3D<FLOAT>(texModel_imag,x1+0.5f,y1+0.5f,z1+0.5f);

					// Set the interpolated value in the 2D output array
					dx00 = d000 + (d001 - d000)*fx;
					dx01 = d100 + (d101 - d100)*fx;
					dx10 = d010 + (d011 - d010)*fx;
					dx11 = d110 + (d111 - d110)*fx;

					dxy0 = dx00 + (dx10 - dx00)*fy;
					dxy1 = dx01 + (dx11 - dx01)*fy;

					val = dxy0 + (dxy1 - dxy0)*fz;
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
//				printf(" diffs = %f, %f \n",ref_real,img_pixel_idx);
//				printf(" diffs = %i, %i ,%i \n",x,y);
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
#else
__global__ void cuda_kernel_PAV_TGE_D2( FLOAT *g_eulers,
		                                FLOAT *g_imgs_real,
		                                FLOAT *g_imgs_imag,
		                                FLOAT *g_model_real,						// note DIFFERENT TYPE input compared to texture-utilising functions
		                                FLOAT *g_model_imag,						// note DIFFERENT TYPE input compared to texture-utilising functions
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
										unsigned long *d_job_num,
										unsigned my_r_max,
										int max_r2,
										long int img_x,
										long int img_y,
										long int mdl_init_y,
										long int mdl_init_z,
										long int mdl_size_x,						// note ADDITIONAL input compared to texture-utilising functions
										long int mdl_size_y,						// note ADDITIONAL input compared to texture-utilising functions
										float padding_factor
										)
{
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;

//    // Specialize BlockReduce for a 1D block of 128 threads on type FLOAT
//    typedef cub::BlockReduce<FLOAT, 128> BlockReduce;
//    // Allocate shared memory for BlockReduce
//    __shared__ typename BlockReduce::TempStorage temp_storage;
	CudaComplex d000, d001, d010, d011, d100, d101, d110, d111;
	CudaComplex dx00, dx01, dx10, dx11, dxy0, dxy1, val;
	FLOAT fx, fy, fz, xp, yp, zp;
	long int r2;
	int x0, x1, y0, y1, z0, z1, pixel;
	bool is_neg_x;
	long int mdl_size_yx=mdl_size_y* mdl_size_x;
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
			s[itrans*BLOCK_SIZE+tid] = 0.0;
		}
		// index of comparison
		unsigned long int ix = d_rot_idx[d_job_idx[bid]];
		unsigned long int iy;
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  ));

		for (unsigned pass = 0; pass < pass_num; pass++) // finish an entire ref image each block
		{
			pixel = (pass * BLOCK_SIZE) + tid;
			if(pixel<image_size)
			{
				int x = pixel % img_x;
				int y = (int)floorf( (float)pixel / (float)img_x);

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
					xp = (__ldg(&g_eulers[ix*9])   * x + __ldg(&g_eulers[ix*9+1]) * y ) * padding_factor;  // FIXME: xp,yp,zp has has accuracy loss
					yp = (__ldg(&g_eulers[ix*9+3]) * x + __ldg(&g_eulers[ix*9+4]) * y ) * padding_factor;  // compared to CPU-based projection. This
					zp = (__ldg(&g_eulers[ix*9+6]) * x + __ldg(&g_eulers[ix*9+7]) * y ) * padding_factor;  // propagates to dx00, dx10, and so on.
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
					// Trilinear interpolation (with physical coords)
					// Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
					// In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
					x0 = floorf(xp);
					fx = xp - x0;
					x1 = x0 + 1;
					xp = fx + x0;

					y0 = floorf(yp);
					fy = yp - y0;
					y0 -= mdl_init_y;
					y1 = y0 + 1;
					yp -= mdl_init_y;

					z0 = floorf(zp);
					fz = zp - z0;
					z0 -= mdl_init_z;
					z1 = z0 + 1;
					zp -= mdl_init_z;

					d000.real = g_model_real[z0*mdl_size_yx+y0*mdl_size_x+x0];
					d001.real = g_model_real[z0*mdl_size_yx+y0*mdl_size_x+x1];
					d010.real = g_model_real[z0*mdl_size_yx+y1*mdl_size_x+x0];
					d011.real = g_model_real[z0*mdl_size_yx+y1*mdl_size_x+x1];
					d100.real = g_model_real[z1*mdl_size_yx+y0*mdl_size_x+x0];
					d101.real = g_model_real[z1*mdl_size_yx+y0*mdl_size_x+x1];
					d110.real = g_model_real[z1*mdl_size_yx+y1*mdl_size_x+x0];
					d111.real = g_model_real[z1*mdl_size_yx+y1*mdl_size_x+x1];

					d000.imag = g_model_imag[z0*mdl_size_yx+y0*mdl_size_x+x0];
					d001.imag = g_model_imag[z0*mdl_size_yx+y0*mdl_size_x+x1];
					d010.imag = g_model_imag[z0*mdl_size_yx+y1*mdl_size_x+x0];
					d011.imag = g_model_imag[z0*mdl_size_yx+y1*mdl_size_x+x1];
					d100.imag = g_model_imag[z1*mdl_size_yx+y0*mdl_size_x+x0];
					d101.imag = g_model_imag[z1*mdl_size_yx+y0*mdl_size_x+x1];
					d110.imag = g_model_imag[z1*mdl_size_yx+y1*mdl_size_x+x0];
					d111.imag = g_model_imag[z1*mdl_size_yx+y1*mdl_size_x+x1];

					// Set the interpolated value in the 2D output array
					dx00 = d000 + (d001 - d000)*fx;
					dx01 = d100 + (d101 - d100)*fx;
					dx10 = d010 + (d011 - d010)*fx;
					dx11 = d110 + (d111 - d110)*fx;
					//-----------------------------
					dxy0 = dx00 + (dx10 - dx00)*fy;
					dxy1 = dx01 + (dx11 - dx01)*fy;
					//-----------------------------
					val = dxy0 + (dxy1 - dxy0)*fz;
					//-----------------------------
					ref_real = val.real;
					ref_imag = val.imag;

					if (is_neg_x)
					{
						ref_imag = -ref_imag;
					}
				}
				else
				{
					ref_real=0.0;
					ref_imag=0.0;
				}

				FLOAT diff_real;
				FLOAT diff_imag;
				for (int itrans=0; itrans<trans_num; itrans++) // finish all translations in each partial pass
				{
					iy=d_trans_idx[d_job_idx[bid]]+itrans;
					unsigned long img_start(iy * image_size);
					unsigned long img_pixel_idx = img_start + pixel;
					diff_real =  ref_real - __ldg(&g_imgs_real[img_pixel_idx]); // TODO  Put g_img_* in texture (in such a way that fetching of next image might hit in cache)
					diff_imag =  ref_imag - __ldg(&g_imgs_imag[img_pixel_idx]);
					s[itrans*BLOCK_SIZE + tid] += (diff_real * diff_real + diff_imag * diff_imag) * 0.5 * __ldg(&g_Minvsigma2[pixel]);
				}
				__syncthreads();
//				printf(" diffs = %f, %f \n",ref_real,img_pixel_idx);
//				printf(" diffs = %i, %i ,%i \n",x,y);
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
#endif // !defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)

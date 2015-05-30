#include "src/gpu_utils/cuda_kernels/projection.cuh"
#include <vector>
#include <iostream>

// uses global memory and explicit interpolation = can do double precision.
__global__ void cuda_kernel_PAV_TGE( FLOAT *g_model_real,
													FLOAT *g_model_imag,
													FLOAT *g_eulers,
													FLOAT *g_Frefs_real,
													FLOAT *g_Frefs_imag,
													int my_r_max,
													int max_r2,
													int min_r2_nn,
													int image_size,
													int orientation_num,
													int XSIZE_img,
													int YSIZE_img,
													int XSIZE_mdl,
													int YSIZE_mdl,
													int STARTINGY_mdl,
													int STARTINGZ_mdl
												   	   )
{
	FLOAT fx, fy, fz, xp, yp, zp;
	int x0, x1, y0, y1, z0, z1; //y2;
	long int r2;
	int pixel;
	int YXSIZE_mdl = YSIZE_mdl * XSIZE_mdl;
	bool is_neg_x;
	CudaComplex d000, d001, d010, d011, d100, d101, d110, d111;
	CudaComplex dx00, dx01, dx10, dx11, dxy0, dxy1, val;
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;
	// inside the padded 2D orientation grid
	if( bid < orientation_num ) // we only need to make
	{
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  ));
		long ref_pixel = bid*(image_size);
		for (unsigned pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
		{
			pixel = (pass * BLOCK_SIZE) + tid;
			if(pixel<image_size)
			{
				int x = pixel % XSIZE_img;
				int y = (int)floorf( (float)pixel / (float)XSIZE_img);

				// Dont search beyond square with side max_r
				if (y > my_r_max)
				{
					if (y >= YSIZE_img - my_r_max)
						y = y - YSIZE_img ;
					else
						x=r2;
				}

				r2 = x*x + y*y;
				if (r2 <= max_r2)
				{
					xp = g_eulers[bid*9]   * x + g_eulers[bid*9+1] * y;  // FIXME: xp,yp,zp has has accuracy loss
					yp = g_eulers[bid*9+3] * x + g_eulers[bid*9+4] * y;  // compared to CPU-based projection. This
					zp = g_eulers[bid*9+6] * x + g_eulers[bid*9+7] * y;  // propagates to dx00, dx10, and so on.
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
					y0 -=  STARTINGY_mdl;
					y1 = y0 + 1;

					z0 = floorf(zp);
					fz = zp - z0;
					z0 -= STARTINGZ_mdl;
					z1 = z0 + 1;

					d000.real = g_model_real[z0*YXSIZE_mdl+y0*XSIZE_mdl+x0];
					d001.real = g_model_real[z0*YXSIZE_mdl+y0*XSIZE_mdl+x1];
					d010.real = g_model_real[z0*YXSIZE_mdl+y1*XSIZE_mdl+x0];
					d011.real = g_model_real[z0*YXSIZE_mdl+y1*XSIZE_mdl+x1];
					d100.real = g_model_real[z1*YXSIZE_mdl+y0*XSIZE_mdl+x0];
					d101.real = g_model_real[z1*YXSIZE_mdl+y0*XSIZE_mdl+x1];
					d110.real = g_model_real[z1*YXSIZE_mdl+y1*XSIZE_mdl+x0];
					d111.real = g_model_real[z1*YXSIZE_mdl+y1*XSIZE_mdl+x1];

					d000.imag = g_model_imag[z0*YXSIZE_mdl+y0*XSIZE_mdl+x0];
					d001.imag = g_model_imag[z0*YXSIZE_mdl+y0*XSIZE_mdl+x1];
					d010.imag = g_model_imag[z0*YXSIZE_mdl+y1*XSIZE_mdl+x0];
					d011.imag = g_model_imag[z0*YXSIZE_mdl+y1*XSIZE_mdl+x1];
					d100.imag = g_model_imag[z1*YXSIZE_mdl+y0*XSIZE_mdl+x0];
					d101.imag = g_model_imag[z1*YXSIZE_mdl+y0*XSIZE_mdl+x1];
					d110.imag = g_model_imag[z1*YXSIZE_mdl+y1*XSIZE_mdl+x0];
					d111.imag = g_model_imag[z1*YXSIZE_mdl+y1*XSIZE_mdl+x1];

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
						val.imag = -val.imag;
					}

				}
				else
				{
					val.real=0.0f;
					val.imag=0.0f;
				}
				g_Frefs_real[ref_pixel+ pixel] = val.real;
				g_Frefs_imag[ref_pixel+ pixel] = val.imag;
			}
		}
	}
}

#if !defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)

// uses texture memory and implicit (texture) interpolation = requires float precision.
__global__ void cuda_kernel_PAV_TTI( FLOAT *g_eulers,
														 FLOAT *g_Frefs_real,
														 FLOAT *g_Frefs_imag,
														 cudaTextureObject_t texModel_real,
														 cudaTextureObject_t texModel_imag,
														 int my_r_max,
														 int max_r2,
														 int min_r2_nn,
														 int image_size,
														 int orientation_num,
														 int XSIZE_img,
														 int YSIZE_img,
														 int STARTINGY_mdl,
														 int STARTINGZ_mdl)
{
	FLOAT xp, yp, zp;
	long int r2;
	int pixel;
	bool is_neg_x;
	CudaComplex val;
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;
	// inside the padded 2D orientation grid
	if( bid < orientation_num )
	{
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  ));
		long ref_pixel = bid*(image_size);
		for (unsigned pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
		{
			pixel = (pass * BLOCK_SIZE) + tid;
			if(pixel<image_size)
			{
				int x = pixel % XSIZE_img;
				int y = (int)floorf( (float)pixel / (float)XSIZE_img);

				// Dont search beyond square with side max_r
				if (y > my_r_max)
				{
					if (y >= YSIZE_img - my_r_max)
						y = y - YSIZE_img ;
					else
						x=r2;
				}

				r2 = x*x + y*y;
				if (r2 <= max_r2)
				{
					xp = __ldg(&g_eulers[bid*9])   * x + __ldg(&g_eulers[bid*9+1]) * y;  // FIXME: xp,yp,zp has has accuracy loss
					yp = __ldg(&g_eulers[bid*9+3]) * x + __ldg(&g_eulers[bid*9+4]) * y;  // compared to CPU-based projection. This
					zp = __ldg(&g_eulers[bid*9+6]) * x + __ldg(&g_eulers[bid*9+7]) * y;  // propagates to dx00, dx10, and so on.
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
					yp -= STARTINGY_mdl;
					zp -= STARTINGZ_mdl;

					val.real=tex3D<FLOAT>(texModel_real,xp+0.5f,yp+0.5f,zp+0.5f);
					val.imag=tex3D<FLOAT>(texModel_imag,xp+0.5f,yp+0.5f,zp+0.5f);


					if (is_neg_x)
					{
						val.imag = -val.imag;
					}

				}
				else
				{
					val.real=0.0f;
					val.imag=0.0f;
				}
				g_Frefs_real[ref_pixel+ pixel] = val.real;
				g_Frefs_imag[ref_pixel+ pixel] = val.imag;
			}
		}
	}
}

#elif !defined(CUDA_DOUBLE_PRECISION)

// uses texture memory and explicit interpolation = requires float precision.
__global__ void cuda_kernel_PAV_TTE( FLOAT *g_eulers,
														 FLOAT *g_Frefs_real,
														 FLOAT *g_Frefs_imag,
														 cudaTextureObject_t texModel_real,
														 cudaTextureObject_t texModel_imag,
														 int my_r_max,
														 int max_r2,
														 int min_r2_nn,
														 int image_size,
														 int orientation_num,
														 int XSIZE_img,
														 int YSIZE_img,
														 int STARTINGY_mdl,
														 int STARTINGZ_mdl)
{
	FLOAT fx, fy, fz, xp, yp, zp;
	int x0, x1, y0, y1, z0, z1; //y2;
	long int r2;
	int pixel;
	bool is_neg_x;
	CudaComplex d000, d001, d010, d011, d100, d101, d110, d111;
	CudaComplex dx00, dx01, dx10, dx11, dxy0, dxy1, val;
	int bid = blockIdx.y * gridDim.x + blockIdx.x;
	int tid = threadIdx.x;
	// inside the padded 2D orientation grid
	if( bid < orientation_num ) // we only need to make
	{
		unsigned pass_num(ceilf(   ((float)image_size) / (float)BLOCK_SIZE  ));
		long ref_pixel = bid*(image_size);
		for (unsigned pass = 0; pass < pass_num; pass++) // finish a reference proj in each block
		{
			pixel = (pass * BLOCK_SIZE) + tid;
			if(pixel<image_size)
			{
				int x = pixel % XSIZE_img;
				int y = (int)floorf( (float)pixel / (float)XSIZE_img);

				// Dont search beyond square with side max_r
				if (y > my_r_max)
				{
					if (y >= YSIZE_img - my_r_max)
						y = y - YSIZE_img ;
					else
						x=r2;
				}

				r2 = x*x + y*y;
				if (r2 <= max_r2)
				{
					xp = g_eulers[bid*9]   * x + g_eulers[bid*9+1] * y;  // FIXME: xp,yp,zp has has accuracy loss
					yp = g_eulers[bid*9+3] * x + g_eulers[bid*9+4] * y;  // compared to CPU-based projection. This
					zp = g_eulers[bid*9+6] * x + g_eulers[bid*9+7] * y;  // propagates to dx00, dx10, and so on.
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
					y0 -=  STARTINGY_mdl;
					y1 = y0 + 1;
					yp -= STARTINGY_mdl;

					z0 = floorf(zp);
					fz = zp - z0;
					z0 -= STARTINGZ_mdl;
					z1 = z0 + 1;
					zp -= STARTINGZ_mdl;

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
						val.imag = -val.imag;
					}

				}
				else
				{
					val.real=0.0f;
					val.imag=0.0f;
				}
				g_Frefs_real[ref_pixel+ pixel] = val.real;
				g_Frefs_imag[ref_pixel+ pixel] = val.imag;
			}
		}
	}
}

#endif //!defined(CUDA_DOUBLE_PRECISION) && defined(USE_TEXINTERP)

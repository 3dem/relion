#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_settings.h"



namespace CudaKernels
{

template <typename T>
__global__ void cuda_kernel_translate2D(	T * g_image_in,
											T * g_image_out,
											int image_size,
											int xdim,
											int ydim,
											int dx,
											int dy)
{
	int tid = threadIdx.x;
	int bid =  blockIdx.x;

	int x,y,xp,yp;
	int pixel=tid + bid*BLOCK_SIZE;
	int new_pixel;

	if(pixel<image_size)
	{
		x = pixel % xdim;
		y = (pixel-x) / (xdim);

		xp = x + dx;
		yp = y + dy;

		if( yp>=0 && xp>=0 && yp<ydim && xp<xdim)
		{
			new_pixel = yp*xdim + xp;
			if(new_pixel>=0 && new_pixel<image_size) // if displacement is negative, new_pixel could be less than 0
				g_image_out[new_pixel] = g_image_in[pixel];
		}
	}
}

template <typename T>
__global__ void cuda_kernel_translate3D(	T * g_image_in,
											T * g_image_out,
											int image_size,
											int xdim,
											int ydim,
											int zdim,
											int dx,
											int dy,
											int dz)
{
	int tid = threadIdx.x;
	int bid =  blockIdx.x;

	int x,y,z,xp,yp,zp,xy;
	int voxel=tid + bid*BLOCK_SIZE;
	int new_voxel;

	int xydim = xdim*ydim;

	if(voxel<image_size)
	{
		z =  voxel / xydim;
		zp = z + dz;

		xy = voxel % xydim;
		y =  xy / xdim;
		yp = y + dy;

		x =  xy % xdim;
		xp = x + dx;

		if( zp>=0 && yp>=0 && xp>=0 && zp<zdim && yp<ydim && xp<xdim)
		{
			new_voxel = zp*xydim +  yp*xdim + xp;
			if(new_voxel>=0 && new_voxel<image_size) // if displacement is negative, new_pixel could be less than 0
				g_image_out[new_voxel] = g_image_in[voxel];
		}
	}
}

} // namespace


template <typename T>
__global__ void cuda_kernel_multi( T *A,
								   T *OUT,
								   T S,
		  	  	  	  	  	  	   int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel] = A[pixel]*S;
}

namespace CudaKernels
{
template <typename T>
__global__ void cuda_kernel_multi(
		T *A,
		T S,
		int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		A[pixel] = A[pixel]*S;
}
}

template <typename T>
__global__ void cuda_kernel_multi( T *A,
								   T *B,
								   T *OUT,
								   T S,
		  	  	  	  	  	  	   int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel] = A[pixel]*B[pixel]*S;
}

__global__ void cuda_kernel_allweights_to_mweights(
		unsigned long * d_iorient,
		XFLOAT * d_allweights,
		XFLOAT * d_mweights,
		unsigned long orientation_num,
		unsigned long translation_num,
        int block_size
		)
{
	size_t idx = blockIdx.x * block_size + threadIdx.x;
	if (idx < orientation_num*translation_num)
		d_mweights[d_iorient[idx/translation_num] * translation_num + idx%translation_num] =
				d_allweights[idx/translation_num * translation_num + idx%translation_num];
}
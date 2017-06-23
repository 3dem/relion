#ifndef ACC_UTILITIES_H_
#define ACC_UTILITIES_H_

#include "src/acc/acc_ptr.h"
#include "src/acc/data_types.h"
#ifdef CUDA
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#else
#include "src/acc/cpu/cpu_kernels/helper.h"
#endif

namespace AccUtilities
{
	template <typename T>
	static void multiply(AccDataTypes::Image<T> &ptr, T value)
	{
#ifdef CUDA
		int BSZ = ( (int) ceilf(( float)ptr.getSize() /(float)BLOCK_SIZE));
		CudaKernels::cuda_kernel_multi<T><<<BSZ,BLOCK_SIZE,0,ptr.getStream()>>>(
			ptr(),
			value,
			ptr.getSize());
#else
		CpuKernels::cpu_kernel_multi<T>(
		ptr(),
		value,
		ptr.getSize());
#endif
	}
	
	template <typename T>
	static void multiply(T *array, T value, size_t size)
	{
#ifdef CUDA
// TODO - How do we figure this out too?
//		int BSZ = ( (int) ceilf(( float)ptr.getSize() /(float)BLOCK_SIZE));
// TODO - how do we get streams?
//		CudaKernels::cuda_kernel_multi<T><<<BSZ,BLOCK_SIZE,0,ptr.getStream()>>>(
	    CudaKernels::cuda_kernel_multi<T><<<128,BLOCK_SIZE,0>>>(
			array,
			value,
			size);
#else
		CpuKernels::cpu_kernel_multi<T>(
		array,
		value,
		size);
#endif
	}

	template <typename T>
	static void translate(
		AccDataTypes::Image<T> &in,
		AccDataTypes::Image<T> &out,
		int dx,
		int dy,
		int dz=0)
	{
#ifdef CUDA
	int BSZ = ( (int) ceilf(( float)in.getxyz() /(float)BLOCK_SIZE));

	if (in.is3D())
	{
		CudaKernels::cuda_kernel_translate3D<T><<<BSZ,BLOCK_SIZE,0,in.getStream()>>>(
			in(),
			out(),
			in.getxyz(),
			in.getx(),
			in.gety(),
			in.getz(),
			dx,
			dy,
			dz);
	}
	else
	{
		CudaKernels::cuda_kernel_translate2D<T><<<BSZ,BLOCK_SIZE,0,in.getStream()>>>(
			in(),
			out(),
			in.getxyz(),
			in.getx(),
			in.gety(),
			dx,
			dy);
	}
#else
	if (in.is3D())
	{
		CpuKernels::cpu_translate3D<T>(
			in(),
			out(),
			in.getxyz(),
			in.getx(),
			in.gety(),
			in.getz(),
			dx,
			dy,
			dz);
	}
	else
	{
		CpuKernels::cpu_translate2D<T>(
			in(),
			out(),
			in.getxyz(),
			in.getx(),
			in.gety(),
			dx,
			dy);
	}
#endif
	}
	
	
static void softMaskBackgroundValue(
		XFLOAT *vol,
		Image<RFLOAT> &img,
		bool     do_Mnoise,
		XFLOAT   radius,
		XFLOAT   radius_p,
		XFLOAT   cosine_width,
		XFLOAT  *g_sum,
		XFLOAT  *g_sum_bg,
		int inblock_dim, 
		int inblock_size)
	{
#ifdef CUDA
		dim3 block_dim = inblock_dim;
		
		cuda_kernel_softMaskBackgroundValue<<<block_dim,inblock_size>>>(	vol,
																	img().nzyxdim,
																	img.data.xdim,
																	img.data.ydim,
																	img.data.zdim,
																	img.data.xdim/2,
																	img.data.ydim/2,
																	img.data.zdim/2, //unused
																	do_Mnoise,
																	radius,
																	radius_p,
																	cosine_width,
																	g_sum,
																	g_sum_bg);
#else
		int block_dim = inblock_dim;
		CpuKernels::SoftMaskBackgroundValue(block_dim, inblock_size,
			vol,
			img().nzyxdim,
			img.data.xdim,
			img.data.ydim,
			img.data.zdim,
			img.data.xdim/2,
			img.data.ydim/2,
			img.data.zdim/2, //unused
			do_Mnoise,
			radius,
			radius_p,
			cosine_width,
			g_sum,
			g_sum_bg);
#endif
	}
};  // namespace 


#endif //ACC_UTILITIES_H_


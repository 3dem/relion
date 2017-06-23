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
	static
	void
	multiply(AccDataTypes::Image<T> &ptr, T value)
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
	static
	void
	translate(
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
};


#endif //ACC_UTILITIES_H_


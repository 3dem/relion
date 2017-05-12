#ifndef ACC_UTILITIES_H_
#define ACC_UTILITIES_H_

#include "src/acc/acc_ptr.h"
#include "src/acc/cuda/utilities.cu"
#include "src/acc/cpu/utilities.cpp"

namespace AccUtilities
{
	template <typename T, int Acc>
	static
	void
	multiply(AccPtr<T,Acc> &ptr, T value)
	{
		if (Acc == ACC_CPU)
		{
			CpuKernels::multiply(
				~ptr,
				value,
				ptr.getSize());
		}
		else
		{
			int BSZ = ( (int) ceilf(( float)ptr.getSize() /(float)BLOCK_SIZE));
			CudaKernels::multiply<<<BSZ,BLOCK_SIZE,0,ptr.getStream()>>>(
				~ptr,
				value,
				ptr.getSize());
		}
	}

	template <typename T, int Acc>
	static
	void
	translate2D(
			AccPtr<T,Acc> &in,
			AccPtr<T,Acc> &out,
			int xdim,
			int ydim,
			int dx,
			int dy)
	{
		if (Acc == ACC_CPU)
		{
			CpuKernels::translate2D(
				~in,
				~out,
				in.getSize(),
				xdim,
				ydim,
				dx,
				dy);
		}
		else
		{
			int BSZ = ( (int) ceilf(( float)in.getSize() /(float)BLOCK_SIZE));
			CudaKernels::translate2D<<<BSZ,BLOCK_SIZE,0,in.getStream()>>>(
				~in,
				~out,
				in.getSize(),
				xdim,
				ydim,
				dx,
				dy);
		}
	}

	template <typename T, int Acc>
	static
	void
	translate3D(
			AccPtr<T,Acc> &in,
			AccPtr<T,Acc> &out,
			int xdim,
			int ydim,
			int zdim,
			int dx,
			int dy,
			int dz)
	{
		if (Acc == ACC_CPU)
		{
			CpuKernels::translate3D(
				~in,
				~out,
				in.getSize(),
				xdim,
				ydim,
				zdim,
				dx,
				dy,
				dz);
		}
		else
		{
			int BSZ = ( (int) ceilf(( float)in.getSize() /(float)BLOCK_SIZE));
			CudaKernels::translate3D<<<BSZ,BLOCK_SIZE,0,in.getStream()>>>(
				~in,
				~out,
				in.getSize(),
				xdim,
				ydim,
				zdim,
				dx,
				dy,
				dz);
		}
	}
};


#endif //ACC_UTILITIES_H_


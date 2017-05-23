#ifndef ACC_UTILITIES_H_
#define ACC_UTILITIES_H_

#include "src/acc/acc_ptr.h"
#include "src/acc/data_types.h"
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
				ptr(),
				value,
				ptr.getSize());
		}
		else
		{
			int BSZ = ( (int) ceilf(( float)ptr.getSize() /(float)BLOCK_SIZE));
			CudaKernels::multiply<<<BSZ,BLOCK_SIZE,0,ptr.getStream()>>>(
				ptr(),
				value,
				ptr.getSize());
		}
	}

	template <typename T, int Acc>
	static
	void
	translate(
		AccDataTypes::Image<T,Acc> &in,
		AccDataTypes::Image<T,Acc> &out,
		int dx,
		int dy,
		int dz=0)
	{

		if (Acc == ACC_CPU)
		{
			if (in.is3D())
			{
				CpuKernels::translate3D(
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
				CpuKernels::translate2D(
					in(),
					out(),
					in.getxyz(),
					in.getx(),
					in.gety(),
					dx,
					dy);
			}
		}
		else
		{
			int BSZ = ( (int) ceilf(( float)in.getxyz() /(float)BLOCK_SIZE));

			if (in.is3D())
			{
				CudaKernels::translate3D<<<BSZ,BLOCK_SIZE,0,in.getStream()>>>(
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
				CudaKernels::translate2D<<<BSZ,BLOCK_SIZE,0,in.getStream()>>>(
					in(),
					out(),
					in.getxyz(),
					in.getx(),
					in.gety(),
					dx,
					dy);
			}
		}
	}
};


#endif //ACC_UTILITIES_H_


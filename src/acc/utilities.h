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
					&ptr[0],
					value,
					ptr.getSize());
		}
		else
		{
			int BSZ = ( (int) ceilf(( float)ptr.getSize() /(float)BLOCK_SIZE));
			CudaKernels::multiply<<<BSZ,BLOCK_SIZE>>>(
				~ptr,
				value,
				ptr.getSize());
		}
	}
};


#endif //ACC_UTILITIES_H_


#ifndef CUDA_DEVICE_UTILS_CUH_
#define CUDA_DEVICE_UTILS_CUH_

#include <cuda_runtime.h>

#ifdef CUDA_DOUBLE_PRECISION
__device__ inline double cuda_atomic_add(double* address, double val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do
	{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
	}
	while (assumed != old); // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	return __longlong_as_double(old);
}
#else
__device__ inline void cuda_atomic_add(float* address, float value)
{
  atomicAdd(address,value);
}
#endif

#endif

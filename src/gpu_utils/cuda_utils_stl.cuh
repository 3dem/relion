#ifndef CUDA_UTILS_STL_CUH_
#define CUDA_UTILS_STL_CUH_

#include <cuda_runtime.h>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <stdio.h>
#include <signal.h>
#include <vector>
#include <algorithm>
#include "src/gpu_utils/cub/device/device_radix_sort.cuh"
#include "src/gpu_utils/cub/device/device_reduce.cuh"

#ifdef CUDA_DOUBLE_PRECISION
#define XFLOAT double
__device__ inline XFLOAT cuda_atomic_add(double* address, double val)
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
#define XFLOAT float
__device__ inline void cuda_atomic_add(float* address, float value)
{
  atomicAdd(address,value);
}
#endif

template <typename T>
static T getMaxOnDevice(CudaGlobalPtr<T> &ptr)
{
#ifdef DEBUG_CUDA
if (ptr.size == 0)
	printf("DEBUG_ERROR: getMinOnDevice called with pointer of zero size.\n");
if (ptr.d_ptr == NULL)
	printf("DEBUG_ERROR: getMinOnDevice called with null device pointer.\n");
if (ptr.getAllocator() == NULL)
	printf("DEBUG_ERROR: getMinOnDevice called with null allocator.\n");
#endif
	ptr.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(0));
	return (T)*std::max_element(ptr.h_ptr,ptr.h_ptr + ptr.size);
}

template <typename T>
static std::pair<int, T> getArgMaxOnDevice(CudaGlobalPtr<T> &ptr)
{
#ifdef DEBUG_CUDA
if (ptr.size == 0)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with pointer of zero size.\n");
if (ptr.d_ptr == NULL)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with null device pointer.\n");
if (ptr.getAllocator() == NULL)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with null allocator.\n");
#endif
	std::pair<int, T> max_pair;
	ptr.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(0));
	max_pair.first  = std::distance(ptr.h_ptr,std::max_element(ptr.h_ptr,ptr.h_ptr + ptr.size));
	max_pair.second = ptr.h_ptr[max_pair.first];

	return max_pair;
}

template <typename T>
static T getMinOnDevice(CudaGlobalPtr<T> &ptr)
{
#ifdef DEBUG_CUDA
if (ptr.size == 0)
	printf("DEBUG_ERROR: getMinOnDevice called with pointer of zero size.\n");
if (ptr.d_ptr == NULL)
	printf("DEBUG_ERROR: getMinOnDevice called with null device pointer.\n");
if (ptr.getAllocator() == NULL)
	printf("DEBUG_ERROR: getMinOnDevice called with null allocator.\n");
#endif
	ptr.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(0));
	return (T)*std::min_element(ptr.h_ptr,ptr.h_ptr + ptr.size);
}


template <typename T>
static std::pair<int, T> getArgMinOnDevice(CudaGlobalPtr<T> &ptr)
{
#ifdef DEBUG_CUDA
if (ptr.size == 0)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with pointer of zero size.\n");
if (ptr.d_ptr == NULL)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with null device pointer.\n");
if (ptr.getAllocator() == NULL)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with null allocator.\n");
#endif
	std::pair<int, T> min_pair;
	ptr.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(0));
	min_pair.first  = std::distance(ptr.h_ptr,std::min_element(ptr.h_ptr,ptr.h_ptr + ptr.size));
	min_pair.second = ptr.h_ptr[min_pair.first];

	return min_pair;
}

template <typename T>
static T getSumOnDevice(CudaGlobalPtr<T> &ptr)
{
#ifdef DEBUG_CUDA
if (ptr.size == 0)
	printf("DEBUG_ERROR: getSumOnDevice called with pointer of zero size.\n");
if (ptr.d_ptr == NULL)
	printf("DEBUG_ERROR: getSumOnDevice called with null device pointer.\n");
if (ptr.getAllocator() == NULL)
	printf("DEBUG_ERROR: getSumOnDevice called with null allocator.\n");
#endif
	ptr.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(0));
	T sum(0);
	for(long int i =0; i<ptr.size; i++)
		sum+=ptr.h_ptr[i];
	return sum;
}

template <typename T>
static void sortOnDevice(CudaGlobalPtr<T> &in, CudaGlobalPtr<T> &out)
{
#ifdef DEBUG_CUDA
if (in.size == 0 || out.size == 0)
	printf("DEBUG_ERROR: sortOnDevice called with pointer of zero size.\n");
if (in.d_ptr == NULL || out.d_ptr == NULL)
	printf("DEBUG_ERROR: sortOnDevice called with null device pointer.\n");
if (in.getAllocator() == NULL)
	printf("DEBUG_ERROR: sortOnDevice called with null allocator.\n");
#endif
	in.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(0));
	std::sort(in.h_ptr,in.h_ptr + in.size);
	for(long int i =0; i<in.size; i++)
		out.h_ptr[i]=in.h_ptr[i];
	out.cp_to_device();
	HANDLE_ERROR(cudaStreamSynchronize(0));
}

#endif


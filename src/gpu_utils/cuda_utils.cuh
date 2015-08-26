#ifndef CUDA_UTILS_CUH_
#define CUDA_UTILS_CUH_

#include <cuda_runtime.h>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <stdio.h>
#include <signal.h>
#include <vector>
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
static cub::KeyValuePair<int, T> getArgMaxOnDevice(CudaGlobalPtr<T> &ptr)
{
#ifdef DEBUG_CUDA
if (ptr.size == 0)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with pointer of zero size.\n");
if (ptr.d_ptr == NULL)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with null device pointer.\n");
if (ptr.getAllocator() == NULL)
	printf("DEBUG_WARNING: getArgMaxOnDevice called with null allocator.\n");
#endif
	CudaGlobalPtr<cub::KeyValuePair<int, T> >  max_pair(1, ptr.getAllocator());
	max_pair.device_alloc();
	size_t temp_storage_size = 0;

	HANDLE_ERROR(cub::DeviceReduce::ArgMax( NULL, temp_storage_size, ~ptr, ~max_pair, ptr.size));

	CudaCustomAllocator::Alloc* alloc = ptr.getAllocator()->alloc(temp_storage_size);

	HANDLE_ERROR(cub::DeviceReduce::ArgMax( alloc->getPtr(), temp_storage_size, ~ptr, ~max_pair, ptr.size));

	max_pair.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(ptr.getStream()));

	ptr.getAllocator()->free(alloc);

	return max_pair[0];
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
	CudaGlobalPtr<T >  min_val(1, ptr.getAllocator());
	min_val.device_alloc();
	size_t temp_storage_size = 0;

	HANDLE_ERROR(cub::DeviceReduce::Min( NULL, temp_storage_size, ~ptr, ~min_val, ptr.size));

	CudaCustomAllocator::Alloc* alloc = ptr.getAllocator()->alloc(temp_storage_size);

	HANDLE_ERROR(cub::DeviceReduce::Min( alloc->getPtr(), temp_storage_size, ~ptr, ~min_val, ptr.size));

	min_val.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(ptr.getStream()));

	ptr.getAllocator()->free(alloc);

	return min_val[0];
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
	CudaGlobalPtr<T >  val(1, ptr.getAllocator());
	val.device_alloc();
	size_t temp_storage_size = 0;

	HANDLE_ERROR(cub::DeviceReduce::Sum( NULL, temp_storage_size, ~ptr, ~val, ptr.size));

	CudaCustomAllocator::Alloc* alloc = ptr.getAllocator()->alloc(temp_storage_size);

	HANDLE_ERROR(cub::DeviceReduce::Sum( alloc->getPtr(), temp_storage_size, ~ptr, ~val, ptr.size));

	val.cp_to_host();
	HANDLE_ERROR(cudaStreamSynchronize(ptr.getStream()));

	ptr.getAllocator()->free(alloc);

	return val[0];
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
	size_t temp_storage_size = 0;

	HANDLE_ERROR(cub::DeviceRadixSort::SortKeys( NULL, temp_storage_size, ~in, ~out, in.size));

	CudaCustomAllocator::Alloc* alloc = in.getAllocator()->alloc(temp_storage_size);

	HANDLE_ERROR(cub::DeviceRadixSort::SortKeys( alloc->getPtr(), temp_storage_size, ~in, ~out, in.size));

	in.getAllocator()->free(alloc);
}

#endif

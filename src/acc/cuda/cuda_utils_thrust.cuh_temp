#ifndef CUDA_UTILS_THRUST_CUH_
#define CUDA_UTILS_THRUST_CUH_

#include <cuda_runtime.h>
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <stdio.h>
#include <signal.h>
#include <vector>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/extrema.h>
#include <thrust/device_vector.h>

template <typename T>
static T getMaxOnDevice(CudaGlobalPtr<T> &ptr)
{
	thrust::device_ptr<T> dp = thrust::device_pointer_cast(~ptr);
	thrust::device_ptr<T> pos = thrust::max_element(dp, dp + ptr.size);
	unsigned int pos_index = thrust::distance(dp, pos);
	T max_val;
	DEBUG_HANDLE_ERROR(cudaMemcpy(&max_val, &ptr.d_ptr[pos_index], sizeof(T), cudaMemcpyDeviceToHost));
	return max_val;
}

template <typename T>
static std::pair<int, T> getArgMaxOnDevice(CudaGlobalPtr<T> &ptr)
{
	std::pair<int, T> pair;
	thrust::device_ptr<T> dp = thrust::device_pointer_cast(~ptr);
	thrust::device_ptr<T> pos = thrust::max_element(dp, dp + ptr.size);
	pair.first = thrust::distance(dp, pos);
	DEBUG_HANDLE_ERROR(cudaMemcpy( &pair.second, &ptr.d_ptr[pair.first], sizeof(T), cudaMemcpyDeviceToHost));
	return pair;
}

template <typename T>
static T getMinOnDevice(CudaGlobalPtr<T> &ptr)
{
	thrust::device_ptr<T> dp = thrust::device_pointer_cast(~ptr);
	thrust::device_ptr<T> pos = thrust::min_element(dp, dp + ptr.size);
	unsigned int pos_index = thrust::distance(dp, pos);
	T min_val;
	DEBUG_HANDLE_ERROR(cudaMemcpy(&min_val, &ptr.d_ptr[pos_index], sizeof(T), cudaMemcpyDeviceToHost));
	return min_val;
}

template <typename T>
static std::pair<int, T> getArgMinOnDevice(CudaGlobalPtr<T> &ptr)
{
	std::pair<int, T> pair;
	thrust::device_ptr<T> dp = thrust::device_pointer_cast(~ptr);
	thrust::device_ptr<T> pos = thrust::min_element(dp, dp + ptr.size);
	pair.first = thrust::distance(dp, pos);
	DEBUG_HANDLE_ERROR(cudaMemcpy( &pair.second, &ptr.d_ptr[pair.first], sizeof(T), cudaMemcpyDeviceToHost));
	return pair;
}

template <typename T>
static T getSumOnDevice(CudaGlobalPtr<T> &ptr)
{
	thrust::device_ptr<T> dp = thrust::device_pointer_cast(~ptr);
	return thrust::reduce(dp, dp + ptr.size);
}

template <typename T>
static void sortOnDevice(CudaGlobalPtr<T> &in, CudaGlobalPtr<T> &out)
{
	//TODO Actually do sorting only on device instead of copying back and forth
	//Copy from "in" to "out" on device and do sorting there
//	DEBUG_HANDLE_ERROR(cudaMemcpy( ~out, ~in, in.size * sizeof(T), cudaMemcpyDeviceToDevice));
//	thrust::device_ptr<T> dp = thrust::device_pointer_cast(~out);
//	thrust::device_vector<T> dv(dp, dp + in.size);
//	thrust::sort(dv.begin(), dv.end() );

	T *h_vec = new T[in.size];
	DEBUG_HANDLE_ERROR(cudaMemcpy( h_vec, ~in, in.size * sizeof(T), cudaMemcpyDeviceToHost));
	thrust::sort(h_vec, h_vec + in.size);
	DEBUG_HANDLE_ERROR(cudaMemcpy( ~out, h_vec, in.size * sizeof(T), cudaMemcpyHostToDevice));
	delete [] h_vec;
}

#endif

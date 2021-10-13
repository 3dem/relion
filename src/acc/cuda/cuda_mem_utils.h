#ifndef CUDA_DEVICE_MEM_UTILS_H_
#define CUDA_DEVICE_MEM_UTILS_H_

#ifdef _CUDA_ENABLED
#include <cuda_runtime.h>
#include <curand.h>
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/custom_allocator.cuh"
#endif


#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "src/complex.h"

// Forward definition
template <typename T>  class AccPtr;

/**
 * Print cuda device memory info
 */
static void cudaPrintMemInfo()
{
	size_t free;
	size_t total;
	DEBUG_HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
	float free_hr(free/(1024.*1024.));
	float total_hr(total/(1024.*1024.));
    printf( "free %.2fMiB, total %.2fMiB, used %.2fMiB\n",
    		free_hr, total_hr, total_hr - free_hr);
}

template< typename T>
static inline
void cudaCpyHostToDevice( T *h_ptr, T *d_ptr, size_t size)
{
	DEBUG_HANDLE_ERROR(cudaMemcpy( d_ptr, h_ptr, size * sizeof(T), cudaMemcpyHostToDevice));
};

template< typename T>
static inline
void cudaCpyHostToDevice( T *h_ptr, T *d_ptr, size_t size, cudaStream_t &stream)
{
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( d_ptr, h_ptr, size * sizeof(T), cudaMemcpyHostToDevice, stream));
};

template< typename T>
static inline
void cudaCpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size)
{
	DEBUG_HANDLE_ERROR(cudaMemcpy( h_ptr, d_ptr, size * sizeof(T), cudaMemcpyDeviceToHost));
};

template< typename T>
static inline
void cudaCpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size, cudaStream_t &stream)
{
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( h_ptr, d_ptr, size * sizeof(T), cudaMemcpyDeviceToHost, stream));
};

template< typename T>
static inline
void cudaCpyDeviceToDevice( T *src, T *des, size_t size, cudaStream_t stream)
{
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( des, src, size * sizeof(T), cudaMemcpyDeviceToDevice, stream));
};

template< typename T>
static inline
void cudaMemInit( T *ptr, T value, size_t size)
{
	DEBUG_HANDLE_ERROR(cudaMemset( ptr, value, size * sizeof(T)));
};

template< typename T>
static inline
void cudaMemInit( T *ptr, T value, size_t size, cudaStream_t &stream)
{
	DEBUG_HANDLE_ERROR(cudaMemsetAsync( ptr, value, size * sizeof(T), stream));
};

#endif

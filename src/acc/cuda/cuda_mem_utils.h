#ifndef CUDA_DEVICE_MEM_UTILS_H_
#define CUDA_DEVICE_MEM_UTILS_H_

#ifdef CUDA
#include <cuda_runtime.h>
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/custom_allocator.cuh"
#endif


#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "src/parallel.h"

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

template <typename T>
class cudaStager
{
public:
	AccPtr<T> AllData;
	size_t size; // size of allocated host-space (AllData.getSize() dictates the amount of memory copied to/from the device)

	/*======================================================
				CONSTRUCTORS WITH ALLOCATORS
	======================================================*/

	inline
	cudaStager(CudaCustomAllocator *allocator):
		AllData(allocator), size(0)
	{};

	inline
	cudaStager(CudaCustomAllocator *allocator, size_t newSize):
		AllData(newSize,allocator), size(newSize)
	{
		AllData.setSize(0);
	};

	/*======================================================
				CONSTRUCTORS WITHOUT ALLOCATORS
	======================================================*/

	inline
	cudaStager():
		AllData(), size(0)
	{};

	inline
	cudaStager(size_t newSize):
		AllData(newSize), size(newSize)
	{
		AllData.setSize(0);
	};

public:

	void prepare_host()
	{

		if(size==0)
		{
			printf("trying to host-alloc a stager with size=0");
			CRITICAL(ERR_STAGEMEM);
		}
		size_t temp_size=AllData.getSize();
		AllData.setSize(size);
		if(AllData.getHostPtr()==NULL)
			AllData.hostAlloc();
		else
			printf("WARNING : host_alloc when host-ptr is non-null");
		AllData.setSize(temp_size);
	}
	void prepare_host(size_t alloc_size)
	{
		if(size==0)
		{
			printf("trying to device-alloc a stager with size=0");
			CRITICAL(ERR_STAGEMEM);
		}
		size_t temp_size=AllData.getSize();
		AllData.setSize(alloc_size);
		if(AllData.getHostPtr()==NULL)
			AllData.hostAlloc();
		else
			printf("WARNING : host_alloc when host-ptr is non-null");
		AllData.setSize(temp_size);
	}
	void prepare_device()
	{
		if(size==0)
		{
			printf("trying to host-alloc a stager with size=0");
			CRITICAL(ERR_STAGEMEM);
		}
		size_t temp_size=AllData.getSize();
		AllData.setSize(size);
		if(AllData.getDevicePtr()==NULL)
			AllData.deviceAlloc();
		else
			printf("WARNING : device_alloc when dev-ptr is non-null");
		AllData.setSize(temp_size);
	}
	void prepare_device(size_t alloc_size)
	{
		if(size==0)
		{
			printf("trying to device-alloc a stager with size=0");
			CRITICAL(ERR_STAGEMEM);
		}
		size_t temp_size=AllData.getSize();
		AllData.setSize(alloc_size);
		if(AllData.getDevicePtr()==NULL)
			AllData.deviceAlloc();
		else
			printf("WARNING : device_alloc when dev-ptr is non-null");
		AllData.setSize(temp_size);
	}
	void prepare()
	{
		 prepare_host();
		 prepare_device();
	}
	void prepare(size_t alloc_size)
	{
		 prepare_host(alloc_size);
		 prepare_device(alloc_size);
	}



	void stage(AccPtr<T> &input)
	{
		if(AllData.getSize()+input.getSize()>size)
		{
			printf("trying to stage more than stager can fit");
			printf(" (attempted to stage %lu addtionally to the allready staged %lu, and total host-allocated capacity is %lu ",input.getSize(),AllData.getSize(),size);
			exit( EXIT_FAILURE );
		}

		for(size_t i=0 ; i<input.getSize(); i++)
			AllData[AllData.getSize()+i] = input[i];

		// reset the staged object to this new position (TODO: disable for pinned mem)
		if(&input[0]!=NULL && input.willFreeHost())
			input.freeHostIfSet();
		input.setHostPtr(&AllData[AllData.getSize()]);
		input.setDevicePtr(&AllData(AllData.getSize()));

		AllData.setSize(AllData.getSize()+input.getSize());
	}

	void cp_to_device()
	{
		AllData.cpToDevice();
	}

	void cp_to_host()
	{
		AllData.cpToHost();
	}
};

#endif

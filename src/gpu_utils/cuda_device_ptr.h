#ifndef CUDA_DEVICE_PTR_H_
#define CUDA_DEVICE_PTR_H_

#include "src/gpu_utils/cuda_settings.h"
#include <cuda_runtime.h>
#include <signal.h>
#include <fstream>
#include <vector>

#ifdef DEBUG_CUDA
#define HANDLE_ERROR2( err ) (HandleError2( err, __FILE__, __LINE__ ))
#else
#define HANDLE_ERROR2( err ) (err) //Do nothing
#endif
static void HandleError2( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "DEBUG_ERROR: %s in %s at line %d\n",
        		cudaGetErrorString( err ), file, line );
		raise(SIGSEGV);
    }
}

template <typename T>
class CudaDevicePtr
{
	T *ptr;
	size_t size;
	size_t allocSize;
	cudaStream_t stream;
	bool free;

public:
	CudaDevicePtr():
		ptr(0), size(0), allocSize(0), stream(0), free(false)
	{};

	CudaDevicePtr(cudaStream_t stream):
		ptr(0), size(0), allocSize(0), stream(stream), free(false)
	{};

	CudaDevicePtr(size_t size):
		ptr(0), size(0), allocSize(0), stream(0), free(false)
	{
		resize(size);
	};

	CudaDevicePtr(size_t size, cudaStream_t stream):
		ptr(0), size(0), allocSize(0), stream(stream), free(false)
	{
		resize(size);
	};

	//Copy constructor
	CudaDevicePtr( const CudaDevicePtr<T>& other ):
		ptr(0), size(0), allocSize(0), stream(other.stream), free(false)
	{
		resize(other.size);
	};

	inline
	operator T*() { return ptr; }

	inline
	T* getPtr() { return ptr; }

	inline
	size_t getSize() { return size; }

	inline
	void resize(size_t newSize)
	{
		if (newSize > allocSize)
		{
			if (free) HANDLE_ERROR2(cudaFree(ptr));

			allocSize = newSize;
			HANDLE_ERROR2(cudaMalloc( (void**) &ptr, allocSize * sizeof(T)));
			free = true;
		}
		size = newSize;
	};

	inline
	void setStream(cudaStream_t s) { stream = s; }

	inline
	cudaStream_t getStream() { return stream; }

	inline
	void init(int value)
	{
		HANDLE_ERROR2(cudaMemsetAsync(ptr, value, size * sizeof(T), stream));
	};
	inline
	void set(T *h_ptr)
	{
		HANDLE_ERROR2(cudaMemcpyAsync( ptr, h_ptr, size * sizeof(T), cudaMemcpyHostToDevice, stream));
	};

	inline
	void set(std::vector<T> &h_ptr)
	{
		resize(h_ptr.size());
		set(&h_ptr[0]);
	};

	inline
	void get(T *h_ptr)
	{
		HANDLE_ERROR2(cudaMemcpyAsync( h_ptr, ptr, size * sizeof(T), cudaMemcpyDeviceToHost, stream));
	};

	inline
	void get(std::vector<T> &h_ptr)
	{
		resize(h_ptr.size());
		get(&h_ptr[0]);
	};

	inline
	T getAt(size_t idx)
	{
		T value;
		HANDLE_ERROR2(cudaMemcpyAsync( &value, ptr, sizeof(T), cudaMemcpyDeviceToHost, stream));
		return value;
	};

	~CudaDevicePtr()
	{
		if (free)
		{
			free = false;
			HANDLE_ERROR2(cudaFree(ptr));
			ptr = 0;
		}
	};
};

#endif

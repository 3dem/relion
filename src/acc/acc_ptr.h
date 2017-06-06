#ifndef ACC_PTR_H_
#define ACC_PTR_H_

#ifdef CUDA
#include "src/acc/cuda/cuda_settings.h"
#include <cuda_runtime.h>
#include "src/acc/cuda/shortcuts.cuh"
#include "src/acc/cuda/custom_allocator.cuh"
#include "src/acc/cuda/cuda_mem_utils.h"
#endif

#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstring>

#include "src/macros.h"
#include "src/error.h"
#include "src/parallel.h"

// Flags for accelerated code paths available
#define ACC_CPU 1
#define ACC_CUDA 2


#define ACC_PTR_DEBUG_FATAL( err ) (HandleAccPtrDebugFatal( err, __FILE__, __LINE__ ))
static void HandleAccPtrDebugFatal( const char *err, const char *file, int line )
{
    	fprintf(stderr, "DEBUG ERROR: %s in %s:%d\n", err, file, line );
		fflush(stdout);
		raise(SIGSEGV);
}

template <typename T, int AccT = ACC_CPU>
class AccPtr
{
	CudaCustomAllocator *allocator;
	CudaCustomAllocator::Alloc *alloc;
	cudaStream_t stream;

	size_t size; //Size used when copying data from and to device
	T *hPtr, *dPtr; //Host and device pointers
	bool doFreeHost, doFreeDevice; //True if host or device needs to be freed

public:

	/*======================================================
				CONSTRUCTORS WITH ALLOCATORS
	======================================================*/

	inline
	AccPtr(CudaCustomAllocator *allocator):
		size(0), hPtr(NULL), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(cudaStream_t stream, CudaCustomAllocator *allocator):
		size(0), hPtr(NULL), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(allocator), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(size_t size, CudaCustomAllocator *allocator):
		size(size), hPtr(NULL), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), hPtr(NULL), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(allocator), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(T * h_start, size_t size, CudaCustomAllocator *allocator):
		size(size), hPtr(h_start), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), hPtr(h_start), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, CudaCustomAllocator *allocator):
		size(size), hPtr(h_start), dPtr(d_start), doFreeHost(false),
		doFreeDevice(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), hPtr(h_start), dPtr(d_start), doFreeHost(false),
		doFreeDevice(false), allocator(allocator), alloc(NULL), stream(stream)
	{}

	/*======================================================
	           CONSTRUCTORS WITHOUT ALLOCATORS
	======================================================*/

	inline
	AccPtr():
		size(0), hPtr(NULL), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(cudaStream_t stream):
		size(0), hPtr(NULL), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(NULL), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(size_t size):
		size(size), hPtr(NULL), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(size_t size, cudaStream_t stream):
		size(size), hPtr(NULL), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(NULL), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(T * h_start, size_t size):
		size(size), hPtr(h_start), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(NULL), alloc(NULL), stream(0)
	{}

	inline
	AccPtr(T * h_start, size_t size, cudaStream_t stream):
		size(size), hPtr(h_start), dPtr(NULL), doFreeHost(false),
		doFreeDevice(false), allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size):
		size(size), hPtr(h_start), dPtr(d_start), doFreeHost(false),
		doFreeDevice(false), allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream):
		size(size), hPtr(h_start), dPtr(d_start), doFreeHost(false),
		doFreeDevice(false), allocator(NULL), alloc(NULL), stream(stream)
	{}

	/*======================================================
	       CONSTRUCTORS WITH OTHER POINTERS
	======================================================*/

	inline
	AccPtr(const AccPtr &ptr):
		size(ptr.size), hPtr(ptr.hPtr), dPtr(ptr.dPtr), doFreeHost(false),
		doFreeDevice(false), allocator(ptr.allocator), alloc(NULL), stream(ptr.stream)
	{}

	inline
	AccPtr(const AccPtr<T> &ptr, size_t start_idx, size_t size):
		size(size), hPtr(&ptr.h_ptr[start_idx]), dPtr(&ptr.dPtr[start_idx]), doFreeHost(false),
		doFreeDevice(false), allocator(ptr.allocator), alloc(NULL), stream(ptr.stream)
	{}


	/*======================================================
	                     METHOD BODY
	======================================================*/

	void markReadyEvent()
	{
		if (AccT == ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (alloc == NULL)
				ACC_PTR_DEBUG_FATAL("markReadyEvent called on null allocation.\n");
#endif
			alloc->markReadyEvent(stream);
		}
	}

	/**
	 * Allocate memory on device
	 */
	inline
	void deviceAlloc()
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if(size==0)
				ACC_PTR_DEBUG_FATAL("deviceAlloc called with size == 0");
			if (d_do_free)
				ACC_PTR_DEBUG_FATAL("Device double allocation.\n");
#endif
			doFreeDevice = true;

			alloc = allocator->alloc(size * sizeof(T));
			dPtr = (T*) alloc->getPtr();

//			DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &dPtr, size * sizeof(T)));
		}
	}

	/**
	 * Allocate memory on device with given size
	 */
	inline
	void deviceAlloc(size_t newSize)
	{
		size = newSize;
		deviceAlloc();
	}

	/**
	 * Allocate memory on host
	 */
	inline
	void hostAlloc()
	{
#ifdef DEBUG_CUDA
		if(size==0)
			ACC_PTR_DEBUG_FATAL("deviceAlloc called with size == 0");
		if (doFreeHost)
			ACC_PTR_DEBUG_FATAL("Host double allocation.\n");
#endif
		doFreeHost = true;
		hPtr = new T[size];
	}

	/**
	 * Allocate memory on host with given size
	 */
	inline
	void hostAlloc(size_t newSize)
	{
		size = newSize;
		hostAlloc();
	}

	inline
	void allAlloc()
	{
		deviceAlloc();
		hostAlloc();
	}

	inline
	void allAlloc(size_t newSize)
	{
		size = newSize;
		deviceAlloc();
		hostAlloc();
	}

	inline
	void accAlloc()
	{
		if (AccT == ACC_CUDA)
			deviceAlloc();
		else
			hostAlloc();
	}

	inline
	void accAlloc(size_t newSize)
	{
		size = newSize;
		if (AccT == ACC_CUDA)
			deviceAlloc(newSize);
		else
			hostAlloc(newSize);
	}

	void resizeHost(size_t newSize)
	{
#ifdef DEBUG_CUDA
		if (size==0)
			ACC_PTR_DEBUG_FATAL("Resizing from size zero (permitted).\n");
#endif
	    T* newArr = new T[newSize];
	    memcpy( newArr, hPtr, newSize * sizeof(T) );

	    size = newSize;
#ifdef DEBUG_CUDA
		if (dPtr!=NULL)
			ACC_PTR_DEBUG_FATAL("Resizing host with present device allocation.\n");
#endif
	    freeHost();
	    setHstPtr(newArr);
	    doFreeHost=true;
	}

	/**
	 * Initiate device memory with provided value
	 */
	inline
	void deviceInit(int value)
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (dPtr == NULL)
				ACC_PTR_DEBUG_FATAL("Memset requested before allocation in deviceInit().\n");
#endif
			cudaMemInit<T>( dPtr, value, size, stream);
		}
	}

	/**
	 * Initiate host memory with provided value
	 */
	inline
	void hostInit(int value)
	{
		memset(hPtr, value, size * sizeof(T));
	}

	/**
	 * Initiate accelerator memory with provided value
	 */
	inline
	void accInit(int value)
	{
		if (AccT == ACC_CUDA)
			deviceInit(value);
		else
			hostInit(value);
	}

	/**
	 * Copy a number (size) of bytes to device stored in the host pointer
	 */
	inline
	void cpToDevice()
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (dPtr == NULL)
				ACC_PTR_DEBUG_FATAL("cpToDevice() called before allocation.\n");
			if (hPtr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cpToDevice().\n");
#endif
			CudaShortcuts::cpyHostToDevice<T>(hPtr, dPtr, size, stream);
		}
	}

	/**
	 * Copy a number (size) of bytes to device stored in the provided host pointer
	 */
	inline
	void cpToDevice(T * hostPtr)
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (hostPtr == NULL)
				ACC_PTR_DEBUG_FATAL("Null-pointer given in cpToDevice(hostPtr).\n");
#endif
			hPtr = hostPtr;
			cpToDevice();
		}
	}

	/**
	 * alloc and copy
	 */
	inline
	void putOnDevice()
	{
		deviceAlloc();
		cpToDevice();
	}

	/**
	 * alloc size and copy
	 */
	inline
	void putOnDevice(size_t newSize)
	{
		size=newSize;
		deviceAlloc();
		cpToDevice();
	}


	/**
	 * Copy a number (size) of bytes from device to the host pointer
	 */
	inline
	void cpToHost()
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (dPtr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_host() called before device allocation.\n");
			if (h_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host().\n");
#endif
			cudaCpyDeviceToHost<T>(dPtr, hPtr, size, stream);
		}
	}

	/**
	 * Copy a number (thisSize) of bytes from device to the host pointer
	 */
	inline
	void cpToHost(size_t thisSize)
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (dPtr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_host(thisSize) called before device allocation.\n");
			if (h_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host(thisSize).\n");
#endif
			cudaCpyDeviceToHost<T>(dPtr, hPtr, thisSize, stream);
		}
	}

	/**
	 * Copy a number (thisSize) of bytes from device to a specific host pointer
	 */
	inline
	void cpToHost(T* hstPtr, size_t thisSize)
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (dPtr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_host(hstPtr, thisSize) called before device allocation.\n");
			if (hstPtr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host(hstPtr, thisSize).\n");
#endif
			cudaCpyDeviceToHost<T>(dPtr, hstPtr, thisSize, stream);
		}
	}

	/**
	 * Copy a number (size) of bytes from device to the host pointer
	 */
	inline
	void cpToHostOnStream(cudaStream_t s)
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (dPtr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_host_on_stream(s) called before device allocation.\n");
			if (h_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host_on_stream(s).\n");
#endif
			cudaCpyDeviceToHost<T>(dPtr, hPtr, size, s);
		}
	}

	/**
	 * Copy a number (size) of bytes from device pointer to the provided new device pointer
	 */
	inline
	void cpOnDevice(T * dstDevPtr)
	{
#ifdef DEBUG_CUDA
		if (dstDevPtr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL-pointer given in cp_on_device(dstDevPtr).\n");
#endif
		CudaShortcuts::cpyDeviceToDevice(dPtr, dstDevPtr, size, stream);
	}

	/**
	 * Copy a number (size) of bytes from host pointer to the provided new host pointer
	 */
	inline
	void cpOnHost(T * dstDevPtr)
	{
#ifdef DEBUG_CUDA
		if (dstDevPtr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL-pointer given in cp_on_host(dstDevPtr).\n");
#endif
		memcpy ( dstDevPtr, hPtr, size );
	}

	inline
	void cpOnAcc(T * dstDevPtr)
	{
		if (AccT == ACC_CUDA)
			cpOnDevice(dstDevPtr);
		else
			cpOnHost(dstDevPtr);
	}

	inline
	void cpOnAcc(AccPtr<T> &devPtr)
	{
		if (AccT == ACC_CUDA)
			cpOnDevice(devPtr.dPtr);
		else
			cpOnHost(devPtr.hPtr);
	}

	/**
	 * Host data quick access
	 */
	inline
	const T& operator[](size_t idx) const
	{
#ifdef DEBUG_CUDA
		if (hPtr == NULL)
			ACC_PTR_DEBUG_FATAL("const operator[] called with NULL host pointer.\n");
#endif
		return hPtr[idx];
	};

	/**
	 * Host data quick access
	 */
	inline
	T& operator[](size_t idx)
	{
#ifdef DEBUG_CUDA
		if (hPtr == NULL)
			ACC_PTR_DEBUG_FATAL("operator[] called with NULL host pointer.\n");
#endif
		return hPtr[idx];
	};

	/**
	 * Acc pointer quick access
	 */
	inline
	T* operator()()
	{
		if (AccT == ACC_CUDA)
		{
#ifdef DEBUG_CUDA
			if (dPtr == NULL)
				ACC_PTR_DEBUG_FATAL("operator~ called with NULL acc pointer.\n");
#endif
			return dPtr;
		}
		else
		{
#ifdef DEBUG_CUDA
			if (hPtr == NULL)
				ACC_PTR_DEBUG_FATAL("operator~ called with NULL acc pointer.\n");
#endif
			return hPtr;
		}
	};

	inline
	void streamSync()
	{
		if (AccT != ACC_CPU)
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
	}

	inline
	T getAccValueAt(size_t idx)
	{
		if (AccT == ACC_CPU)
			return hPtr[idx];
		else
		{
			T value;
			cudaCpyDeviceToHost<T>(&dPtr[idx], &value, 1, stream);
			streamSync();
			return value;
		}
	}


	void dumpDeviceToFile(std::string fileName)
	{
		T *tmp = new T[size];
		cudaCpyDeviceToHost<T>(dPtr, tmp, size, stream);

		std::ofstream f;
		f.open(fileName.c_str());
		streamSync();
		for (unsigned i = 0; i < size; i ++)
			f << tmp[i] << std::endl;
		f.close();
		delete [] tmp;
	}

	void dumpHostToFile(std::string fileName)
	{
		std::ofstream f;
		f.open(fileName.c_str());
		for (unsigned i = 0; i < size; i ++)
			f << hPtr[i] << std::endl;
		f.close();
	}

	/**
	 * Delete device data
	 */
	inline
	void freeDevice()
	{
		if (AccT != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (dPtr == NULL)
				ACC_PTR_DEBUG_FATAL("Free device memory was called on NULL pointer in free_device().\n");
#endif
			doFreeDevice = false;

			if (alloc->getReadyEvent() == 0)
				alloc->markReadyEvent(stream);
			alloc->doFreeWhenReady();
			alloc = NULL;

//			DEBUG_HANDLE_ERROR(cudaFree(dPtr));

			dPtr = NULL;
		}
	}

	/**
	 * Delete host data
	 */
	inline
	void freeHost()
	{
#ifdef DEBUG_CUDA
		if (hPtr == NULL)
			ACC_PTR_DEBUG_FATAL("free_host() called on NULL pointer.\n");
#endif
		doFreeHost = false;
		delete [] hPtr;
		hPtr = NULL;
	}

	inline
	void freeHostIfSet()
	{
		if (doFreeHost)
			freeHost();
	}

	inline
	void freeDeviceIfSet()
	{
		if (doFreeDevice)
			freeDevice();
	}

	/**
	 * Delete both device and host data
	 */
	inline
	void free()
	{
		freeDevice();
		freeHost();
	}

	inline
	void freeIfSet()
	{
		freeHostIfSet();
		freeDeviceIfSet();
	}

	inline
	~AccPtr()
	{
		freeIfSet();
	}


	/*======================================================
	                   GETTERS AND SETTERS
	======================================================*/


	bool willFreeHost()
	{
		return doFreeHost;
	}

	bool willFreeDevice()
	{
		return doFreeDevice;
	}

	void setStream(cudaStream_t s)
	{
		stream = s;
	}

	cudaStream_t getStream()
	{
		return stream;
	}

	void setSize(size_t s)
	{
		size = s;
	}
	size_t getSize()
	{
		return size;
	}

	T *getHostPtr()
	{
		return hPtr;
	}

	void setAllocator(CudaCustomAllocator *a)
	{
		freeDeviceIfSet();
		allocator = a;
	};

	CudaCustomAllocator *getAllocator()
	{
		return allocator;
	}

	void setDevicePtr(T *ptr)
	{
#ifdef DEBUG_CUDA
			if (d_do_free)
				ACC_PTR_DEBUG_FATAL("Device pointer set without freeing the old one.\n");
#endif
		dPtr = ptr;
	}

	T *getDevicePtr()
	{
		return dPtr;
	}

	void setDevicePtr(const AccPtr<T> &ptr)
	{
#ifdef DEBUG_CUDA
		if (ptr.dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("Device pointer is not set.\n");
#endif
		setDevicePtr(ptr.dPtr);
	};

	void setHostPtr(T *ptr)
	{
#ifdef DEBUG_CUDA
		if (doFreeHost)
			ACC_PTR_DEBUG_FATAL("Host pointer set without freeing the old one.\n");
#endif
		hPtr = ptr;
	};

	void setHostPtr(const AccPtr<T> &ptr)
	{
#ifdef DEBUG_CUDA
		if (ptr.hPtr == NULL)
			ACC_PTR_DEBUG_FATAL("Host pointer is not set.\n");
#endif
		setHstPtr(ptr.h_ptr);
	};

};

#endif

#ifndef ACC_PTR_H_
#define ACC_PTR_H_

//#ifdef __NVCC__
#ifdef CUDA
#include "src/acc/cuda/cuda_settings.h"
#include <cuda_runtime.h>
#include "src/acc/cuda/custom_allocator.cuh"
#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/cuda/shortcuts.cuh"
#endif
//#endif

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

// Flags for accelerated code paths available - defined in CMAKE
//#define ACC_CPU 1
//#define ACC_CUDA 2

//#ifndef __NVCC__
//typedef void CudaCustomAllocator;
//typedef void cudaStream_t;
//#else
#include "src/acc/cuda/cuda_mem_utils.h"
//#endif  


#define ACC_PTR_DEBUG_FATAL( err ) (HandleAccPtrDebugFatal( err, __FILE__, __LINE__ ))
static void HandleAccPtrDebugFatal( const char *err, const char *file, int line )
{
    	fprintf(stderr, "DEBUG ERROR: %s in %s:%d\n", err, file, line );
		fflush(stdout);
		raise(SIGSEGV);
}

//template <typename T, int AccT, typename derivedT>
//class AccPtrBase : crtp<derivedT>; // undefined
//{};

//////////////////////////////////////////////////////////////////////////////
// --------------------------- "Generic" Base Implementation ----------------------
//////////////////////////////////////////////////////////////////////////////

template <typename T, int AccT, typename derivedT>
class AccPtrBase 
{
public:
	size_t size; //Size used when copying data from and to device
	T *hPtr, *dPtr; //Host and device pointers
	bool doFreeHost, doFreeDevice; //True if host or device needs to be freed

//	friend derivedT;
	
	/*======================================================
	           CONSTRUCTORS WITHOUT ALLOCATORS
	======================================================*/
	
	inline
	AccPtrBase(size_t sz, T *hptr, T *dptr, bool dfH, bool dfD):
		size(sz), hPtr(hptr), dPtr(dptr), doFreeHost(dfH), doFreeDevice(dfD) 
	{}

	/*======================================================
	                     METHOD BODY
	======================================================*/

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
		// TODO - consider making this std::vector
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
		derivedT& derived = static_cast<derivedT&>(*this);
		derived.deviceAlloc();
		hostAlloc();
	}

	inline
	void allAlloc(size_t newSize)
	{
		size = newSize;
		derivedT& derived = static_cast<derivedT&>(*this);
		derived.deviceAlloc();
		hostAlloc();
	}

	inline
	void accAlloc()
	{
		if (ACC_CUDA==AccT)
		{
			derivedT& derived = static_cast<derivedT&>(*this);
			derived.deviceAlloc();
		}
		else
			hostAlloc();
	}

	inline
	void accAlloc(size_t newSize)
	{
		size = newSize;
		if (ACC_CUDA==AccT)
		{
			derivedT& derived = static_cast<derivedT&>(*this);
			derived.deviceAlloc(newSize);
		}
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
	    setHostPtr(newArr);
	    doFreeHost=true;
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
		if (ACC_CUDA==AccT)
		{
			derivedT& derived = static_cast<derivedT&>(*this);
			derived.deviceInit(value);
		}
		else
			hostInit(value);
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

	/**
	 * Copy a number (size) of bytes from device pointer to the provided new device pointer
	 */
	inline
	void cpOnAcc(T * dstDevPtr)
	{
		if (ACC_CUDA==AccT)
		{
			derivedT& derived = static_cast<derivedT&>(*this);
			derived.cpOnDevice(dstDevPtr);
		}
		else
			cpOnHost(dstDevPtr);
	}

	inline
	void cpOnAcc(AccPtrBase<T, AccT, derivedT> &devPtr)
	{
		if (ACC_CUDA==AccT)
		{
			derivedT& derived = static_cast<derivedT&>(*this);
			derived.cpOnDevice(devPtr.dPtr);
		}
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
	}

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
	}

	/**
	 * Device data quick access
	 */
	inline
	T& operator()(size_t idx) { return dPtr[idx]; };


	/**
	 * Device data quick access
	 */
	inline
	const T& operator()(size_t idx) const { return dPtr[idx]; };
	
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
	}
	
	inline
	T getAccValueAt(size_t idx)
	{
		if (AccT == ACC_CPU)
			return hPtr[idx];
		else
		{
			derivedT& derived = static_cast<derivedT&>(*this);
			return derived.getValueOnDevice(idx);
		}
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

	/**
	 * Delete both device and host data
	 */
	inline
	void free()
	{
		if (AccT == ACC_CUDA)
		{
			derivedT& derived = static_cast<derivedT&>(*this);
			derived.freeDevice();
		}
		freeHost();
	}

	inline
	void freeIfSet()
	{
		freeHostIfSet();
		if (AccT == ACC_CUDA)
		{
			derivedT& derived = static_cast<derivedT&>(*this);
			derived.freeDeviceIfSet();
		}
	}

	inline
	~AccPtrBase()
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

	void setHostPtr(T *ptr)
	{
#ifdef DEBUG_CUDA
		if (doFreeHost)
			ACC_PTR_DEBUG_FATAL("Host pointer set without freeing the old one.\n");
#endif
		hPtr = ptr;
	}

	void setHostPtr(const AccPtrBase<T, AccT, derivedT> &ptr)
	{
#ifdef DEBUG_CUDA
		if (ptr.hPtr == NULL)
			ACC_PTR_DEBUG_FATAL("Host pointer is not set.\n");
#endif
		setHostPtr(ptr.hPtr);
	}

};

//////////////////////////////////////////////////////////////////////////////
// --------------------------- "Generic" Implementation ----------------------
//////////////////////////////////////////////////////////////////////////////

template <typename T, int AccT>
class AccPtr : public AccPtrBase<T, AccT, AccPtr<T,AccT> >
{
//public:
//	~AccPtr()
//	{
	//	AccPtrBase<T, AccT, AccPtr<T,AccT> >::~AccPtrBase();
//	}
};

//////////////////////////////////////////////////////////////////////////////
// --------------------------- ACC_CPU Specialization -----------------------
//////////////////////////////////////////////////////////////////////////////
template <typename T>
class AccPtr<T, ACC_CPU> : public AccPtrBase<T, ACC_CPU, AccPtr<T, ACC_CPU> >
{
public:
	
// TODO - we don't need these - we need versions of the GPU
// constructors that redirect to the to allocator/stream-free CPU
// versions to allow use of the "GPU" code transparently
	inline
	AccPtr():
		AccPtrBase<T,ACC_CPU,AccPtr<T,ACC_CPU> >(0, NULL, NULL, false, false)
	{}
	
	inline
	AccPtr(size_t size):
		AccPtrBase<T,ACC_CPU,AccPtr<T,ACC_CPU> >(size, NULL, NULL, false, false)
	{}

	inline
	AccPtr(T * h_start, size_t size):
		AccPtrBase<T,ACC_CPU,AccPtr<T,ACC_CPU> >(size, h_start, NULL, false, false)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size):
		AccPtrBase<T,ACC_CPU,AccPtr<T,ACC_CPU> >(size, h_start, d_start, false, false)
	{}

	/*======================================================
	       CONSTRUCTORS WITH OTHER POINTERS
	======================================================*/

	inline
	AccPtr(const AccPtr &ptr):
		AccPtrBase<T,ACC_CPU,AccPtr<T,ACC_CPU> >(ptr.size, ptr.hPtr, ptr.dPtr, 
			false, false)
	{}

	inline
	AccPtr(const AccPtr<T, ACC_CPU> &ptr, size_t start_idx, size_t size):
		AccPtrBase<T,ACC_CPU,AccPtr<T,ACC_CPU> >(size, &ptr.h_ptr[start_idx], 
			&ptr.dPtr[start_idx], false, false)
	{}

	/*======================================================
	                     METHOD BODY
	======================================================*/
	void markReadyEvent()
	{}
	
	inline
	void deviceAlloc()
	{}
	
	inline
	void deviceAlloc(size_t newSize)
	{}
	
	inline
	void deviceInit(int value)
	{}
	
	inline
	void cpToDevice()
	{}
	
	inline
	void cpToDevice(T * hostPtr)
	{}
	
	inline
	void putOnDevice()
	{}
	
	inline
	void putOnDevice(size_t newSize)
	{}

	inline
	void cpToHost()
	{}
	
	inline
	void cpToHost(size_t thisSize)
	{}

	inline
	void cpToHost(T* hstPtr, size_t thisSize)
	{}
	
	inline
	void cpToHostOnStream(cudaStream_t s)
	{}
	
	inline
	void cpOnDevice(T * dstDevPtr)
	{}

// TODO - is this appropriate on the host?
	inline
	T* operator~() {
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: \"kernel cast\" on null pointer.\n");
#endif
		return AccPtrBase<T,ACC_CPU,AccPtr<T,ACC_CPU> >::dPtr;
	};
	
	inline
	void streamSync()
	{}
	
	void dumpDeviceToFile(std::string fileName)
	{}
	
	inline
	void freeDevice()
	{}
	
	inline
	void freeDeviceIfSet()
	{}
	
	/*======================================================
	                   GETTERS AND SETTERS
	======================================================*/
	void setStream(cudaStream_t s)
	{}

	cudaStream_t getStream()
	{
		return (cudaStream_t)0;
	}

	void setAllocator(CudaCustomAllocator *a)
	{}

	CudaCustomAllocator *getAllocator()
	{
		return (CudaCustomAllocator*)0;
	}

	void setDevicePtr(T *ptr)
	{}

	T *getDevicePtr()
	{
		return (T *)0;
	}

	void setDevicePtr(const AccPtrBase<T, ACC_CUDA, AccPtr> &ptr)
	{}
	
	inline
	T getValueOnDevice(size_t idx)
	{
		return (T)0;
	}
};

//////////////////////////////////////////////////////////////////////////////
// --------------------------- ACC_CUDA Specialization ----------------------
//////////////////////////////////////////////////////////////////////////////
template <typename T>
class AccPtr<T, ACC_CUDA> : public AccPtrBase<T, ACC_CUDA, AccPtr<T, ACC_CUDA> >
{
public:
        CudaCustomAllocator *allocator;
        CudaCustomAllocator::Alloc *alloc;
        cudaStream_t stream;
		
public:
	
	/*======================================================
				CONSTRUCTORS WITH ALLOCATORS
	======================================================*/

	inline
	AccPtr(CudaCustomAllocator *allocator):
		AccPtrBase<T, ACC_CUDA, AccPtr<T, ACC_CUDA> >(0, NULL, NULL, false, false),
		allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(cudaStream_t stream, CudaCustomAllocator *allocator):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(0, NULL, NULL, false, false),
		allocator(allocator), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(size_t size, CudaCustomAllocator *allocator):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, NULL, NULL, false, false),
		allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, NULL, NULL, false, false),
		allocator(allocator), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(T * h_start, size_t size, CudaCustomAllocator *allocator):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, h_start, NULL, false, false),
		allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, h_start, NULL, false, false),
		allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, CudaCustomAllocator *allocator):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, h_start, d_start, false, false),
		allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, h_start, d_start, false, false),		
		allocator(allocator), alloc(NULL), stream(stream)
	{}

	/*======================================================
	           CONSTRUCTORS WITHOUT ALLOCATORS
	======================================================*/

	inline
	AccPtr(cudaStream_t stream):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(0, NULL, NULL, false, false),
		allocator(NULL), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(size_t size, cudaStream_t stream):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, NULL, NULL, false, false),
		allocator(NULL), alloc(NULL), stream(stream)
	{}
	
	inline
	AccPtr(T * h_start, size_t size, cudaStream_t stream):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, h_start, NULL, false, false),
		allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, h_start, d_start, false, false),
		allocator(NULL), alloc(NULL), stream(stream)
	{}
	
	inline
	AccPtr():
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(0, NULL, NULL, false, false),
		allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{
	}
	
	inline
	AccPtr(size_t size):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, NULL, NULL, false, false),
		allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{
	}

	inline
	AccPtr(T * h_start, size_t size):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, h_start, NULL, false, false),
		allocator(NULL), alloc(NULL), stream(0)
	{
	}

	inline
	AccPtr(T * h_start, T * d_start, size_t size):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, h_start, d_start, false, false),
		allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{
	}
	
	/*======================================================
	       CONSTRUCTORS WITH OTHER POINTERS
	======================================================*/

	inline
	AccPtr(const AccPtr &ptr):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(ptr.size, ptr.hPtr, ptr.dPtr, 
			false, false),
			allocator(ptr.allocator), alloc(NULL), stream(ptr.stream)
	{
	}

	inline
	AccPtr(const AccPtr<T, ACC_CUDA> &ptr, size_t start_idx, size_t size):
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >(size, &ptr.hPtr[start_idx], 
			&ptr.dPtr[start_idx], false, false),
			allocator(ptr.allocator), alloc(NULL), stream(ptr.stream)
	{
	}
	
	/*======================================================
	                     METHOD BODY
	======================================================*/

	void markReadyEvent()
	{
#ifdef DEBUG_CUDA
		if (alloc == NULL)
				ACC_PTR_DEBUG_FATAL("markReadyEvent called on null allocation.\n");
	#endif
		alloc->markReadyEvent(stream);
	}

	/**
	 * Allocate memory on device
	 */
	inline
	void deviceAlloc()
	{
#ifdef DEBUG_CUDA
		if(size==0)
				ACC_PTR_DEBUG_FATAL("deviceAlloc called with size == 0");
		if (doFreeDevice)
				ACC_PTR_DEBUG_FATAL("Device double allocation.\n");
#endif
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::doFreeDevice = true;

		alloc = allocator->alloc(AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size * 
				sizeof(T));
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr = (T*) alloc->getPtr();
//                      DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &dPtr, size * sizeof(T)));
	}

	/**
	 * Allocate memory on device with given size
	 */
	inline
	void deviceAlloc(size_t newSize)
	{
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size = newSize;
		deviceAlloc();
	}

	/**
	 * Initiate device memory with provided value
	 */
	inline
	void deviceInit(int value)
	{
#ifdef DEBUG_CUDA
		if (dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("Memset requested before allocation in deviceInit().\n");
#endif
		cudaMemInit<T>( AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr, 
				value, AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size, stream);
	}
	
	/**
	 * Copy a number (size) of bytes to device stored in the host pointer
	 */
	inline
	void cpToDevice()
	{
#ifdef DEBUG_CUDA
		if (dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("cpToDevice() called before allocation.\n");
		if (hPtr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL host pointer in cpToDevice().\n");
#endif
		CudaShortcuts::cpyHostToDevice<T>(AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::hPtr, 
				AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr, 
				AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size, stream);
	}

	/**
	 * Copy a number (size) of bytes to device stored in the provided host pointer
	 */
	inline
	void cpToDevice(T * hostPtr)
	{
#ifdef DEBUG_CUDA
		if (hostPtr == NULL)
			ACC_PTR_DEBUG_FATAL("Null-pointer given in cpToDevice(hostPtr).\n");
#endif
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::hPtr = hostPtr;
		cpToDevice();
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
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size=newSize;
		deviceAlloc();
		cpToDevice();
	}

	/**
	 * Copy a number (size) of bytes from device to the host pointer
	 */
	inline
	void cpToHost()
	{
#ifdef DEBUG_CUDA
		if (dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("cp_to_host() called before device allocation.\n");
		if (h_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host().\n");
#endif
		cudaCpyDeviceToHost<T>(AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr, 
				AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::hPtr, 
				AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size, stream);
	}

	/**
	 * Copy a number (thisSize) of bytes from device to the host pointer
	 */
	inline
	void cpToHost(size_t thisSize)
	{
#ifdef DEBUG_CUDA
		if (dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("cp_to_host(thisSize) called before device allocation.\n");
		if (h_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host(thisSize).\n");
#endif
		cudaCpyDeviceToHost<T>(AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr, 
				AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::hPtr, thisSize, stream);
	}

	/**
	 * Copy a number (thisSize) of bytes from device to a specific host pointer
	 */
	inline
	void cpToHost(T* hstPtr, size_t thisSize)
	{
#ifdef DEBUG_CUDA
		if (dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("cp_to_host(hstPtr, thisSize) called before device allocation.\n");
		if (hstPtr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host(hstPtr, thisSize).\n");
	#endif
		cudaCpyDeviceToHost<T>(AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr, 
				AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::hstPtr, thisSize, stream);
	}

	/**
	 * Copy a number (size) of bytes from device to the host pointer
	 */
	inline
	void cpToHostOnStream(cudaStream_t s)
	{
#ifdef DEBUG_CUDA
		if (dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("cp_to_host_on_stream(s) called before device allocation.\n");
		if (h_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host_on_stream(s).\n");
#endif
		cudaCpyDeviceToHost<T>(	AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr, 
				AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::hPtr, 
				AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size, s);
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
		CudaShortcuts::cpyDeviceToDevice(AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr, 
				dstDevPtr, AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size, stream);
	}

	/**
	 * Acc pointer quick access
	 */
	inline
	T* operator~() {
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: \"kernel cast\" on null pointer.\n");
#endif
		return AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr;
	};
	
	inline
	void streamSync()
	{
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
	}

	inline
	T getDeviceAt(size_t idx)
	{
		T value;
		cudaCpyDeviceToHost<T>(&AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr[idx],
		&value, 1, stream);
		streamSync();
		return value;
	}
	void dumpDeviceToFile(std::string fileName)
	{
		T *tmp = new T[AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size];
		cudaCpyDeviceToHost<T>(AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr, 
				tmp, AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size, stream);

		std::ofstream f;
		f.open(fileName.c_str());
		streamSync();
		for (unsigned i = 0; i < AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::size; i ++)
			f << tmp[i] << std::endl;
		f.close();
		delete [] tmp;
	}

	/**
	 * Delete device data
	 */
	inline
	void freeDevice()
	{
#ifdef DEBUG_CUDA
		if (AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("Free device memory was called on NULL pointer in free_device().\n");
#endif
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::doFreeDevice = false;

		if (alloc->getReadyEvent() == 0)
			alloc->markReadyEvent(stream);
		alloc->doFreeWhenReady();
		alloc = NULL;

//			DEBUG_HANDLE_ERROR(cudaFree(dPtr));

		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr = NULL;
	}

	inline
	void freeDeviceIfSet()
	{
		if (AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::doFreeDevice)
			freeDevice();
	}
	
	/*======================================================
	                   GETTERS AND SETTERS
	======================================================*/
	void setStream(cudaStream_t s)
	{
		stream = s;
	}

	cudaStream_t getStream()
	{
		return stream;
	}

	void setAllocator(CudaCustomAllocator *a)
	{
		freeDeviceIfSet();
		allocator = a;
	}

	CudaCustomAllocator *getAllocator()
	{
		return allocator;
	}

	void setDevicePtr(T *ptr)
	{
#ifdef DEBUG_CUDA
			if (doFreeDevice)
				ACC_PTR_DEBUG_FATAL("Device pointer set without freeing the old one.\n");
#endif
		AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr = ptr;
	}

	T *getDevicePtr()
	{
		return AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr;
	}

	void setDevicePtr(const AccPtrBase<T, ACC_CUDA, AccPtr> &ptr)
	{
#ifdef DEBUG_CUDA
		if (ptr.dPtr == NULL)
			ACC_PTR_DEBUG_FATAL("Device pointer is not set.\n");
#endif
		setDevicePtr(ptr.dPtr);
	}
	
	inline
	T getValueOnDevice(size_t idx)
	{
		T value;
		cudaCpyDeviceToHost<T>(&AccPtrBase<T,ACC_CUDA,AccPtr<T,ACC_CUDA> >::dPtr[idx], &value, 1, stream);
		streamSync();
		return value;
	}
};

#endif

#ifndef ACC_PTR_H_
#define ACC_PTR_H_

#ifdef CUDA
#include "src/gpu_utils/cuda_settings.h"
#include <cuda_runtime.h>
#include "src/acc/cuda/shortcuts.cuh"
#include "src/acc/cuda/custom_allocator.cuh"
#include "src/gpu_utils/cuda_mem_utils.h"
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

#define ACC_CPU 0
#ifdef CUDA
#define ACC_CUDA 1
#else
#define ACC_CUDA 0
#endif

#define ACC_PTR_DEBUG_FATAL( err ) (HandleAccPtrDebugFatal( err, __FILE__, __LINE__ ))
static void HandleAccPtrDebugFatal( const char *err, const char *file, int line )
{
    	fprintf(stderr, "DEBUG ERROR: %s in %s:%d\n", err, file, line );
		fflush(stdout);
		raise(SIGSEGV);
}

template <typename T, int accType = ACC_CPU, bool CustomAlloc = true>
class AccPtr
{
	CudaCustomAllocator *allocator;
	CudaCustomAllocator::Alloc *alloc;
	cudaStream_t stream;

	size_t size; //Size used when copying data from and to device
	T *h_ptr, *d_ptr; //Host and device pointers
	bool h_do_free, d_do_free; //True if host or device needs to be freed

public:

	/*======================================================
				CONSTRUCTORS WITH ALLOCATORS
	======================================================*/

	inline
	AccPtr(CudaCustomAllocator *allocator):
		size(0), h_ptr(NULL), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(cudaStream_t stream, CudaCustomAllocator *allocator):
		size(0), h_ptr(NULL), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(size_t size, CudaCustomAllocator *allocator):
		size(size), h_ptr(NULL), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), h_ptr(NULL), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(T * h_start, size_t size, CudaCustomAllocator *allocator):
		size(size), h_ptr(h_start), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), h_ptr(h_start), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, CudaCustomAllocator *allocator):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(NULL), stream(stream)
	{}

	/*======================================================
	           CONSTRUCTORS WITHOUT ALLOCATORS
	======================================================*/

	inline
	AccPtr():
		size(0), h_ptr(NULL), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(cudaStream_t stream):
		size(0), h_ptr(NULL), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(NULL), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(size_t size):
		size(size), h_ptr(NULL), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(size_t size, cudaStream_t stream):
		size(size), h_ptr(NULL), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(NULL), alloc(NULL), stream(stream)
	{}

	inline
	AccPtr(T * h_start, size_t size):
		size(size), h_ptr(h_start), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(NULL), alloc(NULL), stream(0)
	{}

	inline
	AccPtr(T * h_start, size_t size, cudaStream_t stream):
		size(size), h_ptr(h_start), d_ptr(NULL), h_do_free(false),
		d_do_free(false), allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(NULL), alloc(NULL), stream(cudaStreamPerThread)
	{}

	inline
	AccPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(NULL), alloc(NULL), stream(stream)
	{}

	/*======================================================
	       CONSTRUCTORS WITH OTHER POINTERS
	======================================================*/

	inline
	AccPtr(const AccPtr &ptr):
		size(ptr.size), h_ptr(ptr.h_ptr), d_ptr(ptr.d_ptr), h_do_free(false),
		d_do_free(false), allocator(ptr.allocator), alloc(NULL), stream(ptr.stream)
	{}

	inline
	AccPtr(const AccPtr<T> &ptr, size_t start_idx, size_t size):
		size(size), h_ptr(&ptr.h_ptr[start_idx]), d_ptr(&ptr.d_ptr[start_idx]), h_do_free(false),
		d_do_free(false), allocator(ptr.allocator), alloc(NULL), stream(ptr.stream)
	{}


	/*======================================================
	                    OTHER STUFF
	======================================================*/

	CudaCustomAllocator *getAllocator() {return allocator; };
	cudaStream_t &getStream() {return stream; };

	void setStream(cudaStream_t s)
	{
		if (accType == ACC_CPU)
			stream = s;
	};

	void setSize(size_t s) { size = s; };
	size_t getSize() { return size; };

	T *get_h_ptr() { return h_ptr; };
	bool get_h_do_free() { return h_do_free; };

	T *get_d_ptr() { return d_ptr; };
	bool get_d_do_free() { return d_do_free; };

	void setAllocator(CudaCustomAllocator *a)
	{
		if (accType == ACC_CPU)
		{
			free_device_if_set();
			allocator = a;
		}
	};

	void markReadyEvent()
	{
		if (accType == ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (alloc == NULL)
				ACC_PTR_DEBUG_FATAL("markReadyEvent called on null allocation.\n");
#endif
			alloc->markReadyEvent(stream);
		}
	}

	void setDevPtr(T *ptr)
	{
		if (accType == ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (d_do_free)
				ACC_PTR_DEBUG_FATAL("Device pointer set without freeing the old one.\n");
#endif
			d_ptr = ptr;
		}
	};

	void setDevPtr(const AccPtr<T> &ptr)
	{
		if (accType == ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (ptr.d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("Device pointer is not set.\n");
#endif
			setHstPtr(ptr.d_ptr);
		}
	};

	void setHstPtr(T *ptr)
	{
#ifdef DEBUG_CUDA
		if (h_do_free)
			ACC_PTR_DEBUG_FATAL("Host pointer set without freeing the old one.\n");
#endif
		h_ptr = ptr;
	};

	void setHstPtr(const AccPtr<T> &ptr)
	{
#ifdef DEBUG_CUDA
		if (ptr.h_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("Host pointer is not set.\n");
#endif
		setHstPtr(ptr.h_ptr);
	};

	/**
	 * Allocate memory on device
	 */
	inline
	void device_alloc()
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if(size==0)
				ACC_PTR_DEBUG_FATAL("device_alloc called with size == 0");
			if (d_do_free)
				ACC_PTR_DEBUG_FATAL("Device double allocation.\n");
#endif
			d_do_free = true;
			if (CustomAlloc)
			{
				alloc = allocator->alloc(size * sizeof(T));
				d_ptr = (T*) alloc->getPtr();
			}
			else
				DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &d_ptr, size * sizeof(T)));
		}
	}

	/**
	 * Allocate memory on device with given size
	 */
	inline
	void device_alloc(size_t newSize)
	{
		size = newSize;
		device_alloc();
	}

	/**
	 * Allocate memory on host
	 */
	inline
	void host_alloc()
	{
#ifdef DEBUG_CUDA
		if(size==0)
			ACC_PTR_DEBUG_FATAL("device_alloc called with size == 0");
		if (h_do_free)
			ACC_PTR_DEBUG_FATAL("Host double allocation.\n");
#endif
		h_do_free = true;
		h_ptr = new T[size];
	}

	/**
	 * Allocate memory on host with given size
	 */
	inline
	void host_alloc(size_t newSize)
	{
		size = newSize;
		host_alloc();
	}

	inline
	void all_alloc()
	{
		device_alloc();
		host_alloc();
	}

	inline
	void all_alloc(size_t newSize)
	{
		size = newSize;
		device_alloc();
		host_alloc();
	}

	inline
	void acc_alloc()
	{
		if (accType == ACC_CUDA)
			device_alloc();
		else
			host_alloc();
	}

	inline
	void acc_alloc(size_t newSize)
	{
		size = newSize;
		if (accType == ACC_CUDA)
			device_alloc(newSize);
		else
			host_alloc(newSize);
	}

	void resize_host(size_t newSize)
	{
#ifdef DEBUG_CUDA
		if (size==0)
			ACC_PTR_DEBUG_FATAL("Resizing from size zero (permitted).\n");
#endif
	    T* newArr = new T[newSize];
	    memcpy( newArr, h_ptr, newSize * sizeof(T) );

	    size = newSize;
#ifdef DEBUG_CUDA
		if (d_ptr!=NULL)
			ACC_PTR_DEBUG_FATAL("Resizing host with present device allocation.\n");
#endif
	    free_host();
	    setHstPtr(newArr);
	    h_do_free=true;
	}

	/**
	 * Initiate device memory with provided value
	 */
	inline
	void device_init(int value)
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("Memset requested before allocation in device_init().\n");
#endif
			cudaMemInit<T>( d_ptr, value, size, stream);
		}
	}

	/**
	 * Initiate host memory with provided value
	 */
	inline
	void host_init(int value)
	{
		memset(h_ptr, value, size * sizeof(T));
	}

	/**
	 * Initiate accelerator memory with provided value
	 */
	inline
	void acc_init(int value)
	{
		if (accType == ACC_CUDA)
			device_init(value);
		else
			host_init(value);
	}

	/**
	 * Copy a number (size) of bytes to device stored in the host pointer
	 */
	inline
	void cp_to_device()
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_device() called before allocation.\n");
			if (h_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_device().\n");
#endif
			CudaShortcuts::cpyHostToDevice<T>(h_ptr, d_ptr, size, stream);
		}
	}

	/**
	 * Copy a number (size) of bytes to device stored in the provided host pointer
	 */
	inline
	void cp_to_device(T * hostPtr)
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (hostPtr == NULL)
				ACC_PTR_DEBUG_FATAL("Null-pointer given in cp_to_device(hostPtr).\n");
#endif
			h_ptr = hostPtr;
			cp_to_device();
		}
	}

	/**
	 * alloc and copy
	 */
	inline
	void put_on_device()
	{
		device_alloc();
		cp_to_device();
	}

	/**
	 * alloc size and copy
	 */
	inline
	void put_on_device(size_t newSize)
	{
		size=newSize;
		device_alloc();
		cp_to_device();
	}


	/**
	 * Copy a number (size) of bytes from device to the host pointer
	 */
	inline
	void cp_to_host()
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_host() called before device allocation.\n");
			if (h_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host().\n");
#endif
			cudaCpyDeviceToHost<T>(d_ptr, h_ptr, size, stream);
		}
	}

	/**
	 * Copy a number (thisSize) of bytes from device to the host pointer
	 */
	inline
	void cp_to_host(size_t thisSize)
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_host(thisSize) called before device allocation.\n");
			if (h_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host(thisSize).\n");
#endif
			cudaCpyDeviceToHost<T>(d_ptr, h_ptr, thisSize, stream);
		}
	}

	/**
	 * Copy a number (thisSize) of bytes from device to a specific host pointer
	 */
	inline
	void cp_to_host(T* hstPtr, size_t thisSize)
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_host(hstPtr, thisSize) called before device allocation.\n");
			if (hstPtr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host(hstPtr, thisSize).\n");
#endif
			cudaCpyDeviceToHost<T>(d_ptr, hstPtr, thisSize, stream);
		}
	}

	/**
	 * Copy a number (size) of bytes from device to the host pointer
	 */
	inline
	void cp_to_host_on_stream(cudaStream_t s)
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("cp_to_host_on_stream(s) called before device allocation.\n");
			if (h_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("NULL host pointer in cp_to_host_on_stream(s).\n");
#endif
			cudaCpyDeviceToHost<T>(d_ptr, h_ptr, size, s);
		}
	}

	/**
	 * Copy a number (size) of bytes from device pointer to the provided new device pointer
	 */
	inline
	void cp_on_device(T * dstDevPtr)
	{
#ifdef DEBUG_CUDA
		if (dstDevPtr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL-pointer given in cp_on_device(dstDevPtr).\n");
#endif
		CudaShortcuts::cpyDeviceToDevice(d_ptr, dstDevPtr, size, stream);
	}

	/**
	 * Copy a number (size) of bytes from host pointer to the provided new host pointer
	 */
	inline
	void cp_on_host(T * dstDevPtr)
	{
#ifdef DEBUG_CUDA
		if (dstDevPtr == NULL)
			ACC_PTR_DEBUG_FATAL("NULL-pointer given in cp_on_host(dstDevPtr).\n");
#endif
		memcpy ( dstDevPtr, h_ptr, size );
	}

	inline
	void cp_on_acc(T * dstDevPtr)
	{
		if (accType == ACC_CUDA)
			cp_on_device(dstDevPtr);
		else
			cp_on_host(dstDevPtr);
	}

	inline
	void cp_on_acc(AccPtr<T> &devPtr)
	{
		if (accType == ACC_CUDA)
			cp_on_device(devPtr.d_ptr);
		else
			cp_on_host(devPtr.h_ptr);
	}

	/**
	 * Host data quick access
	 */
	inline
	T& operator[](size_t idx)
	{
#ifdef DEBUG_CUDA
		if (h_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("operator[] called with NULL host pointer.\n");
#endif
		return h_ptr[idx];
	};


	/**
	 * Host data quick access
	 */
	inline
	const T& operator[](size_t idx) const
	{
#ifdef DEBUG_CUDA
		if (h_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("const operator[] called with NULL host pointer.\n");
#endif
		return h_ptr[idx];
	};

	/**
	 * Device data quick access
	 */
	inline
	T& operator()(size_t idx)
	{
#ifdef DEBUG_CUDA
		if (d_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("operator() called with NULL device pointer.\n");
#endif
		return d_ptr[idx];
	};


	/**
	 * Device data quick access
	 */
	inline
	const T& operator()(size_t idx) const
	{
#ifdef DEBUG_CUDA
		if (d_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("const operator() called with NULL device pointer.\n");
#endif
		return d_ptr[idx];
	};

	/**
	 * Device pointer quick access
	 */
	inline
	T* operator~()
	{
		if (accType == ACC_CUDA)
		{
#ifdef DEBUG_CUDA
			if (d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("operator~ called with NULL acc pointer.\n");
#endif
			return d_ptr;
		}
		else
		{
#ifdef DEBUG_CUDA
			if (h_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("operator~ called with NULL acc pointer.\n");
#endif
			return h_ptr;
		}
	};

	inline
	void streamSync()
	{
		if (accType != ACC_CPU)
			DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
	}

	inline
	T getAccValueAt(size_t idx)
	{
		if (accType == ACC_CPU)
			return h_ptr[idx];
		else
		{
			T value;
			cudaCpyDeviceToHost<T>(&d_ptr[idx], &value, 1, stream);
			streamSync();
			return value;
		}
	}


	void dump_device_to_file(std::string fileName)
	{
		T *tmp = new T[size];
		cudaCpyDeviceToHost<T>(d_ptr, tmp, size, stream);

		std::ofstream f;
		f.open(fileName.c_str());
		streamSync();
		for (unsigned i = 0; i < size; i ++)
			f << tmp[i] << std::endl;
		f.close();
		delete [] tmp;
	}

	void dump_host_to_file(std::string fileName)
	{
		std::ofstream f;
		f.open(fileName.c_str());
		for (unsigned i = 0; i < size; i ++)
			f << h_ptr[i] << std::endl;
		f.close();
	}

	/**
	 * Delete device data
	 */
	inline
	void free_device()
	{
		if (accType != ACC_CPU)
		{
#ifdef DEBUG_CUDA
			if (d_ptr == NULL)
				ACC_PTR_DEBUG_FATAL("Free device memory was called on NULL pointer in free_device().\n");
#endif
			d_do_free = false;
			if (CustomAlloc)
			{
				if (alloc->getReadyEvent() == 0)
					alloc->markReadyEvent(stream);
				alloc->doFreeWhenReady();

				alloc = NULL;
			}
			else
				DEBUG_HANDLE_ERROR(cudaFree(d_ptr));
			d_ptr = NULL;
		}
	}

	/**
	 * Delete host data
	 */
	inline
	void free_host()
	{
#ifdef DEBUG_CUDA
		if (h_ptr == NULL)
			ACC_PTR_DEBUG_FATAL("free_host() called on NULL pointer.\n");
#endif
		h_do_free = false;
		delete [] h_ptr;
		h_ptr = NULL;
	}

	inline
	void free_host_if_set()
	{
		if (h_do_free)
			free_host();
	}

	inline
	void free_device_if_set()
	{
		if (d_do_free)
			free_device();
	}

	/**
	 * Delete both device and host data
	 */
	inline
	void free()
	{
		free_device();
		free_host();
	}

	inline
	void free_if_set()
	{
		free_host_if_set();
		free_device_if_set();
	}

	inline
	~AccPtr()
	{
		free_if_set();
	}
};

#endif

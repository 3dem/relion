#ifndef CUDA_DEVICE_MEM_UTILS_H_
#define CUDA_DEVICE_MEM_UTILS_H_

#include "src/gpu_utils/cuda_settings.h"
#include <cuda_runtime.h>
#include <signal.h>
#include <fstream>
#include <vector>

#ifdef DEBUG_CUDA
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#else
#define HANDLE_ERROR( err ) (err) //Do nothing
#endif
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "DEBUG_ERROR: %s in %s at line %d\n",
        		cudaGetErrorString( err ), file, line );
		raise(SIGSEGV);
    }
}

/**
 * Print cuda device memory info
 */
static void cudaPrintMemInfo()
{
	size_t free;
	size_t total;
	HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
	float free_hr(free/(1024.*1024.));
	float total_hr(total/(1024.*1024.));
    printf( "free %.2fMiB, total %.2fMiB, used %.2fMiB\n",
    		free_hr, total_hr, total_hr - free_hr);
}

template< typename T>
static inline
void cudaCpyHostToDevice( T *h_ptr, T *d_ptr, size_t size)
{
	HANDLE_ERROR(cudaMemcpy( d_ptr, h_ptr, size * sizeof(T), cudaMemcpyHostToDevice));
};

template< typename T>
static inline
void cudaCpyHostToDevice( T *h_ptr, T *d_ptr, size_t size, cudaStream_t stream)
{
	HANDLE_ERROR(cudaMemcpyAsync( d_ptr, h_ptr, size * sizeof(T), cudaMemcpyHostToDevice, stream));
};

template< typename T>
static inline
void cudaCpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size)
{
	HANDLE_ERROR(cudaMemcpy( h_ptr, d_ptr, size * sizeof(T), cudaMemcpyDeviceToHost));
};

template< typename T>
static inline
void cudaCpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size, cudaStream_t stream)
{
	HANDLE_ERROR(cudaMemcpyAsync( h_ptr, d_ptr, size * sizeof(T), cudaMemcpyDeviceToHost, stream));
};

class CudaCustomAllocator
{
	typedef unsigned char BYTE;
	size_t totalSize;

public:
	struct Link
	{
		Link *prev, *next;
		BYTE *ptr;
		size_t size;
		bool free;
	};
private:
	Link *first;
public:

	CudaCustomAllocator(size_t size):
		totalSize(size)
	{
		first = new Link();

		first->prev = NULL;
		first->next = NULL;
		first->size = totalSize;
		first->free = true;

		HANDLE_ERROR(cudaMalloc( (void**) &(first->ptr), size));
	}

	inline
	void* alloc(size_t size)
	{
		Link *link = allocLink(size);
		return (void *) link->ptr;
	};

	inline
	Link* allocLink(size_t size)
	{
		Link *curL = first;

		//Look for the first suited link
		//If not the last and too small or not free go to next link
		while (curL != NULL && ( curL->size <= size || ! curL->free ) )
			curL = curL->next;

		if (curL == NULL)
		{
			printf("ERROR: CudaFixedMalloc out of memory.");
			raise(SIGSEGV);
		}

		if (curL->size == size)
		{
			curL->free = false;

//			printf("ALLOC: ");
//			printState();

			return curL;
		}
		else
		{
			//Setup new pointer
			Link *newL = new Link();
			newL->next = curL;
			newL->ptr = curL->ptr;
			newL->size = size;
			newL->free = false;

			//Modify old pointer
			curL->ptr = &(curL->ptr[size]);
			curL->size -= size;

			//Insert new link into chain
			if(curL->prev == NULL) //If the first link
				first = newL;
			else
				curL->prev->next = newL;
			newL->prev = curL->prev;
			newL->next = curL;
			curL->prev = newL;

//			printf("ALLOC: ");
//			printState();

			return newL;
		}
	};

	inline
	void free(void* ptr)
	{
		Link *curL = first;

		//Find the link holding the pointer
		while (curL != NULL && curL->ptr != ptr )
			curL = curL->next;

		if (curL == NULL)
		{
			printf("ERROR: CudaFixedMalloc could not find provided pointer in list (possible double free).");
			raise(SIGSEGV);
		}

		freeLink(curL);
	};

	inline
	void freeLink(Link* curL)
	{
		curL->free = true;

		//Previous neighbor is free, concatenate
		if ( curL->prev != NULL && curL->prev->free)
		{
			//Resize and set pointer
			curL->size += curL->prev->size;
			curL->ptr = curL->prev->ptr;

			//Fetch secondary neighbor
			Link *ppL = curL->prev->prev;

			//Remove primary neighbor
			if (ppL != NULL)
				ppL->next = curL;
			delete curL->prev;

			//Attach secondary neighbor
			curL->prev = ppL;
		}

		//Next neighbor is free, concatenate
		if ( curL->next != NULL && curL->next->free)
		{
			//Resize and set pointer
			curL->size += curL->next->size;

			//Fetch secondary neighbor
			Link *nnL = curL->next->next;

			//Remove primary neighbor
			if (nnL != NULL)
				nnL->prev = curL;
			delete curL->next;

			//Attach secondary neighbor
			curL->next = nnL;
		}

//		printf("FREE: ");
//		printState();
	};

	void printState()
	{
		Link *curL = first;

		while (curL != NULL)
		{
			if (curL->free)
				printf("[%lu,%lu] ", (unsigned long) curL->size, (unsigned long) curL->ptr);
			else
				printf("(%lu,%lu) ", (unsigned long) curL->size, (unsigned long) curL->ptr);

			curL = curL->next;
		}
		printf("\n");
	}

	~CudaCustomAllocator()
	{
		HANDLE_ERROR(cudaFree( first->ptr ));

		Link *cL = first, *nL;

		while (cL != NULL)
		{
			nL = cL->next;
			delete cL;
			cL = nL;
		}
	}
};


template <typename T, bool CustomAlloc=false>
class CudaGlobalPtr
{
	CudaCustomAllocator *allocator;
	CudaCustomAllocator::Link *link;
	cudaStream_t stream;
public:
	size_t size; //Size used when copying data from and to device
	T *h_ptr, *d_ptr; //Host and device pointers
	bool h_do_free, d_do_free; //True if host or device needs to be freed

	inline
	CudaGlobalPtr():
		size(0), h_ptr(0), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(0), link(0), stream(0)
	{};

	inline
	CudaGlobalPtr(size_t size):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
		d_do_free(false), allocator(0), link(0), stream(0)
	{};

	inline
	CudaGlobalPtr(size_t size, cudaStream_t stream):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
		d_do_free(false), allocator(0), link(0), stream(stream)
	{};

	inline
	CudaGlobalPtr(size_t size, CudaCustomAllocator *allocator):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
		d_do_free(false), allocator(allocator), link(0), stream(0)
	{};

	inline
	CudaGlobalPtr(size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
		d_do_free(false), allocator(allocator), link(0), stream(stream)
	{};

	inline
	CudaGlobalPtr(T * h_start, size_t size):
		size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(0), link(0), stream(0)
	{};

	inline
	CudaGlobalPtr(T * h_start, T * d_start, size_t size):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(0), link(0), stream(0)
	{};

	/**
	 * Allocate memory on device
	 */
	inline
	void device_alloc()
	{
#ifdef DEBUG_CUDA
		if (d_do_free)
			printf("DEBUG_WARNING: Device double allocation.\n");
#endif
		d_do_free = true;
		if (CustomAlloc)
		{
			link = allocator->allocLink(size * sizeof(T));
			d_ptr = (T*) link->ptr;
		}
		else
			HANDLE_ERROR(cudaMalloc( (void**) &d_ptr, size * sizeof(T)));
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
		if (h_do_free)
			printf("DEBUG_WARNING: Host double allocation.\n");
#endif
		h_do_free = true;
		h_ptr = new T[size];
	}

	/**
	 * Initiate device memory with provided value
	 */
	inline
	void device_init(int value)
	{
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: Memset requested before allocation in device_init().\n");
#endif
		HANDLE_ERROR(cudaMemset( d_ptr, value, size * sizeof(T)));
	}

	/**
	 * Copy a number (size) of bytes to device stored in the host pointer
	 */
	inline
	void cp_to_device()
	{
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: cp_to_device() called before allocation.\n");
		if (h_ptr == 0)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_device().\n");
#endif
		cudaCpyHostToDevice<T>(h_ptr, d_ptr, size);
	}

	/**
	 * Copy a number (size) of bytes to device stored in the provided host pointer
	 */
	inline
	void cp_to_device(T * hostPtr)
	{
#ifdef DEBUG_CUDA
		if (h_ptr != 0)
			printf("DEBUG_WARNING: Host pointer already set in call to cp_to_device(hostPtr).\n");
#endif
		h_ptr = hostPtr;
		cp_to_device();
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
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: cp_to_host() called before allocation.\n");
		if (h_ptr == 0)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_host().\n");
#endif
		cudaCpyDeviceToHost<T>(d_ptr, h_ptr, size);
	}

	/**
	 * Host data quick access
	 */
	inline
	T& operator[](size_t idx) { return h_ptr[idx]; };


	/**
	 * Host data quick access
	 */
	inline
	const T& operator[](size_t idx) const { return h_ptr[idx]; };

	/**
	 * Device pointer quick access
	 */
	inline
	T* operator~() {
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: \"kernel cast\" on null pointer.\n");
#endif
		return d_ptr;
	};

	/**
	 * Delete device data
	 */
	inline
	void free_device()
	{
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: Free device memory was called on NULL pointer in free_device().\n");
#endif
		d_do_free = false;
		if (CustomAlloc)
			allocator->freeLink(link);
		else
			HANDLE_ERROR(cudaFree(d_ptr));
		d_ptr = 0;
	}

	/**
	 * Delete host data
	 */
	inline
	void free_host()
	{
#ifdef DEBUG_CUDA
		if (h_ptr == 0)
		{
			printf("DEBUG_ERROR: free_host() called on NULL pointer.\n");
	        exit( EXIT_FAILURE );
		}
#endif
		h_do_free = false;
		delete [] h_ptr;
		h_ptr = 0;
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
	~CudaGlobalPtr()
	{
		if (d_do_free) free_device();
		if (h_do_free) free_host();
	}
};

#endif

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
void cudaCpyHostToDevice( T *h_ptr, T *d_ptr, size_t size, cudaStream_t &stream)
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
void cudaCpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size, cudaStream_t &stream)
{
	HANDLE_ERROR(cudaMemcpyAsync( h_ptr, d_ptr, size * sizeof(T), cudaMemcpyDeviceToHost, stream));
};

class OutOfMemoryHandler
{
public:
	//Returns true if any memory was freed
	virtual void handleOutOfMemory() {};
};

class CudaCustomAllocator
{
	typedef unsigned char BYTE;
	size_t totalSize;
	OutOfMemoryHandler *outOfMemoryHandler;

public:
	class Alloc
	{
		friend class CudaCustomAllocator;

	private:
		Alloc *prev, *next;
		BYTE *ptr;
		size_t size;
		bool free;

	public:
		inline
		BYTE *getPtr() { return ptr; }

		inline
		size_t getSize() { return size; }

		inline
		bool isFree() { return free; }
	};
private:
	Alloc *first;
public:

	CudaCustomAllocator(size_t size):
		totalSize(size),outOfMemoryHandler(NULL), first(0)
	{
#ifndef CUDA_NO_CUSTOM_ALLOCATION
		first = new Alloc();

		first->prev = NULL;
		first->next = NULL;
		first->size = totalSize;
		first->free = true;

		HANDLE_ERROR(cudaMalloc( (void**) &(first->ptr), size));
#endif
	}


	inline
	Alloc* alloc(size_t size)
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		Alloc *nAlloc = new Alloc();
		HANDLE_ERROR(cudaMalloc( (void**) &(nAlloc->ptr), size));
		return nAlloc;
#else
		Alloc *curL = first;

		//Look for the first suited link
		//If not the last and too small or not free go to next link
		while (curL != NULL && ( curL->size <= size || ! curL->free ) )
			curL = curL->next;

		//If out of memory
		if (curL == NULL)
		{
			if (outOfMemoryHandler != NULL) //Is there a handler
			{
				outOfMemoryHandler->handleOutOfMemory(); //Try to handle

				//Is there space now?
				curL = first;
				while (curL != NULL && ( curL->size <= size || ! curL->free ) )
					curL = curL->next;
			}

			//Did we manage to recover?
			if (curL == NULL)
			{
				//Nope... OK I'm out
				printf("ERROR: CudaCustomAllocator out of memory\n [requestedSpace: %lu B]\n [largestContinuousFreeSpace: %lu B]\n [totalFreeSpace: %lu B]\n",
						(unsigned long) size, (unsigned long) getLargestContinuousFreeSpace(), (unsigned long) getTotalFreeSpace());
				fflush(stdout);
				raise(SIGSEGV);
			}
		}

		if (curL->size == size)
		{
			curL->free = false;

			printf("ALLOC: ");
			printState();

			return curL;
		}
		else
		{
			//Setup new pointer
			Alloc *newL = new Alloc();
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
#endif
	};


	inline
	void free(Alloc* curL)
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		HANDLE_ERROR(cudaFree( curL->ptr ));
#else
		curL->free = true;

		//Previous neighbor is free, concatenate
		if ( curL->prev != NULL && curL->prev->free)
		{
			//Resize and set pointer
			curL->size += curL->prev->size;
			curL->ptr = curL->prev->ptr;

			//Fetch secondary neighbor
			Alloc *ppL = curL->prev->prev;

			//Remove primary neighbor
			if (ppL == NULL) //If the previous is first in chain
				first = curL;
			else
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
			Alloc *nnL = curL->next->next;

			//Remove primary neighbor
			if (nnL != NULL)
				nnL->prev = curL;
			delete curL->next;

			//Attach secondary neighbor
			curL->next = nnL;
		}

//		printf("FREE: ");
//		printState();
#endif
	};

	size_t getTotalFreeSpace()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		size_t free, total;
		HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
		return free;
#else
		size_t total = 0;
		Alloc *cL = first;

		while (cL != NULL) //Get last
		{
			if (cL->free)
				total += cL->size;
			cL = cL->next;
		}
		return total;
#endif
	}

	size_t getTotalUsedSpace()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		size_t free, total;
		HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
		return total - free;
#else
		size_t total = 0;
		Alloc *cL = first;

		while (cL != NULL) //Get last
		{
			if (!cL->free)
				total += cL->size;
			cL = cL->next;
		}
		return total;
#endif
	}

	size_t getLargestContinuousFreeSpace()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		return getTotalFreeSpace();
#else
		size_t largest = 0;
		Alloc *cL = first;

		while (cL != NULL) //Get last
		{
			if (cL->free && cL->size > largest)
				largest = cL->size;
			cL = cL->next;
		}
		return largest;
#endif
	}

	void printState()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		printf("Custom allocation is disabled.\n");
#else
		Alloc *curL = first;
		size_t total = 0;

		while (curL != NULL)
		{
			total += curL->size;
			if (curL->free)
				printf("[%luB] ", (unsigned long) curL->size);
			else
				printf("(%luB) ", (unsigned long) curL->size);

			curL = curL->next;
		}
		printf("= (%luB)\n", (unsigned long) total);
#endif
		fflush(stdout);
	}

	void setOutOfMemoryHandler(OutOfMemoryHandler *handler)
	{ outOfMemoryHandler = handler; }

	~CudaCustomAllocator()
	{
		HANDLE_ERROR(cudaFree( first->ptr ));

		Alloc *cL = first, *nL;

		while (cL != NULL)
		{
			nL = cL->next;
			delete cL;
			cL = nL;
		}
	}
};

template <typename T, bool CustomAlloc=true>
class CudaGlobalPtr
{
	CudaCustomAllocator *allocator;
	CudaCustomAllocator::Alloc *alloc;
	cudaStream_t stream;
public:
	size_t size; //Size used when copying data from and to device
	T *h_ptr, *d_ptr; //Host and device pointers
	bool h_do_free, d_do_free; //True if host or device needs to be freed

	/*======================================================
				CONSTRUCTORS WITH ALLOCATORS
	======================================================*/

	inline
	CudaGlobalPtr(CudaCustomAllocator *allocator):
		size(0), h_ptr(0), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(size_t size, CudaCustomAllocator *allocator):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
		d_do_free(false), allocator(allocator), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
		d_do_free(false), allocator(allocator), alloc(0), stream(stream)
	{};

	inline
	CudaGlobalPtr(T * h_start, size_t size, CudaCustomAllocator *allocator):
		size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(T * h_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(T * h_start, T * d_start, size_t size, CudaCustomAllocator *allocator):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(0), stream(stream)
	{};

	/*======================================================
	           CONSTRUCTORS WITHOUT ALLOCATORS
	======================================================*/

	inline
	CudaGlobalPtr():
		size(0), h_ptr(0), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(0), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(size_t size):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
		d_do_free(false), allocator(0), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(size_t size, cudaStream_t stream):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
		d_do_free(false), allocator(0), alloc(0), stream(stream)
	{};

	inline
	CudaGlobalPtr(T * h_start, size_t size):
		size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(0), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(T * h_start, size_t size, cudaStream_t stream):
		size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(0), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(T * h_start, T * d_start, size_t size):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(0), alloc(0), stream(0)
	{};

	inline
	CudaGlobalPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream):
		size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
		d_do_free(false), allocator(0), alloc(0), stream(stream)
	{};

	/*======================================================
	                    OTHER STUFF
	======================================================*/

	CudaCustomAllocator *getAllocator() {return allocator; };
	cudaStream_t &getStream() {return stream; };
	void setStream(cudaStream_t s) { stream = s; };

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
//			HANDLE_ERROR(cudaDeviceSynchronize());
			alloc = allocator->alloc(size * sizeof(T));
//			HANDLE_ERROR(cudaDeviceSynchronize());
			d_ptr = (T*) alloc->getPtr();
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
		HANDLE_ERROR(cudaMemsetAsync( d_ptr, value, size * sizeof(T), stream));
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
		cudaCpyHostToDevice<T>(h_ptr, d_ptr, size, stream);
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
			printf("DEBUG_WARNING: cp_to_host() called before device allocation.\n");
		if (h_ptr == 0)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_host().\n");
#endif
		cudaCpyDeviceToHost<T>(d_ptr, h_ptr, size, stream);
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
		{
//			HANDLE_ERROR(cudaDeviceSynchronize());
			allocator->free(alloc);
//			HANDLE_ERROR(cudaDeviceSynchronize());
		}
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

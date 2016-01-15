#ifndef CUDA_DEVICE_MEM_UTILS_H_
#define CUDA_DEVICE_MEM_UTILS_H_

#include "src/gpu_utils/cuda_settings.h"
#include <cuda_runtime.h>
#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <pthread.h>

#ifdef DEBUG_CUDA
#define DEBUG_HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#else
#define DEBUG_HANDLE_ERROR( err ) (err) //Do nothing
#endif

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
static void HandleError( cudaError_t err, const char *file, int line )
{

    if (err != cudaSuccess)
    {
#ifdef DEBUG_CUDA
        printf( "DEBUG_ERROR: %s in %s at line %d\n",
        		cudaGetErrorString( err ), file, line );
        fflush(stdout);
		raise(SIGSEGV);
#else
        printf( "ERROR: %s in %s at line %d\n",
        		cudaGetErrorString( err ), file, line );
        fflush(stdout);
		raise(SIGSEGV);
#endif
    }

#ifdef DEBUG_CUDA
	cudaError_t peek = cudaPeekAtLastError();
    if (peek != cudaSuccess)
    {
        printf( "DEBUG_ERROR: %s in %s at line %d\n",
        		cudaGetErrorString( peek ), file, line );
        fflush(stdout);
		raise(SIGSEGV);
    }
#endif

}

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
void cudaCpyDeviceToDevice( T *src, T *des, size_t size, cudaStream_t &stream)
{
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( des, src, size * sizeof(T), cudaMemcpyDeviceToDevice, stream));
};

template< typename T>
static inline
void cudaMemInit( T *ptr, int value, size_t size)
{
	DEBUG_HANDLE_ERROR(cudaMemset( ptr, value, size * sizeof(T)));
};

template< typename T>
static inline
void cudaMemInit( T *ptr, int value, size_t size, cudaStream_t &stream)
{
	DEBUG_HANDLE_ERROR(cudaMemsetAsync( ptr, value, size * sizeof(T), stream));
};






class CudaCustomAllocator
{
	typedef unsigned char BYTE;

public:
	class Alloc
	{
		friend class CudaCustomAllocator;

	private:
		Alloc *prev, *next;
		BYTE *ptr;
		size_t size;
		bool free;
		cudaEvent_t readyEvent; //Event record used for auto free
		bool freeWhenReady;

		Alloc():
			prev(NULL), next(NULL),
			ptr(NULL),
			size(0),
			free(0),
			readyEvent(0),
			freeWhenReady(false)
		{}

		~Alloc()
		{
			if (readyEvent != 0)
				DEBUG_HANDLE_ERROR(cudaEventDestroy(readyEvent));
		}

	public:
		inline
		BYTE *getPtr() { return ptr; }

		inline
		size_t getSize() { return size; }

		inline
		bool isFree() { return free; }

		inline
		cudaEvent_t getReadyEvent() { return readyEvent; }

		inline
		void markReadyEvent(cudaStream_t stream = 0)
		{
			//TODO add a debug warning if event already set
			DEBUG_HANDLE_ERROR(cudaEventCreate(&readyEvent));
			DEBUG_HANDLE_ERROR(cudaEventRecord(readyEvent, stream));
		}

		inline
		void doFreeWhenReady() { freeWhenReady = true; }
	};

private:

	Alloc *first;
	size_t totalSize;
	size_t alignmentSize;

	pthread_mutex_t mutex;


	//Look for the first suited space
	inline
	Alloc *_getFirstSuitedFree(size_t size)
	{
		Alloc *a = first;
		//If not the last and too small or not free go to next allocation region
		while (a != NULL && ( a->size <= size || ! a->free ) )
			a = a->next;

		return a;
	}

	//Free allocs with recorded ready events
	inline
	void _syncReadyEvents()
	{
		Alloc *a = first;

		while (a != NULL)
		{
			if (! a->free && a->freeWhenReady && a->readyEvent != 0)
				DEBUG_HANDLE_ERROR(cudaEventSynchronize(a->readyEvent));

			a = a->next;
		}
	}

	//Free allocs with recorded ready events
	inline
	void _freeReadyAllocs()
	{
		Alloc *a = first;

		while (a != NULL)
		{
			if (! a->free && a->freeWhenReady && a->readyEvent != 0)
			{
				cudaError_t e = cudaEventQuery(a->readyEvent);

				if (e == cudaSuccess)
					_free(a);
				else if (e != cudaErrorNotReady)
					HandleError( e, __FILE__, __LINE__ );
			}

			a = a->next;
		}
	}

	inline
	size_t _getTotalFreeSpace()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		size_t free, total;
		DEBUG_HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
		return free;
#else

		size_t total = 0;
		Alloc *a = first;

		while (a != NULL)
		{
			if (a->free)
				total += a->size;
			a = a->next;
		}

		return total;
#endif
	}

	inline
	size_t _getTotalUsedSpace()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		size_t free, total;
		DEBUG_HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
		return total - free;
#else

		size_t total = 0;
		Alloc *a = first;

		while (a != NULL)
		{
			if (!a->free)
				total += a->size;
			a = a->next;
		}

		return total;
#endif
	}

	size_t _getNumberOfAllocs()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		return 0;
#else

		size_t total = 0;
		Alloc *a = first;

		while (a != NULL)
		{
			if (!a->free)
				total ++;
			a = a->next;
		}

		return total;
#endif
	}

	inline
	size_t _getLargestContinuousFreeSpace()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		return getTotalFreeSpace();
#else

		size_t largest = 0;
		Alloc *a = first;

		while (a != NULL)
		{
			if (a->free && a->size > largest)
				largest = a->size;
			a = a->next;
		}

		return largest;
#endif
	}

	inline
	void _printState()
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		printf("Custom allocation is disabled.\n");
#else

		size_t total = 0;
		Alloc *a = first;

		while (a != NULL)
		{
			total += a->size;
			if (a->free)
				printf("[%luB] ", (unsigned long) a->size);
			else if (a->freeWhenReady)
				printf("<%luB> ", (unsigned long) a->size);
			else
				printf("(%luB) ", (unsigned long) a->size);

			a = a->next;
		}

		printf("= %luB\n", (unsigned long) total);
#endif
		fflush(stdout);
	}

	inline
	void _free(Alloc* a)
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		DEBUG_HANDLE_ERROR(cudaFree( a->ptr ));
		a->free = true;
#else

		a->free = true;

		//Previous neighbor is free, concatenate
		if ( a->prev != NULL && a->prev->free)
		{
			//Resize and set pointer
			a->size += a->prev->size;
			a->ptr = a->prev->ptr;

			//Fetch secondary neighbor
			Alloc *ppL = a->prev->prev;

			//Remove primary neighbor
			if (ppL == NULL) //If the previous is first in chain
				first = a;
			else
				ppL->next = a;

			delete a->prev;

			//Attach secondary neighbor
			a->prev = ppL;
		}

		//Next neighbor is free, concatenate
		if ( a->next != NULL && a->next->free)
		{
			//Resize and set pointer
			a->size += a->next->size;

			//Fetch secondary neighbor
			Alloc *nnL = a->next->next;

			//Remove primary neighbor
			if (nnL != NULL)
				nnL->prev = a;
			delete a->next;

			//Attach secondary neighbor
			a->next = nnL;
		}
#endif
	};

public:

	CudaCustomAllocator(size_t size, size_t alignmentSize):
		totalSize(size), alignmentSize(alignmentSize), first(0)
	{
#ifndef CUDA_NO_CUSTOM_ALLOCATION
		first = new Alloc();

		first->prev = NULL;
		first->next = NULL;
		first->size = totalSize;
		first->free = true;

		HANDLE_ERROR(cudaMalloc( (void**) &(first->ptr), size));

//		pthread_mutexattr_t mutex_attr;
//		pthread_mutexattr_init(&mutex_attr);
//		pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_ERRORCHECK);
//		int mutex_error = pthread_mutex_init(&mutex, &mutex_attr);
		int mutex_error = pthread_mutex_init(&mutex, NULL);

		if (mutex_error != 0)
		{
			printf("ERROR: Mutex could not be created for alloactor. CODE: %d.\n", mutex_error);
			fflush(stdout);
			raise(SIGSEGV);
		}
#endif
	}


	inline
	Alloc* alloc(size_t size)
	{
#ifdef CUDA_NO_CUSTOM_ALLOCATION
		Alloc *nAlloc = new Alloc();
		nAlloc->size = size;
		nAlloc->free = false;
		DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &(nAlloc->ptr), size));
		return nAlloc;
#else

		size = alignmentSize*ceilf( (float)size / (float)alignmentSize) ; //To prevent miss-aligned memory

		pthread_mutex_lock(&mutex);

		_freeReadyAllocs();
		Alloc *curAlloc = _getFirstSuitedFree(size);

		//If out of memory
		if (curAlloc == NULL)
		{
#ifdef DEBUG_CUDA
			int spaceDiff = _getTotalFreeSpace();
#endif
			_syncReadyEvents();
			_freeReadyAllocs();
#ifdef DEBUG_CUDA
			spaceDiff = ( (int) _getTotalFreeSpace() ) - spaceDiff;
			printf("DEBUG_INFO: Out of memory handled by waiting for unfinished tasks, which freed %d B.\n", spaceDiff);
#endif

			curAlloc = _getFirstSuitedFree(size); //Is there space now?

			//Did we manage to recover?
			if (curAlloc == NULL)
			{
				printf("ERROR: CudaCustomAllocator out of memory\n [requestedSpace:             %lu B]\n [largestContinuousFreeSpace: %lu B]\n [totalFreeSpace:             %lu B]\n",
						(unsigned long) size, (unsigned long) _getLargestContinuousFreeSpace(), (unsigned long) _getTotalFreeSpace());

				_printState();
				pthread_mutex_unlock(&mutex);

				fflush(stdout);
				raise(SIGSEGV);
			}
		}

		if (curAlloc->size == size)
		{
			curAlloc->free = false;

			pthread_mutex_unlock(&mutex);

			return curAlloc;
		}
		else
		{
			//Setup new pointer
			Alloc *newAlloc = new Alloc();
			newAlloc->next = curAlloc;
			newAlloc->ptr = curAlloc->ptr;
			newAlloc->size = size;
			newAlloc->free = false;

			//Modify old pointer
			curAlloc->ptr = &(curAlloc->ptr[size]);
			curAlloc->size -= size;

			//Insert new allocation region into chain
			if(curAlloc->prev == NULL) //If the first allocation region
				first = newAlloc;
			else
				curAlloc->prev->next = newAlloc;
			newAlloc->prev = curAlloc->prev;
			newAlloc->next = curAlloc;
			curAlloc->prev = newAlloc;

			pthread_mutex_unlock(&mutex);

			return newAlloc;
		}
#endif
	};

	~CudaCustomAllocator()
	{
#ifndef CUDA_NO_CUSTOM_ALLOCATION

		pthread_mutex_lock(&mutex);

		DEBUG_HANDLE_ERROR(cudaFree( first->ptr ));

		Alloc *a = first, *nL;

		while (a != NULL)
		{
			nL = a->next;
			delete a;
			a = nL;
		}

		pthread_mutex_unlock(&mutex);
		pthread_mutex_destroy(&mutex);
#endif
	}

	//Thread-safe wrapper functions

	inline
	void free(Alloc* a)
	{
		pthread_mutex_lock(&mutex);
		_free(a);
		pthread_mutex_unlock(&mutex);
	};

	inline
	void syncReadyEvents()
	{
		pthread_mutex_lock(&mutex);
		_syncReadyEvents();
		pthread_mutex_unlock(&mutex);
	}

	inline
	void freeReadyAllocs()
	{
		pthread_mutex_lock(&mutex);
		_freeReadyAllocs();
		pthread_mutex_unlock(&mutex);
	}

	size_t getTotalFreeSpace()
	{
		pthread_mutex_lock(&mutex);
		size_t size = _getTotalFreeSpace();
		pthread_mutex_unlock(&mutex);
		return size;
	}

	size_t getTotalUsedSpace()
	{
		pthread_mutex_lock(&mutex);
		size_t size = _getTotalUsedSpace();
		pthread_mutex_unlock(&mutex);
		return size;
	}

	size_t getNumberOfAllocs()
	{
		pthread_mutex_lock(&mutex);
		size_t size = _getNumberOfAllocs();
		pthread_mutex_unlock(&mutex);
		return size;
	}

	size_t getLargestContinuousFreeSpace()
	{
		pthread_mutex_lock(&mutex);
		size_t size = _getLargestContinuousFreeSpace();
		pthread_mutex_unlock(&mutex);
		return size;
	}

	void printState()
	{
		pthread_mutex_lock(&mutex);
		_printState();
		pthread_mutex_unlock(&mutex);
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
	CudaGlobalPtr(cudaStream_t stream, CudaCustomAllocator *allocator):
		size(0), h_ptr(0), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(allocator), alloc(0), stream(stream)
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
	CudaGlobalPtr(cudaStream_t stream):
		size(0), h_ptr(0), d_ptr(0), h_do_free(false),
		d_do_free(false), allocator(0), alloc(0), stream(stream)
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
	       CONSTRUCTORS WITH OTHER GLOBAL POINTERS
	======================================================*/

	inline
	CudaGlobalPtr(const CudaGlobalPtr<T> &ptr):
		size(ptr.size), h_ptr(ptr.h_ptr), d_ptr(ptr.d_ptr), h_do_free(false),
		d_do_free(false), allocator(ptr.allocator), alloc(0), stream(ptr.stream)
	{};

	inline
	CudaGlobalPtr(const CudaGlobalPtr<T> &ptr, size_t start_idx, size_t size):
		size(size), h_ptr(&ptr.h_ptr[start_idx]), d_ptr(&ptr.d_ptr[start_idx]), h_do_free(false),
		d_do_free(false), allocator(ptr.allocator), alloc(0), stream(ptr.stream)
	{};



	/*======================================================
	                    OTHER STUFF
	======================================================*/

	CudaCustomAllocator *getAllocator() {return allocator; };
	cudaStream_t &getStream() {return stream; };
	void setStream(cudaStream_t s) { stream = s; };

	void setSize(size_t s) { size = s; };
	size_t getSize() { return size; };

	void markReadyEvent()
	{
#ifdef DEBUG_CUDA
		if (alloc == NULL)
			printf("DEBUG_WARNING: markReadyEvent called on null allocation.\n");
#endif
		alloc->markReadyEvent(stream);
	}

	void setDevPtr(T *ptr)
	{
#ifdef DEBUG_CUDA
		if (d_do_free)
			printf("DEBUG_WARNING: Device pointer set without freeing the old one.\n");
#endif
		d_ptr = ptr;
	};

	void setDevPtr(const CudaGlobalPtr<T> &ptr)
	{
#ifdef DEBUG_CUDA
		if (ptr.d_ptr == NULL)
			printf("DEBUG_WARNING: Device pointer is not set.\n");
#endif
		setHstPtr(ptr.d_ptr);
	};

	void setHstPtr(T *ptr)
	{
#ifdef DEBUG_CUDA
		if (h_do_free)
			printf("DEBUG_WARNING: Host pointer set without freeing the old one.\n");
#endif
		h_ptr = ptr;
	};

	void setHstPtr(const CudaGlobalPtr<T> &ptr)
	{
#ifdef DEBUG_CUDA
		if (ptr.h_ptr == NULL)
			printf("DEBUG_WARNING: Host pointer is not set.\n");
#endif
		setHstPtr(ptr.h_ptr);
	};

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
			alloc = allocator->alloc(size * sizeof(T));
			d_ptr = (T*) alloc->getPtr();
		}
		else
			DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &d_ptr, size * sizeof(T)));
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
	 * Allocate memory on host with given size
	 */
	inline
	void host_alloc(size_t newSize)
	{
		size = newSize;
		host_alloc();
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
		cudaMemInit<T>( d_ptr, value, size, stream);
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
		if (hostPtr == NULL)
			printf("DEBUG_WARNING: Null-pointer given in cp_to_device(hostPtr).\n");
#endif
		h_ptr = hostPtr;
		cp_to_device();
	}

	/**
	 * Copy a number (size) of bytes from device pointer to the provided new device pointer
	 */
	inline
	void cp_on_device(T * dstDevPtr)
	{
#ifdef DEBUG_CUDA
		if (dstDevPtr == NULL)
			printf("DEBUG_WARNING: Null-pointer given in cp_on_device(dstDevPtr).\n");
#endif
		cudaCpyDeviceToDevice(d_ptr, dstDevPtr, size, stream);
	}

	/**
	 * Copy a number (size) of bytes from device pointer to the provided new device pointer
	 */
	inline
	void cp_on_device(CudaGlobalPtr<T> devPtr)
	{
#ifdef DEBUG_CUDA
		if (devPtr.size == 0)
			printf("DEBUG_WARNING: Zero size on provided pointer in cp_on_device.\n");
#endif
		cp_on_device(devPtr.d_ptr);
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
		if (d_ptr == NULL)
			printf("DEBUG_WARNING: cp_to_host() called before device allocation.\n");
		if (h_ptr == NULL)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_host().\n");
#endif
		cudaCpyDeviceToHost<T>(d_ptr, h_ptr, size, stream);
	}

	/**
	 * Copy a number (thisSize) of bytes from device to the host pointer
	 */
	inline
	void cp_to_host(int thisSize)
	{
#ifdef DEBUG_CUDA
		if (d_ptr == NULL)
			printf("DEBUG_WARNING: cp_to_host(thisSize) called before device allocation.\n");
		if (h_ptr == NULL)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_host(thisSize).\n");
#endif
		cudaCpyDeviceToHost<T>(d_ptr, h_ptr, thisSize, stream);
	}

	/**
	 * Copy a number (thisSize) of bytes from device to a specific host pointer
	 */
	inline
	void cp_to_host(T* hstPtr, int thisSize)
	{
#ifdef DEBUG_CUDA
		if (d_ptr == NULL)
			printf("DEBUG_WARNING: cp_to_host(hstPtr, thisSize) called before device allocation.\n");
		if (hstPtr == NULL)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_host(hstPtr, thisSize).\n");
#endif
		cudaCpyDeviceToHost<T>(d_ptr, hstPtr, thisSize, stream);
	}

	/**
	 * Copy a number (size) of bytes from device to the host pointer
	 */
	inline
	void cp_to_host_on_stream(cudaStream_t s)
	{
#ifdef DEBUG_CUDA
		if (d_ptr == NULL)
			printf("DEBUG_WARNING: cp_to_host_on_stream(s) called before device allocation.\n");
		if (h_ptr == NULL)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_host_on_stream(s).\n");
#endif
		cudaCpyDeviceToHost<T>(d_ptr, h_ptr, size, s);
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
	 * Host data quick access
	 */
	inline
	T& operator()(size_t idx) { return d_ptr[idx]; };


	/**
	 * Host data quick access
	 */
	inline
	const T& operator()(size_t idx) const { return d_ptr[idx]; };

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

	inline
	void streamSync()
	{
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
	}

	inline
	T getDeviceAt(size_t idx)
	{
		T value;
		cudaCpyDeviceToHost<T>(&d_ptr[idx], &value, 1, stream);
		streamSync();
		return value;
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
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: Free device memory was called on NULL pointer in free_device().\n");
#endif
		d_do_free = false;

		if (CustomAlloc)
		{
			if (alloc->getReadyEvent() == 0)
				alloc->markReadyEvent();
			alloc->doFreeWhenReady();

			alloc = NULL;
		}
		else
			DEBUG_HANDLE_ERROR(cudaFree(d_ptr));
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
	~CudaGlobalPtr()
	{
		free_if_set();
	}
};

template <typename T>
class cudaStager
{
public:
	CudaGlobalPtr<T> AllData;
	unsigned long size; // size of allocated host-space (AllData.size dictates the amount of memory copied to/from the device)

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
		AllData.size=0;
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
		AllData.size=0;
	};

public:

	void prepare_host()
	{

		if(size==0)
		{
			printf("trying to host-alloc a stager with size=0");
			raise(SIGSEGV);
		}
		int temp_size=AllData.size;
		AllData.size=size;
		if(AllData.h_ptr==NULL)
			AllData.host_alloc();
		else
			printf("WARNING : host_alloc when host-ptr is non-null");
		AllData.size=temp_size;
	}
	void prepare_host(int alloc_size)
	{
		if(size==0)
		{
			printf("trying to device-alloc a stager with size=0");
			raise(SIGSEGV);
		}
		int temp_size=AllData.size;
		AllData.size=alloc_size;
		if(AllData.h_ptr==NULL)
			AllData.host_alloc();
		else
			printf("WARNING : host_alloc when host-ptr is non-null");
		AllData.size=temp_size;
	}
	void prepare_device()
	{
		if(size==0)
		{
			printf("trying to host-alloc a stager with size=0");
			raise(SIGSEGV);
		}
		int temp_size=AllData.size;
		AllData.size=size;
		if(AllData.d_ptr==NULL)
			AllData.device_alloc();
		else
			printf("WARNING : device_alloc when dev-ptr is non-null");
		AllData.size=temp_size;
	}
	void prepare_device(int alloc_size)
	{
		if(size==0)
		{
			printf("trying to device-alloc a stager with size=0");
			raise(SIGSEGV);
		}
		int temp_size=AllData.size;
		AllData.size=alloc_size;
		if(AllData.d_ptr==NULL)
			AllData.device_alloc();
		else
			printf("WARNING : device_alloc when dev-ptr is non-null");
		AllData.size=temp_size;
	}
	void prepare()
	{
		 prepare_host();
		 prepare_device();
	}
	void prepare(int alloc_size)
	{
		 prepare_host(alloc_size);
		 prepare_device(alloc_size);
	}



	void stage(CudaGlobalPtr<T> &input)
	{
		if(AllData.size+input.size>size)
		{
			printf("trying to stage more than stager can fit");
			printf(" (attempted to stage %lu addtionally to the allready staged %lu, and total host-allocated capacity is %lu ",input.size,AllData.size,size);
			exit( EXIT_FAILURE );
		}

		for(int i=0 ; i<input.size; i++)
			AllData.h_ptr[AllData.size+i] = input.h_ptr[i];

		// reset the staged object to this new position (TODO: disable for pinned mem)
		if(input.h_ptr!=NULL && input.h_do_free)
			input.free_host();
		input.h_ptr=&AllData.h_ptr[AllData.size];
		input.d_ptr=&AllData.d_ptr[AllData.size];

		AllData.size+=input.size;
	}



	void cp_to_device()
	{
		AllData.cp_to_device();
	}

	void cp_to_host()
	{
		AllData.cp_to_host();
	}
};

#endif

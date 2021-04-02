#ifndef CUDA_CUSTOM_ALLOCATOR_CUH_
#define CUDA_CUSTOM_ALLOCATOR_CUH_
// This is where custom allocator should be. Commented out for now, to avoid double declaration.

#ifdef _CUDA_ENABLED
#include "src/acc/cuda/cuda_settings.h"
#include <cuda_runtime.h>
#endif

#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "src/macros.h"
#include "src/error.h"
#include "src/parallel.h"

#ifdef CUSTOM_ALLOCATOR_MEMGUARD
#include <execinfo.h>
#include <cxxabi.h>
#endif

#ifdef DUMP_CUSTOM_ALLOCATOR_ACTIVITY
#define CUSTOM_ALLOCATOR_REGION_NAME( name ) (fprintf(stderr, "\n%s", name))
#else
#define CUSTOM_ALLOCATOR_REGION_NAME( name ) //Do nothing
#endif


class CudaCustomAllocator
{

	typedef unsigned char BYTE;

	const static unsigned GUARD_SIZE = 4;
	const static BYTE GUARD_VALUE = 145;
	const static int ALLOC_RETRY = 500;

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


#ifdef CUSTOM_ALLOCATOR_MEMGUARD
		BYTE *guardPtr;
		void *backtrace[20];
		size_t backtraceSize;
#endif

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
			prev = NULL;
			next = NULL;
			ptr = NULL;

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

	bool cache;

	omp_lock_t mutex;


	//Look for the first suited space
	Alloc *_getFirstSuitedFree(size_t size)
	{
		Alloc *a = first;
		//If not the last and too small or not free go to next allocation region
		while (a != NULL && ( a->size <= size || ! a->free ) )
			a = a->next;

		return a;
	}

	//Free allocs with recorded ready events
	bool _syncReadyEvents()
	{
		bool somethingReady(false);
		Alloc *a = first;

		while (a != NULL)
		{
			if (! a->free && a->freeWhenReady && a->readyEvent != 0)
			{
				DEBUG_HANDLE_ERROR(cudaEventSynchronize(a->readyEvent));
				somethingReady = true;
			}

			a = a->next;
		}

		return somethingReady;
	}

	//Free allocs with recorded ready events
	bool _freeReadyAllocs()
	{
		bool somethingFreed(false);
		Alloc *next = first;
		Alloc *curr;

		while (next != NULL)
		{
			curr = next;
			next = curr->next;

			if (! curr->free && curr->freeWhenReady && curr->readyEvent != 0)
			{
				cudaError_t e = cudaEventQuery(curr->readyEvent);

				if (e == cudaSuccess)
				{
					_free(curr);
					next = first; //List modified, restart
					somethingFreed = true;
				}
				else if (e != cudaErrorNotReady)
				{
					_printState();
					HandleError( e, __FILE__, __LINE__ );
				}
			}
		}
		return somethingFreed;
	}

	size_t _getTotalFreeSpace()
	{
		if (cache)
		{
			size_t total = 0;
			Alloc *a = first;

			while (a != NULL)
			{
				if (a->free)
					total += a->size;
				a = a->next;
			}

			return total;
		}
		else
		{
			size_t free, total;
			DEBUG_HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
			return free;
		}
	}

	size_t _getTotalUsedSpace()
	{
		size_t total = 0;
		Alloc *a = first;

		while (a != NULL)
		{
			if (!a->free)
				total += a->size;
			a = a->next;
		}

		return total;
	}

	size_t _getNumberOfAllocs()
	{

		size_t total = 0;
		Alloc *a = first;

		while (a != NULL)
		{
			if (!a->free)
				total ++;
			a = a->next;
		}

		return total;
	}

	size_t _getLargestContinuousFreeSpace()
	{
		if (cache)
		{
			size_t largest = 0;
			Alloc *a = first;

			while (a != NULL)
			{
				if (a->free && a->size > largest)
					largest = a->size;
				a = a->next;
			}

			return largest;
		}
		else
			return _getTotalFreeSpace();
	}

	void _printState()
	{
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
		fflush(stdout);
	}

	void _free(Alloc* a)
	{
//		printf("free: %u ", a->size);
//		_printState();


#ifdef CUSTOM_ALLOCATOR_MEMGUARD
		size_t guardCount = a->size - (a->guardPtr - a->ptr);
		BYTE *guards = new BYTE[guardCount];
		cudaStream_t stream = 0;
		CudaShortcuts::cpyDeviceToHost<BYTE>( a->guardPtr, guards, guardCount, stream);
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
		for (int i = 0; i < guardCount; i ++)
			if (guards[i] != GUARD_VALUE)
			{
				fprintf (stderr, "ERROR: CORRUPTED BYTE GUARDS DETECTED\n");

				char ** messages = backtrace_symbols(a->backtrace, a->backtraceSize);

				// skip first stack frame (points here)
				for (int i = 1; i < a->backtraceSize && messages != NULL; ++i)
				{
					char *mangled_name = 0, *offset_begin = 0, *offset_end = 0;

					// find parantheses and +address offset surrounding mangled name
					for (char *p = messages[i]; *p; ++p)
					{
						if (*p == '(')
						{
							mangled_name = p;
						}
						else if (*p == '+')
						{
							offset_begin = p;
						}
						else if (*p == ')')
						{
							offset_end = p;
							break;
						}
					}

					// if the line could be processed, attempt to demangle the symbol
					if (mangled_name && offset_begin && offset_end &&
						mangled_name < offset_begin)
					{
						*mangled_name++ = '\0';
						*offset_begin++ = '\0';
						*offset_end++ = '\0';

						int status;
						char * real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);

						// if demangling is successful, output the demangled function name
						if (status == 0)
						{
							std::cerr << "[bt]: (" << i << ") " << messages[i] << " : "
									  << real_name << "+" << offset_begin << offset_end
									  << std::endl;

						}
						// otherwise, output the mangled function name
						else
						{
							std::cerr << "[bt]: (" << i << ") " << messages[i] << " : "
									  << mangled_name << "+" << offset_begin << offset_end
									  << std::endl;
						}
//						free(real_name);
					}
					// otherwise, print the whole line
					else
					{
						std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
					}
				}
				std::cerr << std::endl;

//				free(messages);

				exit(EXIT_FAILURE);
			}
		delete[] guards;
#endif

		a->free = true;

		if (cache)
		{
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
		}
		else
		{
			DEBUG_HANDLE_ERROR(cudaFree( a->ptr ));
			a->ptr = NULL;

			if ( a->prev != NULL)
				a->prev->next = a->next;
			else
				first = a->next; //This is the first link

			if ( a->next != NULL)
				a->next->prev = a->prev;

			delete a;
		}
	};

	void _setup()
	{
		first = new Alloc();
		first->prev = NULL;
		first->next = NULL;
		first->size = totalSize;
		first->free = true;

		if (totalSize > 0)
		{
			HANDLE_ERROR(cudaMalloc( (void**) &(first->ptr), totalSize));
			cache = true;
		}
		else
			cache = false;
	}

	void _clear()
	{
		if (first->ptr != NULL)
			DEBUG_HANDLE_ERROR(cudaFree( first->ptr ));

		first->ptr = NULL;

		Alloc *a = first, *nL;

		while (a != NULL)
		{
			nL = a->next;
			delete a;
			a = nL;
		}
	}

public:

	CudaCustomAllocator(size_t size, size_t alignmentSize):
		totalSize(size), alignmentSize(alignmentSize), first(0), cache(true)
	{
		_setup();

		omp_init_lock(&mutex);
	}

	void resize(size_t size)
	{
		Lock ml(&mutex);
		_clear();
		totalSize = size;
		_setup();
	}


	Alloc* alloc(size_t requestedSize)
	{
		Lock ml(&mutex);

		_freeReadyAllocs();

//		printf("alloc: %u ", size);
//		_printState();

		size_t size = requestedSize;

#ifdef CUSTOM_ALLOCATOR_MEMGUARD
		//Ad byte-guards
		size += alignmentSize * GUARD_SIZE; //Ad an integer multiple of alignment size as byte guard size
#endif

#ifdef DUMP_CUSTOM_ALLOCATOR_ACTIVITY
		fprintf(stderr, " %.4f", 100.*(float)size/(float)totalSize);
#endif

		Alloc *newAlloc(NULL);

		if (cache)
		{
			size = alignmentSize*ceilf( (float)size / (float)alignmentSize) ; //To prevent miss-aligned memory

			Alloc *curAlloc = _getFirstSuitedFree(size);

			//If out of memory
			if (curAlloc == NULL)
			{
	#ifdef DEBUG_CUDA
				size_t spaceDiff = _getTotalFreeSpace();
	#endif
				//Try to recover before throwing error
				for (int i = 0; i <= ALLOC_RETRY; i ++)
				{
					if (_syncReadyEvents() && _freeReadyAllocs())
					{
						curAlloc = _getFirstSuitedFree(size); //Is there space now?
						if (curAlloc != NULL)
							break; //Success
					}
					else
						usleep(10000); // 10 ms, Order of magnitude of largest kernels
				}
	#ifdef DEBUG_CUDA
				spaceDiff =  _getTotalFreeSpace() - spaceDiff;
				printf("DEBUG_INFO: Out of memory handled by waiting for unfinished tasks, which freed %lu B.\n", spaceDiff);
	#endif

				//Did we manage to recover?
				if (curAlloc == NULL)
				{
					printf("ERROR: CudaCustomAllocator out of memory\n [requestedSpace:             %lu B]\n [largestContinuousFreeSpace: %lu B]\n [totalFreeSpace:             %lu B]\n",
							(unsigned long) size, (unsigned long) _getLargestContinuousFreeSpace(), (unsigned long) _getTotalFreeSpace());

					_printState();

					fflush(stdout);
					CRITICAL(ERRCUDACAOOM);
				}
			}

			if (curAlloc->size == size)
			{
				curAlloc->free = false;
				newAlloc = curAlloc;
			}
			else //Or curAlloc->size is smaller than size
			{
				//Setup new pointer
				newAlloc = new Alloc();
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
			}
		}
		else
		{
			newAlloc = new Alloc();
			newAlloc->size = size;
			newAlloc->free = false;
			DEBUG_HANDLE_ERROR(cudaMalloc( (void**) &(newAlloc->ptr), size));

			//Just add to start by replacing first
			newAlloc->next = first;
			first->prev = newAlloc;
			first = newAlloc;
		}

#ifdef CUSTOM_ALLOCATOR_MEMGUARD
		newAlloc->backtraceSize = backtrace(newAlloc->backtrace, 20);
		newAlloc->guardPtr = newAlloc->ptr + requestedSize;
		cudaStream_t stream = 0;
		CudaShortcuts::memInit<BYTE>( newAlloc->guardPtr, GUARD_VALUE, size - requestedSize, stream); //TODO switch to specialized stream
		DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
#endif

		return newAlloc;
	};

	~CudaCustomAllocator()
	{
		{
			Lock ml(&mutex);
			_clear();
		}
		omp_destroy_lock(&mutex);
	}

	//Thread-safe wrapper functions

	void free(Alloc* a)
	{
		Lock ml(&mutex);
		_free(a);
	}

	void syncReadyEvents()
	{
		Lock ml(&mutex);
		_syncReadyEvents();
	}

	void freeReadyAllocs()
	{
		Lock ml(&mutex);
		_freeReadyAllocs();
	}

	size_t getTotalFreeSpace()
	{
		Lock ml(&mutex);
		size_t size = _getTotalFreeSpace();
		return size;
	}

	size_t getTotalUsedSpace()
	{
		Lock ml(&mutex);
		size_t size = _getTotalUsedSpace();
		return size;
	}

	size_t getNumberOfAllocs()
	{
		Lock ml(&mutex);
		size_t size = _getNumberOfAllocs();
		return size;
	}

	size_t getLargestContinuousFreeSpace()
	{
		Lock ml(&mutex);
		size_t size = _getLargestContinuousFreeSpace();
		return size;
	}

	void printState()
	{
		Lock ml(&mutex);
		_printState();
	}
};
//

#endif

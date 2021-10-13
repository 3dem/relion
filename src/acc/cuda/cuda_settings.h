#ifndef CUDA_SETTINGS_H_
#define CUDA_SETTINGS_H_

#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "src/macros.h"
#include "src/error.h"

#include <curand.h>

// Required compute capability
#define CUDA_CC_MAJOR 3
#define CUDA_CC_MINOR 5

#define LAUNCH_CHECK
#define CUDA_BENCHMARK_OLD true

// Error handling ----------------------

#ifdef LAUNCH_CHECK
#define LAUNCH_HANDLE_ERROR( err ) (LaunchHandleError( err, __FILE__, __LINE__ ))
#define LAUNCH_PRIVATE_ERROR(func, status) {  \
                       (status) = (func); \
                       LAUNCH_HANDLE_ERROR(status); \
                   }
#else
#define LAUNCH_HANDLE_ERROR( err ) (err) //Do nothing
#define LAUNCH_PRIVATE_ERROR( err ) (err) //Do nothing
#endif

#ifdef DEBUG_CUDA
#define DEBUG_HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define DEBUG_PRIVATE_ERROR(func, status) {  \
                       (status) = (func); \
                       DEBUG_HANDLE_ERROR(status); \
                   }
#else
#define DEBUG_HANDLE_ERROR( err ) (err) //Do nothing
#define DEBUG_PRIVATE_ERROR( err ) (err) //Do nothing
#endif

#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#define PRIVATE_ERROR(func, status) {  \
                       (status) = (func); \
                       HANDLE_ERROR(status); \
                   }

static void HandleError( cudaError_t err, const char *file, int line )
{

    if (err != cudaSuccess)
    {
    	fprintf(stderr, "ERROR: %s in %s at line %d (error-code %d)\n",
						cudaGetErrorString( err ), file, line, err );
		fflush(stdout);
#ifdef DEBUG_CUDA
		raise(SIGSEGV);
#else
		CRITICAL(ERRGPUKERN);
#endif
    }
}

#ifdef LAUNCH_CHECK
static void LaunchHandleError( cudaError_t err, const char *file, int line )
{

    if (err != cudaSuccess)
    {
        printf( "KERNEL_ERROR: %s in %s at line %d (error-code %d)\n",
                        cudaGetErrorString( err ), file, line, err );
        fflush(stdout);
        CRITICAL(ERRGPUKERN);
    }
}
#endif


// GENERAL -----------------------------
#define MAX_RESOL_SHARED_MEM 		32
#define BLOCK_SIZE  				128
// -------------------------------------


// COARSE DIFF -------------------------
#define D2C_BLOCK_SIZE_2D 			512
#define D2C_EULERS_PER_BLOCK_2D 	4

#define D2C_BLOCK_SIZE_REF3D 		128
#define D2C_EULERS_PER_BLOCK_REF3D 	16

#define D2C_BLOCK_SIZE_DATA3D 		64
#define D2C_EULERS_PER_BLOCK_DATA3D 32
// -------------------------------------


// FINE DIFF ---------------------------
#define D2F_BLOCK_SIZE_2D 			256
#define D2F_CHUNK_2D 				7

#define D2F_BLOCK_SIZE_REF3D 		256
#define D2F_CHUNK_REF3D 			7

#define D2F_BLOCK_SIZE_DATA3D 		512
#define D2F_CHUNK_DATA3D			4
// -------------------------------------


// WAVG --------------------------------
#define WAVG_BLOCK_SIZE_DATA3D	 	512
#define WAVG_BLOCK_SIZE 	 		256
// -------------------------------------


// MISC --------------------------------
#define SUMW_BLOCK_SIZE 	  	32
#define SOFTMASK_BLOCK_SIZE 	128
#define CFTT_BLOCK_SIZE 	 	128
#define PROBRATIO_BLOCK_SIZE 	128
#define POWERCLASS_BLOCK_SIZE 	128
#define PROJDIFF_CHUNK_SIZE 	14
// -------------------------------------

// RANDOMIZATION -----------------------
#define RND_BLOCK_NUM                   64
#define RND_BLOCK_SIZE                  32
// -------------------------------------


#define BACKPROJECTION4_BLOCK_SIZE 64
#define BACKPROJECTION4_GROUP_SIZE 16
#define BACKPROJECTION4_PREFETCH_COUNT 3
#define BP_2D_BLOCK_SIZE 128
#define BP_REF3D_BLOCK_SIZE 128
#define BP_DATA3D_BLOCK_SIZE 640


#define REF_GROUP_SIZE 3			// -- Number of references to be treated per block --
									// This applies to wavg and reduces global memory
									// accesses roughly proportionally, but scales shared
									// memory usage by allocating
									// ( 6*REF_GROUP_SIZE + 4 ) * BLOCK_SIZE XFLOATS. // DEPRECATED

#define NR_CLASS_MUTEXES 5

//The approximate minimum amount of memory each process occupies on a device (in MBs)
#define GPU_THREAD_MEMORY_OVERHEAD_MB 200

#endif /* CUDA_SETTINGS_H_ */

#ifndef CUDA_SETTINGS_H_
#define CUDA_SETTINGS_H_

// Required compute capability
#define CUDA_CC_MAJOR 3
#define CUDA_CC_MINOR 5

#define COMPLEXTEXTURE false
#define LAUNCH_CHECK
#define CUDA_BENCHMARK_OLD true

#ifdef CUDA_DOUBLE_PRECISION
	#define XFLOAT double
	#define CUDACOMPLEX double2
#else
	#define XFLOAT float
	#define CUDACOMPLEX float2
#endif

#ifdef RELION_SINGLE_PRECISION
	#define RFLOAT float
#else
	#define RFLOAT double
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

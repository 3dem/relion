#ifndef CUDA_SETTINGS_H_
#define CUDA_SETTINGS_H_

#define COMPLEXTEXTURE false
#define LAUNCH_CHECK

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

#define MAX_RESOL_SHARED_MEM 32
#define BLOCK_SIZE  128         	// -- Number of threads in a block --
									// This is optimally set as big as possible without
									// its ceil:ed multiple exceeding imagesize by too much.


#define D2C_BLOCK_SIZE_3D 		128
#define D2C_BLOCK_SIZE_2D 		512
#define D2C_EULERS_PER_BLOCK_3D 16
#define D2C_EULERS_PER_BLOCK_2D 4


#define WAVG_BLOCK_SIZE 	 	256
#define SUMW_BLOCK_SIZE 	  	32
#define SOFTMASK_BLOCK_SIZE 	1024
#define CFTT_BLOCK_SIZE 	 	128
#define PROBRATIO_BLOCK_SIZE 	128
#define POWERCLASS_BLOCK_SIZE 	128
#define PROJDIFF_CHUNK_SIZE 	14

#define REF_GROUP_SIZE 3			// -- Number of references to be treated per block --
									// This applies to wavg and reduces global memory
									// accesses roughly proportionally, but scales shared
									// memory usage by allocating
									// ( 6*REF_GROUP_SIZE + 4 ) * BLOCK_SIZE XFLOATS. // DEPRECATED

#define NR_CLASS_MUTEXES 5

//The approximate minimum amount of memory each process occupies on a device (in MBs)
#define MEMORY_OVERHEAD_MB 200

#endif /* CUDA_SETTINGS_H_ */

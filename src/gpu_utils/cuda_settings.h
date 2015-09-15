#ifndef CUDA_SETTINGS_H_
#define CUDA_SETTINGS_H_

#ifdef CUDA_DOUBLE_PRECISION
#define XFLOAT double
#else
#define XFLOAT float
#endif

#define MAX_RESOL_SHARED_MEM 32
#define BLOCK_SIZE  128         	// -- Number of threads in a block --
									// This is optimally set as big as possible without
									// its ceil:ed multiple exceeding imagesize by too much.
#define D2C_BLOCK_SIZE  128
#define WAVG_BLOCK_SIZE 256
#define SUMW_BLOCK_SIZE 32

#define PROJDIFF_CHUNK_SIZE 14

#define REF_GROUP_SIZE 3			// -- Number of references to be treated per block --
									// This applies to wavg and reduces global memory
									// accesses roughly proportionally, but scales shared
									// memory usage by allocating
									// ( 6*REF_GROUP_SIZE + 4 ) * BLOCK_SIZE XFLOATS. // DEPRECATED

#define NR_CLASS_MUTEXES 5


#endif /* CUDA_SETTINGS_H_ */

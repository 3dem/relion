/*
 * cuda_settings.h
 *
 *  Created on: May 26, 2015
 *      Author: bjornf
 */

#ifndef CUDA_SETTINGS_H_
#define CUDA_SETTINGS_H_

#ifdef CUDA_DOUBLE_PRECISION
#define FLOAT double
#else
#define FLOAT float
#endif
#define MAX_RESOL_SHARED_MEM 32
#define BLOCK_SIZE  128         	// -- Number of threads in a block --
									// This is optimally set as big as possible without
									// its ceil:ed multiple exceeding imagesize by too much.
#define SUM_BLOCK_SIZE 32

#define REF_GROUP_SIZE 3			// -- Number of references to be treated per block --
									// This applies to wavg and reduces global memory
									// accesses roughly proportionally, but scales shared
									// memory usage by allocating
									// ( 6*REF_GROUP_SIZE + 4 ) * BLOCK_SIZE
									// FLOATS.

#define NR_CLASS_MUTEXES 5


#endif /* CUDA_SETTINGS_H_ */

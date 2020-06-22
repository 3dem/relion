#ifndef ACC_SETTINGS_H_
#define ACC_SETTINGS_H_

#include "src/macros.h"

#ifdef ACC_DOUBLE_PRECISION
	#define XFLOAT double
	#ifndef _CUDA_ENABLED
typedef struct{ XFLOAT x; XFLOAT y;} double2;
	#endif
	#define ACCCOMPLEX double2
#else
	#define XFLOAT float
	#ifndef _CUDA_ENABLED
		typedef struct{ XFLOAT x; XFLOAT y;} float2;
	#endif
	#define ACCCOMPLEX float2
#endif
#ifdef ALTCPU
	#ifndef _CUDA_ENABLED
typedef float cudaStream_t;
		typedef double CudaCustomAllocator;
		#define cudaStreamPerThread 0
	#endif
#endif

#endif /* ACC_SETTINGS_H_ */

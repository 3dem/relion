#ifndef ACC_SETTINGS_H_
#define ACC_SETTINGS_H_

#include "src/macros.h"

#ifdef ACC_DOUBLE_PRECISION
	#define XFLOAT double
	#if !defined(_CUDA_ENABLED) && !defined(_HIP_ENABLED)
		typedef struct{ XFLOAT x; XFLOAT y;} double2;
	#endif
	#define ACCCOMPLEX double2
#else
	#define XFLOAT float
	#if !defined(_CUDA_ENABLED) && !defined(_HIP_ENABLED)
		typedef struct{ XFLOAT x; XFLOAT y;} float2;
	#endif
	#define ACCCOMPLEX float2
#endif

#ifdef _CUDA_ENABLED
	#define accGPUGetDeviceCount cudaGetDeviceCount
	#define accGPUDeviceProp cudaDeviceProp
	#define accGPUGetDeviceProperties cudaGetDeviceProperties
	#define accGPUDeviceSynchronize cudaDeviceSynchronize
	#define accGPUSetDevice cudaSetDevice
	#define accGPUMemGetInfo cudaMemGetInfo
	#define MlOptimiserAccGPU MlOptimiserCuda
	#define AutoPickerAccGPU AutoPickerCuda
#elif _HIP_ENABLED
	#define accGPUGetDeviceCount hipGetDeviceCount
	#define accGPUDeviceProp hipDeviceProp_t
	#define accGPUGetDeviceProperties hipGetDeviceProperties
	#define accGPUDeviceSynchronize hipDeviceSynchronize
	#define accGPUSetDevice hipSetDevice
	#define accGPUMemGetInfo hipMemGetInfo
	#define MlOptimiserAccGPU MlOptimiserHip
	#define AutoPickerAccGPU AutoPickerHip
#elif _SYCL_ENABLED
	#define MlOptimiserAccGPU MlOptimiserSYCL
#endif

#endif /* ACC_SETTINGS_H_ */

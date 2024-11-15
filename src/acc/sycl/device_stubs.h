#ifndef DEVICE_STUBS_H
#define DEVICE_STUBS_H

#undef CUDA
#undef _CUDA_ENABLED
#undef HIP
#undef _HIP_ENABLED

#include <cstddef>
#include "src/acc/sycl/sycl_virtual_dev.h"

using dim3 = int;
using deviceStream_t = virtualSYCL*;
using deviceCustomAllocator = double;

using cudaStream_t = virtualSYCL*;
using CudaCustomAllocator = double;
#define cudaStreamPerThread 0

using hipStream_t = virtualSYCL*;
using HipCustomAllocator = double;
#define hipStreamPerThread 0

#define CUSTOM_ALLOCATOR_REGION_NAME( name ) //Do nothing
#define LAUNCH_PRIVATE_ERROR(func, status)
#define LAUNCH_HANDLE_ERROR( err )
#define DEBUG_HANDLE_ERROR( err )
#define HANDLE_ERROR( err ) 

#endif

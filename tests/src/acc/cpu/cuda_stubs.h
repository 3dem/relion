#ifndef CUDA_STUBS_H
#define CUDA_STUBS_H

#undef CUDA
typedef float cudaStream_t;
typedef double CudaCustomAllocator;
typedef int dim3;
#define cudaStreamPerThread 0
#define CUSTOM_ALLOCATOR_REGION_NAME( name ) //Do nothing
#define LAUNCH_PRIVATE_ERROR(func, status)
#define LAUNCH_HANDLE_ERROR( err )
#define DEBUG_HANDLE_ERROR( err )
#define HANDLE_ERROR( err ) 
#endif
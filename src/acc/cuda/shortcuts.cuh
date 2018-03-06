#ifndef CUDA_SHORTCUTS_CUH_
#define CUDA_SHORTCUTS_CUH_

namespace CudaShortcuts
{

/**
 * Print cuda device memory info
 */
static void printMemInfo()
{
	size_t free;
	size_t total;
	DEBUG_HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
	float free_hr(free/(1024.*1024.));
	float total_hr(total/(1024.*1024.));
    printf( "free %.2fMiB, total %.2fMiB, used %.2fMiB\n",
    		free_hr, total_hr, total_hr - free_hr);
}

template< typename T>
static inline
void cpyHostToDevice( T *h_ptr, T *d_ptr, size_t size)
{
	DEBUG_HANDLE_ERROR(cudaMemcpy( d_ptr, h_ptr, size * sizeof(T), cudaMemcpyHostToDevice));
};

template< typename T>
static inline
void cpyHostToDevice( T *h_ptr, T *d_ptr, size_t size, cudaStream_t stream)
{
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( d_ptr, h_ptr, size * sizeof(T), cudaMemcpyHostToDevice, stream));
};

template< typename T>
static inline
void cpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size)
{
	DEBUG_HANDLE_ERROR(cudaMemcpy( h_ptr, d_ptr, size * sizeof(T), cudaMemcpyDeviceToHost));
};

template< typename T>
static inline
void cpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size, cudaStream_t &stream)
{
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( h_ptr, d_ptr, size * sizeof(T), cudaMemcpyDeviceToHost, stream));
};

template< typename T>
static inline
void cpyDeviceToDevice( T *src, T *des, size_t size, cudaStream_t &stream)
{
	DEBUG_HANDLE_ERROR(cudaMemcpyAsync( des, src, size * sizeof(T), cudaMemcpyDeviceToDevice, stream));
};

template< typename T>
static inline
void memInit( T *ptr, T value, size_t size)
{
	DEBUG_HANDLE_ERROR(cudaMemset( ptr, value, size * sizeof(T)));
};

template< typename T>
static inline
void memInit( T *ptr, T value, size_t size, cudaStream_t &stream)
{
	DEBUG_HANDLE_ERROR(cudaMemsetAsync( ptr, value, size * sizeof(T), stream));
};

}
#endif //CUDA_SHORTCUTS_CUH_

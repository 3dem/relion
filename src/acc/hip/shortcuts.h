/* Portions of this code are under:
   Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/
#ifndef HIP_SHORTCUTS_H_
#define HIP_SHORTCUTS_H_

namespace HipShortcuts
{

/**
 * Print hip device memory info
 */
static void printMemInfo()
{
	size_t free;
	size_t total;
	DEBUG_HANDLE_ERROR(hipMemGetInfo( &free, &total ));
	float free_hr(free/(1024.*1024.));
	float total_hr(total/(1024.*1024.));
    printf( "free %.2fMiB, total %.2fMiB, used %.2fMiB\n",
    		free_hr, total_hr, total_hr - free_hr);
}

template< typename T>
static inline
void cpyHostToDevice( T *h_ptr, T *d_ptr, size_t size)
{
	DEBUG_HANDLE_ERROR(hipMemcpy( d_ptr, h_ptr, size * sizeof(T), hipMemcpyHostToDevice));
};

template< typename T>
static inline
void cpyHostToDevice( T *h_ptr, T *d_ptr, size_t size, hipStream_t stream)
{
	DEBUG_HANDLE_ERROR(hipMemcpyAsync( d_ptr, h_ptr, size * sizeof(T), hipMemcpyHostToDevice, stream));
};

template< typename T>
static inline
void cpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size)
{
	DEBUG_HANDLE_ERROR(hipMemcpy( h_ptr, d_ptr, size * sizeof(T), hipMemcpyDeviceToHost));
};

template< typename T>
static inline
void cpyDeviceToHost( T *d_ptr, T *h_ptr, size_t size, hipStream_t &stream)
{
	DEBUG_HANDLE_ERROR(hipMemcpyAsync( h_ptr, d_ptr, size * sizeof(T), hipMemcpyDeviceToHost, stream));
};

template< typename T>
static inline
void cpyDeviceToDevice( T *src, T *des, size_t size, hipStream_t &stream)
{
	DEBUG_HANDLE_ERROR(hipMemcpyAsync( des, src, size * sizeof(T), hipMemcpyDeviceToDevice, stream));
};

template< typename T>
static inline
void memInit( T *ptr, T value, size_t size)
{
	DEBUG_HANDLE_ERROR(hipMemset( ptr, value, size * sizeof(T)));
};

template< typename T>
static inline
void memInit( T *ptr, T value, size_t size, hipStream_t &stream)
{
	DEBUG_HANDLE_ERROR(hipMemsetAsync( ptr, value, size * sizeof(T), stream));
};

}
#endif //HIP_SHORTCUTS_H_

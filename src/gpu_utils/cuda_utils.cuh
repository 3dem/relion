#include <cuda.h>

static void HandleError( cudaError_t err,
                         const char *file,
                         int line )
{
    if (err != cudaSuccess)
    {
        printf( "CUDA ERROR: %s in %s at line %d\n", cudaGetErrorString( err ),
                file, line );
        exit( EXIT_FAILURE );
    }
}

#ifdef DEBUG_CUDA
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#else
#define HANDLE_ERROR( err ) (err) //Do nothing
#endif

static void cudaPrintMemInfo()
{
	size_t free;
	size_t total;
	HANDLE_ERROR(cudaMemGetInfo( &free, &total ));
	float free_hr(free/(8*1024*1024));
	float total_hr(total/(8*1024*1024));
    printf( "free %.2fMB, total %.2fMB, used %.2fMB\n",
    		free_hr, total_hr, total_hr - free_hr);
}

/*
 * Wraps CUDA functions for memory management on the device only
 */
template <typename T>
class CudaDevicePtr
{
public:
	const size_t sz;
	T *d;
	bool d_free;

	inline
	__host__ CudaDevicePtr<T>(T * d_start, size_t size):
		sz(size), d(d_start), d_free(false)
	{};

	inline
	__host__ CudaDevicePtr<T>(size_t size):
		sz(size), d(0), d_free(true)
	{
		HANDLE_ERROR(cudaMalloc( (void**) &d, size * sizeof(T)));
	};

	inline
	__host__ void device_init(int value)
	{
#ifdef DEBUG_CUDA
		if (d == 0)
		{
			printf("Memset requested before allocation in %s at line %d\n.",
					__FILE__, __LINE__);
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaMemset(d, value, sz * sizeof(T)));
	}

	inline
	__device__ T& operator[](std::size_t idx)       { return d[idx]; };
	inline
	__device__ const T& operator[](std::size_t idx) const { return d[idx]; };

	inline
	__host__ void free()
	{
#ifdef DEBUG_CUDA
		if (d == 0)
		{
			printf("Free device memory was called on NULL pointer in %s at line %d\n.",
					__FILE__, __LINE__);
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaFree(d));
		d = 0;
		d_free = false;
	}

	inline
	__host__ ~CudaDevicePtr()
	{
		if (d_free) free();
	}
};

template <typename T>
class CudaGlobalPtr
{
public:
	const size_t sz;
	T *h, *d;
	bool h_free, d_free;

	inline
	__host__ CudaGlobalPtr<T>(T * h_start, size_t size):
		sz(size), h(h_start), d(0), h_free(false), d_free(false)
	{};

	inline
	__host__ CudaGlobalPtr<T>(size_t size):
		sz(size), h(new T[size]), d(0), h_free(true), d_free(false)
	{};

	inline
	__host__ CudaGlobalPtr<T>(CudaDevicePtr<T> *ptr):
		sz(ptr->size), h(new T[sz]), d(ptr->d), h_free(true), d_free(false)
	{
		cp_to_host();
	};

	inline
	__host__ void device_alloc()
	{
		HANDLE_ERROR(cudaMalloc( (void**) &d, sz * sizeof(T)));
		d_free = true;
	}

	inline
	__host__ void device_init(int value)
	{
#ifdef DEBUG_CUDA
		if (d == 0)
		{
			printf("Memset requested before allocation in %s at line %d\n.",
					__FILE__, __LINE__);
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaMemset( d, value, sz * sizeof(T)));
	}

	inline
	__host__ void cp_to_device()
	{
#ifdef DEBUG_CUDA
		if (d == 0)
		{
			printf("Cpy to device requested before allocation in %s at line %d\n.",
					__FILE__, __LINE__);
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaMemcpy( d, h, sz * sizeof(T), cudaMemcpyHostToDevice));
	}

	inline
	__host__ void cp_to_host()
	{
		HANDLE_ERROR(cudaMemcpy( h, d, sz * sizeof(T), cudaMemcpyDeviceToHost ));
	}

	inline
	__host__ T& operator[](size_t idx)       { return h[idx]; };
	inline
	__host__ const T& operator[](size_t idx) const { return h[idx]; };

	inline
	__host__ void free_device()
	{
#ifdef DEBUG_CUDA
		if (d == 0)
		{
			printf("Free device memory was called on NULL pointer in %s at line %d\n.",
					__FILE__, __LINE__);
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaFree(d));
		d = 0;
		d_free = false;
	}

	inline
	__host__ void free_host()
	{
		delete [] h;
		h = 0;
		h_free = false;
	}

	inline
	__host__ void free()
	{
		free_device();
		free_host();
	}

	inline
	__host__ ~CudaGlobalPtr()
	{
		if (d_free) free_device();
		if (h_free) free_host();
	}
};

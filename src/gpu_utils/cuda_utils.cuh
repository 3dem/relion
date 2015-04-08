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
	float free_hr(free/(1024.*1024.));
	float total_hr(total/(1024.*1024.));
    printf( "free %.2fMiB, total %.2fMiB, used %.2fMiB\n",
    		free_hr, total_hr, total_hr - free_hr);
}

template <typename T>
class CudaGlobalPtr
{
public:
	size_t sz;
	T *h, *d;
	bool h_free, d_free;

	inline
	__host__ CudaGlobalPtr<T>():
		sz(0), h(0), d(0), h_free(false), d_free(false)
	{};

	inline
	__host__ CudaGlobalPtr<T>(T * h_start, size_t size):
		sz(size), h(h_start), d(0), h_free(false), d_free(false)
	{};

	inline
	__host__ CudaGlobalPtr<T>(size_t size):
		sz(size), h(new T[size]), d(0), h_free(true), d_free(false)
	{};

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
			printf("Memset requested before allocation in device_init().\n");
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
			printf("Cpy to device requested before allocation in cp_to_device().\n");
			exit( EXIT_FAILURE );
		}
		if (h == 0)
		{
			printf("NULL host pointer in cp_to_device().\n");
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaMemcpy( d, h, sz * sizeof(T), cudaMemcpyHostToDevice));
	}

	inline
	__host__ void cp_to_device(T * hostPtr)
	{
		h = hostPtr;
		cp_to_device();
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
			printf("Free device memory was called on NULL pointer in free_device().\n");
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

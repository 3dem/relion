
#include <cuda.h>

#ifdef CUDA_DOUBLE_PRECISION
#define FLOAT double
__forceinline__ __device__ FLOAT cuda_atomic_add(FLOAT* address, FLOAT val)
{
	unsigned long long int* address_as_ull = (unsigned long long int*)address;
	unsigned long long int old = *address_as_ull, assumed;
	do
	{
		assumed = old;
		old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)));
	}
	while (assumed != old); // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	return __longlong_as_double(old);
}
#else
#define FLOAT float
__forceinline__ __device__ FLOAT cuda_atomic_add(FLOAT* address, FLOAT val)
{
	unsigned int* address_as_ul = (unsigned int*)address;
	unsigned int old = *address_as_ul, assumed;
	do
	{
		assumed = old;
		old = atomicCAS(address_as_ul, assumed, (int) (val + __longlong_as_double(assumed)));
	}
	while (assumed != old); // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
	return (FLOAT) old;
}
#endif


class CudaComplex
{
public:
	FLOAT real, imag;

	__device__ __host__ CudaComplex(): real(), imag() {};
	__device__ __host__ CudaComplex(FLOAT real, FLOAT imag): real(real), imag(imag) {};

public:
	__device__ CudaComplex operator-(CudaComplex a)  {
	     return CudaComplex(real-a.real,imag-a.imag);
	   }
	__device__ CudaComplex operator+(CudaComplex a)  {
	     return CudaComplex(real+a.real,imag+a.imag);
	   }
	__device__ CudaComplex operator*(FLOAT s)  {
	     return CudaComplex(real*s,imag*s);
	   }
	__device__ void operator +=(CudaComplex s)  {
	     real+=s.real;imag+=s.imag;
	   }

};

#ifdef DEBUG_CUDA
#define HANDLE_ERROR( err ) (HandleError( err, __FILE__, __LINE__ ))
#else
#define HANDLE_ERROR( err ) (err) //Do nothing
#endif
static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "DEBUG_ERROR: %s in %s at line %d\n",
        		cudaGetErrorString( err ), file, line );
        exit( EXIT_FAILURE );
    }
}

/**
 * Print cuda device memory info
 */
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

//Non-concurrent benchmarking tools (only for Linux)
#ifdef CUDA_BENCHMARK
#include <vector>
#include <ctime>
#include <string>

static int cuda_benchmark_find_id(std::string id, std::vector<std::string> v)
{
	for (unsigned i = 0; i < v.size(); i++)
		if (v[i] == id)
			return i;
	return -1;
}

std::vector<std::string> cuda_cpu_identifiers;
std::vector<clock_t>     cuda_cpu_start_times;

#define CUDA_CPU_TIC(ID) (cuda_cpu_tic(ID))
static void cuda_cpu_tic(std::string id)
{
	if (cuda_benchmark_find_id(id, cuda_cpu_identifiers) == -1)
	{
		cuda_cpu_identifiers.push_back(id);
		cuda_cpu_start_times.push_back(clock());
	}
	else
	{
		printf("DEBUG_ERROR: Provided identifier '%s' already exists in call to cuda_cpu_tic.\n", id.c_str());
		exit( EXIT_FAILURE );
	}
}

#define CUDA_CPU_TOC(ID) (cuda_cpu_toc(ID))
static void cuda_cpu_toc(std::string id)
{
	int idx = cuda_benchmark_find_id(id, cuda_cpu_identifiers);
	if (idx == -1)
	{
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to cuda_cpu_toc.\n", id.c_str());
		exit( EXIT_FAILURE );
	}
	else
	{
		clock_t start_time = cuda_cpu_start_times[idx];
		cuda_cpu_identifiers.erase(cuda_cpu_identifiers.begin()+idx);
		cuda_cpu_start_times.erase(cuda_cpu_start_times.begin()+idx);
		FILE *fPtr = fopen("benchmark.dat","a");
		fprintf(fPtr,"CPU: %s \t %.2f ms\n", id.c_str(),
				(((float)clock() - (float)start_time) / CLOCKS_PER_SEC ) * 1000.);
		fclose(fPtr);
	}
}
std::vector<std::string> cuda_gpu_kernel_identifiers;
std::vector<cudaEvent_t> cuda_gpu_kernel_start_times;
std::vector<cudaEvent_t> cuda_gpu_kernel_stop_times;

#define CUDA_GPU_TIC(ID) (cuda_gpu_tic(ID))
static void cuda_gpu_tic(std::string id)
{
	if (cuda_benchmark_find_id(id, cuda_gpu_kernel_identifiers) == -1)
	{
		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);
		cuda_gpu_kernel_identifiers.push_back(id);
		cuda_gpu_kernel_start_times.push_back(start);
		cuda_gpu_kernel_stop_times.push_back(stop);
	}
	else
	{
		printf("DEBUG_ERROR: Provided identifier '%s' already exists in call to cuda_gpu_tic.\n",
				id.c_str());
		exit( EXIT_FAILURE );
	}
}

#define CUDA_GPU_TAC(ID) (cuda_gpu_tac(ID))
static void cuda_gpu_tac(std::string id)
{
	int idx = cuda_benchmark_find_id(id, cuda_gpu_kernel_identifiers);
	if (idx == -1)
	{
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to cuda_gpu_tac.\n",
				id.c_str());
		exit( EXIT_FAILURE );
	}
	else
	{
		cudaEventRecord(cuda_gpu_kernel_stop_times[idx], 0);
		cudaEventSynchronize(cuda_gpu_kernel_stop_times[idx]);
	}
}

#define CUDA_GPU_TOC(ID) (cuda_gpu_toc(ID))
static void cuda_gpu_toc(std::string id)
{
	int idx = cuda_benchmark_find_id(id, cuda_gpu_kernel_identifiers);
	if (idx == -1)
	{
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to cuda_gpu_toc.\n",
				id.c_str());
		exit( EXIT_FAILURE );
	}
	else
	{
		float time;
		cudaEventElapsedTime(&time, cuda_gpu_kernel_start_times[idx],
				cuda_gpu_kernel_stop_times[idx]);
		cudaEventDestroy(cuda_gpu_kernel_start_times[idx]);
		cudaEventDestroy(cuda_gpu_kernel_stop_times[idx]);
		cuda_gpu_kernel_identifiers.erase(cuda_gpu_kernel_identifiers.begin()+idx);
		cuda_gpu_kernel_start_times.erase(cuda_gpu_kernel_start_times.begin()+idx);
		cuda_gpu_kernel_stop_times.erase(cuda_gpu_kernel_stop_times.begin()+idx);

		FILE *fPtr = fopen("benchmark.dat","a");
		fprintf(fPtr,"GPU: %s \t %.2f ms\n", id.c_str(), time);
		fclose(fPtr);
	}
}
#else
#define CUDA_CPU_TIC(ID)
#define CUDA_CPU_TOC(ID)
#define CUDA_GPU_TIC(ID)
#define CUDA_GPU_TAC(ID)
#define CUDA_GPU_TOC(ID)
#endif

template <typename T>
class CudaGlobalPtr
{
public:
	size_t size; //Size used when copying data from and to device
	T *h_ptr, *d_ptr; //Host and device pointers
	bool h_do_free, d_do_free; //True if host or device needs to be freed

	inline
	__host__ CudaGlobalPtr<T>():
		size(0), h_ptr(0), d_ptr(0), h_do_free(false), d_do_free(false)
	{};

	inline
	__host__ CudaGlobalPtr<T>(T * h_start, size_t size):
		size(size), h_ptr(h_start), d_ptr(0), h_do_free(false), d_do_free(false)
	{};

	inline
	__host__ CudaGlobalPtr<T>(size_t size):
		size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true), d_do_free(false)
	{};

	/**
	 * Allocate memory on device
	 */
	inline
	__host__ void device_alloc()
	{
#ifdef DEBUG_CUDA
		if (d_do_free)
			printf("DEBUG_WARNING: Device double allocation.\n");
#endif
		d_do_free = true;
		HANDLE_ERROR(cudaMalloc( (void**) &d_ptr, size * sizeof(T)));
	}

	/**
	 * Allocate memory on host
	 */
	inline
	__host__ void host_alloc()
	{
#ifdef DEBUG_CUDA
		if (h_do_free)
			printf("DEBUG_WARNING: Host double allocation.\n");
#endif
		h_do_free = true;
		h_ptr = new T[size];
	}

	/**
	 * Initiate device memory with provided value
	 */
	inline
	__host__ void device_init(int value)
	{
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: Memset requested before allocation in device_init().\n");
#endif
		HANDLE_ERROR(cudaMemset( d_ptr, value, size * sizeof(T)));
	}

	/**
	 * Copy a number (size) of bytes to device stored in the host pointer
	 */
	inline
	__host__ void cp_to_device()
	{
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: cp_to_device() called before allocation.\n");
		if (h_ptr == 0)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_device().\n");
#endif
		HANDLE_ERROR(cudaMemcpy( d_ptr, h_ptr, size * sizeof(T), cudaMemcpyHostToDevice));
	}

	/**
	 * Copy a number (size) of bytes to device stored in the provided host pointer
	 */
	inline
	__host__ void cp_to_device(T * hostPtr)
	{
#ifdef DEBUG_CUDA
		if (h_ptr != 0)
			printf("DEBUG_WARNING: Host pointer already set in call to cp_to_device(hostPtr).\n");
#endif
		h_ptr = hostPtr;
		cp_to_device();
	}

	/**
	 * Copy a number (size) of bytes from device to the host pointer
	 */
	inline
	__host__ void cp_to_host()
	{
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: cp_to_host() called before allocation.\n");
		if (h_ptr == 0)
			printf("DEBUG_WARNING: NULL host pointer in cp_to_host().\n");
#endif
		HANDLE_ERROR(cudaMemcpy( h_ptr, d_ptr, size * sizeof(T), cudaMemcpyDeviceToHost ));
	}

	/**
	 * Host data quick access
	 */
	inline
	__host__ T& operator[](size_t idx) { return h_ptr[idx]; };


	/**
	 * Host data quick access
	 */
	inline
	__host__ const T& operator[](size_t idx) const { return h_ptr[idx]; };

	/**
	 * Device pointer quick access
	 */
	inline
	__host__ T* operator~() {
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: \"kernel cast\" on null pointer.\n");
#endif
		return d_ptr;
	};

	/**
	 * Delete device data
	 */
	inline
	__host__ void free_device()
	{
#ifdef DEBUG_CUDA
		if (d_ptr == 0)
			printf("DEBUG_WARNING: Free device memory was called on NULL pointer in free_device().\n");
#endif
		d_do_free = false;
		HANDLE_ERROR(cudaFree(d_ptr));
		d_ptr = 0;
	}

	/**
	 * Delete host data
	 */
	inline
	__host__ void free_host()
	{
#ifdef DEBUG_CUDA
		if (h_ptr == 0)
		{
			printf("DEBUG_ERROR: free_host() called on NULL pointer.\n");
	        exit( EXIT_FAILURE );
		}
#endif
		h_do_free = false;
		delete [] h_ptr;
		h_ptr = 0;
	}

	/**
	 * Delete both device and host data
	 */
	inline
	__host__ void free()
	{
		free_device();
		free_host();
	}

	inline
	__host__ ~CudaGlobalPtr()
	{
		if (d_do_free) free_device();
		if (h_do_free) free_host();
	}
};

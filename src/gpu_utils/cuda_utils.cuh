#include <cuda.h>

static void HandleError( cudaError_t err, const char *file, int line )
{
    if (err != cudaSuccess)
    {
        printf( "CUDA ERROR: %s in %s at line %d\n",
        		cudaGetErrorString( err ), file, line );
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

//Non-concurrent benchmarking tools (only for Linux)
#if defined(CUDA_BENCHMARK) && ! defined(__linux__)
#endif
#if defined(CUDA_BENCHMARK) && defined(__linux__)
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
		printf("Provided identifier '%s' already exists in call to cuda_cpu_tic.\n", id.c_str());
		exit( EXIT_FAILURE );
	}
}

#define CUDA_CPU_TOC(ID) (cuda_cpu_toc(ID))
static void cuda_cpu_toc(std::string id)
{
	int idx = cuda_benchmark_find_id(id, cuda_cpu_identifiers);
	if (idx == -1)
	{
		printf("Provided identifier '%s' not found in call to cuda_cpu_toc.\n", id.c_str());
		exit( EXIT_FAILURE );
	}
	else
	{
		clock_t start_time = cuda_cpu_start_times[idx];
		cuda_cpu_identifiers.erase(cuda_cpu_identifiers.begin()+idx);
		cuda_cpu_start_times.erase(cuda_cpu_start_times.begin()+idx);
		FILE *fPtr = fopen("benchmark.dat","a");
		fprintf(fPtr,"CPU: %s \t %.2f \xC2\xB5s\n", id.c_str(),
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
		printf("Provided identifier '%s' already exists in call to cuda_gpu_tic.\n",
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
		printf("Provided identifier '%s' not found in call to cuda_gpu_tac.\n",
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
		printf("Provided identifier '%s' not found in call to cuda_gpu_toc.\n",
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
		fprintf(fPtr,"GPU: %s \t %.2f \xC2\xB5s\n", id.c_str(), time);
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
	size_t size;
	T *hPtr, *dPtr;
	bool h_free, d_free; //True if host or device needs to be freed

	inline
	__host__ CudaGlobalPtr<T>():
		size(0), hPtr(0), dPtr(0), h_free(false), d_free(false)
	{};

	inline
	__host__ CudaGlobalPtr<T>(T * h_start, size_t size):
		size(size), hPtr(h_start), dPtr(0), h_free(false), d_free(false)
	{};

	inline
	__host__ CudaGlobalPtr<T>(size_t size):
		size(size), hPtr(new T[size]), dPtr(0), h_free(true), d_free(false)
	{};

	inline
	__host__ void device_alloc()
	{
		HANDLE_ERROR(cudaMalloc( (void**) &dPtr, size * sizeof(T)));
		d_free = true;
	}

	inline
	__host__ void device_init(int value)
	{
#ifdef DEBUG_CUDA
		if (dPtr == 0)
		{
			printf("Memset requested before allocation in device_init().\n");
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaMemset( dPtr, value, size * sizeof(T)));
	}

	inline
	__host__ void cp_to_device()
	{
#ifdef DEBUG_CUDA
		if (dPtr == 0)
		{
			printf("Cpy to device requested before allocation in cp_to_device().\n");
			exit( EXIT_FAILURE );
		}
		if (hPtr == 0)
		{
			printf("NULL host pointer in cp_to_device().\n");
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaMemcpy( dPtr, hPtr, size * sizeof(T), cudaMemcpyHostToDevice));
	}

	inline
	__host__ void cp_to_device(T * hostPtr)
	{
		hPtr = hostPtr;
		cp_to_device();
	}

	inline
	__host__ void cp_to_host()
	{
		HANDLE_ERROR(cudaMemcpy( hPtr, dPtr, size * sizeof(T), cudaMemcpyDeviceToHost ));
	}

	inline
	__host__ T& operator[](size_t idx) { return hPtr[idx]; };
	inline
	__host__ const T& operator[](size_t idx) const { return hPtr[idx]; };

	inline
	__host__ T* operator~() { return dPtr; };

	inline
	__host__ void free_device()
	{
#ifdef DEBUG_CUDA
		if (dPtr == 0)
		{
			printf("Free device memory was called on NULL pointer in free_device().\n");
			exit( EXIT_FAILURE );
		}
#endif
		HANDLE_ERROR(cudaFree(dPtr));
		dPtr = 0;
		d_free = false;
	}

	inline
	__host__ void free_host()
	{
		delete [] hPtr;
		hPtr = 0;
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

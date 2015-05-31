
#ifndef CUDA_BENCHMARK_UTILS_CUH_
#define CUDA_BENCHMARK_UTILS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/gpu_utils/cuda_utils.cuh"

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

#endif /* CUDA_BENCHMARK_UTILS_CUH_ */

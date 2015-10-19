
#ifndef CUDA_BENCHMARK_UTILS_CUH_
#define CUDA_BENCHMARK_UTILS_CUH_

#include <cuda_runtime.h>
#include <vector>
#include <iostream>
#include <fstream>

//Non-concurrent benchmarking tools (only for Linux)
#if defined CUDA_BENCHMARK_OLD
#include <vector>
#include <time.h>
#include <string>

static int cuda_benchmark_find_id(std::string id, std::vector<std::string> v)
{
	for (unsigned i = 0; i < v.size(); i++)
		if (v[i] == id)
			return i;
	return -1;
}

std::vector<std::string> cuda_cpu_benchmark_identifiers;
std::vector<clock_t>     cuda_cpu_benchmark_start_times;
FILE *cuda_cpu_benchmark_fPtr = fopen("benchmark_cpu.dat","w");

#define CUDA_CPU_TIC(ID) (cuda_cpu_tic(ID))
static void cuda_cpu_tic(std::string id)
{
	if (cuda_benchmark_find_id(id, cuda_cpu_benchmark_identifiers) == -1)
	{
		cuda_cpu_benchmark_identifiers.push_back(id);
		cuda_cpu_benchmark_start_times.push_back(clock());
	}
	else
	{
		printf("DEBUG_ERROR: Provided identifier '%s' already exists in call to cuda_cpu_tic.\n", id.c_str());
		raise(SIGSEGV);
	}
}

#define CUDA_CPU_TOC(ID) (cuda_cpu_toc(ID))
static void cuda_cpu_toc(std::string id)
{
	int idx = cuda_benchmark_find_id(id, cuda_cpu_benchmark_identifiers);
	if (idx == -1)
	{
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to cuda_cpu_toc.\n", id.c_str());
		//exit( EXIT_FAILURE );
	}
	else
	{
		clock_t t = clock() - cuda_cpu_benchmark_start_times[idx];
		cuda_cpu_benchmark_identifiers.erase(cuda_cpu_benchmark_identifiers.begin()+idx);
		cuda_cpu_benchmark_start_times.erase(cuda_cpu_benchmark_start_times.begin()+idx);
		fprintf(cuda_cpu_benchmark_fPtr,"%06.2f ms ......", (float)t / CLOCKS_PER_SEC * 1000.);
		for (int i = 1; i < cuda_cpu_benchmark_identifiers.size(); i++)
			fprintf(cuda_cpu_benchmark_fPtr,"......");
		fprintf(cuda_cpu_benchmark_fPtr," %s\n", id.c_str());
//		printf(,"%s \t %.2f ms\n", id.c_str(), (float)t / CLOCKS_PER_SEC * 1000.);
	}
}
std::vector<std::string> cuda_gpu_benchmark_identifiers;
std::vector<cudaEvent_t> cuda_gpu_benchmark_start_times;
std::vector<cudaEvent_t> cuda_gpu_benchmark_stop_times;
FILE *cuda_gpu_benchmark_fPtr = fopen("benchmark_gpu.dat","w");

#define CUDA_GPU_TIC(ID) (cuda_gpu_tic(ID))
static void cuda_gpu_tic(std::string id)
{
	if (cuda_benchmark_find_id(id, cuda_gpu_benchmark_identifiers) == -1)
	{
		cudaEvent_t start, stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		cudaEventRecord(start, 0);
		cuda_gpu_benchmark_identifiers.push_back(id);
		cuda_gpu_benchmark_start_times.push_back(start);
		cuda_gpu_benchmark_stop_times.push_back(stop);
	}
	else
	{
		printf("DEBUG_ERROR: Provided identifier '%s' already exists in call to cuda_gpu_tic.\n",
				id.c_str());
		raise(SIGSEGV);
	}
}

#define CUDA_GPU_TAC(ID) (cuda_gpu_tac(ID))
static void cuda_gpu_tac(std::string id)
{
	int idx = cuda_benchmark_find_id(id, cuda_gpu_benchmark_identifiers);
	if (idx == -1)
	{
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to cuda_gpu_tac.\n",
				id.c_str());
		raise(SIGSEGV);
	}
	else
	{
		cudaEventRecord(cuda_gpu_benchmark_stop_times[idx], 0);
		cudaEventSynchronize(cuda_gpu_benchmark_stop_times[idx]);
	}
}

#define CUDA_GPU_TOC() (cuda_gpu_toc())
static void cuda_gpu_toc()
{
	if (cuda_gpu_benchmark_identifiers.size() == 0)
	{
		printf("DEBUG_ERROR: There were no identifiers found in the list, on call to cuda_gpu_toc.\n");
		raise(SIGSEGV);
	}
	else
	{
		float time;
		for (int idx = 0; idx < cuda_gpu_benchmark_identifiers.size(); idx ++)
		{
			cudaEventElapsedTime(&time, cuda_gpu_benchmark_start_times[idx],
					cuda_gpu_benchmark_stop_times[idx]);
			cudaEventDestroy(cuda_gpu_benchmark_start_times[idx]);
			cudaEventDestroy(cuda_gpu_benchmark_stop_times[idx]);
			fprintf(cuda_gpu_benchmark_fPtr,"%.2f ms \t %s\n",
					time, cuda_gpu_benchmark_identifiers[idx].c_str());
		}

		cuda_gpu_benchmark_identifiers.clear();
		cuda_gpu_benchmark_start_times.clear();
		cuda_gpu_benchmark_stop_times.clear();
	}
}
#elif defined CUDA_PROFILING
#include <nvToolsExt.h>

#define CUDA_CPU_TIC(ID) (nvtxRangePush(ID))
#define CUDA_CPU_TOC(ID) (nvtxRangePop())
#define CUDA_GPU_TIC(ID)
#define CUDA_GPU_TAC(ID)
#define CUDA_GPU_TOC(ID)
#else
#define CUDA_CPU_TIC(ID)
#define CUDA_CPU_TOC(ID)
#define CUDA_GPU_TIC(ID)
#define CUDA_GPU_TAC(ID)
#define CUDA_GPU_TOC(ID)
#endif

#endif /* CUDA_BENCHMARK_UTILS_CUH_ */

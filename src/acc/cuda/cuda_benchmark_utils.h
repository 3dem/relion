
#ifndef CUDA_BENCHMARK_UTILS_H_
#define CUDA_BENCHMARK_UTILS_H_

//Non-concurrent benchmarking tools (only for Linux)

#include <cuda_runtime.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string>
#include <sstream>

#ifdef TIMING_FILES
#define	CTIC(timer,timing) (timer.cuda_cpu_tic(timing))
#define	CTOC(timer,timing) (timer.cuda_cpu_toc(timing))
#define	GTIC(timer,timing) (timer.cuda_gpu_tic(timing))
#define	GTOC(timer,timing) (timer.cuda_gpu_toc(timing))
#define	GATHERGPUTIMINGS(timer) (timer.cuda_gpu_printtictoc())
#elif defined CUDA_PROFILING
	#include <nvToolsExt.h>
	#define	CTIC(timer,timing) (nvtxRangePush(timing))
	#define	CTOC(timer,timing) (nvtxRangePop())
	#define	GTIC(timer,timing)
	#define	GTOC(timer,timing)
	#define	GATHERGPUTIMINGS(timer)
#else
	#define	CTIC(timer,timing)
	#define	CTOC(timer,timing)
	#define	GTIC(timer,timing)
	#define	GTOC(timer,timing)
	#define	GATHERGPUTIMINGS(timer)
#endif

class relion_timer
{

public:

std::vector<std::string> cuda_cpu_benchmark_identifiers;
std::vector<clock_t>     cuda_cpu_benchmark_start_times;
FILE *cuda_cpu_benchmark_fPtr;

std::vector<std::string> cuda_gpu_benchmark_identifiers;
std::vector<cudaEvent_t> cuda_gpu_benchmark_start_times;
std::vector<cudaEvent_t> cuda_gpu_benchmark_stop_times;
FILE *cuda_gpu_benchmark_fPtr;

relion_timer(std::string fnm)
{
	std::stringstream fnm_cpu, fnm_gpu;
	fnm_cpu << "output/" << fnm << "_cpu.dat";
	cuda_cpu_benchmark_fPtr = fopen(fnm_cpu.str().c_str(),"a");
	fnm_gpu << "output/" << fnm << "_gpu.dat";
	cuda_gpu_benchmark_fPtr = fopen(fnm_gpu.str().c_str(),"a");
}

int cuda_benchmark_find_id(std::string id, std::vector<std::string> v);

void cuda_cpu_tic(std::string id);

void cuda_cpu_toc(std::string id);

void cuda_gpu_tic(std::string id);

void cuda_gpu_toc(std::string id);

void cuda_gpu_printtictoc();

};

#endif /* CUDA_BENCHMARK_UTILS_H_ */

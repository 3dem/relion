/* Portions of this code are under:
   Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/
#ifndef HIP_BENCHMARK_UTILS_H_
#define HIP_BENCHMARK_UTILS_H_

//Non-concurrent benchmarking tools (only for Linux)

#include <hip/hip_runtime.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <time.h>
#include <string>
#include <sstream>

#ifdef TIMING_FILES
#define	CTIC(timer,timing) (timer.hip_cpu_tic(timing))
#define	CTOC(timer,timing) (timer.hip_cpu_toc(timing))
#define	GTIC(timer,timing) (timer.hip_gpu_tic(timing))
#define	GTOC(timer,timing) (timer.hip_gpu_toc(timing))
#define	GATHERGPUTIMINGS(timer) (timer.hip_gpu_printtictoc())
#elif defined HIP_PROFILING
	#include <roctracer/roctx.h>
	#define	CTIC(timer,timing) (roctxRangePush(timing))
	#define	CTOC(timer,timing) (roctxRangePop())
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

std::vector<std::string> hip_cpu_benchmark_identifiers;
std::vector<clock_t>     hip_cpu_benchmark_start_times;
FILE *hip_cpu_benchmark_fPtr;

std::vector<std::string> hip_gpu_benchmark_identifiers;
std::vector<hipEvent_t> hip_gpu_benchmark_start_times;
std::vector<hipEvent_t> hip_gpu_benchmark_stop_times;
FILE *hip_gpu_benchmark_fPtr;

relion_timer(std::string fnm)
{
	std::stringstream fnm_cpu, fnm_gpu;
	fnm_cpu << "output/" << fnm << "_cpu.dat";
	hip_cpu_benchmark_fPtr = fopen(fnm_cpu.str().c_str(),"a");
	fnm_gpu << "output/" << fnm << "_gpu.dat";
	hip_gpu_benchmark_fPtr = fopen(fnm_gpu.str().c_str(),"a");
}

int hip_benchmark_find_id(std::string id, std::vector<std::string> v);

void hip_cpu_tic(std::string id);

void hip_cpu_toc(std::string id);

void hip_gpu_tic(std::string id);

void hip_gpu_toc(std::string id);

void hip_gpu_printtictoc();

};

#endif /* HIP_BENCHMARK_UTILS_H_ */

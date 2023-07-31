/* Portions of this code are under:
   Copyright (c) 2022 Advanced Micro Devices, Inc. All rights reserved.
*/

#include "src/acc/hip/hip_benchmark_utils.h"

//Non-concurrent benchmarking tools (only for Linux)
#include <vector>
#include <time.h>
#include <string>
#include <signal.h>

#include "src/macros.h"
#include "src/error.h"

int relion_timer::hip_benchmark_find_id(std::string id, std::vector<std::string> v)
{
	for (unsigned i = 0; i < v.size(); i++)
		if (v[i] == id)
			return i;
	return -1;
}


void relion_timer::hip_cpu_tic(std::string id)
{
	if (hip_benchmark_find_id(id, hip_cpu_benchmark_identifiers) == -1)
	{
		hip_cpu_benchmark_identifiers.push_back(id);
		hip_cpu_benchmark_start_times.push_back(clock());
	}
	else
	{
		printf("DEBUG_ERROR: Provided identifier '%s' already exists in call to hip_cpu_tic.\n", id.c_str());
		CRITICAL(ERRCTIC);
	}
}

void relion_timer::hip_cpu_toc(std::string id)
{
	int idx = hip_benchmark_find_id(id, hip_cpu_benchmark_identifiers);
	if (idx == -1)
	{
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to hip_cpu_toc.\n", id.c_str());
		//exit( EXIT_FAILURE );
	}
	else
	{
		clock_t t = clock() - hip_cpu_benchmark_start_times[idx];
		hip_cpu_benchmark_identifiers.erase(hip_cpu_benchmark_identifiers.begin()+idx);
		hip_cpu_benchmark_start_times.erase(hip_cpu_benchmark_start_times.begin()+idx);
		fprintf(hip_cpu_benchmark_fPtr,"%06.2f ms ......", (float)t / CLOCKS_PER_SEC * 1000.);
		for (int i = 1; i < hip_cpu_benchmark_identifiers.size(); i++)
			fprintf(hip_cpu_benchmark_fPtr,"......");
		fprintf(hip_cpu_benchmark_fPtr," %s\n", id.c_str());
//		printf(,"%s \t %.2f ms\n", id.c_str(), (float)t / CLOCKS_PER_SEC * 1000.);
	}
}

void relion_timer::hip_gpu_tic(std::string id)
{
	if (hip_benchmark_find_id(id, hip_gpu_benchmark_identifiers) == -1)
	{
		hipEvent_t start, stop;
		hipEventCreate(&start);
		hipEventCreate(&stop);
		hipEventRecord(start, 0);
		hip_gpu_benchmark_identifiers.push_back(id);
		hip_gpu_benchmark_start_times.push_back(start);
		hip_gpu_benchmark_stop_times.push_back(stop);
	}
	else
	{
		printf("DEBUG_ERROR: Provided identifier '%s' already exists in call to hip_gpu_tic.\n",
				id.c_str());
		CRITICAL(ERRGTIC);
	}
}

void relion_timer::hip_gpu_toc(std::string id)
{
	int idx = hip_benchmark_find_id(id, hip_gpu_benchmark_identifiers);
	if (idx == -1)
	{
		printf("DEBUG_ERROR: Provided identifier '%s' not found in call to hip_gpu_tac.\n",
				id.c_str());
		CRITICAL(ERRGTOC);
	}
	else
	{
		hipEventRecord(hip_gpu_benchmark_stop_times[idx], 0);
		hipEventSynchronize(hip_gpu_benchmark_stop_times[idx]);
	}
}

void relion_timer::hip_gpu_printtictoc()
{
	if (hip_gpu_benchmark_identifiers.size() == 0)
	{
		printf("DEBUG_ERROR: There were no identifiers found in the list, on call to hip_gpu_toc.\n");
		CRITICAL(ERRTPC);
	}
	else
	{
		float time;
		for (int idx = 0; idx < hip_gpu_benchmark_identifiers.size(); idx ++)
		{
			hipEventElapsedTime(&time, hip_gpu_benchmark_start_times[idx],
					hip_gpu_benchmark_stop_times[idx]);
			hipEventDestroy(hip_gpu_benchmark_start_times[idx]);
			hipEventDestroy(hip_gpu_benchmark_stop_times[idx]);
			fprintf(hip_gpu_benchmark_fPtr,"%.2f ms \t %s\n",
					time, hip_gpu_benchmark_identifiers[idx].c_str());
		}

		hip_gpu_benchmark_identifiers.clear();
		hip_gpu_benchmark_start_times.clear();
		hip_gpu_benchmark_stop_times.clear();
	}
}

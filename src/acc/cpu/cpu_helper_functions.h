#ifndef CPU_HELPER_FUNCTIONS_H_
#define CPU_HELPER_FUNCTIONS_H_

#include "src/acc/cpu/cpu_ml_optimiser.h"
#include "src/acc/acc_projector.h"
#include "src/acc/cpu/cpu_benchmark_utils.h"
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <signal.h>
#include "src/complex.h"
#include "src/parallel.h"
#include <tbb/parallel_for.h>
#include <tbb/queuing_mutex.h>

namespace CpuKernels
{
#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
void window_fourier_transform(
		size_t grid_dim, 
		size_t Npsi, 
		size_t block_size,
		ACCCOMPLEX *g_in,
		ACCCOMPLEX *g_out,
		size_t iX, size_t iY, size_t iZ, size_t iYX, //Input dimensions
		size_t oX, size_t oY, size_t oZ, size_t oYX, //Output dimensions
		size_t max_idx,
		size_t max_r2 = 0
		)
{
	for(size_t blk=0; blk<grid_dim; blk++) {
		for(size_t psi=0; psi<Npsi; psi++) {
			for(size_t tid=0; tid< block_size; tid++) {
				size_t n = tid + block_size * blk;
				size_t oOFF = oX*oY*oZ*psi;
				size_t iOFF = iX*iY*iZ*psi;
				if (n >= max_idx) return;

				long int k, i, kp, ip, jp;

				if (check_max_r2)
				{
					k = n / (iX * iY);
					i = (n % (iX * iY)) / iX;

					kp = k < iX ? k : k - iZ;
					ip = i < iX ? i : i - iY;
					jp = n % iX;

					if (kp*kp + ip*ip + jp*jp > max_r2)
						return;
				}
				else
				{
					k = n / (oX * oY);
					i = (n % (oX * oY)) / oX;

					kp = k < oX ? k : k - oZ;
					ip = i < oX ? i : i - oY;
					jp = n % oX;
				}

				long int  in_idx = (kp < 0 ? kp + iZ : kp) * iYX + (ip < 0 ? ip + iY : ip)*iX + jp;
				long int out_idx = (kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip)*oX + jp;
				g_out[out_idx + oOFF] =  g_in[in_idx + iOFF];
			} // for tid
		} // for psi
	} // for blk
}

// may need to change to parallel reduce if it becomes the bottle neck.
template <typename T>
static T getMin(T *data, size_t size)
{
	T min = data[0];
	for(size_t i=1; i<size; i++)
		min = data[i] < min ? data[i] : min;
	 
	return min;
}

// may need to change to parallel reduce if it becomes the bottle neck.
template <typename T>
static T getMax(T *data, size_t size)
{
	T max = data[0];
	for(size_t i=1; i<size; i++)
		max = data[i] > max ? data[i] : max;

	return max;
}

template <typename T>
static T getSum(T *data, size_t size)
{
	T sum = data[0];
	for(size_t i=1; i<size; i++)
		sum += data[i];
	 
	return sum;
}

template <typename T>
static std::pair<size_t, T> getArgMin(T *data, size_t size)
{
	std::pair<size_t, T> pair;
	pair.first = 0;
	pair.second = data[0];
	
	for(size_t i=1; i<size; i++)
		if( data[i] < pair.second) {
			pair.first = i;
			pair.second = data[i];
		}
		
	return pair;
}

template <typename T>
static std::pair<size_t, T> getArgMax(T *data, size_t size)
{
	std::pair<size_t, T> pair;
	pair.first = 0;
	pair.second = data[0];
	
	for(size_t i=1; i<size; i++)
		if( data[i] > pair.second) {
			pair.first = i;
			pair.second = data[i];
		}
		
	return pair;
}

} // Namespace CpuKernels

#endif //CPU_HELPER_FUNCTIONS_H_


#ifndef SYCL_HELPER_FUNCTIONS_H_
#define SYCL_HELPER_FUNCTIONS_H_

#include <sys/time.h>
#include <signal.h>
#include <cstdio>
#include <ctime>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <limits>
#include <algorithm>
#include <numeric>

#include "src/complex.h"
#include "src/parallel.h"
#include "src/acc/sycl/sycl_ml_optimiser.h"
#include "src/acc/sycl/sycl_benchmark_utils.h"
#include "src/acc/acc_projector.h"
#include "src/acc/sycl/sycl_kernels/sycl_utils.h"

namespace syclKernels
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
				g_out[out_idx + oOFF] = g_in[in_idx + iOFF];
			} // for tid
		} // for psi
	} // for blk
}

// may need to change to parallel reduce if it becomes the bottle neck.
template <typename T>
static T getMin(T *data, size_t size)
{
//	return *std::min_element(std::execution::unseq, data, data+size);
	T minv = std::numeric_limits<T>::max();
#if _OPENMP >= 201307	// For OpenMP 4.0 and later
	#pragma omp simd reduction(min:minv)
#endif
	for (size_t i = 0; i < size; i++)
		if (data[i] < minv)	minv = data[i];
	return minv;
}

// may need to change to parallel reduce if it becomes the bottle neck.
template <typename T>
static T getMax(T *data, size_t size)
{
//	return *std::max_element(std::execution::unseq, data, data+size);
	T maxv = std::numeric_limits<T>::lowest();
#if _OPENMP >= 201307	// For OpenMP 4.0 and later
	#pragma omp simd reduction(max:maxv)
#endif
	for (size_t i = 0; i < size; i++)
		if (data[i] > maxv)	maxv = data[i];
	return maxv;
}

template <typename T>
static T getSum(T *data, size_t size)
{
//	return std::reduce(std::execution::unseq, data, data+size);
	T sum = static_cast<T>(0);
#if _OPENMP >= 201307	// For OpenMP 4.0 and later
	#pragma omp simd reduction(+:sum)
#endif
	for (size_t i = 0; i < size; i++)
		sum += data[i];
	return sum;
}

template <typename T>
inline void min_loc(std::pair<size_t, T> *out, std::pair<size_t, T> *in)
{
	if (out->second > in->second)
	{
		out->first = in->first;
		out->second = in->second;
	}
}

template <typename T>
inline void max_loc(std::pair<size_t, T> *out, std::pair<size_t, T> *in)
{
	if (out->second < in->second)
	{
		out->first = in->first;
		out->second = in->second;
	}
}

template <typename T>
static std::pair<size_t, T> getArgMin(T *data, size_t size)
{
#if 0
	auto dist = std::min_element(std::execution::unseq, data, data+size);

	std::pair<size_t, T> pair;
	pair.first = std::distance(data, dist);
	pair.second = data[pair.first];
	return pair;
#else
	std::pair<size_t, T> pair {-1, std::numeric_limits<T>::max()};
 #if _OPENMP >= 201307	// For OpenMP 4.0 and later
	#pragma omp declare reduction(minloc: std::pair<size_t, T>: min_loc<T>(&omp_out, &omp_in)) \
		initializer(omp_priv = {-1, std::numeric_limits<T>::max()})
	#pragma omp simd reduction(minloc:pair)
 #endif
	for(size_t i=0; i<size; i++)
		if( data[i] < pair.second)
		{
			pair.first = i;
			pair.second = data[i];
		}

	return pair;
#endif
}

template <typename T>
static std::pair<size_t, T> getArgMax(T *data, size_t size)
{
#if 0
	auto dist = std::max_element(std::execution::unseq, data, data+size);

	std::pair<size_t, T> pair;
	pair.first = std::distance(data, dist);
	pair.second = data[pair.first];
	return pair;
#else
	std::pair<size_t, T> pair {-1, std::numeric_limits<T>::lowest()};
 #if _OPENMP >= 201307	// For OpenMP 4.0 and later
	#pragma omp declare reduction(maxloc: std::pair<size_t, T>: max_loc<T>(&omp_out, &omp_in)) \
		initializer(omp_priv = {-1, std::numeric_limits<T>::lowest()})
	#pragma omp simd reduction(maxloc:pair)
 #endif
	for(size_t i=0; i<size; i++)
		if( data[i] > pair.second)
		{
			pair.first = i;
			pair.second = data[i];
		}

	return pair;
#endif
}

} // Namespace syclKernels

#endif //CPU_HELPER_FUNCTIONS_H_


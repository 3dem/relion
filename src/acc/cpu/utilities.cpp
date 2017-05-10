#include "src/gpu_utils/cuda_settings.h"


namespace CpuKernels
{

template <typename T>
void multiply(
	T *A,
	T *OUT,
	T S,
	int image_size)
{
	for (int i = 0; i < image_size; i ++)
		OUT[i] = A[i]*S;
}

template <typename T>
void multiply(
	T *A,
	T S,
	int image_size)
{
	for (int i = 0; i < image_size; i ++)
		A[i] *= S;
}

template <typename T>
void multiply(
	T *A,
	T *B,
	T *OUT,
	T S,
	int image_size)
{
	for (int i = 0; i < image_size; i ++)
		OUT[i] = A[i]*B[i]*S;
}


} //namespace CpuKernels


#include "src/gpu_utils/cuda_settings.h"


namespace CudaKernels
{
template <typename T>
__global__ void multiply(
	T *A,
	T *OUT,
	T S,
	int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel] = A[pixel]*S;
}
template <typename T>
__global__ void multiply(
	T *A,
	T S,
	int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		A[pixel] = A[pixel]*S;
}

template <typename T>
__global__ void multiply(
	T *A,
	T *B,
	T *OUT,
	T S,
	int image_size)
{
	int pixel = threadIdx.x + blockIdx.x*BLOCK_SIZE;
	if(pixel<image_size)
		OUT[pixel] = A[pixel]*B[pixel]*S;
}

} //namespace CudaKernels


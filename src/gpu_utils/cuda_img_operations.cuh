#ifndef CUDA_IMG_OPERATIONS_CUH_
#define CUDA_IMG_OPERATIONS_CUH_
#include <cuda.h>

__global__ void cuda_kernel_applyAB(double *img, double *myAB, double *shifted_img, int size);

__global__ void cuda_kernel_diff2(double *ref, double *img, double *Minvsigma2, double *partial_sums, int size);

#endif

#include "src/gpu_utils/cuda_img_operations.h"
#include "src/gpu_utils/cuda_img_operations.cuh"
#include <vector>
#include <stdio.h>

#define cuda_block_size 32

__global__ void cuda_kernel_applyAB(double *img, double *myAB, double *shifted_img, int size)
{
	int n = (blockIdx.x * blockDim.x + threadIdx.x)*2;
	if(n < 2*size)
	{
		*(shifted_img + n) = (*(myAB + n)) * (*(img + n))
				- (*(myAB + n + 1)) * (*(img + n + 1));
		*(shifted_img + n + 1) = (*(myAB + n)) * (*(img + n + 1))
				+ (*(myAB + n + 1)) * (*(img + n));
	}
}

__global__ void cuda_kernel_diff2(double *ref, double *img, double *Minvsigma2, double *partial_sums, int size)
{
    int n = (blockIdx.x * blockDim.x + threadIdx.x)*2;
   __shared__ double s[cuda_block_size];

    double diff_real = (*(ref + n)) - (*(img + n));
	double diff_imag = (*(ref + n + 1)) - (*(img + n + 1));

	if(n < 2*size)
		s[threadIdx.x] = (diff_real * diff_real + diff_imag * diff_imag) * 0.5 * (*(Minvsigma2 + n/2));
	else
		s[threadIdx.x] = 0;

	__syncthreads();

	if (threadIdx.x == 0)
	{
		double sum = 0;
		for (int i = 0; i < cuda_block_size; i++)
			sum += s[i];
		partial_sums[blockIdx.x] = sum;
	}
}

void cuda_applyAB(
		int size,
		double *h_exp_local_Fimgs_shifted,
		double *h_myAB,
		double *h_Fimg_otfshift)
{
	int num_blocks=(ceil((float)size/(float)cuda_block_size));

	// declare GPU memory pointers
	double *d_myAB;
	double *d_exp_local_Fimgs_shifted;
	double *d_Fimg_otfshift;

	int N = size * sizeof(double) * 2; // x2 for real and imaginary part
	int N_pad = num_blocks * cuda_block_size * sizeof(double) * 2; // padded number

	//Allocate padding
	cudaMalloc( (void**) &d_myAB, N);
	cudaMalloc( (void**) &d_exp_local_Fimgs_shifted, N);
	cudaMalloc( (void**) &d_Fimg_otfshift, N);

	//Skip padding
	cudaMemcpy( d_myAB, h_myAB, N, cudaMemcpyHostToDevice);
	cudaMemcpy( d_exp_local_Fimgs_shifted, h_exp_local_Fimgs_shifted, N, cudaMemcpyHostToDevice);

	//int Nrest=N_pad-N;
	//for(int i=0; i<Nrest; i++)
	//{
	//	*(d_exp_local_Fimgs_shifted+N+i)=0.;
	//	*(d_myAB+N+i)=1.;
	//}
	cuda_kernel_applyAB<<<num_blocks, cuda_block_size>>>(d_exp_local_Fimgs_shifted, d_myAB, d_Fimg_otfshift, size);

	//Skip padding
	cudaMemcpy( h_Fimg_otfshift, d_Fimg_otfshift, N, cudaMemcpyDeviceToHost );
}

// returns device-allocated scalar diff2 from device-allocated image-triple
double cuda_diff2_deviceimage(
		int size,
		double *d_Frefctf,
		double *d_Fimg_trans,
		double *d_Minvsigma2)
{
	int num_blocks(ceil((float)size/(float)cuda_block_size));

	double sum;
	cudaMalloc( (void**) &sum, sizeof(double));
	double *d_partial_sums;
	int N_par = num_blocks * sizeof(double); // partial sums
	cudaMalloc( (void**) &d_partial_sums, N_par);

	cuda_kernel_diff2<<<num_blocks,cuda_block_size>>>(d_Frefctf, d_Fimg_trans, d_Minvsigma2, d_partial_sums, size);

	for (int i = 0; i < num_blocks; i ++)
		sum += d_partial_sums[i];

	return sum;
}

// returns host-allocated scalar diff2 from host-allocated image-triple
double cuda_diff2_hostimage(
		int size,
		double *h_Frefctf,
		double *h_Fimg_trans,
		double *h_Minvsigma2)
{
	int num_blocks=(ceil((float)size/(float)cuda_block_size));

	double sum(0);
	double *h_partial_sums = new double[num_blocks];

	// declare GPU memory pointers
	double *d_Frefctf;
	double *d_Fimg_trans;
	double *d_Minvsigma2;
	double *d_partial_sums;

	int N = size  * sizeof(double) * 2; // x2 for real and imaginary part
	int N_pad = num_blocks * cuda_block_size * sizeof(double) * 2; // padded number
	int N_par = num_blocks * sizeof(double); // partial sums

	//Allocate padding
	cudaMalloc( (void**) &d_Frefctf, N);
	cudaMalloc( (void**) &d_Fimg_trans, N);
	cudaMalloc( (void**) &d_Minvsigma2, N);
	cudaMalloc( (void**) &d_partial_sums, N_par);

	//Skip padding
	cudaMemcpy( d_Frefctf, h_Frefctf, N, cudaMemcpyHostToDevice);
	cudaMemcpy( d_Fimg_trans, h_Fimg_trans, N, cudaMemcpyHostToDevice);
	cudaMemcpy( d_Minvsigma2, h_Minvsigma2, N, cudaMemcpyHostToDevice);

	cuda_kernel_diff2<<<num_blocks, cuda_block_size>>>(d_Frefctf, d_Fimg_trans, d_Minvsigma2, d_partial_sums, size);
	//Skip padding
	cudaMemcpy( h_partial_sums, d_partial_sums, N_par, cudaMemcpyDeviceToHost);

	for (int i = 0; i < num_blocks; i ++)
		sum += h_partial_sums[i];

	return sum;
}

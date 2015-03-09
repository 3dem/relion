#include "src/gpu_utils/cuda_img_operations.h"
#include "src/gpu_utils/cuda_img_operations.cuh"
#include <vector>
#include <iostream>

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
	//int N_pad = num_blocks * cuda_block_size * sizeof(double) * 2; // padded number

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
		{
			sum += s[i];
		}
		partial_sums[blockIdx.x] = sum;
	}

}

__global__ void cuda_kernel_reduceArray(double * array, int size, double * sum)
{
	int n =	(blockIdx.x * blockDim.x + threadIdx.x);
	int i_max = (int)ceil((float)size/(float)blockDim.x);

	// if there are a number of thread-block-sized parts of the array to go through,
	// sum those parts into the first part of array
	if(i_max>1)
	{
		for(int i=1; i<i_max; i++)
		{
			if((i*blockDim.x+n)<size)
				array[n] += array[i*blockDim.x+n];
			__syncthreads();
		}
	}
	// then sum the work done by all threads
	for(int j=(blockDim.x/2); j>0; j/=2)
	{
		if (n<j and (n+j)<size) { array[n] += array[n+j]; }
			__syncthreads();
	}
	// finally, set the output sum
	if (n == 0 ) { *sum = array[0]; }

}


//       PLAYGROUND
__global__ void cuda_kernel_reduceArray_1(double * Garray, int size, double * sum)
{
	extern __shared__ float Sarray[];

	int tid = threadIdx.x;
	int n =	(blockIdx.x * blockDim.x + threadIdx.x);
	//Since separate blockIdx.x:s have separate shared memory srrays, we can index Sarray by tid only.
	if (n<size)
		Sarray[tid] = Garray[n];
	else
		Sarray[tid] = 0;
	__syncthreads();

//	for(unsigned int s=1; s<blockDim.x; s*=2)
//	{
//		if(tid % (2*s) == 0)
//			Sarray[tid] += Sarray[tid+s];
//		__syncthreads();
//	}
//
//	for(int s=blockDim.x; s>0; s-=1)
//	{
//		if (tid==0 && blockDim.x==s)
//		{
//			*sum += Sarray[tid]+5;
//			__syncthreads();
//		}
//	}
	if (n == 0 )
		*sum = Sarray[0];
}

// returns device-allocated scalar diff2 from device-allocated image-triple
void cuda_diff2_deviceImage(
		int size,
		double *d_Frefctf,
		double *d_Fimg_trans,
		double *d_Minvsigma2,
		double *d_diff2
		)
{
	int num_blocks(ceil((float)size/(float)cuda_block_size));

	double *d_partial_sums;
	int N_par = num_blocks * sizeof(double); // partial sums
	cudaMalloc( (void**) &d_partial_sums, N_par);

	cuda_kernel_diff2<<<num_blocks,cuda_block_size>>>(d_Frefctf, d_Fimg_trans, d_Minvsigma2, d_partial_sums, size);
	cuda_kernel_reduceArray<<<1,cuda_block_size>>>(d_partial_sums, num_blocks, d_diff2);

	cudaFree(d_partial_sums);
}

// returns host-allocated scalar diff2 from host-allocated image-triple
double cuda_diff2_hostImage(
		int size,
		double *h_Frefctf,
		double *h_Fimg_trans,
		double *h_Minvsigma2)
{
	double *h_diff2;
	h_diff2 = (double*) malloc(1);

	// declare GPU memory pointers
	double *d_Frefctf;
	double *d_Fimg_trans;
	double *d_Minvsigma2;
	double *d_diff2;

	//Currently we allocate without padding and instead opt for
	//kernel-based thread control for domain relevance.
	int N = size  * sizeof(double) * 2; // x2 for real and imaginary part

	cudaMalloc( (void**) &d_Frefctf, N);
	cudaMalloc( (void**) &d_Fimg_trans, N);
	cudaMalloc( (void**) &d_Minvsigma2, N);
	cudaMalloc( (void**) &d_diff2, sizeof(double));

	cudaMemcpy( d_Frefctf, h_Frefctf, N, cudaMemcpyHostToDevice);
	cudaMemcpy( d_Fimg_trans, h_Fimg_trans, N, cudaMemcpyHostToDevice);
	cudaMemcpy( d_Minvsigma2, h_Minvsigma2, N, cudaMemcpyHostToDevice);

	cuda_diff2_deviceImage(size, d_Frefctf, d_Fimg_trans, d_Minvsigma2, d_diff2);

	cudaMemcpy( h_diff2, d_diff2, sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(d_diff2);
	return h_diff2[0];
}


// returns a device-allocated pointer to a NxM diff2-matrix, after
// comparing two device-allocated sets of images of sizes N and M respectively.
// THIS DOES NOT WORK YET!!!!!!
void cuda_diff2_deviceImageSets(
		int N, double* d_setN,
		int M, double* d_setM,
		int imageSize,
		double* d_Minvsigma2
		)
{
	int pixelSize = sizeof(double)*2;

	// Allocate device memory for output array
	double** d_diff2Matrix = new double*[M];
	cudaMalloc( (void**) &d_diff2Matrix, M*sizeof(double*));
	for(int i = 0; i < M; ++i)
	 cudaMalloc( (void**) &(d_diff2Matrix[i]), N*sizeof(double));

	// Image pointers
	double* imageM;
	double* imageN;

	for(int j=0; j<M; j++)
	{
		imageM = (d_setM + j*imageSize*pixelSize);
		for(int k=0; k<N; k++)
		{
			imageN = (d_setM + j*imageSize*pixelSize);
			cuda_diff2_deviceImage(imageSize, (double*) imageM, (double*) imageN, (double*) d_Minvsigma2, d_diff2Matrix[j]);
		}
	}
}

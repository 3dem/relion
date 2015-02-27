#include "src/gpu_utils/diff2.h"
#include <vector>

#define cuda_block_size 32

__global__ void cuda_kernel_applyAB(double *img, double *myAB, double *shifted_img)
{
    int n = (blockIdx.x * blockDim.x + threadIdx.x)*2;
    *(shifted_img + n) = (*(myAB + n)) * (*(img + n))
    		- (*(myAB + n + 1)) * (*(img + n + 1));
    *(shifted_img + n + 1) = (*(myAB + n)) * (*(img + n + 1))
			+ (*(myAB + n + 1)) * (*(img + n));
}

//__global__ void cuda_kernel_diff2(Complex *ref, Complex* img, Complex* Minvsigma2, double* diff2)
//{
//    int n = threadIdx.x;
//    double diff_real = (*(ref + n)).real - (*(img + n)).real;
//	double diff_imag = (*(ref + n)).imag - (*(img + n)).imag;
//	// diff2 increment add needs to be atomic
//	diff2 += (diff_real * diff_real + diff_imag * diff_imag) * 0.5 * (*(Minvsigma2 + n));
//}

void cuda_applyAB(
		int img_num,
		double *h_exp_local_Fimgs_shifted,
		double *h_myAB,
		double *h_Fimg_otfshift)
{
	int num_blocks(ceil(img_num/cuda_block_size));

	// declare GPU memory pointers
	double * d_myAB;
	double * d_exp_local_Fimgs_shifted;
	double * d_Fimg_otfshift;

	int N = img_num  * sizeof(double) * 2; // x2 for real and imaginary part

	cudaMalloc( (void**) &d_myAB, N);
	cudaMalloc( (void**) &d_exp_local_Fimgs_shifted, N);
	cudaMalloc( (void**) &d_Fimg_otfshift, N);

	cudaMemcpy( d_myAB, h_myAB, N, cudaMemcpyHostToDevice);
	cudaMemcpy( d_exp_local_Fimgs_shifted, h_exp_local_Fimgs_shifted, N, cudaMemcpyHostToDevice);

	//let's do a simple setup for now; each pixel is a thread, each row is a block
	cuda_kernel_applyAB<<<num_blocks, cuda_block_size>>>(d_exp_local_Fimgs_shifted, d_myAB, d_Fimg_otfshift);

	cudaMemcpy( h_Fimg_otfshift, d_Fimg_otfshift, N, cudaMemcpyDeviceToHost );
}
//
//void cuda_diff2(
//		long int ipart,
//		std::vector<MultidimArray<Complex > > &h_exp_local_Fimgs_shifted,
//		Complex *h_myAB,
//		MultidimArray<Complex > h_Fimg_otfshift)
//{
//	// Size of the image arrays
//	int img_num = NZYXSIZE(h_exp_local_Fimgs_shifted[ipart]);
//	int num_blocks(ceil(img_num/cuda_block_size));
//
//	// declare GPU memory pointers
//	Complex * d_myAB;
//	Complex * d_exp_local_Fimgs_shifted;
//	Complex * d_Fimg_otfshift;
//
//	int N = img_num  * sizeof(Complex);
//
//	cudaMalloc( (void**) &d_myAB, N);
//	cudaMalloc( (void**) &d_exp_local_Fimgs_shifted, N);
//	cudaMalloc( (void**) &d_Fimg_otfshift, N);
//
//	cudaMemcpy( d_myAB, h_myAB, N, cudaMemcpyHostToDevice);
//	cudaMemcpy( d_exp_local_Fimgs_shifted, h_exp_local_Fimgs_shifted[ipart], N, cudaMemcpyHostToDevice);
//
//	//let's do a simple setup for now; each pixel is a thread, each row is a block
//	cuda_kernel_diff2<<<num_blocks, cuda_block_size>>>(d_exp_local_Fimgs_shifted, d_myAB, d_Fimg_otfshift);
//
//	cudaMemcpy( h_Fimg_otfshift, d_Fimg_otfshift, N, cudaMemcpyDeviceToHost );
//}

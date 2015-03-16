#ifndef CUDA_IMG_OPERATIONS_H_
#define CUDA_IMG_OPERATIONS_H_
#include <cuda.h>

void cuda_applyAB(
					int size,
					double *h_exp_local_Fimgs_shifted,
					double *myAB,
					double *shifted_img
				);
				  
				  
double cuda_diff2_hostImage(
					int size,
					double *h_Frefctf,
					double *h_Fimg_trans,
					double *h_Minvsigma2
				);

void cuda_diff2_deviceImage(
					int size,
					double *d_Frefctf,
					double *d_Fimg_trans,
					double *d_Minvsigma2,
					double *d_diff2
				);

#endif

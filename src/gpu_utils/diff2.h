#ifndef CUDA_DIFF2_H_
#define CUDA_DIFF2_H_
#include <cuda.h>

void cuda_applyAB(
					int size,
					double *h_exp_local_Fimgs_shifted,
					double *myAB,
					double *shifted_img
				);
				  
				  
double cuda_diff2(
					int size,
					double *h_Frefctf,
					double *h_Fimg_trans,
					double *h_Minvsigma2
				);


#endif

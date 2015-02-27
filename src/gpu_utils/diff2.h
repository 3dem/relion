#ifndef CUDA_DIFF2_H_
#define CUDA_DIFF2_H_
#include <cuda.h>

void cuda_applyAB(
					int size,
					double *h_exp_local_Fimgs_shifted,
					double *myAB,
					double *shifted_img
				  );
				  
				  
//void cuda_diff2(
//					Complex *ref,
//					Complex* img,
//					Complex* Minvsigma2,
//					double* diff2
//			   );
//
//
#endif

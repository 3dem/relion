
void cuda_kernel_applyAB();

void cuda_kernel_diff2(Complex *ref, Complex* img, Complex* Minvsigma2, double* diff2);

void cuda_applyAB(
					int ipart, 
					std::vector<MultidimArray<Complex > >, 
					Complex *myAB, 
					MultidimArray<Complex > shifted_img
				  );
				  
				  
void cuda_diff2(
					Complex *ref, 
					Complex* img, 
					Complex* Minvsigma2, 
					double* diff2
			   );



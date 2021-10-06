#include <iostream>


#ifdef CUDA
#include <src/jaz/cuda/test00.h>
#endif


int main(int argc, char *argv[])
{
	#ifdef CUDA

		CudaTest00 test;
		test.run();

	#else

		std::cerr << "Not compiled with CUDA." << std::endl;

	#endif

	return 0;
}

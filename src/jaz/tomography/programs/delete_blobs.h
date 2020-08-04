#ifndef DELETE_BLOBS_PROGRAM_H
#define DELETE_BLOBS_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>


class DeleteBlobsProgram
{
	public:
		
		DeleteBlobsProgram(){}
		
			std::string outPath, tomoSetFn, tomoName, spheresFn;
			
			bool diag;
			int SH_bands, num_threads, max_iters;
			double inner_margin, outer_margin;
			
		void readParameters(int argc, char *argv[]);
		void run();
};

#endif

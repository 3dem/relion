#ifndef BACKPROJECT_PROGRAM_H
#define BACKPROJECT_PROGRAM_H

#include <string>

class BackprojectProgram
{
	public:
		
		BackprojectProgram(){}
		
			std::string outTag, tomoSetFn, catFn, motFn, symmName;
			
			bool do_whiten, diag, no_subpix_off, explicit_gridding, no_reconstruction;
			int boxSize, cropSize, num_threads, outer_threads, inner_threads, max_mem_GB;
			double SNR, taper, binning;
			
		void readParameters(int argc, char *argv[]);
		void run();
};

#endif

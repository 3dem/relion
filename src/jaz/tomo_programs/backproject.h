#ifndef BACKPROJECT_PROGRAM_H
#define BACKPROJECT_PROGRAM_H

#include <string>

class BackprojectProgram
{
	public:
		
		BackprojectProgram(){}
		
			std::string outTag, tomoSetFn, catFn, motFn, symmName;
			
			bool diag, no_subpix_off, no_reconstruction;
			int boxSize, cropSize, num_threads, outer_threads, inner_threads, max_mem_GB;
			double WienerFract, taper, binning;
			
		void run();
};

#endif

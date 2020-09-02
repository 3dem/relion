#ifndef SUBSTACK_PROGRAM_H
#define SUBSTACK_PROGRAM_H

#include <string>
#include <vector>

class SubstackProgram
{
	public:
		
		SubstackProgram(){}
		
		
			std::string outTag, tomoListFn, catFn;
			
			int boxSize, num_threads;
			double binning;
			bool diag, do_center;
			
			
		void run();		
};

#endif

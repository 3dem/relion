#ifndef DARK_ERASE_PROGRAM_H
#define DARK_ERASE_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>

class DarkEraseProgram
{
	public:
		
		DarkEraseProgram(){}
		
			int num_threads;
			double thresh, rad;
			
			std::string stackFn, projFn, outFn;
			bool writeNormalized, diag;
			
		void run();		
};

#endif
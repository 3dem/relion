#ifndef POWSPEC_PROGRAM_H
#define POWSPEC_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>

class PowspecProgram
{
	public:
		
		PowspecProgram(){}
		
			int res, ratio;
			double pixSize;	
			std::string stackFn, outFn;
			bool separate;
			
		void run();		
};

#endif

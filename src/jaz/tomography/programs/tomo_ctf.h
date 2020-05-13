#ifndef TOMO_CTF_PROGRAM_H
#define TOMO_CTF_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>

class TomoCtfProgram
{
	public:
		
		TomoCtfProgram(){}
		
			int n_threads;	
			int thickness;	
			double pixSize, voltage, Cs;
			
			std::string stackFn, projFn, outFn;
				
		void run();		
};

#endif

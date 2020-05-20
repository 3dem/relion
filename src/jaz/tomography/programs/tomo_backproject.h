#ifndef TOMO_BP_PROGRAM_H
#define TOMO_BP_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>

class TomoBackprojectProgram
{
	public:
		
		TomoBackprojectProgram(){}
			
			int n_threads;	
			int thickness, w, h, tomoIndex;	
			double spacing, stack_spacing, x0, y0, z0, taperRad;	
			std::string tomoSetFn, outFn;	
			bool weight, zeroDC;	
			double WienerFract;
			
		void run();		
};

#endif

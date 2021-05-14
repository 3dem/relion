#ifndef TOMO_BP_PROGRAM_H
#define TOMO_BP_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/tomography/optimisation_set.h>

class TomoBackprojectProgram
{
	public:
		
		TomoBackprojectProgram(){}
			
			int n_threads;
			int w, h, d;
			double spacing, x0, y0, z0, taperDist, taperFalloff;
			std::string tomoName, outFn;
			bool applyPreWeight, applyWeight, applyCtf, zeroDC, FourierCrop;
			double SNR;

			OptimisationSet optimisationSet;
			
			
		void readParameters(int argc, char *argv[]);
		void run();		
};

#endif

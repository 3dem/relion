#ifndef BFACTOR_FIT_PROGRAM_H
#define BFACTOR_FIT_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <vector>

#include "refinement.h"


class BfactorFitProgram : public RefinementProgram
{
	public:
		
		BfactorFitProgram(int argc, char *argv[]);
		
		int kMin, kMin2;
		bool useL2, useCache;
		
		void readParams();
		void run();
};


#endif

#ifndef ABERRATION_FIT_PROGRAM_H
#define ABERRATION_FIT_PROGRAM_H

#include "refinement.h"
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/optics/aberration_fit.h>
#include <src/jaz/math/tensor2x2.h>


class AberrationFitProgram : public RefinementProgram
{
	public:
			 
		
		AberrationFitProgram(int argc, char *argv[]);
		

			int n_even, n_odd;
			bool do_even, do_odd;
		
		void readParams(IOParser& parser);
		void run();
};

#endif

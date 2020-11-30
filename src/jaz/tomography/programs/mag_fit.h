#ifndef MAG_FIT_PROGRAM_H
#define MAG_FIT_PROGRAM_H

#include "refinement.h"
#include <src/jaz/optics/aberrations_cache.h>
#include <src/jaz/optics/aberration_fit.h>
#include <src/jaz/math/tensor2x2.h>


class MagFitProgram : public RefinementProgram
{
	public:


		MagFitProgram(int argc, char *argv[]);

		void readParams();
		void run();
};

#endif

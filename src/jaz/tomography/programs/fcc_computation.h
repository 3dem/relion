#ifndef FCC_PROGRAM_H
#define FCC_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <vector>

#include "refinement.h"


class FccProgram : public RefinementProgram
{
	public:
		
		FccProgram(int argc, char *argv[]);

		void readParams();
		void run();
};


#endif

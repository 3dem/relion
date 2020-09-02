#ifndef PROTO_POLISH_PROGRAM_H
#define PROTO_POLISH_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <vector>

#include "refinement.h"

class CTF;

class FrameAlignmentProgram : public RefinementProgram
{
	public:
		
		FrameAlignmentProgram(int argc, char *argv[]);
		
			bool shiftOnly, whiten, whiten_abs,  
				const_particles, const_angles, const_shifts;
			
			int range, num_iters;	
			
			double padding, hiPass_px, sig2RampPower;
		
		void readParams(IOParser& parser);
		void run();
};


#endif

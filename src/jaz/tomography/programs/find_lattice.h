#ifndef FIND_LATTICE_PROGRAM_H
#define FIND_LATTICE_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>

class FindLatticeProgram
{
	public:
		
		FindLatticeProgram(){}
		
			int n_threads, overtones;
			
			double spacing_ang, filter_width, angpix, taper, minValue, minDensity;
			bool noise;
			std::string tomoFn, outFn;
			
		void run();		
};

#endif

#ifndef SAMPLE_MANIFOLD_PROGRAM_H
#define SAMPLE_MANIFOLD_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/buffered_image.h>

class SampleManifoldProgram
{
	public:
		
		SampleManifoldProgram(){}
		
		void readParameters(int argc, char *argv[]);
		void run();	
		
};

#endif

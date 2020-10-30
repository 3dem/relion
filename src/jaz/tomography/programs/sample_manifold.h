#ifndef SAMPLE_MANIFOLD_PROGRAM_H
#define SAMPLE_MANIFOLD_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/optimisation_set.h>

class SampleManifoldProgram
{
	public:
		
		SampleManifoldProgram(){}

			std::string output_path;

			OptimisationSet optimisationSet;

			bool avoid_missing_wedge, avoid_present_wedge, store_tilt_series;
			double spacing, depth, max_tilt;
		
		void readParameters(int argc, char *argv[]);
		void run();	
		
};

#endif

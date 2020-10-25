#ifndef REFINEMENT_PROGRAM_H
#define REFINEMENT_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <vector>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/reference_map.h>
#include <src/jaz/optics/optics_data.h>


class CTF;


class RefinementProgram
{
	public:
		
		RefinementProgram(int argc, char *argv[]);
		
			int argc;
			char** argv;
			
			std::string outDir, tomoSetFn, particlesFn, motFn;
			
			bool diag, timing;
			
			int boxSize, num_threads, specified_first_frame, specified_last_frame;

			
			
			ParticleSet dataSet;
			std::vector<std::vector<ParticleIndex>> particles;
			TomogramSet tomogramSet;
			
			TomoReferenceMap referenceMap;
			
			
		void _readParams(IOParser& parser);
		
		void init();
		
		BufferedImage<float> computeFrequencyWeights(
						const Tomogram& tomogram,
						bool whiten, double sig2RampPower, double hiPass_px,
						int num_threads);
};


#endif

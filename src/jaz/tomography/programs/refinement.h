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
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/optics/optics_data.h>


class CTF;


class RefinementProgram
{
	public:
		
		RefinementProgram(int argc, char *argv[]);
		
			int argc;
			char** argv;
			
			std::string outDir;
			OptimisationSet optimisationSet;
			
			bool diag, static_noise, only_do_unfinished, run_from_GUI, run_from_MPI;
			
			int boxSize, num_threads, verbosity;

			
			
			ParticleSet particleSet;
			std::vector<std::vector<ParticleIndex>> particles;
			TomogramSet tomogramSet;
			
			TomoReferenceMap referenceMap;
			
			
		void _readParams(IOParser& parser);
		
		void init();
		
		BufferedImage<float> computeFrequencyWeights(
						const Tomogram& tomogram,
						bool whiten, double sig2RampPower, double hiPass_px,
						bool applyDoseWeight, int num_threads);

		BufferedImage<int> findXRanges(
				const RawImage<float>& freqWeights,
				const RawImage<float>& doseWeights,
				double cutoffFraction);
};


#endif

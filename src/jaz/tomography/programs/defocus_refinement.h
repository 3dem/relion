#ifndef DEFOCUS_REFINE_PROGRAM_H
#define DEFOCUS_REFINE_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/optics/aberration_fit.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <vector>

#include "refinement.h"

class ParticleSet;
class CTF;

class DefocusRefinementProgram : public RefinementProgram
{
	public:
		
		DefocusRefinementProgram(int argc, char *argv[]);
			
			bool clearAstigmatism, scanDefocus, slowScan, refineFast,
				refineAstigmatism, do_plotAstigmatism, regularise;

			int deltaSteps, max_particles, group_count;
			double minDelta, maxDelta, sigma_input;
			
		void run();
		
		
	private:
		
		struct DefocusFit
		{
			double value, stdDev;
			std::vector<double> offsets;
			std::vector<double> totalCost;
			std::vector<std::vector<double>> costByGroup;
		};
		
		static DefocusFit findDefocus(
				int f,  
				double minDelta, 
				double maxDelta,
				int steps, int group_count, double sigma_input,
				const ParticleSet& dataSet,
				std::vector<int>& particles, int max_particles,
				const Tomogram& tomogram,
				std::vector<BufferedImage<fComplex>>& referenceFS,
				const BufferedImage<float>& freqWeights,
				bool flip_value, 
				int num_threads);
		
		static BufferedImage<double> computeOffsetCost(
				int f,
				double z0, double z1, int steps, 
				const ParticleSet& dataSet,
				std::vector<int>& particles, int max_particles,
				const Tomogram& tomogram,
				std::vector<BufferedImage<fComplex>>& referenceFS,
				const BufferedImage<float>& freqWeights,
				bool flip_value, 
				double handedness,
				int num_threads);

		static double scanForDefocus(
				const BufferedImage<aberration::EvenData>& evenData,
				double pixelSize,
				const CTF& ctf0,
				double minDefocus,
				double maxDefocus,
				int steps);

		static gravis::d3Vector findAstigmatism(
				const aberration::EvenSolution& solution,
				const CTF& referenceCtf,
				double initialDeltaZ,
				double pixelSize,
				double initialStep);

		static BufferedImage<double> plotAstigmatism(
				const aberration::EvenSolution& solution,
				const CTF& referenceCtf,
				double initialDeltaZ,
				double range,
				double pixelSize,
				int size);
};


#endif

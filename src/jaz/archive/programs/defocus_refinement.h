#ifndef DEFOCUS_REFINE_PROGRAM_H
#define DEFOCUS_REFINE_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/optics/aberration_fit.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <vector>

#include "refinement.h"

class CTF;

class DefocusRefinementProgram : public RefinementProgram
{
	public:
		
		DefocusRefinementProgram(int argc, char *argv[]);
			
			bool do_clearAstigmatism, do_defocus, do_scanDefocus, do_slowScan, do_refineFast,
				do_refineAstigmatism, do_plotAstigmatism, do_regularise,
				do_slopeFit, do_perTomogramSlope;

			int deltaSteps, max_particles, group_count;
			double minDelta, maxDelta, sigma_input, max_slope_dose;
			
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
				const ParticleSet& particleSet,
				std::vector<ParticleIndex>& particles, int max_particles,
				const Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				std::vector<BufferedImage<fComplex>>& referenceFS,
				const BufferedImage<float>& freqWeights,
				bool flip_value, 
				int num_threads);
		
		static BufferedImage<double> computeOffsetCost(
				int f,
				double z0, double z1, int steps, 
				const ParticleSet& particleSet,
				std::vector<ParticleIndex>& particles, int max_particles,
				const Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				std::vector<BufferedImage<fComplex>>& referenceFS,
				const BufferedImage<float>& freqWeights,
				bool flip_value,
				int num_threads);

		std::vector<gravis::d3Vector> computeSlopeCost(
				double max_dose,
				double m0, double m1, int steps,
				const ParticleSet& particleSet,
				std::vector<ParticleIndex>& particles, int max_particles,
				const Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				std::vector<BufferedImage<fComplex>>& referenceFS,
				const BufferedImage<float>& freqWeights,
				bool flip_value,
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

		void writeSlopeCost(
				const std::vector<gravis::d3Vector>& cost,
				const std::string& filename);
};


#endif

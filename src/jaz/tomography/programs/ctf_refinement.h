#ifndef CTF_REFINE_PROGRAM_H
#define CTF_REFINE_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/optics/aberration_fit.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/optimization/optimization.h>
#include <vector>

#include "refinement.h"

class CTF;

class CtfRefinementProgram : public RefinementProgram
{
	public:
		
		CtfRefinementProgram(int argc, char *argv[]);
			
			bool do_defocus, do_refineFast;

			int deltaSteps;
			double minDelta, maxDelta, lambda_reg;
			
		void run();
		
		
	private:
		
		struct DefocusFit
		{
			double value, stdDev;
			std::vector<double> offsets;
			std::vector<double> totalCost;
			std::vector<std::vector<double>> costByGroup;
		};

		BufferedImage<double> evaluateDefocusRange(
				const BufferedImage<aberration::EvenData>& evenData,
				double pixelSize,
				const std::vector<CTF>& ctfs,
				double minDefocus,
				double maxDefocus,
				int steps);

		static gravis::d3Vector findAstigmatism(
				const aberration::EvenSolution& solution,
				const CTF& referenceCtf,
				double initialDeltaZ,
				double pixelSize,
				double initialStep);

		static std::vector<gravis::d3Vector> findMultiAstigmatism(
				const aberration::EvenSolution& solution,
				const std::vector<CTF>& referenceCtfs,
				double initialDeltaZ,
				double pixelSize,
				double lambda_reg);

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


class GlobalDefocusFit : public Optimization
{
	public:

		GlobalDefocusFit(
			const BufferedImage<double>& dataTerm,
			double defocusStep,
			double lambda);

			const BufferedImage<double>& dataTerm;
			double defocusStep, lambda;

		double f(const std::vector<double>& x, void* tempStorage) const;
};


#endif

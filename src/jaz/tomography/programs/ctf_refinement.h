#ifndef CTF_REFINE_PROGRAM_H
#define CTF_REFINE_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/optics/aberration_fit.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/optics/aberration_fit.h>
#include <vector>

#include "refinement.h"

class CTF;

class CtfRefinementProgram : public RefinementProgram
{
	public:
		
		CtfRefinementProgram(int argc, char *argv[]);
			
			bool do_refine_defocus, do_refine_scale, do_refine_aberrations,
				do_even_aberrations, do_odd_aberrations, do_fit_Beer_Lambert;

			int deltaSteps, n_even, n_odd;
			double minDelta, maxDelta, lambda_reg;
			
		void run();
		
		
	private:

		void processTomograms(
				int first_t,
				int last_t,
				const AberrationsCache& aberrationsCache,
				std::vector<BufferedImage<aberration::EvenData>>& evenData_perGroup,
				std::vector<std::vector<BufferedImage<aberration::EvenData>>>& evenData_perGroup_perThread,
				std::vector<BufferedImage<aberration::OddData>>& oddData_perGroup,
				std::vector<std::vector<BufferedImage<aberration::OddData>>>& oddData_perGroup_perThread,
				int verbosity);

		
		void refineDefocus(
				int t,
				Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				const BufferedImage<float>& freqWeights,
				const BufferedImage<float>& doseWeights);

		void fitScale(
				int t,
				Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				const BufferedImage<float>& freqWeights,
				const BufferedImage<float>& doseWeights);

		void updateAberrations(
				int t,
				const Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				const BufferedImage<float>& freqWeights,
				const BufferedImage<float>& doseWeights,
				std::vector<BufferedImage<aberration::EvenData>>& evenData_perGroup,
				std::vector<std::vector<BufferedImage<aberration::EvenData>>>& evenData_perGroup_perThread,
				std::vector<BufferedImage<aberration::OddData>>& oddData_perGroup,
				std::vector<std::vector<BufferedImage<aberration::OddData>>>& oddData_perGroup_perThread);

		void fitAberrations(
				std::vector<BufferedImage<aberration::EvenData>>& evenData_perGroup,
				std::vector<BufferedImage<aberration::OddData>>& oddData_perGroup,
				double pixelSize);




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
};

class BeerLambertFit : public Optimization
{
	public:

		BeerLambertFit(
			const std::vector<gravis::d4Matrix>& projections,
			const std::vector<double>& sum_prdObs,
			const std::vector<double>& sum_prdSqr);

			const std::vector<gravis::d4Matrix>& projections;
			const std::vector<double>& sum_prdObs;
			const std::vector<double>& sum_prdSqr;
			std::vector<gravis::d3Vector> view_dir;

		double f(const std::vector<double>& x, void* tempStorage) const;

		double getScale(int f, const std::vector<double>& x);

};

#endif

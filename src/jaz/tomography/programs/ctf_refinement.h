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
			
			bool do_refine_defocus,
				do_reset_to_common, do_regularise_defocus,
				do_refine_scale, do_refine_aberrations,
				do_fit_Lambert_per_tomo, do_fit_Lambert_globally,
				do_even_aberrations, do_odd_aberrations;

			int deltaSteps, n_even, n_odd, min_frame, max_frame;
			double minDelta, maxDelta, lambda_reg, k_min_Ang, freqCutoffFract;
			
		void run();
		
		
	protected:

		void parseInput();

		void initTempDirectories();

		void processTomograms(
				const std::vector<int>& tomoIndices,
				const AberrationsCache& aberrationsCache,
				bool per_tomogram_progress);

		void finalise();

	private:
		
		void refineDefocus(
				int t,
				Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				const BufferedImage<float>& freqWeights,
				const BufferedImage<float>& doseWeights,
				const BufferedImage<int>& xRanges,
				double k_min_px,
				int verbosity);

		void updateScale(
				int t,
				Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				const BufferedImage<float>& freqWeights,
				const BufferedImage<float>& doseWeights,
				int verbosity);

		void updateAberrations(
				int t,
				const Tomogram& tomogram,
				const AberrationsCache& aberrationsCache,
				const BufferedImage<float>& freqWeights,
				const BufferedImage<float>& doseWeights,
				const BufferedImage<int>& xRanges,
				int verbosity);


		void collectDefocus();

		void fitGlobalScale();
		void collectScale();

		void fitAberrations(int k_min_px);



		BufferedImage<double> evaluateDefocusRange(
				const BufferedImage<aberration::EvenData>& evenData,
				double pixelSize,
				const std::vector<CTF>& ctfs,
				double minDefocus,
				double maxDefocus,
				int steps,
				double k_min_px);

		static gravis::d3Vector findAstigmatism(
				const aberration::EvenSolution& solution,
				const CTF& referenceCtf,
				double initialDeltaZ,
				double pixelSize,
				double initialStep,
				double k_min_px);

		static std::vector<gravis::d3Vector> findMultiAstigmatism(
				const aberration::EvenSolution& solution,
				const std::vector<CTF>& referenceCtfs,
				double initialDeltaZ,
				double pixelSize,
				double lambda_reg,
				double k_min_px);


		std::string getDefocusTempFilenameRoot(
				const std::string& tomogram_name);

		std::string getScaleTempFilenameRoot(
				const std::string& tomogram_name);

		std::string getEvenAberrationsTempFilename(
				const std::string& tomogram_name, int opticsGroup);

		std::string getOddAberrationsTempFilename(
				const std::string& tomogram_name, int opticsGroup);


		bool defocusAlreadyDone(const std::string& tomogram_name);
		bool scaleAlreadyDone(const std::string& tomogram_name);
		bool aberrationsAlreadyDone(const std::string& tomogram_name, int group_count);


		void writeDefocusEps(const MetaDataTable& table, const std::string& tomo_name);
		void writeScaleEps(const MetaDataTable& table, const std::string& tomo_name);
		void mergeLogFiles();

		void abortIfNeeded();
};


class LambertFit : public Optimization
{
	public:

		LambertFit(
				const std::vector<gravis::d4Matrix>& projections,
				const std::vector<double>& sum_prdObs,
				const std::vector<double>& sum_prdSqr);

			const std::vector<gravis::d4Matrix>& projections;
			const std::vector<double> &sum_prdObs, &sum_prdSqr;
			std::vector<gravis::d3Vector> view_dir;
			gravis::d3Vector tilt_p, tilt_q;

		double f(const std::vector<double>& x, void* tempStorage) const;

		double gradAndValue(const std::vector<double>& x, std::vector<double>& gradDest) const;

		double getScale(int f, const std::vector<double>& x);

};

class MultiLambertFit : public FastDifferentiableOptimization
{
	public:

		MultiLambertFit(
				const std::vector<std::vector<gravis::d4Matrix>>& projections,
				const std::vector<std::vector<double>>& sum_prdObs,
				const std::vector<std::vector<double>>& sum_prdSqr,
				const std::vector<double>& fractional_dose);

			const std::vector<std::vector<gravis::d4Matrix>>& projections;
			const std::vector<std::vector<double>> &sum_prdObs, &sum_prdSqr;
			const std::vector<double>& fractional_dose;
			std::vector<std::vector<gravis::d3Vector>> view_dir;
			std::vector<gravis::d3Vector> tilt_p, tilt_q;


		double gradAndValue(const std::vector<double>& x, std::vector<double>& gradDest) const;

		double getScale(int t, int f, const std::vector<double>& x);

		void report(int iteration, double cost, const std::vector<double>& x) const;

};

#endif

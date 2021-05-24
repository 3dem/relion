#ifndef POLISH_PROGRAM_H
#define POLISH_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/tomography/motion/motion_fit.h>
#include <src/jaz/tomography/motion/modular_alignment/modular_alignment.h>
#include <src/jaz/tomography/motion/modular_alignment/GP_motion_model.h>
#include <src/jaz/tomography/motion/modular_alignment/no_motion_model.h>
#include <src/jaz/tomography/motion/modular_alignment/no_2D_deformation_model.h>
#include <src/jaz/tomography/motion/modular_alignment/spline_2D_deformation_model.h>
#include <src/jaz/optimization/lbfgs.h>
#include <vector>

#include "refinement.h"

class CTF;
class AberrationsCache;

class AlignProgram : public RefinementProgram
{
	public:
		
		
		AlignProgram(int argc, char *argv[]);
		
		
			bool
				do_motion, shiftOnly, whiten, whiten_abs, outputShiftedCCs,
				do_anisotropy, per_tilt_anisotropy;

			double padding, hiPass_px, sig2RampPower;
			int range, num_iters;
						
			ModularAlignmentSettings alignmentSettings;
			GPMotionModel::Settings motionSettings;
			GPMotionModel::MotionParameters motionParameters;
			

		void run();


	protected:

		void parseInput();

		void initialise();

		void finalise();

		void processTomograms(
				const std::vector<int>& tomoIndices,
				const AberrationsCache& aberrationsCache,
				bool per_tomogram_progress);


		std::string getTempFilenameRoot(
				const std::string& tomogram_name);

	private:

		void writeTempData(
				const std::vector<Trajectory>* traj,
				const std::vector<gravis::d4Matrix>& proj,
				const std::vector<gravis::d3Vector>& pos,
				int t);

		void readTempData(int t);

		void mergeLogFiles();
		
		template<class MotionModel>
		void performAlignment(
				MotionModel& motionModel,
				const std::vector<BufferedImage<double>>& CCs,
				const std::vector<gravis::d4Matrix>& projByTime,
				const Tomogram& tomogram,
				int tomo_index,
				int progress_bar_offset,
				bool per_tomogram_progress);
		
};

template<class MotionModel>
void AlignProgram::performAlignment(
		MotionModel& motionModel,
		const std::vector<BufferedImage<double>>& CCs,
		const std::vector<gravis::d4Matrix>& projByTime,
		const Tomogram& tomogram,
		int tomo_index,
		int progress_bar_offset,
		bool per_tomogram_progress)
{
	No2DDeformationModel noDeformationModel;

	ModularAlignment<MotionModel, No2DDeformationModel> alignment(
		CCs, projByTime, particleSet, particles[tomo_index],
		motionModel, noDeformationModel,
		alignmentSettings, tomogram,
		padding,
		progress_bar_offset, num_threads,
		per_tomogram_progress && verbosity > 0);

	std::vector<double> initial(alignment.getParamCount(), 0.0);


	if (verbosity > 0 && per_tomogram_progress)
	{
		Log::beginProgress("Performing optimisation", num_iters);
	}


	std::vector<double> opt = LBFGS::optimize(
		initial, alignment, 1, num_iters, 1e-4, 1e-5);


	if (verbosity > 0 && per_tomogram_progress)
	{
		Log::endProgress();
	}


	std::vector<gravis::d4Matrix> projections = alignment.getProjections(opt, tomogram.frameSequence);
	std::vector<gravis::d3Vector> positions = alignment.getParticlePositions(opt);
	
	if (do_motion)
	{
		std::vector<Trajectory> trajectories = alignment.exportTrajectories(
					opt, particleSet, tomogram.frameSequence);
	
		writeTempData(&trajectories, projections, positions, tomo_index);
	
		alignment.visualiseTrajectories2D(
			opt, 8.0, tomogram.name,
			getTempFilenameRoot(tomogram.name) + "_tracks");
	}
	else
	{
		writeTempData(0, projections, positions, tomo_index);
	}
}




#endif

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
#include <src/jaz/tomography/motion/modular_alignment/Fourier_2D_deformation_model.h>
#include <src/jaz/tomography/motion/modular_alignment/linear_2D_deformation_model.h>
#include <src/jaz/tomography/motion/proto_alignment.h>
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
				do_motion, shiftOnly, globalShift,
				whiten, whiten_abs, outputShiftedCCs,
				do_anisotropy, per_tilt_anisotropy,
				do_deformation, debug;

			double padding, hiPass_px, sig2RampPower, freqCutoffFract;
			int range, num_iters, min_frame, max_frame;
			
			std::string deformationType;

			ModularAlignmentSettings alignmentSettings;
			GPMotionModel::Parameters motionParameters;
			Deformation2D::Parameters deformationParameters;
			

		void run();


	protected:

			std::vector<std::vector<Trajectory>> allTrajectories;

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
		
		void writeTempAlignmentData(
				const std::vector<gravis::d4Matrix>& proj,
				const std::vector<gravis::d3Vector>& pos,
				int t);
		
		void writeTempMotionData(
				const std::vector<Trajectory>& traj,
				int t);
		
		void writeTempDeformationData(
				const std::vector<std::vector<double>>& def,
				int t);

		void readTempData(int t);

		void mergeLogFiles();

		template<class MotionModel>
		void performAlignment(
				MotionModel& motionModel,
				const std::vector<BufferedImage<double>>& CCs,
				const Tomogram& tomogram,
				int tomo_index,
				int progress_bar_offset,
				bool per_tomogram_progress);

		template<class MotionModel, class DeformationModel>
		void performAlignment(
				MotionModel& motionModel,
				DeformationModel& deformationModel,
				const std::vector<BufferedImage<double>>& CCs,
				const Tomogram& tomogram,
				int tomo_index,
				int progress_bar_offset,
				bool per_tomogram_progress);
		
};

template<class MotionModel>
void AlignProgram::performAlignment(
		MotionModel& motionModel,
		const std::vector<BufferedImage<double>>& CCs,
		const Tomogram& tomogram,
		int tomo_index,
		int progress_bar_offset,
		bool per_tomogram_progress)
{
	if (do_deformation)
	{
		if (deformationType == "spline")
		{
			Spline2DDeformationModel deformationModel(
				deformationParameters,
				gravis::i2Vector(tomogram.stack.xdim, tomogram.stack.ydim));
	
			performAlignment(
				motionModel,
				deformationModel,
				CCs,
				tomogram,
				tomo_index,
				progress_bar_offset,
				per_tomogram_progress);
		}
		else if (deformationType == "Fourier")
		{
			Fourier2DDeformationModel deformationModel(
				deformationParameters,
				gravis::i2Vector(tomogram.stack.xdim, tomogram.stack.ydim));
	
			performAlignment(
				motionModel,
				deformationModel,
				CCs,
				tomogram,
				tomo_index,
				progress_bar_offset,
				per_tomogram_progress);
		}
		else if (deformationType == "linear")
		{
			Linear2DDeformationModel deformationModel(
				deformationParameters,
				gravis::i2Vector(tomogram.stack.xdim, tomogram.stack.ydim));
	
			performAlignment(
				motionModel,
				deformationModel,
				CCs,
				tomogram,
				tomo_index,
				progress_bar_offset,
				per_tomogram_progress);
		}
	}
	else
	{
		No2DDeformationModel noDeformationModel;

		performAlignment(
			motionModel,
			noDeformationModel,
			CCs,
			tomogram,
			tomo_index,
			progress_bar_offset,
			per_tomogram_progress);
	}
}

template<class MotionModel, class DeformationModel>
void AlignProgram::performAlignment(
		MotionModel& motionModel,
		DeformationModel& deformationModel,
		const std::vector<BufferedImage<double>>& CCs,
		const Tomogram& tomogram,
		int tomo_index,
		int progress_bar_offset,
		bool per_tomogram_progress)
{
	ModularAlignment<MotionModel, DeformationModel> alignment(
		CCs, particleSet, particles[tomo_index],
		motionModel, deformationModel,
		alignmentSettings, tomogram,
		padding,
		progress_bar_offset, num_threads,
		per_tomogram_progress && verbosity > 0,
		min_frame, max_frame);

	std::vector<double> initial = alignment.originalCoefficients;

	alignment.devMode = debug;


	if (!debug && verbosity > 0 && per_tomogram_progress)
	{
		Log::beginProgress("Performing optimisation", num_iters);
	}

	/*std::cout << "initial info: \n";
	alignment.printDebuggingInfo(initial);*/

	std::vector<double> opt = LBFGS::optimize(
			initial, alignment, 1, num_iters, 1e-6, 1e-4);

	/*std::cout << "final info: \n";
	alignment.printDebuggingInfo(opt);*/

	if (verbosity > 0 && per_tomogram_progress)
	{
		if (!debug) Log::endProgress();

		const int it = alignment.lastIterationNumber;

		if (it >= num_iters)
		{
			Log::warn("Alignment did not converge after " + ZIO::itoa(it) + " iterations");
		}
		else
		{
			Log::print("Alignment converged after " + ZIO::itoa(it) + " iterations");
		}
	}

	std::vector<gravis::d4Matrix> projections = alignment.getProjections(opt, tomogram.frameSequence);
	std::vector<gravis::d3Vector> positions = alignment.getParticlePositions(opt);
	
	writeTempAlignmentData(projections, positions, tomo_index);
	
	if (do_motion)
	{
		std::vector<Trajectory> trajectories = alignment.exportTrajectories(
					opt, tomogram.frameSequence);
	
		writeTempMotionData(trajectories, tomo_index);
	
		alignment.visualiseTrajectories(
			opt, 8.0, tomogram.name,
			getTempFilenameRoot(tomogram.name) + "_tracks");
	}
	
	if (do_deformation)
	{
		std::vector<std::vector<double>> deformations = alignment.get2DDeformations(
					opt, tomogram.frameSequence);
		
		writeTempDeformationData(deformations, tomo_index);
		
		const gravis::i2Vector image_size(
					tomogram.stack.xdim, 
					tomogram.stack.ydim);
		
		const gravis::i2Vector grid_size(
					deformationParameters.grid_width, 
					deformationParameters.grid_height);
		
		alignment.visualise2DDeformations(
					opt, image_size, grid_size, tomogram.name,
					getTempFilenameRoot(tomogram.name) + "_deformations");
	}

	alignment.visualiseShifts(
				opt, tomogram.name,
				getTempFilenameRoot(tomogram.name) + "_shifts");
}




#endif

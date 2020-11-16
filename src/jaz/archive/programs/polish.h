#ifndef POLISH_PROGRAM_H
#define POLISH_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/tomography/motion/motion_fit.h>
#include <vector>

#include "refinement.h"

class CTF;
class AberrationsCache;

class PolishProgram : public RefinementProgram
{
	public:
		
		
		PolishProgram(int argc, char *argv[]);
		
		
			bool whiten, whiten_abs, outputShiftedCCs;
			double padding, hiPass_px, sig2RampPower;
			int range, num_iters;
			
			MotionFit::MotionParameters motParams;
			MotionFit::Settings mfSettings;
		
		void readParams(IOParser& parser);
		void run();


	protected:

		void initialise();
		void finalise();

		void processTomograms(
				const std::vector<int>& tomoIndices,
				const AberrationsCache& aberrationsCache,
				int verbosity,
				bool per_tomogram_progress);


		std::string getTempFilenameRoot(
				const std::string& tomogram_name);

	private:

		void writeTempData(
				const std::vector<Trajectory>& traj,
				const std::vector<gravis::d4Matrix>& proj,
				const std::vector<gravis::d3Vector>& pos,
				int t);

		void readTempData(int t);
};


#endif

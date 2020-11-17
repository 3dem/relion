#ifndef PROTO_POLISH_PROGRAM_H
#define PROTO_POLISH_PROGRAM_H

#include <string>
#include <src/jaz/math/t_complex.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <vector>

#include "refinement.h"

class CTF;
class AberrationsCache;


class FrameAlignmentProgram : public RefinementProgram
{
	public:
		
		FrameAlignmentProgram(int argc, char *argv[]);
		
			bool shiftOnly, whiten, whiten_abs,  
				const_particles, const_angles, const_shifts;
			
			int range, num_iters;
			
			double padding, hiPass_px, sig2RampPower;
		
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
				const std::vector<gravis::d4Matrix>& proj,
				const std::vector<gravis::d3Vector>& pos,
				int t);

		void readTempData(int t);
};


#endif

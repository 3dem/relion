#ifndef TOMO_TEMPLATEPICKER_PROGRAM_H
#define TOMO_TEMPLATEPICKER_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/tomogram_set.h>

class Tomogram;
class ParticleSet;
class AberrationsCache;
class ParticleIndex;


class TemplatePickerProgram
{
	public:

		TemplatePickerProgram(){}

			OptimisationSet optimisation_set;
			double max_freq, max_angle, fiducials_radius_A;
			int num_threads;
			std::string template_filename, out_dir;
			BufferedImage<float> template_map_RS;
			BufferedImage<fComplex> template_map_FS;

			TomogramSet tomogramSet;


		void readBasicParameters(IOParser& parser, int argc, char *argv[]);
		virtual void readParameters(int argc, char *argv[]);
		virtual void run();


	protected:

		void initialise();

		void processTomograms(
				int first_t,
				int last_t,
				const TomogramSet& tomoSet,
				int verbosity);

		void pick(
				double rot, double tilt, double psi,
				const Tomogram& tomogram,
				const BufferedImage<fComplex>& framesFS,
				const BufferedImage<float>& maskedFramesSqRS,
				int num_threads);
};

#endif

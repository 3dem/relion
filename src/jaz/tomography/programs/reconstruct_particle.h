#ifndef BACKPROJECT_PROGRAM_H
#define BACKPROJECT_PROGRAM_H

#include <string>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/optimisation_set.h>

class TomogramSet;
class ParticleSet;
class AberrationsCache;
class ParticleIndex;


class ReconstructParticleProgram
{
	public:
		
		ReconstructParticleProgram();
		

			OptimisationSet optimisationSet;
			std::string outDir, symmName, tmpOutRoot;

			
			bool
				do_whiten, no_reconstruction, only_do_unfinished,
				run_from_GUI, run_from_MPI,
				no_backup, do_circle_crop;

			int boxSize, cropSize, num_threads, outer_threads, inner_threads, max_mem_GB;

			double SNR, taper, binning, freqCutoffFract;

			int nr_helical_asu;
			double helical_rise, helical_twist;


		void readBasicParameters(int argc, char *argv[]);
		virtual void readParameters(int argc, char *argv[]);
		virtual void run();


	protected:

		void processTomograms(
				const std::vector<int>& tomoIndices,
				const TomogramSet& tomoSet,
				const ParticleSet& dataSet,
				const std::vector<std::vector<ParticleIndex>>& particles,
				const AberrationsCache& aberrationsCache,
				std::vector<BufferedImage<dComplex>>& dataImgFS,
				std::vector<BufferedImage<double>>& ctfImgFS,
				const double binnedOutPixelSize,
				int s02D,
				bool do_ctf,
				bool flip_value,
				int verbosity,
				bool per_tomogram_progress);

		void finalise(
				std::vector<BufferedImage<dComplex>>& dataImgFS,
				std::vector<BufferedImage<double>>& ctfImgFS,
				const double binnedOutPixelSize);

		void symmetrise(
				std::vector<BufferedImage<dComplex>>& dataImgFS,
				std::vector<BufferedImage<double>>& ctfImgFS,
				double binnedOutPixelSize);

		void reconstruct(
				BufferedImage<double>& dataImgRS,
				BufferedImage<double>& dataImgDivRS,
				BufferedImage<double>& ctfImgFS,
				BufferedImage<dComplex>& dataImgFS);

		void writeOutput(
				const BufferedImage<double>& corrected,
				const BufferedImage<double>& data,
				const BufferedImage<double>& weight,
				const std::string& tag,
				double pixelSize);
};

#endif

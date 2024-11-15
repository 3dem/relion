#ifndef SUBTOMO_PROGRAM_H
#define SUBTOMO_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/projector.h>



class ParticleSet;
class ParticleIndex;
class TomogramSet;
class Tomogram;
class AberrationsCache;


class SubtomoProgram
{
	public:

		SubtomoProgram();
		
			OptimisationSet optimisationSet;

			std::string outDir;
			
			int 
				boxSize, 
				cropSize,
                min_frames,
                num_threads;
			
			double 
				SNR,
				binning,
				taper,
				env_sigma,
				cone_slope,
				cone_sig0,
				freqCutoffFract,
                maxDose;
			
			bool 
				flip_value, 
				diag,
                do_ctf,
                do_stack2d,
				do_whiten,
				do_center, 
				do_rotate, 
				do_cone_weight,
				do_gridding_precorrection,
				do_circle_crop,
				do_narrow_circle_crop,
				do_circle_precrop,
				do_sum_all,
				write_combined,
				write_ctf,
				write_multiplicity, 
				write_divided, 
				write_normalised,
				do_not_write_any,
				only_do_unfinished,
				apply_offsets,
                apply_orientations,
				write_float16,
				run_from_GUI,
				run_from_MPI,
                do_real_subtomo;

		void readBasicParameters(IOParser& parser);
		virtual void readParameters(int argc, char *argv[]);
		virtual void run();


	protected:

		void initialise(
				ParticleSet& particleSet,
				const std::vector<std::vector<ParticleIndex>>& particles,
				const TomogramSet& tomogramSet,
                bool verbose = true);

        BufferedImage<float> extractSubtomogramsAndReProject(
                ParticleIndex part_id, MultidimArray<RFLOAT> &recTomo,
                const Tomogram& tomogram, const ParticleSet &particleSet,
                const std::vector<bool> &isVisible, RFLOAT tomogram_angpix);

		void processTomograms(
				const std::vector<int>& tomoIndices,
				const TomogramSet& tomogramSet,
				const ParticleSet& particleSet,
				const std::vector<std::vector<ParticleIndex>>& particles,
				const AberrationsCache& aberrationsCache,
				long int s02D,
				long int s2D,
				long int s3D,
				double relative_box_scale,
				int verbose,
				BufferedImage<float>& sum_data,
				BufferedImage<float>& sum_weights );

private:

			bool directoriesPerTomogram;

		std::string getOutputFilename(
				ParticleIndex p,
				int tomogramIndex,
				const ParticleSet& particleSet,
				const TomogramSet& tomogramSet);

		void writeParticleSet(
				const ParticleSet& particleSet,
				const std::vector<std::vector<ParticleIndex>>& particles,
				const TomogramSet& tomogramSet);

		BufferedImage<float> cropAndTaper(
				const BufferedImage<float>& imgFS,
				int boundary,
				int num_threads) const;
};

#endif

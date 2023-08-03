#ifndef SUBTOMO_PROGRAM_H
#define SUBTOMO_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/optimisation_set.h>


class ParticleSet;
class ParticleIndex;
class TomogramSet;
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
				num_threads;
			
			double 
				SNR,
				binning, 
				taper, 
				env_sigma,
				cone_slope,
				cone_sig0,
				freqCutoffFract;
			
			bool 
				flip_value, 
				diag, 
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
				run_from_MPI;

		void readBasicParameters(IOParser& parser);
		virtual void readParameters(int argc, char *argv[]);
		virtual void run();


	protected:

		void initialise(
				const ParticleSet& particleSet,
				const std::vector<std::vector<ParticleIndex>>& particles,
				const TomogramSet& tomogramSet);

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
				bool do_ctf,
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

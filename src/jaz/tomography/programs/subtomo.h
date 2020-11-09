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

		SubtomoProgram(){}
		
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
				cone_sig0;
			
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
				do_sum_all,
				write_combined,
				write_ctf,
				write_multiplicity, 
				write_divided, 
				write_normalised;

		void readBasicParameters(IOParser& parser);
		virtual void readParameters(int argc, char *argv[]);
		virtual void run();


	protected:

		void writeParticleSet(
				const ParticleSet& particleSet,
				const std::vector<std::vector<ParticleIndex>>& particles);

		void processTomograms(
				int first_t,
				int last_t,
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

		BufferedImage<float> cropAndTaper(
				const BufferedImage<float>& imgFS,
				int boundary,
				int num_threads) const;
};

#endif

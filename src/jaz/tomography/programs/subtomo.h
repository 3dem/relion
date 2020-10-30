#ifndef SUBTOMO_PROGRAM_H
#define SUBTOMO_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/optimisation_set.h>

class SubtomoProgram
{
	public:
		
		SubtomoProgram(){}
		
			OptimisationSet optimisationSet;

			std::string outTag;
			
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
		
		void readParameters(int argc, char *argv[]);
		void run();	
		
		BufferedImage<float> cropAndTaper(
				const BufferedImage<float>& imgFS, 
				int boundary, 
				int num_threads) const;
};

#endif

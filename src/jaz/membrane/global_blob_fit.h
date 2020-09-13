#ifndef GLOBAL_BLOB_FIT_3D_H
#define GLOBAL_BLOB_FIT_3D_H

#include "blob_2d.h"
#include <src/jaz/image/buffered_image.h>


class GlobalBlobFit2D
{
	public:
		
		static std::pair<double,std::vector<double>> fit(
				const std::vector<double>& initial_parameters,
				double initial_radius, 
		        int max_frequencies,
		        const RawImage<float>& image,
		        const RawImage<float>& mask);
		
		static BufferedImage<float> drawOutline(
		        const std::vector<double> parameters,
		        double mean_radius,
		        const BufferedImage<float>& image);
};

#endif

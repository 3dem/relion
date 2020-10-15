#ifndef CELLS_H
#define CELLS_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t3Vector.h>
#include <vector>

class Cells
{
	public:
		
		
		static void findCenters(
				const std::vector<gravis::d3Vector>& surfacePoints,
				int bins,
				float maxRadius,
				BufferedImage<float>& outMaxima,
				BufferedImage<float>& outRadii,
				int num_threads);
		
		static void findCenters(
				const BufferedImage<float>& surface,
				int bins,
				float maxRadius,
				float threshVal,
				BufferedImage<float>& outMaxima,
				BufferedImage<float>& outRadii,
				int num_threads);
		
		static BufferedImage<float> findPointSymmetries(
				const BufferedImage<float>& img, 
				double highPass,
				double lowPass);
		
};

#endif

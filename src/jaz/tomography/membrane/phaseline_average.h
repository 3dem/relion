#ifndef PHASE_LINE_AVERAGE_H
#define PHASE_LINE_AVERAGE_H

#include <src/jaz/image/buffered_image.h>
#include <src/jaz/gravis/t2Vector.h>
#include <vector>

class PhaseLineAverage
{
	public:
		
		static std::vector<double> averageNN(
				const RawImage<float>& data, 
				const RawImage<float>& phase,
				int bins);
		
		static BufferedImage<float> expandLIN(
				const RawImage<float>& data, 
				const RawImage<float>& phase, 
				const std::vector<double>& avg);
		
		static gravis::f2Vector findBounds(const RawImage<float>& phase);
};


#endif

#ifndef REFERENCE_MAP_H
#define REFERENCE_MAP_H

#include <src/jaz/image/buffered_image.h>
#include <vector>


class ReferenceMap
{
	public:
		
		ReferenceMap();
		
		ReferenceMap(
			std::string ref1Fn, std::string ref2Fn, int boxSize, 
			std::string maskFn, std::string fscFn,
			bool rcThresh, double threshWidth);
		
		
			std::vector<BufferedImage<float>> image_real;
			BufferedImage<float> mask, freqWeight;
			std::vector<BufferedImage<fComplex>> image_FS;
			
		int getBoxSize() const;
};

#endif

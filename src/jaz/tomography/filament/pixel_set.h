#ifndef PIXEL_SET_H
#define PIXEL_SET_H

#include <vector>

class PixelSet
{
	public:
		
		std::vector<std::vector<float>> splineCoords, pixelValue;
		std::vector<std::vector<long long int>> pixelIndex;
		
		long long int totalNumber;
};

#endif

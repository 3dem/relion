#include "color_helper.h"
#include <cmath>

using namespace gravis;

dRGB ColorHelper::signedToRedBlue(double d, double rbScale)
{		
	return dRGB(
		std::min(1.0, rbScale * std::max(0.0, d) ), 
		std::min(1.0, std::max(0.0, (rbScale * std::abs(d) - 1.0) / (rbScale - 1.0) ) ), 
		std::min(1.0, rbScale * std::max(0.0, -d) ) );
}

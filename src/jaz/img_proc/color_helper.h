#ifndef COLOR_HELPER_H
#define COLOR_HELPER_H

#include <src/jaz/gravis/tRGB.h>

class ColorHelper
{
	public:
		
		static gravis::dRGB signedToRedBlue(double d, double rbScale = 2.0);
};

#endif

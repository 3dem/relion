#include "color_helper.h"
#include <src/macros.h>
#include <cmath>

#ifdef HAVE_PNG
#include <src/jaz/gravis/tImage.h>
#endif

using namespace gravis;

dRGB ColorHelper::signedToRedBlue(double d, double scale, double rbFract)
{		
	const double d_rb = d / (scale * rbFract);
	const double d_g = (std::abs(d)/scale - rbFract) / (1.0 - rbFract);
	
	return dRGB(
		std::min(1.0, std::max(0.0,  d_rb)), 
		std::min(1.0, std::max(0.0,  d_g)), 
		std::min(1.0, std::max(0.0, -d_rb)) );
}

void ColorHelper::writeAngleToPNG(const Image<RFLOAT> &img, std::string filename)
{
	writeSignedToPNG(img, filename+"_[-pi,+pi]", PI);
	writeSignedToPNG(img, filename+"_[-1,+1]", 1.0);
}

void ColorHelper::writeSignedToPNG(const Image<RFLOAT> &img, std::string filename, double scale)
{
	#ifdef HAVE_PNG
	{
		tImage<dRGB> pngOut(img.data.xdim, img.data.ydim);
		pngOut.fill(dRGB(0.f));
		
		for (int y = 0; y < img.data.ydim; y++)
		for (int x = 0; x < img.data.xdim; x++)
		{
			double c = img(y,x);
			pngOut(x,y) = signedToRedBlue(c, scale);
		}
		
		pngOut.writePNG(filename+".png");
	}
	#endif
}

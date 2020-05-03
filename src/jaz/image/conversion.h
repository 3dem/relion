#ifndef IMAGE_CONVERSION_H
#define IMAGE_CONVERSION_H

#include "buffered_image.h"

class Conversion
{
	public:
		
		template <typename T>
		static BufferedImage<double> toDouble(const RawImage<T>& img);
};

template <typename T>
BufferedImage<double> Conversion::toDouble(const RawImage<T>& img)
{
	BufferedImage<double> out(img.xdim, img.ydim, img.zdim);
	
	for (int z = 0; z < out.zdim; z++)
	for (int y = 0; y < out.ydim; y++)
	for (int x = 0; x < out.xdim; x++)
	{
		out(x,y,z) = (double) img(x,y,z);
	}
	
	return out;
}


#endif

#include "image_file_helper.h"
#include <src/image.h>

using namespace gravis;

t3Vector<long int> ImageFileHelper::getSize(const std::string &filename)
{
	Image<RFLOAT> img;
	img.read(filename, false);
	
	return t3Vector<long int>(img.data.xdim, img.data.ydim, img.data.zdim * img.data.ndim);
}

double ImageFileHelper::getSamplingRate(const std::string &filename)
{
	Image<RFLOAT> img;
	img.read(filename, false);

	return (double) img.samplingRateX();
}

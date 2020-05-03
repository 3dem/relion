#ifndef COLOR_HELPER_H
#define COLOR_HELPER_H

#include <src/jaz/gravis/tRGB.h>
#include <src/image.h>
#include "raw_image.h"

class ColorHelper
{
	public:
		
		static gravis::dRGB signedToRedBlue(double d, double scale = 1.0, double rbFract = 0.333);
		
		static void writeAngleToPNG(const Image<RFLOAT>& img, std::string filename);
		
		static void writeSignedToPNG(const Image<RFLOAT>& img, std::string filename, 
									 double scale = 1.0);
		
		static void writeSignedToEPS(std::string filename, int col, const std::vector<Image<RFLOAT> > &imgs,
				const std::vector<double> &scales, const std::vector<std::string> &labels);
		
		template <class T>
		static void writeSignedToPNG(
				const RawImage<T>& img, std::string filename);
		
		template <class T>
		static void writeSignedToPNG(
				const RawImage<T>& img, std::string filename, double scale);
};

template <class T>
void ColorHelper::writeSignedToPNG(
		const RawImage<T>& img, std::string filename)
{
	const int w = img.xdim;
	const int h = img.ydim;
	
	double vmax = 0.0;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double a = std::abs(img(x,y));		
		if (a > vmax) vmax = a;
	}
	
	writeSignedToPNG(img, filename, vmax);
}

template <class T>
void ColorHelper::writeSignedToPNG(
		const RawImage<T>& img, std::string filename, double scale)
{
	gravis::tImage<gravis::dRGB> pngOut(img.xdim, img.ydim);
	pngOut.fill(gravis::dRGB(0.f));
	
	for (int y = 0; y < img.ydim; y++)
	for (int x = 0; x < img.xdim; x++)
	{
		double c = img(x,y);
		pngOut(x,y) = signedToRedBlue(c, scale);
	}
	
	pngOut.writePNG(filename);
}

#endif

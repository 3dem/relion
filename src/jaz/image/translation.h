#ifndef TRANSLATION_H
#define TRANSLATION_H

#include "raw_image.h"

class Translation
{
	public: 
		
		template <typename T>
		static void shiftInFourierSpace2D(RawImage<tComplex<T>>& img, double dx, double dy);
		
		template <typename T>
		static void shiftByHalf(RawImage<tComplex<T>>& img);
};

template <typename T>
void Translation::shiftInFourierSpace2D(RawImage<tComplex<T>>& img, double dx, double dy)
{
	const int w = img.xdim;
	const int h = img.ydim;
	
	const double w0 = 2 * (w - 1);
	const double h0 = h;
	
	const double xshift = -dx / w0;
	const double yshift = -dy / h0;
	
	
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		const double yy = y < h/2? y : y - h;
		const double dotp = 2 * PI * (x * xshift + yy * yshift);

		double a, b;
		SINCOS(dotp, &b, &a);

		const double c = img(x,y).real;
		const double d = img(x,y).imag;
		const double ac = a * c;
		const double bd = b * d;
		const double ab_cd = (a + b) * (c + d);
		
		img(x,y) = tComplex<T>(T(ac - bd), T(ab_cd - ac - bd));
	}
}

template <typename T>
void Translation::shiftByHalf(RawImage<tComplex<T>>& img)
{
	for (long int z = 0; z < img.zdim; z++)
	for (long int y = 0; y < img.ydim; y++)
	for (long int x = 0; x < img.xdim; x++)
	{
		img(x,y,z) *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
	}
}

#endif

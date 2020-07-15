#ifndef IMAGE_NOISE_H
#define IMAGE_NOISE_H

#include "buffered_image.h"
#include <src/complex.h>
#include <src/jaz/math/fft.h>

class ImageNoise
{
	public:
		
		template <typename T>
		static BufferedImage<tComplex<T>> generateSmoothFourierNoise(int size, double falloff_sigma2);
};

template<typename T>
BufferedImage<tComplex<T>> ImageNoise::generateSmoothFourierNoise(int size, double falloff_sigma2)
{
	const int s = size;
	const int sh = s/2 + 1;
	
	BufferedImage<tComplex<T>> out(sh,s);
	BufferedImage<T> noise(s,s);
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double xx = x < s/2? x : x - s;
		const double yy = y < s/2? y : y - s;
		const double r2 = xx*xx + yy*yy;
		
		noise(x,y) = 2 * exp(-0.5*r2/falloff_sigma2) * (rand() / (double)RAND_MAX + 0.5);
	}
	
	FFT::FourierTransform(noise, out, FFT::Both);
	
	return out;
}

#endif

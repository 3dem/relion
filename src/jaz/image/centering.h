#ifndef CENTERING_H
#define CENTERING_H

#include "raw_image.h"
#include "buffered_image.h"

class Centering
{
	public:
		
		template<class T>
		static BufferedImage<T> humanFullToFftwFull(const RawImage<T>& img);
		
		template<class T>
		static BufferedImage<T> humanFullToFftwHalf(const RawImage<T>& img);
		
		template<class T>
		static BufferedImage<T> fftwFullToHumanFull(const RawImage<T>& img);
		
		template<class T>
		static BufferedImage<T> fftwHalfToHumanFull(const RawImage<T>& img);
		
		template<class T>
		static BufferedImage<T> fftwHalfAntisymmetricalToHumanFull(const RawImage<T>& img);
		
		template<class T>
		static BufferedImage<T> fftwHalfToFftwFull(const RawImage<T>& img);
		
		template<class T>
		static BufferedImage<T> fftwHalfAntisymmetricalToFftwFull(const RawImage<T>& img);
		
		template<class T>
		static void shiftInSitu(RawImage<tComplex<T>>& img);
		
};


template<class T>
BufferedImage<T> Centering::humanFullToFftwFull(const RawImage<T>& img)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	const size_t cx = w/2;
	const size_t cy = h/2;
	const size_t cz = d/2;
	
	BufferedImage<T> out(w,h,d);
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		const size_t xx = (x + cx) % w;
		const size_t yy = (y + cy) % h;
		const size_t zz = (z + cz) % d;
		
		out(x,y,z) = img(xx,yy,zz);
	}
	
	return out;
}

template<class T>
BufferedImage<T> Centering::humanFullToFftwHalf(const RawImage<T>& img)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	const size_t wh = w/2 + 1;
	
	const size_t cx = w/2;
	const size_t cy = h/2;
	const size_t cz = d/2;
	
	BufferedImage<T> out(wh,h,d);
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < wh; x++)
	{
		const size_t xx = (x + cx) % w;
		const size_t yy = (y + cy) % h;
		const size_t zz = (z + cz) % d;
		
		out(x,y,z) = img(xx,yy,zz);
	}
	
	return out;
}

template<class T>
BufferedImage<T> Centering::fftwFullToHumanFull(const RawImage<T>& img)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	const size_t cx = w - w/2;
	const size_t cy = h - h/2;
	const size_t cz = d - d/2;
	
	BufferedImage<T> out(w,h,d);
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		const size_t xx = (x + cx) % w;
		const size_t yy = (y + cy) % h;
		const size_t zz = (z + cz) % d;
		
		out(x,y,z) = img(xx,yy,zz);
	}
	
	return out;
}

template<class T>
BufferedImage<T> Centering::fftwHalfToFftwFull(const RawImage<T>& img)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	const int wd = 2 * (w-1);
	
	BufferedImage<T> out(wd,h,d);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		T v = img(x,y,z);
		
		out(x,y,z) = v;
		
		if (x > 0)
		{
			out(wd-x, (h-y)%h, (d-z)%d) = v;
		}
	}
	
	return out;
}

template<class T>
BufferedImage<T> Centering::fftwHalfAntisymmetricalToFftwFull(const RawImage<T>& img)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	const int wd = 2 * (w-1);
	
	BufferedImage<T> out(wd,h,d);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		T v = img(x,y,z);
		
		out(x,y,z) = v;
		
		if (x > 0)
		{
			out(wd-x, (h-y)%h, (d-z)%d) = -v;
		}
	}
	
	return out;
}

template<class T>
BufferedImage<T> Centering::fftwHalfToHumanFull(const RawImage<T>& img)
{
	return fftwFullToHumanFull(fftwHalfToFftwFull(img));
}

template<class T>
BufferedImage<T> Centering::fftwHalfAntisymmetricalToHumanFull(const RawImage<T>& img)
{
	return fftwFullToHumanFull(fftwHalfAntisymmetricalToFftwFull(img));
}


template<class T>
void Centering::shiftInSitu(RawImage<tComplex<T>>& img)
{
	for (long int z = 0; z < img.zdim; z++)
	for (long int y = 0; y < img.ydim; y++)
	for (long int x = 0; x < img.xdim; x++)
	{
		img(x,y,z) *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
	}
}

#endif

#ifndef LOCAL_EXTREMA_H
#define LOCAL_EXTREMA_H

#include "raw_image.h"
#include "buffered_image.h"

class LocalExtrema
{
	public:
		
		template <typename T>
		static BufferedImage<T> boxMaxima(const RawImage<T>& img, int dist);
		
		template <typename T>
		static BufferedImage<T> boxMinima(const RawImage<T>& img, int dist);
		
		template <typename T>
		static BufferedImage<T> boxExtrema(const RawImage<T>& img, int dist, double sign);

		template <typename T>
		static std::vector<gravis::d2Vector> discretePoints2D(
					const RawImage<T>& originalImage,
					const RawImage<T>& boxMaximaImage,
					T threshold);
};

template <typename T>
BufferedImage<T> LocalExtrema::boxMaxima(const RawImage<T>& img, int dist)
{
	return boxExtrema(img, dist, 1.0);
}

template <typename T>
BufferedImage<T> LocalExtrema::boxMinima(const RawImage<T>& img, int dist)
{
	return boxExtrema(img, dist, -1.0);
}

template <typename T>
BufferedImage<T> LocalExtrema::boxExtrema(const RawImage<T>& img, int dist, double sign)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	BufferedImage<T> img1 = img, img2 = img;
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		for (int i = -dist; i <= dist; i++)
		{
			const size_t xx = x + i;
			
			if (xx >= 0 && xx < w && sign * img(xx,y,z) > sign * img1(x,y,z))
			{
				img1(x,y,z) = img(xx,y,z);
			}
		}
	}
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		for (int i = -dist; i <= dist; i++)
		{
			const size_t yy = y + i;
			
			if (yy >= 0 && yy < h && sign * img1(x,yy,z) > sign * img2(x,y,z))
			{
				img2(x,y,z) = img1(x,yy,z);
			}
		}
	}
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		for (int i = -dist; i <= dist; i++)
		{
			const size_t zz = z + i;
			
			if (zz >= 0 && zz < d && sign * img2(x,y,zz) > sign * img1(x,y,z))
			{
				img1(x,y,z) = img2(x,y,zz);
			}
		}
	}
	
	return img1;
}

template<typename T>
std::vector<gravis::d2Vector> LocalExtrema::discretePoints2D(
		const RawImage<T>& originalImage,
		const RawImage<T>& boxMaximaImage,
		T threshold)
{
	std::vector<gravis::d2Vector> out;

	const int w = boxMaximaImage.xdim;
	const int h = boxMaximaImage.ydim;

	out.reserve(w);

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const T v0 = originalImage(x,y);
		const T v1 = boxMaximaImage(x,y);

		if (v0 > threshold && v0 == v1)
		{
			out.push_back(gravis::d2Vector(x,y));
		}
	}

	return out;
}

#endif

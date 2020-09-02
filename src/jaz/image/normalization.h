#ifndef IMG_NORMALIZATION_H
#define IMG_NORMALIZATION_H

#include "raw_image.h"
#include "stack_helper.h"
#include <src/jaz/util/zio.h>

class Normalization
{
    public:
		
		template <class T>
		static BufferedImage<T> toUnitInterval(const RawImage<T> &img);
		
		template <class T>
		static BufferedImage<T> byNormalDist(const RawImage<T> &img);
		
		template <class T>
		static BufferedImage<T> byNormalDistByFrame(const RawImage<T> &img);
		
		template <class T>
		static T computeMean(const RawImage<T>& img);
		
		template <class ImageType, class MaskType>
		static ImageType computeWeightedMean(
				const RawImage<ImageType>& img, const RawImage<MaskType>& mask);
		
		template <class ImageType, class MaskType>
		static ImageType computeWeightedMeanFromWeighted(
				const RawImage<ImageType>& img, const RawImage<MaskType>& mask);
		
		template <class T>
		static T computeVariance(const RawImage<T>& img, T mean);
		
		template <class T>
		static void zeroDC_stack(RawImage<T>& img);
};

template <class T>
BufferedImage<T> Normalization::toUnitInterval(const RawImage<T>& img)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	T minVal = std::numeric_limits<T>::max();
	T maxVal = -std::numeric_limits<T>::max();
			
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		T v = img(x,y,z);
		
		if (v > maxVal) maxVal = v;
		if (v < minVal) minVal = v;
	}
	
	BufferedImage<T> out(w,h,d);
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		T v = img(x,y,z);
		out(x,y,z) = (v - minVal)/(maxVal - minVal);
	}
	
	return out;
}

template <class T>
BufferedImage<T> Normalization::byNormalDist(const RawImage<T>& img)
{
	T mean = computeMean(img);
	T var = computeVariance(img, mean);
	
	BufferedImage<T> out = img - mean;
	out /= sqrt(var);
	
	return out;
}

template <class T>
BufferedImage<T> Normalization::byNormalDistByFrame(const RawImage<T>& img)
{
	BufferedImage<T> out(img.xdim, img.ydim, img.zdim);
	
	const int fc = img.zdim;
	
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<T> slice = NewStackHelper::extractSliceZ(img,f);
		
		T mean = computeMean(slice);
		T var = computeVariance(slice, mean);
		
		slice -= mean;		
		slice /= sqrt(var);
		
		out.getSliceRef(f).copyFrom(slice);
	}
	
	return out;
}

template <class T>
T Normalization::computeMean(const RawImage<T>& img)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	T sum(0);
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		sum += img(x,y,z);
	}
	
	return sum / (w*h*d);
}

template <class DataType, class MaskType>
DataType Normalization::computeWeightedMean(const RawImage<DataType>& img, const RawImage<MaskType>& mask)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	DataType sum(0);
	MaskType wgh(0);
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		sum += mask(x,y,z) * img(x,y,z);
		wgh += mask(x,y,z);
	}
	
	return wgh != MaskType(0)? sum / wgh : DataType(0);
}

template <class DataType, class MaskType>
DataType Normalization::computeWeightedMeanFromWeighted(const RawImage<DataType>& img, const RawImage<MaskType>& mask)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	DataType sum(0);
	MaskType wgh(0);
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		sum += mask(x,y,z) * img(x,y,z);
		wgh += mask(x,y,z) * mask(x,y,z);
	}
	
	return wgh != MaskType(0)? sum / wgh : DataType(0);
}

template <class T>
T Normalization::computeVariance(const RawImage<T>& img, T mean)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t d = img.zdim;
	
	T sum(0);
	
	for (size_t z = 0; z < d; z++)
	for (size_t y = 0; y < h; y++)
	for (size_t x = 0; x < w; x++)
	{
		T t = img(x,y,z) - mean;
		sum += t*t;
	}
	
	return sum / (w*h*d - 1);
}

template <class T>
void Normalization::zeroDC_stack(RawImage<T>& img)
{
	const size_t w = img.xdim;
	const size_t h = img.ydim;
	const size_t fc = img.zdim;
	
	for (size_t f = 0; f < fc; f++)
	{
		double sum = 0.0;
		
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			sum += img(x,y,f);
		}
		
		const double mean = sum / (w*h);
		
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			img(x,y,f) -= mean;
		}
	}
	
}

#endif

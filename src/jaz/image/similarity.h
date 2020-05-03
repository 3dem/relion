#ifndef SIMILARITY_H
#define SIMILARITY_H

#include "raw_image.h"
#include <src/jaz/math/fft.h>
#include <src/macros.h>
#include "padding.h"

class Similarity
{
	public:
		
		template<typename T>
		static BufferedImage<T> convolution(BufferedImage<T>& i0, BufferedImage<T>& i1);
		
		template<typename T>
		static BufferedImage<T> CC(BufferedImage<T>& i0, BufferedImage<T>& i1, bool zeroDC = false);
		
		template<typename T>
		static BufferedImage<T> CC(const RawImage<tComplex<T>>& I0, BufferedImage<T>& i1, bool zeroDC = false);
		
		template<typename T>
		static BufferedImage<T> NCC(BufferedImage<T>& i0, BufferedImage<T>& i1, bool zeroDC = false);
		
		
		template<typename T>
		static BufferedImage<T> CC_2D(BufferedImage<T>& ref, BufferedImage<T>& data, bool centerValue = true);
		
		template<typename T>
		static BufferedImage<T> CC_2D(const RawImage<tComplex<T>>& ref, BufferedImage<T>& data, bool centerValue = true);
		
		template<typename T>
		static BufferedImage<T> NCC_2D(BufferedImage<T>& ref, BufferedImage<T>& data, bool centerValue = true);
		
		template<typename T>
		static BufferedImage<T> weightedNCC_2D(const BufferedImage<T>& ref, const BufferedImage<T>& mask, BufferedImage<T>& data, double eps = 1e-4);
		
		template<typename T>
		static BufferedImage<T> weightedL2_2D(const BufferedImage<T>& ref, const BufferedImage<T>& mask, BufferedImage<T>& data);
		
		
		template<typename T>
		static BufferedImage<T> centerReference(const BufferedImage<T>& ref, const BufferedImage<T>& mask);
};

template<typename T>
BufferedImage<T> Similarity::convolution(BufferedImage<T>& i0, BufferedImage<T>& i1)
{
	const int w = i0.xdim;
	const int h = i0.ydim;
	const int d = i0.zdim;
	
	const int wh = w/2 + 1;
	
	const T scale = w * h * d;
	
	BufferedImage<tComplex<T>> I0, I1;
	
	FFT::FourierTransform(i0, I0);
	FFT::FourierTransform(i1, I1);
			
	BufferedImage<tComplex<T>> CC(wh,h,d);
			
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		DIRECT_NZYX_ELEM(CC.data, 0, z, y, x) = 
			scale *
			DIRECT_NZYX_ELEM(I0.data, 0, z, y, x) * 
			DIRECT_NZYX_ELEM(I1.data, 0, z, y, x);
	}
	
	BufferedImage<T> cc;
	FFT::inverseFourierTransform(CC, cc);
	
	return cc;
}

template<typename T>
BufferedImage<T> Similarity::CC(BufferedImage<T>& i0, BufferedImage<T>& i1, bool zeroDC)
{
	const int w = i0.xdim;
	const int h = i0.ydim;
	const int d = i0.zdim;
	
	const int wh = w/2 + 1;
	
	const T scale = w * h * d;
	
	BufferedImage<tComplex<T>> I0, I1;
	
	FFT::FourierTransform(i0, I0);
	FFT::FourierTransform(i1, I1);
			
	BufferedImage<tComplex<T>> CC(wh,h,d);
			
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		CC(x,y,z) = scale * I0(x,y,z).conj() * I1(x,y,z);
	}
	
	if (zeroDC)
	{
		CC(0,0,0) = tComplex<T>(0,0);
	}
	
	BufferedImage<T> cc;
	FFT::inverseFourierTransform(CC, cc);
	
	return cc;
}

template<typename T>
BufferedImage<T> Similarity::CC(const RawImage<tComplex<T>>& I0, BufferedImage<T>& i1, bool zeroDC)
{
	const int w = i1.xdim;
	const int h = i1.ydim;
	const int d = i1.zdim;
	
	const int wh = w/2 + 1;
	
	const T scale = w * h * d;
	
	BufferedImage<tComplex<T>> I1;
	
	FFT::FourierTransform(i1, I1);
			
	BufferedImage<tComplex<T>> CC(wh,h,d);
			
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		CC(x,y,z) = scale * I0(x,y,z).conj() * I1(x,y,z);
	}
	
	if (zeroDC)
	{
		CC(0,0,0) = tComplex<T>(0,0);
	}
	
	BufferedImage<T> cc;
	FFT::inverseFourierTransform(CC, cc);
	
	return cc;
}

template<typename T>
BufferedImage<T> Similarity::NCC(BufferedImage<T>& i0, BufferedImage<T>& i1, bool zeroDC)
{
	const int w = i0.xdim;
	const int h = i0.ydim;
	const int d = i0.zdim;
	
	const int wh = w/2 + 1;
	
	double p0(0.0), p1(0.0);
	
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		const double v0 = i0(x,y,z);
		const double v1 = i1(x,y,z);
		
		p0 += v0 * v0;
		p1 += v1 * v1;
	}
	
	const T scale = w * h * d / sqrt(p0 * p1);
	
	BufferedImage<tComplex<T>> I0, I1;
	
	FFT::FourierTransform(i0, I0);
	FFT::FourierTransform(i1, I1);
			
	BufferedImage<tComplex<T>> NCC(wh,h,d);
			
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		NCC(x,y,z) = scale * I0(x,y,z).conj() * I1(x,y,z);
	}
	
	if (zeroDC)
	{
		NCC(0,0,0) = tComplex<T>(0,0);
	}
	
	BufferedImage<T> ncc;
	FFT::inverseFourierTransform(NCC, ncc);
	
	return ncc;
}

template<typename T>
BufferedImage<T> Similarity::CC_2D(BufferedImage<T>& ref, BufferedImage<T>& data, bool centerValue)
{
	const int wr = ref.xdim;
	const int hr = ref.ydim;
	const int dr = ref.zdim;
	
	const int wd = data.xdim;
	const int hd = data.ydim;
	const int dd = data.zdim;
	
	if (dr != 1 || dd != 1)
	{
		REPORT_ERROR("Similarity::CC_2D: images are not 2D");
	}
	
	if (wr != wd || hr != hd)
	{
		BufferedImage<T> paddedRef = Padding::padCorner2D_full(ref, wd, hd);	
		return CC(paddedRef, data, centerValue);
	}
	else
	{
		return CC(ref, data, centerValue);
	}
}

template<typename T>
BufferedImage<T> Similarity::CC_2D(const RawImage<tComplex<T>>& ref, BufferedImage<T>& data, bool centerValue)
{
	const int wrh = ref.xdim;
	const int hr = ref.ydim;
	const int dr = ref.zdim;
	
	const int wr = 2*(wrh - 1);
	
	const int wd = data.xdim;
	const int hd = data.ydim;
	const int dd = data.zdim;
	
	if (dr != 1 || dd != 1)
	{
		REPORT_ERROR("Similarity::CC_2D: images are not 2D");
	}
	
	if (wr != wd || hr != hd)
	{
		REPORT_ERROR_STR("Similarity::CC_2D, complex ref.: images are of incorrect size: "
						 << wrh << "[" << wr << "]x" << hr << "x" << dr << " vs. "
						 << wd << "x" << hd << "x" << dd);
	}
	else
	{
		return CC(ref, data, centerValue);
	}
}

template<typename T>
BufferedImage<T> Similarity::NCC_2D(BufferedImage<T>& ref, BufferedImage<T>& data, bool centerValue)
{
	const int wr = ref.xdim;
	const int hr = ref.ydim;
	const int dr = ref.zdim;
	
	const int wd = data.xdim;
	const int hd = data.ydim;
	const int dd = data.zdim;
	
	if (dr != 1 || dd != 1)
	{
		REPORT_ERROR("Similarity::NCC_2D: images are not 2D");
	}
	
	if (wr != wd || hr != hd)
	{
		BufferedImage<T> paddedRef = Padding::padCorner2D_full(ref, wd, hd);	
		return NCC(paddedRef, data, centerValue);
	}
	else
	{
		return NCC(ref, data, centerValue);
	}
}

template<typename T>
BufferedImage<T> Similarity::weightedNCC_2D(
		const BufferedImage<T>& ref, const BufferedImage<T>& mask, BufferedImage<T>& data, double eps)
{
	const int wr = ref.xdim;
	const int hr = ref.ydim;
	const int dr = ref.zdim;
	
	const int wm = mask.xdim;
	const int hm = mask.ydim;
	const int dm = mask.zdim;
	
	const int wd = data.xdim;
	const int hd = data.ydim;
	const int dd = data.zdim;
	
	if (wm != wr || hm != hr)
	{
		REPORT_ERROR("Similarity::weightedNCC_2D: mask and reference have unequal size");
	}
	
	if (dr != 1 || dd != 1 || dm != 1)
	{
		REPORT_ERROR("Similarity::weightedNCC_2D: images are not 2D");
	}
	
	BufferedImage<T> data2;
	data2 = data * data;
	
	BufferedImage<T> centRef = centerReference(ref, mask);
	
	BufferedImage<T> paddedRef = Padding::padCorner2D_full(centRef, wd, hd);
	BufferedImage<T> paddedMask = Padding::padCorner2D_full(mask, wd, hd);
		
	BufferedImage<T> cc = CC(paddedRef, data, true);	
	BufferedImage<T> dataPower = CC(paddedMask, data2);
		
	BufferedImage<T> ncc(wd,hd);
	
	for (int y = 0; y < hd; y++)
	for (int x = 0; x < wd; x++)
	{
		ncc(x,y) = cc(x,y) / sqrt(dataPower(x,y) + eps);
	}
	
	return ncc;
}

template<typename T>
BufferedImage<T> Similarity::weightedL2_2D(
		const BufferedImage<T>& ref, const BufferedImage<T>& mask, BufferedImage<T>& data)
{
	const int wr = ref.xdim;
	const int hr = ref.ydim;
	const int dr = ref.zdim;
	
	const int wm = mask.xdim;
	const int hm = mask.ydim;
	const int dm = mask.zdim;
	
	const int wd = data.xdim;
	const int hd = data.ydim;
	const int dd = data.zdim;
	
	if (wm != wr || hm != hr)
	{
		REPORT_ERROR("Similarity::weightedL2_2D: mask and reference have unequal size");
	}
	
	if (dr != 1 || dd != 1 || dm != 1)
	{
		REPORT_ERROR("Similarity::weightedL2_2D: images are not 2D");
	}
	
	BufferedImage<T> data2;
	data2 = data * data;
	
	BufferedImage<T> paddedRef = Padding::padCorner2D_full(ref, wd, hd);
	BufferedImage<T> paddedMask = Padding::padCorner2D_full(mask, wd, hd);
		
	BufferedImage<T> cc = CC(paddedRef, data, false);	
	BufferedImage<T> dataPower = CC(paddedMask, data2, false);
		
	double refPower = 0.0;
	
	for (int y = 0; y < hr; y++)
	for (int x = 0; x < wr; x++)
	{
		refPower += mask(x,y) * ref(x,y) * ref(x,y);
	}
	
	BufferedImage<T> wl2(wd,hd);
	
	for (int y = 0; y < hd; y++)
	for (int x = 0; x < wd; x++)
	{
		wl2(x,y) = refPower - 2.0 * cc(x,y) + dataPower(x,y);
	}
	
	return wl2;
}

template<typename T>
BufferedImage<T> Similarity::centerReference(const BufferedImage<T>& ref, const BufferedImage<T>& mask)
{
	const int w = ref.xdim;
	const int h = ref.ydim;
	const int d = ref.zdim;
	
	BufferedImage<T> out(w,h,d);
	
	T refAvg = 0.0, maskSum = 0.0;
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		refAvg += mask(x,y,z) * ref(x,y,z);
		maskSum += mask(x,y,z);
	}
	
	refAvg /= maskSum;
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y,z) = mask(x,y,z) * (ref(x,y,z) - refAvg);
	}
	
	return out;
}

#endif

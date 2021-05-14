#ifndef RESAMPLING_H
#define RESAMPLING_H

#include "raw_image.h"
#include "stack_helper.h"
#include "interpolation.h"
#include "cutting.h"
#include "filter.h"
#include <src/jaz/math/fft.h>

class Resampling
{
	public:
		
		template <typename T>
		static BufferedImage<T> upsampleLinear_2D_half(
				const RawImage<T>& img, int wout, int hout);
		
		template <typename T>
		static BufferedImage<T> upsampleLinear_3D_full(
				const RawImage<T>& img, double factor, int wout = -1, int hout = -1, int dout = -1);
		
		template <typename T>
		static BufferedImage<T> upsampleCubic_2D_full(
				const RawImage<T>& img, int z, double factor, int wout = -1, int hout = -1);
		
		template <typename T>
		static BufferedImage<T> upsampleCubic_Stack_full(
				const RawImage<T>& stack, double factor, int wout = -1, int hout = -1);
		
		template <typename T>
		static BufferedImage<T> downsampleFilt_2D_full(
				BufferedImage<T>& img, int wout, int hout, double factor = -1);
		
		template <typename T>
		static BufferedImage<T> downsampleFilt_3D_full(
				BufferedImage<T>& img, int wout, int hout, int dout, double factor = -1);
		
		template <typename T>
		static BufferedImage<T> downsampleGaussFilt_3D_full(
				BufferedImage<T>& img, int wout, int hout, int dout, double factor = -1);
		
		template <typename T>
		static BufferedImage<T> downsampleFiltStack_2D_full(
				BufferedImage<T>& img, double factor, int num_threads, bool keepSizeEven = true);
		
		template <typename T>
		static BufferedImage<T> downsampleMax_2D_full(
				const RawImage<T>& img, int wout, int hout);
		
		template <typename T>
		static BufferedImage<T> lowPassFilter_full(
				BufferedImage<T>& img, double maxFreq0, double maxFreq1);
		
		template <typename T>
		static BufferedImage<T> subsample_2D_full(
				const RawImage<T>& img, int wout, int hout, double factor = -1.0);
		
		template <typename T>
		static BufferedImage<T> subsample_3D_full(
				const BufferedImage<T>& img, int wout, int hout, int dout, double factor = -1.0);
		
		template <typename T>
		static BufferedImage<T> FourierCrop_fullStack(
				const BufferedImage<T>& img, double factor, int num_threads, bool keepSizeEven);
		
		static gravis::i2Vector getFourierCroppedSize2D(
		        int w, int h, double factor, bool keepSizeEven);
		
		template <typename T>
		static BufferedImage<T> FourierCrop_fftwHalfStack(
				const BufferedImage<T>& img, double factor, int num_threads);
		
		template <typename T>
		static BufferedImage<T> FourierCrop_3D(
				const BufferedImage<T>& img, int w, int h, int d, int num_threads);
		
		
		
		
};

template <typename T>
BufferedImage<T> Resampling::upsampleLinear_2D_half(const RawImage<T>& img, int wout, int hout)
{
	BufferedImage<T> out(wout, hout);
	
	const int win = img.xdim;
	const int hin = img.ydim; 
			
	for (int y = 0; y < hout; y++)
	for (int x = 0; x < wout; x++)
	{
		const double xx = win * x / (double) wout;
		const double yy = hin * y / (double) hout;
		
		const int x0 = (int) xx;
		const int y0 = (int) yy;
		
		int x1 = x0+1;
		int y1 = (y0+1) % hin;
		
		if (x1 >= win) 
		{
			x1 = 2*win - 2 - x1;
			y1 = (hin - y1) % hin;
		}
		
		const double xf = xx - x0;
		const double yf = yy - y0;
		
		double vy0 = (1.0 - xf) * img(x0,y0) + xf * img(x0,y1);
		double vy1 = (1.0 - xf) * img(x1,y0) + xf * img(x1,y1);
		
		out(x,y) = (1.0 - yf) * vy0 + yf * vy1;
	}
	
	return out;
}


template <typename T>
BufferedImage<T> Resampling::upsampleLinear_3D_full(
		const RawImage<T>& img, double factor, int wout, int hout, int dout)
{
	const int w0 = img.xdim;
	const int h0 = img.ydim;
	const int d0 = img.zdim;
	
	if (wout < 0) wout = w0 * factor;
	if (hout < 0) hout = h0 * factor;
	if (dout < 0) dout = d0 * factor;
	
	BufferedImage<T> out(wout, hout, dout);
	
	for (int z = 0; z < dout; z++)
	for (int y = 0; y < hout; y++)
	for (int x = 0; x < wout; x++)
	{
		const double xx = (double) x / factor;
		const double yy = (double) y / factor;
		const double zz = (double) z / factor;
		
		out(x,y,z) = Interpolation::linearXYZ_clip(img, xx, yy, zz);
	}
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::upsampleCubic_2D_full(
		const RawImage<T>& img, int z, double factor, int wout, int hout)
{
	const int win = img.xdim;
	const int hin = img.ydim;
	const int din = img.zdim;
	
	if (z < 0 || z >= din) 
	{
		REPORT_ERROR_STR("Resampling::upsampleCubic_2D_full: Illegal Z index: " << z << " / " << din);
	}
	
	if (wout < 0) wout = win * factor;
	if (hout < 0) hout = hin * factor;
			
	BufferedImage<T> out(wout, hout);
	 
	for (int y = 0; y < hout; y++)
	for (int x = 0; x < wout; x++)
	{
		const double xx = (double) x / factor;
		const double yy = (double) y / factor;
		
		out(x,y) = Interpolation::cubicXY_clip(img, xx, yy, z);
	}
	
	return out;
}


template<class T>
BufferedImage<T> Resampling::upsampleCubic_Stack_full(
		const RawImage<T>& stack, double factor, int wout, int hout)
{
	BufferedImage<T> out(wout, hout, stack.zdim);
	
	for (int f = 0; f < stack.zdim; f++)
	{		
		BufferedImage<T> sliceBig = upsampleCubic_2D_full(stack, f, factor, wout, hout);		
		NewStackHelper::insertSliceZ(sliceBig, out, f);
	}
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::downsampleFilt_2D_full(BufferedImage<T>& img, int wout, int hout, double factor)
{
	if (factor == 1.0) return img;
	
	const double q = factor > 0.0? 1.0/factor : wout / (double) img.xdim;
	
    BufferedImage<T> filt = lowPassFilter_full(img, 0.9*q, q);
    BufferedImage<T> out = subsample_2D_full(filt, wout, hout, factor);
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::downsampleFilt_3D_full(BufferedImage<T>& img, int wout, int hout, int dout, double factor)
{
	if (factor == 1.0) return img;
	
	const double q = factor > 0.0? 1.0/factor : wout / (double) img.xdim;
	
    BufferedImage<T> filt = lowPassFilter_full(img, 0.9*q, q);
    BufferedImage<T> out = subsample_3D_full(filt, wout, hout, dout, factor);
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::downsampleGaussFilt_3D_full(BufferedImage<T>& img, int wout, int hout, int dout, double factor)
{
	if (factor == 1.0) return img;
	
	BufferedImage<T> filt = ImageFilter::Gauss3D(img, factor / 2.0, true);
    BufferedImage<T> out = subsample_3D_full(filt, wout, hout, dout, factor);
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::downsampleFiltStack_2D_full(
		BufferedImage<T>& img, double factor, int num_threads, bool keepSizeEven)
{
	if (factor == 1.0) return img;
	
	int wout = (int)(img.xdim / factor);
	int hout = (int)(img.ydim / factor);
	
	if (keepSizeEven)
	{
		wout -= wout % 2;
		hout -= hout % 2;
	}
	
	const double wavelengthPx = 2.0 * factor;
	const double filterEdgeWidthPx = 10.0;
	
	BufferedImage<T> out(wout, hout, img.zdim);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < img.zdim; f++)
	{
		BufferedImage<T> sliceFilt = ImageFilter::lowpass2D(img, f, wavelengthPx, filterEdgeWidthPx, true);
		NewStackHelper::insertSliceZ(subsample_2D_full(sliceFilt, wout, hout, factor), out, f);
	}
			
	return out;
}

template <typename T>
BufferedImage<T> Resampling::downsampleMax_2D_full(const RawImage<T>& img, int wout, int hout)
{	
	const long int w = img.xdim;
	const long int h = img.ydim;
	
	BufferedImage<T> out = subsample_2D_full(img, wout, hout);
	
	double qx = wout / (double) w;
	double qy = hout / (double) h;

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
		T t = img(x,y);
		
		long int xx = qx * x;
		long int yy = qy * y;
		
		if (xx > 0 && xx < wout && yy > 0 && yy < hout && t > out(xx,yy)) 
		{
			out(xx,yy) = t;
		}
    }
	
	for (long int y = hout-1; y >= 0; y--)
    for (long int x = wout-1; x > 0; x--)
    {
		if (out(x-1,y) > out(x,y))
		{
			out(x,y) = out(x-1,y);
		}
	}
	
	for (long int y = hout-1; y > 0; y--)
    for (long int x = wout-1; x >= 0; x--)
    {
		if (out(x,y-1) > out(x,y))
		{
			out(x,y) = out(x,y-1);
		}
	}
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::lowPassFilter_full(BufferedImage<T>& img, double maxFreq0, double maxFreq1)
{
	BufferedImage<T> out(img.xdim, img.ydim, img.zdim);
			
    BufferedImage<tComplex<T>> imgFreq;
	FFT::FourierTransform(img, imgFreq);
	
	for (size_t z = 0; z < imgFreq.zdim; z++)
	for (size_t y = 0; y < imgFreq.ydim; y++)
	for (size_t x = 0; x < imgFreq.xdim; x++)
    {
        double xi = x/(double)imgFreq.xdim;
        double yi = 2.0*y/(double)imgFreq.ydim;
        double zi = 2.0*z/(double)imgFreq.zdim;

        if (yi > 1.0) yi = 2.0 - yi;
        if (zi > 1.0) zi = 2.0 - zi;

        double r = sqrt(xi*xi + yi*yi + zi*zi);

        if (r > maxFreq1)
        {
			imgFreq(x,y,z) = tComplex<T>(0.0, 0.0);
        }
        else if (r > maxFreq0)
        {
            const double t = (r - maxFreq0)/(maxFreq1 - maxFreq0);
            const T q(0.5 * (cos(PI*t) + 1.0));
			
            imgFreq(x,y,z) *= q;
        }
    }

    FFT::inverseFourierTransform(imgFreq, out);
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::subsample_2D_full(const RawImage<T>& img, int wout, int hout, double factor)
{
	if (factor < 0.0) factor = img.xdim / (double) wout;
			
	BufferedImage<T> out(wout, hout);
	
    for (long int y = 0; y < hout; y++)
    for (long int x = 0; x < wout; x++)
    {
		out(x,y) = Interpolation::linearXY_clip(img, factor * x, factor * y);
    }
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::subsample_3D_full(
		const BufferedImage<T>& img, 
		int wout, int hout, int dout, 
		double factor)
{
	if (factor < 0.0) factor = img.xdim / (double) wout;
			
	BufferedImage<T> out(wout, hout, dout);
	
	for (long int z = 0; z < dout; z++)
	for (long int y = 0; y < hout; y++)
    for (long int x = 0; x < wout; x++)
    {
		out(x,y,z) = Interpolation::linearXYZ_clip(img, factor * x, factor * y, factor * z);
    }
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::FourierCrop_fullStack(
		const BufferedImage<T>& img, double factor, int num_threads, bool keepSizeEven)
{
	if (factor == 1.0) return img;
	
	int w_in = img.xdim;
	int h_in = img.ydim;
	const int fc   = img.zdim;
	
	int w_out = (int)(w_in / factor);
	int h_out = (int)(h_in / factor);

	if (keepSizeEven)
	{
		w_out -= w_out % 2;
		h_out -= h_out % 2;
	}
	
	BufferedImage<T> imgCp;
	
	if ((int)(factor * w_out) != w_in || (int)(factor * h_out) != h_in)
	{
		// Maintain the aspect ratio by first cropping 
		// to the exact multiple of the output size: 
		
		Cutting::extract2D(img, imgCp, 0, 0, factor * w_out, factor * h_out);
		w_in = (int) (w_out * factor);
		h_in = (int) (h_out * factor);
	}
	else
	{
		imgCp = img;
	}
	
	const int wh_in = w_in/2 + 1;
	const int wh_out = w_out/2 + 1;
	
	BufferedImage<tComplex<T>> imgFS(wh_in, h_in, fc);
	
	NewStackHelper::FourierTransformStack(imgCp, imgFS, true, num_threads);
	
	BufferedImage<tComplex<T>> imgFS_cropped(wh_out, h_out, fc);

	const double value_scale = 1.0 / factor;
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	for (int y = 0; y < h_out; y++)
	for (int x = 0; x < wh_out; x++)
	{
		int xx = x;
		int yy = y < h_out/2? y : y - h_out + h_in;
		
		imgFS_cropped(x,y,f) = value_scale * imgFS(xx,yy,f);
	}
	
	BufferedImage<float> imgOut(w_out, h_out, fc);
	NewStackHelper::inverseFourierTransformStack(imgFS_cropped, imgOut, true, num_threads);
	
	return imgOut;
}

template <typename T>
BufferedImage<T> Resampling::FourierCrop_fftwHalfStack(const BufferedImage<T>& img, double factor, int num_threads)
{
	const int wh0 = img.xdim;
	const int w0 = (wh0 - 1) * 2;
	const int h0 = img.ydim;
	
	const int w1 = w0 / factor;
	const int wh1 = w1/2 + 1;
	const int h1 = h0 / factor;
	
	const int fc = img.zdim;
		
	BufferedImage<T> out(wh1, h1, fc);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)		
	for (int y = 0; y < h1; y++)
	for (int x = 0; x < wh1; x++)
	{
		int xx = x;
		int yy = y < h1/2? y : y - h1 + h0;
		
		out(x,y,f) = img(xx,yy,f);
	}
	
	return out;
}

template <typename T>
BufferedImage<T> Resampling::FourierCrop_3D(
		const BufferedImage<T>& img, int w, int h, int d, int num_threads)
{
	const int h0 = img.ydim;
	const int d0 = img.zdim;
	
	const int wh = w/2 + 1;
	
	BufferedImage<T> imgCp = img;
		
	BufferedImage<tComplex<T>> IMG;
	FFT::FourierTransform(imgCp, IMG, FFT::Both);
	
	BufferedImage<tComplex<T>> IMG2(wh, h, d);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int z = 0; z < d; z++)		
	for (int y = 0; y < h; y++)
	for (int x = 0; x < wh; x++)
	{
		int yy = y < h/2? y : y - h + h0;
		int zz = z < d/2? z : z - d + d0;
		
		IMG2(x,y,z) = IMG(x,yy,zz);
	}
	
	BufferedImage<T> out;
	
	FFT::inverseFourierTransform(IMG2, out, FFT::Both);
	
	return out;
}

#endif

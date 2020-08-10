#ifndef BANDPASS_FILTER_H
#define BANDPASS_FILTER_H

#include "raw_image.h"
#include "stack_helper.h"
#include "padding.h"
#include <src/jaz/math/fft.h>

class ImageFilter
{
	public:
		
		template<class T>
		static BufferedImage<T> bandpass(
				const RawImage<T>& img, double freqPx, double widthPx, int numOvertones = 0);
		
		
		template<class T>
		static BufferedImage<T> highpass2D(
				const RawImage<T>& img, int z, double freqPx, double widthPx, bool pad);

		template<class T>
		static BufferedImage<tComplex<T>> highpass3D(
				const RawImage<tComplex<T>>& img, double freqPx, double widthPx);

		template<class T>
		static BufferedImage<tComplex<T>> highpassGauss3D(
				const RawImage<tComplex<T>>& img, double sigmaPx);
		
		template<class T>
		static BufferedImage<T> lowpass2D(
				const RawImage<T>& img, int z, double freqPx, double widthPx, bool pad);

		template<class T>
		static BufferedImage<tComplex<T>> lowpass3D(
				const RawImage<tComplex<T>>& img, double freqPx, double widthPx);
		
		template<class T>
		static BufferedImage<T> Gauss2D(
				const RawImage<T>& img, int z, double sigmaRS, bool pad);
		
		template<class T>
		static BufferedImage<T> ramp(const RawImage<T>& img);

		
		
		template<class T>
		static BufferedImage<T> highpassStack(
				const RawImage<T>& stack, double freqPx, double widthPx, bool pad);
		
		template<class T>
		static BufferedImage<T> highpassStackGaussPadded(
				const RawImage<T>& stack, double sigmaRS, int num_threads = 1);
		
		template<class T>
		static BufferedImage<T> lowpassStack(
				const RawImage<T>& stack, double freqPx, double widthPx, bool pad);
		
		template<class T>
		static BufferedImage<T> GaussStack(
				const RawImage<T>& stack, double sigmaRS, bool pad);
		
		template<class T>
		static BufferedImage<T> Gauss3D(
				const RawImage<T>& img, double sigmaRS, bool pad = false);
		
		template<class T>
		static BufferedImage<T> separableGauss3D(
				const RawImage<T>& img, double sigma);
		
};

template<class T>
BufferedImage<T> ImageFilter::bandpass(const RawImage<T>& img, double freqPx, double widthPx, int numOvertones)
{
	const int w = img.xdim;
	const int wh = w/2 + 1;
	const int h = img.ydim;
	const int d = img.zdim;
	
	BufferedImage<T> imgCp = img;
	BufferedImage<tComplex<T>> imgFS;
	
	FFT::FourierTransform(imgCp, imgFS, FFT::Both);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < wh; x++)
	{
		const double xx = x / (double) w;
		const double yy = (((y + h/2) % h) - h/2) / (double) h;
		const double zz = (((z + d/2) % d) - d/2) / (double) d;
		
		const double f = sqrt(xx * xx + yy * yy + zz * zz);
		
		double freqPxAct;
		
		if (numOvertones > 0 && f > freqPx)
		{
			int nearestMultiple = (int) (f / freqPx + 0.5);
			
			if (nearestMultiple > numOvertones + 1) 
			{
				nearestMultiple = numOvertones + 1;
			}
					
			freqPxAct = nearestMultiple * freqPx;
		}
		else
		{
			freqPxAct = freqPx;
		}
		
		const double df = std::abs((f - freqPxAct) / widthPx);
		
		if (df > 1.0) 
		{
			imgFS(x,y,z) = tComplex<T>(0.0, 0.0);
		}
		else
		{
			imgFS(x,y,z) *= (T) (0.5 * (cos(PI * df) + 1.0));
		}
	}
	
	BufferedImage<T> out;
	
	FFT::inverseFourierTransform(imgFS, out, FFT::Both);
	
	return out;
}

template<class T>
BufferedImage<T> ImageFilter::highpass2D(const RawImage<T>& img, int z, double freqPx, double widthPx, bool pad)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int wh = w/2 + 1;
	
	int padding = pad? (int) (1.0 / freqPx + 0.5) : 0;
	
	const int wp = w + 2 * padding;
	const int hp = h + 2 * padding;
	const int wph = wp/2 + 1;
	
	BufferedImage<T> imgPadded(wp, hp);
	
	for (int y = 0; y < hp; y++)
	for (int x = 0; x < wp; x++)
	{
		int xx = x - padding;
		int yy = y - padding;
		
		if (xx < 0) xx = 0;
		else if (xx >= w) xx = w - 1;
		
		if (yy < 0) yy = 0;
		else if (yy >= h) yy = h - 1;
		
		imgPadded(x,y) = img(xx,yy,z);
	}

	BufferedImage<tComplex<T>> imgFS;
	
	FFT::FourierTransform(imgPadded, imgFS, FFT::Both);
			
	for (int y = 0; y < hp; y++)
	for (int x = 0; x < wph; x++)
	{
		const double xx = x / (double) wp;
		const double yy = (((y + hp/2) % hp) - hp/2) / (double) hp;
		
		const double f = sqrt(xx * xx + yy * yy);
				
		const double df = (f - freqPx) / widthPx + 0.5;
		
		if (df < 0.0) 
		{
			imgFS(x,y) = tComplex<T>(0.0, 0.0);
		}
		else if (df < 1.0)
		{
			imgFS(x,y) *= (T) (1.0 - 0.5 * (cos(PI * df) + 1.0));
		}	
	}
	
	BufferedImage<T> out, outCropped(w,h);
	
	FFT::inverseFourierTransform(imgFS, out, FFT::Both);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		outCropped(x,y) = out(x + padding, y + padding);
	}
	
	return outCropped;
}

template<class T>
BufferedImage<tComplex<T>> ImageFilter::highpass3D(
	const RawImage<tComplex<T>>& img, double freqPx, double widthPx)
{
	const int wh = img.xdim;
	const int w = (wh - 1) * 2;
	const int h = img.ydim;
	const int d = img.zdim;

	BufferedImage<tComplex<T>> out(wh, h, d);

	for (int z = 0; z < d;  z++)
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < wh; x++)
	{
		const double xx = x;
		const double yy = (((y + h/2) % h) - h/2);
		const double zz = (((z + d/2) % d) - d/2);

		const double f = sqrt(xx * xx + yy * yy + zz * zz);

		const double df = (f - freqPx) / widthPx + 0.5;

		if (df < 0.0)
		{
			out(x,y,z) = tComplex<T>(0.0, 0.0);
		}
		else if (df < 1.0)
		{
			out(x,y,z) = (T) (1.0 - 0.5 * (cos(PI * df) + 1.0)) * img(x,y,z);
		}
		else
		{
			out(x,y,z) = img(x,y,z);
		}
	}

	return out;
}

template<class T>
BufferedImage<tComplex<T>> ImageFilter::highpassGauss3D(
	const RawImage<tComplex<T>>& img, double sigmaPx)
{
	const int wh = img.xdim;
	const int w = (wh - 1) * 2;
	const int h = img.ydim;
	const int d = img.zdim;

	BufferedImage<tComplex<T>> out(wh, h, d);

	for (int z = 0; z < d;  z++)
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < wh; x++)
	{
		const double xx = x;
		const double yy = (((y + h/2) % h) - h/2);
		const double zz = (((z + d/2) % d) - d/2);

		const double f = sqrt(xx * xx + yy * yy + zz * zz);
		const double q = f / sigmaPx;

		out(x,y,z) = (1 - exp(-0.5 * q * q)) * img(x,y,z);
	}

	return out;
}

template<class T>
BufferedImage<T> ImageFilter::lowpass2D(const RawImage<T>& img, int z, double freqPx, double widthPx, bool pad)
{
	const int w = img.xdim;
	const int h = img.ydim;
	
	int padding = pad? (int) (1.0 / freqPx + 0.5) : 0;
	
	const int wp = w + 2 * padding;
	const int hp = h + 2 * padding;
	const int wph = wp/2 + 1;
	
	BufferedImage<T> imgPadded(wp, hp);
	
	for (int y = 0; y < hp; y++)
	for (int x = 0; x < wp; x++)
	{
		int xx = x - padding;
		int yy = y - padding;
		
		if (xx < 0) xx = 0;
		else if (xx >= w) xx = w - 1;
		
		if (yy < 0) yy = 0;
		else if (yy >= h) yy = h - 1;
		
		imgPadded(x,y) = img(xx,yy,z);
	}
	
	BufferedImage<tComplex<T>> imgFS;
	
	FFT::FourierTransform(imgPadded, imgFS, FFT::Both);
		
	for (int y = 0; y < hp; y++)
	for (int x = 0; x < wph; x++)
	{
		const double xx = x / (double) wp;
		const double yy = (((y + hp/2) % hp) - hp/2) / (double) hp;
		
		const double f = sqrt(xx * xx + yy * yy);
				
		const double df = (f - freqPx) / widthPx + 0.5;
		
		if (df > 1.0) 
		{
			imgFS(x,y) = tComplex<T>(0.0, 0.0);
		}
		else if (df > 0.0)
		{
			imgFS(x,y) *= (T) (0.5 * (cos(PI * df) + 1.0));
		}	
	}
	
	BufferedImage<T> out, outCropped(w,h);
	
	FFT::inverseFourierTransform(imgFS, out, FFT::Both);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		outCropped(x,y) = out(x + padding, y + padding);
	}
	
	return outCropped;
}

template<class T>
BufferedImage<tComplex<T>> ImageFilter::lowpass3D(
	const RawImage<tComplex<T>>& img, double freqPx, double widthPx)
{
	const int wh = img.xdim;
	const int w = (wh - 1) * 2;
	const int h = img.ydim;
	const int d = img.zdim;

	BufferedImage<tComplex<T>> out(wh, h, d);

	for (int z = 0; z < d;  z++)
	for (int y = 0; y < h;  y++)
	for (int x = 0; x < wh; x++)
	{
		const double xx = x;
		const double yy = (((y + h/2) % h) - h/2);
		const double zz = (((z + d/2) % d) - d/2);

		const double f = sqrt(xx * xx + yy * yy + zz * zz);

		const double df = (f - freqPx) / widthPx + 0.5;

		if (df < 0.0)
		{
			out(x,y,z) = img(x,y,z);
		}
		else if (df < 1.0)
		{
			out(x,y,z) *= (T) (0.5 * (cos(PI * df) + 1.0));
		}
		else
		{
			out(x,y,z) = tComplex<T>(0.0, 0.0);
		}

	}

	return out;
}

template<class T>
BufferedImage<T> ImageFilter::Gauss2D(const RawImage<T>& img, int z, double sigmaRS, bool pad)
{
	const int w = img.xdim;
	const int h = img.ydim;
	
	int padding = pad? (int) (2.0 * sigmaRS + 0.5) : 0;
	
	const double sigmaFSx = 0.5 / (PI * sigmaRS);
	const double sigmaFSy = 0.5 / (PI * sigmaRS);
		
	const int wp = w + 2 * padding;
	const int hp = h + 2 * padding;
	const int wph = wp/2 + 1;
	
	BufferedImage<T> imgPadded(wp, hp);
	
	for (int y = 0; y < hp; y++)
	for (int x = 0; x < wp; x++)
	{
		int xx = x - padding;
		int yy = y - padding;
		
		if (xx < 0) xx = 0;
		else if (xx >= w) xx = w - 1;
		
		if (yy < 0) yy = 0;
		else if (yy >= h) yy = h - 1;
		
		imgPadded(x,y) = img(xx,yy,z);
	}
	
	BufferedImage<tComplex<T>> imgFS;
	
	FFT::FourierTransform(imgPadded, imgFS, FFT::Both);
		
	for (int y = 0; y < hp; y++)
	for (int x = 0; x < wph; x++)
	{
		const double xx = x / (double) wp;
		const double yy = (((y + hp/2) % hp) - hp/2) / (double) hp;
		
		const double d = 
			  xx * xx / (sigmaFSx * sigmaFSx) 
			+ yy * yy / (sigmaFSy * sigmaFSy);
				
		imgFS(x,y) *= exp(-d / 2.0);
	}
	
	BufferedImage<T> out, outCropped(w,h);
	
	FFT::inverseFourierTransform(imgFS, out, FFT::Both);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		outCropped(x,y) = out(x + padding, y + padding);
	}
	
	return outCropped;
}

template<class T>
BufferedImage<T> ImageFilter::ramp(const RawImage<T>& img)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	const int wh = w/2 + 1;
	
	BufferedImage<T> imgCp = img;
	BufferedImage<tComplex<T>> imgFS;
	
	FFT::FourierTransform(imgCp, imgFS, FFT::Both);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < wh; x++)
	{
		double xx = x;
		double yy = y < h/2? y : y - h;
		double zz = z < d/2? z : z - d;
		
		imgFS(x,y,z) *= sqrt(xx*xx + yy*yy + zz*zz) / wh;
	}
	
	FFT::inverseFourierTransform(imgFS, imgCp, FFT::Both);
	
	return imgCp;
}

template<class T>
BufferedImage<T> ImageFilter::highpassStack(
		const RawImage<T>& stack, double freqPx, double widthPx, bool pad)
{
	BufferedImage<T> out(stack);
	
	const int fc = stack.zdim;
	
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<T> sliceFilt = highpass2D(stack, f, freqPx, widthPx, pad);
		NewStackHelper::insertSliceZ(sliceFilt, out, f);
	}
	
	return out;
}

template<class T>
BufferedImage<T> ImageFilter::highpassStackGaussPadded(
		const RawImage<T>& stack, double sigmaRS, int num_threads)
{
	BufferedImage<T> out(stack.xdim, stack.ydim, stack.zdim);
	
	const int fc = stack.zdim;
	
	BufferedImage<T> mask(stack.xdim, stack.ydim, 1);
	mask.fill(T(1));
	mask = Padding::padCenter2D_full(mask, 2*sigmaRS);
	mask = Gauss2D(mask, 0, sigmaRS, false);
	mask = Padding::unpadCenter2D_full(mask, 2*sigmaRS);

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{	
		BufferedImage<T> slice0 = NewStackHelper::extractSliceZ(stack, f);
		BufferedImage<T> slice = Padding::padCenter2D_full(slice0, 2*sigmaRS);
		
		slice = Gauss2D(slice, 0, sigmaRS, false);
		slice = Padding::unpadCenter2D_full(slice, 2*sigmaRS);
		slice /= mask;
		slice = slice0 - slice;
		
		out.getSliceRef(f).copyFrom(slice);
	}
	
	return out;
}


template<class T>
BufferedImage<T> ImageFilter::lowpassStack(
		const RawImage<T>& stack, double freqPx, double widthPx, bool pad)
{
	BufferedImage<T> out(stack);
	
	const int fc = stack.zdim;
	
	for (int f = 0; f < fc; f++)
	{		
		BufferedImage<T> sliceFilt = lowpass2D(stack, f, freqPx, widthPx, pad);		
		NewStackHelper::insertSliceZ(sliceFilt, out, f);
	}
	
	return out;
}

template<class T>
BufferedImage<T> ImageFilter::GaussStack(
		const RawImage<T>& stack, double sigmaRS, bool pad)
{
	BufferedImage<T> out(stack);
	
	const int fc = stack.zdim;
	
	for (int f = 0; f < fc; f++)
	{		
		BufferedImage<T> sliceFilt = Gauss2D(stack, f, sigmaRS, pad);		
		NewStackHelper::insertSliceZ(sliceFilt, out, f);
	}
	
	return out;
}

template<class T>
BufferedImage<T> ImageFilter::Gauss3D(const RawImage<T>& img, double sigmaRS, bool pad)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	int padding = pad? (int) (2.0 * sigmaRS + 0.5) : 0;
	
	const double sigmaFSx = 0.5 / (PI * sigmaRS);
		
	const int wp = w + 2 * padding;
	const int hp = h + 2 * padding;
	const int dp = d + 2 * padding;
	const int wph = wp/2 + 1;
	
	BufferedImage<T> imgPadded(wp, hp, dp);
	
	for (int z = 0; z < dp; z++)
	for (int y = 0; y < hp; y++)
	for (int x = 0; x < wp; x++)
	{
		int xx = x - padding;
		int yy = y - padding;
		int zz = z - padding;
		
		if (xx < 0) xx = 0;
		else if (xx >= w) xx = w - 1;
		
		if (yy < 0) yy = 0;
		else if (yy >= h) yy = h - 1;
		
		if (zz < 0) zz = 0;
		else if (zz >= d) zz = d - 1;
		
		if (xx >= 0 && xx < w &&
			yy >= 0 && yy < h &&
			zz >= 0 && zz < d) 
		{
			imgPadded(x,y,z) = img(xx,yy,zz);
		}
		else
		{
			imgPadded(x,y,z) = 0;
		}
	}
	
	BufferedImage<tComplex<T>> imgFS;
	
	FFT::FourierTransform(imgPadded, imgFS, FFT::Both);
	
	for (int z = 0; z < dp; z++)
	for (int y = 0; y < hp; y++)
	for (int x = 0; x < wph; x++)
	{
		const double xx = x / (double) wp;
		const double yy = (((y + hp/2) % hp) - hp/2) / (double) hp;
		const double zz = (((z + dp/2) % dp) - dp/2) / (double) dp;
		
		const double dd = (xx * xx + yy * yy + zz * zz) / (sigmaFSx * sigmaFSx);
				
		imgFS(x,y,z) *= exp(-dd / 2.0);
	}
	
	BufferedImage<T> out, outCropped(w,h,d);
	
	FFT::inverseFourierTransform(imgFS, out, FFT::Both);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		outCropped(x,y,z) = out(x + padding, y + padding, z + padding);
	}
	
	return outCropped;
}

template<typename T>
BufferedImage<T> ImageFilter::separableGauss3D(const RawImage<T>& img, double sigma)
{
    const int k = (int)(3 * sigma + 0.5);
	
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;

    BufferedImage<T> out(img);

    std::vector<double> kernel(2 * k + 1);
    const double s2 = sigma * sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-i*i/s2);
    }

    BufferedImage<T> temp(w, h, d);

    for (size_t z = 0; z < d; z++)
    for (size_t y = 0; y < h; y++)
    for (size_t x = 0; x < w; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int xx = x + i;
            if (xx < 0 || xx >= w) continue;

            v += kernel[i+k] * img(xx,y,z);
            m += kernel[i+k];
        }

        out(x,y,z) = v/m;
    }

    for (size_t z = 0; z < d; z++)
    for (size_t y = 0; y < h; y++)
    for (size_t x = 0; x < w; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = y + i;
            if (yy < 0 || yy >= h) continue;

            v += kernel[i+k] * out(x,yy,z);
            m += kernel[i+k];
        }

        temp(x,y,z) = v/m;
    }

    for (size_t z = 0; z < d; z++)
    for (size_t y = 0; y < h; y++)
    for (size_t x = 0; x < w; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int zz = z + i;
            if (zz < 0 || zz >= d) continue;

            v += kernel[i+k] * temp(x,y,zz);
            m += kernel[i+k];
        }

        out(x,y,z) = v/m;
    }
	
	return out;
}

#endif


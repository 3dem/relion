/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * C-CINA Basel
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef RECONSTRUCTION_H
#define RECONSTRUCTION_H

#include <string>
#include <omp.h>
#include <src/jaz/gravis/t3Vector.h>

#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/interpolation.h>
#include "extraction.h"
#include "tomo_stack.h"


class Reconstruction
{
    public:
		
		template <typename T>
		static void griddingCorrect3D(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<T>& psfImgFS,
				BufferedImage<T>& out,
				bool center, 
				int num_threads);
		
		template <typename T>
		static void ctfCorrect3D(
				BufferedImage<T>& img,
				BufferedImage<T>& ctfImgFS,
				BufferedImage<T>& out,
				double WienerOffset,
				int num_threads);
		
		template <typename T>
		static void taper(
				BufferedImage<T>& img,
				double dist,
				bool do_center,
				int num_threads);
		
		template <typename T>
		static void GaussEnvelope(
				BufferedImage<T>& img,
				double sigma,
				bool do_center,
				int num_threads);		
		
		
		template <typename T>
		static void correct3D_RS(
				BufferedImage<T>& dataImgRS,
				BufferedImage<T>& psfImgRS,
				BufferedImage<T>& out,
				double WienerOffset,
				int num_threads);
		
		// used after a forward projection:
		template <typename T>
		static void correctStack(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<tComplex<T>>& psfImgFS,
				RawImage<T>& out,
				bool center, 
				int num_threads);
		
		template <typename T>
		static void correctStackFS(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<tComplex<T>>& psfImgFS,
				BufferedImage<tComplex<T>>& out,
				bool center, 
				int num_threads);
		
		// used after a cross projection:
		template <typename T>
		static void correctStack(
				BufferedImage<tComplex<T>>& dataImgFS,
				BufferedImage<T>& psfImgFS,
				BufferedImage<T>& ctfImgFS,
				RawImage<T>& out,
				double WienerOffset,
				bool center, 
				int num_threads);

};

template <typename T>
void Reconstruction :: griddingCorrect3D(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<T>& psfImgFS,
		BufferedImage<T>& out,
		bool center, 
		int num_threads)
{
	const long int wh = dataImgFS.xdim;
	const long int w = 2 * (wh - 1);
	const long int h = dataImgFS.ydim;
	const long int d = dataImgFS.zdim;
	
	if (center)
	{
		#pragma omp parallel for num_threads(num_threads)
		for (long int z = 0; z < d;  z++)
		for (long int y = 0; y < h;  y++)
		for (long int x = 0; x < wh; x++)
		{
			dataImgFS(x,y,z) *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		}
	}
	
	// gridding correction
	
	BufferedImage<T> dataImg(w,h,d);	
	FFT::inverseFourierTransform(dataImgFS, dataImg, FFT::Both);

	BufferedImage<tComplex<T>> psfImgFS_complex(wh,h,d);
	
	#pragma omp parallel for num_threads(num_threads)	
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		psfImgFS_complex(x,y,z) = psfImgFS(x,y,z);
		
		if (center) 
		{
			psfImgFS_complex(x,y,z) *= (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		}
	}
	
	BufferedImage<T> dataEnv(w,h,d);	
	FFT::inverseFourierTransform(psfImgFS_complex, dataEnv, FFT::Both);
		
	const double eps = 1e-10;
	
	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		if (dataEnv(x,y,z) > eps)
		{
			out(x,y,z) = dataImg(x,y,z) / dataEnv(x,y,z);
		}
	}
}
	
template <typename T>
void Reconstruction :: ctfCorrect3D(
		BufferedImage<T>& img,
		BufferedImage<T>& ctfImgFS,
		BufferedImage<T>& out,
		double WienerOffset,
		int num_threads)
{
	const long int w = img.xdim;
	const long int wh = w/2 + 1;
	const long int h = img.ydim;
	const long int d = img.zdim;
	
	// CTF correction	

	BufferedImage<tComplex<T>> dataImgCorrFS(wh,h,d);	
	FFT::FourierTransform(img, dataImgCorrFS, FFT::Both);
	
	#pragma omp parallel for num_threads(num_threads)
	for (long int z = 0; z < d;  z++)
	for (long int y = 0; y < h;  y++)
	for (long int x = 0; x < wh; x++)
	{
		dataImgCorrFS(x,y,z) /= ctfImgFS(x,y,z) + WienerOffset;
	}
	
	FFT::inverseFourierTransform(dataImgCorrFS, out, FFT::Both);
}

template <typename T>
void Reconstruction :: taper(
		BufferedImage<T>& img,
		double dist,
		bool do_center,
		int num_threads)
{
	double innerAvg = 0.0;
	double innerWgh = 0.0;
	
	const int s = img.xdim;
	const double m = s/2.0;
	
	if (do_center)
	{
		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			const double dx = x - m;
			const double dy = y - m;
			const double dz = z - m;
			
			const double r = sqrt(dx*dx + dy*dy + dz*dz);
			
			if (r < m)
			{
				double c = r < m - dist? 1.0 : 0.5 - 0.5 * cos(PI * (m - r) / dist);
				
				innerAvg += c * img(x,y,z);
				innerWgh += c;
			}
		}
		
		innerAvg /= innerWgh;
	}
	else
	{
		innerAvg = 0.0;
	}
	
	
	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double dx = x - m;
		const double dy = y - m;
		const double dz = z - m;
		
		const double r = sqrt(dx*dx + dy*dy + dz*dz);
		
		if (r < m)
		{
			double c = r < m - dist? 1.0 : 0.5 - 0.5 * cos(PI * (m - r) / dist);
			
			img(x,y,z) = c * (img(x,y,z) - innerAvg);
		}
		else
		{
			img(x,y,z) = 0;
		}
	}
}


template <typename T>
void Reconstruction :: GaussEnvelope(
		BufferedImage<T>& img,
		double sigma,
		bool do_center,
		int num_threads)
{
	if (sigma < 0.0) return;
	
	double innerAvg = 0.0;
	double innerWgh = 0.0;
	
	const int s = img.xdim;
	const double m = s/2.0;
	const double beta = -1.0 / (2.0 * sigma * sigma);
	
	if (do_center)
	{
		for (int z = 0; z < s; z++)
		for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			const double dx = x - m;
			const double dy = y - m;
			const double dz = z - m;
			
			const double c = exp(beta * (dx*dx + dy*dy + dz*dz));
			
			innerAvg += c * img(x,y,z);
			innerWgh += c;
		}
		
		innerAvg /= innerWgh;
	}
	else
	{
		innerAvg = 0.0;
	}
	
	
	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double dx = x - m;
		const double dy = y - m;
		const double dz = z - m;
		
		const double c = exp(beta * (dx*dx + dy*dy + dz*dz));
		
		img(x,y,z) = c * (img(x,y,z) - innerAvg);
	}
}

template <typename T>
void Reconstruction :: correct3D_RS(
		BufferedImage<T>& dataImgRS,
		BufferedImage<T>& psfImgRS,
		BufferedImage<T>& out,
		double WienerOffset,
		int num_threads)
{
	BufferedImage<fComplex> dataFS, psfFS;
	
	FFT::FourierTransform(dataImgRS, dataFS, FFT::Both);
	FFT::FourierTransform(psfImgRS, psfFS, FFT::Both);
	
	for (int z = 0; z < dataFS.zdim; z++)
	for (int y = 0; y < dataFS.ydim; y++)
	for (int x = 0; x < dataFS.xdim; x++)
	{
		const double wa = psfFS(x,y,z).abs() + WienerOffset;
		
		if (wa > 0.0) dataFS(x,y,z) /= wa;
	}
	
	FFT::inverseFourierTransform(dataFS, out, FFT::Both);
}

template <typename T>
void Reconstruction :: correctStack(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<tComplex<T>>& psfImgFS,
		RawImage<T>& out,
		bool center, 
		int num_threads)
{
	const int wh = dataImgFS.xdim;
	const int w  = 2 * (wh - 1);
	const int h  = dataImgFS.ydim;
	const int fc = dataImgFS.zdim;
	
	BufferedImage<T> data(w,h,fc), psf(w,h,fc);
	
	StackHelper::inverseFourierTransformStack(dataImgFS, data, center);	
	StackHelper::inverseFourierTransformStack(psfImgFS, psf, center);
		
	#pragma omp parallel for num_threads(num_threads)			
	for (long int f = 0; f < fc; f++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < w; x++)
	{
		out(x,y,f) = data(x,y,f) / psf(x,y,f);
	}
}

template <typename T>
void Reconstruction :: correctStackFS(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<tComplex<T>>& psfImgFS,
		BufferedImage<tComplex<T>>& out,
		bool center, 
		int num_threads)
{
	const int wh = dataImgFS.xdim;
	const int w  = 2 * (wh - 1);
	const int h  = dataImgFS.ydim;
	const int fc = dataImgFS.zdim;
	
	BufferedImage<T> temp(w,h,fc);
	
	correctStack(dataImgFS, psfImgFS, temp, center, num_threads);
	
	out = StackHelper::FourierTransformStack(temp, center, num_threads);
}


template <typename T>
void Reconstruction :: correctStack(
		BufferedImage<tComplex<T>>& dataImgFS,
		BufferedImage<T>& psfImgFS,
		BufferedImage<T>& ctfImgFS,
		RawImage<T>& out,
		double WienerOffset,
		bool center, 
		int num_threads)
{
	const long int wh = dataImgFS.xdim;
	const long int w = 2 * (wh - 1);
	const long int h = dataImgFS.ydim;
	const long int fc = dataImgFS.zdim;
	
	// deconvolve by Fourier-space 2D-PSF:
	
	std::vector<BufferedImage<T>> 
			tempDataRS(num_threads, BufferedImage<T>(w,h)),
			tempEnvRS(num_threads, BufferedImage<T>(w,h));
	
	std::vector<BufferedImage<tComplex<T>>> 
			tempDataFS(num_threads, BufferedImage<tComplex<T>>(wh,h)),
			tempEnvFS(num_threads, BufferedImage<tComplex<T>>(wh,h));
	
	#pragma omp parallel for num_threads(num_threads)			
	for (long int f = 0; f < fc;  f++)
	{
		const int t = omp_get_thread_num();
		
		if (center)
		{
			for (long int y = 0; y < h;  y++)
			for (long int x = 0; x < wh; x++)
			{
				tempDataFS[t](x,y) = dataImgFS(x,y,f) * (1 - 2*(x%2)) * (1 - 2*(y%2));
				tempEnvFS[t](x,y) = psfImgFS(x,y,f) * (1 - 2*(x%2)) * (1 - 2*(y%2));
			}
		}
		else
		{
			for (long int y = 0; y < h;  y++)
			for (long int x = 0; x < wh; x++)
			{
				tempDataFS[t](x,y) = dataImgFS(x,y,f);
				tempEnvFS[t](x,y) = psfImgFS(x,y,f);
			}
		}
		
		FFT::inverseFourierTransform(tempDataFS[t], tempDataRS[t], FFT::Both);
		FFT::inverseFourierTransform(tempEnvFS[t], tempEnvRS[t], FFT::Both);
		
		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < w; x++)
		{
			out(x,y,f) = tempDataRS[t](x,y) / tempEnvRS[t](x,y);
		}
	}
	
	// divide by ctf in Fourier space
	
	#pragma omp parallel for num_threads(num_threads)			
	for (long int f = 0; f < fc;  f++)
	{
		const int t = omp_get_thread_num();
		
		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < w; x++)
		{
			tempDataRS[t](x,y) = out(x,y,f);
		}
		
		FFT::FourierTransform(tempDataRS[t], tempDataFS[t], FFT::Both);
		
		for (long int y = 0; y < h;  y++)
		for (long int x = 0; x < wh; x++)
		{
			tempDataFS[t](x,y) /= (ctfImgFS(x,y,f) + WienerOffset);
		}
		
		FFT::inverseFourierTransform(tempDataFS[t], tempDataRS[t], FFT::Both);
		
		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < w; x++)
		{
			out(x,y,f) = tempDataRS[t](x,y);
		}
	}

}

#endif

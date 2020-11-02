/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
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

#ifndef TAPERING_H
#define TAPERING_H

#include <src/jaz/image/raw_image.h>
#include <string>
#include <src/jaz/gravis/t3Vector.h>
	
class Tapering
{
	public:

		inline static double getRadialTaperWeight(
			double r,
			double r0,
			double r1);

		inline static double getRadialTaperWeight2D(
			double x,
			double y,
			int w,
			int h,
			double r0,
			double r1);

		inline static double getTaperWeight2D(
			double x,
			double y,
			int w,
			int h,
			double falloff,
			double dist);
		
		inline static double getTaperWeight3D(
			double x, 
			double y,
			double z, 
			int w, 
			int h,
			int d,
			double r);

		template <typename T>
		static void taper(
			RawImage<T>& img,
			double d,
			int num_threads = 1,
			bool exactMean = false);

		template <typename T>
		static void taper2D(
			RawImage<T>& img,
			double d,
			int num_threads = 1,
			bool exactMean = false);

		template <typename T>
		static void taperCircularly2D(
			RawImage<T>& img,
			double r0, double r1,
			int num_threads = 1,
			bool exactMean = false);
};

inline double Tapering::getRadialTaperWeight(double r, double r0, double r1)
{
	if (r < r0) return 1.0;
	else if (r > r1) return 0.0;
	else return (cos(PI * (r - r0) / (r1 - r0)) + 1.0)/2.0;
}

inline double Tapering::getRadialTaperWeight2D(
		double x,
		double y,
		int w,
		int h,
		double r0,
		double r1)
{
	const double dx = x - w/2;
	const double dy = y - h/2;
	const double d = sqrt(dx*dx + dy*dy);
	const double d0 = w/2;

	return getRadialTaperWeight(d, r0, r1);
}

inline double Tapering::getTaperWeight2D(
		double x, double y,
		int w, int h,
		double falloff, double dist)
{
	double wx(1.0), wy(1.0);

	if (x < dist)
	{
		return 0.0;
	}
	else if (x < dist + falloff)
	{
		wx *= (1.0 - cos(PI * (x + 1 - dist) / falloff)) / 2.0;
	}

	if (x >= w - dist)
	{
		return 0.0;
	}
	if (x >= w - dist - falloff)
	{
		wx *= (1.0 - cos(PI * (w - dist - x) / falloff)) / 2.0;
	}

	if (y < dist)
	{
		return 0.0;
	}
	else if (y < dist + falloff)
	{
		wy *= (1.0 - cos(PI * (y + 1 - dist) / falloff)) / 2.0;
	}

	if (y >= h - dist)
	{
		return 0.0;
	}
	if (y >= h - dist - falloff)
	{
		wy *= (1.0 - cos(PI * (h - dist - y) / falloff)) / 2.0;
	}

	return wx * wy;
}

inline double Tapering::getTaperWeight3D(
		double x, double y, double z, 
		int w, int h, int d, 
		double r)
{
	double wx(1.0), wy(1.0), wz(1.0);

	if (x < r) 
	{
		wx *= (1.0 - cos(PI * (x+1) / r))/2.0;
	}

	if (x >= w - r) 
	{
		wx *= (1.0 - cos(PI * (w - x) / r))/2.0;
	}

	if (y < r) 
	{
		wy *= (1.0 - cos(PI * (y+1) / r))/2.0;
	}

	if (y >= h - r) 
	{
		wy *= (1.0 - cos(PI * (h - y) / r))/2.0;
	}
	
	if (z < r) 
	{
		wz *= (1.0 - cos(PI * (z+1) / r))/2.0;
	}
	
	if (z >= d - r) 
	{
		wz *= (1.0 - cos(PI * (d - z) / r))/2.0;
	}

	return wx * wy * wz;	
}

template <typename T>
void Tapering::taper(RawImage<T>& img, double r, int num_threads, bool exactMean)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;

	if (exactMean)
	{
		double dSum = 0.0;
		double wgSum = 0.0;

		for (size_t z = 0; z < d; z++)
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			const double t = getTaperWeight3D(x,y,z,w,h,d,r);

			dSum += t * img(x,y,z);
			wgSum += t;
		}

		const double wgMean = dSum / wgSum;

		#pragma omp parallel for num_threads(num_threads)
		for (size_t z = 0; z < d; z++)
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			const double t = getTaperWeight3D(x,y,z,w,h,d,r);

			img(x,y,z) = t * img(x,y,z) + (1.0 - t) * wgMean;
		}
	}
	else
	{
		#pragma omp parallel for num_threads(num_threads)
		for (size_t z = 0; z < d; z++)
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			img(x,y,z) *= getTaperWeight3D(x,y,z,w,h,d,r);
		}
	}
}

template <typename T>
void Tapering::taper2D(RawImage<T>& img, double r, int num_threads, bool exactMean)
{
	const int w = img.xdim;
	const int h = img.ydim;

	if (exactMean)
	{
		double dSum = 0.0;
		double wgSum = 0.0;

		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			const double t = getTaperWeight2D(x,y,w,h,r,0);

			dSum += t * img(x,y);
			wgSum += t;
		}

		const double wgMean = dSum / wgSum;

		#pragma omp parallel for num_threads(num_threads)
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			const double t = getTaperWeight2D(x,y,w,h,r,0);

			img(x,y) = t * img(x,y) + (1.0 - t) * wgMean;
		}
	}
	else
	{
		#pragma omp parallel for num_threads(num_threads)
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			img(x,y) *= getTaperWeight2D(x,y,w,h,r,0);
		}
	}
}

template <typename T>
void Tapering::taperCircularly2D(RawImage<T>& img, double r0, double r1, int num_threads, bool exactMean)
{
	const int w = img.xdim;
	const int h = img.ydim;

	if (exactMean)
	{
		double dSum = 0.0;
		double wgSum = 0.0;

		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			const double t = getRadialTaperWeight2D(x,y,w,h,r0,r1);

			dSum += t * img(x,y);
			wgSum += t;
		}

		const double wgMean = dSum / wgSum;

		#pragma omp parallel for num_threads(num_threads)
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			const double t = getRadialTaperWeight2D(x,y,w,h,r0,r1);

			img(x,y) = t * img(x,y) + (1.0 - t) * wgMean;
		}
	}
	else
	{
		#pragma omp parallel for num_threads(num_threads)
		for (size_t y = 0; y < h; y++)
		for (size_t x = 0; x < w; x++)
		{
			img(x,y) *= getRadialTaperWeight2D(x,y,w,h,r0,r1);
		}
	}
}

#endif

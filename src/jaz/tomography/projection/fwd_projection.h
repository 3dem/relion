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

#ifndef FORWARD_PROJECTION_H
#define FORWARD_PROJECTION_H

#include <string>
#include <omp.h>
#include <src/jaz/gravis/t3Vector.h>

#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/interpolation.h>
#include <iostream>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/tomo_stack.h>


class ForwardProjection
{
	public:

		template <typename T>
		static void forwardProject(
				const RawImage<tComplex<T>>& reference,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<tComplex<T>>& dest,
				int num_threads = 1);

		template <typename T>
		static void forwardProjectWithinRange(
				const int* xRanges,
				const RawImage<tComplex<T>>& reference,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<tComplex<T>>& dest,
				int num_threads = 1);

		template <typename T>
		static void forwardProject3DGradient(
				const RawImage<tComplex<T>>& reference,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<gravis::t3Vector<tComplex<T>>>& dest,
				int num_threads = 1);

		template <typename T>
		static void forwardProject2DGradient(
				const RawImage<tComplex<T>>& reference,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<gravis::t2Vector<tComplex<T>>>& dest,
				int num_threads = 1);

		template <typename T>
		static void forwardProject_withPSF(
				const RawImage<tComplex<T>>& reference,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<tComplex<T>>& destData,
				RawImage<tComplex<T>>& destPsf,
				int num_threads = 1);
};


template <typename T>
void ForwardProjection::forwardProject(
		const RawImage<tComplex<T>>& reference,
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<tComplex<T>>& dest,
		int num_threads)
{
	const int wh2 = dest.xdim;
	const int h2 = dest.ydim;
	const int fc = dest.zdim;

	std::vector<gravis::d3Matrix> projTransp(fc);

	for (int f = 0; f < fc; f++)
	{
		projTransp[f] = gravis::d3Matrix(
				proj[f](0,0), proj[f](1,0), proj[f](2,0),
				proj[f](0,1), proj[f](1,1), proj[f](2,1),
				proj[f](0,2), proj[f](1,2), proj[f](2,2) );
	}

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		for (long int yi = 0; yi < h2;  yi++)
		for (long int xi = 0; xi < wh2; xi++)
		{
			const double xd = xi;
			const double yd = yi >= h2/2? yi - h2 : yi;

			const gravis::d3Vector pi(xd, yd, 0.0);
			gravis::d3Vector pw = projTransp[f] * pi;

			dest(xi, yi, f) = Interpolation::linearXYZ_FftwHalf_complex(
						reference, pw.x, pw.y, pw.z);
		}
	}
}

template <typename T>
void ForwardProjection::forwardProjectWithinRange(
		const int* xRanges,
		const RawImage<tComplex<T>>& reference,
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<tComplex<T>>& dest,
		int num_threads)
{
	const int wh2 = dest.xdim;
	const int h2 = dest.ydim;
	const int fc = dest.zdim;

	std::vector<gravis::d3Matrix> projTransp(fc);

	for (int f = 0; f < fc; f++)
	{
		projTransp[f] = gravis::d3Matrix(
				proj[f](0,0), proj[f](1,0), proj[f](2,0),
				proj[f](0,1), proj[f](1,1), proj[f](2,1),
				proj[f](0,2), proj[f](1,2), proj[f](2,2) );
	}

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		for (long int yi = 0; yi < h2;  yi++)
		{
			for (long int xi = 0; xi < xRanges[yi]; xi++)
			{
				const double xd = xi;
				const double yd = yi >= h2/2? yi - h2 : yi;

				const gravis::d3Vector pi(xd, yd, 0.0);
				gravis::d3Vector pw = projTransp[f] * pi;

				dest(xi, yi, f) = Interpolation::linearXYZ_FftwHalf_complex(
							reference, pw.x, pw.y, pw.z);
			}

			for (long int xi = xRanges[yi]; xi < wh2; xi++)
			{
				dest(xi, yi, f) = fComplex(0.f, 0.f);
			}
		}
	}
}

template <typename T>
void ForwardProjection::forwardProject3DGradient(
		const RawImage<tComplex<T>>& reference,
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<gravis::t3Vector<tComplex<T>>>& dest,
		int num_threads)
{
	const int wh2 = dest.xdim;
	const int h2 = dest.ydim;
	const int fc = dest.zdim;

	std::vector<gravis::d3Matrix> projTransp(fc);

	for (int f = 0; f < fc; f++)
	{
		projTransp[f] = gravis::d3Matrix(
				proj[f](0,0), proj[f](1,0), proj[f](2,0),
				proj[f](0,1), proj[f](1,1), proj[f](2,1),
				proj[f](0,2), proj[f](1,2), proj[f](2,2) );
	}

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		for (long int yi = 0; yi < h2;  yi++)
		for (long int xi = 0; xi < wh2; xi++)
		{
			const double xd = xi;
			const double yd = yi >= h2/2? yi - h2 : yi;

			const gravis::d3Vector pi(xd, yd, 0.0);
			gravis::d3Vector pw = projTransp[f] * pi;

			dest(xi, yi, f) = Interpolation::linearXYZGradient_FftwHalf_complex(
						reference, pw.x, pw.y, pw.z);
		}
	}
}

template <typename T>
void ForwardProjection::forwardProject2DGradient(
		const RawImage<tComplex<T>>& reference,
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<gravis::t2Vector<tComplex<T>>>& dest,
		int num_threads)
{
	const int w = dest.xdim;
	const int h = dest.ydim;
	const int fc = dest.zdim;

	BufferedImage<gravis::t3Vector<tComplex<T>>> grad3D(w,h,fc);

	forwardProject3DGradient(reference, proj, grad3D, num_threads);

	for (int f = 0; f < fc; f++)
	{
		for (int yi = 0; yi < h; yi++)
		for (int xi = 0; xi < w; xi++)
		{
			const gravis::t3Vector<tComplex<T>> g = grad3D(xi,yi,f);

			dest(xi,yi).x = proj[f](0,0) * g.x + proj[f](0,1) * g.y + proj[f](0,2) * g.z;
			dest(xi,yi).y = proj[f](1,0) * g.x + proj[f](1,1) * g.y + proj[f](1,2) * g.z;
		}
	}
}

template <typename T>
void ForwardProjection::forwardProject_withPSF(
		const RawImage<tComplex<T>>& reference,
		const std::vector<gravis::d4Matrix>& proj,
		RawImage<tComplex<T>>& destData,
		RawImage<tComplex<T>>& destPsf,
		int num_threads)
{
	const int wh2 = destData.xdim;
	const int h2 = destData.ydim;
	const int fc = destData.zdim;

	std::vector<gravis::d3Matrix> projTransp(fc);

	for (int f = 0; f < fc; f++)
	{
		projTransp[f] = gravis::d3Matrix(
				proj[f](0,0), proj[f](1,0), proj[f](2,0),
				proj[f](0,1), proj[f](1,1), proj[f](2,1),
				proj[f](0,2), proj[f](1,2), proj[f](2,2) );
	}

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		for (long int yi = 0; yi < h2;  yi++)
		for (long int xi = 0; xi < wh2; xi++)
		{
			const double xd = xi;
			const double yd = yi >= h2/2? yi - h2 : yi;

			const gravis::d3Vector pi(xd, yd, 0.0);
			gravis::d3Vector pw = projTransp[f] * pi;

			destData(xi, yi, f) = Interpolation::linearXYZ_FftwHalf_complex(
						reference, pw.x, pw.y, pw.z);

			const double ax = std::abs(pw.x);
			const double ay = std::abs(pw.y);
			const double az = std::abs(pw.z);

			if (ax < 1.0 && ay < 1.0 && az < 1.0)
			{
				destPsf(xi, yi, f) = (1.0 - ax) * (1.0 - ay) * (1.0 - az);
			}
			else
			{
				destPsf(xi, yi, f) = 0.0;
			}
		}
	}
}

#endif

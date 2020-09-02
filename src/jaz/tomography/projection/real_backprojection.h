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

#ifndef REAL_BACKPROJECTION_H
#define REAL_BACKPROJECTION_H

#include <string>
#include <omp.h>
#include <src/jaz/gravis/t3Vector.h>

#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/interpolation.h>
#include <iostream>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/image/tapering.h>


class RealSpaceBackprojection
{
	public:
		
		enum InterpolationType {Linear, Cubic};
		
		template <typename SrcType, typename DestType>
		static void backproject(
			const RawImage<SrcType>& stack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<DestType>& dest,
			int num_threads = 1,
			gravis::d3Vector origin = gravis::d3Vector(0.0, 0.0, 0.0),
			double spacing = 1.0,
			InterpolationType interpolation = Linear,
			double taper = 0.0);
		
		template <typename SrcType, typename DestType>
		static void backprojectCoverage(
			const RawImage<SrcType>& stack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<DestType>& dest,
			int num_threads = 1,
			gravis::d3Vector origin = gravis::d3Vector(0.0, 0.0, 0.0),
			double spacing = 1.0,
			InterpolationType interpolation = Linear,
			double taper = 0.0);
		
		template <typename SrcType, typename DestType>
		static void backprojectPsf(
			const RawImage<SrcType>& stack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<DestType>& dest,
			int num_threads = 1,
			gravis::d3Vector origin = gravis::d3Vector(0.0, 0.0, 0.0),
			double spacing = 1.0,
			InterpolationType interpolation = Linear);
				
		template <typename T>
		static void backprojectRaw(
			const TomoStack<T>& stack,
			RawImage<T>& dest, 
			RawImage<T>& maskDest,
			gravis::d3Vector origin, 
			gravis::d3Vector spacing = 1.0, 
			int num_threads = 1,
			InterpolationType interpolation = Linear,
			double taperX = 20, 
			double taperY = 20, 
			double wMin = 3.0, 
			int frame0 = 0, 
			int frames = -1);
		
		template <typename T>
		static void backprojectRaw(
			const TomoStack<T>& stack,
			const std::vector<BufferedImage<T>>& images,
			RawImage<T>& dest, 
			RawImage<T>& maskDest,
			gravis::d3Vector origin, 
			gravis::d3Vector spacing = 1.0, 
			int num_threads = 1,
			InterpolationType interpolation = Linear,
			double taperX = 20, 
			double taperY = 20, 
			double wMin = 3.0, 
			int frame0 = 0, 
			int frames = -1);
};
						
						
template <typename SrcType, typename DestType>
void RealSpaceBackprojection::backproject(
				const RawImage<SrcType>& stack,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<DestType>& dest,
				int num_threads,
				gravis::d3Vector origin,
				double spacing,
				InterpolationType interpolation,
				double taper)
{
	const int fc = stack.zdim;
	
	const bool doTaper = taper != 0.0;
	
	#pragma omp parallel for num_threads(num_threads)	
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
	{
		double sum = 0.0;
		double wgh = 0.0;
		double taperMax = 0.0;
	
		gravis::d4Vector pw(
			origin.x + x * spacing, 
			origin.y + y * spacing, 
			origin.z + z * spacing, 
			1.0);
	
		for (int f = 0; f < fc; f++)
		{
			gravis::d4Vector pi = proj[f] * pw;
	
			if (pi.x >= 0.0 && pi.x < stack.xdim && pi.y >= 0.0 && pi.y < stack.ydim)
			{
				if (doTaper)
				{
					const double t = Tapering::getTaperWeight2D(
								pi.x, pi.y, stack.xdim, stack.ydim, taper);
					
					if (t > taperMax) taperMax = t;
				}
						
				if (interpolation == Linear)
				{
					sum += Interpolation::linearXY_clip(stack, pi.x, pi.y, f);
				}
				else
				{
					sum += Interpolation::cubicXY_clip(stack, pi.x, pi.y, f);
				}
				
				wgh += 1.0;
			}
		}
		
		if (doTaper) sum *= taperMax;
	
		if (wgh > 0.0) dest(x,y,z) += sum / wgh;
	}
}

template <typename SrcType, typename DestType>
void RealSpaceBackprojection::backprojectCoverage(
				const RawImage<SrcType>& stack,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<DestType>& dest,
				int num_threads,
				gravis::d3Vector origin,
				double spacing,
				InterpolationType interpolation,
				double taper)
{
	const int fc = stack.zdim;
	
	const bool doTaper = taper != 0.0;
	
	#pragma omp parallel for num_threads(num_threads)	
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
	{
		double wgh = 0.0;
		double taperMax = 0.0;
	
		gravis::d4Vector pw(
			origin.x + x * spacing, 
			origin.y + y * spacing, 
			origin.z + z * spacing, 
			1.0);
	
		for (int f = 0; f < fc; f++)
		{
			gravis::d4Vector pi = proj[f] * pw;
	
			if (pi.x >= 0.0 && pi.x < stack.xdim && pi.y >= 0.0 && pi.y < stack.ydim)
			{
				if (doTaper)
				{
					const double t = Tapering::getTaperWeight2D(
								pi.x, pi.y, stack.xdim, stack.ydim, taper);
					
					if (t > taperMax) taperMax = t;
				}
				
				wgh += 1.0;
			}
		}
		
		if (doTaper) 
		{
			dest(x,y,z) += wgh;
		}
		else 
		{
			dest(x,y,z) += taperMax * wgh;
		}
	}
}

template <typename SrcType, typename DestType>
void RealSpaceBackprojection::backprojectPsf(		
				const RawImage<SrcType>& stack,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<DestType>& dest,
				int num_threads,
				gravis::d3Vector origin,
				double spacing,
				InterpolationType interpolation)
{
	const int fc = proj.size();
	std::vector<gravis::d2Vector> originProj(fc);

    for (int f = 0; f < fc; f++)
    {
        gravis::d4Vector orig4 = proj[f] * gravis::d4Vector(dest.xdim/2, dest.ydim/2, dest.zdim/2, 1.0);
		originProj[f] = gravis::d2Vector(orig4.x, orig4.y);
    }
	
	#pragma omp parallel for num_threads(num_threads)	
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
	{
		double sum = 0.0;
	
		gravis::d4Vector pw(x, y, z, 1.0);
	
		for (int f = 0; f < fc; f++)
		{
			gravis::d4Vector pi = proj[f] * pw;
	
			double dx = (pi.x - originProj[f].x) * spacing;
			double dy = (pi.y - originProj[f].y) * spacing;
	
			if (interpolation == Linear)
			{
				sum += Interpolation::linearXYkernel(dx, dy);
			}
			else
			{
				sum += Interpolation::cubicXYkernel(dx, dy);
			}
		}
	
		dest(x,y,z) += sum;
	}
}

template <typename T>
void RealSpaceBackprojection::backprojectRaw(
				const TomoStack<T>& stack,
				RawImage<T>& dest, 
				RawImage<T>& maskDest,
				gravis::d3Vector origin, 
				gravis::d3Vector spacing, 
				int num_threads,
				InterpolationType interpolation,
				double taperX, 
				double taperY, 
				double wMin, 
				int frame0, 
				int frames)
{
	backprojectRaw(stack, stack.images, dest, maskDest, origin, spacing, num_threads,
				   interpolation, taperX, taperY, wMin, frame0, frames);
}

template <typename T>
void RealSpaceBackprojection::backprojectRaw(
				const TomoStack<T>& stack,
				const std::vector<BufferedImage<T>>& images,
				RawImage<T>& dest, 
				RawImage<T>& maskDest,
				gravis::d3Vector origin, 
				gravis::d3Vector spacing,
				int num_threads, 
				InterpolationType interpolation,
				double taperX, 
				double taperY, 
				double wMin, 
				int frame0, 
				int frames)
{
	gravis::d4Matrix vol2world;

    vol2world(0,0) = spacing.x;
    vol2world(1,1) = spacing.y;
    vol2world(2,2) = spacing.z;
    vol2world(0,3) = origin.x;
    vol2world(1,3) = origin.y;
    vol2world(2,3) = origin.z;
	
    const int ic = frames > 0? frames + frame0 : images.size();

    std::vector<gravis::d4Matrix> vol2img(ic);

    for (int im = 0; im < ic; im++)
    {
        vol2img[im] = stack.worldToImage[im] * vol2world;
    }
	
	#pragma omp parallel for num_threads(num_threads)	
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
    {
        double sum = 0.0;
        double wgh = 0.0;

        gravis::d4Vector pw(x,y,z,1.0);

        for (int im = frame0; im < ic; im++)
        {
            gravis::d4Vector pi = vol2img[im] * pw;

			if (pi.x >= 0.0 && pi.x < images[im].xdim-1 
				&& pi.y >= 0.0 && pi.y < images[im].ydim-1)
            {
                double wghi = Tapering::getTaperWeight2D(
							pi.x, pi.y, images[im].xdim, images[im].ydim, taperX);

                if (interpolation == Linear)
                {
                    sum += wghi * Interpolation::linearXY_clip(images[im], pi.x, pi.y, 0);
                }
                else
                {
                    sum += wghi * Interpolation::cubicXY_clip(images[im], pi.x, pi.y, 0);
                }

                wgh += wghi;
            }
        }

        if (wgh > 0.0)
        {
            sum /= wgh;
        }

        dest(x,y,z) = sum;
        maskDest(x,y,z) = wgh;
    }

    double mean = 0.0, sum = 0.0;
	
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
    {
        mean += maskDest(x,y,z) * dest(x,y,z);
        sum += maskDest(x,y,z);
    }

	if (sum > 0.0)
	{
		mean /= sum;
	}
	
	#pragma omp parallel for num_threads(num_threads)	
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
    {
        double t = maskDest(x,y,z) / wMin;

        if (t < 1.0)
        {
            dest(x,y,z) = t * dest(x,y,z) + (1.0 - t) * mean;
        }
    }
}

#endif

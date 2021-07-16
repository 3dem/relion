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
#include <src/jaz/tomography/tomogram.h>
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
			double taperFalloff = 20,
			double taperDist = 0);

		template <typename SrcType, typename DestType>
		static void backprojectSmooth(
			const RawImage<SrcType>& stack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<DestType>& dest,
			int num_threads = 1,
			gravis::d3Vector origin = gravis::d3Vector(0.0, 0.0, 0.0),
			double spacing = 1.0,
			InterpolationType interpolation = Linear,
			double taperFalloff = 20,
			double taperDist = 0);

		template <typename SrcType, typename DestType>
		static void backprojectCoverage(
			const RawImage<SrcType>& stack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<DestType>& dest,
			int num_threads = 1,
			gravis::d3Vector origin = gravis::d3Vector(0.0, 0.0, 0.0),
			double spacing = 1.0,
			InterpolationType interpolation = Linear,
			double taperFalloff = 20,
			double taperDist = 0);

		template <typename SrcType, typename DestType>
		static void backprojectSmoothCoverage(
			const RawImage<SrcType>& stack,
			const std::vector<gravis::d4Matrix>& proj,
			RawImage<DestType>& dest,
			int num_threads = 1,
			gravis::d3Vector origin = gravis::d3Vector(0.0, 0.0, 0.0),
			double spacing = 1.0,
			InterpolationType interpolation = Linear,
			double taperFalloff = 20,
			double taperDist = 0);
		
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
			const Tomogram& tomogram,
			RawImage<T>& dest,
			RawImage<T>& maskDest,
			gravis::d3Vector origin,
			gravis::d3Vector spacing = gravis::d3Vector(1.0,1.0,1.0),
			int num_threads = 1,
			InterpolationType interpolation = Linear,
			double taperFalloff = 20,
			double taperDist = 0,
			double wMin = 3.0);
		
		template <typename T>
		static void backprojectRaw(
			const std::vector<gravis::d4Matrix>& projections,
			const RawImage<T>& imageStack,
			RawImage<T>& dest, 
			RawImage<T>& maskDest,
			gravis::d3Vector origin, 
			gravis::d3Vector spacing = gravis::d3Vector(1.0,1.0,1.0),
			int num_threads = 1,
			InterpolationType interpolation = Linear,
			double taperFalloff = 20,
			double taperDist = 0,
			double wMin = 3.0);
		
		template <typename SrcType>
		static BufferedImage<SrcType> preWeight(
			const RawImage<SrcType>& stack,
			const std::vector<gravis::d4Matrix>& proj, 
			int num_threads = 1);
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
				double taperFalloff,
				double taperDist)
{
	const int fc = stack.zdim;

	const bool doTaper = taperFalloff != 0.0 || taperDist != 0.0;

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
								pi.x, pi.y, stack.xdim, stack.ydim, taperFalloff, taperDist);

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
void RealSpaceBackprojection::backprojectSmooth(
				const RawImage<SrcType>& stack,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<DestType>& dest,
				int num_threads,
				gravis::d3Vector origin,
				double spacing,
				InterpolationType interpolation,
				double taperFalloff,
				double taperDist)
{
	const int fc = stack.zdim;

	#pragma omp parallel for num_threads(num_threads)
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
	{
		double sum = 0.0;
		double wgh = 0.0;

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
				const double t = Tapering::getTaperWeight2D(
						pi.x, pi.y, stack.xdim, stack.ydim, taperFalloff, taperDist);

				if (interpolation == Linear)
				{
					sum += t * Interpolation::linearXY_clip(stack, pi.x, pi.y, f);
				}
				else
				{
					sum += t * Interpolation::cubicXY_clip(stack, pi.x, pi.y, f);
				}

				wgh += t;
			}
		}

		if (wgh > 1e-6)
		{
			dest(x,y,z) += sum / wgh;
		}
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
				double taperFalloff,
				double taperDist)
{
	const int fc = stack.zdim;

	const bool doTaper = taperFalloff != 0.0 || taperDist != 0.0;
	
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
								pi.x, pi.y, stack.xdim, stack.ydim, taperFalloff, taperDist);
					
					if (t > taperMax) taperMax = t;
				}
				
				wgh += 1.0;
			}
		}
		
		if (doTaper) 
		{
			dest(x,y,z) += taperMax * wgh;
		}
		else 
		{
			dest(x,y,z) += wgh;
		}
	}
}

template <typename SrcType, typename DestType>
void RealSpaceBackprojection::backprojectSmoothCoverage(
				const RawImage<SrcType>& stack,
				const std::vector<gravis::d4Matrix>& proj,
				RawImage<DestType>& dest,
				int num_threads,
				gravis::d3Vector origin,
				double spacing,
				InterpolationType interpolation,
				double taperFalloff,
				double taperDist)
{
	const int fc = stack.zdim;

	#pragma omp parallel for num_threads(num_threads)
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
	{
		double sum = 0.0;
		double wgh = 0.0;

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
				const double t = Tapering::getTaperWeight2D(
						pi.x, pi.y, stack.xdim, stack.ydim, taperFalloff, taperDist);

				if (interpolation == Linear)
				{
					sum += t * Interpolation::linearXY_clip(stack, pi.x, pi.y, f);
				}
				else
				{
					sum += t * Interpolation::cubicXY_clip(stack, pi.x, pi.y, f);
				}

				wgh += t;
			}
		}

		if (wgh > 1e-6)
		{
			dest(x,y,z) += wgh;
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
				const Tomogram& tomogram,
				RawImage<T>& dest, 
				RawImage<T>& maskDest,
				gravis::d3Vector origin, 
				gravis::d3Vector spacing, 
				int num_threads,
				InterpolationType interpolation,
				double taperFalloff,
				double taperDist,
				double wMin)
{
	backprojectRaw(
		tomogram.projectionMatrices, tomogram.stack, dest, maskDest,
		origin, spacing, num_threads,
		interpolation, taperFalloff, taperDist, wMin);
}

template <typename T>
void RealSpaceBackprojection::backprojectRaw(
				const std::vector<gravis::d4Matrix>& projections,
				const RawImage<T>& imageStack,
				RawImage<T>& dest, 
				RawImage<T>& maskDest,
				gravis::d3Vector origin, 
				gravis::d3Vector spacing,
				int num_threads, 
				InterpolationType interpolation,
				double taperFalloff,
				double taperDist,
				double wMin)
{
	gravis::d4Matrix vol2world;

	vol2world(0,0) = spacing.x;
	vol2world(1,1) = spacing.y;
	vol2world(2,2) = spacing.z;
	vol2world(0,3) = origin.x;
	vol2world(1,3) = origin.y;
	vol2world(2,3) = origin.z;

	const int w = imageStack.xdim;
	const int h = imageStack.ydim;
	const int ic = imageStack.zdim;

	std::vector<gravis::d4Matrix> vol2img(ic);

	for (int im = 0; im < ic; im++)
	{
		vol2img[im] = projections[im] * vol2world;
	}
	
	#pragma omp parallel for num_threads(num_threads)	
	for (size_t z = 0; z < dest.zdim; z++)
	for (size_t y = 0; y < dest.ydim; y++)
	for (size_t x = 0; x < dest.xdim; x++)
	{
		double sum = 0.0;
		double wgh = 0.0;

		gravis::d4Vector pw(x,y,z,1.0);

		for (int im = 0; im < ic; im++)
		{
			gravis::d4Vector pi = vol2img[im] * pw;

			if (pi.x >= 0.0 && pi.x < w-1
					&& pi.y >= 0.0 && pi.y < h-1)
			{
				double wghi = Tapering::getTaperWeight2D(
							pi.x, pi.y, w, h, taperFalloff, taperDist);

				if (interpolation == Linear)
				{
					sum += wghi * Interpolation::linearXY_clip(imageStack, pi.x, pi.y, im);
				}
				else
				{
					sum += wghi * Interpolation::cubicXY_clip(imageStack, pi.x, pi.y, im);
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

template<typename SrcType>
BufferedImage<SrcType> RealSpaceBackprojection::preWeight(
        const RawImage<SrcType>& stack, 
        const std::vector<gravis::d4Matrix>& proj, 
        int num_threads)
{
	const int w  = stack.xdim;
	const int h  = stack.ydim;
	const int fc = stack.zdim;
	const int wh = w/2 + 1;
	
	BufferedImage<SrcType> out(w,h,fc);
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		std::vector<gravis::d3Vector> f_to_others_z(fc);
		
		gravis::d3Matrix proj_f = gravis::d3Matrix::extract(proj[f]);
		proj_f.transpose();
		        
		for (int ff = 0; ff < fc; ff++)	
		{
			gravis::d3Matrix proj_ff = gravis::d3Matrix::extract(proj[ff]);
			proj_ff.transpose();
			proj_ff.invert();
			
			gravis::d3Matrix A = proj_ff * proj_f;
			
			f_to_others_z[ff][0] = A(2,0);
			f_to_others_z[ff][1] = A(2,1);
			f_to_others_z[ff][2] = A(2,2);
		}
				        
		BufferedImage<SrcType> frameRS = stack.getConstSliceRef(f);
		BufferedImage<tComplex<SrcType>> frameFS;
		
		FFT::FourierTransform(frameRS, frameFS);
		
		for (int xi = 0; xi < wh; xi++)
		for (int yi = 0; yi < h; yi++)
		{
			const double xx = xi;
			const double yy = yi < h/2? yi : yi - h;
			
			const gravis::d3Vector r(xx,yy,0.0);
			
			double sum = 0.0;
			
			for (int ff = 0; ff < fc; ff++)	
			{
				const double z_f = r.dot(f_to_others_z[ff]);
				const double weight = 1.0 - std::abs(z_f);
				
				if (weight > 0.0)
				{
					sum += weight;
				}
			}
			
			frameFS(xi,yi) *= (SrcType)(1.0 / (sum));
		}
		
		FFT::inverseFourierTransform(frameFS, frameRS);
		
		out.copySliceFrom(f, frameRS);
	}
	
	return out;
}

#endif

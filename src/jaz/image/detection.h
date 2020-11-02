#ifndef DETECTION_H
#define DETECTION_H

#include "raw_image.h"
#include "normalization.h"
#include "local_extrema.h"
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/projection/projection.h>

class Detection
{
	public:
		
		template <typename T>
		static BufferedImage<T> smallCircleKernel(double avgRad, int w, int h);
		
		template <typename T>
		static BufferedImage<T> circleKernel(double minRad, double maxRad, int w, int h);
		
		template <typename T>
		static std::vector<gravis::d3Vector> findLocalMaxima(
				const Tomogram& tomogram,
				const RawImage<T>& similarities2D,
				gravis::d3Vector origin, 
				gravis::d3Vector spacing,
				gravis::d3Vector diagonal,
				T minValue, 
				int maxPeaks,
				double minDist,
				int num_threads,
				int binning,
				std::string diagFn);
};

template <typename T>
BufferedImage<T> Detection::smallCircleKernel(double avgRad, int w, int h)
{
	BufferedImage<T> out(w,h);
	
	double sum = 0.0;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		double xx = x < w/2? x : w - x;
		double yy = y < h/2? y : h - y;
		
		const double r = sqrt(xx*xx + yy*yy);
				
		if (r < avgRad)
		{
			out(x,y) = -1.0;
		}
		else
		{
			out(x,y) = exp(-r/avgRad);
		}
		
		sum += out(x,y);
	}
	
	const double mean = sum / (w*h);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y) -= mean;
	}
	
	double sum2 = 0.0;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		sum2 += out(x,y) * out(x,y);
	}
	
	const double nrm = sqrt(sum2/(w*h));
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y) /= nrm;
	}
	
	return out;
}

template <typename T>
BufferedImage<T> Detection::circleKernel(double minRad, double maxRad, int w, int h)
{
	BufferedImage<T> out(w,h);
	
	double sum = 0.0;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		double xx = x < w/2? x : w - x;
		double yy = y < h/2? y : h - y;
		
		const double r = sqrt(xx*xx + yy*yy);
				
		if (r < minRad)
		{
			out(x,y) = 0.0;
		}
		else if (r < maxRad)
		{
			const double t = (r - minRad) / (maxRad - minRad);
			out(x,y) = -sin(2.0 * PI * t) / r;
		}
		else
		{			
			out(x,y) = 0.0;
		}
		
		sum += out(x,y);
	}
	
	const double mean = sum / (w*h);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y) -= mean;
	}
	
	double sum2 = 0.0;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		sum2 += out(x,y) * out(x,y);
	}
	
	const double nrm = sqrt(sum2/(w*h));
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y) /= nrm;
	}
	
	return out;
}

template <typename T>
std::vector<gravis::d3Vector> Detection::findLocalMaxima(
	const Tomogram& tomogram,
	const RawImage<T>& similarities2D,
	gravis::d3Vector origin, 
	gravis::d3Vector spacing,
	gravis::d3Vector diagonal,
	T minValue, 
	int maxPeaks, 
	double minDist,
	int num_threads,
	int binning,
	std::string diagFn)
{
	const int w2D = similarities2D.xdim;
	const int h2D = similarities2D.ydim;
	
	const int wt2D = tomogram.stack.xdim;
	const int ht2D = tomogram.stack.ydim;
	const int fc = tomogram.stack.zdim;


	
	if (w2D != wt2D || h2D != ht2D)
	{
		REPORT_ERROR_STR("Detection::findLocalMaxima: images are of incorrect size: "
						 << w2D << "x" << h2D << " vs. " << wt2D << "x" << ht2D);
	}
		
	const int w2Db = w2D / binning;
	const int h2Db = h2D / binning;
	
	BufferedImage<T> binnedSimilarity(w2Db, h2Db, fc);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<float> binned = Resampling::downsampleMax_2D_full(
					similarities2D.getConstSliceRef(f), w2Db, h2Db);

		binnedSimilarity.getSliceRef(f).copyFrom(binned);
	}

	Tomogram binnedTomogram = tomogram.FourierCrop(binning, num_threads, false);
	
	const int w3D = diagonal.x;
	const int h3D = diagonal.y;
	const int d3D = diagonal.z;
	
	const int w3Db = w3D / binning;
	const int h3Db = h3D / binning;
	const int d3Db = d3D / binning;


	BufferedImage<T> coarseVol(w3Db, h3Db, d3Db), coarseMask(w3Db, h3Db, d3Db);
			
	RealSpaceBackprojection::backprojectRaw(
		binnedTomogram.projectionMatrices, binnedSimilarity,
		coarseVol, coarseMask,
		origin, spacing * binning, num_threads, 
		RealSpaceBackprojection::Linear,
		20, 0, 10.0);
	
	coarseVol = Normalization::byNormalDist(coarseVol);
	
	if (diagFn != "")
	{
		coarseVol.write(
			diagFn + "detection_binned-CC.mrc",
			binnedTomogram.optics.pixelSize);
	}
	
	
	BufferedImage<T> boxMax = LocalExtrema::boxMaxima(coarseVol, minDist/binning);
			
	
	if (diagFn != "")
	{
		boxMax.write(
			diagFn + "detection_coarse-max.mrc",
			binnedTomogram.optics.pixelSize);
	}
	
	std::vector<gravis::d3Vector> coarseMaxima, fineMaxima;
	coarseMaxima.reserve(1024);
	fineMaxima.reserve(1024);

	const double localSearchRadius = 3.0 * binning;
	const int localSearchRegion = 2 * localSearchRadius;
	
	BufferedImage<T> peaks(w3Db,h3Db,d3Db);
	BufferedImage<T> subvolume(localSearchRegion, localSearchRegion, localSearchRegion);
	BufferedImage<T> subvolumeMask(localSearchRegion, localSearchRegion, localSearchRegion);
	
	
	for (int z = 0; z < d3Db; z++)
	for (int y = 0; y < h3Db; y++)
	for (int x = 0; x < w3Db; x++)
	{
		if (coarseVol(x,y,z) > minValue && coarseVol(x,y,z) == boxMax(x,y,z))
		{
			const gravis::d3Vector p = 
				origin + (double)binning * gravis::d3Vector(
						spacing.x * x,
						spacing.y * y,
						spacing.z * z);
			
			coarseMaxima.push_back(p);
			
			peaks(x,y,z) = coarseVol(x,y,z);
			
			if (coarseMaxima.size() > maxPeaks)
			{
				REPORT_ERROR_STR("Detection::findLocalMaxima: maximum number of local maxima exceeded ("
								 << maxPeaks << ") - please increase the cutoff threshold (minValue) "
								 << "or the maximum allowed number of peaks (maxPeaks)" );
			}
			
			RealSpaceBackprojection::backprojectRaw(
				tomogram.projectionMatrices,
				similarities2D,
				subvolume, subvolumeMask,
				p - gravis::d3Vector(localSearchRadius),
				gravis::d3Vector(1.0), 
				num_threads, 
				RealSpaceBackprojection::Linear,
				20, 0, 10.0);
			
			double maxVal = subvolume(0,0,0);
			gravis::d3Vector bestInd(0,0,0);
					
			for (int zz = 0; zz < localSearchRegion; zz++)
			for (int yy = 0; yy < localSearchRegion; yy++)
			for (int xx = 0; xx < localSearchRegion; xx++)
			{
				if (subvolume(xx,yy,zz) > maxVal)
				{
					maxVal = subvolume(xx,yy,zz);
					bestInd = gravis::d3Vector(xx,yy,zz);
				}
			}
			
			gravis::d3Vector p2 = p - gravis::d3Vector(localSearchRadius) + bestInd;

			fineMaxima.push_back(p2);
		}
		else
		{
			peaks(x,y,z) = 0;
		}		
	}
	
	return fineMaxima;
	
}

#endif

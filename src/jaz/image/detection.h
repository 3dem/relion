#ifndef DETECTION_H
#define DETECTION_H

#include "raw_image.h"
#include "normalization.h"
#include "local_extrema.h"
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomography/tomo_stack.h>
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
				const TomoStack<T>& tomogram, 
				const std::vector<BufferedImage<T>>& similarities2D,
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
	const TomoStack<T>& tomogram, 
	const std::vector<BufferedImage<T>>& similarities2D,
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
	const int w2D = similarities2D[0].xdim;
	const int h2D = similarities2D[0].ydim;
	
	const int wt2D = tomogram.images[0].xdim;
	const int ht2D = tomogram.images[0].ydim;
	
	const int fc = tomogram.images.size();
	
	
	if (w2D != wt2D || h2D != ht2D)
	{
		REPORT_ERROR_STR("Detection::findLocalMaxima: images are of incorrect size: "
						 << w2D << "x" << h2D << " vs. " << wt2D << "x" << ht2D);
	}
		
	const int w2Db = w2D / binning;
	const int h2Db = h2D / binning;
	
	std::vector<BufferedImage<T>> simBin(fc);
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		simBin[f] = Resampling::downsampleMax_2D_full(similarities2D[f], w2Db, h2Db);
	}
					
	TomoStack<float> tsb = tomogram.downsample(binning, num_threads, false);
	
	const int w3D = diagonal.x;
	const int h3D = diagonal.y;
	const int d3D = diagonal.z;
	
	const int w3Db = w3D / binning;
	const int h3Db = h3D / binning;
	const int d3Db = d3D / binning;
	
	BufferedImage<T> coarseVol(w3Db, h3Db, d3Db), coarseMask(w3Db, h3Db, d3Db);
			
	RealSpaceBackprojection::backprojectRaw(
		tsb, simBin, coarseVol, coarseMask, 
		origin, spacing * binning, num_threads, 
		RealSpaceBackprojection::Linear,
		20, 20, 10.0);
	
	coarseVol = Normalization::byNormalDist(coarseVol);
	
	if (diagFn != "")
	{
		coarseVol.writeVtk(diagFn + "detection_binned-CC.vtk", origin, spacing * binning);
	}
	
	
	BufferedImage<T> boxMax = LocalExtrema::boxMaxima(coarseVol, minDist/binning);
			
	
	if (diagFn != "")
	{
		boxMax.writeVtk(diagFn + "detection_coarse-max.vtk", origin, spacing * binning);
	}
	
	std::vector<gravis::d3Vector> coarseMaxima, fineMaxima;
	coarseMaxima.reserve(1024);
	fineMaxima.reserve(1024);
	
	BufferedImage<T> peaks(w3Db,h3Db,d3Db);
	BufferedImage<T> subvolume(2*binning, 2*binning, 2*binning);	
	BufferedImage<T> subvolumeMask(2*binning, 2*binning, 2*binning);	
	
	
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
				tomogram, similarities2D, subvolume, subvolumeMask, 
				p - gravis::d3Vector(binning), 
				gravis::d3Vector(1.0), 
				num_threads, 
				RealSpaceBackprojection::Linear,
				20, 20, 10.0);
			
			double maxVal = subvolume(0,0,0);
			gravis::d3Vector bestInd(0,0,0);
					
			for (int zz = 0; zz < 2*binning; zz++)
			for (int yy = 0; yy < 2*binning; yy++)
			for (int xx = 0; xx < 2*binning; xx++)
			{
				if (subvolume(xx,yy,zz) > maxVal)
				{
					maxVal = subvolume(xx,yy,zz);
					bestInd = gravis::d3Vector(xx,yy,zz);
				}
			}
			
			gravis::d3Vector p2 = p - gravis::d3Vector(binning) + bestInd;
						
			fineMaxima.push_back(p2);
		}
		else
		{
			peaks(x,y,z) = 0;
		}		
	}
	
	/*{
		peaks.writeVtk("dev/peaks.vtk", origin, spacing * binning);
		
		const gravis::d3Vector cent = coarseMaxima[0];
		const int s = 300;
		const gravis::d3Vector halfDiag(s/2.0);
				
		TomoStack<float> tsVes = tomogram.extractSubStack(cent, 4*s/3, 4*s/3);
		
		tsVes.saveImages("dev/ves.vtk");
	
		Image<float> ves1(s,s,s), maskVes1(s,s,s);
		
		
		BackprojectionHelper::backprojectRaw(
			tsVes, ves1, maskVes1, 
			cent - halfDiag, spacing, num_threads, 
			BackprojectionHelper::Linear,
			20, 20, 10.0);
		
		ves1.writeVtk("dev/ves0.vtk", cent - halfDiag, spacing);
		
		
		
		BackprojectionHelper::backprojectRaw(
			tomogram, ves1, maskVes1, 
			cent - halfDiag, spacing, num_threads, 
			BackprojectionHelper::Linear,
			20, 20, 10.0);
		
		ves1.writeVtk("dev/ves0_0.vtk", cent - halfDiag, spacing);
	}*/
	
	return fineMaxima;
	
}

#endif

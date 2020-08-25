#include "cells.h"
#include <src/jaz/image/padding.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/math/fft.h>

using namespace gravis;

void Cells::findCenters(
	const std::vector<d3Vector>& surfacePoints, 
	int bins,
	float maxRadius,
	BufferedImage<float>& outMaxima, 
	BufferedImage<float>& outRadii,
	int num_threads)
{
	const int w = outMaxima.xdim;
	const int h = outMaxima.ydim;
	const int d = outMaxima.zdim;
	
	double dist2bin = (double) bins / (double) maxRadius;
	
	#pragma omp parallel for num_threads(num_threads)		
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		d3Vector p0(x,y,z);
		
		std::vector<float> hist(bins, 0.f);
		
		for (int i = 0; i < surfacePoints.size(); i++)
		{
			double d =  dist2bin * (surfacePoints[i] - p0).length() ;
			
			int d0 = (int) d;
			int d1 = d0 + 1;
			double df = d - d0;
			
			if (d0 < bins)
			{
				hist[d0] += 1.0 - df;
			}
			
			if (d1 < bins)
			{
				hist[d1] += df;
			}
		}
		
		float maxHist = 0, maxRad = 0;
		
		for (int i = 0; i < bins; i++)
		{
			if (hist[i] > maxHist)
			{
				maxHist = hist[i];
				maxRad = i;
			}
		}
		
		outMaxima(x,y,z) = maxHist;
		outRadii(x,y,z) = maxRad / dist2bin;
	}
}

void Cells::findCenters(
	const BufferedImage<float> &surface, 
	int bins, 
	float maxRadius, 
	float threshVal,
	BufferedImage<float> &outMaxima, 
	BufferedImage<float> &outRadii, 
	int num_threads)
{
	const int w = surface.xdim;
	const int h = surface.ydim;
	const int d = surface.zdim;
	const int wh = w/2 + 1;
	
	double dist2bin = (double) bins / (double) maxRadius;
	
	int padding = (int) (maxRadius + 0.5f);
	
	const int w2 = w + 2 * padding;
	const int h2 = h + 2 * padding;
	const int d2 = d + 2 * padding;
	const int w2h = w2/2 + 1;
			
	
	BufferedImage<float> surfPad = Padding::padCenter3D_full(surface, padding);
	
	for (int z = 0; z < d2; z++)
	for (int y = 0; y < h2; y++)
	for (int x = 0; x < w2; x++)
	{
		surfPad(x,y,z) = surfPad(x,y,z) > threshVal? 1 : 0;
	}
	
	
	BufferedImage<fComplex> surfPadFS;
	FFT::FourierTransform(surfPad, surfPadFS, FFT::Both);
	
	outMaxima.fill(0.f);
	outRadii.fill(0.f);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int b = 0; b < bins; b++)
	{
		const double r0 = (b + 0.5) / dist2bin;
		std::cout << b << " / " << bins << ": " << r0 << std::endl;
		
		BufferedImage<float> byDim[3];
		
		for (int dim = 0; dim < 3; dim++)
		{
			std::cout << "  " << dim << std::endl;
			
			byDim[dim] = BufferedImage<float>(w,h,d);
			
			BufferedImage<float> bySign[2];
			
			for (int sg = 0; sg < 2; sg++)
			{
				const double sign = 2*sg - 1;
				
				std::cout << "    " << sign << std::endl;
				
				d3Vector dir;
				dir[dim] = sign;
				dir[(dim+1)%3] = 0;
				dir[(dim+2)%3] = 0;		
				
				BufferedImage<float> kernelRS(w2,h2,d2);		
				
				for (int z = 0; z < d2; z++)
				for (int y = 0; y < h2; y++)
				for (int x = 0; x < w2; x++)
				{
					const double xx = x < w2/2? x : x - w2;
					const double yy = y < h2/2? y : y - h2;
					const double zz = z < d2/2? z : z - d2;
					
					d3Vector v(xx,yy,zz);
										
					const double r = sqrt(xx*xx + yy*yy + zz*zz);
					
					double cs = v.dot(dir) / r;
					
					const double db = dist2bin * std::abs(r - r0);
								
					kernelRS(x,y,z) = (db < 1.0 && cs > 0.0)? cs : 0.0;
				}
				
				BufferedImage<fComplex> kernelFS;
				FFT::FourierTransform(kernelRS, kernelFS, FFT::Both);
				
				for (int z = 0; z < d2; z++)
				for (int y = 0; y < h2; y++)
				for (int x = 0; x < w2h; x++)
				{
					kernelFS(x,y,z) = surfPadFS(x,y,z) * kernelFS(x,y,z).conj();
				}
				
				FFT::inverseFourierTransform(kernelFS, bySign[sg], FFT::Both);
				
				bySign[sg].write("debug/bySign_R"
								 +ZIO::itoa(b)
								 +"_"
								 +ZIO::itoa(sign*(dim+1))
								 +".mrc");
			}
			
			for (int z = 0; z < d; z++)
			for (int y = 0; y < h; y++)
			for (int x = 0; x < w; x++)
			{
				const double a = bySign[0](x+padding, y+padding, z+padding);
				const double b = bySign[1](x+padding, y+padding, z+padding);
				
				byDim[dim](x,y,z) = a > b? b : a;
			}
			
			byDim[dim].write("debug/byDim_R"
							 +ZIO::itoa(b)
							 +"_"
							 +ZIO::itoa((dim+1))
							 +".mrc");
		}
		
		BufferedImage<float> byRad(w,h,d);
		
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			const double a = byDim[0](x,y,z);
			const double b = byDim[1](x,y,z);
			const double c = byDim[2](x,y,z);
			
			byRad(x,y,z) = a > b? (b > c? c : b) : (a > c? c : a);
		}
		
		for (int z = 0; z < d; z++)
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			if (outMaxima(x,y,z) < byRad(x,y,z))
			{
				outMaxima(x,y,z) = byRad(x,y,z);
				outRadii(x,y,z) = r0;
			}
		}
	}	
}

BufferedImage<float> Cells::findPointSymmetries(
		const BufferedImage<float> &img, 
		double highPass,
		double lowPass)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	const int wh = w/2 + 1;
	
	BufferedImage<float> imgCp = img;
	BufferedImage<fComplex> imgFS;
	
	FFT::FourierTransform(imgCp, imgFS, FFT::Both);
			
	BufferedImage<fComplex> ac(wh,h,d);
	
	double power = 0.0;
			
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		const fComplex zz = imgFS(x,y,z);
		ac(x,y,z) = zz*zz;
		power += zz.norm();
	}
	
	const float scale = (float) (1.0 / power);
	const double hp2 = highPass * highPass;
	const double lp2 = lowPass * lowPass;
	
	if (hp2 > 0.0)
	{
		for (long int z = 0; z < d; z++)
		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < wh; x++)
		{
			double xx = x;
			double yy = y < h/2? y : y - h;
			double zz = z < d/2? z : z - d;
			
			double r2 = xx*xx + yy*yy + zz*zz;
			
			ac(x,y,z) *= (1 - exp(-0.5 * r2 / hp2)) * exp(-0.5 * r2 / lp2);
		}
	}
	
	BufferedImage<fComplex> ac_half(wh,h,d);
	
	for (long int z = 0; z < d; z++)
	for (long int y = 0; y < h; y++)
	for (long int x = 0; x < wh; x++)
	{
		double center = (1 - 2*(x%2)) * (1 - 2*(y%2)) * (1 - 2*(z%2));
		
		double xx = x;
		double yy = y < h/2? y : y - h;
		double zz = z < d/2? z : z - d;
		
		fComplex half_val = Interpolation::linearXYZ_FftwHalf_complex(
					ac, xx/2.0, yy/2.0, zz/2.0);
		
		ac_half(x,y,z) = center * half_val;
	}
	
	BufferedImage<float> acRS;
	FFT::inverseFourierTransform(ac_half, acRS, FFT::Both);
	
	return acRS;
}

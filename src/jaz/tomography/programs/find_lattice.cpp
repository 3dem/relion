#include "find_lattice.h"
#include <src/jaz/image/filter.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/structure_tensor.h>
#include <src/jaz/image/local_extrema.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/single_particle/vtk_helper.h>

#include <omp.h>

using namespace gravis;


void FindLatticeProgram::run()
{
	BufferedImage<float> tomo;
	tomo.read(tomoFn);
	
	if (noise)
	{
		#pragma omp parallel for num_threads(n_threads)	
		for (size_t z = 0; z < tomo.zdim; z++)
		for (size_t y = 0; y < tomo.ydim; y++)
		for (size_t x = 0; x < tomo.xdim; x++)
		{
			tomo(x,y,z) = (float)(2.0 * rand() / (double)RAND_MAX - 1.0);
		}
	}
	
	if (taper > 0.0)
	{
		Tapering::taper(tomo, taper, n_threads, false);
	}
	
	const double spacing_px = spacing_ang / angpix;
		
	BufferedImage<float> tomoBP = ImageFilter::bandpass(tomo, spacing_px, filter_width, overtones);
	
	tomoBP.write(outFn+"_bandpass.mrc");
	
	BufferedImage<float> latticeDensity(tomoBP.xdim, tomoBP.ydim, tomoBP.zdim);
	
	{
		BufferedImage<float> tomoBP_pos = ImageFilter::bandpass(
					tomo, spacing_px + 2, filter_width, overtones);
		
		BufferedImage<float> tomoBP_neg = ImageFilter::bandpass(
					tomo, spacing_px - 2, filter_width, overtones);
		
		BufferedImage<float> tomoBP_env = tomoBP_pos + tomoBP_neg;
		
		tomoBP_env = ImageFilter::Gauss3D(tomoBP_env * tomoBP_env, spacing_px/2.0);
		BufferedImage<float> tomoBP_gauss = ImageFilter::Gauss3D(tomoBP * tomoBP, spacing_px/2.0);
	
				
		#pragma omp parallel for num_threads(n_threads)	
		for (size_t z = 0; z < tomoBP.zdim; z++)
		for (size_t y = 0; y < tomoBP.ydim; y++)
		for (size_t x = 0; x < tomoBP.xdim; x++)
		{
			const double v = tomoBP_gauss(x,y,z);
			const double ve = tomoBP_env(x,y,z);
			
			latticeDensity(x,y,z) = v - ve;
		}
		
		latticeDensity.write(outFn+"_density.mrc");
	}
	
	
	
	BufferedImage<float> boxMax = LocalExtrema::boxMaxima(tomoBP, (int)(spacing_px/sqrt(3.0)));
	BufferedImage<float> boxMin = LocalExtrema::boxMinima(tomoBP, (int)(spacing_px/sqrt(3.0)));
	
	std::vector<d3Vector> maxima(0), safeMaxima(0), minima(0), safeMinima(0);
	maxima.reserve(1024);
	safeMaxima.reserve(1024);
	minima.reserve(1024);
	safeMinima.reserve(1024);
	
	for (size_t z = 0; z < tomoBP.zdim; z++)
	for (size_t y = 0; y < tomoBP.ydim; y++)
	for (size_t x = 0; x < tomoBP.xdim; x++)
	{
		if (tomoBP(x,y,z) > minValue && tomoBP(x,y,z) == boxMax(x,y,z))
		{
			const d3Vector p = Interpolation::localQuadraticMaxXYZ(tomoBP, x, y, z);
			
			maxima.push_back(p);
					
			if (latticeDensity(x,y,z) > minDensity)
			{
				safeMaxima.push_back(p);
			}
		}
		
		if (tomoBP(x,y,z) < -minValue && tomoBP(x,y,z) == boxMin(x,y,z))
		{
			const d3Vector p = Interpolation::localQuadraticMinXYZ(tomoBP, x, y, z);
			
			minima.push_back(p);
					
			if (latticeDensity(x,y,z) > minDensity)
			{
				safeMinima.push_back(p);
			}
		}
	}
	
	const int xc = maxima.size();
	const int sxc = safeMaxima.size();
	
	const int nc = minima.size();
	const int snc = safeMinima.size();
	
	std::cout << xc << " maxima found, " << sxc << " safe(r).\n";
	std::cout << nc << " minima found, " << snc << " safe(r).\n";
	
	{
		std::ofstream ptstream(outFn+"_maxima.csv");
		
		for (int i = 0; i < xc; i++)
		{
			const d3Vector p = maxima[i];
			ptstream << p.x << "," << p.y << "," << p.z << "\n";
		}
	}
	
	{
		std::ofstream sptstream(outFn+"_safe-maxima.csv");
					
		for (int i = 0; i < sxc; i++)
		{
			const d3Vector p = safeMaxima[i];
			sptstream << p.x << "," << p.y << "," << p.z << "\n";
		}
	}
	
	{
		std::ofstream ptstream(outFn+"_minima.csv");
		
		for (int i = 0; i < nc; i++)
		{
			const d3Vector p = minima[i];
			ptstream << p.x << "," << p.y << "," << p.z << "\n";
		}
	}
	
	{
		std::ofstream sptstream(outFn+"_safe-minima.csv");
					
		for (int i = 0; i < snc; i++)
		{
			const d3Vector p = safeMinima[i];
			sptstream << p.x << "," << p.y << "," << p.z << "\n";
		}
	}
}

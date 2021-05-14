#include "dark_erase.h"
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/conversion.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/volume_converter.h>
#include <src/jaz/util/image_file_helper.h>

#include <omp.h>

void DarkEraseProgram::run()
{
	std::cout << "loading images..." << std::endl;
	
	BufferedImage<float> stack0;
	stack0.read(stackFn);
	
	const int w0  = stack0.xdim;
	const int h0  = stack0.ydim;
	const int fc = stack0.zdim;
	
	
	const int bin = 8;
	const double minFreq0 = 500.0;
	const double maxFreq0 = 64.0;
	
	const double minFreq = minFreq0 / bin;
	const double maxFreq = maxFreq0 / bin;
	
	std::cout << "resampling images..." << std::endl;
		
	BufferedImage<float> stackBin0 = Resampling::downsampleFiltStack_2D_full(stack0, bin, num_threads);
	
	const int wb = stackBin0.xdim;
	const int hb = stackBin0.ydim;
	
	
	BufferedImage<float> stackBin = stackBin0;
	
	std::cout << "highpass filtering images at " << minFreq << " px^-1..." << std::endl;
	
	stackBin = ImageFilter::highpassStack(stackBin, minFreq, 10.0, true);
	stackBin = ImageFilter::lowpassStack(stackBin, maxFreq, 10.0, true);
	
	std::cout << "normalizing images..." << std::endl;
	
	BufferedImage<float> nrm(wb,hb,fc);
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		double avg = 0.0;
				
		for (int y = 0; y < hb; y++)
		for (int x = 0; x < wb; x++)
		{
			avg += stackBin(x,y,f);
		}
			
		avg /= wb * hb;
		
		double var = 0.0;
		
		for (int y = 0; y < hb; y++)
		for (int x = 0; x < wb; x++)
		{
			const double d = stackBin(x,y,f) - avg;
			var += d * d;
		}
			
		var /= wb * hb - 1;
		
		const double sd = sqrt(var);
		
		for (int y = 0; y < hb; y++)
		for (int x = 0; x < wb; x++)
		{
			nrm(x,y,f) = (stackBin(x,y,f) - avg) / sd;
		}
	}
	
	if (writeNormalized)
	{
		nrm.write(outFn);		
		return;
	}
	
	std::cout << "constructing masks..." << std::endl;
	
	BufferedImage<float> mask(wb,hb,fc);
		
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < hb; y++)
		for (int x = 0; x < wb; x++)
		{
			if (nrm(x,y,f) > thresh)
			{
				mask(x,y,f) = 1.0;
			}
			else
			{
				mask(x,y,f) = 0.0;
			}
		}
	}
	
	if (diag)
	{
		mask.write("mask0.mrc");
	}
	
	mask = ImageFilter::GaussStack(mask, rad / bin, true);
	
	
	if (diag)
	{
		mask.write("maskLP0.mrc");
	}
		
	const double t90 = 0.9;
			
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < hb; y++)
		for (int x = 0; x < wb; x++)
		{
			if (mask(x,y,f) < t90) mask(x,y,f) = 0.0;
			else mask(x,y,f) = (mask(x,y,f) - t90) / (1.0 - t90);
		}
	}
	
	if (diag)
	{
		mask.write("maskLP1.mrc");
	}
	
	mask = ImageFilter::GaussStack(mask, rad / bin, true);
	
	if (diag)
	{
		mask.write("maskLP2.mrc");
	}
	
	BufferedImage<float> fillMask(wb,hb,fc);
	
	const double t0 = 0.3;
	const double t1 = 0.95;
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < hb; y++)
		for (int x = 0; x < wb; x++)
		{
			const double m = mask(x,y,f);
			
			if (m > t1)
			{
				fillMask(x,y,f) = 1.0;
				mask(x,y,f) = 1.0;
			}
			else if (m > t0)
			{
				fillMask(x,y,f) = 0.0;
				mask(x,y,f) = (m - t0) / (t1 - t0);
			}
			else
			{
				fillMask(x,y,f) = 0.0;
				mask(x,y,f) = 0.0;
			}
		}
	}
	
	if (diag)
	{
		mask.write("mask_fin.mrc");
		fillMask.write("fillMask.mrc");
	}
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < hb; y++)
		for (int x = 0; x < wb; x++)
		{
			const double m = mask(x,y,f);
			
			if (m > t1)
			{
				fillMask(x,y,f) = 1.0;
				mask(x,y,f) = 1.0;
			}
			else if (m > t0)
			{
				fillMask(x,y,f) = 0.0;
				mask(x,y,f) = (m - t0) / (t1 - t0);
			}
			else
			{
				fillMask(x,y,f) = 0.0;
				mask(x,y,f) = 0.0;
			}
		}
	}
	
	std::cout << "filling gaps..." << std::endl;
	
	BufferedImage<float> maskedStack = fillMask * stackBin0;
	
	if (diag)
	{
		maskedStack.write("maskedStack.mrc");
	}
	
	BufferedImage<double> extendedMaskedStack = ImageFilter::GaussStack(
			Conversion::toDouble(maskedStack), 2.0 * rad / bin, true);
	
	BufferedImage<double> extendedMask = ImageFilter::GaussStack(
			Conversion::toDouble(fillMask), 2.0 * rad / bin, true);
	
	if (diag)
	{
		extendedMaskedStack.write("extendedMaskedStack.mrc");
		extendedMask.write("extendedMask.mrc");
	}
	
	const double eps = 1e-16;
	
	BufferedImage<float> filledStack(wb,hb,fc);
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < hb; y++)
		for (int x = 0; x < wb; x++)
		{
			const double n = extendedMaskedStack(x,y,f);
			const double d = extendedMask(x,y,f);
			
			if (d > eps) filledStack(x,y,f) = (float)(n / d);
			else filledStack(x,y,f) = (float)(n / eps);
		}
	}
	
	filledStack = ImageFilter::GaussStack(filledStack, 2.0 * rad / bin, true);
	
	if (diag)
	{
		filledStack.write("filledStack.mrc");
	}
	
	std::cout << "applying masks..." << std::endl;
		
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		BufferedImage<float> mask_f = Resampling::upsampleCubic_2D_full(mask, f, bin, w0, h0);	
		BufferedImage<float> fill_f = Resampling::upsampleCubic_2D_full(filledStack, f, bin, w0, h0);
		
		for (int y = 0; y < h0; y++)
		for (int x = 0; x < w0; x++)
		{
			const float a = mask_f(x,y);			
			stack0(x,y,f) = a * stack0(x,y,f) + (1.0 - a) * fill_f(x,y);
		}
	}
	
	stack0.write(outFn);
}

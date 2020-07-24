#include "tomo_backproject.h"
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/args.h>
#include <src/jaz/gravis/t4Matrix.h>

#include <omp.h>

using namespace gravis;


void TomoBackprojectProgram::run()
{
	TomogramSet tomogramSet(tomoSetFn);
	Tomogram tomogram = tomogramSet.loadTomogram(tomoIndex, true);
	
	const int w0 = tomogram.w0;
	const int h0 = tomogram.h0;
	
	std::cout << w0 << " x " << h0 << std::endl;
	
	if (thickness < 0)
	{
		thickness = tomogram.d0;
		std::cout << "Using thickness from '" << tomoSetFn << "': " << thickness << std::endl;
	}
	
	if (zeroDC) Normalization::zeroDC_stack(tomogram.stack);
		
	
	const int fc = tomogram.frameCount;
		
	BufferedImage<float> stackAct;
	std::vector<d4Matrix> projAct(fc);
	
	
	const int w1 = w > 0? w : w0 / spacing + 0.5;
	const int h1 = h > 0? h : h0 / spacing + 0.5;
	const int t1 = (int)(thickness/spacing);
	
	std::cout << w1 << "x" << h1 << "x" << t1 << std::endl; 
		
	if (std::abs(spacing - 1.0) < 1e-2)
	{
		projAct = tomogram.proj;
		stackAct = tomogram.stack;
	}
	else
	{
		for (int f = 0; f < fc; f++)
		{
			projAct[f] = tomogram.proj[f] / spacing;
			projAct[f](3,3) = 1.0;
		}
		
		if (std::abs(spacing / stack_spacing - 1.0) > 1e-2)
		{
			std::cout << "resampling image stack... " << std::endl;
			
			stackAct = Resampling::downsampleFiltStack_2D_full(tomogram.stack, spacing / stack_spacing, n_threads);
		}
		else
		{
			stackAct = tomogram.stack;
		}
	}
	
	
	d3Vector orig(x0, y0, z0);
	BufferedImage<float> out(w1, h1, t1);
	
	std::cout << "backprojecting... " << std::endl;
	
	RealSpaceBackprojection::backproject(
			stackAct, projAct, out, n_threads, 
			orig, spacing, RealSpaceBackprojection::Linear, taperRad);
	
	if (weight)
	{
		BufferedImage<float> psf(w1, h1, t1);
		
		RealSpaceBackprojection::backprojectPsf(stackAct, projAct, psf, n_threads, orig, spacing);
		
		Reconstruction::correct3D_RS(out, psf, out, WienerFract, n_threads);
	}
	
	std::cout << "writing output..." << std::endl;
	
	const double samplingRate = tomogram.optics.pixelSize * spacing;

	out.write(outFn, samplingRate);
}

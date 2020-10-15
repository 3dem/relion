
#include <src/jaz/tomography/tomo_stack.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/membrane/blob_3d.h>
#include <src/jaz/membrane/blob_fit_3d.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/volume_converter.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	if (false)
	{
		Image<double> stack;
		stack.read("imod_good/b001ts001.ali:mrc");
		
		BufferedImage<double> stack2(stack.data.xdim, stack.data.ydim, stack.data.zdim);
		
		for (int z = 0; z < stack.data.zdim; z++)
		for (int y = 0; y < stack.data.ydim; y++)
		for (int x = 0; x < stack.data.xdim; x++)
		{
			stack2(x,y,z) = stack(z,y,x);
		}
		
		stack2.write("imod_good/b001ts001_nohead.ali:mrc");
		
		for (int z = 0; z < stack.data.zdim; z++)
		for (int y = 0; y < stack.data.ydim; y++)
		for (int x = 0; x < stack.data.xdim; x++)
		{
			stack(z,y,x) = (x == y)? 1.0 : 0.0;
		}
		
		stack.write("diagonal.ali:mrc");
		
		for (int z = 0; z < stack.data.zdim; z++)
		for (int y = 0; y < stack.data.ydim; y++)
		for (int x = 0; x < stack.data.xdim; x++)
		{
			stack2(x,y,z) = (x == y)? 1.0 : 0.0;
		}
		
		stack2.write("diagonal_nohead.ali:mrc");
		
		return 0;
	}
	
	/*if (false)
	{
		std::string stack_fn = "gridStack.mrc";
		std::string tlt_fn = "test_off2.tlt";
		std::string xf_fn = "testC.xf";
		
		TomoStack<float> ts0 = TomoStack<float>(stack_fn, -1, tlt_fn, xf_fn, "", 1.0, 1.0);
		
		const int w = ts0.images[0].xdim;
		const int h = ts0.images[0].ydim;
		const int d = 81; // field THICKNESS in .com file
					
		BufferedImage<float> dest(w,d,h);
		BufferedImage<float> maskDest(w,d,h);
		
		dest.fill(0.0);
		maskDest.fill(0.0);
		
		d3Vector origin(0.0, -d/2.0 + 1.0, 0.0);
		d3Vector spacing(1.0);
		
		const int n_threads = 6;
		
		RealSpaceBackprojection::backprojectRaw(
			ts0, dest, maskDest, origin, spacing, n_threads,
			RealSpaceBackprojection::Linear,
			2, 2, 1.0);
	
		dest.writeVtk("gridStack_dyn_testC_off2.vtk", d3Vector(0.0, d/2, 0.0), d3Vector(1.0));
	}*/
	
	//if (false)
	{
		const int w = 240;
		const int h = 160;
		const int s = 240;
		const int fc = 5;
		
		BufferedImage<float> testImg(w,h,fc);
		testImg.fill(0.f);
		
		for (int i = 0; i < s; i++)
		{
			for (int d = 1; d <= 3; d++)
			{
				const double e = s / (2 << d);
				const double intensity = 1.0 / d;
				const int qMax = (int) FLOOR(s / e);
				
				for (int q = 0; q <= qMax; q++)
				{
					const int r = ((int)(q * e + 0.5) + s/2) % s;
					
					for (int f = 0; f < fc; f++)
					for (int t = 0; t < s; t++)
					{
						if (r < w && t < h) testImg(r,t,f) += intensity;
						if (r < h && t < w) testImg(t,r,f) += intensity;
					}
				}
			}
		}
		
		testImg.write("gridStack.mrc");
		testImg.write("gridStack.vtk");
		
		return 1;
	}
			
	
	
	
}

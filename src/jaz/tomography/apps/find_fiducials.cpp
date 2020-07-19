#include <src/args.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/tomography/projection/real_backprojection.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/image/resampling.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	std::string tomoSetFn, optFn, outDir;
	double thresh, binning, handedness;
	int num_threads;
	
	IOParser parser;

	try
	{	
		int gen_section = parser.addSection("General refinement options");
		
		tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		thresh = textToDouble(parser.getOption("--d", "Detection threshold", "3.5"));
		binning = textToDouble(parser.getOption("--bin", "Binning level", "16"));
			
		optFn = parser.getOption("--ctf", "Consider a CTF using parameters from the supplied file", "");
		handedness = textToDouble(parser.getOption("--handedness", "Handedness", "1"));
		
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		outDir = parser.getOption("--o", "Output filename pattern");
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	TomogramSet tomogramSet(tomoSetFn);
	
	const int tc = tomogramSet.size();
	
	
	for (int t = 0; t < tc; t++)
	{
		Tomogram tomogram = tomogramSet.loadTomogram(t, true);
		const int fc = tomogram.frameCount;
		
		std::vector<d4Matrix> projAct(fc);
				
		for (int f = 0; f < fc; f++)
		{
			projAct[f] = tomogram.projectionMatrices[f] / binning;
			projAct[f](3,3) = 1.0;
		}
		
		BufferedImage<float> stackAct = Resampling::downsampleFiltStack_2D_full(
					tomogram.stack, binning, num_threads);
		
		const int w1 = tomogram.w0 / binning;
		const int h1 = tomogram.h0 / binning;
		const int d1 = tomogram.d0 / binning;
		
		d3Vector orig(0, 0, 0);
		BufferedImage<float> out(w1, h1, d1);
		
		std::cout << "backprojecting... " << std::endl;
		
		RealSpaceBackprojection::backproject(
				stackAct, projAct, out, num_threads, 
				orig, binning, RealSpaceBackprojection::Linear, 10);
		
		out.write("debug_bp.mrc");
		std::exit(0);
	}
	
	return 0;
	
	/*
	Tomogram tomo
	
	for (int f = 0; f < fc; f++)
	{
		projAct[f] = proj[f] / spacing;
		projAct[f](3,3) = 1.0;
	}
	
	if (std::abs(spacing / stack_spacing - 1.0) > 1e-2)
	{
		std::cout << "resampling image stack... " << std::endl;
		
		stackAct = Resampling::downsampleFiltStack_2D_full(stack, spacing / stack_spacing, n_threads);
	}
	else
	{
		stackAct = stack;
	}
}

	RealSpaceBackprojection::backprojectRaw(
		tomogram, similarities2D, subvolume, subvolumeMask, 
		p - gravis::d3Vector(binning), 
		gravis::d3Vector(1.0), 
		num_threads, 
		RealSpaceBackprojection::Linear,
		20, 20, 10.0);*/
}

#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/args.h>
#include <src/jaz/gravis/t4Matrix.h>

#include <omp.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	IOParser parser;
	
	int n_threads = 1;
	
	double spacing;
	std::string stackFn, projFn, outStackFn, outProjFn;	
	bool zeroDC, FourierCrop;


	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");

	stackFn = parser.getOption("--i", "Tilt sequence image stack (e.g. *.st:mrc)");
	zeroDC = parser.checkOption("--0dc", "Zero the DC component of each frame");
	spacing = textToDouble(parser.getOption("--bin", "Binning (pixel spacing)", "8.0"));
	FourierCrop = parser.checkOption("--F", "Use Fourier cropping to downsample (warning: can produce an inexact aspect ratio)");
	n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
	outStackFn = parser.getOption("--o", "Output stack filename");

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	
	BufferedImage<float> stack;
	stack.read(stackFn);
	
	if (zeroDC) Normalization::zeroDC_stack(stack);

	std::cout << "resampling image stack... " << std::endl;
	
	BufferedImage<float> stackAct;
	
	if (FourierCrop)
	{
		stackAct = Resampling::FourierCrop_fullStack(stack, spacing, n_threads, true);
	}
	else
	{
		stackAct = Resampling::downsampleFiltStack_2D_full(stack, spacing, n_threads);
	}
	
	stackAct.write(outStackFn);
	
	return 0;
}

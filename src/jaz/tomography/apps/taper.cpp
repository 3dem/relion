#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/args.h>
#include <src/jaz/gravis/t4Matrix.h>

#include <omp.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	IOParser parser;
	
	int n_threads = 1;
	
	double width;
	std::string inFn, outFn;

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");

	inFn = parser.getOption("--i", "Image to be tapered");
	width = textToDouble(parser.getOption("--w", "Filter edge width", "10"));
	n_threads = textToInteger(parser.getOption("--j", "Number of threads", "5"));
	outFn = parser.getOption("--o", "Output filename");

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	BufferedImage<float> imgIn;
	imgIn.read(inFn);
	
	Reconstruction::taper(imgIn, width, true, n_threads);
	
	imgIn.write(outFn);
	
	return 0;
}

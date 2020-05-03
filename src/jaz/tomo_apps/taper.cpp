#include <src/jaz/tomo/projection/projection.h>
#include <src/jaz/tomo/extraction.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/tomo/reconstruction.h>
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

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		inFn = parser.getOption("--i", "Image to be tapered");
		width = textToDouble(parser.getOption("--w", "Filter edge width", "10"));
		n_threads = textToInteger(parser.getOption("--j", "Number of threads", "5"));
		outFn = parser.getOption("--o", "Output filename");
			
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	BufferedImage<float> imgIn;
	imgIn.read(inFn);
	
	Reconstruction::taper(imgIn, width, true, n_threads);
	
	imgIn.write(outFn);
	
	return 0;
}

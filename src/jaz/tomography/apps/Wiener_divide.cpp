#include <src/jaz/image/centering.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/image/power_spectrum.h>
#include <src/jaz/tomography/reconstruction.h>
#include <iostream>

using namespace gravis;


int main(int argc, char *argv[])
{
	IOParser parser;
	
	int n_threads;
	
	double SNR, taperWidth;
	std::string dataFn, wghFn, outFn;


	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");

	dataFn = parser.getOption("--i", "Data term");
	wghFn = parser.getOption("--w", "Data term");
	SNR = textToDouble(parser.getOption("--SNR", "Assumed SNR", "0.0001"));
	taperWidth = textToDouble(parser.getOption("--tw", "Taper edge width", "10"));
	n_threads = textToInteger(parser.getOption("--j", "Number of threads", "5"));
	outFn = parser.getOption("--o", "Output filename");

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	BufferedImage<float> data;
	data.read(dataFn);
	
	BufferedImage<float> weight;
	weight.read(wghFn);
	
	const int wd = data.xdim;
	const int hd = data.ydim;
	const int dd = data.zdim;
	
	const int ww = weight.xdim;
	const int hw = weight.ydim;
	const int dw = weight.zdim;
	
	if (hd != hw || dd != dw || !(wd == ww || wd/2 + 1 == ww))
	{
		REPORT_ERROR_STR("incompatible image sizes: " << data.getSizeString() 
						 << " vs. " << weight.getSizeString());
	}
	
	BufferedImage<float> weightHalf;
	
	if (wd == ww)
	{
		weightHalf = Centering::humanFullToFftwHalf(weight);
	}
	else
	{
		weightHalf = weight;
	}
	
	
	BufferedImage<fComplex> dataFS;
	
	FFT::FourierTransform(data, dataFS, FFT::Both);
	
	const int wh = wd/2 + 1;
	
	const double off = 1.0 / SNR;
	
	std::cout << "offset: " << off << std::endl;
	
	dataFS.write("debug_dataFS_0.vtk");
	weight.write("debug_weight.mrc");
	
	for (int z = 0; z < dd; z++)
	for (int y = 0; y < hd; y++)
	for (int x = 0; x < wh; x++)
	{
		dataFS(x,y,z) /= weightHalf(x,y,z) + off;
	}
	
	dataFS.write("debug_dataFS_1.vtk");
	
	FFT::inverseFourierTransform(dataFS, data, FFT::Both);
	
	Reconstruction::taper(data, taperWidth, true, n_threads);
	
	
	data.write(outFn);
	
	return 0;
}

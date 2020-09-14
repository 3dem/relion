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
	
	double slope, d0;
	std::string imgFn, outFn;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		imgFn = parser.getOption("--i", "map");
		const double alpha = 0.5 * textToDouble(parser.getOption("--angle", "opening angle [°]", "2"));
		slope = sin(DEG2RAD(alpha));
		d0 = textToDouble(parser.getOption("--d0", "d0", "1"));		
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
	
	std::cout << "slope = " << slope << std::endl;
	
	BufferedImage<float> img;
	img.read(imgFn);
	
	const int s = img.xdim;
	const int sh = s/2 + 1;
	
	BufferedImage<fComplex> imgFS;
	
	FFT::FourierTransform(img, imgFS);	
	imgFS.write("debug_0.vtk");
	
	
	BufferedImage<float> mask(sh,s,s);
	
	for (int z = 0; z < s; z++)
	for (int y = 0; y < s; y++)
	for (int x = 0; x < sh; x++)
	{
		const double xx = x;
		const double yy = y < s/2? y : y - s;
		const double zz = z < s/2? z : z - s;
		
		const double rho = sqrt(xx*xx + yy*yy);
		const double t = rho / (std::abs(zz) * slope + d0);
		
		const double m = 1.0 - exp(-0.5*t*t);
		
		imgFS(x,y,z) *= m;
		
		mask(x,y,z) = m;
	}
	
	imgFS.write("debug_1.vtk");
	FFT::inverseFourierTransform(imgFS, img);
	
	mask.write("debug_mask.mrc");
	img.write(outFn);
}

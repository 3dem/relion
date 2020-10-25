#include <src/jaz/math/Euler_angles_dynamo.h>
#include <src/macros.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/math/fft.h>
#include <omp.h>

using namespace gravis;


int main(int argc, char *argv[])
{
	if (argc < 3)
	{
		std::cerr << "usage: ramp <in> <out>\n";
		return 1;
	}
	
	std::string fnIn = argv[1];
	std::string fnOut = argv[2];
	
	BufferedImage<double> img;
	img.read(fnIn);
	
	BufferedImage<double> out = ImageFilter::ramp(img);
	
	out.write(fnOut);
	
	return 0;
}

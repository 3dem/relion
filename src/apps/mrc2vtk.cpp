#include <src/image.h>
#include <src/jaz/vtk_helper.h>

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		std::cout << "usage: relion_mrc2vtk X.(mrc/mrcs/tiff/spi)\n		  -> X.vtk\n";
		return 666;
	}

	std::string fn(argv[1]), fn2;

	if (fn.find_last_of('.') != std::string::npos)
	{
		fn2 = fn.substr(0, fn.find_last_of('.')) + ".vtk";
	}
	else
	{
		fn2 = fn + ".vtk";
	}

	Image<RFLOAT> img;
	img.read(fn);

	Image<RFLOAT> imgZ = VtkHelper::allToZ(img);
	VtkHelper::writeVTK(imgZ, fn2);
}

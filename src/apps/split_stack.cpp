
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/volume_converter.h>

#include <src/jaz/gravis/t4Matrix.h>
using namespace gravis;

int main(int argc, char *argv[])
{
	std::string inPath = "./";
	std::string inName = "TS_03";

	Image<RFLOAT> img0;
	img0.read(inPath+inName+".st:mrcs", true);

	Image<RFLOAT> img1(img0.data.xdim, img0.data.ydim, 1, 1);

	for (int i = 0; i < img0.data.ndim; i++)
	{
		SliceHelper::extractStackSlice(img0, img1, i);

		std::stringstream sts;
		sts << i;
		std::string fn;
		sts >> fn;

		img1.write(inPath+"frames/"+inName+"_f"+fn+".mrc");
	}

	return RELION_EXIT_SUCCESS;
}

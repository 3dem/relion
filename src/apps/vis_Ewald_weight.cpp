
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/ctf.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/new_ft.h>
#include <src/jaz/single_particle/noise_helper.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/ctf/delocalisation_helper.h>
#include <src/jaz/image/color_helper.h>

#ifdef HAVE_PNG
#include <src/jaz/gravis/tImage.h>
#endif

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	std::string starFn, outPath;
	int threads, mg_index, part_index, s_cl;
	double mask_rad;
	
	IOParser parser;
	
	try
	{
		parser.setCommandLine(argc, argv);
		
		parser.addSection("General options");
		
		starFn = parser.getOption("--i", "Input particle *.star file");
		mask_rad = textToDouble(parser.getOption("--rad", "Mask radius [A]", "50"));
		
		mg_index = textToInteger(parser.getOption("--m", "Micrograph index", "0"));
		part_index = textToInteger(parser.getOption("--p", "Particle index", "0"));
		
		s_cl = textToInteger(parser.getOption("--s", "Box size (overrides particle file)", "-1"));
						
		threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		
		outPath = parser.getOption("--o", "Output path");

		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

	ObservationModel obsModel;	
	MetaDataTable mdt0;
	
	ObservationModel::loadSafely(starFn, obsModel, mdt0);
	
	std::vector<MetaDataTable> allMdts = StackHelper::splitByMicrographName(mdt0);
	
	const int optGroup = obsModel.getOpticsGroup(allMdts[mg_index], part_index);
	
	const int s = s_cl > 0? s_cl : obsModel.getBoxSize(optGroup);
	const int sh = s/2 + 1;
	
	const double angpix = obsModel.getPixelSize(optGroup);
	
	CTF ctf;
	ctf.readByGroup(allMdts[mg_index], &obsModel, part_index);
	
	std::cout << "drawing A...\n";
	
	Image<RFLOAT> img0(sh,s), img1(sh,s);	
	ctf.applyWeightEwaldSphereCurvature(img0.data, s, s, angpix, 2*mask_rad);
	ctf.applyWeightEwaldSphereCurvature_new(img1.data, s, s, angpix, 2*mask_rad);
	
	Image<RFLOAT> img0Full, img1Full;
	FftwHelper::decenterDouble2D(img0(), img0Full());
	FftwHelper::decenterDouble2D(img1(), img1Full());
	
	VtkHelper::writeVTK(img0Full(), outPath + "_Ewald_weight_0.vtk");
	VtkHelper::writeVTK(img1Full(), outPath + "_Ewald_weight_1.vtk");
	
	return RELION_EXIT_SUCCESS;
}

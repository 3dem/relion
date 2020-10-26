
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

Image<RFLOAT> visualiseBand(
	CTF& ctf, int s, double r0, double r1, double flank, double angpix, 
	std::string outPath, std::string tag, 
	bool delocSupp, double mask_rad, bool writeVtks);

void pasteImage(const Image<RFLOAT>& src, Image<RFLOAT>& dest, int x0, int y0);

int main(int argc, char *argv[])
{
	std::string starFn, outPath;
	int threads, mg_index, part_index, s_cl;
	double band_width, band_step, mask_rad;
	bool writeVtks, delocSupp;
	
	IOParser parser;
	
	try
	{
		parser.setCommandLine(argc, argv);
		
		parser.addSection("General options");
		
		starFn = parser.getOption("--i", "Input particle *.star file");
		delocSupp = parser.checkOption("--ds", "Apply delocalisation suppression");
		mask_rad = textToDouble(parser.getOption("--rad", "Mask radius [Px]", "50"));
		
		mg_index = textToInteger(parser.getOption("--m", "Micrograph index", "0"));
		part_index = textToInteger(parser.getOption("--p", "Particle index", "0"));
		
		s_cl = textToInteger(parser.getOption("--s", "Box size (overrides particle file)", "-1"));
		
		band_width = textToDouble(parser.getOption("--bw", "Width of each frequency band [Px]", "50"));
		band_step = textToDouble(parser.getOption("--bs", "Falloff of frequency bands [Px]", "15"));
		
		writeVtks = parser.checkOption("--write_vtk", "Write VTK files for individual images");
				
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
	
	int bandCount = sh / band_width + 1;
	
	Image<RFLOAT> output(3 * s, (2 + bandCount)*s);
	
	Image<RFLOAT> square = visualiseBand(
		ctf, s, 0, 2*s, band_step, angpix, outPath, "full-square", 
		delocSupp, mask_rad, writeVtks);
	
	Image<RFLOAT> circle = visualiseBand(
		ctf, s, 0, sh, band_step, angpix, outPath, "full-circle", 
		delocSupp, mask_rad, writeVtks);
	
	pasteImage(square, output, 0, 0*s);
	pasteImage(circle, output, 0, 1*s);
	
	for (int b = 0; b < bandCount; b++)
	{
		double r0 = b * band_width;
		double r1 = (b + 1) * band_width;
		
		std::stringstream sts;
		sts << b;
	
		Image<RFLOAT> band = visualiseBand(
			ctf, s, r0, r1, band_step, angpix, outPath, "band_"+sts.str(), 
			delocSupp, mask_rad, writeVtks);
		
		pasteImage(band, output, 0, (b+2)*s);
	}
	
	#ifdef HAVE_PNG
	{
		tImage<dRGB> pngOut(output.data.xdim, output.data.ydim);
		pngOut.fill(dRGB(0.f));
		
		for (int y = 0; y < output.data.ydim; y++)
		for (int x = 0; x < output.data.xdim; x++)
		{
			double c = output(y,x);
			//pngOut(x,y) = fRGB(std::max(c,0.0), c*c, std::max(-c,0.0));
			pngOut(x,y) = ColorHelper::signedToRedBlue(c);
		}
		
		pngOut.writePNG(outPath + "_all.png");
	}
	#endif
	
	return RELION_EXIT_SUCCESS;
}

Image<RFLOAT> visualiseBand(
	CTF& ctf, int s, double r0, double r1, double flank, double angpix, 
	std::string outPath, std::string tag, 
	bool delocSupp, double mask_rad, bool writeVtks)
{
	const int sh = s/2 + 1;
	
	Image<RFLOAT> one(sh,s);
	one.data.initConstant(1);
	
	Image<RFLOAT> mask = FilterHelper::raisedCosEnvRingFreq2D(one, r0, r1, flank);
	
	Image<RFLOAT> ctfImg(sh,s), ctfImgFull(s,s);
	ctf.getFftwImage(ctfImg(), s, s, angpix);
	
	if (delocSupp)
	{
		DelocalisationHelper::maskOutsideBox(ctf, mask_rad, angpix, s, ctfImg(), 0.0, 0.0);
	}
	
	ctfImg.data *= mask.data;
			
	FftwHelper::decenterDouble2D(ctfImg(), ctfImgFull());
	if (writeVtks) VtkHelper::writeVTK(ctfImgFull, outPath + "_" + tag + "_ctf.vtk");
	
	Image<Complex> ctfImgComplex(sh,s);
	
	for (int y = 0; y < s;  y++)
	for (int x = 0; x < sh; x++)
	{
		ctfImgComplex(y,x) = ctfImg(y,x);
	}
	
	Image<RFLOAT> psf, psfFull;
	NewFFT::inverseFourierTransform(ctfImgComplex(), psf(), NewFFT::Both);
	
	FftwHelper::decenterFull(psf(), psfFull());
	if (writeVtks) VtkHelper::writeVTK(psfFull, outPath + "_" + tag + "_psf.vtk");
	
	Image<RFLOAT> deloc = DelocalisationHelper::plotDelocalisation(ctf, mask, angpix);
	
	Image<RFLOAT> delocFull;
	FftwHelper::decenterFull(deloc(), delocFull());
	if (writeVtks) VtkHelper::writeVTK(delocFull, outPath + "_" + tag + "_deloc.vtk");
	
	Image<RFLOAT> out(3*s, s);
	
	pasteImage(ctfImgFull, out, 0, 0);
	pasteImage(FilterHelper::normaliseToUnitIntervalSigned(psfFull), out, s, 0);
	pasteImage(FilterHelper::normaliseToUnitIntervalSigned(delocFull), out, 2*s, 0);
	
	return out;
}

void pasteImage(const Image<RFLOAT>& src, Image<RFLOAT>& dest, int x0, int y0)
{
	const int ws = src.data.xdim;
	const int hs = src.data.ydim;
	
	const int wd = dest.data.xdim;
	const int hd = dest.data.ydim;
	
	for (int y = 0; y < hs; y++)
	for (int x = 0; x < ws; x++)
	{
		const int xx = x0 + x;
		const int yy = y0 + y;
		
		if (xx >= 0 && xx < wd && yy >= 0 && yy < hd)
		{
			dest(yy,xx) = src(y,x);
		}
	}
}

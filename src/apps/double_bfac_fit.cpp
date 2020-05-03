
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/jaz/single_particle/damage_helper.h>
#include <src/jaz/single_particle/image_log.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;
	
	parser.setCommandLine(argc, argv);
	parser.addSection("General options");
	
	std::string fccFn0 = parser.getOption("--fcc0", "FCC 0");
	std::string fccFn1 = parser.getOption("--fcc1", "FCC 1");
	
	int k0 = textToInteger(parser.getOption("--k0", "k0"));
	int k1 = textToInteger(parser.getOption("--k1", "k1"));
	
	double angpix = textToDouble(parser.getOption("--angpix", "pixel size"));
	
	std::string outPath = parser.getOption("--o", "output path");
	
	
	if (parser.checkForErrors()) return RELION_EXIT_FAILURE;
	
	Image<RFLOAT> fcc0, fcc1;
	fcc0.read(fccFn0);
	fcc1.read(fccFn1);
	
	const int w = fcc0.data.xdim;
	const int h0 = fcc0.data.ydim;
	const int h1 = fcc1.data.ydim;	
	const int h = h0 + h1;
	
	Image<RFLOAT> fcc(w,h);
	
	for (int y = 0; y < h0; y++)
	{
		for (int x = 0; x < w; x++)
		{
			fcc(y,x) = fcc0(y,x);
		}
	}
	
	for (int y = 0; y < h1; y++)
	{
		for (int x = 0; x < w; x++)
		{
			fcc(y+h0,x) = fcc1(y,x);
		}
	}
	
	const int sh = w;
	const int s = 2*(sh-1);
	const int fc = h;
	
	std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs
            = DamageHelper::fitBkFactors(fcc, k0, k1);
	
	mktree(outPath + "/");

	Image<RFLOAT> bfacFit = DamageHelper::renderBkFit(bkFacs, sh, fc);
	Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs, sh, fc, true);
	std::vector<Image<RFLOAT>> pixelWeights(2, Image<RFLOAT>(sh,h0));
	
	for (int x = 0; x < sh; x++)
	{
		double sum = 0.0;
		
		for (int i = 0; i < h0; i++)
		{
			sum += bfacFitNoScale(i,x);			
		}
		
		for (int i = 0; i < h0; i++)
		{
			pixelWeights[0](i,x) = bfacFitNoScale(i,x) / sum;
		}
	}
	
	ImageLog::write(pixelWeights[0], outPath + "/pixel_weights_0");
	
	for (int x = 0; x < sh; x++)
	{
		double sum = 0.0;
		
		for (int i = 0; i < h1; i++)
		{
			sum += bfacFitNoScale(i+h0,x);			
		}
		
		for (int i = 0; i < h1; i++)
		{
			pixelWeights[1](i,x) = bfacFitNoScale(i+h0,x) / sum;
		}
	}
	
	ImageLog::write(pixelWeights[1], outPath + "/pixel_weights_1");
			
	ImageLog::write(bfacFit, outPath + "/glob_Bk-fit");
	ImageLog::write(bfacFitNoScale, outPath + "/glob_Bk-fit_noScale");
	ImageLog::write(fcc, outPath + "/glob_Bk-data");

	std::ofstream bfacsDat(outPath + "/Bfac.dat");
	std::ofstream kfacsDat(outPath + "/kfac.dat");
	
	const double cf = 8.0 * angpix * angpix * sh * sh;
	
	double avg0 = 0.0;
	
	for (int i = 0; i < h0; i++)
	{
		double s = bkFacs.first[i].x;
		double b = -cf/(s*s);

		bfacsDat << i << " " << b << std::endl;
		kfacsDat << i << " " << log(bkFacs.first[i].y) << std::endl;
		
		avg0 += b/h0;
	}
	
	bfacsDat << "\n";
	kfacsDat << "\n";
	
	std::ofstream dfacsDat(outPath + "/Dfac.dat");
	
	for (int i = 0; i < sh; i++)
	{
		dfacsDat << i << " " << bkFacs.second[i] << "\n";
		
	}	
	
	dfacsDat << "\n";
	
	double avg1 = 0.0;
	
	for (int i = 0; i < h1; i++)
	{
		double s = bkFacs.first[i+h0].x;
		double b = -cf/(s*s);

		bfacsDat << i << " " << b << std::endl;
		kfacsDat << i << " " << log(bkFacs.first[i+h0].y) << std::endl;
		
		avg1 += b/h1;
	}	

	bfacsDat.close();
	kfacsDat.close();
	
	std::cout << "avg0 = " << avg0 << "\n";
	std::cout << "avg1 = " << avg1 << "\n";
	std::cout << "avg diff = " << (avg1 - avg0) << "\n";
	
	return RELION_EXIT_SUCCESS;
}

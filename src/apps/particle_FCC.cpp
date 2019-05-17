
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/obs_model.h>
#include <src/jaz/image_log.h>
#include <src/jaz/new_ft.h>
#include <src/jaz/noise_helper.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/reference_map.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
    std::string starFn, outPath;
	double minFreqPx;
	bool oppositeHalf, predictCTF;
	int minMG, maxMG, threads;
	
	ReferenceMap reference;

    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");
		
		starFn = parser.getOption("--i", "Input particle *.star file");
		reference.read(parser, argc, argv);
		
		minFreqPx = textToDouble(parser.getOption("--min_freq", "Min. image frequency [px]", "0"));
		
		oppositeHalf = parser.checkOption("--opposite_half", "Correlate with opposite half-set");
		predictCTF = parser.checkOption("--predict_CTF", "Modulate prediction by CTF");
		
		minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
		maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));

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
	
	reference.load(1, false);
	
	const int s = reference.s;
	const int sh = s/2 + 1;
	
	if (maxMG < 0) maxMG = allMdts.size() - 1;
	
	std::vector<double> num(sh, 0.0), denom0(sh, 0.0), denom1(sh, 0.0);
	
	for (int m = 0; m <= maxMG; m++)
	{
		std::vector<Image<Complex>> obs, pred;
		
		int opticsGroup;
		allMdts[m].getValue(EMDL_IMAGE_OPTICS_GROUP, opticsGroup, 0);
		opticsGroup--;
		
		// both defocus_tit and tilt_fit need the same observations
		obs = StackHelper::loadStackFS(allMdts[m], "", threads, true, &obsModel);
		
		pred = reference.predictAll(
			allMdts[m], obsModel, 
			oppositeHalf? ReferenceMap::Opposite : ReferenceMap::Own, 
			threads, predictCTF, true, false);
		
		const int pc = obs.size();
		
		for (int p = 0; p < pc; p++)
		{
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				const double xx = x;
				const double yy = (y + s/2) % s - s/2;
				
				const int r = (int)(sqrt(xx*xx + yy*yy) + 0.5);
				
				if (r >= sh) continue;
				
				const Complex z_pred = pred[p](y,x);
				const Complex z_obs = obs[p](y,x);
				
				num[r] += z_pred.real * z_obs.real + z_pred.imag * z_obs.imag;
				denom0[r] += z_pred.norm();
				denom1[r] += z_obs.norm();
			}
		}
	}
	
	std::ofstream os(outPath+"_FCC.dat");
	
	for (int r = minFreqPx; r < sh; r++)
	{
		double wgh = denom0[r] * denom1[r];
		
		if (wgh > 0.0)
		{
			double fcc = num[r] / sqrt(wgh);
			
			os << r << " " << fcc << "\n";
		}
	}		

	return RELION_EXIT_SUCCESS;
}

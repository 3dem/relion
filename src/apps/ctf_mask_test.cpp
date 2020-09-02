
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/new_ft.h>
#include <src/jaz/single_particle/noise_helper.h>
#include <src/jaz/single_particle/fftw_helper.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
    std::string starFn, outPath;
    int s, threads, mg;
	double rad, step, flankWidth;

    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");
		
		starFn = parser.getOption("--i", "Input particle *.star file");
		s = textToInteger(parser.getOption("--s", "Image size"));
		rad = textToDouble(parser.getOption("--r", "Particle radius"));
		step = textToDouble(parser.getOption("--t", "Frequency step"));
		flankWidth = textToInteger(parser.getOption("--tw", "Filter step width"));
		threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		mg = textToInteger(parser.getOption("--mg", "Micrograph index", "0"));
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
	int opticsGroup = obsModel.getOpticsGroup(allMdts[mg], 0);
	
	const int sh = s/2 + 1;
	
	Image<RFLOAT> ctfImg(sh,s), one(sh,s);
	one.data.initConstant(1.0);
	
	const double angpix = obsModel.getPixelSize(opticsGroup);
	
	CTF ctf;
	ctf.readByGroup(allMdts[mg], &obsModel, 0);
	ctf.getFftwImage(ctfImg(), s, s, angpix);
		
	const int tc = (int)(sh/step + 1.0);
	
	const int maxBin = sh;
	
	std::vector<Image<RFLOAT>> mask(tc+1), psf(tc+1), maskedCTF(tc+1), slopeHistRad(tc+1);
	std::vector<std::vector<double>> slopeHist(tc+1, std::vector<double>(maxBin));
	
	for (int t = 0; t < tc+1; t++)
	{
		const double k0 = t * step;
		const double k1 = (t+1) * step;
		
		if (t < tc)
		{
			mask[t] = FilterHelper::raisedCosEnvRingFreq2D(one, k0, k1, flankWidth);
		}
		else
		{
			mask[t] = one;
		}
		
		Image<Complex> ctfZ(sh,s);		
		Image<RFLOAT> maskedCTF_half(sh,s);
				
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			maskedCTF_half(y,x) = ctfImg(y,x) * mask[t](y,x);
			ctfZ(y,x) = ctfImg(y,x) * mask[t](y,x);
		}
		
		FftwHelper::decenterDouble2D(maskedCTF_half.data, maskedCTF[t].data);
				
		NewFFT::inverseFourierTransform(ctfZ.data, psf[t].data);
		
		const double as = s * angpix;
		
		double minSlope = 100, maxSlope = -100;
		
		for (int y = 0; y < s;  y++)
		for (int x = 0; x < sh; x++)
		{
			double xx = x/as;
			double yy = y < sh? y/as : (y - s)/as;
			
			double slope = ctf.getGammaGrad(xx,yy).length()/(as*PI);
			
			if (mask[t](y,x) >= 0.5)
			{
				if (slope < minSlope) minSlope = slope;
				if (slope > maxSlope) maxSlope = slope;
			}
			
			int si = slope * sh;
			
			if (si < maxBin)
			{
				slopeHist[t][si] += mask[t](y,x);
			}
		}
		
		double maxHist = 0.0;
		
		for (int b = 0; b < maxBin; b++)
		{
			if (slopeHist[t][b] > maxHist)
			{
				maxHist = slopeHist[t][b];
			}
		}
		
		if (maxHist > 0.0)
		{
			for (int b = 0; b < maxBin; b++)
			{
				slopeHist[t][b] /= maxHist;
			}
		}
		
		std::cout << t << ": " << minSlope << " - " << maxSlope << "\n";
		
		slopeHistRad[t] = NoiseHelper::radialMap(slopeHist[t], false);
	}
	
	JazConfig::writeMrc = false;
	JazConfig::writeVtk = true;
	
	ImageLog::write(maskedCTF, outPath+"_maskedCTF");
	ImageLog::write(mask, outPath+"_mask");
	ImageLog::write(psf, outPath+"_psf", ImageLog::CenterXY);
	ImageLog::write(slopeHistRad, outPath+"_slopeHist", ImageLog::CenterXY);
}


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
    std::string starFn, outPath, name;
    int threads, s, optGroup;
	double rad, maxFreqAng, minFreqAng;
	bool allParts;

    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");
		
		starFn = parser.getOption("--i", "Input particle *.star file");
		rad = textToDouble(parser.getOption("--rad", "Particle radius [Å]"));
		optGroup = textToInteger(parser.getOption("--og", "Optics group", "1")) - 1;
		maxFreqAng = textToDouble(parser.getOption("--max_freq", "Max. image frequency [Å] (default is Nyquist)", "-1"));
		minFreqAng = textToDouble(parser.getOption("--min_freq", "Min. image frequency [Å]", "0"));
		name = parser.getOption("--name", "Name of dataset (for the plot)", "");
		allParts = parser.checkOption("--all_part", "Consider all particles, instead of only the first one in each micrograph");
		s = textToInteger(parser.getOption("--s", "Square size for estimation", "256"));
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
	
	const int sh = s/2 + 1;
	
	const double angpix = obsModel.getPixelSize(optGroup);
	
	if (maxFreqAng < 0) maxFreqAng = 2*angpix;
	
	const double r2max = 1.0 / (maxFreqAng * maxFreqAng);
	const double r2min = minFreqAng > 0? 1.0 / (minFreqAng * minFreqAng) : -1;
			
	const int radPx = (int)(rad / angpix + 0.5);
	
	const int maxBin = 5*s;
	
	const double as = s * angpix;
		
	std::vector<double> histCent(maxBin, 0.0);
	std::vector<double> histWorst(maxBin, 0.0);
	
	for (int m = 0; m < allMdts.size(); m++)
	{
		const int pc = allMdts[m].numberOfObjects();
		
		const double mgContrib = allParts? 1.0 : pc;
		const int p_max = allParts? pc : 1;
		
		for (int p = 0; p < p_max; p++)
		{
			int ogp = obsModel.getOpticsGroup(allMdts[m], p);
			
			if (ogp != optGroup) continue;
					
			CTF ctf;
			ctf.readByGroup(allMdts[m], &obsModel, p);
			
			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double xx = x/as;
				double yy = y < sh? y/as : (y - s)/as;
				
				const double r2 = xx*xx + yy*yy;
				
				if (r2 > r2max || r2 < r2min) continue;
				
				
				t2Vector<RFLOAT> delocCent = RFLOAT(1.0 / (2 * angpix * PI)) * ctf.getGammaGrad(xx,yy);
				
				double delocCentVal = delocCent.normLInf();
				
				int sic = (int)(delocCentVal + 0.5);
				if (sic >= maxBin) sic = maxBin - 1;
				
				histCent[sic] += mgContrib;
				
				
				d2Vector delocWorst(std::abs(delocCent.x) + radPx, std::abs(delocCent.y) + radPx);
				
				double delocWorstVal = delocWorst.normLInf();
				
				int siw = (int)(delocWorstVal + 0.5);			
				if (siw >= maxBin) siw = maxBin - 1;
				
				histWorst[siw] += mgContrib;
			}
		}
	}
	
	std::vector<double> histCentCumul(maxBin, 0.0), histWorstCumul(maxBin, 0.0);
	double cumulC = 0.0, cumulW = 0.0;
	
	int first = -1;
	
	for (int b = maxBin-1; b >= 0; b--)
	{
		cumulC += histCent[b];
		histCentCumul[b] = cumulC;
		
		cumulW += histWorst[b];
		histWorstCumul[b] = cumulW;
		
		if (first < 0 && cumulW > 0.0)
		{
			first = b;
		}
	}
	
	if (first < 0)
	{
		std::cerr << "No data found!\n";
		return RELION_EXIT_FAILURE;
	}
	
	
	CPlot2D plot2D("");
	
	std::stringstream ogsts;
	ogsts << (optGroup + 1);
	
	std::string title = "Delocalisation";
	if (name != "") title = title + " for " + name + " (opt. gr. " + ogsts.str() + ")";
	
	std::stringstream pssts;
	pssts << angpix;
	
	std::stringstream frq0sts;
	frq0sts << minFreqAng;
	
	std::stringstream frq1sts;
	frq1sts << maxFreqAng;
	
	title = title + " at " + pssts.str() + " A/px";
	
	if (minFreqAng <= 0)
	{
		title = title + " (up to " + frq1sts.str() + " A)";
	}
	else
	{
		title = title + " (" + frq0sts.str() + " A - " + frq1sts.str() + " A)";
	}
	
	plot2D.SetTitle(title);
	plot2D.SetDrawLegend(true);
	
	CDataSet center;
	center.SetDrawMarker(false);
	center.SetDatasetColor(0.0,0.0,0.0);
	center.SetLineWidth(0.5);
	center.SetDatasetTitle("particle center");
	
	std::stringstream radsts;
	radsts << rad;
	
	CDataSet edge;
	edge.SetDrawMarker(false);
	edge.SetDatasetColor(0.3,0.3,0.6);
	edge.SetLineWidth(0.5);
	edge.SetDatasetTitle("worst periphery point (radius " + radsts.str() + " A)");
	
	for (int i = 0; i < first + radPx + 1; i++)
	{
		if (i < maxBin && i <= first)
		{
			CDataPoint point(2*i, histCentCumul[i]/histCentCumul[0]);
			center.AddDataPoint(point);
		}
		
		if (i < maxBin && i <= first)
		{
			CDataPoint point(2*i, histWorstCumul[i]/histWorstCumul[0]);
			edge.AddDataPoint(point);
		}
	}
	
	plot2D.AddDataSet(center);
	plot2D.AddDataSet(edge);
	
	plot2D.SetXAxisTitle("box size (pixels)");
	plot2D.SetYAxisTitle("fraction of pixels outside of box");
	
	plot2D.OutputPostScriptPlot(outPath+".eps");
	
	return RELION_EXIT_SUCCESS;
}

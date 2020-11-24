#include <src/jaz/tomography/dynamo/catalogue.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/image/power_spectrum.h>
#include <iostream>

using namespace gravis;


int main(int argc, char *argv[])
{
	IOParser parser;
	
	std::string outTag, tomoListFn;
	
	int num_threads;
    double binning;
		

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");

	tomoListFn = parser.getOption("--t", "Tomogram list");
	num_threads = textToInteger(parser.getOption("--j", "Number of threads", "6"));
	binning = textToDouble(parser.getOption("--bin", "Binning factor", "8"));

	outTag = parser.getOption("--o", "Output filename pattern");

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	TomoList tomoList(tomoListFn);
	
	std::vector<long int> tss = tomoList.getKnownIndices();
	
	BufferedImage<float> ts0;
	ts0.read(tomoList.getTiltSeriesFilename(tss[0]));
	
	const int w0 = ts0.xdim;
	const int h0 = ts0.ydim;
	const int fc0 = ts0.zdim;
	
	const int wb = w0 / binning;
	const int hb = h0 / binning;
	
	const long int tc = tss.size();
	
	BufferedImage<float> cat(wb,hb,tc);
	
	std::ofstream indexFile(outTag+"_indices.dat");
	
	for (long int t = 0; t < tc; t++)
	{
		const int tt = tss[t];
		indexFile << t << ' ' << tt << '\n';
	}
	
	indexFile.close();
	
	#pragma omp parallel for num_threads(num_threads)
	for (long int t = 0; t < tc; t++)
	{
		const int tt = tss[t];
		
		std::cout << t << '/' << tc << " (" << tt << ')' << std::endl;
				
		BufferedImage<float> ts;
		ts.read(tomoList.getTiltSeriesFilename(tt));
		
		const int fc = ts.zdim;
		
		BufferedImage<float> slice = NewStackHelper::extractSliceZ(ts, fc/2);
		
		BufferedImage<float> sliceBinned = Resampling::downsampleFilt_2D_full(
				slice, wb, hb);
		
		NewStackHelper::insertSliceZ(sliceBinned, cat, t);
	}
	
	cat.write(outTag+"_catalogue.mrc");
}

#include <src/args.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/stack_helper.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/tomography/tomogram_set.h>



int main(int argc, char *argv[])
{
	IOParser parser;
	parser.setCommandLine(argc, argv);
	
	std::string tomoSetFn, optFn, outDir;
	int boxSize;
	bool normalise;
	
	try
	{
		int gen_section = parser.addSection("General options");
		
		tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		optFn = parser.getOption("--ctf", "Consider a CTF using parameters from the supplied file");
		boxSize = textToInteger(parser.getOption("--b", "Box size"));
		normalise = parser.checkOption("--nrm", "Scale each column to make its sum 1");
		outDir = parser.getOption("--o", "Output directory");
		
		if (parser.checkForErrors()) 
		{
			parser.writeUsage(std::cout);
			std::exit(-1);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	TomogramSet tomogramSet(tomoSetFn);
	const int tc = tomogramSet.size();
	
	if (outDir[outDir.length()-1] != '/')
	{
		outDir = outDir + "/";
	}
	
	int res = system(("mkdir -p "+outDir).c_str());
		
	for (int t = 0; t < tc; t++)
	{		
		Tomogram tomogram = tomogramSet.loadTomogram(t, false);
		const int fc = tomogram.frameCount;
		const int s = boxSize;
		const int sh = s/2 + 1;
		
		BufferedImage<float> doseWeights = tomogram.computeDoseWeight(s, 1.0);
		
		BufferedImage<float> out(sh, fc);
		
		for (int f = 0; f < fc; f++)
		for (int x = 0; x < sh; x++)
		{
			out(x,f) = doseWeights(x,0,f);
		}
		
		if (normalise)
		{
			for (int x = 0; x < sh; x++)
			{
				double sum(0.0);
				
				for (int f = 0; f < fc; f++)
				{
					sum += out(x,f);
				}
				
				if (sum != 0.0)
				{
					for (int f = 0; f < fc; f++)
					{
						out(x,f) /= sum;
					}
				}
			}
		}
		
		out.write(outDir + ZIO::itoa(t) + ".mrc");
	}
	
	return 0;
}

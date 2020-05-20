#include <src/args.h>
#include <src/jaz/tomography/programs/substack.h>

int main(int argc, char *argv[])
{
	IOParser parser;
	
	SubstackProgram sp;
		
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		sp.catFn = parser.getOption("--i", "Input particle set");
		sp.tomoListFn = parser.getOption("--t", "Tomogram list", "tomolist.txt");
		sp.boxSize = textToInteger(parser.getOption("--b", "Box size", "512"));
        sp.binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));
		
		sp.do_center = !parser.checkOption("--no_center", "Do not subtract the mean from the voxel values");		
		
		sp.diag = parser.checkOption("--diag", "Write out diagnostic information");
		
		sp.num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		sp.outTag = parser.getOption("--o", "Output filename pattern");
		
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	sp.run();
		
	return 0;
}

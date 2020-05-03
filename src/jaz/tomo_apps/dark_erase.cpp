#include <src/args.h>
#include <src/jaz/tomo_programs/dark_erase.h>

int main(int argc, char *argv[])
{
	IOParser parser;
	
	DarkEraseProgram dep;	
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		dep.stackFn = parser.getOption("--ts", "Tilt sequence image stack (e.g. *.st:mrc)");
		
		dep.thresh = textToDouble(parser.getOption("--th", "Threshold value (in sigma)", "-2.5"));
		dep.rad = textToDouble(parser.getOption("--r", "Fill radius (in pixels)", "32.0"));
		
		dep.writeNormalized = parser.checkOption("--vis", "Only write out a normalized stack to assist in finding an optimal threshold value");
		dep.diag = parser.checkOption("--diag", "Write out diagnostic data");
		
		dep.num_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		
		dep.outFn = parser.getOption("--o", "Output filename");
				
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	dep.run();
	
	return 0;
}

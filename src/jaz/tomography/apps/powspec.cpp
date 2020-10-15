#include <src/args.h>
#include <src/jaz/tomography/programs/powspec.h>


int main(int argc, char *argv[])
{
	IOParser parser;
	
	PowspecProgram pp;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		pp.stackFn = parser.getOption("--ts", "Tilt stack(s) (e.g. 'a.mrc,b.mrc')");	
		pp.separate = parser.checkOption("--sep", "Compute separate spectra for all frames");
		pp.res = textToInteger(parser.getOption("--r", "Resolution (in pixels)", "512"));
		pp.pixSize = textToDouble(parser.getOption("--angpix", "Pixel size in Å", "1.0"));
		pp.ratio = textToInteger(parser.getOption("--ratio", "Compute the ratio against this frame", "-1"));
		
		pp.outFn = parser.getOption("--o", "Output filename");
				
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	pp.run();
		
	return 0;
}

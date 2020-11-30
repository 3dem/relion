#include <src/args.h>
#include <src/jaz/tomography/programs/substack.h>

int main(int argc, char *argv[])
{
	IOParser parser;
	
	SubstackProgram sp;

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");

	sp.particlesFn = parser.getOption("--i", "Input particle set");
	sp.tomoListFn = parser.getOption("--t", "Tomogram list", "tomolist.txt");
	sp.boxSize = textToInteger(parser.getOption("--b", "Box size"));
	sp.binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));

	sp.do_center = !parser.checkOption("--no_center", "Do not subtract the mean from the voxel values");

	sp.diag = parser.checkOption("--diag", "Write out diagnostic information");

	sp.num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
	sp.outTag = parser.getOption("--o", "Output filename pattern");

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	sp.run();
		
	return 0;
}

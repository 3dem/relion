#include <src/args.h>
#include <src/jaz/tomography/programs/tomo_ctf.h>



int main(int argc, char *argv[])
{
	IOParser parser;
	
	TomoCtfProgram tcp;	

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");

	tcp.stackFn = parser.getOption("--ts", "Tilt sequence image stack (e.g. *.st:mrc)");
	tcp.projFn = parser.getOption("--p", "Tilt projections file (*.tlt)");

	tcp.pixSize = textToDouble(parser.getOption("--angpix", "Pixel size in Å", "1.0"));
	tcp.voltage = textToDouble(parser.getOption("--kV", "Voltage in kV", "1.0"));
	tcp.Cs = textToDouble(parser.getOption("--Cs", "Spherical aberration in mm", "2.7"));

	tcp.n_threads = textToInteger(parser.getOption("--j", "Number of threads", "5"));

	tcp.outFn = parser.getOption("--o", "Output filename");

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	tcp.run();
	
	return 0;
}

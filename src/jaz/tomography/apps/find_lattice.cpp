#include <src/args.h>
#include <src/jaz/tomography/programs/find_lattice.h>


int main(int argc, char *argv[])
{
	IOParser parser;
	
	FindLatticeProgram flp;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		flp.tomoFn = parser.getOption("--t", "Tomogram");
		flp.spacing_ang = textToDouble(parser.getOption("--d", "Lattice spacing [Å]"));
		flp.angpix = textToDouble(parser.getOption("--angpix", "Pixel size [Å]"));
		flp.filter_width = textToDouble(parser.getOption("--f", "Filter width [~Px]", "10.0"));
		flp.minValue = textToDouble(parser.getOption("--mv", "Minimal abs. value of bandpass filtered image", "0.333"));
		flp.minDensity = textToDouble(parser.getOption("--md", "Minimal lattice density", "0.0"));
		flp.taper = textToDouble(parser.getOption("--td", "Tapering distance [Px]", "0.0"));
		flp.overtones = textToInteger(parser.getOption("--ot", "Number of overtones", "0"));
		flp.noise = parser.checkOption("--noise", "Replace tomogram contents by white noise");
				
		flp.n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		
		flp.outFn = parser.getOption("--o", "Output filename");
				
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	flp.run();
		
	return 0;
}

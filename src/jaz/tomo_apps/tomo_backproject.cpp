#include <src/args.h>
#include <src/jaz/tomo_programs/tomo_backproject.h>


int main(int argc, char *argv[])
{
	IOParser parser;	
	
	TomoBackprojectProgram tbp;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		tbp.tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		tbp.tomoIndex = textToInteger(parser.getOption("--ti", "Tomogram index", "0"));
		
		tbp.weight = parser.checkOption("--wg", "Perform weighting in Fourier space (using a Wiener filter)");
		tbp.WienerFract = textToDouble(parser.getOption("--SNR", "SNR assumed by the Wiener filter", "0.1"));
		
		tbp.zeroDC = parser.checkOption("--0dc", "Zero the DC component of each frame");
		
		tbp.taperRad = textToDouble(parser.getOption("--td", "Tapering distance", "0.0"));
	
        tbp.thickness = textToInteger(parser.getOption("--th", "Thickness (read from .tlt file by default)", "-1"));
		
		tbp.x0 = textToDouble(parser.getOption("--x0", "X origin", "1.0"));
		tbp.y0 = textToDouble(parser.getOption("--y0", "Y origin", "1.0"));
		tbp.z0 = textToDouble(parser.getOption("--z0", "Z origin", "1.0"));
		
		tbp.w = textToInteger(parser.getOption("--w", "Width",  "-1.0"));
		tbp.h = textToInteger(parser.getOption("--h", "Height", "-1.0"));
		
		tbp.spacing = textToDouble(parser.getOption("--bin", "Binning (pixel spacing)", "8.0"));
		tbp.stack_spacing = textToDouble(parser.getOption("--stack_bin", "Binning level of the stack", "1.0"));	
		tbp.n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));
		
		tbp.outFn = parser.getOption("--o", "Output filename");
				
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	tbp.run();
	
	return 0;
}

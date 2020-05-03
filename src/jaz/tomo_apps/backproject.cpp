#include <src/args.h>
#include <src/jaz/tomo_programs/backproject.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
	IOParser parser;
	
	BackprojectProgram bp;
	
	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");
		
		bp.catFn = parser.getOption("--i", "Input particle set");
		bp.tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		bp.boxSize = textToInteger(parser.getOption("--b", "Box size", "100"));
		bp.cropSize = textToInteger(parser.getOption("--crop", "Size of (additionally output) cropped image", "-1"));
		bp.binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));
		bp.taper = textToDouble(parser.getOption("--taper", "Taper against the sphere by this number of pixels (only if cropping)", "10"));
		bp.WienerFract = textToDouble(parser.getOption("--SNR", "SNR assumed by the Wiener filter", "0.1"));
		bp.symmName = parser.getOption("--sym", "Symmetry group", "C1");
				
		bp.max_mem_GB = textToInteger(parser.getOption("--mem", "Max. amount of memory to use for accumulation (--j_out will be reduced)", "-1"));
				
		bp.diag = parser.checkOption("--diag", "Write out diagnostic information");
		bp.no_subpix_off = parser.checkOption("--nso", "No subpixel offset (debugging)");
		
		bp.motFn = parser.getOption("--mot", "Particle trajectories", "");
		
		bp.num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		bp.inner_threads = textToInteger(parser.getOption("--j_in", "Number of inner threads (slower, needs less memory)", "3"));
		bp.outer_threads = textToInteger(parser.getOption("--j_out", "Number of outer threads (faster, needs more memory)", "2"));
		
		bp.no_reconstruction = parser.checkOption("--no_recon", "Do not reconstruct the volume, only backproject (for benchmarking purposes)");
		
		bp.outTag = parser.getOption("--o", "Output filename pattern");
		
		Log::readParams(parser);
		
		parser.checkForErrors();
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	bp.run();
	
	return 0;
}

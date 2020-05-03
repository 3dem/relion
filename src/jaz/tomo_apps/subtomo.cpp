#include <src/args.h>
#include <src/jaz/tomo_programs/subtomo.h>
#include <src/jaz/util/log.h>


int main(int argc, char *argv[])
{
    IOParser parser;

    SubtomoProgram sp;

    try
    {
        parser.setCommandLine(argc, argv);
        int gen_section = parser.addSection("General options");

        sp.catFn = parser.getOption("--i", "Catalogue .tbl or .star file");
        sp.tomoSetFn = parser.getOption("--t", "Tomogram set", "tomograms.star");
		sp.boxSize = textToInteger(parser.getOption("--b", "Projection box size", "100"));
		sp.cropSize = textToInteger(parser.getOption("--crop", "Output box size", "-1"));
        sp.binning = textToDouble(parser.getOption("--bin", "Binning factor", "1"));
		sp.write_multiplicity = parser.checkOption("--multi", "Write out multiplicity volumes");		
        sp.do_rotate = parser.checkOption("--rot", "Rotate the particles according to current angles");
        sp.WienerFract = textToDouble(parser.getOption("--SNR", "SNR assumed by the Wiener filter (exaggerated for subtomograms!)", "1.0"));
		
		sp.do_cone_weight = parser.checkOption("--cone_weight", "Weight down a double cone along Z");		
		const double alpha = 0.5 * textToDouble(parser.getOption("--cone_angle", "Opening angle of the cone in degrees", "10"));
		sp.cone_slope = sin(DEG2RAD(alpha));
		sp.cone_sig0 = textToDouble(parser.getOption("--cone_sig0", "Cone width at Z = 0", "2"));	
				
        sp.taper = textToDouble(parser.getOption("--taper", "Taper against the sphere by this number of pixels", "5"));
		sp.env_sigma = textToDouble(parser.getOption("--env", "Sigma of a Gaussian envelope applied before cropping", "-1"));
		
        sp.do_center = !parser.checkOption("--no_center", "Do not subtract the mean from the voxel values");
		
		sp.flip_value = !parser.checkOption("--no_ic", "Do not invert contrast (keep particles dark)");
		sp.write_combined = !parser.checkOption("--no_comb", "Do not write the concatenated CTF-multiplicity image");
		sp.write_ctf = parser.checkOption("--ctf", "Write 3D CTFs");
		sp.write_divided = parser.checkOption("--div", "Write CTF-corrected subtomograms");
		sp.write_normalised = parser.checkOption("--nrm", "Write multiplicity-normalised subtomograms");
		
		sp.motFn = parser.getOption("--mot", "Particle trajectories", "");
		
        sp.diag = parser.checkOption("--diag", "Write out diagnostic information");

        sp.num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
        sp.outTag = parser.getOption("--o", "Output filename pattern");
		
		Log::readParams(parser);

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

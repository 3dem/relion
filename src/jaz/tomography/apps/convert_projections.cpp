#include <src/args.h>
#include <src/jaz/tomography/programs/convert_projections.h>

int main(int argc, char *argv[])
{
	IOParser parser;
	
	ConvertProjectionsProgram cpp;
			

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");

	cpp.newstComFn = parser.getOption("--nc", "Input command to IMOD's newstack", "newst.com");
	cpp.tltComFn = parser.getOption("--tc", "Input command to IMOD's tilt", "tilt.com");

	cpp.tsFn = parser.getOption("--ts", "Original tilt series (only size needed, unless views have to be excluded)");

	cpp.thicknessCmd = textToDouble(parser.getOption(
		"--thick",
		"Thickness of original IMOD tomogram (overrides the value in tilt.com)",
		"-1.0"));

	cpp.offset3Dx = textToDouble(parser.getOption("--offx", "3D offset, X", "0.0"));
	cpp.offset3Dy = textToDouble(parser.getOption("--offy", "3D offset, Y", "0.0"));
	cpp.offset3Dz = textToDouble(parser.getOption("--offz", "3D offset, Z", "0.0"));

	cpp.flipYZ = parser.checkOption("--flipYZ", "Interchange the Y and Z coordinates");
	cpp.flipZ = parser.checkOption("--flipZ", "Change the sign of the Z coordinate");
	cpp.flipAngles = parser.checkOption("--flipAng", "Change the sign of all tilt angles");

	cpp.ali = parser.checkOption("--ali", "Map to aligned stack (.ali)");
	cpp.aliSize = parser.checkOption("--ali_size", "Use the size indicated in newst.com");

	cpp.n_threads = textToInteger(parser.getOption("--j", "Number of threads", "1"));

	cpp.outFn = parser.getOption("--o", "Output filename");
	cpp.outFnCrop = parser.getOption("--oc",
		"Output filename for culled stack (with frames excluded)",
		"");

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}
	
	cpp.run();
	
	return 0;	
}

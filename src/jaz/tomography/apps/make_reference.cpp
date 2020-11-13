#include <src/args.h>
#include <src/jaz/tomography/tomolist.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/optimisation_set.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	IOParser parser;

	OptimisationSet os;
	std::string reconstructionPath, evalMaskFn, refMaskFn, outDir;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		reconstructionPath = parser.getOption("--rec", "Reconstruction path (--o argument to reconstruct_particle)");
		evalMaskFn = parser.getOption("--mask", "Mask to be used for FSC computation");
		refMaskFn = parser.getOption("--opt_mask", "Mask to be used for subsequent optimisation procedures (optional)", "");

		outDir = parser.getOption("--o", "Output directory");

		os.read(
			parser,
			true,           // optimisation set
			true,   false,  // particles
			true,   false,  // tomograms
			true,   false,  // trajectories
			true,   false,  // manifolds
			false,  false); // reference

		if (parser.checkForErrors()) std::exit(-1);
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	if (outDir[outDir.length()-1] != '/')
	{
		outDir = outDir + "/";
	}

	int res = system(("mkdir -p "+outDir).c_str());

	{
		std::ofstream ofs(outDir+"/note.txt");

		ofs << "Command:\n\n";

		for (int i = 0; i < argc; i++)
		{
			ofs << argv[i] << ' ';
		}

		ofs << '\n';
	}

	res = system(("mkdir -p "+outDir+"PostProcess").c_str());

	res = system(("relion_postprocess --i " + reconstructionPath + "_half1.mrc --mask "
				  + evalMaskFn+" --o " + outDir + "PostProcess/post").c_str());

	os.refMap1 = reconstructionPath + "_half1.mrc";
	os.refMap2 = reconstructionPath + "_half2.mrc";
	os.refFSC = outDir + "PostProcess/post.star";

	if (refMaskFn != "")
	{
		os.refMask = refMaskFn;
	}

	os.write(outDir+"optimisation_set.star");

	return 0;
}

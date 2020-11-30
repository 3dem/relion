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
	try
	{
		IOParser parser;

		OptimisationSet os;
		std::string reconstructionPath, evalMaskFn, refMaskFn, outDir;


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

		if (parser.checkForErrors())
		{
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
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

		os.refMap1 = reconstructionPath + "half1.mrc";
		os.refMap2 = reconstructionPath + "half2.mrc";
		os.refFSC = outDir + "PostProcess/postprocess.star";
		if (refMaskFn != "")
		{
			os.refMask = refMaskFn;
		}

		res = system(("mkdir -p "+outDir+"PostProcess").c_str());

		res = system(("relion_postprocess --i " + os.refMap1 + " --mask "
					  + evalMaskFn+" --o " + outDir + "PostProcess/postprocess").c_str());

		os.write(outDir+"optimisation_set.star");
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

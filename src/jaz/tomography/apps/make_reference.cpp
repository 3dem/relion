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

		reconstructionPath = parser.getOption("--rec", "Reconstruction path (--o argument to reconstruct_particle. Optional)", "");
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

		outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);
		ZIO::makeDir(outDir+"PostProcess");

		if (reconstructionPath.length() > 0)
		{
			if (reconstructionPath[reconstructionPath.length()-1] != '/')
			{
				reconstructionPath = reconstructionPath + "/";
			}

			os.refMap1 = reconstructionPath + "half1.mrc";
			os.refMap2 = reconstructionPath + "half2.mrc";
		}

		int res = system(("relion_postprocess --i " + os.refMap1 + " --mask "
					  + evalMaskFn+" --o " + outDir + "PostProcess/postprocess").c_str());

		if (res != RELION_EXIT_SUCCESS) return res;

		os.refFSC = outDir + "PostProcess/postprocess.star";

		if (refMaskFn != "")
		{
			os.refMask = refMaskFn;
		}

		os.write(outDir+"optimisation_set.star");
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

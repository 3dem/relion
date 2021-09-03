#include <src/args.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/optimisation_set.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/util/log.h>
#include <src/args.h>

#include <omp.h>

using namespace gravis;


void run(int argc, char *argv[]);

int main(int argc, char *argv[])
{
	try
	{
		run(argc, argv);
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

void run(int argc, char *argv[])
{
	OptimisationSet optimisationSet;

	IOParser parser;

	parser.setCommandLine(argc, argv);

	optimisationSet.read(
		parser,
		true,             // optimisation set
		false,   false,   // particles
		true,    true,    // tomograms
		true,	 false,   // trajectories
		false,   false,   // manifolds
		false,   false);  // reference

	const int gen_section = parser.addSection("General options");
	const double slope = textToDouble(parser.getOption("--ds", "Defocus slope", "1.0"));
	std::string outDir = parser.getOption("--o", "Output directory");

	Log::readParams(parser);

	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	outDir = ZIO::prepareTomoOutputDirectory(outDir, argc, argv);


	TomogramSet tomogram_set(optimisationSet.tomograms);

	const int tc = tomogram_set.size();

	for (int t = 0; t < tc; t++)
	{
		tomogram_set.globalTable.setValue(EMDL_TOMO_DEFOCUS_SLOPE, slope, t);
	}

	tomogram_set.write(outDir + "tomograms.star");
	optimisationSet.tomograms = outDir + "tomograms.star";

	optimisationSet.write(outDir + "optimisation_set.star");
}


#include <string>
#include <src/metadata_table.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/args.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <omp.h>


using namespace gravis;


int main(int argc, char *argv[])
{
	std::string particlesFn, outDir;	
	
	IOParser parser;
	
	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input file (e.g. run_it023_data.star)");
		outDir = parser.getOption("--o", "Output directory");

		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	outDir = ZIO::makeOutputDir(outDir);

	ObservationModel obs_model;
	MetaDataTable particles_table;
	ObservationModel::loadSafely(particlesFn, obs_model, particles_table);
	
	
	std::vector<MetaDataTable> particles_by_micrograph = StackHelper::splitByMicrographName(particles_table);
	
	const int micrograph_count = particles_by_micrograph.size();
	
	
	MetaDataTable output_table;
	
	for (int m = 0; m < micrograph_count; m++)
	{
		output_table.addObject();
		
		output_table.setValue(
			EMDL_MICROGRAPH_ID, 
			m,
			m);
		
		output_table.setValue(
			EMDL_MICROGRAPH_NAME, 
			particles_by_micrograph[m].getString(EMDL_MICROGRAPH_NAME, 0),
			m);
	}

	output_table.write(outDir + "micrographs.star");
	
	return 0;
}

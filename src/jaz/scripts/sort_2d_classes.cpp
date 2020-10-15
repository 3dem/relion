#include <string>
#include <src/metadata_table.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/class_helper.h>
#include <src/args.h>
#include <src/metadata_table.h>
#include <src/ctf.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/index_sort.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/translation.h>
#include <src/jaz/image/tapering.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/fft.h>
#include <omp.h>


using namespace gravis;


int main(int argc, char *argv[])
{
	std::string particles_filename, class_averages_filename, outDir;
	int num_threads;
	
	IOParser parser;
	
	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particles_filename = parser.getOption("--i", "Input file (e.g. run_it023_data.star)");
		class_averages_filename = parser.getOption("--ca", "Class averages stack");
		
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
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
	ObservationModel::loadSafely(particles_filename, obs_model, particles_table);

	const int class_count = ClassHelper::countClasses(particles_table);
	std::vector<int> particle_count = ClassHelper::getClassSizes(particles_table, class_count);
	std::vector<int> order = ClassHelper::sortByAscendingFrequency(particle_count);
	        
	BufferedImage<float> class_averages;
	class_averages.read(class_averages_filename);
	
	std::cout << class_averages.getSizeString() << std::endl;
	
	const int s = class_averages.xdim;
	
	BufferedImage<float> output_stack(s, s, class_count);
	
	std::ofstream class_size_file(outDir+"class_sizes.txt");
	
	int total = 0;
	
	for (int i = 0; i < class_count; i++)
	{
		const int j = order[i];
		
		output_stack.getSliceRef(i).copyFrom(class_averages.getSliceRef(j));
		class_size_file << i << " (old " << j << ") " << particle_count[j] << '\n';
		
		total += particle_count[j];
	}
	
	class_size_file << "\ntotal: " << total << '\n';
	
	output_stack.write(outDir+"sorted_classes.mrc");
	
	return 0;
}


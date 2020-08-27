#include <string>
#include <src/metadata_table.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/args.h>
#include <src/metadata_table.h>
#include <src/ctf.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/index_sort.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/translation.h>
#include <src/jaz/image/tapering.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/image/local_extrema.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/membrane/blob_fit_2d.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <omp.h>


using namespace gravis;


int main(int argc, char *argv[])
{
	std::string particlesFn, class_averages_filename, outDir;
	int num_threads, num_frequencies, max_class, min_MG, max_MG, max_iterations;
	double radius, edge_padding;
	bool diag;


	IOParser parser;

	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input STAR file output by align_2d_classes");
		class_averages_filename = parser.getOption("--ca", "Class averages stack");
		radius = textToDouble(parser.getOption("--r", "Average blob radius (in bin-1 pixels)", "600"));
		edge_padding = textToDouble(parser.getOption("--pad", "Edge padding (in bin-1 pixels)", "20"));
		num_frequencies = textToInteger(parser.getOption("--f", "Number of 1D frequencies used to model the blob", "10"));

		max_class = textToInteger(parser.getOption("--max_class", "Last 2D class to consider", "25"));
		min_MG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
		max_MG = textToInteger(parser.getOption("--max_MG", "Last micrograph index", "-1"));

		max_iterations = textToInteger(parser.getOption("--it", "Maximal number of iterations", "1000"));
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		diag = parser.checkOption("--diag", "Write out diagnostic information");
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

	BufferedImage<float> class_averages;
	class_averages.read(class_averages_filename);


	const int w = class_averages.xdim;
	const int h = class_averages.ydim;
	const int num_classes = class_averages.zdim;

	if (max_class >= num_classes)
	{
		Log::warn(
			"Only "+ZIO::itoa(num_classes)+" 2D classes found in "+class_averages_filename);

		max_class = num_classes - 1;
	}


	BufferedImage<float> blob_fits(w,h,max_class+1);
	std::vector<std::vector<double>> all_optimal_parameters(max_class+1);

	for (int class_id = 0; class_id <= max_class; class_id++)
	{
		RawImage<float> slice = class_averages.getSliceRef(class_id);
		const d2Vector initial_centre(w/2, h/2 - radius);

		BlobFit2D blob_fit(
			slice, initial_centre, num_frequencies, radius, h, 0, num_threads);

		for (int y = 0; y < blob_fit.weight.ydim; y++)
		for (int x = 0; x < blob_fit.weight.xdim; x++)
		{
			const double dx = x - w/2;
			const double dy = y - h/2;
			const double r = sqrt(dx*dx + dy*dy);

			if (r < w/2 - edge_padding)
			{
				blob_fit.weight(x,y) = 1;
			}
			else
			{
				blob_fit.weight(x,y) = 0;
			}
		}


		std::vector<double> intial_parameters(num_frequencies + 2, 0.0);
		intial_parameters[0] = initial_centre.x;
		intial_parameters[1] = initial_centre.y;

		std::vector<double> optimal_parameters = NelderMead::optimize(
				intial_parameters, blob_fit, 2, 0.001, max_iterations, 1.0, 2.0, 0.5, 0.5, true);

		all_optimal_parameters[class_id] = optimal_parameters;

		std::cout << class_id << ":  "
			<< intial_parameters[0] << ", " << intial_parameters[1] << " => "
			<< optimal_parameters[0] << ", " << optimal_parameters[1] << std::endl;


		Blob2D blob(optimal_parameters, radius + 2*h);

		std::vector<double> radial_average = blob.radialAverage(slice, blob_fit.weight);
		BufferedImage<float> projection = blob.radialAverageProjection(slice, radial_average);

		blob_fits.copySliceFrom(class_id, projection);
	}

	blob_fits.write(outDir+"blob_fits.mrc");

	return 0;
}

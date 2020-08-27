#include <src/args.h>
#include <src/metadata_table.h>
#include <src/ctf.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/translation.h>
#include <src/jaz/image/tapering.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <omp.h>


using namespace gravis;


int main(int argc, char *argv[])
{
	std::string particlesFn, outDir;
	double margin, SNR;
	int num_threads;


	IOParser parser;

	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input file (e.g. run_it023_data.star)", "");
		SNR = textToDouble(parser.getOption("--SNR", "Assumed signal-to-noise ratio", "0.1"));
		margin = textToDouble(parser.getOption("--m", "Margin around the particle [Px]", "20"));
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
	ObservationModel::loadSafely(particlesFn, obs_model, particles_table);

	int max_class = -1;

	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(EMDL_PARTICLE_CLASS, p);

		if (class_id > max_class) max_class = class_id;
	}

	const int class_count = max_class + 1;

	if (class_count == 1)
	{
		Log::print("1 class found");
	}
	else
	{
		Log::print(ZIO::itoa(class_count)+" classes found");
	}


	std::vector<int> class_size(class_count, 0);

	for (long int p = 0; p < particles_table.numberOfObjects(); p++)
	{
		const int class_id = particles_table.getIntMinusOne(
					EMDL_PARTICLE_CLASS, p);

		class_size[class_id]++;
	}



	const int box_size = obs_model.getBoxSize(0);
	const double pixel_size = obs_model.getPixelSize(0);



	BufferedImage<dComplex> data(box_size / 2 + 1, box_size, num_threads * class_count);
	data.fill(dComplex(0.0, 0.0));

	BufferedImage<double> weight(box_size / 2 + 1, box_size, num_threads * class_count);
	weight.fill(0.0);


	std::vector<MetaDataTable> particles_by_micrograph = StackHelper::splitByMicrographName(particles_table);

	const int micrograph_count = particles_by_micrograph.size();

	for (int micrograph_id = 0; micrograph_id < micrograph_count; micrograph_id++)
	{
		Log::print("Micrograph " + ZIO::itoa(micrograph_id+1));

		const MetaDataTable& particles_in_micrograph = particles_by_micrograph[micrograph_id];
		const int particle_count = particles_in_micrograph.numberOfObjects();


		BufferedImage<float> micrograph;
		float mean_value = 0, std_dev = 1;
		double micrograph_pixel_size = pixel_size;
		double extraction_scale = 1;


		const std::string micrograph_filename = particles_in_micrograph.getString(EMDL_MICROGRAPH_NAME, 0);
		micrograph_pixel_size = ImageFileHelper::getSamplingRate(micrograph_filename);
		extraction_scale = pixel_size / micrograph_pixel_size;
		extraction_scale = std::round(1000 * extraction_scale)/1000.0;


		micrograph.read(micrograph_filename);

		mean_value = Normalization::computeMean(micrograph);
		std_dev = sqrt(Normalization::computeVariance(micrograph, mean_value));

		const i2Vector micrograph_size(micrograph.xdim, micrograph.ydim);


		BufferedImage<float>
				out_loaded(box_size, box_size, particle_count),
				out_extracted(box_size, box_size, particle_count),
				out_centres = micrograph;

		for (long int p = 0; p < particle_count; p++)
		{
			BufferedImage<float> particle_image_RS;

			std::string img_fn = particles_in_micrograph.getString(EMDL_IMAGE_NAME, p);
			particle_image_RS.read(img_fn);

			out_loaded.copySliceFrom(p, particle_image_RS);

			const int extraction_box_size = std::round(extraction_scale * box_size);

			const d2Vector global_position_0(
				particles_in_micrograph.getDouble(EMDL_IMAGE_COORD_X, p),
				particles_in_micrograph.getDouble(EMDL_IMAGE_COORD_Y, p));

			const d2Vector global_position = global_position_0;

			i2Vector integral_position(std::round(global_position.x), std::round(global_position.y));

			BufferedImage<float> extraction_buffer(extraction_box_size, extraction_box_size);

			const int x0 = integral_position.x - extraction_box_size / 2;
			const int y0 = integral_position.y - extraction_box_size / 2;

			for (int y = 0; y < extraction_box_size; y++)
			for (int x = 0; x < extraction_box_size; x++)
			{
				int xx = x0 + x;
				int yy = y0 + y;

				if (xx < 0) xx = 0;
				else if (xx >= micrograph_size.x) xx = micrograph_size.x - 1;

				if (yy < 0) yy = 0;
				else if (yy >= micrograph_size.y) yy = micrograph_size.y - 1;

				extraction_buffer(x,y) = -micrograph(xx, yy);

				out_centres(xx,yy) = mean_value + 6 * std_dev;
			}

			const double mean_value_p = Normalization::computeMean(extraction_buffer);
			extraction_buffer -= mean_value_p;

			extraction_buffer /= std_dev;

			if (std::abs(micrograph_pixel_size - pixel_size) > 0.001)
			{
				particle_image_RS = Resampling::FourierCrop_fullStack(
							extraction_buffer, extraction_scale, num_threads, true);
			}
			else
			{
				particle_image_RS = extraction_buffer;
			}

			out_extracted.copySliceFrom(p, particle_image_RS);
		}

		out_extracted.write(outDir+"micrograph_"+ZIO::itoa(micrograph_id)+"_extracted.mrc");
		out_loaded.write(outDir+"micrograph_"+ZIO::itoa(micrograph_id)+"_loaded.mrc");
		out_centres.write(outDir+"micrograph_"+ZIO::itoa(micrograph_id)+"_centres.mrc");
	}

	return 0;
}

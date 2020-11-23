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
	std::string particlesFn, class_averages_filename, classes_filename, outDir;
	int num_threads, num_frequencies, min_MG, max_MG, max_iterations;
	double radius, edge_padding;
	bool diag, flip_y;


	IOParser parser;

	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input STAR file output by align_2d_classes");
		class_averages_filename = parser.getOption("--ca", "Class averages stack");
		classes_filename = parser.getOption("--classes", "File with a list of 2D classes to consider");
		radius = textToDouble(parser.getOption("--r", "Average blob radius (in bin-1 pixels)", "600"));
		flip_y = parser.checkOption("--flip_y", "Centre of blob lies in positive y direction");
		edge_padding = textToDouble(parser.getOption("--pad", "Edge padding (in bin-1 pixels)", "20"));
		num_frequencies = textToInteger(parser.getOption("--f", "Number of 1D frequencies used to model the blob", "10"));

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


	const int box_size = class_averages.xdim;
	const int num_classes = class_averages.zdim;

	std::vector<int> classes_to_consider = ZIO::readInts(classes_filename);
	const int relevant_class_count = classes_to_consider.size();

	std::vector<int> class_to_subset(num_classes, -1);

	for (int i = 0; i < relevant_class_count; i++)
	{
		class_to_subset[classes_to_consider[i]] = i;
	}


	const double pixel_size = obs_model.getPixelSize(0);


	BufferedImage<float> blob_fits(box_size, box_size, relevant_class_count);
	std::vector<std::vector<double>> all_optimal_parameters(relevant_class_count);



	Log::beginProgress("Fitting blobs", relevant_class_count/num_threads);

	#pragma omp parallel for num_threads(num_threads)
	for (int cc = 0; cc < relevant_class_count; cc++)
	{
		const int thread_id = omp_get_thread_num();

		if (thread_id == 0)
		{
			Log::updateProgress(cc);
		}

		const int class_id = classes_to_consider[cc];

		RawImage<float> slice = class_averages.getSliceRef(class_id);
		const d2Vector initial_centre(box_size/2, box_size/2 + (flip_y? radius : -radius));

		BlobFit2D blob_fit(
			slice, initial_centre, radius/2, 0, 0, num_threads);

		for (int y = 0; y < blob_fit.weight.ydim; y++)
		for (int x = 0; x < blob_fit.weight.xdim; x++)
		{
			const double dx = x - box_size/2;
			const double dy = y - box_size/2;
			const double r = sqrt(dx*dx + dy*dy);

			if (r < box_size/2 - edge_padding)
			{
				blob_fit.weight(x,y) = 1;
			}
			else
			{
				blob_fit.weight(x,y) = 0;
			}
		}


		std::vector<double> initial_parameters(num_frequencies + 2, 0.0);
		initial_parameters[0] = initial_centre.x;
		initial_parameters[1] = initial_centre.y;

		std::vector<double> optimal_parameters = NelderMead::optimize(
				initial_parameters, blob_fit, 2, 0.001, max_iterations, 1.0, 2.0, 0.5, 0.5, false);

		all_optimal_parameters[cc] = optimal_parameters;


		Blob2D blob(optimal_parameters, radius + 2*box_size);

		std::vector<double> radial_average = blob.radialAverage(slice, blob_fit.weight);
		BufferedImage<float> projection = blob.radialAverageProjection(slice, radial_average);

		blob_fits.copySliceFrom(cc, projection);

		// DEBUG:
		/*{
			const double c = box_size/2;
			std::vector<double> coeffs0 = {c, 0.0, 20.0, 0.0, 20.0, 0.0};
			Blob2D b0(coeffs0, box_size);

			std::vector<double> radial_average(box_size+1);

			for (int i = 0; i < box_size; i++)
			{
				radial_average[i] = sin(i);
			}

			BufferedImage<float> projection;

			projection = b0.radialAverageProjection(slice, radial_average);

			projection.write("DEBUG_projection_0.mrc", pixel_size);

			std::vector<double> coeffs1 = Blob2D::rotate(coeffs0, DEG2RAD(30), d2Vector(c,c));
			Blob2D b1(coeffs1, box_size);

			projection = b1.radialAverageProjection(slice, radial_average);

			projection.write("DEBUG_projection_1.mrc", pixel_size);

			std::exit(0);
		}*/
	}

	Log::endProgress();

	blob_fits.write(outDir+"blob_fits_by_class.mrc", pixel_size);



	ZIO::makeOutputDir(outDir+"Frames");

	MetaDataTable output_particles;

	std::vector<MetaDataTable> particles_by_micrograph = StackHelper::splitByMicrographName(particles_table);
	const int micrograph_count = particles_by_micrograph.size();

	Log::beginProgress("Deleting Blobs", micrograph_count);

	for (int m = 0; m < micrograph_count; m++)
	{
		Log::updateProgress(m);

		const MetaDataTable& particles = particles_by_micrograph[m];
		const int particle_count = particles.numberOfObjects();


		std::vector<int> compressed_particle_indices(particle_count, -1);

		int current_index = 0;

		for (int p = 0; p < particle_count; p++)
		{
			const int global_class_id = particles.getIntMinusOne(EMDL_PARTICLE_CLASS, p);
			const int selected_class_id = class_to_subset[global_class_id];

			if (selected_class_id >= 0)
			{
				compressed_particle_indices[p] = current_index;
				current_index++;
			}
		}

		const int relevant_particle_count = current_index;

		if (relevant_particle_count > 0)
		{
			BufferedImage<float> all_particle_images(box_size, box_size, relevant_particle_count);
			BufferedImage<float> all_blob_images(box_size, box_size, relevant_particle_count);

			#pragma omp parallel for num_threads(num_threads)
			for (int p = 0; p < particle_count; p++)
			{
				const int global_class_id = particles.getIntMinusOne(EMDL_PARTICLE_CLASS, p);
				const int selected_class_id = class_to_subset[global_class_id];

				if (selected_class_id >= 0)
				{
					const double dx_A = particles.getDouble(EMDL_ORIENT_ORIGIN_X_ANGSTROM, p);
					const double dy_A = particles.getDouble(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, p);

					std::vector<double> raw_blob_parameters = all_optimal_parameters[selected_class_id];
					raw_blob_parameters[0] += dx_A / pixel_size;
					raw_blob_parameters[1] += dy_A / pixel_size;

					const d2Vector raw_blob_centre(raw_blob_parameters[0], raw_blob_parameters[1]);
					const d2Vector image_centre((double)(box_size/2), (double)(box_size/2));
					const double max_dist = (image_centre - raw_blob_centre).length() + box_size / sqrt(2.0) + 1;

					const double phi = DEG2RAD(particles.getDouble(EMDL_ORIENT_PSI, p));

					const double centre = box_size / 2;
					const d2Vector axis(centre,centre);

					std::vector<double> rotated_blob_parameters = Blob2D::rotate(raw_blob_parameters, phi, axis);
					Blob2D blob(rotated_blob_parameters, max_dist);

					BufferedImage<float> weight(box_size,box_size);
					weight.fill(1.f);


					std::string img_fn = particles.getString(EMDL_IMAGE_NAME, p);
					BufferedImage<float> particle_image;
					particle_image.read(img_fn);

					std::vector<double> rad_average = blob.radialAverage(particle_image, weight);
					BufferedImage<float> projection = blob.radialAverageProjection(particle_image, rad_average);

					particle_image -= projection;

					const int index = compressed_particle_indices[p];
					all_particle_images.copySliceFrom(index, particle_image);
					all_blob_images.copySliceFrom(index, projection);
				}
			}

			std::string img_fn = particles.getString(EMDL_IMAGE_NAME, 0);
			img_fn = outDir + "Frames/" + img_fn.substr(img_fn.find_last_of('/')+1);

			all_particle_images.write(img_fn);
			all_blob_images.write(img_fn.substr(0,img_fn.find_last_of('.'))+"_blobs.mrcs");

			for (int p = 0; p < particle_count; p++)
			{
				const int index = compressed_particle_indices[p];

				if (index >= 0)
				{
					output_particles.addObject();
					output_particles.setObject(particles.getObject(p));
					output_particles.setValue(EMDL_IMAGE_NAME, ZIO::itoa(index+1) + "@" + img_fn);
				}
			}
		}
	}

	Log::endProgress();

	ObservationModel::saveNew(output_particles, obs_model.opticsMdt, outDir + "particles.star");

	return 0;
}

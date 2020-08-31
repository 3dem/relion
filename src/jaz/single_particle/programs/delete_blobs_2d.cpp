#include "delete_blobs_2d.h"

#include <src/jaz/membrane/blob_fit_2d.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/util/drawing.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/local_extrema.h>


using namespace gravis;


void DeleteBlobs2DProgram::readParameters(int argc, char *argv[])
{
	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		micrographs_list_filename = parser.getOption("--i", "Micrograph lists filename");
		micrographs_dir = parser.getOption("--md", "Micrographs directory");
		blobs_dir = parser.getOption("--bd", "Initial blobs directory");

		prior_sigma_A = textToDouble(parser.getOption("--sig", "Weight of initial position", "0"));
		max_binning = textToDouble(parser.getOption("--bin0", "Initial (maximal) binning factor", "8"));
		min_binning = textToDouble(parser.getOption("--bin1", "Final (minimal) binning factor", "4"));

		diag = parser.checkOption("--diag", "Write out diagnostic information");

		max_frequencies = textToInteger(parser.getOption("--n", "Number of frequencies", "12"));
		blob_thickness = textToDouble(parser.getOption("--th", "Blob thickness [fraction of radius]", "0.5"));
		highpass_sigma_real_A = textToDouble(parser.getOption("--hp", "High-pass sigma [Å, real space]", "300"));
		max_iters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "1000"));
		convergence_threshold = textToDouble(parser.getOption("--cth", "Convergence threshold", "0.02"));
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));

		outPath = parser.getOption("--o", "Output filename pattern");

		Log::readParams(parser);

		if (parser.checkForErrors())
		{
			parser.writeUsage(std::cout);
			exit(1);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}

	outPath = ZIO::makeOutputDir(outPath);
	
	ZIO::makeOutputDir(outPath + "Frames");
	ZIO::makeOutputDir(outPath + "Blobs");

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag");
	}
	
	micrographs_dir = ZIO::ensureEndingSlash(micrographs_dir);
	blobs_dir = ZIO::ensureEndingSlash(blobs_dir);
}

void DeleteBlobs2DProgram::run()
{
	std::vector<std::string> all_micrograph_names;
	all_micrograph_names.reserve(1000);
	
	std::ifstream file(micrographs_list_filename);
	
	if (!file)
	{
		REPORT_ERROR("DelineatedBlob2D::run: unable to read " + micrographs_list_filename);
	}
	
	std::string line;
	
	while (std::getline(file, line))
	{
		std::stringstream sts;
		sts << line;
		
		std::string name;
		sts >> name;
		
		all_micrograph_names.push_back(name);
	}
	
	const int micrograph_count = all_micrograph_names.size();
	
	t3Vector<long int> micrograph_size = ImageFileHelper::getSize(
	            micrographs_dir + all_micrograph_names[0] + ".mrc");
	
	const double visualisation_binning = 32.0;
	
	i2Vector visualisation_size = Resampling::getFourierCroppedSize2D(
	            micrograph_size.x, micrograph_size.y, visualisation_binning, true);
	
	BufferedImage<float> visualisation(
	            visualisation_size.x, 
	            visualisation_size.y, 
	            3 * micrograph_count);
	
	std::vector<std::string> failed_micrographs;
	
	
	#pragma omp parallel for num_threads(num_threads)
	for (int m = 0; m < micrograph_count; m++)
	{
		const int thread_id = omp_get_thread_num();
		const bool verbose = thread_id == 0;
		
		if (verbose)
		{
			Log::beginSection(
				"Micrograph "+ZIO::itoa(m+1)+"/"
			    +ZIO::itoa((int)std::ceil(micrograph_count / (double)num_threads))+" (out of "
			    +ZIO::itoa(micrograph_count)+")");
		}
		
		try 
		{
			processMicrograph(
				m,
				micrographs_dir + all_micrograph_names[m] + ".mrc",
				blobs_dir + all_micrograph_names[m] + ".blobs",
				visualisation,
				visualisation_binning,
				verbose);
		}
		catch (Exception e)
		{
			#pragma omp critical
			{
				failed_micrographs.push_back(all_micrograph_names[m]);
			}
		}
		
		if (verbose)
		{
			Log::endSection();
		}
	}
	
	visualisation.write(outPath + "diagnostic.mrc");
	
	if (failed_micrographs.size() > 0)
	{
		ZIO::writeToFile(failed_micrographs, outPath + "failed_micrographs.txt");
	}
}

void DeleteBlobs2DProgram::processMicrograph(
        int micrograph_index,
        const std::string& micrograph_filename,
		const std::string& blobs_filename,
        RawImage<float>& visualisation,
        double visualisation_binning,
        bool verbose)
{	        
	if (verbose)
	{
		Log::print("Loading");
	}

	BufferedImage<float> micrograph;
	micrograph.read(micrograph_filename);

	const int w_full = micrograph.xdim;
	const int h_full = micrograph.ydim;

	const double pixel_size = ImageFileHelper::getSamplingRate(micrograph_filename);

	const double highpass_sigma_real = highpass_sigma_real_A / pixel_size;
	
	std::vector<DelineatedBlob2D> delineated_blobs = DelineatedBlob2D::read(blobs_filename);
	
	
	const int blob_count = delineated_blobs.size();

	if (verbose)
	{
		Log::print("Filtering");
	}

	BufferedImage<float> micrograph_filtered = ImageFilter::highpassStackGaussPadded(
				micrograph, highpass_sigma_real, 1);
	
	BufferedImage<float> blobs_image(w_full, h_full, 1);
	blobs_image.fill(0.f);
	
	BufferedImage<float> erased_image = micrograph;
	
	BufferedImage<float> dummy_weight(w_full, h_full, 1);
	dummy_weight.fill(1.f);
	
	std::string micrograph_name = micrograph_filename.substr(
	            micrograph_filename.find_last_of('/')+1);
	
	micrograph_name = micrograph_name.substr(
	            0, micrograph_filename.find_last_of('.'));


	for (int blob_id = 0; blob_id < blob_count; blob_id++)
	{
		if (verbose)
		{
			Log::beginSection("Blob " + ZIO::itoa(blob_id + 1)+"/"+ZIO::itoa(blob_count));
		}
		
		const double radius = delineated_blobs[blob_id].radius;
		Blob2D initial_blob = delineated_blobs[blob_id].blob;
		
		std::vector<double> initial_parameters_fullsize = initial_blob.toVector();

		std::vector<double> blob_coeffs = initial_parameters_fullsize;

		i2Vector window_origin_full;
		BufferedImage<float> blob_region_full;

		{
			const double blob_padding_full = radius / 2.0;
			const double smoothing_radius_full = 0.5 * radius;

			// determine bounding box at bin 1

			Blob2D blob_full(initial_parameters_fullsize, smoothing_radius_full);

			std::pair<d2Vector, d2Vector> bounding_box_full = blob_full.scanForBoundingBox(
						radius, blob_padding_full, (int)(2 * PI * radius));

			window_origin_full = i2Vector(
					std::floor(bounding_box_full.first.x),
					std::floor(bounding_box_full.first.y));

			i2Vector window_size_full(
					std::ceil(bounding_box_full.second.x) - window_origin_full.x,
					std::ceil(bounding_box_full.second.y) - window_origin_full.y);

			i2Vector window_size_binned_max(
						std::round(window_size_full.x / max_binning),
						std::round(window_size_full.y / max_binning));

			if (window_size_binned_max.x < 1 || window_size_binned_max.y < 1)
			{
				Log::endSection();
				continue;
			}


			// ensure an even box size at the greatest binning level

			if (window_size_binned_max.x % 2 == 1) window_size_binned_max.x += 1;
			if (window_size_binned_max.y % 2 == 1) window_size_binned_max.y += 1;


			// ensure an integral ratio between window sizes

			window_size_full.x = max_binning * window_size_binned_max.x;
			window_size_full.y = max_binning * window_size_binned_max.y;

			if (window_size_full.x > micrograph.xdim) window_size_full.x = micrograph.xdim;
			if (window_size_full.y > micrograph.ydim) window_size_full.y = micrograph.ydim;


			// shift full-size window to lie inside the image

			if (window_origin_full.x + window_size_full.x > micrograph.xdim)
			{
				window_origin_full.x = micrograph.xdim - window_size_full.x;
			}
			else if (window_origin_full.x < 0)
			{
				window_origin_full.x = 0;
			}

			if (window_origin_full.y + window_size_full.y > micrograph.ydim)
			{
				window_origin_full.y = micrograph.ydim - window_size_full.y;
			}
			else if (window_origin_full.y < 0)
			{
				window_origin_full.y = 0;
			}

			blob_region_full.resize(window_size_full.x, window_size_full.y);

			// extract full-size image

			for (int y = 0; y < window_size_full.y; y++)
			for (int x = 0; x < window_size_full.x; x++)
			{
				const int x0 = window_origin_full.x;
				const int y0 = window_origin_full.y;

				blob_region_full(x,y) = micrograph_filtered(x0 + x, y0 + y);
			}
		}


		std::vector<double> initial_parameters_cropped = initial_parameters_fullsize;
		initial_parameters_cropped[0] -= window_origin_full.x;
		initial_parameters_cropped[1] -= window_origin_full.y;

		std::vector<double> blob_parameters_cropped = initial_parameters_cropped;



		double current_binning = max_binning;

		while (current_binning > min_binning - 1e-6)
		{
			if (verbose)
			{
				if (current_binning == max_binning)
				{
					Log::beginSection("Fitting at bin " + ZIO::itoa((int)current_binning));
				}
				else
				{
					Log::beginSection("Refining at bin " + ZIO::itoa((int)current_binning));
				}
			}
			
			blob_parameters_cropped = fitBlob(
				blob_id, blob_parameters_cropped,
				radius, pixel_size, current_binning,
				blob_region_full,
				micrograph_name,
				verbose);

			current_binning /= 2;
			
			if (verbose)
			{
				Log::endSection();
			}
		}
		
		if (verbose)
		{
			Log::print("Erasing");
		}

		std::vector<double> final_parameters_fullsize = blob_parameters_cropped;
		final_parameters_fullsize[0] += window_origin_full.x;
		final_parameters_fullsize[1] += window_origin_full.y;

		Blob2D final_blob(blob_coeffs, radius/2);
		final_blob.erase(micrograph, erased_image, blobs_image, dummy_weight, 1.5 * radius, radius);
		
		if (verbose)
		{
			Log::endSection();
		}
	}	  
	
	if (verbose)
	{
		Log::print("Writing output");
	}

	erased_image.write(outPath + "Frames/" + micrograph_name + ".mrc", pixel_size);
	blobs_image.write(outPath + "Blobs/" + micrograph_name + ".mrc", pixel_size);
	
	BufferedImage<float> cropped_original = Resampling::FourierCrop_fullStack(
	            micrograph, visualisation_binning, 1, true);
	
	BufferedImage<float> cropped_erased = Resampling::FourierCrop_fullStack(
	            erased_image, visualisation_binning, 1, true);
	
	BufferedImage<float> cropped_blobs = Resampling::FourierCrop_fullStack(
	            blobs_image, visualisation_binning, 1, true);
	
	const float mu = Normalization::computeMean(cropped_original);
	cropped_original -= mu;
	cropped_erased -= mu;
	
	visualisation.copySliceFrom(3 * micrograph_index,     cropped_original);
	visualisation.copySliceFrom(3 * micrograph_index + 1, cropped_erased);
	visualisation.copySliceFrom(3 * micrograph_index + 2, cropped_blobs);
}

std::vector<double> DeleteBlobs2DProgram::fitBlob(
		int blob_id,
		const std::vector<double>& initial_blob_params_cropped,
		double radius_full,
		double pixel_size_full,
		double binning_factor,
		BufferedImage<float>& blob_region_full,
		const std::string& image_name,
		bool verbose)
{
	const double radius_binned = radius_full / binning_factor;
	const double prior_sigma = prior_sigma_A / pixel_size_full;
	const double initial_step = 2.0;
	const double smoothing_radius_binned = 0.5 * radius_full / binning_factor;

	std::string blobTag = "blob_" + ZIO::itoa(blob_id);
	std::string binTag = "bin_" + ZIO::itoa((int)binning_factor);

	std::string outTag = outPath + "diag/" + image_name + "/" + blobTag + "/" + binTag;
	
		

	const d2Vector initial_position_cropped = d2Vector(initial_blob_params_cropped[0], initial_blob_params_cropped[1]);

	BufferedImage<float> blob_region_binned = Resampling::FourierCrop_fullStack(
				blob_region_full, binning_factor, 1, true);
	
	//blob_region_binned = ImageFilter::Gauss2D(blob_region_binned, 0, 3, true);

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag/" + image_name + "/" + blobTag);
		blob_region_binned.write(outTag + "_data.mrc");
	}


	std::vector<double> last_optimum = fromBin1(initial_blob_params_cropped, binning_factor);
	int initial_frequencies = initial_blob_params_cropped.size() < 2? 0 : (initial_blob_params_cropped.size() - 2) / 2;


	for (int current_frequencies = initial_frequencies;
		 current_frequencies <= max_frequencies; current_frequencies++)
	{
		if (verbose)
		{
			Log::print(ZIO::itoa(current_frequencies) + " frequencies ("
				   + ZIO::itoa(2*current_frequencies + 2) + " parameters)");
		}

		std::string tag = outTag + "_N_" + ZIO::itoa(current_frequencies);

		
		const d2Vector initial_position_cropped_binned = initial_position_cropped / binning_factor;

		BlobFit2D blob_fit(
			blob_region_binned, initial_position_cropped_binned, current_frequencies,
			smoothing_radius_binned, prior_sigma, 1);
		
		Blob2D blob0(last_optimum, smoothing_radius_binned);
		
		{
			const double r0 = radius_binned - blob_thickness * radius_binned / 2;
			const double r1 = radius_binned + blob_thickness * radius_binned / 2;
			
			blob_fit.computeWeight(blob0, r0, r1);
		}
		
		std::vector<double> current_optimum(2*current_frequencies + 2, 0.0);

		for (int i = 0; i < last_optimum.size(); i++)
		{
			current_optimum[i] = last_optimum[i];
		}

		std::vector<double> new_optimum;

		new_optimum = NelderMead::optimize(
						current_optimum, blob_fit, 
						initial_step, convergence_threshold, max_iters, 
						1.0, 2.0, 0.5, 0.5, false);

		if (diag)
		{
			Blob2D blob2(new_optimum, smoothing_radius_binned);
			
			const double E0 = blob_fit.f(current_optimum, 0);
			const double E1 = blob_fit.f(new_optimum, 0);

			if (verbose)
			{
				Log::print("  E: "+ZIO::itoa(E0)+" -> "+ZIO::itoa(E1));
			}

			drawTestStack(blob2, blob_region_binned, blob_fit.weight).write(tag+"_fit.mrc");
			drawFit(blob2, blob_region_binned, blob_fit.weight).write(tag+"_residual.mrc");
			blob_fit.weight.write(tag+"_weight.mrc");
		}

		last_optimum = new_optimum;
	}

	return toBin1(last_optimum, binning_factor);
}


BufferedImage<float> DeleteBlobs2DProgram::drawFit(
		Blob2D& blob,
		const BufferedImage<float>& image,
		const BufferedImage<float>& realWeight)
{
	std::vector<double> radAvg = blob.radialAverage(image, realWeight);
	BufferedImage<float> diff = blob.drawError(image, realWeight, radAvg);

	return diff;
}


BufferedImage<float> DeleteBlobs2DProgram::drawTestStack(
		Blob2D& blob,
		const BufferedImage<float>& image,
		const BufferedImage<float>& realWeight)
{
	std::vector<double> radAvg = blob.radialAverage(image, realWeight);
	BufferedImage<float> proj = blob.radialAverageProjection(image, radAvg);

	return proj;
}

std::vector<double> DeleteBlobs2DProgram::toBin1(
		const std::vector<double> &parameters,
		double binning_factor)
{
	std::vector<double> upscaled_optimum = parameters;

	for (int i = 0; i < upscaled_optimum.size(); i++)
	{
		upscaled_optimum[i] *= binning_factor;
	}

	return upscaled_optimum;
}

std::vector<double> DeleteBlobs2DProgram::fromBin1(
		const std::vector<double> &parameters,
		double binning_factor)
{
	const int xc = parameters.size();
	std::vector<double> downscaled_optimum(xc);

	for (int i = 0; i < downscaled_optimum.size(); i++)
	{
		downscaled_optimum[i] = parameters[i] / binning_factor;
	}

	return downscaled_optimum;
}


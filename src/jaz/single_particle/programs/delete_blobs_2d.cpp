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

		blobs_filename = parser.getOption("--i", "Initial blob locations");
		micrograph_filename = parser.getOption("--m", "Micrograph file name");

		prior_sigma_A = textToDouble(parser.getOption("--sig", "Weight of initial position", "0"));
		max_binning = textToDouble(parser.getOption("--bin0", "Initial (maximal) binning factor", "8"));
		min_binning = textToDouble(parser.getOption("--bin1", "Final (minimal) binning factor", "2"));

		diag = parser.checkOption("--diag", "Write out diagnostic information");

		max_frequencies = textToInteger(parser.getOption("--n", "Number of frequencies", "50"));
		highpass_sigma_real_A = textToDouble(parser.getOption("--hp", "High-pass sigma [Å, real space]", "300"));
		max_iters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "1000"));
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

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag");
	}
}

void DeleteBlobs2DProgram::run()
{
	Log::print("Loading");

	BufferedImage<float> micrograph;
	micrograph.read(micrograph_filename);

	const int w_full = micrograph.xdim;
	const int h_full = micrograph.ydim;

	const double pixel_size = ImageFileHelper::getSamplingRate(micrograph_filename);

	const double highpass_sigma_real = highpass_sigma_real_A / pixel_size;
	
	std::vector<DelineatedBlob2D> delineated_blobs = DelineatedBlob2D::read(blobs_filename);
	
	
	
	const int blob_count = delineated_blobs.size();

	Log::print("Filtering");


	micrograph.write(outPath+"DEBUG_micrograph.mrc");

	BufferedImage<float> micrograph_filtered = ImageFilter::highpassStackGaussPadded(
				micrograph, highpass_sigma_real, num_threads);

	micrograph_filtered.write(outPath+"DEBUG_micrograph_filtered.mrc");
	
	BufferedImage<float> blobs_image(w_full, h_full, 1);
	blobs_image.fill(0.f);
	
	
	std::string micrograph_name = micrograph_filename.substr(
	            micrograph_filename.find_last_of('/')+1);
	
	micrograph_name = micrograph_name.substr(
	            0, micrograph_filename.find_last_of('.'));


	for (int blob_id = 0; blob_id < blob_count; blob_id++)
	{
		Log::beginSection("Blob #" + ZIO::itoa(blob_id + 1));
		
		Blob2D blob = delineated_blobs[blob_id].blob;
		const double radius = delineated_blobs[blob_id].radius;
		
		std::vector<double> initial = blob.toVector();

		std::vector<double> blob_coeffs = initial;


		double current_binning = max_binning;

		while (current_binning > min_binning - 1e-6)
		{
			if (current_binning == max_binning)
			{
				Log::beginSection("Fitting at bin " + ZIO::itoa((int)current_binning));
			}
			else
			{
				Log::beginSection("Refining at bin " + ZIO::itoa((int)current_binning));
			}
			
			blob_coeffs = fitBlob(
						blob_id, blob_coeffs, 
			            radius, pixel_size, current_binning, 
			            micrograph_filtered,
						micrograph_name);

			current_binning /= 2;
			
			Log::endSection();
		}

		//Blob2D blob(blob_coeffs, outer_radius_full);


			/*RawImage<float> stackSlice = original_stack.getSliceRef(f);
			RawImage<float> blobsSlice = blobs_stack.getSliceRef(f);

			blob.decompose(
				stackSlice, blobsSlice,
				tomogram0.projectionMatrices[f],
				fiducialsMask.getConstSliceRef(f),
				50.0);*/

		
		Log::endSection();
	}

	/*const std::string tag = outPath + tomogram0.name;

	original_stack.write(tag + "_subtracted.mrc", tomogram0.optics.pixelSize);
	blobs_stack.write(tag + "_blobs.mrc", tomogram0.optics.pixelSize);

	subtracted_tomogram_set.setTiltSeriesFile(tomo_index, tag + "_subtracted.mrc");
	blobs_tomogram_set.setTiltSeriesFile(tomo_index, tag + "_blobs.mrc");

	if (diag)
	{
		visualisation.getSliceRef(3 * tomo_batch_index + 1).copyFrom(
				original_stack.getConstSliceRef(frame_count / 2));

		visualisation.getSliceRef(3 * tomo_batch_index + 2).copyFrom(
				blobs_stack.getConstSliceRef(frame_count / 2));
	}*/
}

std::vector<double> DeleteBlobs2DProgram::fitBlob(
		int blob_id,
		const std::vector<double>& initial_parameters,
        double radius_full,
		double pixel_size_full,
		double binning_factor,
		const RawImage<float>& image_full,
		const std::string& image_name)
{
	const double radius_binned = radius_full / binning_factor;
	const double prior_sigma = prior_sigma_A / pixel_size_full;
	const double initial_step = 2.0;
	const double blob_padding_full = radius_full / 2.0;
	const double smoothing_radius_full = 0.5 * radius_full;
	const double smoothing_radius_binned = 0.5 * radius_full / binning_factor;

	std::string blobTag = "blob_" + ZIO::itoa(blob_id);
	std::string binTag = "bin_" + ZIO::itoa((int)binning_factor);

	std::string outTag = outPath + "diag/" + image_name + "/" + blobTag + "/" + binTag;
	
	
	// determine bounding box at bin 1
	
	Blob2D blob_full(initial_parameters, smoothing_radius_full);
	
	std::pair<d2Vector, d2Vector> bounding_box_full = blob_full.scanForBoundingBox(
	            radius_full, blob_padding_full, (int)(4 * PI * radius_binned));
	
	i2Vector window_origin_full(
			std::floor(bounding_box_full.first.x),
			std::floor(bounding_box_full.first.y));
	          
	i2Vector window_size_full(
			std::ceil(bounding_box_full.second.x) - window_origin_full.x,
			std::ceil(bounding_box_full.second.y) - window_origin_full.y);	
	
	i2Vector window_size_binned_max(
	            std::round(window_size_full.x / max_binning),
	            std::round(window_size_full.y / max_binning));
	
	
	// ensure an even box size at greatest binning level
	
	if (window_size_binned_max.x % 2 == 1) window_size_binned_max.x += 1;
	if (window_size_binned_max.y % 2 == 1) window_size_binned_max.y += 1;
	
	
	// ensure an integral ratio between window sizes
	
	window_size_full.x = max_binning * window_size_binned_max.x;
	window_size_full.y = max_binning * window_size_binned_max.y;
	
	if (window_size_full.x > image_full.xdim) window_size_full.x = image_full.xdim;
	if (window_size_full.y > image_full.ydim) window_size_full.y = image_full.ydim;
	
	       
	// shift full-size window to lie inside the image
	
	if (window_origin_full.x + window_size_full.x > image_full.xdim)
	{
		window_origin_full.x = image_full.xdim - window_size_full.x;
	}
	else if (window_origin_full.x < 0)
	{
		window_origin_full.x = 0;
	}
	
	if (window_origin_full.y + window_size_full.y > image_full.ydim)
	{
		window_origin_full.y = image_full.ydim - window_size_full.y;
	}
	else if (window_origin_full.y < 0)
	{
		window_origin_full.y = 0;
	}
	
	// extract full-size image
	
	BufferedImage<float> blob_region_full(window_size_full.x, window_size_full.y);

	for (int y = 0; y < window_size_full.y; y++)
	for (int x = 0; x < window_size_full.x; x++)
	{
		const int x0 = window_origin_full.x;
		const int y0 = window_origin_full.y;

		blob_region_full(x,y) = image_full(x0 + x, y0 + y);
	}

	std::vector<double> initial_cropped = initial_parameters;
	initial_cropped[0] -= window_origin_full.x;
	initial_cropped[1] -= window_origin_full.y;
	

	const d2Vector initial_position_cropped = d2Vector(initial_cropped[0], initial_cropped[1]);

	BufferedImage<float> blob_region_binned = Resampling::FourierCrop_fullStack(
				blob_region_full, binning_factor, num_threads, true);
	
	blob_region_full.write(outPath + "blob_region_full.mrc");
	blob_region_binned.write(outPath + "blob_region_binned.mrc");
	
	//blob_region_binned = ImageFilter::Gauss2D(blob_region_binned, 0, 3, true);

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag/" + image_name + "/" + blobTag);
		blob_region_binned.write(outTag + "_data.mrc");
	}


	std::vector<double> last_optimum = fromBin1(initial_cropped, binning_factor);
	int initial_frequencies = initial_cropped.size() < 2? 0 : (initial_cropped.size() - 2) / 2;


	for (int current_frequencies = initial_frequencies;
		 current_frequencies <= max_frequencies; current_frequencies++)
	{
		Log::print(ZIO::itoa(current_frequencies) + " frequencies ("
				   + ZIO::itoa(2*current_frequencies + 2) + " parameters)");

		std::string tag = outTag + "_N_" + ZIO::itoa(current_frequencies);

		
		const d2Vector initial_position_cropped_binned = initial_position_cropped / binning_factor;

		BlobFit2D blob_fit(
			blob_region_binned, initial_position_cropped_binned, current_frequencies,
			smoothing_radius_binned, prior_sigma, num_threads);
		
		Blob2D blob0(last_optimum, smoothing_radius_binned);
		
		blob_fit.computeWeight(blob0, 0.5 * radius_binned, 1.5 * radius_binned);
		
		std::vector<double> current_optimum(2*current_frequencies + 2, 0.0);

		for (int i = 0; i < last_optimum.size(); i++)
		{
			current_optimum[i] = last_optimum[i];
		}

		std::vector<double> new_optimum;

		new_optimum = NelderMead::optimize(
				current_optimum, blob_fit, initial_step, 0.001, max_iters, 1.0, 2.0, 0.5, 0.5, false);

		if (diag)
		{
			Blob2D blob2(new_optimum, smoothing_radius_binned);
			
			const double E0 = blob_fit.f(current_optimum, 0);
			const double E1 = blob_fit.f(new_optimum, 0);

			Log::print("  E: "+ZIO::itoa(E0)+" -> "+ZIO::itoa(E1));

			drawTestStack(blob2, blob_region_binned, blob_fit.weight).write(tag+"_fit.mrc");
			drawFit(blob2, blob_region_binned, blob_fit.weight).write(tag+"_residual.mrc");
			blob_fit.weight.write(tag+"_weight.mrc");
		}

		last_optimum = new_optimum;
	}


	std::vector<double> upscaled_optimum = toBin1(last_optimum, binning_factor);

	upscaled_optimum[0] += window_origin_full.x;
	upscaled_optimum[1] += window_origin_full.y;

	return upscaled_optimum;
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


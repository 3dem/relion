#include "delete_blobs_2d.h"

#include <src/jaz/membrane/blob_fit_2d.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/image/filter.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/image/normalization.h>


using namespace gravis;


void DeleteBlobs2DProgram::readParameters(int argc, char *argv[])
{
	double sphere_thickness_0;

	IOParser parser;

	try
	{
		parser.setCommandLine(argc, argv);
		int gen_section = parser.addSection("General options");

		blob_list_filename = parser.getOption("--i", "Initial blob locations", "blobs.star");
		micrographs_list_filename = parser.getOption("--m", "Micrographs list", "micrographs.star");

		blob_radius_A = textToDouble(parser.getOption("--r", "Blob radius [Å]"));
		blob_thickness_A = textToDouble(parser.getOption("--th", "Blob membrane thickness [Å]"));

		prior_sigma_A = textToDouble(parser.getOption("--sig", "Uncertainty std. dev. of initial position [Å]", "10"));
		max_binning = textToDouble(parser.getOption("--bin0", "Initial (maximal) binning factor", "8"));
		min_binning = textToDouble(parser.getOption("--bin1", "Final (minimal) binning factor", "2"));

		diag = parser.checkOption("--diag", "Write out diagnostic information");

		min_MG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
		max_MG = textToInteger(parser.getOption("--max_MG", "Last micrograph index", "-1"));

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
	MetaDataTable micrographs_table;
	micrographs_table.read(micrographs_list_filename);
	const int micrograph_count = micrographs_table.numberOfObjects();

	std::ifstream blob_list_file(blob_list_filename);

	if (!blob_list_file)
	{
		REPORT_ERROR("Unable to read: "+blob_list_filename);
	}

	std::vector<MetaDataTable> blob_tables = MetaDataTable::readAll(blob_list_file, micrograph_count);


	BufferedImage<float> visualisation(1,1,1);

	if (max_MG < 0)
	{
		max_MG = micrograph_count - 1;
	}

	for (int m = min_MG; m <= max_MG; m++)
	{
		std::string micrograph_name = micrographs_table.getString(EMDL_MICROGRAPH_NAME, m);

		Log::beginSection("Micrograph " + micrograph_name);

		processMicrograph(
			blob_tables[m],
			micrograph_name,
			visualisation,
			m, micrograph_count);
		
		Log::endSection();
	}

	if (diag)
	{
		visualisation.write(outPath + "diagnostic.mrc");
	}
}

void DeleteBlobs2DProgram::processMicrograph(
		MetaDataTable& detected_blobs,
		const std::string& micrograph_name,
		BufferedImage<float>& visualisation,
		int micrograph_index,
		int micrograph_count)
{
	Log::print("Loading");

	BufferedImage<float> micrograph;
	micrograph.read(micrograph_name);

	const int w_full = micrograph.xdim;
	const int h_full = micrograph.ydim;

	const double pixel_size = ImageFileHelper::getSamplingRate(micrograph_name);

	const double highpass_sigma_real = highpass_sigma_real_A / pixel_size;
	const double outer_radius_A = blob_radius_A + blob_thickness_A / 2;
	const double outer_radius_full = outer_radius_A / pixel_size;

	const int blob_count = detected_blobs.numberOfObjects();

	Log::print("Filtering");


	micrograph.write("DEBUG_micrograph.mrc");

	BufferedImage<float> micrograph_filtered = ImageFilter::highpassStackGaussPadded(
				micrograph, highpass_sigma_real, num_threads);


	micrograph_filtered.write("DEBUG_micrograph_filtered.mrc");


	/*{
		const double bin = 32;

		BufferedImage<float> micrograph_binned = Resampling::FourierCrop_fullStack(
					micrograph_filtered, bin, num_threads, true);

		micrograph_binned.write("DEBUG_micrograph_bin32.mrc");
		BufferedImage<float> symm = evaluateRotationalSymmetry(micrograph_binned, 17, 23, 5);

		float mean = Normalization::computeMean(symm);

		symm -= mean;

		symm.write("DEBUG_symmetry_bin32.mrc");


		BufferedImage<float> box_maxima = LocalExtrema::boxMaxima(detections, radius / binning);



		for (int blob_id = 0; blob_id < blob_count; blob_id++)
		{
			const d2Vector pos = d2Vector(
				detected_blobs.getDouble(EMDL_IMAGE_COORD_X, blob_id),
				detected_blobs.getDouble(EMDL_IMAGE_COORD_Y, blob_id)) / bin;

			int px = (int)(pos.x + 0.5);
			int py = (int)(pos.y + 0.5);

			symm(px,py) = 150;
		}

		symm.write("DEBUG_symmetry+dots_bin32.mrc");


		std::exit(0);
	}*/

	/*if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag/" + micrograph_name);

		if (visualisation.xdim == 0)
		{
			visualisation.resize(w_full, h_full, 3 * micrograph_count);
		}

		visualisation.getSliceRef(3 * tomo_batch_index).copyFrom(
				original_stack.getConstSliceRef(frame_count / 2));
	}*/

	
	BufferedImage<float> blobs_image(w_full, h_full, 1);
	blobs_image.fill(0.f);



	//for (int blob_id = 0; blob_id < blob_count; blob_id++)
	{
		int blob_id = 1;

		Log::beginSection("Blob #" + ZIO::itoa(blob_id + 1));
		
		const d2Vector initial_blob_position(
				detected_blobs.getDouble(EMDL_IMAGE_COORD_X, blob_id),
				detected_blobs.getDouble(EMDL_IMAGE_COORD_Y, blob_id));

		std::vector<double> initial = {initial_blob_position.x, initial_blob_position.y};

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
						blob_id, blob_coeffs, outer_radius_full,
						pixel_size, current_binning, micrograph_filtered,
						micrograph_index);

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
		const std::vector<double>& initial,
		double outer_radius_full,
		double pixel_size_full,
		double binning_factor,
		RawImage<float>& image_full,
		int micrograph_index)
{
	const double pixel_size_binned = pixel_size_full * binning_factor;
	const double outer_radius_binned = outer_radius_full / binning_factor;
	const double prior_sigma = prior_sigma_A / pixel_size_full;
	const double initial_step = 2.0;


	std::string blobTag = "blob_" + ZIO::itoa(blob_id);
	std::string binTag = "bin_" + ZIO::itoa((int)binning_factor);

	std::string outTag = outPath + "diag/micrograph_" + ZIO::itoa(micrograph_index) + "/" + blobTag + "/" + binTag;


	const int window_size_binned = 2 * std::ceil(outer_radius_full / max_binning);
	const int window_size = (int)(window_size_binned * max_binning);

	i2Vector window_origin(
				std::round(initial[0]) - window_size/2,
				std::round(initial[1]) - window_size/2);

	if (window_origin.x + window_size > image_full.xdim)
	{
		window_origin.x = image_full.xdim - window_size;
	}
	else if (window_origin.x < 0)
	{
		window_origin.x = 0;
	}

	if (window_origin.y + window_size > image_full.ydim)
	{
		window_origin.y = image_full.ydim - window_size;
	}
	else if (window_origin.y < 0)
	{
		window_origin.y = 0;
	}

	BufferedImage<float> blob_region_full(window_size, window_size);

	for (int y = 0; y < window_size; y++)
	for (int x = 0; x < window_size; x++)
	{
		const int x0 = window_origin.x;
		const int y0 = window_origin.y;

		blob_region_full(x,y) = image_full(x0 + x, y0 + y);
	}

	std::vector<double> initial_cropped = initial;
	initial_cropped[0] -= window_origin.x;
	initial_cropped[1] -= window_origin.y;

	const d2Vector initial_position_cropped = d2Vector(initial_cropped[0], initial_cropped[1]);

	BufferedImage<float> blob_region_binned = Resampling::FourierCrop_fullStack(
				blob_region_full, binning_factor, num_threads, true);

	blob_region_binned = ImageFilter::Gauss2D(blob_region_binned, 0, 3, true);

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag/micrograph_" + ZIO::itoa(micrograph_index) + "/" + blobTag);
		blob_region_binned.write(outTag + "_data.mrc");
	}


	std::vector<double> last_optimum = fromBin1(initial_cropped, binning_factor);
	int initial_frequencies = initial_cropped.size() < 2? 0 : initial_cropped.size() - 2;


	for (int current_frequencies = initial_frequencies;
		 current_frequencies <= max_frequencies; current_frequencies++)
	{
		Log::print(ZIO::itoa(current_frequencies) + " frequencies ("
				   + ZIO::itoa(current_frequencies + 2) + " parameters)");

		std::string tag = outTag + "_N_" + ZIO::itoa(current_frequencies);

		const d2Vector initial_position_cropped_binned = initial_position_cropped / binning_factor;
		const double radius_binned = blob_radius_A / pixel_size_binned;
		const double thickness_binned = blob_thickness_A / pixel_size_binned;

		BlobFit2D bf(
			blob_region_binned, initial_position_cropped_binned, current_frequencies,
			radius_binned, thickness_binned, prior_sigma, num_threads);

		std::vector<double> current_optimum(current_frequencies + 2, 0.0);

		for (int i = 0; i < last_optimum.size(); i++)
		{
			current_optimum[i] = last_optimum[i];
		}

		std::vector<double> new_optimum;

		new_optimum = NelderMead::optimize(
				current_optimum, bf, initial_step, 0.001, max_iters, 1.0, 2.0, 0.5, 0.5, false);

		if (diag)
		{
			Blob2D blob2(new_optimum, outer_radius_binned);
			
			const double E0 = bf.f(current_optimum, 0);
			const double E1 = bf.f(new_optimum, 0);

			Log::print("  E: "+ZIO::itoa(E0)+" -> "+ZIO::itoa(E1));

			drawTestStack(blob2, blob_region_binned, bf.weight).write(tag+"_fit.mrc");
			drawFit(blob2, blob_region_binned, bf.weight).write(tag+"_residual.mrc");
		}

		last_optimum = new_optimum;
	}


	std::vector<double> upscaled_optimum = toBin1(last_optimum, binning_factor);

	upscaled_optimum[0] += window_origin.x;
	upscaled_optimum[1] += window_origin.y;

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

BufferedImage<float> DeleteBlobs2DProgram::evaluateRotationalSymmetry(
		const RawImage<float> &image, double radius, double max_radius, double sigma)
{
	const int w = image.xdim;
	const int h = image.ydim;

	BufferedImage<float> out(w,h);

	BufferedImage<float> weight(w,h);
	weight.fill(1.f);

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		std::vector<double> params{(double)x, (double)y};
		Blob2D blob(params, max_radius);

		std::vector<double> radAvg = blob.radialAverage(image, weight);

		double sum = 0, sum_wgh = 0;

		for (int i = 0; i < radAvg.size(); i++)
		{
			sum += radAvg[i] * i;
			sum_wgh += i;
		}

		const double avg = sum / sum_wgh;

		double power = 0;

		for (int i = 0; i < radAvg.size(); i++)
		{
			const double d = radAvg[i] - avg;
			const double dr = i - radius;
			const double g = exp(-dr/(2*sigma*sigma));

			power += g * d * d * i;
		}

		out(x,y) = power;
	}

	return out;
}


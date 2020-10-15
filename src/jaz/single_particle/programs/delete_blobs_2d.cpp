#include "delete_blobs_2d.h"

#include <src/jaz/membrane/blob_fit_2d.h>
#include <src/jaz/membrane/global_blob_fit.h>
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
#include <src/jaz/single_particle/stack_helper.h>


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
		particles_file = parser.getOption("--ptc", "Optional particles file for phase flipping", "");
		
		ring_min = textToDouble(parser.getOption("--ring_min", "Inner radius of ring to isolate (relative to blob surface)", "0"));
		ring_max = textToDouble(parser.getOption("--ring_max", "Outer radius of ring to isolate", "0"));
		ring_edge_sigma = textToDouble(parser.getOption("--ring_edge_sigma", "Smoothness sigma of isolation ring", "32"));
		ring_filter_sigma = textToDouble(parser.getOption("--ring_edge_sigma", "Filtering sigma of isolation ring", "128"));
		
		global_prefit = !parser.checkOption("--c2f", "Apply a coarse-to-fine progression instead of initialising with a global pre-fit");
		roundedness = textToDouble(parser.getOption("--rnd", "Roundedness prior", "0"));
		smoothness = textToDouble(parser.getOption("--smt", "Outline average smoothness (negative means full average)", "-1"));
		prior_sigma_A = textToDouble(parser.getOption("--sig", "Weight of initial position", "0"));
		max_binning = textToDouble(parser.getOption("--bin0", "Initial (maximal) binning factor", "8"));
		min_binning = textToDouble(parser.getOption("--bin1", "Final (minimal) binning factor", "4"));
		mask_other_blobs = !parser.checkOption("--nomask", "Do not mask out neighbouring blobs");
		mask_smooth_sigma = textToDouble(parser.getOption("--sig", "Mask smoothing sigma [bin-1 pixels]", "100"));
		        
		diag = parser.checkOption("--diag", "Write out diagnostic information");

		max_frequencies = textToInteger(parser.getOption("--n", "Number of frequencies", "12"));
		blob_thickness = textToDouble(parser.getOption("--th", "Blob thickness [fraction of radius]", "0.25"));
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
	
	do_isolate_ring = ring_min != ring_max;
}

void DeleteBlobs2DProgram::run()
{
	std::vector<std::string> all_micrograph_names;
	all_micrograph_names.reserve(1000);
	
	std::ifstream file(micrographs_list_filename);
	
	if (!file)
	{
		REPORT_ERROR("DeleteBlobs2DProgram::run: unable to read " + micrographs_list_filename);
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
	
	const bool do_phase_flipping = particles_file != "";
	
	
	std::map<std::string, CTF> micrograph_to_ctf;
	
	if (do_phase_flipping)
	{
		Log::print("Will flip the phases of all micrographs according to: "+particles_file);
		
		ObservationModel obs_model;
		MetaDataTable particles_table;
		ObservationModel::loadSafely(particles_file, obs_model, particles_table);
		
		
		std::vector<MetaDataTable> all_tables = StackHelper::splitByMicrographName(particles_table);
		
		for (int i = 0; i < all_tables.size(); i++)
		{
			CTF ctf;
			ctf.readByGroup(all_tables[i], &obs_model, 0);
			
			const std::string micrograph_filename = all_tables[i].getString(EMDL_MICROGRAPH_NAME, 0);
			
			std::string micrograph_name = micrograph_filename.substr(
			            micrograph_filename.find_last_of('/')+1);
			
			micrograph_name = micrograph_name.substr(
			            0, micrograph_name.find_last_of('.'));
			
			micrograph_to_ctf[micrograph_name] = ctf;
		}
	}
	
	
	#pragma omp parallel for num_threads(num_threads)
	for (int m = 0; m < micrograph_count; m++)
	{
		const int thread_id = omp_get_thread_num();
		const bool verbose = thread_id == 0;
		
		if (verbose)
		{
			Log::beginSection(
				"Micrograph "+ZIO::itoa(m+1)+"/"
			    +ZIO::itoa((int)std::ceil(micrograph_count / (double)num_threads))+" on this thread (out of "
			    +ZIO::itoa(micrograph_count)+" on all threads)");
		}
		
		try 
		{
			const std::string micrograph_name = all_micrograph_names[m];
			
			CTF* ctf = 0;
			
			if (do_phase_flipping)
			{
				std::map<std::string, CTF>::iterator it = micrograph_to_ctf.find(micrograph_name);
				
				if (it == micrograph_to_ctf.end())
				{
					Log::warn("Unable to find a CTF for "+micrograph_name+" in "+particles_file);
				}
				else
				{
					ctf = &it->second;
				}
			}
			
			processMicrograph(
				m,
				micrographs_dir + micrograph_name + ".mrc",
				blobs_dir + micrograph_name + ".blobs",
				visualisation,
				visualisation_binning,
				verbose,
				ctf);
		}
		catch (...)
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
        bool verbose,
        CTF* ctf)
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
	
	
	BufferedImage<int> closest_blob = findClosestBlob(delineated_blobs, w_full, h_full);
	        

	if (verbose)
	{
		Log::print("Filtering");
	}
	
	BufferedImage<float> micrograph_filtered;
	
	if (ctf != 0)
	{
		micrograph_filtered = ImageFilter::phaseFlip(micrograph, *ctf, pixel_size);
		
		micrograph_filtered = ImageFilter::highpassStackGaussPadded(
				micrograph_filtered, highpass_sigma_real, 1);
	}
	else
	{
		micrograph_filtered = ImageFilter::highpassStackGaussPadded(
				micrograph, highpass_sigma_real, 1);
	}
	
	BufferedImage<float> blobs_image(w_full, h_full, 1);
	blobs_image.fill(0.f);
	
	BufferedImage<float> erased_image = micrograph;
	
	
	std::string micrograph_name = micrograph_filename.substr(
	            micrograph_filename.find_last_of('/')+1);
	
	micrograph_name = micrograph_name.substr(
	            0, micrograph_name.find_last_of('.'));
	
	
	std::vector<std::vector<double>> all_blob_parameters;


	for (int blob_id = 0; blob_id < blob_count; blob_id++)
	{
		if (verbose)
		{
			Log::beginSection("Blob " + ZIO::itoa(blob_id + 1)+"/"+ZIO::itoa(blob_count));
		}
		
		const double radius = delineated_blobs[blob_id].radius;
		Blob2D initial_blob = delineated_blobs[blob_id].getBlob2D();
		
		std::vector<double> initial_parameters_fullsize = initial_blob.toVector();

		i2Vector window_origin_full;
		BufferedImage<float> blob_region_full;
		BufferedImage<float> blob_mask_full;

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
				if (verbose)
				{
					Log::endSection();
				}
				
				all_blob_parameters.push_back(std::vector<double>(0));
				            
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
			blob_mask_full.resize(window_size_full.x, window_size_full.y);

			// extract full-size image

			for (int y = 0; y < window_size_full.y; y++)
			for (int x = 0; x < window_size_full.x; x++)
			{
				const int x0 = window_origin_full.x;
				const int y0 = window_origin_full.y;

				blob_region_full(x,y) = micrograph_filtered(x0 + x, y0 + y);
			}
			
			if (mask_other_blobs)
			{
				for (int y = 0; y < window_size_full.y; y++)
				for (int x = 0; x < window_size_full.x; x++)
				{
					const int x0 = window_origin_full.x;
					const int y0 = window_origin_full.y;
	
					blob_mask_full(x,y) = closest_blob(x0 + x, y0 + y) == blob_id;
				}
			}
			else
			{
				blob_mask_full.fill(1.f);
			}
		}


		std::vector<double> initial_parameters_cropped = initial_parameters_fullsize;
		initial_parameters_cropped[0] -= window_origin_full.x;
		initial_parameters_cropped[1] -= window_origin_full.y;

		std::vector<double> blob_parameters_cropped = initial_parameters_cropped;

		double estimated_radius = radius;

		if (global_prefit)
		{
			BufferedImage<float> dummy_weight = blob_mask_full;
			dummy_weight.fill(1.f);
			
			if (verbose)
			{
				Log::print("Pre-fitting globally");
			}
			
			std::pair<double, std::vector<double>> global_estimate =
				GlobalBlobFit2D::fit(
					blob_parameters_cropped, 
					radius, max_frequencies,
					blob_region_full,
					dummy_weight);
			
			estimated_radius = global_estimate.first;
			blob_parameters_cropped = global_estimate.second;
			
			if (diag)
			{
				BufferedImage<float> outline = GlobalBlobFit2D::drawOutline(
						blob_parameters_cropped,
						estimated_radius,
						blob_region_full);
				
				
				std::string blobTag = "blob_" + ZIO::itoa(blob_id);
				ZIO::makeOutputDir(outPath + "diag/" + micrograph_name + "/" + blobTag);
				
				outline.write(
					outPath + "diag/" + micrograph_name + "/" +
					blobTag + "/global_initial.mrc");
			}
			
			if (verbose)
			{
				Log::print("Refining at bin " + ZIO::itoa((int)min_binning));
			}
			
			// try much thinner blobs
			blob_parameters_cropped = fitBlob(
					blob_id, blob_parameters_cropped,
					radius, pixel_size, min_binning,
					blob_region_full,
					blob_mask_full,
					micrograph_name,
					verbose);
		}
		else
		{
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
					blob_mask_full,
					micrograph_name,
					verbose);
	
				current_binning /= 2;
				
				if (verbose)
				{
					Log::endSection();
				}
			}
		}
		
		
		if (verbose)
		{
			Log::print("Erasing");
		}
		
		
		std::vector<double> final_parameters_fullsize = blob_parameters_cropped;
		final_parameters_fullsize[0] += window_origin_full.x;
		final_parameters_fullsize[1] += window_origin_full.y;
		
		BufferedImage<float> blob_weight(w_full, h_full, 1);
		
		for (int y = 0; y < h_full; y++)
		for (int x = 0; x < w_full; x++)
		{
			if (closest_blob(x,y) == blob_id)
			{
				blob_weight(x,y) = 1.f;
			}
			else
			{
				blob_weight(x,y) = 0.f;
			}
		}
		
		blob_weight = ImageFilter::Gauss2D(blob_weight, 0, mask_smooth_sigma, true);
		
		Blob2D final_blob(final_parameters_fullsize, radius/2);
		
		if (smoothness < 0)
		{
			final_blob.erase(micrograph, erased_image, blobs_image, blob_weight, 1.5 * radius, radius);
		}
		else
		{
			final_blob.eraseLocally(micrograph, erased_image, blobs_image, blob_weight, 1.5 * radius, radius, smoothness);
		}
		
		const int param_count = final_parameters_fullsize.size();
		std::vector<double> final_parameters_fullsize_with_radius(param_count + 1);
		
		final_parameters_fullsize_with_radius[0] = estimated_radius;
		        
		for (int i = 0; i < param_count; i++)
		{
			final_parameters_fullsize_with_radius[i+1] = final_parameters_fullsize[i];
		}
		
		all_blob_parameters.push_back(final_parameters_fullsize_with_radius);
		        
		if (verbose)
		{
			Log::endSection();
		}
	}
	
	if (do_isolate_ring)
	{
		BufferedImage<float> closest_blob_vis(closest_blob.xdim, closest_blob.ydim);
		
		for (long int i = 0; i < closest_blob_vis.getSize(); i++)
		{
			closest_blob_vis[i] = (float) closest_blob[i];
		}
		
		closest_blob_vis.write("closest_blob_vis.mrc");
		
		erased_image = isolateRing(all_blob_parameters, closest_blob, erased_image, micrograph_name);
	}
	
	if (verbose)
	{
		Log::print("Writing output");
	}

	erased_image.write(outPath + "Frames/" + micrograph_name + ".mrc", pixel_size);
	blobs_image.write(outPath + "Blobs/" + micrograph_name + ".mrc", pixel_size);
	
	ZIO::writeToFile(all_blob_parameters, outPath + "Blobs/" + micrograph_name + ".blobs");
	
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
        BufferedImage<float>& blob_mask_full,
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
	

	const d2Vector initial_position_cropped = d2Vector(
	            initial_blob_params_cropped[0], 
				initial_blob_params_cropped[1]);

	BufferedImage<float> blob_region_binned = Resampling::FourierCrop_fullStack(
				blob_region_full, binning_factor, 1, true);
	
	BufferedImage<float> blob_mask_binned(blob_region_binned.xdim, blob_region_binned.ydim);
	
	for (int y = 0; y < blob_region_binned.ydim; y++)
	for (int x = 0; x < blob_region_binned.xdim; x++)
	{
		blob_mask_binned(x,y) = blob_mask_full(
		            (int)(x*binning_factor), 
		            (int)(y*binning_factor));
	}
	
	
	//blob_region_binned = ImageFilter::Gauss2D(blob_region_binned, 0, 3, true);

	if (diag)
	{
		ZIO::makeOutputDir(outPath + "diag/" + image_name + "/" + blobTag);
		blob_region_binned.write(outTag + "_data.mrc");
	}


	std::vector<double> last_optimum = fromBin1(initial_blob_params_cropped, binning_factor);
	int initial_frequencies = initial_blob_params_cropped.size() < 2? 0 : (initial_blob_params_cropped.size() - 2) / 2;
	
	
	
	const d2Vector initial_position_cropped_binned = initial_position_cropped / binning_factor;

	BlobFit2D blob_fit(
		blob_region_binned, initial_position_cropped_binned, 
		smoothing_radius_binned, prior_sigma, roundedness, 1);
	
	{
		Blob2D blob0(last_optimum, smoothing_radius_binned);
		
		const double r0 = radius_binned - blob_thickness * radius_binned / 2;
		const double r1 = radius_binned + blob_thickness * radius_binned / 2;
		
		blob_fit.computeWeight(blob0, r0, r1, blob_mask_binned);
	}

	for (int current_frequencies = initial_frequencies;
		 current_frequencies <= max_frequencies; current_frequencies++)
	{
		if (verbose)
		{
			Log::print(ZIO::itoa(current_frequencies) + " frequencies ("
				   + ZIO::itoa(2*current_frequencies + 2) + " parameters)");
		}

		std::string tag = outTag + "_N_" + ZIO::itoa(current_frequencies);
				
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
			Blob2D blob0(current_optimum, smoothing_radius_binned);
			drawTestStack(blob0, blob_region_binned, blob_fit.weight).write(tag+"_initial_fit.mrc");
			
			Blob2D blob2(new_optimum, smoothing_radius_binned);
			
			const double E0 = blob_fit.f(current_optimum, 0);
			const double E1 = blob_fit.f(new_optimum, 0);

			if (verbose)
			{
				Log::print("  E: "+ZIO::itoa(E0)+" -> "+ZIO::itoa(E1));
			}

			drawTestStack(blob2, blob_region_binned, blob_fit.weight).write(tag+"_final_fit.mrc");
			drawFit(blob2, blob_region_binned, blob_fit.weight).write(tag+"_residual.mrc");
			blob_fit.weight.write(tag+"_weight.mrc");
		}

		last_optimum = new_optimum;
	}

	return toBin1(last_optimum, binning_factor);
}

BufferedImage<float> DeleteBlobs2DProgram::isolateRing(
        const std::vector<std::vector<double> >& blob_parameters, 
        const BufferedImage<int>& closest_blob,
        const BufferedImage<float>& image,
        const std::string& micrograph_name)
{
	const int w = image.xdim;
	const int h = image.ydim;
	
	
	BufferedImage<float> mask(w,h);
	
	bool warned_already = false;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const int b = closest_blob(x,y);
		
		if (b >= 0)
		{
			if (!warned_already && (b >= blob_parameters.size() || blob_parameters[b].size() == 0))
			{
				std::cout << "micrograph_name: " << micrograph_name << '\n';
				std::cout << "blob_index: " << b << " out of " << blob_parameters.size() << '\n';
				std::cout << "x, y: " << x << ", " << y << '\n' << std::endl;
				
				warned_already = true;
			}
		}
		
		if (b >= 0 && b < blob_parameters.size() && blob_parameters[b].size() > 0)
		{
			DelineatedBlob2D blob(blob_parameters[b]);
			
			const double r = blob.getSignedDistance(d2Vector(x,y));
			
			mask(x,y) = (r >= ring_min && r <= ring_max)? 1.f : 0.f;
		}
		else
		{
			mask(x,y) = 0.f;
		}
	}
	
	BufferedImage<float> masked_image = image * mask;	
	BufferedImage<float> local_average = ImageFilter::Gauss2D(masked_image, 0, ring_filter_sigma, true);
	BufferedImage<float> smooth_mask = ImageFilter::Gauss2D(mask, 0, ring_filter_sigma, true);
	
	smooth_mask = ImageFilter::thresholdAbove(smooth_mask, 0.001f);
	
	local_average /= smooth_mask;
	
	BufferedImage<float> sharp_mask = ImageFilter::Gauss2D(mask, 0, ring_edge_sigma, true);
	
	BufferedImage<float> out = image - local_average;
	
	out *= sharp_mask;
		
	return out;
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

BufferedImage<int> DeleteBlobs2DProgram::findClosestBlob(const std::vector<DelineatedBlob2D> &blobs, int w, int h)
{
	double max_rad = 0;
	
	for (int b = 0; b < blobs.size(); b++)
	{
		const DelineatedBlob2D& blob = blobs[b];
		
		if (blob.radius > max_rad) max_rad = blob.radius;
	}
	
	BufferedImage<int> out(w,h);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const d2Vector pos(x,y);
		
		double min_dist = sqrt(w*w + h*h);
		int best_blob = -1;
		
		for (int b = 0; b < blobs.size(); b++)
		{
			const DelineatedBlob2D& blob = blobs[b];
			
			if ((blob.center - pos).length() < 2.0 * max_rad)
			{
				const double dist = blob.getSignedDistance(pos);
				
				if (dist < min_dist)
				{
					min_dist = dist;
					best_blob = b;
				}
			}
		}
		
		out(x,y) = best_blob;
	}
	
	return out;
}


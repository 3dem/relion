#include <src/jaz/util/topaz_helper.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/util/drawing.h>
#include <src/jaz/util/index_sort.h>
#include <src/jaz/image/local_extrema.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/membrane/blob_2d.h>
#include <src/jaz/membrane/point_blob_fit_2d.h>
#include <src/jaz/membrane/area_point_blob_fit.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	std::string points_file_name, image_directory, outDir;
	double binning, score_threshold, particle_spacing, min_radius_bin1, max_radius_bin1, tolerance_bin1, threshold, tethering, 
	        aspect_cost, contrast_cost, acceptance_threshold;
	int radius_steps, max_iterations, max_frequencies, num_threads;
	bool diag;
		
	
	IOParser parser;
	
	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		points_file_name = parser.getOption("--i", "Topaz membrane detections");
		image_directory = parser.getOption("--id", "Directory containing the micrographs");
		score_threshold = textToDouble(parser.getOption("--pt", "Score threshold for input particles", "-4"));
		particle_spacing = textToDouble(parser.getOption("--sp", "Input particle spacing (bin-1 pixels)", "25"));
		binning = textToDouble(parser.getOption("--bin", "Binning level", "32"));
		min_radius_bin1 = textToDouble(parser.getOption("--r0", "Minimal blob radius (bin-1 pixels)", "300"));
		max_radius_bin1 = textToDouble(parser.getOption("--r1", "Maximal blob radius (bin-1 pixels)", "400"));
		radius_steps = textToInteger(parser.getOption("--rs", "Radius steps", "5"));
		tolerance_bin1 = textToDouble(parser.getOption("--rt", "Radius tolerance (bin-1 pixels)", "50"));
		threshold = textToDouble(parser.getOption("--dt", "Blob centre detection threshold (#particles)", "15"));
		
		max_iterations = textToInteger(parser.getOption("--it", "Max. number of iterations", "1000"));
		max_frequencies = textToInteger(parser.getOption("--frq", "Max. number of blob frequencies", "6"));
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		
		tethering = textToDouble(parser.getOption("--tth", 
			"Blob tethering to its initial position", 
		    "0.0"));
		
		aspect_cost = textToDouble(parser.getOption("--ac", 
			"Cost of deviating from a circular shape", 
		    "0.02"));
		
		contrast_cost = textToDouble(parser.getOption("--cc", 
			"Cost of contrast mismatch (dark pixels outside or bright pixels inside)", 
		    "2.0"));
		
		acceptance_threshold = textToDouble(parser.getOption("--at", 
			"Acceptance threshold for final blobs", 
		    "0.5"));
		
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
	
	if (image_directory.length() > 0)
	{
		if (image_directory[image_directory.length()-1] != '/')
		{
			image_directory = image_directory + "/";
		}
	}
	
	        
	TopazParticleMap particles_by_image = TopazHelper::read(points_file_name, score_threshold);

	const std::string first_image_name = particles_by_image.begin()->first;
	
	t3Vector<long int> full_image_size = ImageFileHelper::getSize(image_directory + first_image_name + ".mrc");
	
	const i2Vector binned_image_size(
	            full_image_size.x / binning,
	            full_image_size.y / binning);
	
	const double avg_radius_bin1 = (min_radius_bin1 + max_radius_bin1) / 2.0;
	const double binned_avg_radius = avg_radius_bin1 / binning;
	const double binned_tolerance = tolerance_bin1 / binning;
	
	
	const int micrograph_count = particles_by_image.size();
	
	std::vector<std::string> all_image_names(micrograph_count);
	std::vector<std::vector<TopazHelper::Particle>> all_particles(micrograph_count);
	
	{
		int index = 0;
		
		for (TopazParticleMap::iterator it = particles_by_image.begin();
		     it != particles_by_image.end(); it++)
		{
			all_image_names[index] = it->first;
			all_particles[index] = it->second;
			
			index++;
		}
	}
	
	{
		std::ofstream names_list(outDir+"micrographs.txt");
		
		for (int m = 0; m < micrograph_count; m++)
		{
			names_list << all_image_names[m] << '\n';
		}
	}
	
	
	BufferedImage<float> diagnostic(binned_image_size.x, binned_image_size.y, micrograph_count);
	
	int res;
	res = system(("mkdir -p "+outDir+"Frames").c_str());
	res = system(("mkdir -p "+outDir+"diag").c_str());
	

	Log::beginProgress("Finding blobs", micrograph_count / num_threads);

	#pragma omp parallel for num_threads(num_threads)	
	for (int m = 0; m < micrograph_count; m++)
	{
		const int th = omp_get_thread_num();

		if (th == 0)
		{
			Log::updateProgress(m);
		}

		const std::string image_name = all_image_names[m];
		const std::vector<TopazHelper::Particle> particles = all_particles[m];
		
		BufferedImage<float> centre_quality(binned_image_size.x, binned_image_size.y);
		centre_quality.fill(0.f);
		
		BufferedImage<float> blob_radius(binned_image_size.x, binned_image_size.y);
		blob_radius.fill(0.f);
		
		BufferedImage<float> micrograph;
		micrograph.read(image_directory + image_name + ".mrc");
		micrograph = Resampling::FourierCrop_fullStack(micrograph, binning, 1, true);		
		micrograph = Normalization::byNormalDist(micrograph);
		
		BufferedImage<float> lowpass0 = ImageFilter::Gauss2D(
					micrograph, 0, 0.125 * binned_avg_radius, true);

		BufferedImage<float> lowpass1 = ImageFilter::Gauss2D(
					micrograph, 0, 1.500 * binned_avg_radius, true);

		BufferedImage<float> dog = lowpass1 - lowpass0;
		
		if (diag)
		{
			dog.write(outDir+"diag/"+image_name+"_DoG.mrc");
			micrograph.write(outDir+"diag/"+image_name+"_micrograph.mrc");
		}

		
		for (int y = 0; y < binned_image_size.y; y++)
		for (int x = 0; x < binned_image_size.x; x++)
		{
			const d2Vector point_pos(x,y);
			
			double max_score = 0.0;
			double best_radius = 0.0;
			
			for (int rr = 0; rr < radius_steps; rr++)
			{
				const double t = rr / (double)(radius_steps - 1);
				const double rad_bin1 = (1 - t) * min_radius_bin1 + t * max_radius_bin1;
				const double rad = rad_bin1 / binning;

				double score = 0.0;
				
				for (int p = 0; p < particles.size(); p++)
				{
					const TopazHelper::Particle& particle = particles[p];
					
					const d2Vector binned_pos(
						particle.coordinates.x / binning,
						particle.coordinates.y / binning);
					
					const double d = (binned_pos - point_pos).length();
					const double dr = (d - rad) / binned_tolerance;
					
					score += exp(-dr*dr/2);
				}
				
				if (score > max_score)
				{
					best_radius = rad;
					max_score = score;
				}
				
			}
			
			centre_quality(x,y) = (float) max_score;
			blob_radius(x,y) = (float) best_radius;
		}
		
		if (diag)
		{
			centre_quality.write(outDir + "centre_quality.mrc");
		}
		
		BufferedImage<float> box_maxima = LocalExtrema::boxMaxima(centre_quality, (int)(0.72 * binned_avg_radius));
		std::vector<d2Vector> detections = LocalExtrema::discretePoints2D(centre_quality, box_maxima, (float)threshold);
		
		
		
		std::vector<d2Vector> all_particle_positions(particles.size()); 
		
		for (int i = 0; i < particles.size(); i++)
		{
			i2Vector coords = particles[i].coordinates;
			all_particle_positions[i] = d2Vector(coords.x, coords.y) / binning;
		}
		
		
		if (diag)
		{
			BufferedImage<float> detections_image(binned_image_size.x, binned_image_size.y);
			detections_image.fill(0.f);
			Drawing::drawCrosses(detections, 1.f, 10, detections_image);		
			detections_image.write(outDir + "topaz_blob_detections.mrc");
			detections_image.fill(0.f);
			Drawing::drawCrosses(all_particle_positions, 1.f, 3, detections_image);		
			detections_image.write(outDir + "topaz_input.mrc");
		}
		
		
		BufferedImage<float> density_map(binned_image_size.x, binned_image_size.y);
		
		for (int y = 0; y < binned_image_size.y; y++)
		for (int x = 0; x < binned_image_size.x; x++)
		{
			density_map(x,y) = 0.f;
			
			const double sigma = 2 * particle_spacing / binning;
			
			for (int p = 0; p < all_particle_positions.size(); p++)
			{
				const double d = (all_particle_positions[p] - d2Vector(x,y)).length() / sigma;
				const double rho = exp(-d*d / 2);
				
				if (density_map(x,y) < rho)
				{
					density_map(x,y) = rho;
				}
			}
		}
		
		if (diag)
		{
			density_map.write(outDir + "density_map.mrc");
		}
		
		
		BufferedImage<float> final_plots_sum(binned_image_size.x, binned_image_size.y);		
		final_plots_sum.fill(0.f);
		
		BufferedImage<float> initial_plots_sum;	
		
		if (diag)
		{
			initial_plots_sum.resize(binned_image_size.x, binned_image_size.y);
			initial_plots_sum.fill(0.f);
		}
		
		std::vector<double> all_costs(detections.size());
		std::vector<std::vector<double>> all_optimal_parameters(detections.size());
		std::vector<bool> is_accepted(detections.size());
		
		for (int detection_id = 0; detection_id < detections.size(); detection_id++)
		{
			const d2Vector d = detections[detection_id];
			
			const double radius = blob_radius((int)std::round(d.x), (int)std::round(d.y));
			
			AreaPointBlobFit point_blob_fit(
				density_map, dog, d, radius, binned_tolerance, tethering, aspect_cost, contrast_cost);
			            
			/*PointBlobFit2D point_blob_fit(
				all_particle_positions, d, radius, binned_tolerance, tethering);*/
			
			std::vector<double> initial_parameters(3+2*max_frequencies, 0.0);
			initial_parameters[0] = radius;
			initial_parameters[1] = d.x;
			initial_parameters[2] = d.y;
			
			if (diag)
			{
				BufferedImage<float> final_plot = point_blob_fit.visualise(
			            initial_parameters, binned_image_size.x, binned_image_size.y);
				
				initial_plots_sum += final_plot;
			}
			
						
			std::vector<double> optimal_parameters = NelderMead::optimize(
			            initial_parameters, point_blob_fit, 0.5, 0.0001, max_iterations, 
			            1.0, 2.0, 0.5, 0.5, false);
			
			const double cost = point_blob_fit.f(optimal_parameters,0);
			
			
			is_accepted[detection_id] = cost < acceptance_threshold;
			
			if (!(cost == cost))
			{
				std::cout << "nan: " << d << ", " << radius << "\n";
				std::cout << "optimal_parameters.size() = " << optimal_parameters.size() << "\n";
				
				for (int i = 0; i < optimal_parameters.size(); i++)
				{
					std::cout << optimal_parameters[i] << " ";
				}
				std::cout << "\n";
			}
			
			if (is_accepted[detection_id])
			{
				BufferedImage<float> final_plot = point_blob_fit.visualise(
			            optimal_parameters, binned_image_size.x, binned_image_size.y);
				
				final_plots_sum += final_plot;
			}
			
			all_costs[detection_id] = cost;
			all_optimal_parameters[detection_id] = optimal_parameters;
		}
		
		
		BufferedImage<float> all_particle_dots(binned_image_size.x, binned_image_size.y);	
		Drawing::drawCrosses(all_particle_positions, 0.1f, 3, all_particle_dots);
		final_plots_sum += all_particle_dots;
		
		
		std::vector<int> sorted_blobs = IndexSort<double>::sortIndices(all_costs);
		
		
		std::ofstream blob_file(outDir+"Frames/"+image_name+".blobs");

		for (int i = 0; i < sorted_blobs.size(); i++)
		{
			if (is_accepted[sorted_blobs[i]])
			{
				const std::vector<double> parameters = all_optimal_parameters[sorted_blobs[i]];
				
				for (int j = 0; j < parameters.size(); j++)
				{
					blob_file << binning * parameters[j] << " ";
				}
				
				blob_file << '\n';
				
				Drawing::drawCross(
					d2Vector(parameters[1], parameters[2]),
					(float)all_costs[sorted_blobs[i]], 5, final_plots_sum);
			}
		}
		
		diagnostic.getSliceRef(m).copyFrom(final_plots_sum);
		
		if (diag)
		{
			initial_plots_sum.write(outDir+"diag/"+image_name+"_initial.mrc");
			final_plots_sum.write(outDir+"diag/"+image_name+"_final.mrc");
		}
	}

	Log::endProgress();
	
	diagnostic.write(outDir + "diagnostic.mrc");
	
	return 0;
}

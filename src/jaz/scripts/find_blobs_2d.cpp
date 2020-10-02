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
#include <omp.h>


using namespace gravis;


int main(int argc, char *argv[])
{
	std::string particlesFn, class_averages_filename, outDir;
	int num_threads, max_class, min_MG, max_MG;
	double radius, binning, additional_binning;
	bool flip_y, diag;
	
	
	IOParser parser;
	
	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input file (e.g. run_it023_data.star)");
		class_averages_filename = parser.getOption("--ca", "Class averages stack");
		radius = textToDouble(parser.getOption("--r", "Average blob radius (in bin-1 pixels)", "300"));
		flip_y = parser.checkOption("--flip_y", "Centre of blob lies in positive y direction");

		
		max_class = textToInteger(parser.getOption("--max_class", "Last 2D class to consider", "25"));
		min_MG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
		max_MG = textToInteger(parser.getOption("--max_MG", "Last micrograph index", "-1"));
		
		binning = textToDouble(parser.getOption("--bin", "Binning level to work at", "8"));
		additional_binning = textToDouble(parser.getOption("--bin2", "Additional binning for the log file", "4"));
		
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
	
	
	std::vector<MetaDataTable> particles_by_micrograph = StackHelper::splitByMicrographName(particles_table);
	
	const int micrograph_count = particles_by_micrograph.size();
	
	if (max_MG < 0) 
	{
		max_MG = micrograph_count - 1;
	}
	
	i2Vector full_size, binned_size, log_size;	
	{
		
		const std::string micrograph_fn = particles_by_micrograph[0].getString(EMDL_MICROGRAPH_NAME, 0);
		
		BufferedImage<float> micrograph;
		micrograph.read(micrograph_fn);
		
		full_size.x = micrograph.xdim;
		full_size.y = micrograph.ydim;
		
		BufferedImage<float> micrograph_binned = Resampling::FourierCrop_fullStack(
					micrograph, binning, num_threads, true);
		
		binned_size.x = micrograph_binned.xdim;
		binned_size.y = micrograph_binned.ydim;
		
		BufferedImage<float> micrograph_binned_2 = Resampling::FourierCrop_fullStack(
					micrograph_binned, additional_binning, num_threads, true);
		
		log_size.x = micrograph_binned_2.xdim;
		log_size.y = micrograph_binned_2.ydim;
		
	}
	
	const int w_log = log_size.x;
	const int h_log = log_size.y;
	
	const int w_binned = binned_size.x;
	const int h_binned = binned_size.y;
	
	BufferedImage<float> log_image(w_log, h_log, max_MG + 1);
	log_image.fill(0.f);
	
	

	Log::beginProgress("Finding blobs", (max_MG - min_MG + 1)/num_threads);

	std::vector<std::vector<d2Vector>> blob_coordinates(max_MG + 1);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int m = min_MG; m <= max_MG; m++)
	{
		const int thread_id = omp_get_thread_num();

		if (thread_id == 0)
		{
			Log::updateProgress(m);
		}
		
		MetaDataTable& particles = particles_by_micrograph[m];
		
		const std::string micrograph_fn = particles.getString(EMDL_MICROGRAPH_NAME, 0);
		
		BufferedImage<float> micrograph;
		micrograph.read(micrograph_fn);
		
		BufferedImage<float> micrograph_binned = Resampling::FourierCrop_fullStack(
					micrograph, binning, 1, true);
		
		if (diag)
		{
			micrograph_binned.write(outDir + "mg_" + ZIO::itoa(m) + "_micrograph_binned.mrc");
		}
		
		BufferedImage<float> accumulation_image(w_binned, h_binned);
		
		accumulation_image.fill(0.f);
		
		
		const int particle_count = particles.numberOfObjects();
		
		for (int p = 0; p < particle_count; p++)
		{
			const int class_id = particles.getIntMinusOne(EMDL_PARTICLE_CLASS, p);
			
			if (class_id <= max_class)
			{
				const d2Vector image_coord(
					particles.getDouble(EMDL_IMAGE_COORD_X, p),
					particles.getDouble(EMDL_IMAGE_COORD_Y, p));
			
				const double phi = DEG2RAD(particles.getDouble(EMDL_ORIENT_PSI, p));
				     
				const d2Vector predicted_centre = flip_y?
						(image_coord + radius * d2Vector(sin(phi), cos(phi))) / binning :
						(image_coord - radius * d2Vector(sin(phi), cos(phi))) / binning;
				
				const int xi = (int)(predicted_centre.x + 0.5);
				const int yi = (int)(predicted_centre.y + 0.5);
				
				if (xi >= 0 && xi < w_binned && yi >= 0 && yi < h_binned)
				{
					accumulation_image(xi, yi) += 1;
				}                      
			}
		}
		
		BufferedImage<float> blurred_accumulation_image = ImageFilter::Gauss2D(
					accumulation_image, 0, 0.5 * radius / binning, true);
		
		if (diag)
		{
			blurred_accumulation_image.write(outDir + "mg_" + ZIO::itoa(m) + "_likelihood_by_particles.mrc");
		}
		
		BufferedImage<float> blurred_micrograph_0 = ImageFilter::Gauss2D(
					micrograph_binned, 0, 1.00 * radius / binning, true);
		
		BufferedImage<float> blurred_micrograph_1 = ImageFilter::Gauss2D(
					micrograph_binned, 0, 0.25 * radius / binning, true);
		
		BufferedImage<float> dog_map = blurred_micrograph_0 - blurred_micrograph_1;
		
		if (diag)
		{
			dog_map.write(outDir + "mg_" + ZIO::itoa(m) + "_likelihood_by_appearance.mrc");
		}
		
		BufferedImage<float> detections(w_binned, h_binned);
		
		for (int y = 0; y < h_binned; y++)
		for (int x = 0; x < w_binned; x++)
		{
			const float dog = dog_map(x,y);
			const float centre = blurred_accumulation_image(x,y);
			
			if (dog > 0.f && centre > 0.f)
			{
				detections(x,y) = dog * centre;
			}
		}
		
		if (diag)
		{
			detections.write(outDir + "mg_" + ZIO::itoa(m) + "_likelihood_total.mrc");
		}
		
		BufferedImage<float> box_maxima = LocalExtrema::boxMaxima(detections, radius / binning);
		
		const float mean_DOG = Normalization::computeMean(detections);
		const float threshold = mean_DOG;
		
		BufferedImage<float> sparse_detections_image(w_binned, h_binned);
		std::vector<d2Vector> sparse_detections;
		sparse_detections.reserve(100);
		        
		for (int y = 0; y < h_binned; y++)
		for (int x = 0; x < w_binned; x++)
		{
			const float bm = box_maxima(x,y);
			const float dt = detections(x,y);
			
			if (bm == dt && dt > threshold)
			{
				sparse_detections_image(x,y) = dt;
				sparse_detections.push_back(binning * d2Vector(x,y));
			}
			else
			{
				sparse_detections_image(x,y) = 0.f;
			}
		}
		
		blob_coordinates[m] = sparse_detections;
		
		BufferedImage<float> micrograph_binned_2 = Resampling::FourierCrop_fullStack(
		            micrograph_binned, additional_binning, 1, true);
		
		const float mean = Normalization::computeMean(micrograph_binned_2);
		const float variance = Normalization::computeVariance(micrograph_binned_2, mean);
		const float dot_value = mean + 7.f * sqrt(variance);
		
		RawImage<float> log_slice = log_image.getSliceRef(m);
		log_slice.copyFrom(micrograph_binned_2);
		
		for (int i = 0; i < sparse_detections.size(); i++)
		{
			const d2Vector d = sparse_detections[i] / (binning * additional_binning);
			const int dx = (int)(d.x + 0.5);
			const int dy = (int)(d.y + 0.5);
			
			if (dx >= 0 && dx < w_log && dy >= 0 && dy < h_log)
			{
				log_slice(dx,dy) = dot_value;
			}
		}
	}

	Log::endProgress();
	
	log_image.write(outDir+"diagnostic.mrc");
	
	{
		const int mc = blob_coordinates.size();
		
		std::ofstream ofs(outDir+"blobs.star");
	    
	    for (int m = 0; m < mc; m++)
	    {
			MetaDataTable mdt;
	        mdt.setName(ZIO::itoa(m));
			
			const int bc = blob_coordinates[m].size();
	
	        for (int b = 0; b < bc; b++)
	        {
				const d2Vector p = blob_coordinates[m][b];
				
	            mdt.addObject();
				
	            mdt.setValue(EMDL_IMAGE_COORD_X, p.x);
				mdt.setValue(EMDL_IMAGE_COORD_Y, p.y);
	        }
	
	        mdt.write(ofs);
		}
	}
	
	return 0;
	
	
	/*
	  
	{
		const double bin = 32;
		const double radius_bin1 = 350;
		const double radius_binned = radius_bin1 / bin;
		const double threshold = 0.1;

		BufferedImage<float> micrograph_binned = Resampling::FourierCrop_fullStack(
					micrograph_filtered, bin, num_threads, true);

		micrograph_binned.write(outPath+"DEBUG_micrograph_binned.mrc");

		BufferedImage<float> symm = evaluateRotationalSymmetry(
					micrograph_binned, radius_binned, 1.5*radius_binned, 0.25*radius_binned);

		float mean = Normalization::computeMean(symm);

		symm -= mean;

		symm.write(outPath+"DEBUG_symmetry.mrc");


		BufferedImage<float> lowpass0 = ImageFilter::Gauss2D(
					micrograph_binned, 0, 0.5 * radius_binned, true);

		BufferedImage<float> lowpass1 = ImageFilter::Gauss2D(
					micrograph_binned, 0, 1.500 * radius_binned, true);

		BufferedImage<float> dog = lowpass1 - lowpass0;

		dog.write(outPath+"DEBUG_dog.mrc");

		double var_dog = Normalization::computeVariance(dog, 0.f);
		double var_symm = Normalization::computeVariance(symm, 0.f);

		dog /= sqrt(var_dog);
		symm /= sqrt(var_symm);

		dog = ImageFilter::thresholdAbove(dog, 0.f);
		symm = ImageFilter::thresholdAbove(symm, 0.f);

		BufferedImage<float> product = dog * symm;

		product.write(outPath+"DEBUG_product.mrc");

		BufferedImage<float> box_maxima = LocalExtrema::boxMaxima(
					product, (int)(0.72 * radius_binned));

		box_maxima.write(outPath+"DEBUG_product_maxima.mrc");

		std::vector<d2Vector> peaks = LocalExtrema::discretePoints2D(
				product, box_maxima, (float)threshold);

		const float mean_mg = Normalization::computeMean(micrograph_binned);
		const float var_mg = Normalization::computeVariance(micrograph_binned, mean_mg);
		const float drawing_value = mean_mg + 5 * sqrt(var_mg);

		Drawing::drawCrosses(peaks, drawing_value, 5, micrograph_binned);

		micrograph_binned.write(outPath+"DEBUG_detections.mrc");

		std::exit(0);
	}  

	*/
}

/*
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
}*/

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
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/translation.h>
#include <src/jaz/image/tapering.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/image/resampling.h> 
#include <src/jaz/image/local_extrema.h>
#include <src/jaz/math/fft.h>
#include <omp.h>


using namespace gravis;


BufferedImage<float> rotate(const RawImage<float>& img, double phi);
BufferedImage<float> flip_x(const RawImage<float>& img);
BufferedImage<float> correlate(const RawImage<float>& img0, const RawImage<float>& img1);
float max_correlation(const RawImage<float>& img0, const RawImage<float>& img1);

int main(int argc, char *argv[])
{
	std::string particlesFn, class_averages_filename, outDir;
	int num_threads;
	
	
	IOParser parser;
	
	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input file (e.g. run_it023_data.star)");
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
	ObservationModel::loadSafely(particlesFn, obs_model, particles_table);
	
	const double binning = 8.0;
	const double radius = 300.0;
	const int max_class = 100;
	
	std::vector<MetaDataTable> particles_by_micrograph = StackHelper::splitByMicrographName(particles_table);
	
	const int micrograph_count = particles_by_micrograph.size();
	
	for (int m = 0; m < 1; m++)
	{
		MetaDataTable& particles = particles_by_micrograph[m];
		
		const std::string micrograph_fn = particles.getString(EMDL_MICROGRAPH_NAME, 0);
		
		BufferedImage<float> micrograph;
		micrograph.read(micrograph_fn);
		
		BufferedImage<float> micrograph_binned = Resampling::FourierCrop_fullStack(
		            micrograph, binning, num_threads, true);
		
		
		micrograph_binned.write(outDir + "mg_" + ZIO::itoa(m) + "_micrograph_binned.mrc");
		
		const int w_binned = micrograph_binned.xdim;
		const int h_binned = micrograph_binned.ydim;
		
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
				     
				const d2Vector predicted_centre = 
						(image_coord - radius * d2Vector(sin(phi), cos(phi))) / binning;
				
				const int xi = (int)(predicted_centre.x + 0.5);
				const int yi = (int)(predicted_centre.y + 0.5);
				
				if (xi >= 0 && xi < w_binned && yi >= 0 && yi < h_binned)
				{
					accumulation_image(xi, yi) += 1;
				}
				
				/*
				for (int r = 0; r < 101; r++)
				{
					const d2Vector predicted_centre = 
							(image_coord - 0.01 * r * radius * d2Vector(sin(phi), cos(phi))) / binning;
					
					const int xi = (int)(predicted_centre.x + 0.5);
					const int yi = (int)(predicted_centre.y + 0.5);
					
					if (xi >= 0 && xi < w_binned && yi >= 0 && yi < h_binned)
					{
						accumulation_image(xi, yi) += 1;
					}
				}*/
				                        
			}
		}
		
		accumulation_image.write(outDir + "mg_" + ZIO::itoa(m) + "_accumulation_image.mrc");
		
		BufferedImage<float> blurred_accumulation_image = ImageFilter::Gauss2D(
					accumulation_image, 0, 0.5 * radius / binning, true);
		
		blurred_accumulation_image.write(outDir + "mg_" + ZIO::itoa(m) + "_blurred_accumulation_image.mrc");
		
		BufferedImage<float> blurred_micrograph_0 = ImageFilter::Gauss2D(
					micrograph_binned, 0, 1.0 * radius / binning, true);
		
		BufferedImage<float> blurred_micrograph_1 = ImageFilter::Gauss2D(
					micrograph_binned, 0, 0.5 * radius / binning, true);
		
		blurred_micrograph_0.write(outDir + "mg_" + ZIO::itoa(m) + "_blurred_micrograph_0.mrc");
		blurred_micrograph_1.write(outDir + "mg_" + ZIO::itoa(m) + "_blurred_micrograph_1.mrc");
		
		BufferedImage<float> dog_map = blurred_micrograph_0 - blurred_micrograph_1;
		dog_map.write(outDir + "mg_" + ZIO::itoa(m) + "_DOG.mrc");
		
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
		
		detections.write(outDir + "mg_" + ZIO::itoa(m) + "_detections.mrc");
		
		BufferedImage<float> box_maxima = LocalExtrema::boxMaxima(detections, radius / binning);
		
		BufferedImage<float> sparse_detections(w_binned, h_binned);
		        
		for (int y = 0; y < h_binned; y++)
		for (int x = 0; x < w_binned; x++)
		{
			const float bm = box_maxima(x,y);
			const float dt = detections(x,y);
			
			if (bm == dt)
			{
				sparse_detections(x,y) = dt;
			}
			else
			{
				sparse_detections(x,y) = 0.f;
			}
		}
		
		sparse_detections.write(outDir + "mg_" + ZIO::itoa(m) + "_sparse_detections.mrc");
		
	}
	
}

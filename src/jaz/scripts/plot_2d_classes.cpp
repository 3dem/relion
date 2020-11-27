#include <src/jaz/image/buffered_image.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/class_helper.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/args.h>
#include <src/metadata_table.h>
#include <omp.h>
#include <string>


using namespace gravis;


int main(int argc, char *argv[])
{
	std::string particlesFn, class_averages_filename, classes_filename, micrograph_dir, outDir;
	int num_threads, num_MG, num_best_classes, pad;
	bool flip_contrast, diag;


	IOParser parser;

	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		particlesFn = parser.getOption("--i", "Input STAR file output by align_2d_classes");
		class_averages_filename = parser.getOption("--ca", "Class averages stack");
		micrograph_dir = parser.getOption("--mgdir", "Micrographs directory", "");
		classes_filename = parser.getOption("--classes", "File with a list of 2D classes to consider", "");
		num_best_classes = textToInteger(parser.getOption("--bc", "Number of best 2D classes to consider otherwise", "20"));
		num_MG = textToInteger(parser.getOption("--mgs", "Number of micrographs to use", "12"));
		pad = textToInteger(parser.getOption("--pad", "Image padding (pixels)", "20"));
		num_threads = textToInteger(parser.getOption("--j", "Number of OMP threads", "6"));
		flip_contrast = parser.checkOption("--keep_contrast", "Do not flip the contrast");
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
	const int midbox = box_size/2;
	const double radius = midbox - pad;
	const int num_classes = class_averages.zdim;

	std::vector<int> classes_to_consider;
	
	if (classes_filename != "")
	{
		classes_to_consider = ZIO::readInts(classes_filename);
	}
	else	
	{
		const int class_count = ClassHelper::countClasses(particles_table);
		
		if (class_count < num_best_classes)
		{
			REPORT_ERROR_STR("Cannot consider " << num_best_classes 
			                 << " 2D classes, only " << class_count << " present");
		}
		
		std::vector<int> particle_count = ClassHelper::getClassSizes(particles_table, class_count);
		std::vector<int> order = ClassHelper::sortByAscendingFrequency(particle_count);
		
		classes_to_consider = std::vector<int>(num_best_classes);
		
		for (int i = 0; i < num_best_classes; i++)
		{
			classes_to_consider[i] = order[i];
		}
	}
	
	const int relevant_class_count = classes_to_consider.size();

	std::vector<int> class_to_subset(num_classes, -1);

	for (int i = 0; i < relevant_class_count; i++)
	{
		class_to_subset[classes_to_consider[i]] = i;
	}

	
	const double particle_pixel_size = obs_model.getPixelSize(0);
	

	std::vector<MetaDataTable> particles_by_micrograph = StackHelper::splitByMicrographName(particles_table);
	const int micrograph_count = std::min(num_MG, (int)particles_by_micrograph.size());
	
	
	t3Vector<long int> mg_size = ImageFileHelper::getSize(
	            particles_by_micrograph[0].getString(EMDL_MICROGRAPH_NAME, 0));
		
	BufferedImage<float> output(mg_size.x, mg_size.y, 2 * micrograph_count);
	
	const float scale = flip_contrast? -1.f : 1.f;

	
	Log::beginProgress("Plotting 2D class averages", micrograph_count);

	for (int m = 0; m < micrograph_count; m++)
	{
		Log::updateProgress(m);
		
		if (particles_by_micrograph[m].numberOfObjects() == 0) continue;
		
		MetaDataTable& particles = particles_by_micrograph[m];
		        
		std::string micrograph_filename = particles.getString(EMDL_MICROGRAPH_NAME, 0);
		        
		BufferedImage<float> micrograph;
		micrograph.read(micrograph_filename);
		
		const float average_value = Normalization::computeMean(micrograph);
		
		output.getSliceRef(2*m).copyFrom(micrograph);
		
		const int w = micrograph.xdim;
		const int h = micrograph.ydim;
		
		const double micrograph_pixel_size = ImageFileHelper::getSamplingRate(micrograph_filename);
		const double particle_scale = particle_pixel_size / micrograph_pixel_size;
		
		const int pc = particles.numberOfObjects();
		
		for (int p = 0; p < pc; p++)
		{
			const int global_class_id = particles.getIntMinusOne(EMDL_PARTICLE_CLASS, p);
			const int selected_class_id = class_to_subset[global_class_id];

			if (selected_class_id >= 0)
			{
				const double dx_A = particles.getDouble(EMDL_ORIENT_ORIGIN_X_ANGSTROM, p);
				const double dy_A = particles.getDouble(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, p);
				const double phi = DEG2RAD(particles.getDouble(EMDL_ORIENT_PSI, p));
				
				const double xc = particles.getDouble(EMDL_IMAGE_COORD_X, p) - dx_A / micrograph_pixel_size;
				const double yc = particles.getDouble(EMDL_IMAGE_COORD_Y, p) - dy_A / micrograph_pixel_size;
				
				const int x0 = std::max(std::ceil(xc - particle_scale * box_size), 0.0);
				const int y0 = std::max(std::ceil(yc - particle_scale * box_size), 0.0);
				const int x1 = std::min(std::ceil(xc + particle_scale * box_size), (double) w);
				const int y1 = std::min(std::ceil(yc + particle_scale * box_size), (double) h);
				
				for (int y = y0; y < y1; y++)
				for (int x = x0; x < x1; x++)
				{
					const double xm = x - xc;
					const double ym = y - yc;
										
					const double xp = (cos(phi) * xm - sin(phi) * ym) / particle_scale;
					const double yp = (sin(phi) * xm + cos(phi) * ym) / particle_scale;
					
					if (xp * xp + yp * yp < radius * radius)
					{
						micrograph(x,y) = average_value + scale * Interpolation::linearXY_clip(
								class_averages.getSliceRef(global_class_id),
								midbox + xp,
								midbox + yp);
					}
				}
			}
		}
		
		output.getSliceRef(2*m+1).copyFrom(micrograph);
	}
	
	Log::endProgress();
	
	output.write(outDir+"class_average_plots.mrc");
	
	return 0;
}

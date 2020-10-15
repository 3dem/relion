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
	double binning;
	int num_threads;
		
	
	IOParser parser;
	
	try
	{
		IOParser parser;

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		points_file_name = parser.getOption("--i", "Topaz membrane detections");
		image_directory = parser.getOption("--id", "Directory containing the micrographs");
		binning = textToDouble(parser.getOption("--bin", "Binning level", "32"));
		
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
	
	if (image_directory.length() > 0)
	{
		if (image_directory[image_directory.length()-1] != '/')
		{
			image_directory = image_directory + "/";
		}
	}
	
	        
	TopazParticleMap particles_by_image = TopazHelper::read(points_file_name, -1000);

	const std::string first_image_name = particles_by_image.begin()->first;
	
	
	BufferedImage<float> micrograph;
	micrograph.read(image_directory + first_image_name + ".mrc");
	micrograph = Resampling::FourierCrop_fullStack(micrograph, binning, 1, true);	
	
	/*const i2Vector binned_image_size(
	            full_image_size.x / binning,
	            full_image_size.y / binning);*/
		
	const i2Vector binned_image_size(micrograph.xdim, micrograph.ydim);
	        
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
			names_list << m << " " << all_image_names[m] << '\n';
		}
	}
	
	
	BufferedImage<float> all_particle_images(binned_image_size.x, binned_image_size.y, micrograph_count);
	BufferedImage<float> all_micrographs(binned_image_size.x, binned_image_size.y, micrograph_count);
	
	
	int res = system(("mkdir -p "+outDir+"Frames").c_str());
	

	Log::beginProgress("Plotting micrographs and particles", micrograph_count / num_threads);

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
		
		all_micrographs.getSliceRef(m).copyFrom(micrograph);
		
		
		BufferedImage<float> particles_image(binned_image_size.x, binned_image_size.y);
		particles_image.fill(0.f);
		
		
		for (int i = 0; i < particles.size(); i++)
		{
			TopazHelper::Particle p = particles[i];
			i2Vector c = p.coordinates;
			
			Drawing::drawCross(d2Vector(c.x, c.y) / binning, (float)p.score, 2, particles_image);
		}
		
		all_particle_images.getSliceRef(m).copyFrom(particles_image);
	}
	
	Log::endProgress();
	
	all_particle_images.write(outDir + "particles.mrc");
	all_micrographs.write(outDir + "micrographs.mrc");
	
	return 0;
}

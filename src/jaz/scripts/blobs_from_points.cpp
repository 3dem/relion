#include <src/jaz/util/topaz_helper.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/image_file_helper.h>
#include <src/jaz/membrane/blob_2d.h>

using namespace gravis;

int main(int argc, char *argv[])
{
	std::string points_file_name = "model_dilute_mg0.txt";
	std::string image_directory = "MotionCorr/job002/Frames/";
	const double binning = 32.0;
	const double radius = 350.0;
	const double tolerance = 100.0;
	        
	TopazParticleMap particles_by_image = 
	        TopazHelper::read(points_file_name);

	const std::string first_image_name = particles_by_image.begin()->first;
	
	t3Vector<long int> full_image_size = ImageFileHelper::getSize(image_directory + first_image_name + ".mrc");
	
	const i2Vector binned_image_size(
	            full_image_size.x / binning,
	            full_image_size.y / binning);
	
	const double binned_radius = radius / binning;
	const double binned_tolerance = tolerance / binning;
	
	
	for (TopazParticleMap::iterator it = particles_by_image.begin();
	     it != particles_by_image.end(); it++)
	{
		const std::string image_name = it->first;
		const std::vector<TopazHelper::Particle> particles = it->second;
		
		BufferedImage<float> map(binned_image_size.x, binned_image_size.y);
		map.fill(0.f);
		
		for (int y = 0; y < binned_image_size.y; y++)
		for (int x = 0; x < binned_image_size.x; x++)
		{
			const d2Vector point_pos(x,y);
			
			double score = 0.0;
			
			for (int p = 0; p < particles.size(); p++)
			{
				const TopazHelper::Particle& particle = particles[p];
				
				const d2Vector binned_pos(
					particle.coordinates.x / binning,
					particle.coordinates.y / binning);
				
				const double d = ((binned_pos - point_pos).length() - binned_radius) / binned_tolerance;
				
				score += exp(-d*d/2);
			}
			
			map(x,y) = (float) score;
		}
		
		map.write("DEV_topaz_map.mrc");
		std::exit(0);
	}
	
	
	
	return 0;
}

#ifndef FILAMENT_MAPPING_H
#define FILAMENT_MAPPING_H

#include <src/jaz/image/buffered_image.h>
#include <string>

class FilamentMapping
{
	public:
		
		BufferedImage<float> index, signed_dist, mask;	
		
		void write(std::string filename);
		static FilamentMapping read(std::string filename);
};

#endif

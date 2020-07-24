#ifndef IMAGE_FILE_HELPER_H
#define IMAGE_FILE_HELPER_H

#include <src/jaz/gravis/t3Vector.h>

class ImageFileHelper
{
	public:
		
		static gravis::t3Vector<long int> getSize(const std::string& filename);

		static double getSamplingRate(const std::string& filename);
};

#endif

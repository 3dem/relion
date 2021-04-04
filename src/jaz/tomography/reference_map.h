#ifndef TOMO_REFERENCE_MAP_H
#define TOMO_REFERENCE_MAP_H

#include <src/jaz/image/buffered_image.h>
#include <vector>

class OptimisationSet;

class TomoReferenceMap
{
	public:

		static void presharpen(BufferedImage<float>& map_RS, double padding);
		

		TomoReferenceMap();

		void read(IOParser& parser);
		void read(const OptimisationSet& optimisationSet);
		
			std::string mapFilenames[2], maskFilename, fscFilename;
			bool useFscThreshold;
			double fscThresholdWidth;
			int lastShell;

			std::vector<BufferedImage<float>> image_real;
			BufferedImage<float> mask, freqWeight;
			std::vector<BufferedImage<fComplex>> image_FS;


		void load(int boxSize = -1, int verbosity = 1);
			
		int getBoxSize() const;

};

#endif

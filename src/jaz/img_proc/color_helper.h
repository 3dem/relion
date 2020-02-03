#ifndef COLOR_HELPER_H
#define COLOR_HELPER_H

#include <src/jaz/gravis/tRGB.h>
#include <src/image.h>

class ColorHelper
{
	public:

		static gravis::dRGB signedToRedBlue(double d, double scale = 1.0, double rbFract = 0.333);
		static void writeAngleToPNG(const Image<RFLOAT>& img, std::string filename);
		static void writeSignedToPNG(const Image<RFLOAT>& img, std::string filename, double scale = 1.0);
		static void writeSignedToEPS(std::string filename, int col, const std::vector<Image<RFLOAT> > &imgs,
				const std::vector<double> &scales, const std::vector<std::string> &labels);
};

#endif

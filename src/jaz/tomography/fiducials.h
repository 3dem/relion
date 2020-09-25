#ifndef FIDUCIALS_H
#define FIDUCIALS_H

#include <src/jaz/gravis/t3Vector.h>
#include <vector>
#include "tomogram.h"
#include <src/jaz/image/detection.h>
#include <src/jaz/image/similarity.h>
#include <src/jaz/math/fft.h>

class Fiducials
{
	public:

		static std::vector<gravis::d3Vector> read(
		        const std::string &filename,
				double pixelSize);

		static std::string write(
				const std::vector<gravis::d3Vector>& positions,
				double pixelSize,
				const std::string& tomoName,
				const std::string &path);

		static void drawMask(
				const std::vector<gravis::d3Vector>& positions,
				const gravis::d4Matrix& proj,
				double radius,
				RawImage<float>& destination,
				double value);
		
		static void erase(
		        const std::vector<gravis::d3Vector>& positions,
		        double radius,
		        Tomogram& tomogram,
		        int num_threads);
};


#endif

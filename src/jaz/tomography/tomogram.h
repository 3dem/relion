#ifndef DYN_TOMOGRAM_H
#define DYN_TOMOGRAM_H

#include <src/jaz/image/buffered_image.h>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/optics/optics_data.h>
#include <src/ctf.h>
#include "tomolist.h"



class Tomogram
{
	public:
		
		Tomogram();
			
			
			bool hasOptics, hasImage;
			OpticsData optics;
			int frameCount;
			double handedness;
			
			BufferedImage<float> stack;
			std::vector<gravis::d4Matrix> projectionMatrices;
			
			std::vector<CTF> centralCTFs;
			std::vector<double> cumulativeDose;
			gravis::d3Vector centre;
			int w0, h0, d0;
			std::vector<int> frameSequence;
			std::string name;
		
			
		double getFrameDose() const;

		// - Dissolve OpticsData?
		
		// gravis::d2Vector project(gravis::d3Vector pos, int frame);
		// Image<float> drawCtf(gravis::d3Vector pos, int size);
		// Image<float> drawDoseWeights(int size);

		BufferedImage<float> computeDoseWeight(int boxSize, double binning) const;
		BufferedImage<float> computeNoiseWeight(int boxSize, double binning, double overlap = 2.0) const;
		
		CTF getCtf(int frame, gravis::d3Vector position, double zOffset = 0.0) const;

		Tomogram extractSubstack(gravis::d3Vector position, int width, int height) const;
		Tomogram FourierCrop(double factor, int num_threads, bool downsampleData = true) const;
};


#endif

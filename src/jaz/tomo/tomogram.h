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
			std::vector<gravis::d4Matrix> proj;
			
			std::vector<CTF> centralCTFs;
			std::vector<double> cumulativeDose;
			gravis::d3Vector centre;
			int w0, h0, d0;
			std::vector<int> frameSequence;
		
			
		double getFrameDose() const;
		
		// - Get rid of boxSize ctor parameter -> losing doseWeights
		// - Dissolve OpticsData?
		
		// gravis::d2Vector project(gravis::d3Vector pos, int frame);
		// Image<float> drawCtf(gravis::d3Vector pos, int size);
		// Image<float> drawDoseWeights(int size);
		
		BufferedImage<float> computeDoseWeight(int boxSize, double binning) const;
		
		
};


#endif

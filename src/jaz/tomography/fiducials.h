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
		template <typename T>
		std::vector<gravis::d3Vector> detect(
				const Tomogram& tomo, 
				gravis::d3Vector origin, 
				gravis::d3Vector spacing,
				gravis::d3Vector diagonal,
				double avgRad, double exclusionRad, 
				int maxNum, double minCC = 10.0,
				int num_threads = 1, int coarseBin = 16);
};

template <typename T>
std::vector<gravis::d3Vector> Fiducials::detect(
		const Tomogram& tomo, 
		gravis::d3Vector origin, 
		gravis::d3Vector spacing,
		gravis::d3Vector diagonal,
		double avgRad, double exclusionRad, 
		int maxNum, double minCC, 
		int num_threads, int coarseBin)
{
	const int w = tomo.stack.xdim;
	const int h = tomo.stack.ydim;
	
	
}

#endif

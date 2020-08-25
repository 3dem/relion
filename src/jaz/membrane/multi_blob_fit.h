#ifndef MULTI_BLOB_FIT_H
#define MULTI_BLOB_FIT_H

#include <vector>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/tomography/tomogram.h>

class MultiBlobFit : public Optimization
{
	public:
		
		MultiBlobFit(
			const Tomogram& tomogram,
			const std::vector<gravis::d4Vector>& spheres,
			double marginIn, double marginOut);
		
		// clustering is independent of scale
		std::vector<std::vector<int>> cluster(
			const Tomogram& tomogram,
			const std::vector<gravis::d4Vector>& spheres,
			double marginOut);
		
		
		
		
		
		
		
};


#endif

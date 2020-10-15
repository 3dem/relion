#ifndef RECURSIVE_LINE_SEARCH_H
#define RECURSIVE_LINE_SEARCH_H

#include "optimization.h"

class RecursiveLineSearch
{
	public:
		
		static double optimize(
			const Optimization& opt,
			double x0, 
			double x1, 
			int steps, 
			int recursions, 
			bool verbose = false);
};

#endif

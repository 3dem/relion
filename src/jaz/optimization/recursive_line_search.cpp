#include "recursive_line_search.h"
#include <limits>
#include <iostream>


double RecursiveLineSearch::optimize(
		const Optimization &opt, double x0, double x1, int steps, int recursions, bool verbose)
{
	double bestX = x0;
	double minCost = std::numeric_limits<double>::max();
	
	const double step = (x1 - x0) / (steps - 1);
	
	if (verbose)
	{
		std::cout << "testing: " 
				  << x0 << ", " 
				  << x0 + step << ", " 
				  << x0 + 2*step 
				  << ", ... " 
				  << x1 << std::endl;
	}
	
	for (int i = 0; i < steps; i++)
	{
		const double x = x0 + i * step;		
		const double fx = opt.f({x}, 0);
		
		if (fx < minCost)
		{
			minCost = fx;
			bestX = x;
		}
	}
	
	if (verbose)
	{
		std::cout << "   best: " << bestX << " at " << minCost << std::endl;
	}
	
	if (recursions > 0)
	{
		return optimize(opt, bestX - step, bestX + step, steps, recursions - 1, verbose);
	}
	else
	{
		return bestX;
	}
}

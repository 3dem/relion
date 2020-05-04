#ifndef CLOSEST_SPLINE_PT_PROBLEM_H
#define CLOSEST_SPLINE_PT_PROBLEM_H

#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/math/spline.h>


class ClosestSplinePointProblem : public Optimization
{
	public:
		
		ClosestSplinePointProblem(
				const Spline<gravis::d2Vector,double>& spline2D, 
				gravis::d2Vector pos);
		
		const Spline<gravis::d2Vector,double>& spline2D;
		gravis::d2Vector pos;
		
		double f(const std::vector<double>& x, void* tempStorage) const;
		
		double findRecursively(double t0, double t1, int iterations);
};

#endif

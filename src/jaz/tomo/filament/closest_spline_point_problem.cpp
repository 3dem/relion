#include "closest_spline_point_problem.h"

using namespace gravis;

ClosestSplinePointProblem::ClosestSplinePointProblem(
		const Spline<d2Vector,double> &spline2D, d2Vector pos
		)
: spline2D(spline2D),
  pos(pos)
{}

double ClosestSplinePointProblem::f(const std::vector<double> &x, void *tempStorage) const
{
	return (spline2D.getValue(x[0]) - pos).norm2();
}

double ClosestSplinePointProblem::findRecursively(double t0, double t1, int iterations)
{
	const d2Vector s0 = spline2D.getValue(t0);
	const d2Vector v0 = spline2D.getTangent(t0);
	
	double f0 = (pos - s0).dot(v0);
	
	const d2Vector s1 = spline2D.getValue(t1);
	const d2Vector v1 = spline2D.getTangent(t1);
	
	double f1 = (pos - s1).dot(v1);
	
	
	if (f0 > 0.0 && f1 < 0.0)
	{
		for (int i = 0; i < iterations; i++)
		{
			double t = f0 * (t1 - t0) / (f0 - f1);
					
			const d2Vector s = spline2D.getValue(t);
			const d2Vector v = spline2D.getTangent(t);
			const double f = (pos - s).dot(v);
			
			if (f > 0.0)
			{
				f0 = f;
				t0 = t;
			}
			else if (f == 0.0) 
			{
				return t;
			}
			else 
			{
				f1 = f;
				t1 = t;
			}
		}
	}
	else
	{
		std::cerr << "warning: f0, f1 = " << f0 << ", " << f1 << std::endl;
	}
	
	return (t0 + t1) / 2.0;
}

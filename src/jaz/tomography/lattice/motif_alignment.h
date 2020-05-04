#ifndef MOTIF_ALIGNMENT_H
#define MOTIF_ALIGNMENT_H

#include <src/jaz/optimization/optimization.h>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>

class MotifAlignment : public DifferentiableOptimization
{
	public:
		
		MotifAlignment(
			const std::vector<gravis::d3Vector>& motifPts,
			const std::vector<gravis::d3Vector>& targetPts,
			double tolerance);
		
		
			std::vector<gravis::d3Vector> motifPts, targetPts;
			double tolerance;
					
		
		double f(const std::vector<double> &x, void *tempStorage) const;
		
		void grad(const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const;
		
		
		static gravis::d4Matrix getMatrix(const std::vector<double> &x);
		
		static std::vector<double> getParams(gravis::d4Matrix A);
				
};

#endif

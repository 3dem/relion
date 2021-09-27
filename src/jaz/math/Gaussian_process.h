#ifndef GAUSS_PROCESS_H
#define GAUSS_PROCESS_H

#include <src/matrix2d.h>
#include <src/jaz/gravis/t3Vector.h>

class GpKernel;

class GaussianProcess
{
	public:
		
		struct Basis
		{
			std::vector<double> eigenvectors, eigenvalues;
		};
		
		static Matrix2D<double> computeCovariance(
				const std::vector<gravis::d3Vector>& points,
				GpKernel* kernel);
		
		static Basis getBasis(const Matrix2D<double>& C, int maxDims, double eps = 1e-10);
};

class GpKernel
{
	public:
			
		virtual double covariance(gravis::d3Vector a, gravis::d3Vector b) = 0;
		
};

class SquareExponentialKernel : public GpKernel
{
	public:
		
		SquareExponentialKernel(double sig_vel, double sig_div);
		
			double sig_vel, sig_div;
		
		double covariance(gravis::d3Vector a, gravis::d3Vector b);
};

class ExponentialKernel : public GpKernel
{
	public:
		
		ExponentialKernel(double sig_vel, double sig_div);
		
			double sig_vel, sig_div;
		
		double covariance(gravis::d3Vector a, gravis::d3Vector b);
};


#endif

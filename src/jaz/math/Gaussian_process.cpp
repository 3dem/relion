#include "Gaussian_process.h"

#include <src/Eigen/Dense>
#include <src/Eigen/SVD>

using namespace Eigen;


Matrix2D<double> GaussianProcess::computeCovariance(
		const std::vector<gravis::d3Vector> &points, 
		GpKernel *kernel)
{
	const int pc = points.size();
	
	Matrix2D<double> A(pc,pc);

    for (int i = 0; i < pc; i++)
    for (int j = i; j < pc; j++)
    {
		const double k = kernel->covariance(points[i], points[j]);
        A(i,j) = k;
        A(j,i) = k;
    }
	
	return A;
}

GaussianProcess::Basis GaussianProcess::getBasis(
		const Matrix2D<double> &C, int maxDims, double eps)
{
	const int N = C.Xdim();
	
	int dc = (maxDims < 0 || maxDims > N)? N : maxDims;

	MatrixXd A0(N,N);

	for (int c = 0; c < N; c++)
	for (int r = 0; r < N; r++)
	{
		A0(r,c) = C(r,c);
	}

	BDCSVD<MatrixXd> svd(A0, ComputeFullV);

	// remove eigendeformations with too small eigenvalues

	for (int d = 0; d < dc; d++)
	{
		if (svd.singularValues()[d] < eps)
		{
			dc = d;
			break;
		}
	}

	GaussianProcess::Basis out;
	out.eigenvalues.resize(dc);
	out.eigenvectors.resize(N*dc);

	for (int d = 0; d < dc; d++)
	{
		const double l = sqrt(svd.singularValues()[d]);

		for (int p = 0; p < N; p++)
		{
			out.eigenvectors[dc*p + d] = l * svd.matrixV()(p,d);
		}
	}

	for (int d = 0; d < dc; d++)
	{
		out.eigenvalues[d] = svd.singularValues()[d];
	}
	
	return out;
}

SquareExponentialKernel::SquareExponentialKernel(double sig_vel, double sig_div)
:	sig_vel(sig_vel),
	sig_div(sig_div)
{}

double SquareExponentialKernel::covariance(gravis::d3Vector a, gravis::d3Vector b)
{
	double d2 = (a - b).norm2() / (sig_div * sig_div);
	
	return sig_vel * sig_vel * exp(-d2);
}

ExponentialKernel::ExponentialKernel(double sig_vel, double sig_div)
:	sig_vel(sig_vel),
	sig_div(sig_div)
{}

double ExponentialKernel::covariance(gravis::d3Vector a, gravis::d3Vector b)
{
	double d = (a - b).length() / sig_div;
	
	return sig_vel * sig_vel * exp(-d);
}

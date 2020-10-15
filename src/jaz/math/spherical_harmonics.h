#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H

#include <vector>
#include <src/jaz/gravis/t3Vector.h>

class SphericalHarmonics
{
	public:
		
		static double getValue(gravis::d3Vector v, int bands, double* coeffs);
		static double getValue(double phi, double theta, int bands, double* coeffs);
		
		
	private:
		
		static std::vector<std::vector<std::vector<double>>> Legendre_coeffs;
		
		static void initializeCoeffs(int bands);
		static double binomialCoeff(int N, int k);
};
		
#endif

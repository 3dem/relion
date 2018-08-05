#ifndef ZERNIKE_H
#define ZERNIKE_H

#include <vector>

class Zernike
{
	public:
		
		static double Z(int m, int n, double rho, double phi);
		static double R(int m, int n, double rho);
		
	private:
		
		static std::vector<std::vector<std::vector<double>>> R_coeffs;
		
		static double factorial(int k);
		static void prepCoeffs(int n);
};

#endif

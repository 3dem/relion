#ifndef ZERNIKE_H
#define ZERNIKE_H

#include <vector>

class Zernike
{
	public:
		
		static double Z(int m, int n, double rho, double phi);
		static double Z_cart(int m, int n, double x, double y);
		static double R(int m, int n, double rho);
		
		static void evenIndexToMN(int i, int& m, int& n);
		static int numberOfEvenCoeffs(int n_max);
		
		static void oddIndexToMN(int i, int& m, int& n);
		static int numberOfOddCoeffs(int n_max);
		
	private:
		
		static std::vector<std::vector<std::vector<double>>> R_coeffs;
		
		static double factorial(int k);
		static void prepCoeffs(int n);
};

#endif

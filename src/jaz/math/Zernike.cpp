#include "Zernike.h"
#include <src/error.h>
#include <sstream>
#include <cmath>

std::vector<std::vector<std::vector<double>>> Zernike::R_coeffs 
	= std::vector<std::vector<std::vector<double>>>(0);

double Zernike::Z(int m, int n, double rho, double phi)
{
	if (m >= 0)
	{
		return R(m, n, rho) * cos(m * phi);
	}
	else
	{
		return R(-m, n, rho) * sin(-m * phi);
	}
}

double Zernike::Z_cart(int m, int n, double x, double y)
{
	if (x == 0 && y == 0)
	{
		return Z(m, n, sqrt(x*x + y*y), 0.0);
	}
	else
	{
		return Z(m, n, sqrt(x*x + y*y), atan2(y,x));
	}
}

double Zernike::R(int m, int n, double rho)
{
	if (m > n)
	{
		REPORT_ERROR_STR("Zernike::R: illegal argument: m = " << m << ", n = " << n << ".\n");
	}
	
	if ((n - m) % 2 == 1) return 0.0;
	
	if (R_coeffs.size() <= n)
	{
		prepCoeffs(n);
	}
	
	double out = 0.0;
	
	for (int k = 0; k <= (n-m)/2; k++)
	{
		out += R_coeffs[n][m][k] * pow(rho, n - 2*k);
	}
	
	return out;
}

void Zernike::evenIndexToMN(int i, int &m, int &n)
{
	const int k = (int)sqrt((double)i);
	
	m = 2*(i - k*k - k);
	n = 2*k;
}

int Zernike::numberOfEvenCoeffs(int n_max)
{
	const int l = n_max / 2;	
	return l*l + 2*l + 1;
}

void Zernike::oddIndexToMN(int i, int& m, int& n)
{
	const int k = (int)((sqrt(1 + 4 * i) - 1.0) / 2.0);
	const int i0 = k*k + k;
			
	n = 2 * k + 1;
	m = 2 * (i - i0) - n;
}

int Zernike::numberOfOddCoeffs(int n_max)
{
	const int l = (n_max - 1) / 2 + 1;
	return l * l + l;
}

double Zernike::factorial(int k)
{
	// @TODO_C++11: replace by tgamma(k+1) once C++11 becomes available
	
	double out = 1.0;
	
	for (int i = 2; i <= k; i++)
	{
		out *= (double)i;
	}
	
	return out;
}

void Zernike::prepCoeffs(int n)
{
	std::vector<std::vector<std::vector<double>>> newCoeffs(n+1);
	
	for (int nn = 0; nn < R_coeffs.size(); nn++)
	{
		newCoeffs[nn] = R_coeffs[nn];
	}
	
	for (int nn = R_coeffs.size(); nn <= n; nn++)
	{
		newCoeffs[nn] = std::vector<std::vector<double>>(nn+1);
				
		for (int m = 0; m <= nn; m++)
		{			
			if ((nn - m) % 2 == 1) continue;
			
			newCoeffs[nn][m] = std::vector<double>((nn-m)/2 + 1);
			
			for (int k = 0; k <= (nn-m)/2; k++)
			{
				newCoeffs[nn][m][k] = 
					  (1 - 2*(k%2)) * factorial(nn-k) 
					/ (factorial(k) * factorial((nn+m)/2 - k) * factorial((nn-m)/2 - k));				
			}
		}
	}

	R_coeffs = newCoeffs;	
}

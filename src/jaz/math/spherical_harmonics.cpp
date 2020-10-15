#include "spherical_harmonics.h"

using namespace gravis;


std::vector<std::vector<std::vector<double>>> SphericalHarmonics::Legendre_coeffs 
	= std::vector<std::vector<std::vector<double>>>(0);


double SphericalHarmonics::getValue(gravis::d3Vector v, int bands, double *coeffs)
{
	const double phi = atan2(v.y, v.x);
	const double rho = sqrt(v.x * v.x + v.y * v.y);
	const double theta = atan2(v.z, rho);
	
	return getValue(phi, theta, bands, coeffs);
}

double SphericalHarmonics::getValue(double phi, double theta, int bands, double *coeffs)
{
	if (Legendre_coeffs.size() < bands)
	{
		initializeCoeffs(bands);
	}
}

void SphericalHarmonics::initializeCoeffs(int bands)
{
	Legendre_coeffs.resize(bands);
	
	for (int l = 0; l < bands; l++)
	{
		Legendre_coeffs[l] = std::vector<std::vector<double>>(l+1);
		
		for (int m = 0; m <= l; m++)
		{
			Legendre_coeffs[l][m] = std::vector<double>(l+1);
			
			for (int k = 0; k <= l; k++)
			{
				if (m > k)
				{
					Legendre_coeffs[l][m][k] = 0.0;
				}
				else
				{
					double derivFact = 1.0;
					
					for (int j = k - m + 1; j <= k; j++)
					{
						derivFact *= (double) j;
					}
							
					//Legendre_coeffs[l][m][k] = binomialCoeff(l,k) * binomialCoeff(
				}
			}
		}
	}
}

double SphericalHarmonics::binomialCoeff(int N, int k)
{
	double out = 1.0;
	
	for (int i = k+1; i < N-k; i++)
	{
		out *= (double) i;
	}
	
	return out;
}

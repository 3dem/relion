#include "blob.h"
#include <src/spherical-harmonics/SphericalHarmonics.h>

#include <omp.h>

Blob::Blob()
:	center(0.0), 
	outer_radius(0), 
	sphericalHarmonics(0),
	shCoeffs(0), shBands(0)
{}

Blob::Blob(gravis::d3Vector center, int outer_radius)
:	center(center), 
	outer_radius(outer_radius), 
	sphericalHarmonics(0),
	shCoeffs(0), shBands(0)
{}
	
Blob::Blob(const std::vector<double>& params, int outer_radius,
		   au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics)
:	center(params[0], params[1], params[2]),
	outer_radius(outer_radius), 
	sphericalHarmonics(sphericalHarmonics),
	shCoeffs(params.size() - 3),
	shBands(sqrt(params.size() - 3) + 0.5 - 1)
{
	for (int i = 0; i < shCoeffs.size(); i++)
	{
		shCoeffs[i] = params[i+3];
	}
}

std::vector<double> Blob::toVector()
{
	std::vector<double> out(shCoeffs.size() + 3);
	
	for (int i = 0; i < 3; i++)
	{
		out[i] = center[i];
	}
	
	for (int i = 0; i < shCoeffs.size(); i++)
	{
		out[i+3] = shCoeffs[i];
	}
	
	return out;
}

double Blob::getOffset(gravis::d3Vector v)
{
	const int cc = shCoeffs.size();
	
	if (cc < 1) return 0.0;
	
	std::vector<double> basis(cc);
	getBasis(v, &basis[0]);
	
	double out(0.0);
	
	for (int i = 1; i < cc; i++)
	{
		out += shCoeffs[i] * basis[i];
	}
	
	return out;
}

void Blob::getBasis(gravis::d3Vector v, double *dest)
{
	const int cc = shCoeffs.size();
	
	if (cc < 1) return;
	
	v.normalize();
	
	const double phi = atan2(v.y, v.x);
	
	std::vector<double> Y(cc);
	
	#pragma omp critical
	{
		sphericalHarmonics->computeY(shBands, v.z, phi, &Y[0]);
	}
	
	for (int i = 1; i < cc; i++)
	{
		dest[i] = Y[i];
	}
}

std::vector<double> Blob::accelerate(gravis::d3Vector ux, gravis::d3Vector uy, int bins)
{
	std::vector<double> accSH(bins);
	
	for (int i = 0; i < bins; i++)
	{
		const double phi = 2.0 * PI * i / (double)bins;
		const double dx = cos(phi);
		const double dy = sin(phi);
		
		accSH[i] = getOffset(dx * ux + dy * uy);
	}
	
	return accSH;
}

std::vector<double> Blob::accelerateBasis(gravis::d3Vector ux, gravis::d3Vector uy, int bins)
{
	const int cc = shCoeffs.size();	
	
	if (cc < 1) return std::vector<double>(0);
	
	std::vector<double> accSHbasis(bins*cc);
		
	for (int i = 0; i < bins; i++)
	{
		const double phi = 2.0 * PI * i / (double)bins;
		const double dx = cos(phi);
		const double dy = sin(phi);
		
		getBasis(dx * ux + dy * uy, &accSHbasis[i*cc]);
	}
	
	return accSHbasis;
}

double Blob::getOffsetAcc(double dx, double dy, const std::vector<double>& accSH)
{
	if (dx == 0.0 && dy == 0.0) return 0.0;
	
	const int bins = accSH.size();
	
	double phi = atan2(dy,dx);
	if (phi < 0.0) phi += 2.0 * PI;
	
	const double id = bins * phi / (2.0 * PI);
	
	const int i0 = (int)id;
	const int i1 = (i0+1) % bins;
	const double di = id - i0;
	
	return (1.0 - di) * accSH[i0] + di * accSH[i1];
}


double Blob::getBasisAcc(double dx, double dy, int b, const std::vector<double>& accSHbasis)
{
	if (dx == 0.0 && dy == 0.0) return 0.0;
	
	const int cc = shCoeffs.size();
	const int bins = accSHbasis.size() / cc;
	
	double phi = atan2(dy,dx);
	if (phi < 0.0) phi += 2.0 * PI;
	
	const double id = bins * phi / (2.0 * PI);
	
	const int i0 = (int)id;
	const int i1 = (i0+1) % bins;
	const double di = id - i0;
	
	return (1.0 - di) * accSHbasis[i0*cc + b] + di * accSHbasis[i1*cc + b];
	
}

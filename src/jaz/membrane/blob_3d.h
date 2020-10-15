#ifndef BLOB_2D_H
#define BLOB_2D_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/spherical-harmonics/SphericalHarmonics.h>


class Blob3D
{
	public:
		
		Blob3D();
		
		Blob3D(gravis::d3Vector center, int outer_radius);
		
		Blob3D(const std::vector<double>& params, int outer_radius,
			 au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics);
		
		
			gravis::d3Vector center;
			int outer_radius;
			au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics;
			
			std::vector<double> shCoeffs;//, accSH, accSHbasis;
			int shBands;
		

		std::vector<double> radialAverage(
				const RawImage<float>& frame,
				const gravis::d4Matrix& proj,
				const RawImage<float>& weight,
				int radius = -1);

		double radialAverageError(
				const RawImage<float>& frame,
				const gravis::d4Matrix& proj,
				const RawImage<float>& weight,
				const std::vector<double>& radAvg);

		BufferedImage<float> drawError(
				const RawImage<float>& frame,
				const gravis::d4Matrix& proj,
				const RawImage<float>& weight,
				const std::vector<double>& radAvg);

		BufferedImage<float> radialAverageProjection(
				const RawImage<float>& frame,
				const gravis::d4Matrix& proj,
				const std::vector<double>& radAvg);

		void decompose(
				RawImage<float>& frame,
				RawImage<float>& blob,
				const gravis::d4Matrix& proj,
				const RawImage<float>& weight,
				double taper);
		
		
		inline std::vector<double> toVector();
		
		inline double getOffset(gravis::d3Vector v);
		inline void getBasis(gravis::d3Vector v, double* dest);
		
		inline std::vector<double> accelerate(gravis::d3Vector ux, gravis::d3Vector uy, int bins);
		inline std::vector<double> accelerateBasis(gravis::d3Vector ux, gravis::d3Vector uy, int bins);
		
		inline double getOffsetAcc(double dx, double dy, const std::vector<double>& accSH);
		inline double getBasisAcc(double dx, double dy, int b, const std::vector<double>& accSHbasis);

		inline double smoothOrigin(double r, double radius);
		
};


inline std::vector<double> Blob3D::toVector()
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

inline double Blob3D::getOffset(gravis::d3Vector v)
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

inline void Blob3D::getBasis(gravis::d3Vector v, double *dest)
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

inline std::vector<double> Blob3D::accelerate(gravis::d3Vector ux, gravis::d3Vector uy, int bins)
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

inline std::vector<double> Blob3D::accelerateBasis(gravis::d3Vector ux, gravis::d3Vector uy, int bins)
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

inline double Blob3D::getOffsetAcc(double dx, double dy, const std::vector<double>& accSH)
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


inline double Blob3D::getBasisAcc(double dx, double dy, int b, const std::vector<double>& accSHbasis)
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

inline double Blob3D::smoothOrigin(double r, double radius)
{
	if (r < radius / 2)
	{
		return r * r / radius + radius / 4;
	}
	else
	{
		return r;
	}
}

#endif

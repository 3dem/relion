#ifndef BLOB_2D_H
#define BLOB_2D_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/spherical-harmonics/SphericalHarmonics.h>


class Blob2D
{
	public:
		
		Blob2D();
		
		Blob2D(gravis::d2Vector center, int outer_radius);
		
		Blob2D(const std::vector<double>& params, int outer_radius,
			 au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics);
		
		
			gravis::d2Vector center;
			int outer_radius;
			
			std::vector<dComplex> amplitudes;
		

		std::vector<double> radialAverage(
				const RawImage<float>& frame,
				const RawImage<float>& weight,
				int radius = -1);

		double radialAverageError(
				const RawImage<float>& frame,
				const RawImage<float>& weight,
				const std::vector<double>& radAvg);

		BufferedImage<float> drawError(
				const RawImage<float>& frame,
				const RawImage<float>& weight,
				const std::vector<double>& radAvg);

		BufferedImage<float> radialAverageProjection(
				const RawImage<float>& frame,
				const std::vector<double>& radAvg);

		void decompose(
				RawImage<float>& frame,
				RawImage<float>& blob,
				const RawImage<float>& weight,
				double taper);
		
		
		inline std::vector<double> toVector();
		
		inline double getOffset(gravis::d2Vector v);

		inline double smoothOrigin(double r, double radius);
		
};


inline std::vector<double> Blob2D::toVector()
{
	std::vector<double> out(amplitudes.size() + 2);
	
	for (int i = 0; i < 2; i++)
	{
		out[i] = center[i];
	}
	
	for (int i = 0; i < amplitudes.size(); i++)
	{
		out[i+2] = amplitudes[i];
	}
	
	return out;
}

inline double Blob2D::getOffset(gravis::d2Vector v)
{
	const int cc = amplitudes.size();
	
	if (cc < 1) return 0.0;
	
	const double phi = atan2(v.y, v.x);

	double out = 0.0;
	
	for (int i = 0; i < cc; i++)
	{
		out += amplitudes[i].real * cos((i+1) * phi);
		out += amplitudes[i].imag * sin((i+1) * phi);
	}
	
	return out;
}

inline double Blob2D::smoothOrigin(double r, double radius)
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

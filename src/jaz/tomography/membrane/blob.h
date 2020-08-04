#ifndef BLOB_H
#define BLOB_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomography/tomogram.h>

namespace au { namespace edu { namespace anu { namespace qm { namespace ro 
{	class SphericalHarmonics;	}}}}}


class Blob
{
	public:
		
		Blob();		
		
		Blob(gravis::d3Vector center, int outer_radius);
		
		Blob(const std::vector<double>& params, int outer_radius, 
			 au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics);
		
		
			gravis::d3Vector center;
			int outer_radius;
			au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics;
			
			std::vector<double> shCoeffs;//, accSH, accSHbasis;
			int shBands;
		

		std::vector<double> radialAverage(
				const Tomogram& tomogram, int f, int radius = -1,
				const RawImage<float>* mask = 0);

		double radialAverageError(
				const Tomogram& tomogram, int f, const std::vector<double>& radAvg,
				const RawImage<float>* mask = 0);

		std::vector<double> radialAverageErrorGrad(
				const Tomogram& tomogram, int f, const std::vector<double>& radAvg);

		BufferedImage<float> radialAverageProjection(
				const Tomogram& tomogram, int f, const std::vector<double>& radAvg);

		// change to pure image stack
		void subtract(Tomogram& tomogram, int f, double taper);
		
		std::vector<double> toVector();
		
		double getOffset(gravis::d3Vector v);
		void getBasis(gravis::d3Vector v, double* dest);
		
		std::vector<double> accelerate(gravis::d3Vector ux, gravis::d3Vector uy, int bins);
		std::vector<double> accelerateBasis(gravis::d3Vector ux, gravis::d3Vector uy, int bins);
		
		double getOffsetAcc(double dx, double dy, const std::vector<double>& accSH);
		double getBasisAcc(double dx, double dy, int b, const std::vector<double>& accSHbasis);
		
};


#endif

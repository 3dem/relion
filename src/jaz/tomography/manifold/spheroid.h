#ifndef SPHEROID_H
#define SPHEROID_H

#include "manifold.h"
#include <src/spherical-harmonics/SphericalHarmonics.h>


class Spheroid : public Manifold
{
	public:
		
		
		Spheroid(const std::vector<double>& parameters, int index);
		
		
			gravis::d3Vector position;
			std::vector<double> SH_coefficients;
			
		
		std::vector<RigidAlignment> sampleParticles(
				double spacing,
				double depth,
				double minTilt,
				double maxTilt,
				bool sample_present_wedge,
				bool sample_missing_wedge,
				const std::vector<gravis::d4Matrix>& projections) const;

		double getRadius(au::edu::anu::qm::ro::SphericalHarmonics& SH,
				  double phi, double theta) const;

		gravis::d3Vector getSurfacePoint(au::edu::anu::qm::ro::SphericalHarmonics& SH,
				  double phi, double theta) const;
		
		std::vector<double> getParameters() const;

		static std::string getTypeName();
};

#endif

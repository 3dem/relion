#ifndef SPHERE_H
#define SPHERE_H

#include "manifold.h"
#include <src/spherical-harmonics/SphericalHarmonics.h>


class Sphere : public Manifold
{
	public:


		Sphere(const std::vector<double>& parameters, int index);


			gravis::d3Vector position;
			double radius;


		std::vector<RigidAlignment> sampleParticles(
				double spacing,
				double depth,
				double minTilt,
				double maxTilt,
				bool sample_present_wedge,
				bool sample_missing_wedge,
				const std::vector<gravis::d4Matrix>& projections) const;

		gravis::d3Vector getSurfacePoint(
				double phi, double theta) const;

		std::vector<double> getParameters() const;

		static std::string getTypeName();
};

#endif

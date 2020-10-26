#include "spheroid.h"
#include <src/euler.h>
#include <src/spherical-harmonics/SphericalHarmonics.h>
#include <src/jaz/tomography/tilt_geometry.h>
#include <src/jaz/math/spherical_coordinates.h>

using namespace gravis;
using namespace au::edu::anu::qm::ro;


Spheroid::Spheroid(const std::vector<double> &parameters, int index)
: Manifold(index, getTypeName())
{
	for (int i = 0; i < 3; i++)
	{
		position[i] = parameters[i];
	}
	
	SH_coefficients.resize(parameters.size() - 3);
	
	for (int i = 3; i < parameters.size(); i++)
	{
		SH_coefficients[i-3] = parameters[i];
	}
}

std::vector<RigidAlignment> Spheroid::sampleParticles(
		double spacing,
		double depth,
		double minTilt,
		double maxTilt,
		bool sample_present_wedge,
		bool sample_missing_wedge,
		const std::vector<d4Matrix>& projections) const
{
	const double sphere_radius = SH_coefficients[0] / (2 * sqrt(PI));
	const double wedge = (maxTilt - minTilt) * sphere_radius;
	const int latitudes = std::round(wedge / spacing);
	const int num_SH_coeffs = SH_coefficients.size();
	const int max_L = sqrt(num_SH_coeffs) - 1;
	const double angular_spacing = spacing / sphere_radius;

	SphericalHarmonics SH(max_L);

	std::vector<RigidAlignment> out;
	out.reserve((latitudes + 1) * (std::round(2 * PI * sphere_radius / spacing)));

	d3Matrix world_to_tilt = TiltGeometry::worldToTiltSpace(projections);

	d3Matrix tilt_to_world = world_to_tilt;
	tilt_to_world.invert();

	std::pair<d3Vector,d3Vector> present_wedge = TiltGeometry::findPresentWedge(
			projections);

	for (int lat = 0; lat <= latitudes; lat++)
	{
		const double tilt0 = minTilt + lat * (maxTilt - minTilt) / latitudes;

		const double perimeter = 2 * PI * sphere_radius * cos(tilt0);
		const int longitudes = std::round(perimeter / spacing);

		for (int lon = 0; lon < longitudes; lon++)
		{
			const double phi0 = 2 * PI * lon / (double) longitudes;

			SphericalCoordinates tilt_space_spherical(phi0, tilt0);
			d3Vector world_space_Cartesian = tilt_to_world * tilt_space_spherical.toCartesian();
			SphericalCoordinates world_space_spherical(world_space_Cartesian);

			const double phi = world_space_spherical.phi;
			const double tilt = world_space_spherical.theta;
			const double d = angular_spacing;

			const d3Vector surface = getSurfacePoint(SH, phi, tilt);

			const d3Vector surface_phi_0 = getSurfacePoint(SH, phi - d, tilt);
			const d3Vector surface_phi_1 = getSurfacePoint(SH, phi + d, tilt);

			const d3Vector surface_tilt_0 = getSurfacePoint(SH, phi, tilt - d);
			const d3Vector surface_tilt_1 = getSurfacePoint(SH, phi, tilt + d);

			d3Vector normal = (surface_phi_1 - surface_phi_0).cross(
						surface_tilt_1 - surface_tilt_0).normalize();

			if (normal.dot(surface - position) < 0.0)
			{
				normal = -normal;
			}

			double wedge_0 = present_wedge.first.dot(normal);
			double wedge_1 = present_wedge.second.dot(normal);

			const bool in_present_wedge = wedge_0 * wedge_1 < 0.0;

			if (sample_present_wedge && in_present_wedge
				|| sample_missing_wedge && !in_present_wedge)
			{
				out.push_back(RigidAlignment(
						surface - depth * normal,
						0,
						-RAD2DEG(acos(normal.z)),
						-RAD2DEG(atan2(normal.y, normal.x))));
			}
		}
	}

	return out;
}

double Spheroid::getRadius(SphericalHarmonics& SH, double phi, double theta) const
{
	const int num_SH_coeffs = SH_coefficients.size();
	const int max_L = sqrt(num_SH_coeffs) - 1;

	std::vector<double> Y(num_SH_coeffs);
	SH.computeY(max_L, sin(theta), phi, &Y[0]);

	double radius = 0.0;

	for (int i = 0; i < num_SH_coeffs; i++)
	{
		radius += SH_coefficients[i] * Y[i];
	}

	return radius;
}

d3Vector Spheroid::getSurfacePoint(SphericalHarmonics &SH, double phi, double theta) const
{
	const double r = getRadius(SH, phi, theta);

	const d3Vector direction(
			cos(phi) * cos(theta),
			sin(phi) * cos(theta),
			sin(theta));

	return position + r * direction;
}

std::vector<double> Spheroid::getParameters() const
{
	std::vector<double> parameters(SH_coefficients.size() + 3);

	for (int i = 0; i < 3; i++)
	{
		parameters[i] = position[i];
	}
	
	for (int i = 3; i < parameters.size(); i++)
	{
		parameters[i] = SH_coefficients[i-3];
	}
	
	return parameters;
}

std::string Spheroid::getTypeName()
{
	return "spheroid";
}

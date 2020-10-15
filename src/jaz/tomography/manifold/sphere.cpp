#include "sphere.h"
#include <src/euler.h>
#include <src/spherical-harmonics/SphericalHarmonics.h>
#include <src/jaz/tomography/tilt_geometry.h>
#include <src/jaz/math/spherical_coordinates.h>

using namespace gravis;


Sphere::Sphere(const std::vector<double> &parameters, int index)
: Manifold(index, getTypeName())
{
	for (int i = 0; i < 3; i++)
	{
		position[i] = parameters[i];
	}

	radius = parameters[3];
}

std::vector<RigidAlignment> Sphere::sampleParticles(
		double spacing,
		double depth,
		double minTilt,
		double maxTilt,
		bool sample_present_wedge,
		bool sample_missing_wedge,
		const std::vector<d4Matrix>& projections) const
{
	const double wedge = (maxTilt - minTilt) * radius;
	const int latitudes = std::round(wedge / spacing);
	const double angular_spacing = spacing / radius;


	std::vector<RigidAlignment> out;
	out.reserve((latitudes + 1) * (std::round(2 * PI * radius / spacing)));

	d3Matrix world_to_tilt = TiltGeometry::worldToTiltSpace(projections);

	d3Matrix tilt_to_world = world_to_tilt;
	tilt_to_world.invert();

	std::pair<d3Vector,d3Vector> present_wedge = TiltGeometry::findPresentWedge(
			projections);

	for (int lat = 0; lat <= latitudes; lat++)
	{
		const double tilt0 = minTilt + lat * (maxTilt - minTilt) / latitudes;

		const double perimeter = 2 * PI * radius * cos(tilt0);
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

			const d3Vector surface = getSurfacePoint(phi, tilt);

			const d3Vector surface_phi_0 = getSurfacePoint(phi - d, tilt);
			const d3Vector surface_phi_1 = getSurfacePoint(phi + d, tilt);

			const d3Vector surface_tilt_0 = getSurfacePoint(phi, tilt - d);
			const d3Vector surface_tilt_1 = getSurfacePoint(phi, tilt + d);

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

d3Vector Sphere::getSurfacePoint(double phi, double theta) const
{
	const d3Vector direction(
			cos(phi) * cos(theta),
			sin(phi) * cos(theta),
			sin(theta));

	return position + radius * direction;
}

std::vector<double> Sphere::getParameters() const
{
	return std::vector<double> {
		position.x,
		position.y,
		position.z,
		radius};
}

std::string Sphere::getTypeName()
{
	return "sphere";
}

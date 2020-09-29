#include "spherical_coordinates.h"

using namespace gravis;

SphericalCoordinates::SphericalCoordinates()
{
}

SphericalCoordinates::SphericalCoordinates(double phi, double theta)
:	phi(phi), theta(theta)
{
}

SphericalCoordinates::SphericalCoordinates(const d3Vector& cartesian)
{
	const double len = cartesian.length();

	if (len == 0.0)
	{
		phi = 0.0;
		theta = 0.0;
	}
	else
	{
		const d3Vector v = cartesian / len;
		const double r = sqrt(v.x * v.x + v.y * v.y);

		theta = atan2(v.z, r);
		phi = atan2(v.y, v.x);
	}
}

d3Vector SphericalCoordinates::toCartesian() const
{
	return d3Vector(
		cos(phi) * cos(theta),
		sin(phi) * cos(theta),
		sin(theta));
}

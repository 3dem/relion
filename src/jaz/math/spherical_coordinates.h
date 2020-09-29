#ifndef SPHERICAL_COORDINATES_H
#define SPHERICAL_COORDINATES_H

#include <src/jaz/gravis/t3Vector.h>


class SphericalCoordinates
{
	public:

		SphericalCoordinates();
		SphericalCoordinates(double phi, double theta);
		SphericalCoordinates(const gravis::d3Vector& cartesian);

			double phi, theta;

		gravis::d3Vector toCartesian() const;
};

#endif

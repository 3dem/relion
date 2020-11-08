#ifndef TILT_GEOMETRY_H
#define TILT_GEOMETRY_H

#include <src/jaz/gravis/t4Matrix.h>
#include <vector>


class TiltGeometry
{
	public:

		static gravis::d3Vector estimateTiltAxis(
				const std::vector<gravis::d4Matrix>& projections);

		/* returns a matrix with the following properties:
		  - z axis (third row) points parallel to the tilt axis
		  - y axis (second row) points forward (average view-z direction)
		  - all axes orthogonal */
		static gravis::d3Matrix worldToTiltSpace(
			const std::vector<gravis::d4Matrix>& projections);

		static std::pair<gravis::d3Vector,gravis::d3Vector> findPresentWedge(
			const std::vector<gravis::d4Matrix>& projections);
};

#endif

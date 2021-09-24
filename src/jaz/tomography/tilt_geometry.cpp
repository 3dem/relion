#include <stdexcept>
#include <limits>

#include "tilt_geometry.h"
#include <src/jaz/util/index_sort.h>

using namespace gravis;

extern "C"
{
		#include <src/jaz/single_particle/d3x3/dsyevh3.h>
		#include <src/jaz/single_particle/d3x3/dsyevc3.h>
}

d3Vector TiltGeometry::estimateTiltAxis(const std::vector<d4Matrix> &projections)
{
	const int fc = projections.size();

	double A[3][3];

	for (int r = 0; r < 3; r++)
	for (int c = 0; c < 3; c++)
	{
		A[r][c] = 0.0;
	}

	for (int f = 0; f < fc; f++)
	{
		const d3Matrix Af = d3Matrix::extract(projections[f]);
		const d3Vector z(Af(2,0), Af(2,1), Af(2,2));

		for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
		{
			A[r][c] += z[r]*z[c];
		}
	}

	double Q[3][3];
	std::vector<double> w(3);

	dsyevh3(A, Q, &w[0]);

	std::vector<int> inds = IndexSort<double>::sortIndices(w);

	return d3Vector(
		Q[0][inds[0]],
		Q[1][inds[0]],
		Q[2][inds[0]]);
}

d3Matrix TiltGeometry::worldToTiltSpace(
		const std::vector<d4Matrix>& projections)
{
	const int fc = projections.size();

	const d3Vector a = estimateTiltAxis(projections);

	d3Vector mean_z(0.0, 0.0, 0.0);

	for (int f = 0; f < fc; f++)
	{
		const d4Matrix& A = projections[f];
		const d3Vector z = d3Vector(A(2,0), A(2,1), A(2,2));

		mean_z += z;
	}

	mean_z.normalize();

	d3Vector p = mean_z.cross(a).normalize();
	const d3Vector q = a.cross(p).normalize();

	return d3Matrix(
		p.x, p.y, p.z,
		q.x, q.y, q.z,
		a.x, a.y, a.z );
}


std::pair<d3Vector, d3Vector> TiltGeometry::findPresentWedge(
		const std::vector<d4Matrix> &projections)
{
	const int fc = projections.size();

	std::vector<d3Vector> all_z(fc);
	d3Vector mean_z(0.0, 0.0, 0.0);


	double B[3][3];

	for (int r = 0; r < 3; r++)
	for (int c = 0; c < 3; c++)
	{
		B[r][c] = 0.0;
	}

	for (int f = 0; f < fc; f++)
	{
		const d4Matrix& A = projections[f];
		const d3Vector z = d3Vector(A(2,0), A(2,1), A(2,2));

		all_z[f] = z;
		mean_z += z;

		for (int r = 0; r < 3; r++)
		for (int c = 0; c < 3; c++)
		{
			B[r][c] += z[r] * z[c];
		}
	}

	double Q[3][3];
	std::vector<double> w(3);
	dsyevh3(B, Q, &w[0]);
	const int d_min = IndexSort<double>::sortIndices(w)[0];

	d3Vector tilt_axis(Q[0][d_min], Q[1][d_min], Q[2][d_min]);

	double min_cross = std::numeric_limits<double>::max();
	double max_cross = -std::numeric_limits<double>::max();
	int max_cross_f = fc/2;
	int min_cross_f = fc/2;

	for (int f = 0; f < fc; f++)
	{
		const double cross = mean_z.cross(all_z[f]).dot(tilt_axis);

		if (cross > max_cross)
		{
			max_cross = cross;
			max_cross_f = f;
		}

		if (cross < min_cross)
		{
			min_cross = cross;
			min_cross_f = f;
		}
	}

	return std::make_pair(all_z[min_cross_f], all_z[max_cross_f]);
}

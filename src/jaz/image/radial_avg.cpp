#include "radial_avg.h"

double RadialAvg::get1DIndex(
		double x, double y, double z,
		int w, int h, int d)
{
	const double xx = x;
	const double yy = (y >= h/2)? y - h: y;
	const double zz = (z >= d/2)? z - d: z;

	return sqrt(xx * xx + yy * yy + zz * zz);
}

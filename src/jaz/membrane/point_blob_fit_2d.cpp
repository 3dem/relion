#include "point_blob_fit_2d.h"

using namespace gravis;



PointBlobFit2D::PointBlobFit2D(
		const std::vector<d2Vector> &allPoints, 
		d2Vector origin, 
		double radius, 
		double tolerance)
:	
  origin(origin),
  radius(radius),
  tolerance(tolerance)
{
	points.reserve(allPoints.size());
	
	for (int i = 0; i < allPoints.size(); i++)
	{
		double distance = (allPoints[i] - origin).length();
		
		if (std::abs(distance - radius) > 3 * tolerance)
		{
			points.push_back(allPoints[i]);
		}
	}
}

double PointBlobFit2D::f(const std::vector<double> &x, void *tempStorage) const
{
	Blob2D blob(x, radius);

	double out = 0.0;
	
	for (int i = 0; i < points.size(); i++)
	{
		const d2Vector d = points[i] - blob.center;
		const double r = blob.getOffset(d);
		const double dr = (r - radius) / tolerance;
		
		out -= exp(-dr*dr / 2);
	}
	
	return out;
}

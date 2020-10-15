#include "point_blob_fit_2d.h"
#include <src/jaz/util/drawing.h>

using namespace gravis;



PointBlobFit2D::PointBlobFit2D(
		const std::vector<d2Vector> &allPoints, 
		d2Vector origin, 
		double initial_radius, 
        double tolerance, 
        double tethering)
:	
  origin(origin),
  initial_radius(initial_radius),
  tolerance(tolerance),
  tethering(tethering)
{
	points.reserve(allPoints.size());
	
	for (int i = 0; i < allPoints.size(); i++)
	{
		double distance = (allPoints[i] - origin).length();
		
		if (distance < 2 * initial_radius)
		{
			points.push_back(allPoints[i]);
		}
	}
}

double PointBlobFit2D::f(const std::vector<double> &x, void *tempStorage) const
{
	std::vector<double> x_blob(x.size()-1);
	const double radius = x[0];
	
	for (int i = 1; i < x.size(); i++)
	{
		x_blob[i-1] = x[i];
	}
	        
	Blob2D blob(x_blob, radius);

	double out = 0.0;
	
	for (int i = 0; i < points.size(); i++)
	{
		const d2Vector d = points[i] - blob.center;
		const double r = d.length() + blob.getOffset(d);
		const double dr = (r - radius) / tolerance;
		
		out -= exp(-dr*dr / 2);
	}
	
	out += tethering * (blob.center - origin).norm2();
	
	return out;
}

BufferedImage<float> PointBlobFit2D::visualise(const std::vector<double>& x, int w, int h)
{
	std::vector<double> x_blob(x.size()-1);
	const double radius = x[0];
	
	for (int i = 1; i < x.size(); i++)
	{
		x_blob[i-1] = x[i];
	}
	        
	Blob2D blob(x_blob, radius);
	
	BufferedImage<float> out(w,h);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const d2Vector d = d2Vector(x,y) - blob.center;
		const double r = d.length() + blob.getOffset(d);
		const double dr = (r - radius) / tolerance;
		
		out(x,y) = exp(-dr*dr / 2);
	}
	
	//Drawing::drawCrosses(points, 1.5f, 5, out);
	
	return out;
}

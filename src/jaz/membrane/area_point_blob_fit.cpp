#include "area_point_blob_fit.h"
#include <src/jaz/util/drawing.h>
#include <src/jaz/image/interpolation.h>

using namespace gravis;



AreaPointBlobFit::AreaPointBlobFit(
        const RawImage<float>& pointDensity,
		d2Vector origin, 
		double initial_radius, 
        double tolerance, 
        double tethering,
        double aspectCost)
:	
  origin(origin),
  initial_radius(initial_radius),
  tolerance(tolerance),
  tethering(tethering),
  pointDensity(pointDensity),
  aspectCost(aspectCost)
{
	radialSamples = (int) std::round(2 * PI * initial_radius);
}

double AreaPointBlobFit::f(const std::vector<double> &x, void *tempStorage) const
{
	std::vector<double> x_blob(x.size()-1);
	const double radius = x[0];
	
	for (int i = 1; i < x.size(); i++)
	{
		x_blob[i-1] = x[i];
	}
	        
	Blob2D blob(x_blob, radius);

	double out = 0.0;
	
	for (int i = 0; i < radialSamples; i++)
	{
		const double phi = 2 * PI * i / (double)radialSamples;
		const double offset = blob.getOffset(phi);
		const double r = radius + offset;
		
		const d2Vector p = blob.center + r * d2Vector(cos(phi), sin(phi));
		
		const double rho = Interpolation::linearXY_clip(pointDensity, p.x, p.y);
		
		out += 1.0 - rho;
	}
	
	out /= radialSamples;
	
	out += tethering * (blob.center - origin).norm2();
	
	for (int i = 3; i < x.size(); i++)
	{
		out += aspectCost * x[i] * x[i];
	}
	
	return out;
}

BufferedImage<float> AreaPointBlobFit::visualise(const std::vector<double>& x, int w, int h)
{
	std::vector<double> x_blob(x.size()-1);
	const double radius = x[0];
	
	for (int i = 1; i < x.size(); i++)
	{
		x_blob[i-1] = x[i];
	}
	        
	Blob2D blob(x_blob, radius);
	
	BufferedImage<float> out(w,h);
	
	for (int i = 0; i < radialSamples; i++)
	{
		const double phi = 2 * PI * i / (double)radialSamples;
		const double offset = blob.getOffset(phi);
		const double r = radius + offset;
		
		const d2Vector p = blob.center + r * d2Vector(cos(phi), sin(phi));
		
		Drawing::addPoint(p, 1.f, out);
	}
	
	return out;
}

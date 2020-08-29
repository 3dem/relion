#include "area_point_blob_fit.h"
#include <src/jaz/util/drawing.h>
#include <src/jaz/image/interpolation.h>

using namespace gravis;



AreaPointBlobFit::AreaPointBlobFit(
        const RawImage<float>& pointDensity,
		const RawImage<float>& regionalTerm,
		gravis::d2Vector origin,
		double initial_radius, 
		double tolerance,
	    double tethering,
		double aspectCost,
	    double contrastCost)
:	
  origin(origin),
  initial_radius(initial_radius),
  tolerance(tolerance),
  tethering(tethering),
  pointDensity(pointDensity),
  regionalTerm(regionalTerm),
  aspectCost(aspectCost),
  contrastCost(contrastCost)
{
	radialSamples = (int) std::round(2 * PI * initial_radius);
	
	const int r = (int) std::round(initial_radius);
	        
	boxSize = (2*r + 1) * i2Vector(1,1);
	
	int x0 = std::ceil(origin.x - r);
	int y0 = std::ceil(origin.y - r);
	
	if (x0 < 0) 
	{
		x0 = 0;
	}
	else if (x0 + 2*r + 1 >= regionalTerm.xdim) 
	{
		x0 -= x0 + 2*r - regionalTerm.xdim;
	}
	
	if (y0 < 0) 
	{
		y0 = 0;
	}
	else if (y0 + 2*r + 1 >= regionalTerm.ydim) 
	{
		y0 -= y0 + 2*r - regionalTerm.ydim;
	}
	
	boxOrigin = i2Vector(x0, y0);
	
	boxArea = boxSize.x * boxSize.y;
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
	
	if (contrastCost > 0.0)
	{
		double contrastTerm = 0.0;
		
		for (int y = 0; y < boxSize.y; y++)
		for (int x = 0; x < boxSize.x; x++)
		{
			const int xx = boxOrigin.x + x;
			const int yy = boxOrigin.y + y;
			
			const d2Vector pp(xx - blob.center.x, yy - blob.center.y);
			
			const double rad = radius + blob.getOffset(pp);
			const double rad0 = pp.length();
			
			double indicator;
			
			if (rad > rad0 + 1) indicator = 0;
			else if (rad > rad0 - 1) indicator = (rad0 + 1 - rad) / 2.0;
			else indicator = 1.0;
			
			contrastTerm += indicator * regionalTerm(xx,yy);
		}
		
		out += contrastCost * contrastTerm / boxArea;
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

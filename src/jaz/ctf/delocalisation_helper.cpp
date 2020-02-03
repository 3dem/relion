#include "delocalisation_helper.h"
#include <src/jaz/gravis/t2Vector.h>

using namespace gravis;

void DelocalisationHelper::maskOutsideBox(
		const CTF &ctf, double radius, 
		double angpix, int s_orig, 
		MultidimArray<RFLOAT> &fftwCtfImg,
		double offsetx, double offsety)
{
	const int s = fftwCtfImg.ydim;
	const int sh = fftwCtfImg.xdim;
	const double r2 = radius * radius;
	const double as = angpix * s;
	
	const t2Vector<RFLOAT> origin(offsetx, offsety);
			
	for (int y = 0; y < s; y++)
	for (int x = 0; x < sh; x++)
	{
		const double xx = x/as;
		const double yy = y < sh? y/as : (y - s)/as;
		
		t2Vector<RFLOAT> delocCent = RFLOAT(1.0 / (2 * angpix * PI)) * ctf.getGammaGrad(xx,yy);
		
		double out = 0.0;
		double cnt = 0.0;
		
		for (RFLOAT sign = (x == 0? 1.0 : -1.0); sign <= 1.0; sign += 2.0)
		{
			t2Vector<RFLOAT> p = origin + sign * delocCent;
			
			double dx, dy;
			
			if (p.x > 0) dx = s_orig/2 - p.x;
			else dx = s_orig/2 + p.x;
			
			if (p.y > 0) dy = s_orig/2 - p.y;
			else dy = s_orig/2 +p.y;
			
			double fx, fy;
			
			if (dx < -radius) fx = 0.0;
			else if (dx < radius) fx = 1.0 - acos(dx/radius)/PI + dx * sqrt(r2 - dx*dx)/(PI * r2);
			else fx = 1.0;
			
			if (dy < -radius) fy = 0.0;
			else if (dy < radius) fy = 1.0 - acos(dy/radius)/PI + dy * sqrt(r2 - dy*dy)/(PI * r2);
			else fy = 1.0;
					
			const double f = fx * fy;
			
			out += f;
			cnt += 1;
		}
		
		DIRECT_A2D_ELEM(fftwCtfImg, y, x) *= out/cnt;
	}
}

Image<RFLOAT> DelocalisationHelper::plotDelocalisation(
		const CTF &ctf, Image<RFLOAT> &mask, double angpix)
{
	const int s = mask.data.ydim;
	const int sh = mask.data.xdim;
	
	Image<RFLOAT> out(s,s);
	
	const double as = s * angpix;
	
	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		double xx = x < sh? x/as : (x - s)/as;
		double yy = y < sh? y/as : (y - s)/as;
		
		t2Vector<RFLOAT> delocCent = RFLOAT(1.0 / (2 * angpix * PI)) * ctf.getGammaGrad(xx,yy);
		
		if (delocCent.x > -sh && delocCent.x < sh 
		 && delocCent.y > -sh && delocCent.y < sh)
		{
			int rx0 = (int)(delocCent.x + s)%s;
			int ry0 = (int)(delocCent.y + s)%s;
			int rx1 = (int)(delocCent.x + s + 1)%s;
			int ry1 = (int)(delocCent.y + s + 1)%s;
			
			double rxf = delocCent.x - std::floor(delocCent.x);
			double ryf = delocCent.y - std::floor(delocCent.y);
			
			int mx = x < sh? x : s - x;
			int my = x < sh? y : (s - y)%s;
			double mv = mask(my, mx);
			
			out(ry0,rx0) += mv * (1 - rxf) * (1 - ryf);
			out(ry0,rx1) += mv * rxf * (1 - ryf);
			out(ry1,rx0) += mv * (1 - rxf) * ryf;
			out(ry1,rx1) += mv * rxf * ryf;
		}
	}
	
	return out;
}

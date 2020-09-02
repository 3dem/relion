#include "quadric_disc_fit.h"
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/reconstruction.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/filter.h>

using namespace gravis;


QuadricDiscFit::QuadricDiscFit(
		BufferedImage<float>& tiltSeries, 
		const std::vector<d4Matrix> proj, 
		d3Vector surface, 
		d3Vector inside,
		double bin,
		double diam,
		double thickness,
		int num_threads)
{
	const int fc = proj.size();
		
	const d3Vector eh = (surface - inside).normalize();      // north
	const d3Vector ep = (std::abs(eh.x) > std::abs(eh.y))?   // left
				eh.cross(d3Vector(0,1,0)).normalize() :
				eh.cross(d3Vector(1,0,0)).normalize();
	const d3Vector eq = ep.cross(eh).normalize();            // back
	
	const int s = (int) diam;
	const double r = s/2.0;
	
	d3Vector shift = surface - r * bin * (ep + eq + eh);
	d3Vector zero = d3Vector(0,0,0);
	
	localToTomo = d4Matrix
		(	bin * ep.x, bin * eq.x, bin * eh.x, shift.x,
			bin * ep.y, bin * eq.y, bin * eh.y, shift.y,
			bin * ep.z, bin * eq.z, bin * eh.z, shift.z,
			 0.0,  0.0,  0.0,   1.0  );
	
	std::vector<d4Matrix> projAct(fc);
	
	for (int f = 0; f < fc; f++)
	{
		projAct[f] = proj[f] / bin;
		projAct[f](3,3) = 1.0;
		
		projAct[f] = projAct[f] * localToTomo;
	}	
	
	BufferedImage<float> stackAct = Resampling::downsampleFiltStack_2D_full(tiltSeries, bin, num_threads);
			
	for (int f = 0; f < fc; f++)
	{
		RawImage<float> imgf0 = stackAct.getSliceRef(f);
		
		BufferedImage<float> imgLP = ImageFilter::Gauss2D(imgf0, 0, 2.0 * thickness, true);
		
		BufferedImage<float> imgf = imgf0 - imgLP;
		imgf = ImageFilter::Gauss2D(imgf, 0, thickness, true);
		
		stackAct.getSliceRef(f).copyFrom(imgf);
	}
	
	localTomo = BufferedImage<float>(s,s,s);
	
	BufferedImage<float> data(s,s,s), psf(s,s,s), coverage(s,s,s);
	data.fill(0.f);
	psf.fill(0.f);
	coverage.fill(0.f);
	
	const double taperDist = 20.0;
	
	RealSpaceBackprojection::backproject(
		stackAct, projAct, data, num_threads, zero, 1, RealSpaceBackprojection::Linear, taperDist);
		
	RealSpaceBackprojection::backprojectPsf(stackAct, projAct, psf, num_threads, zero, 1);
		
	RealSpaceBackprojection::backprojectCoverage(
		stackAct, projAct, coverage, num_threads, zero, 1, RealSpaceBackprojection::Linear, taperDist);
		
	float wghMean = Normalization::computeWeightedMeanFromWeighted(data, coverage);
	
	data -= wghMean * coverage;
		
	Reconstruction::correct3D_RS(data, psf, localTomo, 1, num_threads);	
	
	float avg = Normalization::computeMean(localTomo);
	float var = Normalization::computeVariance(localTomo, avg);
	
	localTomo /= sqrt(var);
}

double QuadricDiscFit::f(const std::vector<double> &x, void *tempStorage) const
{
	const double Qxx = x[0];
	const double Qxy = x[1];
	const double Qxz = x[2];
	const double Qxw = x[3];
	
	const double Qyy = x[4];
	const double Qyz = x[5];
	const double Qyw = x[6];
	
	const double Qzz = x[7];
	const double Qzw = x[8];
	
	const double Qww = x[9];	
	
	const int s = localTomo.xdim;
	const double m = s/2;
	
	double cost = 0.0;
	
	const double initial_stdev = 2.0;
		
	for (int yi = 0; yi < s; yi++)
	for (int xi = 0; xi < s; xi++)
	{		
		const double x = xi - m;
		const double y = yi - m;
		
		const double a = Qzz;
		const double b = 2*(Qxz*x + Qyz*y + Qzw);
		
		const double c = Qxx*x*x + 2*Qxy*x*y + 2*Qxw*x
		                         +   Qyy*y*y + 2*Qyw*y
		                                     +   Qww  ;		
		
		const double discr = b*b - 4*a*c;
		
		if (discr > 0.0)
		{
			const double dh = b > 0.0? sqrt(discr) : -sqrt(discr);
			const double h = m + (-b + dh) / (2*a);
			
			if (h >= 0 && h < s-1)
			{				
				const int hi = (int) h;
				const double hf = h - hi;
				
				const double v0 = (double) localTomo(xi,yi,hi);
				const double v1 = (double) localTomo(xi,yi,hi+1);
				
				cost += (1-hf)*v0 + hf*v1;
			}
			// else: see below
		}
		// else: implicit zero cost (avg. value)
	}
	
	cost /= s * s;
	
	{
		const double a = Qzz;
		const double b = 2 * Qzw;		
		const double c = Qww;		
		
		const double discr = b*b - 4*a*c;
		
		if (discr > 0.0)
		{
			const double dh = b > 0.0? sqrt(discr) : -sqrt(discr);
			const double eh = (-b + dh) / (2*a);
		
			cost += eh * eh / (initial_stdev * initial_stdev);
		}
		else
		{
			cost += 1000.0;
		}
	}
	
	return cost;
}

std::vector<double> QuadricDiscFit::getInitial(double dist)
{
	const double rn2 = 1.0 / (dist * dist);
	const double rn1 = 1.0 / dist;
	
	return {rn2, 0.0, 0.0, 0.0, 
	             rn2, 0.0, 0.0, 
	                  rn2, rn1,
	                       0.0 };
}

d4Matrix QuadricDiscFit::getMatrix(const std::vector<double> &x)
{
	return d4Matrix(x[0], x[1], x[2], x[3],
	                x[1], x[4], x[5], x[6],
	                x[2], x[5], x[7], x[8],
	                x[3], x[6], x[8], x[9] );
}

BufferedImage<float> QuadricDiscFit::visualize(const std::vector<double> &x)
{
	const double Qxx = x[0];
	const double Qxy = x[1];
	const double Qxz = x[2];
	const double Qxw = x[3];
	
	const double Qyy = x[4];
	const double Qyz = x[5];
	const double Qyw = x[6];
	
	const double Qzz = x[7];
	const double Qzw = x[8];
	
	const double Qww = x[9];	
	
	const int s = localTomo.xdim;
	const double m = s/2;
	
	
	BufferedImage<float> out(s,s,s);
	out.fill(0.f);
			
	for (int yi = 0; yi < s; yi++)
	for (int xi = 0; xi < s; xi++)
	{		
		const double x = xi - m;
		const double y = yi - m;
		
		const double a = Qzz;
		const double b = 2*(Qxz*x + Qyz*y + Qzw);
		
		const double c = Qxx*x*x + 2*Qxy*x*y + 2*Qxw*x
		                         +   Qyy*y*y + 2*Qyw*y
		                                     +   Qww  ;		
		
		const double discr = b*b - 4*a*c;
		
		if (discr > 0.0)
		{
			const double dh = b > 0.0? sqrt(discr) : -sqrt(discr);
			const double h = m + (-b + dh) / (2*a);
						
			const int hi = (int) h;
			const double hf = h - hi;
						
			if (hi >= 0 && hi < s)
			{
				out(xi,yi,hi) += (1.0 - hf);
			}
			
			if (hi >= -1 && hi < s-1)
			{
				out(xi,yi,hi+1) += hf;
			}
		}
	}	
	
	return out;
}

BufferedImage<float> QuadricDiscFit::visualize2D(
		const std::vector<double> &x,
		const BufferedImage<float>& stack,
		const std::vector<d4Matrix>& proj )
{
	const double Qxx = x[0];
	const double Qxy = x[1];
	const double Qxz = x[2];
	const double Qxw = x[3];
	
	const double Qyy = x[4];
	const double Qyz = x[5];
	const double Qyw = x[6];
	
	const double Qzz = x[7];
	const double Qzw = x[8];
	
	const double Qww = x[9];	
	
	const int s = localTomo.xdim;
	const double m = s/2;
	
	const int w = stack.xdim;
	const int h = stack.ydim;
	const int fc = stack.zdim;
	
	
	d4Matrix quadricToLocal;
	quadricToLocal.loadIdentity();
		
	quadricToLocal(0,3) = m;
	quadricToLocal(1,3) = m;
	quadricToLocal(2,3) = m;
	
	std::vector<d4Matrix> projAct(fc);
	
	for (int f = 0; f < fc; f++)
	{
		projAct[f] = proj[f] * localToTomo * quadricToLocal;
	}
	
	BufferedImage<float> out(stack);
	
	for (int f = 0; f < fc; f++)
	{
		for (int yi = 0; yi < s; yi++)
		for (int xi = 0; xi < s; xi++)
		{		
			const double x = xi - m;
			const double y = yi - m;
			
			const double a = Qzz;
			const double b = 2*(Qxz*x + Qyz*y + Qzw);
			
			const double c = Qxx*x*x + 2*Qxy*x*y + 2*Qxw*x
									 +   Qyy*y*y + 2*Qyw*y
												 +   Qww  ;		
			
			const double discr = b*b - 4*a*c;
			
			if (discr > 0.0)
			{
				const double zi = m + (-b + sqrt(discr)) / (2*a);
				
				const double z = zi - m;
				
				d4Vector j4 = projAct[f] * d4Vector(x,y,z,1);
				
				const int jx = (int)(j4.x);
				const int jy = (int)(j4.y);
				
				if (jx >= 0 && jx < w && jy >= 0 && jy < h)
				{
					out(jx, jy, f) += 1000.;
				}
				
				const double z2 = (-b - sqrt(discr)) / (2*a);
				
				d4Vector j42 = projAct[f] * d4Vector(x,y,z2,1);
				
				const int jx2 = (int)(j42.x);
				const int jy2 = (int)(j42.y);
				
				if (jx2 >= 0 && jx2 < w && jy2 >= 0 && jy2 < h)
				{
					out(jx2, jy2, f) -= 1000.;
				}
			}
		}
	}	
	
	return out;
}

BufferedImage<float> QuadricDiscFit::writePhaselines(
		const std::vector<double> &x, int w, int h, 
		const std::vector<d4Matrix>& proj)
{
	const int s = localTomo.xdim;
	const double m = s/2;
	const int fc = proj.size();
	
	d4Matrix quadricToLocal;
	quadricToLocal.loadIdentity();
		
	quadricToLocal(0,3) = m;
	quadricToLocal(1,3) = m;
	quadricToLocal(2,3) = m;
	
	std::vector<d4Matrix> projAct(fc);
	
	for (int f = 0; f < fc; f++)
	{
		projAct[f] = proj[f] * localToTomo * quadricToLocal;
	}
	
	d4Matrix Q = getMatrix(x);
	
	d4Matrix Qa = Q.adjugate();
	
	BufferedImage<float> out(w,h,fc);
		
	
	for (int f = 0; f < fc; f++)
	{
		d4Matrix P = projAct[f];
						
		d4Matrix Pt = P;
		Pt.transpose();
		
		d4Matrix Ca_4 = P * Qa * Pt;
		
		d3Matrix Ca(
			Ca_4(0,0), Ca_4(0,1), Ca_4(0,3), 
			Ca_4(1,0), Ca_4(1,1), Ca_4(1,3), 
			Ca_4(3,0), Ca_4(3,1), Ca_4(3,3) );
		
		d3Matrix C = Ca.adjugate();
		
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			d3Vector p(x,y,1);
			
			out(x,y,f) = p.dot(C*p);
		}
	}
	
	return out;
}

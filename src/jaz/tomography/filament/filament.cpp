#include "filament.h"
#include "closest_spline_point_problem.h"
#include <src/jaz/optimization/recursive_line_search.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <map>
#include <omp.h>

using namespace gravis;


void Filament::computeArcLength(
		const Spline<d3Vector, double>& spline, 
		int samples,
		std::map<double,double>& arc_to_index_out,
		std::map<double,double>& index_to_arc_out,
		double& length_out)
{
	double currLen = 0.0;
	
	const double totLen = spline.length();
	const double dt_di = totLen / (samples - 1);
	
	for (int i = 0; i < samples; i++)
	{
		double t = i * dt_di;
		double v = spline.getTangent(t).length() * dt_di;
		
		arc_to_index_out[currLen] = t;
		index_to_arc_out[t] = currLen;
		
		if (i < samples-1) currLen += v;
	}
	
	length_out = currLen;
}

double Filament::translateIndex(double a, const std::map<double, double> &table)
{
	std::map<double,double>::const_iterator it1 = table.upper_bound(a);
	
	if (it1 == table.begin())
	{
		return it1->second;
	}
	
	std::map<double,double>::const_iterator it0 = it1;
	it0--;
	
	if (it1 == table.end()) return it0->second;
	else if (it0 == table.end()) return it1->second;
	
	double t0, t1, a0, a1;
	
	a0 = it0->first;
	t0 = it0->second;
	a1 = it1->first;
	t1 = it1->second;
	
	const double f = (a - a0) / (a1 - a0);
	
	return f * t1 + (1 - f) * t0;
}

void Filament::computeSplineCoords(
		int w, int h,
		double binning,
		const Spline<gravis::d3Vector, double>& spline3D, 
		const std::map<double,double>& index_to_arc,
		const std::map<double,double>& arc_to_index,
		double arcLen,
		const gravis::d4Matrix& proj, 
		double maxDistBinned,
		RawImage<float>& index_out, 
		RawImage<float>& signed_dist_out,
		RawImage<float>& mask_out)
{
	const int len = spline3D.coeffs_a.size();
		
	Spline<gravis::d2Vector, double> spline2D = projectSpline(spline3D, proj);
	
	const int taps0 = 20;
	std::vector<double> initialTaps(taps0);
	
	for (int i = 0; i < taps0; i++)
	{
		initialTaps[i] = translateIndex(arcLen * i / (double) (taps0-1), arc_to_index);
	}
	
	const double halfStep = 0.5 * arcLen / (binning * (taps0 - 1));
	const double safeDist = sqrt(maxDistBinned*maxDistBinned + halfStep*halfStep);
	
	index_out.fill(0.f);
	signed_dist_out.fill(0.f);
	mask_out.fill(0.f);
		
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const d2Vector pos0(x,y);
		ClosestSplinePointProblem cspp(spline2D, pos0);
		
		double t0 = 0, minDist2 = std::numeric_limits<double>::max();
		int bestTap = 0;
		
		for (int i = 0; i < taps0; i++)
		{
			const double dist2 = (spline2D.getValue(initialTaps[i]) - pos0).norm2();
			
			if (dist2 < minDist2)
			{
				minDist2 = dist2;
				t0 = initialTaps[i];
				bestTap = i;
			}
		}
		
		const double t0_pre = (bestTap == 0)? t0 : initialTaps[bestTap-1];
		const double t0_post = (bestTap == taps0-1)? t0 : initialTaps[bestTap+1];
		
		//const double t0 = RecursiveLineSearch::optimize(cspp, 0.0, len, 20, 0);
				
		bool is_relevant = t0 >= 0.0 && t0 <= len;
		
		if (is_relevant)
		{
			double dist0 = sqrt(minDist2);
			
			if (dist0 > safeDist + 2.0)
			{
				is_relevant = false;
			}
			else
			{
				double t = RecursiveLineSearch::optimize(cspp, t0_pre, t0_post, 5, 4);
				
				double dist = sqrt(cspp.f({t}, 0));
				
				if (dist > maxDistBinned || t < 0.0 || t >= len)
				{
					is_relevant = false;
				}
				else
				{
					index_out(x,y) = translateIndex(t, index_to_arc);
					
					d2Vector pos  = spline2D.getValue(t);
					d2Vector tang = spline2D.getTangent(t).normalize();
					
					d2Vector q(tang.y, -tang.x);
					
					signed_dist_out(x,y) = q.dot(pos0 - pos) * binning;
					mask_out(x,y) = 1.f;
				}
			}
		}	
	}
}

Spline<gravis::d2Vector,double> Filament::projectSpline(
		const Spline<gravis::d3Vector,double>& spline3D,
		const gravis::d4Matrix& proj)
{
	Spline<gravis::d2Vector, double> spline2D;
	
	const int len = spline3D.coeffs_a.size();
	
	spline2D.coeffs_a.resize(len);
	spline2D.coeffs_b.resize(len);
	spline2D.coeffs_c.resize(len);
	spline2D.coeffs_d.resize(len);
	
	for (int s = 0; s < len; s++)
	{
		gravis::d4Vector a4 = proj * gravis::d4Vector(spline3D.coeffs_a[s], 0.0);
		gravis::d4Vector b4 = proj * gravis::d4Vector(spline3D.coeffs_b[s], 0.0);
		gravis::d4Vector c4 = proj * gravis::d4Vector(spline3D.coeffs_c[s], 0.0);
		gravis::d4Vector d4 = proj * gravis::d4Vector(spline3D.coeffs_d[s], 1.0);
		
		spline2D.coeffs_a[s] = gravis::d2Vector(a4.x, a4.y);
		spline2D.coeffs_b[s] = gravis::d2Vector(b4.x, b4.y);
		spline2D.coeffs_c[s] = gravis::d2Vector(c4.x, c4.y);
		spline2D.coeffs_d[s] = gravis::d2Vector(d4.x, d4.y);
	}
	
	return spline2D;
}




Filament::Filament()
{
	
}

Filament::Filament(const Spline<gravis::d3Vector,double>& spline, double maxRadius)
:	spline(spline),
	maxRadius(maxRadius),
	arcLen(0.0)
{
	computeArcLength(spline, 1000, arc_to_index, index_to_arc, arcLen);
}

FilamentMapping Filament::rasteriseCoordinates(
		int w, int h, 
		const std::vector<d4Matrix>& proj,
		double binning,
		int num_threads) const
{
	FilamentMapping out;
	
	const int fc = proj.size();
	
	out.index = BufferedImage<float>(w,h,fc);
	out.signed_dist = BufferedImage<float>(w,h,fc);
	out.mask = BufferedImage<float>(w,h,fc);
			
	out.mask.fill(1.f);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		RawImage<float> splineIndex_f = out.index.getSliceRef(f);
		RawImage<float> signedDist_f = out.signed_dist.getSliceRef(f);
		RawImage<float> splineMask_f = out.mask.getSliceRef(f);
		
		computeSplineCoords(
			w, h, binning,
			spline, index_to_arc, arc_to_index, arcLen, proj[f], maxRadius / binning,
			splineIndex_f, signedDist_f, splineMask_f);
	}
	
	return out;
}

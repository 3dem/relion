#ifndef BLOB_H
#define BLOB_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomo/tomo_stack.h>

namespace au { namespace edu { namespace anu { namespace qm { namespace ro 
{	class SphericalHarmonics;	}}}}}


class Blob
{
	public:
		
		Blob();		
		
		Blob(gravis::d3Vector center, int outer_radius);
		
		Blob(const std::vector<double>& params, int outer_radius, 
			 au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics);
		
		
			gravis::d3Vector center;			
			int outer_radius;			
			au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics;
			
			std::vector<double> shCoeffs;//, accSH, accSHbasis;
			int shBands;
		
			
		template <typename T>
		std::vector<double> radialAverage(
				const TomoStack<T>& ts, int f, int radius = -1,
				const RawImage<float>* mask = 0);
		
		template <typename T>
		double radialAverageError(
				const TomoStack<T>& ts, int f, const std::vector<double>& radAvg,
				const RawImage<float>* mask = 0);
		
		template <typename T>
		std::vector<double> radialAverageErrorGrad(
				const TomoStack<T>& ts, int f, const std::vector<double>& radAvg);
		
		template <typename T>
		BufferedImage<T> radialAverageProjection(
				const TomoStack<T>& ts, int f, const std::vector<double>& radAvg);
		
		template <typename T>
		void subtract(TomoStack<T>& ts, int f, double taper);
		
		std::vector<double> toVector();
		
		double getOffset(gravis::d3Vector v);
		void getBasis(gravis::d3Vector v, double* dest);
		
		std::vector<double> accelerate(gravis::d3Vector ux, gravis::d3Vector uy, int bins);
		std::vector<double> accelerateBasis(gravis::d3Vector ux, gravis::d3Vector uy, int bins);
		
		double getOffsetAcc(double dx, double dy, const std::vector<double>& accSH);
		double getBasisAcc(double dx, double dy, int b, const std::vector<double>& accSHbasis);
		
};


template <typename T>
std::vector<double> Blob::radialAverage(const TomoStack<T>& ts, int f, int radius,
										const RawImage<float>* mask)
{
	if (radius < 0) radius = outer_radius+1;
	
	std::vector<double> radAvg(radius, 0.0), radCnt(radius, 0.0);
	
	const int w = ts.images[f].xdim;
	const int h = ts.images[f].ydim;
	
	gravis::d4Vector pw(center);	
	gravis::d4Vector pi = ts.worldToImage[f] * pw;
	
	gravis::d3Vector dwdx(ts.worldToImage[f](0,0), ts.worldToImage[f](0,1), ts.worldToImage[f](0,2));
	gravis::d3Vector dwdy(ts.worldToImage[f](1,0), ts.worldToImage[f](1,1), ts.worldToImage[f](1,2));
			
	std::vector<double> accSH = accelerate(dwdx, dwdy, 2.0 * PI * radius);
	
	const int x0 = (int) (pi.x + 0.5);
	const int y0 = (int) (pi.y + 0.5);
	
	
	for (int y = y0 - radius + 1; y <= y0 + radius - 1; y++)
	for (int x = x0 - radius + 1; x <= x0 + radius - 1; x++)
	{
		if (x < 0 || x >= w || y < 0 || y >= h) continue;
		
		const double dx = x - pi.x;
		const double dy = y - pi.y;
		
		const double r = sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH);
		
		if (r >= radius - 1 || r < 0) continue;
		
		const double m = mask == 0? 1.0 : (*mask)(x,y);
		
		const int r0 = (int)r;
		const int r1 = r0 + 1;
		
		const double dr = r - r0;
		
		radAvg[r0] += m * (1 - dr) * ts.images[f](x,y); 
		radCnt[r0] += m * (1 - dr);
				
		radAvg[r1] += m * dr * ts.images[f](x,y);
		radCnt[r1] += m * dr;
	}
	
	for (int i = 0; i < radius; i++)
	{
		if (radCnt[i] > 0.0) radAvg[i] /= radCnt[i];
	}
	
	return radAvg;
}


template <typename T>
double Blob::radialAverageError(
		const TomoStack<T>& ts, int f, const std::vector<double>& radAvg, const RawImage<float>* mask)
{
	double out(0.0), totWgh(0.0);
	
	const int radius = radAvg.size();
	
	const int w = ts.images[f].xdim;
	const int h = ts.images[f].ydim;
	
	gravis::d4Vector pw(center);	
	gravis::d4Vector pi = ts.worldToImage[f] * pw;
	
	gravis::d3Vector dwdx(ts.worldToImage[f](0,0), ts.worldToImage[f](0,1), ts.worldToImage[f](0,2));
	gravis::d3Vector dwdy(ts.worldToImage[f](1,0), ts.worldToImage[f](1,1), ts.worldToImage[f](1,2));
			
	std::vector<double> accSH = accelerate(dwdx, dwdy, 2.0 * PI * radius);
	
	const int x0 = (int) (pi.x + 0.5);
	const int y0 = (int) (pi.y + 0.5);
	
	for (int y = y0 - radius + 1; y <= y0 + radius - 1; y++)
	for (int x = x0 - radius + 1; x <= x0 + radius - 1; x++)
	{
		double obs(0.0), m(1.0);
		
		if (x < 0 || x >= w || y < 0 || y >= h)
		{
			int xa, ya;
			
			if (x < 0) xa = 0;
			else if (x >= w) xa = w-1;
			else xa = x;
			
			if (y < 0) ya = 0;
			else if (y >= h) ya = h-1;
			else ya = y;
			
			obs = ts.images[f](xa,ya);
			m = 1e-4;
		}
		else
		{
			obs = ts.images[f](x,y);
			if (mask != 0) m = (*mask)(x,y);
		}
		
		const double dx = x - pi.x;
		const double dy = y - pi.y;
		
		double r = sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH);
		
		if (r >= radius - 1 - 1e-6) r = radius - 1 - 1e-6;
		if (r < 0) r = 0;
		
		const int r0 = (int)r;
		const int r1 = r0 + 1;
		
		const double dr = r - r0;
		
		const double pred = (1 - dr) * radAvg[r0] + dr * radAvg[r1];
		const double err = pred - obs;
		
		out += m * err * err;
		totWgh += m;
	}
	
	if (totWgh > 0.0) return out / totWgh;
	return out;
}

template <typename T>
std::vector<double> Blob::radialAverageErrorGrad(
		const TomoStack<T>& ts, int f, const std::vector<double>& radAvg)
{
/*  
	C = SUM_p (e_p^2)
	e_p = pred_p - obs_p
	pred_p = (1 - dr) * radAvg[r0] + dr * radAvg[r1];
		
	radAvg[r] = (SUM_p m_p * obs_p) / SUM_p m_p
	m_p = 1 - max(r(p) - r, 1)
	
	r(p) = |p - c| + SUM_i q_i Y_i
	
	
	=> dC/dq = SUM_p ( 2 e_p d(e_p)/dq )
	
	d(e_p)/dq = d(pred_p)/dq
	 =	d(1 - dr)/dq * radAvg[r0] + (1 - dr) * d(radAvg[r0])/dq
		+ d(dr)/dq * radAvg[r1] + dr * d(radAvg[r1])/dq
	 =	- d(dr)/dq * radAvg[r0] + (1 - dr) * d(radAvg[r0])/dq
		+ d(dr)/dq * radAvg[r1] + dr * d(radAvg[r1])/dq
	 =	d(dr)/dq * (radAvg[r1] - radAvg[r0]) 
		+ (1 - dr) * d(radAvg[r0])/dq
		dr * d(radAvg[r1])/dq
		
	d(r)/dc = (c - p)/|c - p| (EXCEPT for c == p!) + SUM_i q_i dY_i/dc
	d(r)/dq = q
	
	d(radAvg[r])/dq = ( d(SUM_p m_p * obs_p)/dq * SUM_p m_p - (SUM_p m_p * obs_p) * d(SUM_p m_p)/dq )
					  / (SUM_p m_p)^2
					  
	d(SUM_p m_p * obs_p)/dq = SUM_p obs_p * d(m_p)/dq
	d(SUM_p m_p)/dq = SUM_p d(m_p)/dq
		
	            { 1 for r - 1 <= r(p) < r
	d(m_p)/dq =	{ -1 for r <= r(p) < r + 1
				{ 0 else
				
*/
	const int radius = radAvg.size();
	
	const int w = ts.images[f].xdim;
	const int h = ts.images[f].ydim;
	
	const int qc = shCoeffs.size() + 3;
	
	
	gravis::d4Vector pw(center);	
	gravis::d4Vector pi = ts.worldToImage[f] * pw;
	
	gravis::d3Vector dwdx(ts.worldToImage[f](0,0), ts.worldToImage[f](0,1), ts.worldToImage[f](0,2));
	gravis::d3Vector dwdy(ts.worldToImage[f](1,0), ts.worldToImage[f](1,1), ts.worldToImage[f](1,2));
			
	std::vector<double> accSH = accelerate(dwdx, dwdy, 2.0 * PI * radius);
	std::vector<double> accSHbasis = accelerateBasis(dwdx, dwdy, 2.0 * PI * radius);
	
	const int x0 = (int) (pi.x + 0.5);
	const int y0 = (int) (pi.y + 0.5);
	
	
	std::vector<double> radAvgGrad(radius*qc, 0.0), radCntGrad(radius*qc, 0.0);
	std::vector<double> radAvgVal(radius, 0.0), radCntVal(radius, 0.0);
	
	for (int y = y0 - radius + 1; y <= y0 + radius - 1; y++)
	for (int x = x0 - radius + 1; x <= x0 + radius - 1; x++)
	{
		if (x < 0 || x >= w || y < 0 || y >= h) continue;
		
		const double dx = x - pi.x;
		const double dy = y - pi.y;
		
		double r = sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH);
		
		if (r < radius - 1 - 1e-6 && r >= 0) 
		{
			const int r0 = (int)r;
			const int r1 = r0 + 1;
			
			const T val = ts.images[f](x,y);
			
			gravis::d2Vector cp(pi.x - x, pi.y - y);
			
			if (cp.length() > 1e-20) cp.normalize();
			
			for (int i = 0; i < 3; i++)
			{
				const double drdq = dwdx[i] * cp.x + dwdy[i] * cp.y;
				
				radAvgGrad[qc*r0 + i] -= drdq * val;
				radAvgGrad[qc*r1 + i] += drdq * val;
				
				radCntGrad[qc*r0 + i] -= drdq;
				radCntGrad[qc*r1 + i] += drdq;
			}
			
			for (int q = 4; q < qc; q++)
			{
				const double drdq = getBasisAcc(dx, dy, q-3, accSHbasis);
				
				radAvgGrad[qc*r0 + q] -= drdq * val;
				radAvgGrad[qc*r1 + q] += drdq * val;
				
				radCntGrad[qc*r0 + q] -= drdq;
				radCntGrad[qc*r1 + q] += drdq;
			}
			
			const double dr = r - r0;
			
			radAvgVal[r0] += (1 - dr) * val; 
			radCntVal[r0] += 1 - dr;
					
			radAvgVal[r1] += dr * val;
			radCntVal[r1] += dr;
		}
	}
	
	for (int r = 0; r < radius; r++)
	for (int q = 0; q < qc; q++)
	{
		const size_t i = qc*r + q;
		
		if (radCntVal[r] > 0.0)
		{
			radAvgGrad[i] = (radAvgGrad[i] * radCntVal[r] - radAvgVal[r] * radCntGrad[i]) 
		                / (radCntVal[r] * radCntVal[r]);
		}
		else
		{
			radAvgGrad[i] = 0.0;
		}
	}
	
	
	std::vector<double> out(shCoeffs.size() + 3, 0.0);
	
	
	for (int y = y0 - radius + 1; y <= y0 + radius - 1; y++)
	for (int x = x0 - radius + 1; x <= x0 + radius - 1; x++)
	{
		if (x < 0 || x >= w || y < 0 || y >= h) continue;
		
		const double dx = x - pi.x;
		const double dy = y - pi.y;
		
		double r = sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH);
		
		if (r < radius - 1 - 1e-6 && r >= 0) 
		{
			const int r0 = (int)r;
			const int r1 = r0 + 1;
			
			const double dr = r - r0;
			
			const double pred = (1 - dr) * radAvg[r0] + dr * radAvg[r1];
			const double obs = ts.images[f](x,y);
			const double err = pred - obs;
			
			/*
				d(e_p)/dq = d(pred_p)/dq
				 =	d(1 - dr)/dq * radAvg[r0] + (1 - dr) * d(radAvg[r0])/dq
					+ d(dr)/dq * radAvg[r1] + dr * d(radAvg[r1])/dq
				 =	- d(dr)/dq * radAvg[r0] + (1 - dr) * d(radAvg[r0])/dq
					+ d(dr)/dq * radAvg[r1] + dr * d(radAvg[r1])/dq
				 =	d(dr)/dq * (radAvg[r1] - radAvg[r0]) 
					+ (1 - dr) * d(radAvg[r0])/dq
					dr * d(radAvg[r1])/dq
			*/
			
			gravis::d2Vector cp(pi.x - x, pi.y - y);
			
			if (cp.length() > 1e-20) cp.normalize();
				
			for (int q = 0; q < qc; q++)
			{
				double drdq = 0.0;
				
				if (q < 3)
				{
					drdq = dwdx[q] * cp.x + dwdy[q] * cp.y;
				}
				else if (q > 3)
				{
					drdq = getBasisAcc(dx, dy, q-3, accSHbasis);
				}
				
				const double dpred_dq = 
						drdq * (radAvg[r1] - radAvg[r0])
						+ (1 - dr) * radAvgGrad[qc*r0 + q] 
						+ dr * radAvgGrad[qc*r1 + q];
						
				//const double dpred_dq = drdq * (radAvg[r1] - radAvg[r0]);
				
				out[q] += 2.0 * err * dpred_dq;
			}
		}
	}
	
	return out;
}

template <typename T>
BufferedImage<T> Blob::radialAverageProjection(
		const TomoStack<T>& ts, int f, const std::vector<double>& radAvg)
{
	BufferedImage<T> out(ts.images[f].xdim, ts.images[f].ydim);
	
	const int radius = radAvg.size();
	
	const int w = ts.images[f].xdim;
	const int h = ts.images[f].ydim;
	
	gravis::d4Vector pw(center);	
	gravis::d4Vector pi = ts.worldToImage[f] * pw;
	
	gravis::d3Vector dwdx(ts.worldToImage[f](0,0), ts.worldToImage[f](0,1), ts.worldToImage[f](0,2));
	gravis::d3Vector dwdy(ts.worldToImage[f](1,0), ts.worldToImage[f](1,1), ts.worldToImage[f](1,2));
	
	std::vector<double> accSH = accelerate(dwdx, dwdy, 2.0 * PI * radius);
	
	const int x0 = (int) (pi.x + 0.5);
	const int y0 = (int) (pi.y + 0.5);
	
	for (int y = y0 - radius + 1; y <= y0 + radius - 1; y++)
	for (int x = x0 - radius + 1; x <= x0 + radius - 1; x++)
	{
		if (x < 0 || x >= w || y < 0 || y >= h) continue;
		
		const double dx = x - pi.x;
		const double dy = y - pi.y;
		
		const double r = sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH);
		
		if (r >= radius - 1 || r < 0)
		{
			out(x,y) = ts.images[f](x,y);
			continue;
		}
		
		const int r0 = (int)r;
		const int r1 = r0 + 1;
		
		const double dr = r - r0;
		
		out(x,y) = (1 - dr) * radAvg[r0] + dr * radAvg[r1];
	}
	
	return out;
}


template <typename T>
void Blob::subtract(
	TomoStack<T>& ts, int f, double taper)
{
	const int radius = outer_radius + taper;
	
	std::vector<double> radAvg = radialAverage(ts, f, radius);	
	
	const int w = ts.images[f].xdim;
	const int h = ts.images[f].ydim;
	
	gravis::d4Vector pw(center);	
	gravis::d4Vector pi = ts.worldToImage[f] * pw;
	
	gravis::d3Vector dwdx(ts.worldToImage[f](0,0), ts.worldToImage[f](0,1), ts.worldToImage[f](0,2));
	gravis::d3Vector dwdy(ts.worldToImage[f](1,0), ts.worldToImage[f](1,1), ts.worldToImage[f](1,2));
			
	std::vector<double> accSH = accelerate(dwdx, dwdy, 2.0 * PI * radius);
	
	const int x0 = (int) (pi.x + 0.5);
	const int y0 = (int) (pi.y + 0.5);
	
	double outside_val(0.0), outside_wgh(0.0);
	
	for (int y = y0 - radius + 1; y <= y0 + radius - 1; y++)
	for (int x = x0 - radius + 1; x <= x0 + radius - 1; x++)
	{
		if (x < 0 || x >= w || y < 0 || y >= h) continue;
		
		const double dx = x - pi.x;
		const double dy = y - pi.y;
		
		const double ru = sqrt(dx*dx + dy*dy);
		const double r = ru + getOffsetAcc(dx, dy, accSH);
		
		if (r >= radius - 1 || r < 0)
		{
			continue;
		}
		
		double wgh;
		
		if (ru < outer_radius) wgh = 1.0;
		else if (ru < outer_radius + taper) wgh = cos(0.5 * PI * (ru - outer_radius) / taper);
		else wgh = 0.0;
		
		outside_val += (1 - wgh) * ts.images[f](x,y);
		outside_wgh += (1 - wgh);
	}
	
	if (outside_wgh > 0.0) outside_val /= outside_wgh;
	
	for (int y = y0 - radius + 1; y <= y0 + radius - 1; y++)
	for (int x = x0 - radius + 1; x <= x0 + radius - 1; x++)
	{
		if (x < 0 || x >= w || y < 0 || y >= h) continue;
		
		const double dx = x - pi.x;
		const double dy = y - pi.y;
		
		const double ru = sqrt(dx*dx + dy*dy);
		const double r = ru + getOffsetAcc(dx, dy, accSH);
		
		if (r >= radius - 1 || r < 0)
		{
			continue;
		}
		
		const int r0 = (int)r;
		const int r1 = r0 + 1;
		
		const double dr = r - r0;
		
		const double pred = (1 - dr) * radAvg[r0] + dr * radAvg[r1] - outside_val;
		
		double wgh;
		
		if (ru < outer_radius) wgh = 1.0;
		else if (ru < outer_radius + taper) wgh = cos(0.5 * PI * (ru - outer_radius) / taper);
		else wgh = 0.0;
		
		ts.images[f](x,y) -= wgh * pred;
	}
}

#endif

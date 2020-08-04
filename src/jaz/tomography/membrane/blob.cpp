#include "blob.h"
#include <src/spherical-harmonics/SphericalHarmonics.h>

#include <omp.h>

Blob::Blob()
:	center(0.0), 
	outer_radius(0), 
	sphericalHarmonics(0),
	shCoeffs(0), shBands(0)
{}

Blob::Blob(gravis::d3Vector center, int outer_radius)
:	center(center), 
	outer_radius(outer_radius), 
	sphericalHarmonics(0),
	shCoeffs(0), shBands(0)
{}
	
Blob::Blob(const std::vector<double>& params, int outer_radius,
		   au::edu::anu::qm::ro::SphericalHarmonics* sphericalHarmonics)
:	center(params[0], params[1], params[2]),
	outer_radius(outer_radius), 
	sphericalHarmonics(sphericalHarmonics),
	shCoeffs(params.size() - 3),
	shBands(sqrt(params.size() - 3) + 0.5 - 1)
{
	for (int i = 0; i < shCoeffs.size(); i++)
	{
		shCoeffs[i] = params[i+3];
	}
}

std::vector<double> Blob::radialAverage(
		const Tomogram& tomogram, int f, int radius,
		const RawImage<float>* mask)
{
	if (radius < 0) radius = outer_radius+1;

	std::vector<double> radAvg(radius, 0.0), radCnt(radius, 0.0);

	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;

	const gravis::d4Matrix& proj = tomogram.projectionMatrices[f];

	gravis::d4Vector pw(center);
	gravis::d4Vector pi = proj * pw;

	gravis::d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	gravis::d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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

		radAvg[r0] += m * (1 - dr) * tomogram.stack(x,y,f);
		radCnt[r0] += m * (1 - dr);

		radAvg[r1] += m * dr * tomogram.stack(x,y,f);
		radCnt[r1] += m * dr;
	}

	for (int i = 0; i < radius; i++)
	{
		if (radCnt[i] > 0.0) radAvg[i] /= radCnt[i];
	}

	return radAvg;
}

double Blob::radialAverageError(
		const Tomogram& tomogram, int f, const std::vector<double>& radAvg, const RawImage<float>* mask)
{
	double out(0.0), totWgh(0.0);

	const int radius = radAvg.size();

	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;


	const gravis::d4Matrix& proj = tomogram.projectionMatrices[f];

	gravis::d4Vector pw(center);
	gravis::d4Vector pi = proj * pw;

	gravis::d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	gravis::d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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

			obs = tomogram.stack(xa,ya,f);
			m = 1e-4;
		}
		else
		{
			obs = tomogram.stack(x,y,f);
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

std::vector<double> Blob::radialAverageErrorGrad(
		const Tomogram& tomogram, int f, const std::vector<double>& radAvg)
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

	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;

	const int qc = shCoeffs.size() + 3;


	const gravis::d4Matrix& proj = tomogram.projectionMatrices[f];

	gravis::d4Vector pw(center);
	gravis::d4Vector pi = proj * pw;

	gravis::d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	gravis::d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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

			const float val = tomogram.stack(x,y,f);

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
			const double obs = tomogram.stack(x,y,f);
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

BufferedImage<float> Blob::radialAverageProjection(
		const Tomogram& tomogram, int f, const std::vector<double>& radAvg)
{
	BufferedImage<float> out(tomogram.stack.xdim, tomogram.stack.ydim);

	const int radius = radAvg.size();

	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;

	const gravis::d4Matrix& proj = tomogram.projectionMatrices[f];

	gravis::d4Vector pw(center);
	gravis::d4Vector pi = proj * pw;

	gravis::d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	gravis::d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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
			out(x,y) = tomogram.stack(x,y,f);
			continue;
		}

		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		out(x,y) = (1 - dr) * radAvg[r0] + dr * radAvg[r1];
	}

	return out;
}

void Blob::subtract(
	Tomogram& tomogram, int f, double taper)
{
	const int radius = outer_radius + taper;

	std::vector<double> radAvg = radialAverage(tomogram, f, radius);

	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;

	const gravis::d4Matrix& proj = tomogram.projectionMatrices[f];

	gravis::d4Vector pw(center);
	gravis::d4Vector pi = proj * pw;

	gravis::d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	gravis::d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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

		outside_val += (1 - wgh) * tomogram.stack(x,y,f);
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

		tomogram.stack(x,y,f) -= wgh * pred;
	}
}


std::vector<double> Blob::toVector()
{
	std::vector<double> out(shCoeffs.size() + 3);
	
	for (int i = 0; i < 3; i++)
	{
		out[i] = center[i];
	}
	
	for (int i = 0; i < shCoeffs.size(); i++)
	{
		out[i+3] = shCoeffs[i];
	}
	
	return out;
}

double Blob::getOffset(gravis::d3Vector v)
{
	const int cc = shCoeffs.size();
	
	if (cc < 1) return 0.0;
	
	std::vector<double> basis(cc);
	getBasis(v, &basis[0]);
	
	double out(0.0);
	
	for (int i = 1; i < cc; i++)
	{
		out += shCoeffs[i] * basis[i];
	}
	
	return out;
}

void Blob::getBasis(gravis::d3Vector v, double *dest)
{
	const int cc = shCoeffs.size();
	
	if (cc < 1) return;
	
	v.normalize();
	
	const double phi = atan2(v.y, v.x);
	
	std::vector<double> Y(cc);
	
	#pragma omp critical
	{
		sphericalHarmonics->computeY(shBands, v.z, phi, &Y[0]);
	}
	
	for (int i = 1; i < cc; i++)
	{
		dest[i] = Y[i];
	}
}

std::vector<double> Blob::accelerate(gravis::d3Vector ux, gravis::d3Vector uy, int bins)
{
	std::vector<double> accSH(bins);
	
	for (int i = 0; i < bins; i++)
	{
		const double phi = 2.0 * PI * i / (double)bins;
		const double dx = cos(phi);
		const double dy = sin(phi);
		
		accSH[i] = getOffset(dx * ux + dy * uy);
	}
	
	return accSH;
}

std::vector<double> Blob::accelerateBasis(gravis::d3Vector ux, gravis::d3Vector uy, int bins)
{
	const int cc = shCoeffs.size();	
	
	if (cc < 1) return std::vector<double>(0);
	
	std::vector<double> accSHbasis(bins*cc);
		
	for (int i = 0; i < bins; i++)
	{
		const double phi = 2.0 * PI * i / (double)bins;
		const double dx = cos(phi);
		const double dy = sin(phi);
		
		getBasis(dx * ux + dy * uy, &accSHbasis[i*cc]);
	}
	
	return accSHbasis;
}

double Blob::getOffsetAcc(double dx, double dy, const std::vector<double>& accSH)
{
	if (dx == 0.0 && dy == 0.0) return 0.0;
	
	const int bins = accSH.size();
	
	double phi = atan2(dy,dx);
	if (phi < 0.0) phi += 2.0 * PI;
	
	const double id = bins * phi / (2.0 * PI);
	
	const int i0 = (int)id;
	const int i1 = (i0+1) % bins;
	const double di = id - i0;
	
	return (1.0 - di) * accSH[i0] + di * accSH[i1];
}


double Blob::getBasisAcc(double dx, double dy, int b, const std::vector<double>& accSHbasis)
{
	if (dx == 0.0 && dy == 0.0) return 0.0;
	
	const int cc = shCoeffs.size();
	const int bins = accSHbasis.size() / cc;
	
	double phi = atan2(dy,dx);
	if (phi < 0.0) phi += 2.0 * PI;
	
	const double id = bins * phi / (2.0 * PI);
	
	const int i0 = (int)id;
	const int i1 = (i0+1) % bins;
	const double di = id - i0;
	
	return (1.0 - di) * accSHbasis[i0*cc + b] + di * accSHbasis[i1*cc + b];
	
}

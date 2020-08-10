#include "blob.h"
#include <src/spherical-harmonics/SphericalHarmonics.h>

#include <omp.h>

using namespace gravis;


Blob::Blob()
:	center(0.0), 
	outer_radius(0), 
	sphericalHarmonics(0),
	shCoeffs(0), shBands(0)
{}

Blob::Blob(d3Vector center, int outer_radius)
:	center(center), 
	outer_radius(outer_radius), 
	sphericalHarmonics(0),
	shCoeffs(0), shBands(0)
{
}
	
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
		const RawImage<float>& frame,
		const d4Matrix& proj,
		const RawImage<float>& weight,
		int radius)
{
	if (radius < 0) radius = outer_radius + 1;

	std::vector<double> radAvg(radius, 0.0), radCnt(radius, 0.0);

	const int w = frame.xdim;
	const int h = frame.ydim;

	d4Vector pw(center);
	d4Vector pi = proj * pw;

	d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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


		const double m = weight(x,y);

		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		radAvg[r0] += m * (1 - dr) * frame(x,y);
		radCnt[r0] += m * (1 - dr);

		radAvg[r1] += m * dr * frame(x,y);
		radCnt[r1] += m * dr;
	}

	for (int i = 0; i < radius; i++)
	{
		if (radCnt[i] > 0.0) radAvg[i] /= radCnt[i];
	}

	return radAvg;
}

double Blob::radialAverageError(
		const RawImage<float>& frame,
		const d4Matrix& proj,
		const RawImage<float>& weight,
		const std::vector<double>& radAvg)
{
	double out(0.0), totWgh(0.0);

	const int radius = radAvg.size();

	const int w = frame.xdim;
	const int h = frame.ydim;


	d4Vector pw(center);
	d4Vector pi = proj * pw;

	d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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

			obs = frame(xa,ya);
			m = 1e-4;
		}
		else
		{
			obs = frame(x,y);
			m = weight(x,y);
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

BufferedImage<float> Blob::radialAverageProjection(
		const RawImage<float>& frame,
		const gravis::d4Matrix& proj,
		const std::vector<double>& radAvg)
{
	BufferedImage<float> out(frame.xdim, frame.ydim);

	const int radius = radAvg.size();

	const int w = frame.xdim;
	const int h = frame.ydim;


	d4Vector pw(center);
	d4Vector pi = proj * pw;

	d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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
			out(x,y) = frame(x,y);
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
	RawImage<float>& frame,
	const d4Matrix& proj,
	const RawImage<float>& weight,
	double taper)
{
	const int radius = outer_radius + taper;

	std::vector<double> radAvg = radialAverage(frame, proj, weight, radius);

	const int w = frame.xdim;
	const int h = frame.ydim;

	d4Vector pw(center);
	d4Vector pi = proj * pw;

	d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

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

		outside_val += (1 - wgh) * frame(x,y);
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

		frame(x,y) -= wgh * pred;
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

double Blob::getOffset(d3Vector v)
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

void Blob::getBasis(d3Vector v, double *dest)
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

std::vector<double> Blob::accelerate(d3Vector ux, d3Vector uy, int bins)
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

std::vector<double> Blob::accelerateBasis(d3Vector ux, d3Vector uy, int bins)
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

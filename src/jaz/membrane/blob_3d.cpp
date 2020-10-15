#include "blob_3d.h"

#include <omp.h>

using namespace gravis;


Blob3D::Blob3D()
:	center(0.0), 
	outer_radius(0), 
	sphericalHarmonics(0),
	shCoeffs(0), shBands(0)
{}

Blob3D::Blob3D(d3Vector center, int outer_radius)
:	center(center), 
	outer_radius(outer_radius), 
	sphericalHarmonics(0),
	shCoeffs(0), shBands(0)
{
}
	
Blob3D::Blob3D(const std::vector<double>& params, int outer_radius,
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

	if (shCoeffs.size() > 0) shCoeffs[0] = 0.0;
}

std::vector<double> Blob3D::radialAverage(
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

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double dx = x - pi.x;
		const double dy = y - pi.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH), radius);

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

double Blob3D::radialAverageError(
		const RawImage<float>& frame,
		const d4Matrix& proj,
		const RawImage<float>& weight,
		const std::vector<double>& radAvg)
{
	double out = 0.0;

	const int radius = radAvg.size();

	const int w = frame.xdim;
	const int h = frame.ydim;


	d4Vector pw(center);
	d4Vector pi = proj * pw;

	d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

	std::vector<double> accSH = accelerate(dwdx, dwdy, 2.0 * PI * radius);

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double obs = frame(x,y);
		const double m = weight(x,y);

		const double dx = x - pi.x;
		const double dy = y - pi.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH), radius);

		if (r >= radius - 1 - 1e-6) r = radius - 1 - 1e-6;
		if (r < 0) r = 0;

		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		const double pred = (1 - dr) * radAvg[r0] + dr * radAvg[r1];
		const double err = pred - obs;

		out += m * err * err;
	}

	return out;
}

BufferedImage<float> Blob3D::drawError(
		const RawImage<float>& frame,
		const d4Matrix& proj,
		const RawImage<float>& weight,
		const std::vector<double>& radAvg)
{
	const int radius = radAvg.size();

	const int w = frame.xdim;
	const int h = frame.ydim;

	BufferedImage<float> out(w,h);

	d4Vector pw(center);
	d4Vector pi = proj * pw;

	d3Vector dwdx(proj(0,0), proj(0,1), proj(0,2));
	d3Vector dwdy(proj(1,0), proj(1,1), proj(1,2));

	std::vector<double> accSH = accelerate(dwdx, dwdy, 2.0 * PI * radius);

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double obs = frame(x,y);
		const double m = weight(x,y);

		const double dx = x - pi.x;
		const double dy = y - pi.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH), radius);

		if (r >= radius - 1 - 1e-6) r = radius - 1 - 1e-6;
		if (r < 0) r = 0;

		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		const double pred = (1 - dr) * radAvg[r0] + dr * radAvg[r1];

		out(x,y) = sqrt(m) * (pred - obs);
	}

	return out;
}

BufferedImage<float> Blob3D::radialAverageProjection(
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

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double dx = x - pi.x;
		const double dy = y - pi.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffsetAcc(dx, dy, accSH), radius);

		if (r >= radius - 1 - 1e-6) r = radius - 1 - 1e-6;
		if (r < 0) r = 0;

		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		out(x,y) = (1 - dr) * radAvg[r0] + dr * radAvg[r1];
	}

	return out;
}

void Blob3D::decompose(
	RawImage<float>& frame,
	RawImage<float>& blob,
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
		const double r = smoothOrigin(ru + getOffsetAcc(dx, dy, accSH), radius);

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
		const double r = smoothOrigin(ru + getOffsetAcc(dx, dy, accSH), radius);

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
		blob(x,y)  += wgh * pred;
	}
}


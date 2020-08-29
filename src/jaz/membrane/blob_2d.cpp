#include "blob_2d.h"

#include <omp.h>

using namespace gravis;


Blob2D::Blob2D()
:	center(0.0), 
	outer_radius(0), 
	amplitudes(0)
{}

Blob2D::Blob2D(d2Vector center, int outer_radius)
:	center(center), 
	outer_radius(outer_radius), 
	amplitudes(0)
{
}
	
Blob2D::Blob2D(const std::vector<double>& params, int outer_radius)
:	center(params[0], params[1]),
	outer_radius(outer_radius),
	amplitudes((params.size() - 2)/2)
{
	for (int i = 0; i < amplitudes.size(); i++)
	{
		amplitudes[i].real = params[2*i+2];
		amplitudes[i].imag = params[2*i+3];
	}
}

std::vector<double> Blob2D::radialAverage(
		const RawImage<float>& frame,
		const RawImage<float>& weight,
		int radius)
{
	if (radius < 0) radius = outer_radius + 1;

	std::vector<double> radAvg(radius, 0.0), radCnt(radius, 0.0);

	const int w = frame.xdim;
	const int h = frame.ydim;


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double dx = x - center.x;
		const double dy = y - center.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)), radius);

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

double Blob2D::radialAverageError(
		const RawImage<float>& frame,
		const RawImage<float>& weight,
		const std::vector<double>& radAvg)
{
	double out = 0.0;

	const int radius = radAvg.size();

	const int w = frame.xdim;
	const int h = frame.ydim;


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double obs = frame(x,y);
		const double m = weight(x,y);

		const double dx = x - center.x;
		const double dy = y - center.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)), radius);

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

BufferedImage<float> Blob2D::drawError(
		const RawImage<float>& frame,
		const RawImage<float>& weight,
		const std::vector<double>& radAvg)
{
	const int radius = radAvg.size();

	const int w = frame.xdim;
	const int h = frame.ydim;

	BufferedImage<float> out(w,h);


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double obs = frame(x,y);
		const double m = weight(x,y);

		const double dx = x - center.x;
		const double dy = y - center.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)), radius);

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

BufferedImage<float> Blob2D::radialAverageProjection(
		const RawImage<float>& frame,
		const std::vector<double>& radAvg)
{
	BufferedImage<float> out(frame.xdim, frame.ydim);

	const int radius = radAvg.size();

	const int w = frame.xdim;
	const int h = frame.ydim;

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double dx = x - center.x;
		const double dy = y - center.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)), radius);

		if (r >= radius - 1 - 1e-6) r = radius - 1 - 1e-6;
		if (r < 0) r = 0;

		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		out(x,y) = (1 - dr) * radAvg[r0] + dr * radAvg[r1];
	}

	return out;
}

void Blob2D::decompose(
	RawImage<float>& frame,
	RawImage<float>& blob,
	const RawImage<float>& weight,
	double taper)
{
	const int radius = outer_radius + taper;

	std::vector<double> radAvg = radialAverage(frame, weight, radius);

	const int w = frame.xdim;
	const int h = frame.ydim;


	const int x0 = (int) (center.x + 0.5);
	const int y0 = (int) (center.y + 0.5);

	double outside_val(0.0), outside_wgh(0.0);

	for (int y = y0 - radius + 1; y <= y0 + radius - 1; y++)
	for (int x = x0 - radius + 1; x <= x0 + radius - 1; x++)
	{
		if (x < 0 || x >= w || y < 0 || y >= h) continue;

		const double dx = x - center.x;
		const double dy = y - center.y;

		const double ru = sqrt(dx*dx + dy*dy);
		double r = smoothOrigin(ru + getOffset(d2Vector(dx, dy)), radius);

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

		const double dx = x - center.x;
		const double dy = y - center.y;

		const double ru = sqrt(dx*dx + dy*dy);
		double r = smoothOrigin(ru + getOffset(d2Vector(dx, dy)), radius);

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

double Blob2D::scanForMinimalRadius(int samples)
{
	double min = std::numeric_limits<double>::max();
	
	for (int i = 0; i < samples; i++)
	{
		double phi = 2 * PI * i / (double) samples;
		
		double value = getOffset(phi);
		
		if (value < min) min = value;
	}
	
	return min;
}

double Blob2D::scanForMaximalRadius(int samples)
{
	double max = -std::numeric_limits<double>::max();
	
	for (int i = 0; i < samples; i++)
	{
		double phi = 2 * PI * i / (double) samples;
		
		double value = getOffset(phi);
		
		if (value > max) max = value;
	}
	
	return max;
}

std::vector<double> Blob2D::rotate(const std::vector<double> &params, double angle, d2Vector axis)
{
	std::vector<double> out = params;
	const double cos_psi = cos(angle);
	const double sin_psi = sin(angle);

	out[0] = axis.x + cos_psi * (params[0] - axis.x) + sin_psi * (params[1] - axis.y);
	out[1] = axis.y - sin_psi * (params[0] - axis.x) + cos_psi * (params[1] - axis.y);

	const int Fourier_params = params.size() - 2;

	for (int i = 0; i < Fourier_params/2; i++)
	{
		const int n = i + FIRST_BLOB_FREQUENCY;
		const double cos_n_psi = cos(n * angle);
		const double sin_n_psi = sin(n * angle);

		const int j = 2 + 2 * i;

		out[j]   = cos_n_psi * params[j]   + sin_n_psi * params[j+1];
		out[j+1] = cos_n_psi * params[j+1] - sin_n_psi * params[j];
	}

	return out;
}


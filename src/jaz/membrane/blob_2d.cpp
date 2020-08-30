#include "blob_2d.h"

#include <omp.h>

using namespace gravis;


Blob2D::Blob2D()
:	center(0.0), 
	amplitudes(0)
{}

Blob2D::Blob2D(d2Vector center, double smoothingRadius)
:	center(center), 
	amplitudes(0),
    smoothingRadius(smoothingRadius)
{
}
	
Blob2D::Blob2D(const std::vector<double>& params, double smoothingRadius)
:	center(params[0], params[1]),
	amplitudes((params.size() - 2)/2),
    smoothingRadius(smoothingRadius)
{
	for (int i = 0; i < amplitudes.size(); i++)
	{
		amplitudes[i].real = params[2*i+2];
		amplitudes[i].imag = params[2*i+3];
	}
}

std::vector<double> Blob2D::radialAverage(
		const RawImage<float>& frame,
		const RawImage<float>& weight) const
{
	int maxRadius = findMaxRadius(i2Vector(frame.xdim, frame.ydim));

	std::vector<double> radAvg(maxRadius, 0.0), radCnt(maxRadius, 0.0);

	const int w = frame.xdim;
	const int h = frame.ydim;


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double dx = x - center.x;
		const double dy = y - center.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)));

		if (r >= maxRadius - 1 || r < 0) continue;

		const double m = weight(x,y);

		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		radAvg[r0] += m * (1 - dr) * frame(x,y);
		radCnt[r0] += m * (1 - dr);

		radAvg[r1] += m * dr * frame(x,y);
		radCnt[r1] += m * dr;
	}

	for (int i = 0; i < maxRadius; i++)
	{
		if (radCnt[i] > 0.0) radAvg[i] /= radCnt[i];
	}

	return radAvg;
}

double Blob2D::radialAverageError(
		const RawImage<float>& frame,
		const RawImage<float>& weight,
		const std::vector<double>& radAvg) const
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

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)));

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
		const std::vector<double>& radAvg) const
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

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)));

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
		const std::vector<double>& radAvg) const
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

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)));

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
	double outerRadius, double taper) const
{

	std::vector<double> radAvg = radialAverage(frame, weight);
	const int radius = radAvg.size();

	const int w = frame.xdim;
	const int h = frame.ydim;


	const int x0 = (int) (center.x + 0.5);
	const int y0 = (int) (center.y + 0.5);

	double outside_val(0.0), outside_wgh(0.0);

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double dx = x - center.x;
		const double dy = y - center.y;

		const double ru = sqrt(dx*dx + dy*dy);
		double r = smoothOrigin(ru + getOffset(d2Vector(dx, dy)));

		if (r >= radius - 1 || r < 0)
		{
			continue;
		}

		double wgh;

		if (ru < outerRadius) wgh = 1.0;
		else if (ru < outerRadius + taper) wgh = cos(0.5 * PI * (ru - outerRadius) / taper);
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
		double r = smoothOrigin(ru + getOffset(d2Vector(dx, dy)));

		if (r >= radius - 1 || r < 0)
		{
			continue;
		}

		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		const double pred = (1 - dr) * radAvg[r0] + dr * radAvg[r1] - outside_val;

		double wgh;

		if (ru < outerRadius) wgh = 1.0;
		else if (ru < outerRadius + taper) wgh = cos(0.5 * PI * (ru - outerRadius) / taper);
		else wgh = 0.0;

		frame(x,y) -= wgh * pred;
		blob(x,y)  += wgh * pred;
	}
}

double Blob2D::scanForMinimalRadius(int samples) const
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

double Blob2D::scanForMaximalRadius(int samples) const
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

std::pair<d2Vector, d2Vector> Blob2D::scanForBoundingBox(double radius, double padding, int samples) const
{
	d2Vector p0 = center;
	d2Vector p1 = center;
	
	for (int i = 0; i < samples; i++)
	{
		const double phi = 2 * PI * i / (double) samples;		
		const double rad = radius + getOffset(phi);
		const d2Vector r = center + rad * d2Vector(cos(phi), sin(phi));
		
		if (p0.x > r.x) p0.x = r.x;
		if (p0.y > r.y) p0.y = r.y;
		if (p1.x < r.x) p1.x = r.x;
		if (p1.y < r.y) p1.y = r.y;
	}
	
	const d2Vector v0 = d2Vector(
	            p0.x - padding, 
				p0.y - padding);
	
	const d2Vector v1 = d2Vector(
	            p1.x + padding, 
				p1.y + padding );
	
	return std::make_pair(v0,v1);
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


DelineatedBlob2D::DelineatedBlob2D(d2Vector center, double radius, double smoothing_radius)
:	blob(center, smoothing_radius),
	radius(radius)	
{
	
}

DelineatedBlob2D::DelineatedBlob2D(const std::vector<double>& params)
:	radius(params[0]),
	blob(stripRadius(params), 2 * params[0])
{
}

std::vector<DelineatedBlob2D> DelineatedBlob2D::read(const std::string &filename)
{
	std::vector<DelineatedBlob2D> out;
	
	std::ifstream file(filename);
	
	if (!file)
	{
		REPORT_ERROR("DelineatedBlob2D::read: unable to read " + filename);
	}
	
	std::string line;
	
	while (std::getline(file, line))
	{
		std::stringstream sts;
		sts << line;
		
		std::vector<double> values;
		
		while (sts)
		{
			double value;
			sts >> value;
			values.push_back(value);
		}
		
		out.push_back(DelineatedBlob2D(values));
	}
	
	return out;
}

std::vector<double> DelineatedBlob2D::stripRadius(const std::vector<double> &params)
{
	std::vector<double> blob_params(params.size()-1);
	
	for (int i = 0; i < params.size()-1; i++)
	{
		blob_params[i] = params[i+1];
	}
	
	return blob_params;
}


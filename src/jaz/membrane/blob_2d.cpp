#include "blob_2d.h"

#include <omp.h>

using namespace gravis;


Blob2D::Blob2D()
:	center(0.0), 
	amplitudes(0)
{
}

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
		const RawImage<float>& weight,
        double relevantRadius) const
{
	int maxRadius = findMaxRadius(i2Vector(frame.xdim, frame.ydim));
	
	if (relevantRadius > 0 && maxRadius > relevantRadius) maxRadius = relevantRadius;

	std::vector<double> radAvg(maxRadius, 0.0), radCnt(maxRadius, 0.0);

	const int w = frame.xdim;
	const int h = frame.ydim;


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double m = weight(x,y);
		
		if (m == 0.f) continue;
		
		const double dx = x - center.x;
		const double dy = y - center.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)));

		if (r >= maxRadius - 1 || r < 0) continue;


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


std::pair<std::vector<double>,std::vector<double>> Blob2D::radialAverageAndWeight(
		const RawImage<float>& frame,
		const RawImage<float>& weight,
        double relevantRadius) const
{
	int maxRadius = findMaxRadius(i2Vector(frame.xdim, frame.ydim));
	
	if (relevantRadius > 0 && maxRadius > relevantRadius) maxRadius = relevantRadius;

	std::vector<double> radAvg(maxRadius, 0.0), radCnt(maxRadius, 0.0);

	const int w = frame.xdim;
	const int h = frame.ydim;


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double m = weight(x,y);
		
		if (m == 0.f) continue;
		
		const double dx = x - center.x;
		const double dy = y - center.y;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(d2Vector(dx, dy)));

		if (r >= maxRadius - 1 || r < 0) continue;


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

	return std::make_pair(radAvg, radCnt);
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
		const double m = weight(x,y);
		
		if (m == 0.f) continue;
		
		const double obs = frame(x,y);
		

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
	out.fill(0.f);


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double m = weight(x,y);
		
		if (m == 0.f) continue;
		
		const double obs = frame(x,y);

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

void Blob2D::erase(
	const RawImage<float>& micrographs,
	RawImage<float>& erased_out,
	RawImage<float>& blob_out,
	const RawImage<float>& weight,
	double radius, 
    double taper) const
{
	std::pair<std::vector<double>,std::vector<double>> radAvgAndWgh = radialAverageAndWeight(
	            micrographs, weight, radius + taper);
	
	const double cappingRadius = 1.5 * smoothingRadius;
	
	std::vector<double> radAvg = radAvgAndWgh.first;
	std::vector<double> radWgh = radAvgAndWgh.second;
	
	int first_r = (int) cappingRadius;
	
	for (int r = 0; r < cappingRadius; r++)
	{
		if (radWgh[r] > 0)
		{
			first_r = r;
			break;
		}
	}
	
	const double cappingRange = cappingRadius - first_r;

	double tipAvg = 0.0;
	double tipWgh = 0.0;
	
	for (int r = first_r; r < cappingRadius; r++)
	{
		const int rr = r - first_r;
		const double t = radWgh[r] * (cos(PI * rr / cappingRange) + 1.0) / 2;
		
		tipAvg += t * radAvg[r];
		tipWgh += t;
	}
	
	tipAvg /= tipWgh;
	
	for (int r = first_r; r < cappingRadius; r++)
	{
		const int rr = r - first_r;
		const double t = (cos(PI * rr / cappingRange) + 1.0) / 2;
		
		radAvg[r] = (1 - t) * radAvg[r] + t * tipAvg;
	}

	
	const int w = micrographs.xdim;
	const int h = micrographs.ydim;

	double outside_val(0.0), outside_wgh(0.0);
	
	std::pair<gravis::d2Vector,gravis::d2Vector> boundingBox = scanForBoundingBox(
	        radius + taper, 2, 2 * PI * (radius+taper));
	
	const d2Vector minPos = boundingBox.first;
	const d2Vector maxPos = boundingBox.second;
	
	const int x0 = minPos.x < 0?       0 : std::ceil(minPos.x);
	const int x1 = maxPos.x >= w-1?  w-1 : std::floor(maxPos.x);
	const int y0 = minPos.y < 0?       0 : std::ceil(minPos.y);
	const int y1 = maxPos.y >= h-1?  h-1 : std::floor(maxPos.y);

	for (int y = y0; y <= y1; y++)
	for (int x = x0; x <= x1; x++)
	{
		const double dx = x - center.x;
		const double dy = y - center.y;

		const double ru = sqrt(dx*dx + dy*dy);
		double r = smoothOrigin(ru + getOffset(d2Vector(dx, dy)));

		if (r >= radius + taper || r < 0)
		{
			continue;
		}

		double wgh;

		if (ru < radius) wgh = 1.0;
		else if (ru < radius + taper) wgh = (cos(PI * (ru - radius) / taper) + 1) / 2;
		else wgh = 0.0;

		if (wgh > 0.0)
		{
			outside_val += (1 - wgh) * micrographs(x,y);
			outside_wgh += (1 - wgh);
		}
	}

	if (outside_wgh > 0.0) outside_val /= outside_wgh;
	
	for (int y = y0; y <= y1; y++)
	for (int x = x0; x <= x1; x++)
	{
		if (x < 0 || x >= w || y < 0 || y >= h) continue;

		const double dx = x - center.x;
		const double dy = y - center.y;

		const double ru = sqrt(dx*dx + dy*dy);
		double r = smoothOrigin(ru + getOffset(d2Vector(dx, dy)));

		if (r >= radius + taper - 1 || r < 0)
		{
			continue;
		}

		const int r0 = (int)r;
		const int r1 = r0 + 1;
		
		if (radWgh[r0] == 0 || radWgh[r1] == 0)
		{
			continue;
		}

		const double dr = r - r0;

		const double pred = (1 - dr) * radAvg[r0] + dr * radAvg[r1] - outside_val;

		double wgh;

		if (ru < radius) wgh = 1.0;
		else if (ru < radius + taper) wgh = (cos(PI * (ru - radius) / taper) + 1) / 2;
		else wgh = 0.0;
		
		const double blob_value = wgh * pred;

		erased_out(x,y) -= blob_value;
		blob_out(x,y)   += blob_value;
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


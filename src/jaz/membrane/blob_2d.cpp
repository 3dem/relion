#include "blob_2d.h"
#include <src/jaz/image/interpolation.h>
#include <src/jaz/math/fft.h>

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
	
	
	double maxAmp = 0.0;
	
	for (int i = 0; i < amplitudes.size(); i++)
	{
		maxAmp += amplitudes[i].abs();
	}
	
	int x0 = std::ceil( center.x - maxRadius - maxAmp);
	int x1 = std::floor(center.x + maxRadius + maxAmp)+1;
	int y0 = std::ceil( center.y - maxRadius - maxAmp);
	int y1 = std::floor(center.y + maxRadius + maxAmp)+1;
	
	if (x0 < 0) x0 = 0;
	if (x1 > w) x1 = w;
	
	if (y0 < 0) y0 = 0;
	if (y1 > h) y1 = h;
	
	for (int y = y0; y < y1; y++)
	for (int x = x0; x < x1; x++)
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



std::pair<std::vector<double>,std::vector<double>> Blob2D::radialAverageAndWeightInSectors(
		const RawImage<float>& frame,
		const RawImage<float>& weight,
        int sectors,
        double relevantRadius) const
{
	int maxRadius = findMaxRadius(i2Vector(frame.xdim, frame.ydim));
	
	if (relevantRadius > 0 && maxRadius > relevantRadius) maxRadius = relevantRadius;

	std::vector<double> radAvg(sectors*maxRadius, 0.0), radCnt(sectors*maxRadius, 0.0);

	const int w = frame.xdim;
	const int h = frame.ydim;


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double m = weight(x,y);
		
		if (m == 0.f) continue;
		
		const double dx = x - center.x;
		const double dy = y - center.y;
		
		const double phi = atan2(dx, dy);
		
		const double tau = sectors * (phi + PI) / (2*PI);
		const int sector = ((int)std::round(tau)) % sectors;

		double r = smoothOrigin(sqrt(dx*dx + dy*dy) + getOffset(phi));

		if (r >= maxRadius - 1 || r < 0) continue;


		const int r0 = (int)r;
		const int r1 = r0 + 1;

		const double dr = r - r0;

		radAvg[sector * maxRadius + r0] += m * (1 - dr) * frame(x,y);
		radCnt[sector * maxRadius + r0] += m * (1 - dr);

		radAvg[sector * maxRadius + r1] += m * dr * frame(x,y);
		radCnt[sector * maxRadius + r1] += m * dr;
	}

	for (int i = 0; i < sectors * maxRadius; i++)
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
	const RawImage<float>& micrograph,
	RawImage<float>& erased_out,
	RawImage<float>& blob_out,
	const RawImage<float>& weight,
	double radius, 
    double taper) const
{
	std::pair<std::vector<double>,std::vector<double>> radAvgAndWgh = radialAverageAndWeight(
	            micrograph, weight, radius + taper);
	
	const int max_radius = radAvgAndWgh.first.size();
	
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

	
	const int w = micrograph.xdim;
	const int h = micrograph.ydim;

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

		if (r >= radius + taper || r >= max_radius - 1 || r < 0)
		{
			continue;
		}

		double wgh;

		if (ru < radius) wgh = 1.0;
		else if (ru < radius + taper) wgh = (cos(PI * (ru - radius) / taper) + 1) / 2;
		else wgh = 0.0;

		if (wgh > 0.0)
		{
			outside_val += (1 - wgh) * micrograph(x,y);
			outside_wgh += (1 - wgh);
		}
	}

	if (outside_wgh > 0.0) outside_val /= outside_wgh;
	
	for (int y = y0; y <= y1; y++)
	for (int x = x0; x <= x1; x++)
	{
		const double dx = x - center.x;
		const double dy = y - center.y;

		const double ru = sqrt(dx*dx + dy*dy);
		double r = smoothOrigin(ru + getOffset(d2Vector(dx, dy)));

		if (r >= radius + taper - 1 || r >= max_radius - 1 || r < 0)
		{
			continue;
		}

		const int r0 = (int)r;
		const int r1 = r0 + 1;
		
		const double eps = 1e-6;
		
		if (radWgh[r0] < eps || radWgh[r1] < eps)
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


void Blob2D::eraseLocally(
	const RawImage<float>& micrograph,
	RawImage<float>& erased_out,
	RawImage<float>& blob_out,
	const RawImage<float>& weight,
	double radius, 
    double taper,
	double smoothness) const
{
	std::pair<BufferedImage<float>, BufferedImage<float>> 
		polarAndMask =	transformToPolar(micrograph, weight, radius + taper);
		
	BufferedImage<float> polarImage = polarAndMask.first;
	BufferedImage<float> polarWeight = polarAndMask.second;
		
	const int phi_samples = polarImage.xdim;
	const int max_radius = polarImage.ydim;
	
	const bool debug = false;
	
	if (debug)
	{
		polarImage.write("DEBUG_polarImage_0.mrc");
		polarWeight.write("DEBUG_polarWeight_0.mrc");
	}
	
	BufferedImage<float> polarImageNormalized;
	        
	const bool no_shear = true;
	
	if (no_shear)
	{
		polarImageNormalized.resize(phi_samples, max_radius);
		
		BufferedImage<float> image_curve_RS(phi_samples, 1);
		BufferedImage<float> weight_curve_RS(phi_samples, 1);
		BufferedImage<fComplex> image_curve_FS(phi_samples/2 + 1, 1);
		BufferedImage<fComplex> weight_curve_FS(phi_samples/2 + 1, 1);
		
		for (int y = 0; y < max_radius; y++)
		{
			const double sigma_RS = smoothness * (radius+1) / (double) (y+1);
			const double sigma_FS = PI * phi_samples / sigma_RS;			
			const double sig22 = 2.0 * sigma_FS * sigma_FS;
			        
			for (int x = 0; x < phi_samples; x++)
			{
				image_curve_RS(x,0) = polarImage(x,y);
				weight_curve_RS(x,0) = polarWeight(x,y);
			}
			
			FFT::FourierTransform(image_curve_RS, image_curve_FS);
			FFT::FourierTransform(weight_curve_RS, weight_curve_FS);
			
			for (int x = 0; x < phi_samples/2 + 1; x++)
			{
				image_curve_FS(x,0) *= exp(-x*x/sig22);
				weight_curve_FS(x,0) *= exp(-x*x/sig22);
			}
			
			FFT::inverseFourierTransform(image_curve_FS, image_curve_RS);
			FFT::inverseFourierTransform(weight_curve_FS, weight_curve_RS);
			
			const double eps = 1e-3;
			
			for (int x = 0; x < phi_samples; x++)
			{
				polarImage(x,y) = image_curve_RS(x,0);
				polarWeight(x,y) = weight_curve_RS[x] > eps? weight_curve_RS[x] : eps;
			}
		}
	}
	else
	{
		BufferedImage<fComplex> polarImageFS, polarWeightFS;
		
		FFT::FourierTransform(polarImage, polarImageFS);
		FFT::FourierTransform(polarWeight, polarWeightFS);
		
		const double sigma = PI * phi_samples / smoothness;
		const double sig22 = 2.0 * sigma * sigma;
				
		for (int y = 0; y < max_radius; y++)
		for (int x = 0; x < phi_samples/2 + 1; x++)
		{
			polarImageFS(x,y)  *= exp(-x*x/sig22);
			polarWeightFS(x,y) *= exp(-x*x/sig22);
		}
		
		FFT::inverseFourierTransform(polarImageFS, polarImage);
		FFT::inverseFourierTransform(polarWeightFS, polarWeight);
		
		const double eps = 1e-3;
		
		for (int i = 0; i < polarWeight.getSize(); i++)
		{
			if (polarWeight[i] < eps)
			{
				polarWeight[i] = eps;
			}
		}
	}
		
	if (debug)
	{
		polarImage.write("DEBUG_polarImage_1.mrc");
		polarWeight.write("DEBUG_polarWeight_1.mrc");
	}
	
	polarImageNormalized = polarImage / polarWeight;
	
	
	if (debug)
	{
		polarImageNormalized.write("DEBUG_polarImageNormalized.mrc");
	}
	
	const double cappingRadius = 1.75 * smoothingRadius;
	
	double tipAvg = 0.0;
	double tipWgh = 0.0;
	
	for (int p = 0; p < phi_samples; p++)
	for (int r = 0; r < cappingRadius && r < max_radius; r++)
	{
		const double t = (cos(PI * r / cappingRadius) + 1.0) / 2;
		
		tipAvg += t * polarImage(p,r);
		tipWgh += t * polarWeight(p,r);
	}
	
	tipAvg /= tipWgh;
	
	for (int p = 0; p < phi_samples; p++)
	for (int r = 0; r < cappingRadius && r < max_radius; r++)
	{
		const double t = (cos(PI * r / cappingRadius) + 1.0) / 2;
		
		polarImageNormalized(p,r) = t * tipAvg + (1-t) * polarImageNormalized(p,r);
	}
	
	if (debug)
	{
		polarImageNormalized.write("DEBUG_polarImageNormalized_2.mrc");
	}
	
	const int w = micrograph.xdim;
	const int h = micrograph.ydim;

	double outside_val(0.0), outside_wgh(0.0);
	
	std::pair<gravis::d2Vector,gravis::d2Vector> boundingBox = scanForBoundingBox(
	        radius + taper, 2, 2 * PI * (radius+taper));
	
	const d2Vector minPos = boundingBox.first;
	const d2Vector maxPos = boundingBox.second;
	
	const int x0 = minPos.x < 0?       0 : std::ceil(minPos.x);
	const int x1 = maxPos.x >= w-1?  w-1 : std::floor(maxPos.x);
	const int y0 = minPos.y < 0?       0 : std::ceil(minPos.y);
	const int y1 = maxPos.y >= h-1?  h-1 : std::floor(maxPos.y);
	
	const double eps2 = 0.001;

	for (int y = y0; y <= y1; y++)
	for (int x = x0; x <= x1; x++)
	{
		const double dx = x - center.x;
		const double dy = y - center.y;
		
		const double ru = sqrt(dx*dx + dy*dy);
		const double phi = atan2(dy, dx);
		double r = smoothOrigin(ru + getOffset(phi));

		if (r >= radius + taper || r >= max_radius - 1  || r < 0)
		{
			continue;
		}
		
		double wgh = weight(x,y) * Interpolation::cubicXY_wrap(
		            polarWeight, 
		            phi_samples * phi / (2 * PI), 
		            r);
		
		if (wgh < eps2) continue;
		
		if (ru > radius)
		{
			if (ru < radius + taper) wgh *= (cos(PI * (ru - radius) / taper) + 1) / 2;
			else wgh = 0.0;
		}

		if (wgh > 0.0)
		{
			outside_val += (1 - wgh) * micrograph(x,y);
			outside_wgh += (1 - wgh);
		}
	}

	if (outside_wgh > 0.0) outside_val /= outside_wgh;
	
	for (int y = y0; y <= y1; y++)
	for (int x = x0; x <= x1; x++)
	{
		const double dx = x - center.x;
		const double dy = y - center.y;

		double phi = atan2(dy, dx) + 2 * PI;
		if (phi >= 2*PI) phi -= 2*PI;
		
		const double ru = sqrt(dx*dx + dy*dy);
		double r = smoothOrigin(ru + getOffset(phi));

		if (r >= radius + taper - 1 || r >= max_radius - 1 || r < 0)
		{
			continue;
		}
		
		const double pred = Interpolation::cubicXY_wrap(
		            polarImageNormalized, 
		            phi_samples * phi / (2 * PI), 
		            r);
		
		double wgh = weight(x,y) * Interpolation::cubicXY_wrap(
		            polarWeight, 
		            phi_samples * phi / (2 * PI), 
		            r);
		
		if (wgh < eps2) continue;
		
		if (ru > radius)
		{
			if (ru < radius + taper) wgh *= (cos(PI * (ru - radius) / taper) + 1) / 2;
			else wgh = 0.0;
		}
		
		const double blob_value = wgh * (pred - outside_val);

		if (blob_value == blob_value)
		{
			erased_out(x,y) -= blob_value;
			blob_out(x,y)   += blob_value;
		}
	}
}

std::pair<BufferedImage<float>, BufferedImage<float>> 
Blob2D::transformToPolar(
        const RawImage<float> &frame,
        const RawImage<float> &mask, 
        double maxRadius) const
{
	const double len = 2 * PI * maxRadius;
	        
	const int h_out = 2 * (int) (maxRadius / 2.0) + 2;
	const int w_out = 2 * (int) (len / 2.0);
	
	BufferedImage<float> out_data(w_out, h_out);
	BufferedImage<float> out_weight(w_out, h_out);
	
	for (int y = 0; y < h_out; y++)
	for (int x = 0; x < w_out; x++)
	{
		const double phi = 2 * PI * x / (double) w_out;
		const double r = y - getOffset(phi);		
		
		const d2Vector pos = center + r * d2Vector(cos(phi), sin(phi));
		
		if (frame.containsPoint(pos))
		{
			const float wgh = Interpolation::cubicXY_clip(mask, pos.x, pos.y);
			out_data(x,y) = wgh * Interpolation::cubicXY_clip(frame, pos.x, pos.y);
			out_weight(x,y) = wgh;
		}
		else
		{
			out_data(x,y) = 0;
			out_weight(x,y) = 0;
		}
	}
	
	return std::make_pair(out_data, out_weight);
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
		const double rad = radius - getOffset(phi);
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


DelineatedBlob2D::DelineatedBlob2D()
{	
}

DelineatedBlob2D::DelineatedBlob2D(const Blob2D& blob, double radius)
:	Blob2D(blob),
	radius(radius)	
{
}

DelineatedBlob2D::DelineatedBlob2D(d2Vector center, double radius, double smoothing_radius)
:	Blob2D(center, smoothing_radius),
	radius(radius)	
{
}

DelineatedBlob2D::DelineatedBlob2D(const std::vector<double>& params)
:	Blob2D(stripRadius(params), 2 * params[0]),
	radius(params[0])	
{
}

double DelineatedBlob2D::getRadius(double phi) const
{
	return radius - getOffset(phi);
}

double DelineatedBlob2D::getSignedDistance(d2Vector imgPos) const
{
	return getDistance(imgPos) - radius;
}

double DelineatedBlob2D::getRelativeSignedDistance(d2Vector imgPos) const
{
	return (getDistance(imgPos) - radius) / radius;
}

d2Vector DelineatedBlob2D::getOutlinePoint(double phi) const
{
	const double r = getRadius(phi);
	return center + r * d2Vector(cos(phi), sin(phi));
}

d2Vector DelineatedBlob2D::estimateNormal(double phi, double scale) const
{
	d2Vector p0 = getOutlinePoint(phi - scale / 2.0);
	d2Vector p1 = getOutlinePoint(phi + scale / 2.0);
	
	return d2Vector(p1.y - p0.y, p0.x - p1.x).normalize();
}

double DelineatedBlob2D::perimeter() const
{
	const int samples = std::round(2 * PI * radius);
	
	double sum = 0.0;
	        
	for (int i = 0; i < samples; i++)
	{
		const double phi0 = 2 * PI * (i  ) / (double) samples;
		const double phi1 = 2 * PI * (i+1) / (double) samples;
		
		const double rad0 = radius - Blob2D::getOffset(phi0);
		const double rad1 = radius - Blob2D::getOffset(phi1);
		
		const d2Vector r0 = Blob2D::center + rad0 * d2Vector(cos(phi0), sin(phi0));
		const d2Vector r1 = Blob2D::center + rad1 * d2Vector(cos(phi1), sin(phi1));
		
		sum += (r1 - r0).length();
	}
	
	return sum;
}

Blob2D DelineatedBlob2D::getBlob2D() const
{
	Blob2D out = *this;
	return out; 
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
		
		if (values.size() >= 3 && values[0] > 0.0)
		{
			out.push_back(DelineatedBlob2D(values));
		}
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

std::vector<double> DelineatedBlob2D::addRadius(double radius, const std::vector<double> &params)
{
	std::vector<double> blob_params(params.size()+1);
	
	blob_params[0] = radius;
	
	for (int i = 0; i < params.size(); i++)
	{
		blob_params[i+1] = params[i];
	}
	
	return blob_params;
}


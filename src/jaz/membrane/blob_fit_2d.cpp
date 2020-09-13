#include "blob_fit_2d.h"
#include <src/jaz/tomography/fiducials.h>

using namespace gravis;


BlobFit2D::BlobFit2D(
	const RawImage<float>& image,
	d2Vector position,
	double smoothingRadius,
	double priorSigma,
	double roundedness,
	int num_threads)
:
		image(image),
		initialPos(position),
		smoothingRadius(smoothingRadius),
		num_threads(num_threads),
		priorSigma2(priorSigma*priorSigma),
        roundedness(roundedness)
{
	const int w = image.xdim;
	const int h = image.ydim;

	weight.resize(w,h);
	weight.fill(1.f);
}

double BlobFit2D::f(const std::vector<double>& x, void* tempStorage) const
{
	Blob2D blob(x, smoothingRadius);

	double out = 0.0;

	std::vector<double> radAvg = blob.radialAverage(image, weight);

	out += blob.radialAverageError(image, weight, radAvg) / (image.xdim * (double)image.ydim);

	if (priorSigma2 > 0.0)
	{

		out += (blob.center - initialPos).norm2() / priorSigma2;
	}
	
	if (roundedness > 0.0)
	{
		for (int i = 2; i < x.size(); i++)
		{
			const double f = (i - 2) / 2;
			out += roundedness * x[i] * x[i] * f * f;
		}
	}
	
	return out;
}

void BlobFit2D::computeWeight(const Blob2D& blob, double minRadius, double maxRadius,
                              const RawImage<float>& mask)
{
	const int w = weight.xdim;
	const int h = weight.ydim;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double r = blob.getDistance(d2Vector(x,y));
		
		if (r < minRadius || r > maxRadius)
		{
			weight(x,y) = 0.f;
		}
		else
		{
			weight(x,y) = mask(x,y);
		}
	}
}

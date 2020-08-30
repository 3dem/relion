#include "blob_fit_2d.h"
#include <src/jaz/tomography/fiducials.h>

using namespace gravis;


BlobFit2D::BlobFit2D(
	const RawImage<float>& image,
	d2Vector position,
	int frequencies,
	double smoothingRadius,
	double priorSigma,
	int num_threads)
:
		image(image),
		initialPos(position),
		smoothingRadius(smoothingRadius),
		frequencies(frequencies),
		num_threads(num_threads),
		priorSigma2(priorSigma*priorSigma)
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

	out += blob.radialAverageError(image, weight, radAvg);

	if (priorSigma2 > 0.0)
	{

		out += (blob.center - initialPos).norm2() / priorSigma2;
	}
	
	return out;
}

void BlobFit2D::computeWeight(const Blob2D& blob, double minRadius, double maxRadius)
{
	const int w = weight.xdim;
	const int h = weight.ydim;
	
	const float eps = 1e-6;
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double r = blob.getDistance(d2Vector(x,y));
		
		if (r > maxRadius || r < minRadius)
		{
			weight(x,y) = eps;
		}
		else
		{
			weight(x,y) = 1.f;
		}
	}
}

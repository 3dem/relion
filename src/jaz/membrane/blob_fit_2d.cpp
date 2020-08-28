#include "blob_fit_2d.h"
#include <src/jaz/tomography/fiducials.h>

using namespace gravis;


BlobFit2D::BlobFit2D(
	const RawImage<float>& image,
	d2Vector position,
	int frequencies,
	double radius,
	double thickness,
	double priorSigma,
	int num_threads)
:
		image(image),
		initialPos(position),
		outer_radius(radius + thickness/2),
		frequencies(frequencies),
		num_threads(num_threads),
		priorSigma2(priorSigma*priorSigma)
{
	const int w = image.xdim;
	const int h = image.ydim;

	const double maskedVal = 1e-6;

	weight.resize(w,h);
	weight.fill(1.f);

	const double outer_rad = radius + thickness / 2;


	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		const double dxb = position.x - x;
		const double dyb = position.y - y;
		const double distB = sqrt(dxb * dxb + dyb * dyb);

		const double d = 2 * std::abs(distB - outer_rad) / thickness;

		if (d > 1)
		{
			weight(x,y) *= maskedVal;
		}
		else
		{
			const double a = 1 - d;
			weight(x,y) *= maskedVal + (1 - maskedVal) * a * a;
		}
	}
}

double BlobFit2D::f(const std::vector<double>& x, void* tempStorage) const
{
	Blob2D blob(x, outer_radius);

	double out = 0.0;

	std::vector<double> radAvg = blob.radialAverage(image, weight);

	out += blob.radialAverageError(image, weight, radAvg);

	if (priorSigma2 > 0.0)
	{

		out += (blob.center - initialPos).norm2() / priorSigma2;
	}
	
	return out;
}

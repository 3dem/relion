#include "blob_fit_3d.h"
#include <src/jaz/tomography/fiducials.h>

using namespace gravis;


BlobFit3D::BlobFit3D(
	const Tomogram &tomogram,
	d3Vector position,
	int sh_bands,
	double radius,
	double thickness,
	const std::vector<d4Vector> &allSpheres,
	const std::vector<d3Vector> &fiducials,
	double fiducialRadius,
	double priorSigma,
	int num_threads)
:
		tomogram(tomogram),
		initialPos(position),
		outer_radius(radius + thickness/2),
		sh_bands(sh_bands),
		num_threads(num_threads),
		priorSigma2(priorSigma*priorSigma)
{
	const int fc = tomogram.stack.zdim;
	const int w = tomogram.stack.xdim;
	const int h = tomogram.stack.ydim;

	const double maskedVal = 1e-6;

	weight.resize(w,h,fc);
	weight.fill(1.f);

	const double outer_rad_2 = (radius + thickness / 2)*(radius + thickness / 2);

	for (int f = 0; f < fc; f++)
	{
		const d4Vector blob_img = tomogram.projectionMatrices[f]
				* d4Vector(position);

		RawImage<float> weightSlice = weight.getSliceRef(f);

		Fiducials::drawMask(fiducials, tomogram.projectionMatrices[f], fiducialRadius, weightSlice, maskedVal);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			const double dxb = blob_img.x - x;
			const double dyb = blob_img.y - y;
			const double distB = sqrt(dxb * dxb + dyb * dyb);

			const double d = 2 * std::abs(distB - outer_rad_2) / thickness;

			if (d > 1)
			{
				weight(x,y,f) *= maskedVal;
			}
			else
			{
				const double a = 1 - d;
				weight(x,y,f) *= maskedVal + (1 - maskedVal) * a * a;
			}
		}
	}
}

double BlobFit3D::f(const std::vector<double>& x, void* tempStorage) const
{
	SphericalHarmonics* sh = (SphericalHarmonics*) tempStorage;
	Blob3D blob(x, outer_radius, sh);

	const int fc = tomogram.stack.zdim;
	
	int padding = 4069;
	std::vector<double> per_thread(padding * num_threads, 0.0);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();
		
		std::vector<double> radAvg = blob.radialAverage(
				tomogram.stack.getConstSliceRef(f),
				tomogram.projectionMatrices[f],
				weight);

		per_thread[padding*t] += blob.radialAverageError(
				tomogram.stack.getConstSliceRef(f),
				tomogram.projectionMatrices[f],
				weight,
				radAvg);
	}
	
	double out = 0.0;
	
	for (int t = 0; t < num_threads; t++)
	{
		out += per_thread[padding*t];
	}
	
	out += (blob.center - initialPos).norm2() / priorSigma2;
	
	return out;
}

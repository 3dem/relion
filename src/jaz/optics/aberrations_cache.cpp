#include "aberrations_cache.h"
#include <src/jaz/math/Zernike.h>

using namespace gravis;


AberrationsCache::AberrationsCache(const MetaDataTable &opticsTable, const int size, const double pixelSize)
:	size(size), pixelSize(pixelSize)
{
	const int s = size;
	const int sh = s/2 + 1;
	const int gc = opticsTable.numberOfObjects();

	const d2Matrix M(1.0, 0.0, 0.0, 1.0);

	hasSymmetrical = opticsTable.containsLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS);
	hasAntisymmetrical = opticsTable.containsLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS);

	if (hasSymmetrical)
	{
		symmetrical.resize(gc);
	}

	if (hasAntisymmetrical)
	{
		antisymmetrical.resize(gc);
		phaseShift.resize(gc);
	}

	for (int g = 0; g < gc; g++)
	{
		const double as = pixelSize * s;

		if (hasSymmetrical)
		{
			std::vector<double> coeffs = opticsTable.getDoubleVector(
					EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, g);

			symmetrical[g] = BufferedImage<double>(sh,s);

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double phase = 0.0;

				for (int i = 0; i < coeffs.size(); i++)
				{
					int m, n;
					Zernike::evenIndexToMN(i, m, n);

					const double xx0 = x/as;
					const double yy0 = y < s/2? y/as : (y-s)/as;

					const double xx = M(0,0) * xx0 + M(0,1) * yy0;
					const double yy = M(1,0) * xx0 + M(1,1) * yy0;

					phase += coeffs[i] * Zernike::Z_cart(m,n,xx,yy);
				}

				symmetrical[g](x,y) = phase;
			}
		}

		if (hasAntisymmetrical)
		{
			std::vector<double> coeffs = opticsTable.getDoubleVector(
					EMDL_IMAGE_ODD_ZERNIKE_COEFFS, g);

			antisymmetrical[g] = BufferedImage<double>(sh,s);
			phaseShift[g] = BufferedImage<fComplex>(sh,s);

			for (int y = 0; y < s;  y++)
			for (int x = 0; x < sh; x++)
			{
				double phase = 0.0;

				for (int i = 0; i < coeffs.size(); i++)
				{
					int m, n;
					Zernike::oddIndexToMN(i, m, n);

					const double xx0 = x/as;
					const double yy0 = y < s/2? y/as : (y-s)/as;

					const double xx = M(0,0) * xx0 + M(0,1) * yy0;
					const double yy = M(1,0) * xx0 + M(1,1) * yy0;

					phase += coeffs[i] * Zernike::Z_cart(m,n,xx,yy);
				}

				antisymmetrical[g](x,y) = phase;
				phaseShift[g](x,y) = fComplex(cos(phase), sin(phase));
			}
		}
	}
}

void AberrationsCache::correctObservations(RawImage<fComplex> &observations, int optics_group) const
{
	if (!hasAntisymmetrical) return;

	const int sh = observations.xdim;
	const int s  = observations.ydim;
	const int fc = observations.zdim;

	if (phaseShift[optics_group].ydim != s)
	{
		REPORT_ERROR_STR(
			"reconstruct_particle: wrong cached phase-shift size. Box size: "
			<< s << ", cache size: " << phaseShift[optics_group].ydim);
	}

	for (int f  = 0; f  < fc; f++)
	for (int yi = 0; yi < s;  yi++)
	for (int xi = 0; xi < sh; xi++)
	{
		const fComplex r = phaseShift[optics_group](xi,yi);
		const fComplex z = observations(xi,yi,f);

		observations(xi,yi,f).real = z.real * r.real + z.imag * r.imag;
		observations(xi,yi,f).imag = z.imag * r.real - z.real * r.imag;
	}
}

#include "aberrations_cache.h"
#include <src/jaz/math/Zernike.h>

using namespace gravis;


AberrationsCache::AberrationsCache(const MetaDataTable &opticsTable, const int size)
:	size(size)
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
		const double pixelSize = opticsTable.getDouble(EMDL_IMAGE_PIXEL_SIZE, g);
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

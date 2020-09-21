#include "tilt_space_blob_fit.h"
#include <src/spherical-harmonics/SphericalHarmonics.h>

using namespace gravis;
using namespace au::edu::anu::qm::ro;


TiltSpaceBlobFit::TiltSpaceBlobFit(
	int sh_bands, 
	const RawImage<float>& correlation, 
	const RawImage<d3Vector>& directions_xz)
	:	
	sh_bands(sh_bands),
	correlation(correlation)
{
	SphericalHarmonics* SH = new SphericalHarmonics(sh_bands);
	
	const int w  = correlation.xdim;
	const int h  = correlation.ydim;
	const int fc = correlation.zdim;
	const int SH_params = sh_bands * sh_bands;
	
	basis.resize(SH_params, w, fc);
	
	for (int f = 0; f < fc; f++)
	for (int x = 0; x < w; x++)
	{
		const d3Vector v = directions_xz(x,0,f);
		const double phi = atan2(v.y, v.x);
		
		SH->computeY(sh_bands, v.z, phi, &basis(0,x,f));
	}
}

double TiltSpaceBlobFit::f(const std::vector<double> &x, void *tempStorage) const
{
	
}

void TiltSpaceBlobFit::grad(const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
	
}

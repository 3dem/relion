#include "tomogram.h"
#include "projection_IO.h"
#include "tomo_ctf_helper.h"
#include <src/jaz/optics/damage.h>
#include <src/jaz/image/power_spectrum.h>


using namespace gravis;


Tomogram::Tomogram()
{
	
}

double Tomogram::getFrameDose() const
{
	return cumulativeDose[frameSequence[1]] - cumulativeDose[frameSequence[0]];
}

BufferedImage<float> Tomogram::computeDoseWeight(int boxSize, double binning) const
{
	// @TODO: add support for B/k factors
	
	return Damage::weightStack_GG(cumulativeDose, optics.pixelSize * binning, boxSize);
}

BufferedImage<float> Tomogram::computeNoiseWeight(int boxSize, double binning, double overlap) const
{
	const int s0 = (int)(binning * boxSize + 0.5);
	const int s = boxSize;
	const int sh = s/2 + 1;
	const int fc = stack.zdim;

	BufferedImage<float> out(sh, s, fc);

	for (int f = 0; f < fc; f++)
	{
		BufferedImage<double> powSpec = PowerSpectrum::periodogramAverage2D(
			stack, s0, s0, overlap, f, false);

		std::vector<double> powSpec1D = RadialAvg::fftwHalf_2D_lin(powSpec);

		std::vector<float> frqWghts1D(powSpec1D.size());

		for (int i = 0; i < powSpec1D.size(); i++)
		{
			if (powSpec1D[i] > 0.0)
			{
				frqWghts1D[i] = (float)(1.0 / sqrt(powSpec1D[i]));
			}
			else
			{
				frqWghts1D[i] = 0.f;
			}
		}

		RawImage<float> outSlice = out.getSliceRef(f);
		RadialAvg::toFftwHalf_2D_lin(frqWghts1D, sh, s, outSlice, binning);
	}

	return out;
}

CTF Tomogram::getCtf(int frame, d3Vector position, double zOffset) const
{
	const d4Matrix& projFrame = projectionMatrices[frame];
	d4Vector pos2D = projFrame * d4Vector(position);
	d4Vector cent2D = projFrame * d4Vector(centre);
	
	double dz_pos = pos2D.z - cent2D.z;
	double dz = handedness * optics.pixelSize * dz_pos + zOffset;
	
	CTF ctf = centralCTFs[frame];
	
	ctf.DeltafU += dz;
	ctf.DeltafV += dz;
	
	ctf.initialise();
	
	return ctf;
}

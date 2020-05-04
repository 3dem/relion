#include "tomogram.h"
#include "projection_IO.h"
#include "tomo_ctf_helper.h"
#include <src/jaz/optics/damage.h>


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

#ifndef DELOCALISATION_HELPER_H
#define DELOCALISATION_HELPER_H

#include <src/ctf.h>

class DelocalisationHelper
{
	public:
		
		static void maskOutsideBox(
				const CTF& ctf, double radius, 
				double angpix, int s_orig, 
				MultidimArray<RFLOAT>& fftwCtfImg,
				double offsetx, double offsety);
		
		static Image<RFLOAT> plotDelocalisation(
				const CTF& ctf, Image<RFLOAT>& mask,
				double angpix);
};

#endif

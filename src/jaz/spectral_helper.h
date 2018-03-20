#ifndef SPECTRAL_HELPER_H
#define SPECTRAL_HELPER_H

#include <src/image.h>

class SpectralHelper
{
    public:

        static void computePhase(const Image<Complex>& src, Image<RFLOAT>& dest);
        static void computeAbs(const Image<Complex>& src, Image<RFLOAT>& dest);
};

#endif

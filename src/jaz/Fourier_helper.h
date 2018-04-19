#ifndef FOURIER_HELPER_H
#define FOURIER_HELPER_H

#include <src/image.h>
#include <src/complex.h>

class FourierHelper
{
    public:

        static void FourierShift2D(MultidimArray<Complex>& img, RFLOAT xshift, RFLOAT yshift);
        static void FourierShift2D(MultidimArray<RFLOAT>& img, RFLOAT xshift, RFLOAT yshift);
};

#endif

#ifndef SVD_HELPER_H
#define SVD_HELPER_H

#include <src/matrix2d.h>

class SvdHelper
{
    public:

        static void decompose(
            const Matrix2D<RFLOAT>& A,
            Matrix2D<RFLOAT>& U,
            Matrix1D<RFLOAT>& S,
            Matrix2D<RFLOAT>& Vt);
};

#endif

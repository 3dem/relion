#include <src/jaz/spectral_helper.h>

void SpectralHelper :: computePhase(const Image<Complex>& src, Image<RFLOAT>& dest)
{
    dest = Image<RFLOAT>(src.data.xdim, src.data.ydim, src.data.zdim, src.data.ndim);

    FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(src.data)
    {
        Complex z = DIRECT_NZYX_ELEM(src.data, l, k, i, j);
        DIRECT_NZYX_ELEM(dest.data, l, k, i, j) = atan2(z.imag, z.real);
    }
}

void SpectralHelper::computeAbs(const Image<Complex>& src, Image<RFLOAT>& dest)
{
    dest = Image<RFLOAT>(src.data.xdim, src.data.ydim, src.data.zdim, src.data.ndim);

    FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(src.data)
    {
        Complex z = DIRECT_NZYX_ELEM(src.data, l, k, i, j);
        DIRECT_NZYX_ELEM(dest.data, l, k, i, j) = z.abs();
    }
}

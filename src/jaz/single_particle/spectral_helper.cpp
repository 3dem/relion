/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include <src/jaz/single_particle/spectral_helper.h>

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

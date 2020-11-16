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

#include <src/jaz/single_particle/fftw_helper.h>

void FftwHelper::decenterUnflip2D(const MultidimArray<RFLOAT> &src, MultidimArray<RFLOAT> &dest)
{
    const long int w = src.xdim;
    const long int h = src.ydim;

    if (h != 2 * (w - 1))
    {
        std::cerr << "w = " << w << " h = " << h << std::endl;
        REPORT_ERROR("FftwHelper::decenterUnflip2D: illegal input size");
    }
    dest.reshape(h, h);

    const long int origin = w - 1;
    const long int nyquist = origin;

    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        // Logical coordinates of destination
        long int xl = x - origin;
        long int yl = y - origin;
        int sign = 1;

        if (xl < 0)
        {
            sign = -1;
            xl = -xl;
            yl = -yl;
        }

	// Cannot trust Nyquist
        if (xl == nyquist)
            xl--;
        if (yl == nyquist)
            yl--;
        if (xl == -nyquist)
            yl++;

        DIRECT_A2D_ELEM(dest, y, x) = sign * DIRECT_A2D_ELEM(src, (yl + h) % h, xl); 
    }
}

void FftwHelper::decenterDouble2D(const MultidimArray<RFLOAT> &src, MultidimArray<RFLOAT> &dest)
{
    const long int w = src.xdim;
    const long int h = src.ydim;

    if (h != 2 * (w - 1))
    {
        std::cerr << "w = " << w << " h = " << h << std::endl;
        REPORT_ERROR("FftwHelper::decenterDouble2D: illegal input size");
    }
    dest.reshape(h, h);

    const long int origin = w - 1;
    const long int nyquist = origin;

    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        // Logical coordinates of destination
        long int xl = x - origin;
        long int yl = y - origin;

        if (xl < 0)
        {
            xl = -xl;
            yl = -yl;
        }

	// Cannot trust Nyquist
        if (xl == nyquist)
            xl--;
        if (yl == nyquist)
            yl--;
        if (xl == -nyquist)
            yl++;

        DIRECT_A2D_ELEM(dest, y, x) = DIRECT_A2D_ELEM(src, (yl + h) % h, xl); 
    }
}

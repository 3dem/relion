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

    dest.reshape(h, 2*(w - 1));

    const long int yc = dest.ydim/2;

    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        long int xs = x - w;

        if (xs < 0)
        {
            long int ys = (y + yc - 1) % dest.ydim;
            DIRECT_A2D_ELEM(dest, y, x) = -DIRECT_A2D_ELEM(src, h-ys-1, -xs-1);
        }
        else
        {
            long int ys = (y + yc) % dest.ydim;
            DIRECT_A2D_ELEM(dest, y, x) =  DIRECT_A2D_ELEM(src, ys, xs);
        }
    }
}

void FftwHelper::decenterDouble2D(const MultidimArray<RFLOAT> &src, MultidimArray<RFLOAT> &dest)
{
    const long int w = src.xdim;
    const long int h = src.ydim;

    dest.reshape(h, 2*(w - 1));

    const long int yc = dest.ydim/2;

    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        long int xs = x - w;

        if (xs < 0)
        {
            long int ys = (y + yc - 1) % dest.ydim;
            DIRECT_A2D_ELEM(dest, y, x) = DIRECT_A2D_ELEM(src, h-ys-1, -xs-1);
        }
        else
        {
            long int ys = (y + yc) % dest.ydim;
            DIRECT_A2D_ELEM(dest, y, x) =  DIRECT_A2D_ELEM(src, ys, xs+1);
        }
    }
}

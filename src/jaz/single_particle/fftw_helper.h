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

#ifndef FFTW_HELPER_H
#define FFTW_HELPER_H

#include <src/multidim_array.h>

class FftwHelper
{
    public:

        template <typename T>
        static void decenterHalf(const MultidimArray<T> &src, MultidimArray<T> &dest);

        template <typename T>
        static void recenterHalf(const MultidimArray<T> &src, MultidimArray<T> &dest);

        template <typename T>
        static void decenterFull(const MultidimArray<T> &src, MultidimArray<T> &dest);

        template <typename T>
        static void recenterFull(const MultidimArray<T> &src, MultidimArray<T> &dest);

	/*
 	    decenterUnflip2D and decenterDouble2D uses Jasenko's convention of
	    taking the X Nyquist at positive and the Y Nyquist negative.
	    For N = 6, the input is treated as kx = [0, 1, 2, 3] and ky = [0, 1, 2, -3, -2, -1].

	    The output obeys the XMIPP convention, so both X and Y are [-3, -2, -1, 0, 1, 2].
  
	    We have two problems:
	        - Outside the jaz folder, the Y Nyquist is positive. So ky = [0, 1, 2, 3, -2, -1].
	          See comments in fftw.h and ctf.h.
	        - For aberration functions, f(kx = 3) is not always the same as f(kx = -3), because it is not periodic.

	    Because the Nyquist component is unreliable due to these issues, values from neighbouring row/column
            are taken instead.
  	*/
        static void decenterUnflip2D(const MultidimArray<RFLOAT> &src, MultidimArray<RFLOAT> &dest);
        static void decenterDouble2D(const MultidimArray<RFLOAT> &src, MultidimArray<RFLOAT> &dest);
};

template <typename T>
void FftwHelper::decenterHalf(const MultidimArray<T> &src, MultidimArray<T> &dest)
{
    dest.reshape(src);

    for (long int z = 0; z < dest.zdim; z++)
    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        long int zp = z < dest.xdim? z : z - dest.zdim;
        long int yp = y < dest.xdim? y : y - dest.ydim;

        DIRECT_A3D_ELEM(dest, z, y, x) = DIRECT_A3D_ELEM(src, zp, yp, x);
    }
}

template <typename T>
void FftwHelper::recenterHalf(const MultidimArray<T> &src, MultidimArray<T> &dest)
{
    dest.reshape(src);

    const long int zc = dest.zdim - dest.zdim/2;
    const long int yc = dest.ydim - dest.ydim/2;

    for (long int z = 0; z < dest.zdim; z++)
    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        long int zs = (z + zc) % dest.zdim;
        long int ys = (y + yc) % dest.ydim;

        DIRECT_A3D_ELEM(dest, z, y, x) = DIRECT_A3D_ELEM(src, zs, ys, x);
    }
}

template <typename T>
void FftwHelper::decenterFull(const MultidimArray<T> &src, MultidimArray<T> &dest)
{
    dest.reshape(src);

    const long int zc = dest.zdim/2;
    const long int yc = dest.ydim/2;
    const long int xc = dest.xdim/2;

    for (long int z = 0; z < dest.zdim; z++)
    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        long int zs = (z + zc) % dest.zdim;
        long int ys = (y + yc) % dest.ydim;
        long int xs = (x + xc) % dest.xdim;

        DIRECT_A3D_ELEM(dest, z, y, x) = DIRECT_A3D_ELEM(src, zs, ys, xs);
    }
}

template <typename T>
void FftwHelper::recenterFull(const MultidimArray<T> &src, MultidimArray<T> &dest)
{
    dest.reshape(src);

    const long int zc = dest.zdim - dest.zdim/2;
    const long int yc = dest.ydim - dest.ydim/2;
    const long int xc = dest.xdim - dest.xdim/2;

    for (long int z = 0; z < dest.zdim; z++)
    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        long int zs = (z + zc) % dest.zdim;
        long int ys = (y + yc) % dest.ydim;
        long int xs = (x + xc) % dest.xdim;

        DIRECT_A3D_ELEM(dest, z, y, x) = DIRECT_A3D_ELEM(src, zs, ys, xs);
    }
}

#endif

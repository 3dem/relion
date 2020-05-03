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

#ifndef RESAMPLING_HELPER_H
#define RESAMPLING_HELPER_H

#include <src/multidim_array.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/interpolation.h>

class ResamplingHelper
{
    public:

        // low-pass filter in real-space, then subsample by factor n
        template <typename T>
        static void downsampleGauss2D(const Image<T>& src, double n, Image<T>& dest);

        template <typename T>
        static void downsampleBox2D(const Image<T>& src, double n, Image<T>& dest);

        template <typename T>
        static void subsample2D(const Image<T>& src, double n, Image<T>& dest);

        template <typename T>
        static void subsample2D_cubic(const Image<T>& src, double n, Image<T>& dest, bool wrap = false);

        template <typename T>
        static void subsample3D(const Image<T>& src, int n, Image<T>& dest);

        template <typename T>
        static void upsample2D_linear(const Image<T>& src, int n, Image<T>& dest, bool wrap = false);

        template <typename T>
        static void upsample2D_cubic(const Image<T>& src, int n, Image<T>& dest,
                                 bool wrap = false, int w = -1, int h = -1);
};

template <typename T>
void ResamplingHelper::downsampleGauss2D(const Image<T>& src, double n, Image<T>& dest)
{
    Image<T> temp(src.data.ydim, src.data.xdim);
    FilterHelper::separableGaussianXY(src, temp, 0.5*(n-1), n-1);
    subsample2D(temp, n, dest);
}

template <typename T>
void ResamplingHelper::downsampleBox2D(const Image<T>& src, double n, Image<T>& dest)
{
    const int w0 = src.data.xdim;
    const int h0 = src.data.ydim;

    const int w1 = (int)(w0/n);
    const int h1 = (int)(h0/n);

    if (dest.data.xdim != w1 || dest.data.ydim != h1
        || dest.data.zdim != 1 || dest.data.ndim != 1)
    {
        dest = Image<RFLOAT>(w1,h1);
    }

    for (int y = 0; y < h1; y++)
    for (int x = 0; x < w1; x++)
    {
        T val = 0.0;
        RFLOAT wgh = 0.0;

        for (int yy = 0; yy < n; yy++)
        for (int xx = 0; xx < n; xx++)
        {
            int xin = (int)(x*n) + xx;
            int yin = (int)(y*n) + yy;

            if (xin < w0 && yin < h0)
            {
                val += DIRECT_NZYX_ELEM(src(), 0, 0, yin, xin);
                wgh += 1.0;
            }
        }

        if (wgh > 0.0)
        {
            val /= wgh;
        }

        DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = val;
    }
}

template <typename T>
void ResamplingHelper::subsample2D(const Image<T>& src, double n, Image<T>& dest)
{
    if (dest.data.xdim != (int)(src.data.xdim/n)
        || dest.data.ydim != (int)(src.data.ydim/n))
    {
        dest = Image<T>((int)(src.data.xdim/n), (int)(src.data.ydim/n));
    }

    for (int y = 0; y < dest.data.ydim; y++)
    for (int x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x)
                = DIRECT_NZYX_ELEM(src.data, 0, 0, (int)(n*y), (int)(n*x));
    }
}

template <typename T>
void ResamplingHelper::subsample2D_cubic(const Image<T>& src, double n, Image<T>& dest, bool wrap)
{
    dest.data.reshape(src.data.zdim, src.data.ydim/n, src.data.xdim/n);

    for (size_t z = 0; z < dest.data.zdim; z++)
    for (size_t y = 0; y < dest.data.ydim; y++)
    for (size_t x = 0; x < dest.data.xdim; x++)
    {
        double xx = x * (double)n;
        double yy = y * (double)n;

        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = Interpolation::cubicXY(src, xx, yy, z, 0, wrap);
    }
}

template <typename T>
void ResamplingHelper::subsample3D(const Image<T>& src, int n, Image<T>& dest)
{
    dest.data.reshape(src.data.zdim/n, src.data.ydim/n, src.data.xdim/n);

    for (size_t z = 0; z < dest.data.zdim; z++)
    for (size_t y = 0; y < dest.data.ydim; y++)
    for (size_t x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(src.data, 0, n*z, n*y, n*x);
    }
}

template <typename T>
void ResamplingHelper::upsample2D_linear(const Image<T>& src, int n, Image<T>& dest, bool wrap)
{
    dest.data.reshape(src.data.zdim, src.data.ydim*n, src.data.xdim*n);

    for (size_t z = 0; z < dest.data.zdim; z++)
    for (size_t y = 0; y < dest.data.ydim; y++)
    for (size_t x = 0; x < dest.data.xdim; x++)
    {
        int x0 = x/n;
        int y0 = y/n;
        int x1 = x0 + 1;
        int y1 = y0 + 1;

        double xf = (x / (double)n) - x0;
        double yf = (y / (double)n) - y0;

        if (wrap)
        {
            x1 = (x1 + src.data.xdim) % src.data.xdim;
            y1 = (y1 + src.data.ydim) % src.data.ydim;
        }
        else
        {
            if (x1 >= src.data.xdim) x1 = src.data.xdim - 1;
            if (y1 >= src.data.ydim) y1 = src.data.ydim - 1;
        }

        T v00 = DIRECT_NZYX_ELEM(src.data, 0, z, y0, x0);
        T v01 = DIRECT_NZYX_ELEM(src.data, 0, z, y0, x1);
        T v10 = DIRECT_NZYX_ELEM(src.data, 0, z, y1, x0);
        T v11 = DIRECT_NZYX_ELEM(src.data, 0, z, y1, x1);

        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) =
                xf*(yf*v11 + (1-yf)*v01) + (1-xf)*(yf*v10 + (1-yf)*v00);
    }
}

template <typename T>
void ResamplingHelper::upsample2D_cubic(const Image<T>& src, int n, Image<T>& dest, bool wrap, int w, int h)
{
    if (w < 0) w = src.data.xdim*n;
    if (h < 0) h = src.data.ydim*n;

    dest.data.reshape(src.data.zdim, h, w);

    for (size_t z = 0; z < dest.data.zdim; z++)
    for (size_t y = 0; y < h; y++)
    for (size_t x = 0; x < w; x++)
    {
        double xx = x / (double)n;
        double yy = y / (double)n;

        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = Interpolation::cubicXY(src, xx, yy, z, 0, wrap);
    }
}

#endif

#ifndef RESAMPLING_HELPER_H
#define RESAMPLING_HELPER_H

#include <src/multidim_array.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/interpolation.h>

class ResamplingHelper
{
    public:

    // low-pass filter in real-space, then subsample by factor n
    template <typename T>
    static void downsample(const Image<T>& src, int n, Image<T>& dest);

    template <typename T>
    static void subsample2D(const Image<T>& src, int n, Image<T>& dest);

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
void ResamplingHelper::downsample(const Image<T>& src, int n, Image<T>& dest)
{
    Image<T> temp(src.data.zdim, src.data.ydim, src.data.xdim);
    FilterHelper::separableGaussian(src.data, temp.data, 0.5*(n-1), n-1);
    subsample(temp, n, dest);
}

template <typename T>
void ResamplingHelper::subsample2D(const Image<T>& src, int n, Image<T>& dest)
{
    dest.data.reshape(src.data.zdim, src.data.ydim/n, src.data.xdim/n);

    for (size_t z = 0; z < dest.data.zdim; z++)
    for (size_t y = 0; y < dest.data.ydim; y++)
    for (size_t x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(src.data, 0, z, n*y, n*x);
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

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

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

#ifndef IMAGE_OPS_H
#define IMAGE_OPS_H

#include <src/image.h>

class ImageOp
{
    public:

    template<typename T1, typename T2>
    static void linearCombination(const Image<T1>& src0, const Image<T1>& src1, T2 a0, T2 a1, Image<T1>& dest);

    template<typename T1, typename T2, typename T3>
    static void multiply(const Image<T1>& i0, const Image<T2>& i1, Image<T3>& dest);

    template<typename T1, typename T2>
    static void multiplyBy(Image<T1>& dest, const Image<T2>& i1);

    template<typename T1, typename T2>
    static void linearCombination(const Image<T1>& src0, T1 src1, T2 a0, T2 a1, Image<T1>& dest);


    template<typename T1, typename T2>
    static void linearCombination(const MultidimArray<T1>& src0, const MultidimArray<T1>& src1, T2 a0, T2 a1, MultidimArray<T1>& dest);

    template<typename T1, typename T2, typename T3>
    static void multiply(const MultidimArray<T1>& i0, const MultidimArray<T2>& i1, MultidimArray<T3>& dest);

    template<typename T1, typename T2>
    static void linearCombination(const MultidimArray<T1>& src0, T1 src1, T2 a0, T2 a1, MultidimArray<T1>& dest);

    template<typename T1>
    static void flipX(const MultidimArray<T1>& src0, MultidimArray<T1>& dest);

    template<typename T1>
    static void flipY(const MultidimArray<T1>& src0, MultidimArray<T1>& dest);

    template<typename T1>
    static void rotate90(const MultidimArray<T1>& src0, MultidimArray<T1>& dest);

    template<typename T1>
    static void rotate180(const MultidimArray<T1>& src0, MultidimArray<T1>& dest);

    template<typename T1>
    static void rotate270(const MultidimArray<T1>& src0, MultidimArray<T1>& dest);
};

template<typename T1, typename T2>
void ImageOp :: linearCombination(const Image<T1>& src0, const Image<T1>& src1, T2 a0, T2 a1, Image<T1>& dest)
{
    for (long int n = 0; n < src0.data.ndim; n++)
    for (long int z = 0; z < src0.data.zdim; z++)
    for (long int y = 0; y < src0.data.ydim; y++)
    for (long int x = 0; x < src0.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = a0 * DIRECT_NZYX_ELEM(src0.data, n, z, y, x) + a1 * DIRECT_NZYX_ELEM(src1.data, n, z, y, x);
    }
}

template<typename T1, typename T2, typename T3>
void ImageOp :: multiply(const Image<T1>& i0, const Image<T2>& i1, Image<T3>& dest)
{
    dest = Image<T3>(i0.data.xdim, i0.data.ydim, i0.data.zdim, i0.data.ndim);

    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = 0; y < i0.data.ydim; y++)
    for (long int x = 0; x < i0.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = DIRECT_NZYX_ELEM(i0.data, n, z, y, x) * DIRECT_NZYX_ELEM(i1.data, n, z, y, x);
    }
}

template<typename T1, typename T2>
void ImageOp :: multiplyBy(Image<T1>& dest, const Image<T2>& i1)
{
    for (long int n = 0; n < dest.data.ndim; n++)
    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) *= DIRECT_NZYX_ELEM(i1.data, n, z, y, x);
    }
}

template<typename T1, typename T2>
void ImageOp :: linearCombination(const Image<T1>& src0, T1 src1, T2 a0, T2 a1, Image<T1>& dest)
{
    for (long int n = 0; n < src0.data.ndim; n++)
    for (long int z = 0; z < src0.data.zdim; z++)
    for (long int y = 0; y < src0.data.ydim; y++)
    for (long int x = 0; x < src0.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = a0 * DIRECT_NZYX_ELEM(src0.data, n, z, y, x) + a1 * src1;
    }
}


template<typename T1, typename T2>
void ImageOp :: linearCombination(const MultidimArray<T1>& src0, const MultidimArray<T1>& src1, T2 a0, T2 a1, MultidimArray<T1>& dest)
{
    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest, n, z, y, x) = a0 * DIRECT_NZYX_ELEM(src0, n, z, y, x) + a1 * DIRECT_NZYX_ELEM(src1, n, z, y, x);
    }
}

template<typename T1, typename T2, typename T3>
void ImageOp :: multiply(const MultidimArray<T1>& i0, const MultidimArray<T2>& i1, MultidimArray<T3>& dest)
{
    dest = MultidimArray<T3>(i0.xdim, i0.ydim, i0.zdim, i0.ndim);

    for (long int n = 0; n < i0.ndim; n++)
    for (long int z = 0; z < i0.zdim; z++)
    for (long int y = 0; y < i0.ydim; y++)
    for (long int x = 0; x < i0.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest, n, z, y, x) = DIRECT_NZYX_ELEM(i0, n, z, y, x) * DIRECT_NZYX_ELEM(i1, n, z, y, x);
    }
}

template<typename T1, typename T2>
void ImageOp :: linearCombination(const MultidimArray<T1>& src0, T1 src1, T2 a0, T2 a1, MultidimArray<T1>& dest)
{
    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest, n, z, y, x) = a0 * DIRECT_NZYX_ELEM(src0, n, z, y, x) + a1 * src1;
    }
}

// flip 'left-to-right'
// MotionCor2's FlipGain 2
template<typename T1>
void ImageOp :: flipX(const MultidimArray<T1>& src0, MultidimArray<T1>& dest)
{
    dest.reshape(src0);

    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest, n, z, y, x) = DIRECT_NZYX_ELEM(src0, n, z, y, src0.xdim - 1 - x);
    }
}

// flip 'upside down'
// MotionCor2's FlipGain 1
template<typename T1>
void ImageOp :: flipY(const MultidimArray<T1>& src0, MultidimArray<T1>& dest)
{
    dest.reshape(src0);

    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest, n, z, y, x) = DIRECT_NZYX_ELEM(src0, n, z, src0.ydim - 1 - y, x);
    }
}

// This is equivalent to MotionCor2's -RotGain 1.
// In relion_display, this looks clock-wise.
template<typename T1>
void ImageOp :: rotate90(const MultidimArray<T1>& src0, MultidimArray<T1>& dest)
{
    dest.reshape(src0.ndim, src0.zdim, src0.xdim, src0.ydim);

    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest, n, z, x, src0.ydim - 1 - y) = DIRECT_NZYX_ELEM(src0, n, z, y, x);
    }
}

// MotionCor2's RotGain 2
template<typename T1>
void ImageOp :: rotate180(const MultidimArray<T1>& src0, MultidimArray<T1>& dest)
{
    dest.reshape(src0);

    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest, n, z, src0.ydim - 1 - y, src0.xdim - 1 - x) = DIRECT_NZYX_ELEM(src0, n, z, y, x);
    }
}

// MotionCor2's RotGain 3
template<typename T1>
void ImageOp :: rotate270(const MultidimArray<T1>& src0, MultidimArray<T1>& dest)
{
    dest.reshape(src0.ndim, src0.zdim, src0.xdim, src0.ydim);

    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest, n, z, src0.xdim - 1 - x, y) = DIRECT_NZYX_ELEM(src0, n, z, y, x);
    }
}

#endif

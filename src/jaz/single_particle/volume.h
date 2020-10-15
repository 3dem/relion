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

#ifndef JAZ_VOLUME_H
#define JAZ_VOLUME_H

#include <vector>
#include <stddef.h>
#include <src/jaz/single_particle/config.h>

/* class Volume: represents a grid of voxels of arbitrary type.
   Only provides methods for direct access at integral coordinates.
   Everything else is handled by external classes. */

#define FOR_ALL_VOXELS(V) \
    for (size_t z = 0; z < (V).dimz; z++) \
    for (size_t y = 0; y < (V).dimy; y++) \
    for (size_t x = 0; x < (V).dimx; x++)

template <typename T>
class Volume
{
    public:

        Volume(){}

        /* The constructor only allocates memory, it does not initialize the values.*/

        Volume(size_t dimx, size_t dimy, size_t dimz);

            long int dimx, dimy, dimz;
            std::vector<T> voxels;


        /* operator (x,y,z): returns a reference to the indicated voxel.
           The correct version (const or non-const) will be chosen by the compiler,
           depending on whether the instance is declared as const or not.*/

        const T& operator() (size_t, size_t, size_t) const;
        T& operator() (size_t, size_t, size_t);


        /* data(): returns a pointer to the first data element.*/

        const T* data() const;
        T* data();

        void resize(size_t dimx, size_t dimy, size_t dimz);
        void resize(const Volume& example);
        void fill(T t);
        
        Volume& operator += (const Volume& v)
        {
            for (int i = 0; i < voxels.size(); i++)
            {
                voxels[i] += v.voxels[i];
            }

            return *this;
        }
};

template <class T>
Volume<T>::Volume(size_t dimx, size_t dimy, size_t dimz)
    :   dimx(dimx), dimy(dimy), dimz(dimz)
{
    voxels.resize(dimx*dimy*dimz);
}

template <class T>
inline const T& Volume<T>::operator() (size_t x, size_t y, size_t z) const
{
    return voxels[(z*dimy + y)*dimx + x];
}

template <class T>
inline T& Volume<T>::operator() (size_t x, size_t y, size_t z)
{
    return voxels[(z*dimy + y)*dimx + x];
}

template <class T>
inline const T* Volume<T>::data() const
{
    return &voxels[0];
}

template <class T>
inline T* Volume<T>::data()
{
    return &voxels[0];
}

template <class T>
inline void Volume<T>::resize(size_t dimx, size_t dimy, size_t dimz)
{
    this->dimx = dimx;
    this->dimy = dimy;
    this->dimz = dimz;

    voxels.resize(dimx*dimy*dimz);
}

template <class T>
inline void Volume<T>::resize(const Volume& example)
{
    dimx = example.dimx;
    dimy = example.dimy;
    dimz = example.dimz;

    voxels.resize(dimx*dimy*dimz);
}

template <class T>
inline void Volume<T>::fill(T t)
{
    for (long int i = 0; i < voxels.size(); i++)
    {
        voxels[i] = t;
    }
}

#endif

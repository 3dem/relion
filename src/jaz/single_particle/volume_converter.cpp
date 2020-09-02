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

#include <src/jaz/single_particle/volume_converter.h>

void VolumeConverter::convert(const Image<RFLOAT>& src, Volume<RFLOAT>& dest)
{
    dest.resize(src.data.xdim, src.data.ydim, src.data.zdim);

    FOR_ALL_VOXELS(dest)
    {
        dest(x,y,z) = DIRECT_A3D_ELEM(src.data, z, y, x);
    }
}

void VolumeConverter::convertStack(const Image<RFLOAT>& src, Volume<RFLOAT>& dest)
{
    dest.resize(src.data.xdim, src.data.ydim, src.data.ndim);

    FOR_ALL_VOXELS(dest)
    {
        dest(x,y,z) = DIRECT_NZYX_ELEM(src.data, z, 1, y, x);
    }
}

void VolumeConverter::convert(const Volume<RFLOAT>& src, Image<RFLOAT>& dest)
{
    dest.data.resize(src.dimz, src.dimy, src.dimx);

    FOR_ALL_VOXELS(src)
    {
        DIRECT_A3D_ELEM(dest.data, z, y, x) = src(x,y,z);
    }
}


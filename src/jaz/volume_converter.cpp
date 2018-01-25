#include <src/jaz/volume_converter.h>

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


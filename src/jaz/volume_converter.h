#ifndef VOLUME_CONVERTER_H
#define VOLUME_CONVERTER_H

#include <src/jaz/volume.h>
#include <src/image.h>
#include <src/multidim_array.h>

/* class VolumeConverter: facilitates conversion between 'Volume' and Relion's 'Image' and 'MultidimArray' classes. */

class VolumeConverter
{
    public:

        static void convert(const Image<RFLOAT>& src, Volume<RFLOAT>& dest);
        static void convertStack(const Image<RFLOAT>& src, Volume<RFLOAT>& dest);
        static void convert(const Volume<RFLOAT>& src, Image<RFLOAT>& dest);
};

#endif

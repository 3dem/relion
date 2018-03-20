#ifndef VOLUME_INTEGRATION_H
#define VOLUME_INTEGRATION_H

#include "src/jaz/volume.h"
#include "src/jaz/gravis/t4Matrix.h"
#include "src/image.h"

class VolumeIntegration
{
    public:

    static void integrateAlongZ(const Volume<RFLOAT>& vol, gravis::d4Matrix vol2img, Image<RFLOAT>& dest);
};

#endif

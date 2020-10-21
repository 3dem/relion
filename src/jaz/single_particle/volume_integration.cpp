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

#include <src/jaz/single_particle/volume_integration.h>
#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/gravis/t2Matrix.h>

using namespace gravis;

void VolumeIntegration :: integrateAlongZ(const Volume<RFLOAT>& vol, gravis::d4Matrix vol2img, Image<RFLOAT>& dest)
{
    const int xsv = vol.dimx;
    const int ysv = vol.dimy;
    const int zsv = vol.dimz;

    const int xsi = dest.data.xdim;
    const int ysi = dest.data.ydim;

    d2Matrix A2;
    A2(0,0) = vol2img(0,0);
    A2(0,1) = vol2img(0,1);
    A2(1,0) = vol2img(1,0);
    A2(1,1) = vol2img(1,1);
    A2.invert();

    #if JAZ_USE_OPENMP
    #pragma omp parallel for
    #endif
    for (int yi = 0; yi < ysi; yi++)
    for (int xi = 0; xi < xsi; xi++)
    {
        DIRECT_A2D_ELEM(dest.data, yi, xi) = 0.0;
    }

    for (int zv = 0; zv < zsv; zv++)
    {
        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (int yi = 0; yi < ysi; yi++)
        for (int xi = 0; xi < xsi; xi++)
        {
            d2Vector b(xi - zv*vol2img(0,2) - vol2img(0,3), yi - zv*vol2img(1,2) - vol2img(1,3));
            d2Vector v = A2 * b;

            int xvi = (int) v.x;
            int yvi = (int) v.y;

            double xvf = v.x - xvi;
            double yvf = v.y - yvi;

            if (xvi < 0 || yvi < 0 || xvi >= xsv-1 || yvi >= ysv-1) continue;

            RFLOAT vv00 = vol(xvi,yvi,zv);
            RFLOAT vv10 = vol(xvi+1,yvi,zv);
            RFLOAT vv01 = vol(xvi,yvi+1,zv);
            RFLOAT vv11 = vol(xvi+1,yvi+1,zv);

            RFLOAT vv0 = yvf * vv01 + (1.0 - yvf) * vv00;
            RFLOAT vv1 = yvf * vv11 + (1.0 - yvf) * vv10;

            RFLOAT vv = xvf * vv1 + (1.0 - xvf) * vv0;

            DIRECT_A2D_ELEM(dest.data, yi, xi) += vv;
        }
    }
}

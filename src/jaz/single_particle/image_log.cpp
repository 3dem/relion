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

#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>


void ImageLog::write(
    Image<Complex> &img, std::string fn,
    bool polar, Centering center,
    double originX, double originY, double originZ,
    double spacingX, double spacingY, double spacingZ)
{
    if (polar)
    {
        Image<RFLOAT> argImg, absImg;
        FilterHelper::getPhase(img, argImg);
        FilterHelper::getAbs(img, absImg);

        write(argImg, fn+"_arg", center, originX, originY, originZ, spacingX, spacingY, spacingZ);
        write(absImg, fn+"_abs", center, originX, originY, originZ, spacingX, spacingY, spacingZ);
    }
    else
    {
        Image<RFLOAT> realImg, imagImg;
        FilterHelper::getReal(img, realImg);
        FilterHelper::getImag(img, imagImg);

        write(realImg, fn+"_re", center, originX, originY, originZ, spacingX, spacingY, spacingZ);
        write(imagImg, fn+"_im", center, originX, originY, originZ, spacingX, spacingY, spacingZ);
    }
}

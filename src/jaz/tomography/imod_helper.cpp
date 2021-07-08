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

#include "imod_helper.h"
#include <src/error.h>
#include <src/macros.h>
#include <fstream>

using namespace gravis;

std::vector<d4Matrix> ImodHelper::readTiltTransforms(std::string fn, d4Matrix vol2world, double cix, double ciy)
{
    std::ifstream anglesFile(fn.c_str());

    if (!anglesFile.is_open())
    {
        REPORT_ERROR("ImodHelper::readTiltTransforms: failed to open "+fn+".");
    }

    std::vector<d4Matrix> vol2img;

    const double deg2rad = PI/180.0;

    while (anglesFile.good())
    {
        double a;
        anglesFile >> a;
        a *= deg2rad;

        d4Matrix w2i;
        w2i(0,0) =  cos(a);
        w2i(0,2) =  sin(a);
        w2i(2,0) = -sin(a);
        w2i(2,2) =  cos(a);
        w2i(0,3) =  cix;
        w2i(1,3) =  ciy;

        vol2img.push_back(w2i * vol2world);
    }

    return vol2img;
}

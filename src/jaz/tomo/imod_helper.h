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

#ifndef IMOD_HELPER_H
#define IMOD_HELPER_H

#include <vector>
#include <string>
#include <src/jaz/gravis/t4Matrix.h>

class ImodHelper
{
    public:

    static std::vector<gravis::d4Matrix> readTiltTransforms(std::string fn, gravis::d4Matrix vol2world, double cix, double ciy);
    static std::vector<gravis::d4Matrix> readAffineTransforms(std::string fn);
};

#endif

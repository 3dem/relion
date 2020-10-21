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

#ifndef PROJECTION_HELPER_H
#define PROJECTION_HELPER_H

#include <vector>
#include <src/jaz/gravis/t4Matrix.h>

class ProjectionHelper
{
    public:

    /* loads a sequence of tilt angles and returns 4x4 matrices that map *world space* coordinates to image coordinates*/
    static std::vector<gravis::d4Matrix> loadTiltProjections(
            std::string tiltFile,               // file containing the tilt angles in ASCII
            double centerX, double centerY);    // world origin projected into the images (usually, the image center)

    /* loads a sequence of tilt angles and returns 4x4 matrices that map *voxel* coordinates to image coordinates*/
    static std::vector<gravis::d4Matrix> loadTiltProjectionsVol(
            std::string tiltFile,               // file containing the tilt angles in ASCII
            double centerX, double centerY,     // world origin projected into the images (usually, the image center)
            double X0, double Y0, double Z0,    // origin of the volume in world coordinates
            double spacing = 1.0);              // volume resolution

    /* loads a sequence of affine transforms*/
    static std::vector<gravis::d4Matrix> loadAffineTransforms(
            std::string xformFile,              // file containing the affine transforms in ASCII (ie. the .xf-file from imod)
            double cx, double cy,               // coordinates of image center
            bool square_result = true);         // coordinate origin in the output image is at (cx,cx), not (cx,cy)
};

#endif

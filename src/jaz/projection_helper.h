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

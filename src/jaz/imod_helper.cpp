#include <src/jaz/imod_helper.h>
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

std::vector<d4Matrix> ImodHelper::readAffineTransforms(std::string fn)
{
    /*std::ifstream mapFile(fn.c_str());

    if (!mapFile.is_open())
    {
        REPORT_ERROR("ImodHelper::readAffineTransforms: failed to open "+fn+".");
    }

    std::vector<double> angles;
    std::vector<d4Matrix> vol2img;

    while (anglesFile.good())
    {
        double a;
        anglesFile >> a;
        a *= deg2rad;

        angles.push_back(a);

        d4Matrix w2i;
        w2i(0,0) =  cos(a);
        w2i(0,2) =  sin(a);
        w2i(2,0) = -sin(a);
        w2i(2,2) =  cos(a);
        w2i(0,3) =  cix;
        w2i(1,3) =  ciy;

        vol2img.push_back(w2i * vol2world);
    }

    return vol2img;*/
}

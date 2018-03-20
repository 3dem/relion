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

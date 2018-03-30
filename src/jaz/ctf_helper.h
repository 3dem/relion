#ifndef CTF_HELPER_H
#define CTF_HELPER_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class CtfHelper
{
    public:

        static std::vector<CTF> loadCtffind4(std::string path, int imageCount,
                                             double voltage = 300.0, double Cs = 2.2,
                                             double Q0 = 0.1, double Bfac = 0.0,
                                             double scale = 1.0);
};

#endif

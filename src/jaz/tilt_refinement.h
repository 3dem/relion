#ifndef TILT_REFINEMENT_H
#define TILT_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class TiltOptimization : public Optimization
{
    public:

        TiltOptimization(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                double angpix,
                bool L1 = false,
                bool anisotropic = false);

        double f(const std::vector<double>& x) const;

    private:

        const Image<Complex>& xy;
        const Image<RFLOAT>& weight;
        const double angpix;
        const bool L1, anisotropic;
};

class TiltRefinement
{
    public:

        static void updateTiltShift(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                CTF& ctf, RFLOAT angpix,
                Image<Complex>& xyDest,
                Image<RFLOAT>& wDest);

        static void updateTiltShiftPar(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                CTF& ctf, RFLOAT angpix,
                Image<Complex>& xyDest,
                Image<RFLOAT>& wDest);

        static void fitTiltShift(
                const Image<RFLOAT>& phase,
                const Image<RFLOAT>& weight,
                RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
                RFLOAT* shift_x, RFLOAT* shift_y,
                RFLOAT* tilt_x, RFLOAT* tilt_y,
                Image<RFLOAT>* fit,
                gravis::d2Matrix magCorr = gravis::d2Matrix());

        static void optimizeTilt(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
                bool L1,
                RFLOAT shift0_x, RFLOAT shift0_y,
                RFLOAT tilt0_x, RFLOAT tilt0_y,
                RFLOAT* shift_x, RFLOAT* shift_y,
                RFLOAT* tilt_x, RFLOAT* tilt_y,
                Image<RFLOAT>* fit);

        static void drawPhaseShift(
                RFLOAT shift_x, RFLOAT shift_y,
                RFLOAT tilt_x, RFLOAT tilt_y,
                int w, int h, double as,
                gravis::d2Matrix magCorr,
                Image<RFLOAT>* tgt);

};

#endif

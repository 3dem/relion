#ifndef TILT_REFINEMENT_H
#define TILT_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class TiltRefinement
{
    public:

        static void updateTiltShift(const Image<Complex>& prediction,
                                 const Image<Complex>& observation,
                                 CTF& ctf, RFLOAT angpix,
                                 Image<Complex>& xyDest,
                                 Image<RFLOAT>& wDest);

        static void updateTiltShiftPar(const Image<Complex>& prediction,
                                 const Image<Complex>& observation,
                                 CTF& ctf, RFLOAT angpix,
                                 Image<Complex>& xyDest,
                                 Image<RFLOAT>& wDest);

        static void fitTiltShift(const Image<RFLOAT>& phase,
                                const Image<RFLOAT>& weight,
                                RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
                                RFLOAT* shift_x, RFLOAT* shift_y,
                                RFLOAT* tilt_x, RFLOAT* tilt_y,
                                Image<RFLOAT>* fit,
                                gravis::d2Matrix magCorr = gravis::d2Matrix());
};

#endif

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

#ifndef TILT_REFINEMENT_H
#define TILT_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
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

        double f(const std::vector<double>& x, void* tempStorage) const;

    private:

        const Image<Complex>& xy;
        const Image<RFLOAT>& weight;
        const double angpix;
        const bool L1, anisotropic;
};

class TiltHelper
{
    public:

        static void updateTiltShift(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                const CTF& ctf, double angpix,
                Image<Complex>& xyDest,
                Image<RFLOAT>& wDest);

        static void updateTiltShiftPar(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                const CTF& ctf, double angpix,
                Image<Complex>& xyDest,
                Image<RFLOAT>& wDest);

        static void fitTiltShift(
                const Image<RFLOAT>& phase,
                const Image<RFLOAT>& weight,
                double Cs, double lambda, double angpix,
                double* shift_x, double* shift_y,
                double* tilt_x, double* tilt_y,
                Image<RFLOAT>* fit,
                gravis::d2Matrix magCorr = gravis::d2Matrix());

        static void optimizeTilt(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                double Cs, double lambda, double angpix,
                bool L1,
                double shift0_x, double shift0_y,
                double tilt0_x, double tilt0_y,
                double* shift_x, double* shift_y,
                double* tilt_x, double* tilt_y,
                Image<RFLOAT>* fit);

        static void optimizeAnisoTilt(
                const Image<Complex>& xy,
                const Image<RFLOAT>& weight,
                double Cs, double lambda, double angpix,
                bool L1,
                double shift0_x, double shift0_y,
                double tilt0_x, double tilt0_y,
                double* shift_x, double* shift_y,
                double* tilt_x, double* tilt_y,
                double* tilt_xx, double* tilt_xy, double* tilt_yy,
                Image<RFLOAT>* fit);

        static void drawPhaseShift(
                double shift_x, double shift_y,
                double tilt_x, double tilt_y,
                int w, int h, double as,
                gravis::d2Matrix magCorr,
                Image<RFLOAT>* tgt);

        static void drawPhaseShift(
                double shift_x, double shift_y,
                double tilt_x, double tilt_y,
                double tilt_xx, double tilt_xy, double tilt_yy,
                int w, int h, double as,
                gravis::d2Matrix magCorr,
                Image<RFLOAT>* tgt);

};

#endif

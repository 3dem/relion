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

#ifndef DEFOCUS_REFINEMENT_H
#define DEFOCUS_REFINEMENT_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/single_particle/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <vector>

class AstigmatismOptimization : public Optimization
{
    public:

        AstigmatismOptimization(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                const Image<RFLOAT>& weight,
                const CTF& ctf0,
                RFLOAT angpix);

        double f(const std::vector<double>& x) const;

    private:

        const Image<Complex>& prediction, observation;
        const Image<RFLOAT>& weight;
        const CTF& ctf0;
        RFLOAT angpix;
};

class AstigmatismOptimizationAcc : public Optimization
{
    public:

        AstigmatismOptimizationAcc(
                const Image<Complex>& prediction,
                const Image<Complex>& observation,
                const Image<RFLOAT>& weight,
                const CTF& ctf0,
                bool phaseShift,
                bool spherAberr,
                RFLOAT angpix,
                RFLOAT phiScale = 10,
                RFLOAT csScale = 10);

        AstigmatismOptimizationAcc(
                const std::vector<Image<Complex>>& prediction,
                const std::vector<Image<Complex>>& observation,
                const Image<RFLOAT>& weight,
                const CTF& ctf0,
                bool phaseShift,
                bool spherAberr,
                RFLOAT angpix,
                RFLOAT phiScale = 10,
                RFLOAT csScale = 100);

        double f(const std::vector<double>& x, void* tempStorage) const;

        double getU(const std::vector<double> &x);
        double getV(const std::vector<double> &x);
        double getPhi(const std::vector<double> &x);
        double getPhase(const std::vector<double> &x);
        double getCs(const std::vector<double> &x);

        std::vector<double> getInitialParams();

    private:

        Image<Complex> data;
        const CTF& ctf0;
        bool phaseShift, spherAberr;
        RFLOAT angpix, phiScale, csScale;
};

class DefocusHelper
{
    public:

        static RFLOAT findDefocus1D(const Image<Complex>& prediction,
            const Image<Complex>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, RFLOAT angpix,
            double* destU, double* destV,
            RFLOAT range = 1000.0, int steps = 11,
            int recDepth = 2, RFLOAT recScale = 10.0);

        static void findAstigmatismNM(
            const Image<Complex>& prediction,
            const Image<Complex>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, RFLOAT angpix,
            double* destU, double* destV, double* destPhi);

        static void findAstigmatismAndPhaseNM(
            const std::vector<Image<Complex>>& prediction,
            const std::vector<Image<Complex>>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, RFLOAT angpix,
            double* destU, double* destV, double* destPhi, double* destPhase);

        static void findAstigmatismPhaseAndCsNM(
            const std::vector<Image<Complex>>& prediction,
            const std::vector<Image<Complex>>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, RFLOAT angpix,
            double* destU, double* destV,
            double* destPhi, double* destPhase, double* destCs);

        static void findAstigmatismNM(
            const std::vector<Image<Complex>>& prediction,
            const std::vector<Image<Complex>>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, RFLOAT angpix,
            double* destU, double* destV, double* destPhi);

        static std::vector<gravis::d2Vector> diagnoseDefocus(const Image<Complex>& prediction,
            const Image<Complex>& observation,
            const Image<RFLOAT>& weight,
            const CTF& ctf0, RFLOAT angpix,
            double range, int steps, int threads);
};

#endif

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

#ifndef GP_MOTION_FIT
#define GP_MOTION_FIT

#include <src/image.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/gravis/t2Vector.h>
#include <vector>

class GpMotionFit : public DifferentiableOptimization
{
    public:

        GpMotionFit(
                const std::vector<std::vector<Image<double>>>& correlation,
                double cc_pad,
				double sig_vel_px, double sig_div_px, double sig_acc_px,
                int maxDims,
                const std::vector<gravis::d2Vector>& positions,
                const std::vector<gravis::d2Vector>& perFrameOffsets,
                int threads, bool expKer);

        double f(const std::vector<double>& x) const;
        double f(const std::vector<double>& x, void* tempStorage) const;

        void grad(const std::vector<double>& x, std::vector<double>& gradDest) const;
        void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;

        void* allocateTempStorage() const;
        void deallocateTempStorage(void* ts) const;

        void paramsToPos(const std::vector<double>& x,
                         std::vector<std::vector<gravis::d2Vector>>& pos) const;

        void posToParams(const std::vector<std::vector<gravis::d2Vector>>& pos,
                         std::vector<double>& x) const;

        class TempStorage
        {
            public:

                int pad;

                std::vector<std::vector<gravis::d2Vector>> pos, ccg_pf;
                std::vector<std::vector<double>> gradDestT;
                std::vector<double> e_t;
        };

    private:

        bool expKer;
        int pc, fc, dc, threads;
        double cc_pad, sig_vel_px, sig_div_px, sig_acc_px;

        Matrix2D<RFLOAT> basis;
        std::vector<double> eigenVals;

        const std::vector<std::vector<Image<double>>>& correlation;
        const std::vector<gravis::d2Vector>& positions;
        const std::vector<gravis::d2Vector>& perFrameOffsets;
};

#endif

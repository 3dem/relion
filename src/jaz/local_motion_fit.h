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

#ifndef LOCAL_MOTION_FIT
#define LOCAL_MOTION_FIT

#include <src/image.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/gravis/t2Vector.h>
#include <vector>

class LocalMotionFit : public DifferentiableOptimization
{
    public:

        LocalMotionFit(
                const std::vector<std::vector<Image<RFLOAT>>>& correlation,
                const std::vector<double>& velWgh,
                const std::vector<double>& accWgh,
                const std::vector<std::vector<std::vector<double>>>& divWgh,
                const std::vector<gravis::d2Vector>& offsets,
                int threads);

        double f(const std::vector<double>& x, void* tempStorage) const;
        void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;

    private:

        int pc, fc, threads;
        const std::vector<std::vector<Image<RFLOAT>>>& correlation;
        const std::vector<double>& velWgh;
        const std::vector<double>& accWgh;
        const std::vector<std::vector<std::vector<double>>>& divWgh;
        const std::vector<gravis::d2Vector>& offsets;
};

#endif

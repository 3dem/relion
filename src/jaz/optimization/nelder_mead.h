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

#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

#include <vector>
#include "optimization.h"

class NelderMead
{
    public:

        static std::vector<double> optimize(
                const std::vector<double>& initial,
                const Optimization& opt,
                double initialStep, double tolerance, long maxIters,
                double alpha = 1.0, double gamma = 2.0,
                double rho = 0.5, double sigma = 0.5,
                bool verbose = false, double* minCost = 0);

        static void test();
};

#endif

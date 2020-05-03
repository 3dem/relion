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

#include "gradient_descent.h"
#include <iostream>

std::vector<double> GradientDescent::optimize(
    const std::vector<double> &initial,
    const DifferentiableOptimization &opt,
    double step, double minStep, double minDiff, long maxIters,
    double inertia, bool verbose)
{
    std::vector<double> x = initial;
    std::vector<double> last_x = x;

    const long n = initial.size();

    std::vector<double> g(n, 0.0), v(n, 0.0);

    void* tempStorage = opt.allocateTempStorage();

    double last_f = opt.f(initial, tempStorage);

    double act_step = step;
    int goodSince = 0, accAfter = 5;

    if (verbose)
    {
        std::cout << "initial: " << last_f << "\n";
    }

    for (int i = 0; i < maxIters; i++)
    {
        opt.grad(x, g, tempStorage);

        for (int j = 0; j < n; j++)
        {
            v[j] = inertia * v[j] - (1.0 - inertia) * act_step * g[j];
            x[j] += v[j];
        }

        double f = opt.f(x, tempStorage);

        if (verbose)
        {
			opt.report(i, f, x);
            std::cout << i << ": " << f << " (" << act_step << ") [" << (f - last_f) << "]";
        }

        if (f > last_f)
        {
            if (verbose)
            {
                std::cout << " *\n";
            }

            for (int j = 0; j < n; j++)
            {
                x[j] = last_x[j];
                v[j] = 0.0;
            }

            act_step /= 2.0;

            if (act_step < minStep) break;
        }
        else
        {
            if (verbose)
            {
                std::cout << "\n";
            }

            if (last_f - f < minDiff) break;

            for (int j = 0; j < n; j++)
            {
                last_x[j] = x[j];
            }

            if (act_step < step/2.0)
            {
                goodSince++;

                if (goodSince >= accAfter)
                {
                    goodSince = 0;
                    act_step *= 2.0;
                }
            }
            else
            {
                act_step = step;
            }

            last_f = f;
        }

    }

    opt.deallocateTempStorage(tempStorage);

    return x;
}

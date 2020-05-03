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

#include <src/jaz/single_particle/local_motion_fit.h>
#include <src/jaz/single_particle/interpolation.h>
#include <omp.h>

using namespace gravis;

LocalMotionFit::LocalMotionFit(const std::vector<std::vector<Image<RFLOAT>>>& correlation,
        const std::vector<double>& velWgh,
        const std::vector<double>& accWgh,
        const std::vector<std::vector<std::vector<double>>> &divWgh,
        const std::vector<d2Vector>& offsets,
        int threads)
:   pc(correlation.size()),
    fc(correlation[0].size()),
    threads(threads),
    correlation(correlation),
    velWgh(velWgh),
    accWgh(accWgh),
    divWgh(divWgh),
    offsets(offsets)
{
}

double LocalMotionFit::f(const std::vector<double> &x, void* tempStorage) const
{
    double e_tot = 0.0;

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    {
        double e = 0.0;

        for (int f = 0; f < fc; f++)
        {
            const double xpf = x[2*(p*fc + f) + 0];
            const double ypf = x[2*(p*fc + f) + 1];

            e -= Interpolation::cubicXY(correlation[p][f],
                    xpf+offsets[f].x, ypf+offsets[f].y, 0, 0, true);

            if (f > 0 && f < fc-1)
            {
                const double xpfn = x[2*(p*fc + f - 1) + 0];
                const double ypfn = x[2*(p*fc + f - 1) + 1];

                const double xpfp = x[2*(p*fc + f + 1) + 0];
                const double ypfp = x[2*(p*fc + f + 1) + 1];

                const double ax = xpfn + xpfp - 2.0 * xpf;
                const double ay = ypfn + ypfp - 2.0 * ypf;

                e += accWgh[f] * (ax * ax + ay * ay);
            }

            if (f > 0)
            {
                const double xpfn = x[2*(p*fc + f - 1) + 0];
                const double ypfn = x[2*(p*fc + f - 1) + 1];

                const double vx = xpf - xpfn;
                const double vy = ypf - ypfn;

                e += velWgh[f-1] * (vx * vx + vy * vy);

                for (int q = p+1; q < pc; q++)
                {
                    if (divWgh[f-1][p][q] <= 0.0) continue;

                    const double xqf = x[2*(q*fc + f) + 0];
                    const double yqf = x[2*(q*fc + f) + 1];

                    const double xqfn = x[2*(q*fc + f - 1) + 0];
                    const double yqfn = x[2*(q*fc + f - 1) + 1];

                    const double cx = (xpf - xpfn) - (xqf - xqfn);
                    const double cy = (ypf - ypfn) - (yqf - yqfn);

                    e += divWgh[f-1][p][q] * (cx * cx + cy * cy);
                }
            }
        }

        #pragma omp atomic
        e_tot += e;
    }

    return e_tot;
}


void LocalMotionFit::grad(const std::vector<double> &x, std::vector<double> &gradDest, void* tempStorage) const
{
    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        gradDest[2*(p*fc + f) + 0] = 0.0;
        gradDest[2*(p*fc + f) + 1] = 0.0;
    }

    std::vector<std::vector<double>> tempGrad(threads);

    for (int i = 0; i < threads; i++)
    {
        tempGrad[i] = std::vector<double>(2*pc*fc, 0.0);
    }

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++)
    {
        int th = omp_get_thread_num();

        for (int f = 0; f < fc; f++)
        {
            const double xpf = x[2*(p*fc + f) + 0];
            const double ypf = x[2*(p*fc + f) + 1];

            gravis::t2Vector<RFLOAT> g = Interpolation::cubicXYgrad(
                    correlation[p][f], xpf+offsets[f].x, ypf+offsets[f].y, 0, 0, true);

            tempGrad[th][2*(p*fc + f) + 0] -= g.x;
            tempGrad[th][2*(p*fc + f) + 1] -= g.y;

            if (f > 0 && f < fc-1)
            {
                const double xpfn = x[2*(p*fc + f - 1) + 0];
                const double ypfn = x[2*(p*fc + f - 1) + 1];

                const double xpfp = x[2*(p*fc + f + 1) + 0];
                const double ypfp = x[2*(p*fc + f + 1) + 1];

                const double ax = xpfn + xpfp - 2.0 * xpf;
                const double ay = ypfn + ypfp - 2.0 * ypf;

                tempGrad[th][2*(p*fc + f - 1) + 0] += 2.0 * accWgh[f] * ax;
                tempGrad[th][2*(p*fc + f - 1) + 1] += 2.0 * accWgh[f] * ay;

                tempGrad[th][2*(p*fc + f) + 0] -= 4.0 * accWgh[f] * ax;
                tempGrad[th][2*(p*fc + f) + 1] -= 4.0 * accWgh[f] * ay;

                tempGrad[th][2*(p*fc + f + 1) + 0] += 2.0 * accWgh[f] * ax;
                tempGrad[th][2*(p*fc + f + 1) + 1] += 2.0 * accWgh[f] * ay;
            }

            if (f > 0)
            {
                const double xpfn = x[2*(p*fc + f - 1) + 0];
                const double ypfn = x[2*(p*fc + f - 1) + 1];

                const double vx = xpf - xpfn;
                const double vy = ypf - ypfn;

                tempGrad[th][2*(p*fc + f) + 0] += 2.0 * velWgh[f-1] * vx;
                tempGrad[th][2*(p*fc + f) + 1] += 2.0 * velWgh[f-1] * vy;
                tempGrad[th][2*(p*fc + f - 1) + 0] -= 2.0 * velWgh[f-1] * vx;
                tempGrad[th][2*(p*fc + f - 1) + 1] -= 2.0 * velWgh[f-1] * vy;

                for (int q = p+1; q < pc; q++)
                {
                    if (divWgh[f-1][p][q] <= 0.0) continue;

                    const double xqf = x[2*(q*fc + f) + 0];
                    const double yqf = x[2*(q*fc + f) + 1];

                    const double xqfn = x[2*(q*fc + f - 1) + 0];
                    const double yqfn = x[2*(q*fc + f - 1) + 1];

                    const double cx = (xpf - xpfn) - (xqf - xqfn);
                    const double cy = (ypf - ypfn) - (yqf - yqfn);

                    const double wgh = divWgh[f-1][p][q];

                    tempGrad[th][2*(p*fc + f - 1) + 0] -= 2.0 * wgh * cx;
                    tempGrad[th][2*(p*fc + f - 1) + 1] -= 2.0 * wgh * cy;

                    tempGrad[th][2*(p*fc + f) + 0] += 2.0 * wgh * cx;
                    tempGrad[th][2*(p*fc + f) + 1] += 2.0 * wgh * cy;

                    tempGrad[th][2*(q*fc + f - 1) + 0] += 2.0 * wgh * cx;
                    tempGrad[th][2*(q*fc + f - 1) + 1] += 2.0 * wgh * cy;

                    tempGrad[th][2*(q*fc + f) + 0] -= 2.0 * wgh * cx;
                    tempGrad[th][2*(q*fc + f) + 1] -= 2.0 * wgh * cy;

                }
            }
        }
    }

    for (int i = 0; i < threads; i++)
    {
        for (int j = 0; j < 2*pc*fc; j++)
        {
            gradDest[j] += tempGrad[i][j];
        }
    }
    /*std::cout << "x:\n";

    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        std::cout << x[2*(p*fc + f) + 0] << ", " << x[2*(p*fc + f) + 1] << "\n";
    }

    std::cout << "\n";

    std::cout << "grad:\n";

    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        std::cout << gradDest[2*(p*fc + f) + 0] << ", " << gradDest[2*(p*fc + f) + 1] << "\n";
    }

    std::cout << "\n";*/
}

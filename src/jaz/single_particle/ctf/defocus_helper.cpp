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

#include "defocus_helper.h"

#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/refinement_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>

#include <src/projector.h>

using namespace gravis;

RFLOAT DefocusHelper::findDefocus1D(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        const Image<RFLOAT>& weight,
        const CTF &ctf0, RFLOAT angpix,
        double *destU, double *destV,
        RFLOAT range, int steps,
        int recDepth, RFLOAT recScale)
{
    const RFLOAT delta = ctf0.DeltafV - ctf0.DeltafU;
    CTF ctf(ctf0);

    double minErr = std::numeric_limits<double>::max();
    double bestU = ctf0.DeltafU;
    double bestV = ctf0.DeltafV;

    for (int s = -steps/2; s <= steps/2; s++)
    {
        const RFLOAT u = ctf0.DeltafU + s * range / (steps/2);
        const RFLOAT v = u + delta;

        ctf.DeltafU = u;
        ctf.DeltafV = v;
        ctf.initialise();

        double err = RefinementHelper::squaredDiff(prediction, observation, ctf, angpix, weight);

        if (err < minErr)
        {
            minErr = err;
            bestU = u;
            bestV = v;
        }
    }

    if (recDepth > 0)
    {
        ctf.DeltafU = bestU;
        ctf.DeltafV = bestV;
        ctf.initialise();

        findDefocus1D(prediction, observation, weight,
                ctf, angpix, destU, destV,
                range/recScale, steps,
                recDepth-1, recScale);
    }
    else
    {
        *destU = bestU;
        *destV = bestV;
    }

    return minErr;
}

void DefocusHelper::findAstigmatismNM(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, RFLOAT angpix,
        double *destU, double *destV, double *destPhi)
{
    AstigmatismOptimizationAcc opt(prediction, observation, weight, ctf0, false, false, angpix);

    std::vector<double> initial = opt.getInitialParams();

    std::vector<double> params = NelderMead::optimize(initial, opt, 50.0, 1.0, 1000);

    *destU = opt.getU(params);
    *destV = opt.getV(params);
    *destPhi = opt.getPhi(params);
}

void DefocusHelper::findAstigmatismAndPhaseNM(
        const std::vector<Image<Complex>>& prediction,
        const std::vector<Image<Complex>>& observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, RFLOAT angpix,
        double *destU, double *destV, double *destPhi, double *destPhase)
{
    AstigmatismOptimizationAcc opt(prediction, observation, weight, ctf0, true, false, angpix);

    std::vector<double> initial = opt.getInitialParams();

    std::vector<double> params = NelderMead::optimize(initial, opt, 5.0, 0.01, 100000);

    *destU = opt.getU(params);
    *destV = opt.getV(params);
    *destPhi = opt.getPhi(params);
    *destPhase = opt.getPhase(params);
}

void DefocusHelper::findAstigmatismPhaseAndCsNM(
        const std::vector<Image<Complex>>& prediction,
        const std::vector<Image<Complex>>& observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, RFLOAT angpix,
        double *destU, double *destV, double *destPhi, double *destPhase, double* destCs)
{
    AstigmatismOptimizationAcc opt(prediction, observation, weight, ctf0, true, true, angpix);

    std::vector<double> initial = opt.getInitialParams();

    std::vector<double> params = NelderMead::optimize(initial, opt, 5.0, 0.01, 100000);

    *destU = opt.getU(params);
    *destV = opt.getV(params);
    *destPhi = opt.getPhi(params);
    *destPhase = opt.getPhase(params);
    *destCs = opt.getCs(params);
}

void DefocusHelper::findAstigmatismNM(
        const std::vector<Image<Complex>>& prediction,
        const std::vector<Image<Complex>>& observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, RFLOAT angpix,
        double *destU, double *destV, double *destPhi)
{
    AstigmatismOptimizationAcc opt(prediction, observation, weight, ctf0, false, false, angpix);

    std::vector<double> initial = opt.getInitialParams();

    std::vector<double> params = NelderMead::optimize(initial, opt, 50.0, 1.0, 1000);

    *destU = opt.getU(params);
    *destV = opt.getV(params);
    *destPhi = opt.getPhi(params);
}

std::vector<d2Vector> DefocusHelper::diagnoseDefocus(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, RFLOAT angpix,
        double range, int steps, int threads)
{
    const RFLOAT delta = ctf0.DeltafV - ctf0.DeltafU;

    std::vector<d2Vector> out(steps);

    #pragma omp parallel for num_threads(threads)
    for (int s = 0; s < steps; s++)
    {
        CTF ctf(ctf0);
        const RFLOAT u = ctf0.DeltafU + (s - steps/2) * range / (double)steps;
        const RFLOAT v = u + delta;

        ctf.DeltafU = u;
        ctf.DeltafV = v;
        ctf.initialise();

        out[s][0] = u;
        out[s][1] = RefinementHelper::squaredDiff(prediction, observation, ctf, angpix, weight);
    }

    return out;
}

AstigmatismOptimization::AstigmatismOptimization(
        const Image<Complex>& prediction, const Image<Complex>& observation,
        const Image<RFLOAT>& weight, const CTF& ctf0, RFLOAT angpix)
:   prediction(prediction),
    observation(observation),
    weight(weight),
    ctf0(ctf0),
    angpix(angpix)
{
}

double AstigmatismOptimization::f(const std::vector<double>& x) const
{
    CTF ctf(ctf0);

    ctf.DeltafU = x[0];
    ctf.DeltafV = x[1];
    ctf.azimuthal_angle = x[2];
    ctf.initialise();

    return RefinementHelper::squaredDiff(prediction, observation, ctf, angpix, weight);
}

AstigmatismOptimizationAcc::AstigmatismOptimizationAcc(const Image<Complex>& prediction,
        const Image<Complex>& observation,
        const Image<RFLOAT>& weight,
        const CTF& ctf0,
        bool phaseShift,
        bool spherAberr,
        RFLOAT angpix,
        RFLOAT phiScale,
        RFLOAT csScale)
:   ctf0(ctf0),
    phaseShift(phaseShift),
    spherAberr(spherAberr),
    angpix(angpix),
    phiScale(phiScale),
    csScale(csScale)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    data = Image<Complex>(w,h);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        const Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);
        const RFLOAT vw = DIRECT_A2D_ELEM(weight.data, y, x);

        const RFLOAT x2 = vx.real*vx.real + vx.imag*vx.imag;

        const RFLOAT yxb = x2 > 0.0? (vy.real*vx.real + vy.imag*vx.imag)/x2 : 0.0;
        const RFLOAT wp = vw * x2;

        DIRECT_A2D_ELEM(data.data, y, x) = Complex(yxb, wp);
    }
}

AstigmatismOptimizationAcc::AstigmatismOptimizationAcc(
        const std::vector<Image<Complex>>& prediction,
        const std::vector<Image<Complex>>& observation,
        const Image<RFLOAT>& weight,
        const CTF& ctf0,
        bool phaseShift,
        bool spherAberr,
        RFLOAT angpix,
        RFLOAT phiScale,
        RFLOAT csScale)
:   ctf0(ctf0),
    phaseShift(phaseShift),
    spherAberr(spherAberr),
    angpix(angpix),
    phiScale(phiScale),
    csScale(csScale)
{
    const long w = prediction[0].data.xdim;
    const long h = prediction[0].data.ydim;

    data = Image<Complex>(w,h);
    data.data.initZeros();

    for (int i = 0; i < prediction.size(); i++)
    {
        for (long y = 0; y < h; y++)
        for (long x = 0; x < w; x++)
        {
            Complex vx = DIRECT_A2D_ELEM(prediction[i].data, y, x);
            const Complex vy = DIRECT_A2D_ELEM(observation[i].data, y, x);
            const RFLOAT vw = DIRECT_A2D_ELEM(weight.data, y, x);

            const RFLOAT x2 = vx.real*vx.real + vx.imag*vx.imag;

            const RFLOAT yxb = x2 > 0.0? (vy.real*vx.real + vy.imag*vx.imag)/x2 : 0.0;
            const RFLOAT wp = vw * x2;

            DIRECT_A2D_ELEM(data.data, y, x) += Complex(yxb, wp);
        }
    }
}

double AstigmatismOptimizationAcc::f(const std::vector<double> &x, void* tempStorage) const
{
    CTF ctf(ctf0);

    ctf.DeltafU = x[0];
    ctf.DeltafV = x[1];
    ctf.azimuthal_angle = x[2]/phiScale;

    if (phaseShift) ctf.phase_shift = x[3]/phiScale;
    if (spherAberr) ctf.Cs = x[phaseShift?4:3]/csScale;

    ctf.initialise();

    const long w = data.data.xdim;
    const long h = data.data.ydim;

    double out = 0.0;

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        Complex vd = DIRECT_A2D_ELEM(data.data, y, x);

        RFLOAT vm = ctfImg(y,x);
        RFLOAT dx = vd.real - vm;
        out += vd.imag * dx * dx;
    }

    return out;
}

double AstigmatismOptimizationAcc::getU(const std::vector<double> &x)
{
    return x[0];
}

double AstigmatismOptimizationAcc::getV(const std::vector<double> &x)
{
    return x[1];
}

double AstigmatismOptimizationAcc::getPhi(const std::vector<double> &x)
{
    return x[2] / phiScale;
}

double AstigmatismOptimizationAcc::getPhase(const std::vector<double> &x)
{
    return x[3] / phiScale;
}

double AstigmatismOptimizationAcc::getCs(const std::vector<double> &x)
{
    return x[phaseShift?4:3] / csScale;
}

std::vector<double> AstigmatismOptimizationAcc::getInitialParams()
{
    int num = 3;
    if (phaseShift) num++;
    if (spherAberr) num++;

    std::vector<double> initial(num);

    initial[0] = ctf0.DeltafU;
    initial[1] = ctf0.DeltafU;
    initial[2] = phiScale * ctf0.azimuthal_angle;

    if (phaseShift) initial[3] = phiScale * ctf0.phase_shift;
    if (spherAberr) initial[phaseShift?4:3] = csScale * ctf0.Cs;

    return initial;
}

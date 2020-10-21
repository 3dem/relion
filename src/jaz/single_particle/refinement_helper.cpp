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

#include <src/jaz/single_particle/refinement_helper.h>
#include <src/jaz/single_particle/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/single_particle/vtk_helper.h>


using namespace gravis;

void RefinementHelper::drawFSC(const MetaDataTable *mdt, std::vector<double>& dest1D,
                               Image<RFLOAT> &dest, double thresh)
{
    const int n = mdt->numberOfObjects();
    const int w = 2*(n-1);
    const int h = 2*(n-1);

    dest1D = std::vector<double>(n);

    for (int i = 0; i < n; i++)
    {
        int idx;
        mdt->getValue(EMDL_SPECTRAL_IDX, idx, i);
        mdt->getValue(EMDL_POSTPROCESS_FSC_TRUE, dest1D[i], i);

        if (dest1D[i] < thresh) dest1D[i] = 0.0;
    }

    dest = Image<RFLOAT>(n,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < n; x++)
    {
        double xx = x;
        double yy = y > h/2? y - h : y;

        double r = sqrt(xx*xx + yy*yy);
        int ri = (int)(r+0.5);
        if (ri > w/2) ri = w/2;

        DIRECT_A2D_ELEM(dest.data, y, x) = dest1D[ri];
    }
}

void RefinementHelper::computeSNR(const MetaDataTable *mdt, Image<RFLOAT> &dest, double eps)
{
    const int n = mdt->numberOfObjects();
    const int w = 2*(n-1);
    const int h = 2*(n-1);

    std::vector<double> snr(n);

    for (int i = 0; i < n; i++)
    {
        int idx;
        double fsc;
        mdt->getValue(EMDL_SPECTRAL_IDX, idx, i);
        mdt->getValue(EMDL_POSTPROCESS_FSC_TRUE, fsc, i);

        if (fsc > 1.0 - eps) fsc = 1.0 - eps;
        //else if (fsc < eps) fsc = 0.0;

        snr[i] = fsc / (1.0 - fsc);
    }

    dest = Image<RFLOAT>(n,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < n; x++)
    {
        double xx = x;
        double yy = y > h/2? y - h : y;

        double r = sqrt(xx*xx + yy*yy);
        int ri = (int)(r+0.5);
        if (ri > w/2) ri = w/2;

        DIRECT_A2D_ELEM(dest.data, y, x) = snr[ri];
    }
}

void RefinementHelper::computeSigInvSq(const MetaDataTable *mdt, const std::vector<double>& signalPow,
                                       Image<RFLOAT> &dest, double eps)
{
    const int n = mdt->numberOfObjects();
    const int w = 2*(n-1);
    const int h = 2*(n-1);

    std::vector<double> sigInvSq(n);

    for (int i = 0; i < n; i++)
    {
        int idx;
        double fsc;
        mdt->getValue(EMDL_SPECTRAL_IDX, idx, i);
        mdt->getValue(EMDL_POSTPROCESS_FSC_TRUE, fsc, i);

        if (fsc > 1.0 - eps) fsc = 1.0 - eps;
        //else if (fsc < eps) fsc = 0.0;

        double snr = fsc / (1.0 - fsc);
        double sigPow = signalPow[i];

        if (sigPow < eps) sigPow = eps;

        sigInvSq[i] = snr / sigPow;
    }

    dest = Image<RFLOAT>(n,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < n; x++)
    {
        double xx = x;
        double yy = y > h/2? y - h : y;

        double r = sqrt(xx*xx + yy*yy);
        int ri = (int)(r+0.5);
        if (ri > w/2) ri = w/2;

        DIRECT_A2D_ELEM(dest.data, y, x) = sigInvSq[ri];
    }
}

Image<RFLOAT> RefinementHelper::correlation(const Image<Complex> &prediction, const Image<Complex> &observation)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    Image<RFLOAT> out(w,h);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);

        DIRECT_A2D_ELEM(out.data, y, x) = (vy.real * vx.real + vy.imag * vx.imag);
    }

    return out;
}

Image<RFLOAT> RefinementHelper::correlation(
        const std::vector<Image<Complex> >& predictions,
        const std::vector<Image<Complex> >& observations)
{
    const long w = predictions[0].data.xdim;
    const long h = predictions[0].data.ydim;
    const long c = predictions.size();

    Image<RFLOAT> out(w,h);
    out.data.initZeros();

    for (long i = 0; i < c; i++)
    {
        for (long y = 0; y < h; y++)
        for (long x = 0; x < w; x++)
        {
            Complex vx = DIRECT_A2D_ELEM(predictions[i].data, y, x);
            Complex vy = DIRECT_A2D_ELEM(observations[i].data, y, x);

            DIRECT_A2D_ELEM(out.data, y, x) += (vy.real * vx.real + vy.imag * vx.imag);
        }
    }

    return out;
}

void RefinementHelper::addToQR(const Image<Complex>& prediction,
                        const Image<Complex>& observation,
                        Image<Complex>& q, Image<RFLOAT>& r)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);

        DIRECT_A2D_ELEM(q.data, y, x) += vy.conj() * vx;
        DIRECT_A2D_ELEM(r.data, y, x) += vx.norm();
    }
}

void RefinementHelper::addToPQR(const Image<Complex>& prediction,
                        const Image<Complex>& observation,
                        Image<RFLOAT>& p, Image<Complex>& q, Image<RFLOAT>& r)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);

        DIRECT_A2D_ELEM(p.data, y, x) += vx.norm();
        DIRECT_A2D_ELEM(q.data, y, x) += vy.conj() * vx;
        DIRECT_A2D_ELEM(r.data, y, x) += vx.norm();
    }
}

double RefinementHelper::squaredDiff(
        const Image<Complex>& prediction,
        const Image<Complex>& observation,
        CTF& ctf, RFLOAT angpix, const Image<RFLOAT>& weight)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    double out = 0.0;
	
	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        const Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);
        const RFLOAT vw = DIRECT_A2D_ELEM(weight.data, y, x);

        RFLOAT vm = ctfImg(y,x);
		
        out += vw * (vy - vm * vx).norm();
    }

    return out;
}

double RefinementHelper::squaredDiff(
        const std::vector<Image<Complex> > &predictions,
        const std::vector<Image<Complex> > &observations,
        CTF &ctf, RFLOAT angpix, const Image<RFLOAT> &weight)
{
    double out = 0.0;

    for (long i = 0; i < predictions.size(); i++)
    {
        out += squaredDiff(predictions[i], observations[i], ctf, angpix, weight);
    }

    return out;
}

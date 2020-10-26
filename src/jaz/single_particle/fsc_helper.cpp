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

#include <src/jaz/single_particle/fsc_helper.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/optimization/nelder_mead.h>

using namespace gravis;

void FscHelper::computeFscTable(
        const std::vector<std::vector<Image<Complex> > >& frames,
        const std::vector<Image<Complex> >& predictions,
        Image<RFLOAT>& table, Image<RFLOAT>& weight)
{
    const int w = predictions[0].data.xdim;
    const int pc = frames.size();
    const int fc = frames[0].size();

    table = Image<RFLOAT>(w,fc);
    weight = Image<RFLOAT>(w,fc);

    table.data.initZeros();
    weight.data.initZeros();

    Image<RFLOAT> weight1 = Image<RFLOAT>(w,fc);
    Image<RFLOAT> weight2 = Image<RFLOAT>(w,fc);

    weight1.data.initZeros();
    weight2.data.initZeros();

    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frames[p][f]())
        {
            int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));

            if (idx >= w)
            {
                continue;
            }

            Complex z1 = DIRECT_A3D_ELEM(frames[p][f](), k, i, j);
            Complex z2 = DIRECT_A3D_ELEM(predictions[p](), k, i, j);

            DIRECT_A2D_ELEM(table.data, f, idx) += z1.real * z2.real + z1.imag * z2.imag;
            DIRECT_A2D_ELEM(weight1.data, f, idx) += z1.norm();
            DIRECT_A2D_ELEM(weight2.data, f, idx) += z2.norm();
        }
    }

    for (int f = 0; f < fc; f++)
    for (int x = 0; x < w; x++)
    {
        RFLOAT w1 = DIRECT_A2D_ELEM(weight1.data, f, x);
        RFLOAT w2 = DIRECT_A2D_ELEM(weight2.data, f, x);
        RFLOAT ww = sqrt(w1 * w2);

        DIRECT_A2D_ELEM(weight.data, f, x) = ww;
        DIRECT_A2D_ELEM(table.data, f, x) /= ww;
    }
}

void FscHelper::computeFscRow(
        const MultidimArray<Complex> &data0, const MultidimArray<Complex> &data1,
        int row, Image<RFLOAT> &table, Image<RFLOAT> &weight)
{
    const int w = data0.xdim;

    std::vector<double> weight1(w, 0.0), weight2(w, 0.0), data(w, 0.0);

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(data0)
    {
        int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));

        if (idx >= w)
        {
            continue;
        }

        Complex z1 = DIRECT_A3D_ELEM(data0, k, i, j);
        Complex z2 = DIRECT_A3D_ELEM(data1, k, i, j);

        data[idx] += z1.real * z2.real + z1.imag * z2.imag;
        weight1[idx] += z1.norm();
        weight2[idx] += z2.norm();
    }

    for (int x = 0; x < w; x++)
    {
        if (x >= table.data.xdim) continue;

        RFLOAT ww = sqrt(weight1[x] * weight2[x]);

        DIRECT_A2D_ELEM(table.data, row, x) = ww > 0.0? data[x] / ww : 0.0;
        DIRECT_A2D_ELEM(weight.data, row, x) = ww;
    }
}

void FscHelper::initFscTable(int kc, int tc, Image<RFLOAT> &table, Image<RFLOAT> &weight0, Image<RFLOAT> &weight1)
{
    table = Image<RFLOAT>(kc,tc);
    weight0 = Image<RFLOAT>(kc,tc);
    weight1 = Image<RFLOAT>(kc,tc);

    table.data.initZeros();
    weight0.data.initZeros();
    weight1.data.initZeros();
}

void FscHelper::updateFscTable(const std::vector<Image<Complex> > &frames,
                               const Image<Complex> &prediction, double scale,
                               Image<RFLOAT> &table,
                               Image<RFLOAT> &weight0,
                               Image<RFLOAT> &weight1)
{
    const int fc = frames.size();

    for (int f = 0; f < fc; f++)
    {
        updateFscTable(frames[f], f, prediction, scale, table, weight0, weight1);
    }
}

void FscHelper::updateFscTable(const Image<Complex> &frame,
                               int f,
                               const Image<Complex> &prediction, double scale,
                               Image<RFLOAT> &table,
                               Image<RFLOAT> &weight0,
                               Image<RFLOAT> &weight1)
{
    const int w = prediction.data.xdim;

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frame())
    {
        int idx = ROUND(scale * sqrt(kp*kp + ip*ip + jp*jp));

        if (idx >= w)
        {
            continue;
        }

        Complex z1 = DIRECT_A3D_ELEM(frame(), k, i, j);
        Complex z2 = DIRECT_A3D_ELEM(prediction(), k, i, j);

        DIRECT_A2D_ELEM(table.data, f, idx) += z1.real * z2.real + z1.imag * z2.imag;
        DIRECT_A2D_ELEM(weight0.data, f, idx) += z1.norm();
        DIRECT_A2D_ELEM(weight1.data, f, idx) += z2.norm();
    }
}

void FscHelper::updateFscTableVelWgh(
        const std::vector<Image<Complex> > &frames,
        const std::vector<d2Vector> &velocities,
        const Image<Complex> &prediction,
        Image<RFLOAT> &table,
        Image<RFLOAT> &weight0,
        Image<RFLOAT> &weight1)
{
    const int w = prediction.data.xdim;
    const int fc = frames.size();

    for (int f = 0; f < fc; f++)
    {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frames[f]())
        {
            int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));

            if (idx >= w)
            {
                continue;
            }

            double kv = (ip * velocities[f].y + jp * velocities[f].x)/(double)w;
            double wgh = kv < 1e-20? 1.0 : sin(PI*kv) / (PI*kv);
            //double wgh = exp(-0.5*kv*kv/0.5);

            Complex z1 = DIRECT_A3D_ELEM(frames[f](), k, i, j);
            Complex z2 = wgh * DIRECT_A3D_ELEM(prediction(), k, i, j);

            DIRECT_A2D_ELEM(table.data, f, idx) += z1.real * z2.real + z1.imag * z2.imag;
            DIRECT_A2D_ELEM(weight0.data, f, idx) += z1.norm();
            DIRECT_A2D_ELEM(weight1.data, f, idx) += z2.norm();
        }
    }
}

void FscHelper::updateVelFscTable(
        const std::vector<Image<Complex> > &frames,
        const std::vector<d2Vector> &velocities,
        const Image<Complex> &prediction,
        Image<RFLOAT> &table,
        Image<RFLOAT> &weight0,
        Image<RFLOAT> &weight1,
        int kmin, int kmax)
{
    const int w = prediction.data.xdim;
    const int fc = frames.size();

    if (kmax < 0) kmax = w;

    for (int f = 0; f < fc; f++)
    {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frames[f]())
        {
            int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));

            if (idx >= w || idx < kmin || idx > kmax)
            {
                continue;
            }

            double kv = ip * velocities[f].y + jp * velocities[f].x;
            int kvidx = ROUND(kv);

            if (kvidx < 0) kvidx = -kvidx;
            if (kvidx >= table().xdim) continue;

            Complex z1 = DIRECT_A3D_ELEM(frames[f](), k, i, j);
            Complex z2 = DIRECT_A3D_ELEM(prediction(), k, i, j);

            DIRECT_A2D_ELEM(table.data, f, kvidx) += z1.real * z2.real + z1.imag * z2.imag;
            DIRECT_A2D_ELEM(weight0.data, f, kvidx) += z1.norm();
            DIRECT_A2D_ELEM(weight1.data, f, kvidx) += z2.norm();
        }
    }
}

void FscHelper::mergeFscTables(const std::vector<Image<RFLOAT> > &tables,
                               const std::vector<Image<RFLOAT> > &weights0,
                               const std::vector<Image<RFLOAT> > &weights1,
                               Image<RFLOAT> &table, Image<RFLOAT> &weight)
{
    const int w = tables[0].data.xdim;
    const int fc = tables[0].data.ydim;
    const int mgc = tables.size();

    Image<RFLOAT> tableSum = Image<RFLOAT>(w,fc);
    Image<RFLOAT> weightSum0 = Image<RFLOAT>(w,fc);
    Image<RFLOAT> weightSum1 = Image<RFLOAT>(w,fc);

    tableSum.data.initZeros();
    weightSum0.data.initZeros();
    weightSum1.data.initZeros();

    table = Image<RFLOAT>(w,fc);
    weight = Image<RFLOAT>(w,fc);

    for (int m = 0; m < mgc; m++)
    {
        for (int f = 0; f < fc; f++)
        for (int x = 0; x < w; x++)
        {
            DIRECT_A2D_ELEM(tableSum.data, f, x) += DIRECT_A2D_ELEM(tables[m].data, f, x);
            DIRECT_A2D_ELEM(weightSum0.data, f, x) += DIRECT_A2D_ELEM(weights0[m].data, f, x);
            DIRECT_A2D_ELEM(weightSum1.data, f, x) += DIRECT_A2D_ELEM(weights1[m].data, f, x);
        }
    }

    for (int f = 0; f < fc; f++)
    for (int x = 0; x < w; x++)
    {
        RFLOAT w1 = DIRECT_A2D_ELEM(weightSum0.data, f, x);
        RFLOAT w2 = DIRECT_A2D_ELEM(weightSum1.data, f, x);
        RFLOAT ww = sqrt(w1 * w2);

        DIRECT_A2D_ELEM(weight.data, f, x) = ww;

        if (ww > 0.0)
        {
            DIRECT_A2D_ELEM(table.data, f, x) = DIRECT_A2D_ELEM(tableSum.data, f, x) / ww;
        }
        else
        {
            DIRECT_A2D_ELEM(table.data, f, x) = 0.0;
        }
    }
}

double FscHelper::computeTsc(
        const std::vector<Image<RFLOAT>> &tables,
        const std::vector<Image<RFLOAT>> &weights0,
        const std::vector<Image<RFLOAT>> &weights1,
        int k0, int k1)
{
    //const int kc = tables[0].data.xdim;
    const int fc = tables[0].data.ydim;
    const int tc = tables.size();

    double tSum = 0.0, w0Sum = 0.0, w1Sum = 0.0;

    for (int t = 0; t < tc; t++)
    for (int f = 0; f < fc; f++)
    for (int k = k0; k < k1; k++)
    {
        tSum += tables[t](f,k);
        w0Sum += weights0[t](f,k);
        w1Sum += weights1[t](f,k);
    }

    double ww = sqrt(w0Sum*w1Sum);

    if (ww > 0.0)
    {
        return tSum / ww;
    }
    else
    {
        return 0.0;
    }
}

void FscHelper::computeNoiseSq(
        std::vector<std::vector<Image<Complex> > > frames,
        std::vector<Image<Complex> > predictions, Image<RFLOAT> &sigma2)
{
    const int w = predictions[0].data.xdim;
    const int pc = frames.size();
    const int fc = frames[0].size();

    sigma2 = Image<RFLOAT>(w,fc);
    Image<RFLOAT> count(w,fc);
    sigma2.data.initZeros();
    count.data.initZeros();

    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frames[p][f]())
        {
            int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));

            if (idx >= w)
            {
                continue;
            }

            Complex z1 = DIRECT_A3D_ELEM(frames[p][f](), k, i, j);
            Complex z2 = DIRECT_A3D_ELEM(predictions[p](), k, i, j);

            DIRECT_A2D_ELEM(sigma2.data, f, idx) += (z2 - z1).norm();
            DIRECT_A2D_ELEM(count.data, f, idx) += 1.0;
        }
    }

    for (int f = 0; f < fc; f++)
    for (int x = 0; x < w; x++)
    {
        if (DIRECT_A2D_ELEM(count.data, f, x) > 0.0)
        {
            DIRECT_A2D_ELEM(sigma2.data, f, x) /= DIRECT_A2D_ELEM(count.data, f, x);
        }
    }
}

Image<RFLOAT> FscHelper::computeSignalSq(const Image<RFLOAT> &sigma2, const Image<RFLOAT> &frc)
{
    const int kc = sigma2.data.xdim;
    const int fc = sigma2.data.ydim;

    Image<RFLOAT> out(kc,fc);

    const RFLOAT eps = 1e-3;

    for (int f = 0; f < fc; f++)
    for (int k = 0; k < kc; k++)
    {
        RFLOAT c = DIRECT_A2D_ELEM(frc.data, f, k);
        RFLOAT s2 = DIRECT_A2D_ELEM(sigma2.data, f, k);

        if (c < eps)
        {
            DIRECT_A2D_ELEM(out.data, f, k) = 0.0;
        }
        else
        {
            if (c > 1.0 - eps) c = 1.0 - eps;
            RFLOAT snr2 = c / (1.0 - c);
            DIRECT_A2D_ELEM(out.data, f, k) = snr2*s2;
        }
    }

    return out;
}

std::vector<d2Vector> FscHelper::fitBfactorsNM(const Image<RFLOAT> &tau2, const Image<RFLOAT> &weight, int cutoff)
{
    const int kc = tau2.data.xdim;
    const int ic = tau2.data.ydim;
    std::vector<d2Vector> out(ic);

    const double Bw = 0.001, Cw = 10.0;

    Image<RFLOAT> tauRel(kc,ic);
    std::vector<RFLOAT> tauAvg(kc,0.0);

    for (int k = 0; k < kc; k++)
    {
        double ta = 0.0;
        double ts = 0.0;

        for (int f = 0; f < ic; f++)
        {
            double t2 = DIRECT_A2D_ELEM(tau2.data, f, k);

            if (t2 >= 0.0)
            {
                ta += t2;
                ts += 1.0;
            }
        }

        if (ts > 0.0)
        {
            tauAvg[k] = ta/ts;
        }

        for (int f = 0; f < ic; f++)
        {
            double t2 = DIRECT_A2D_ELEM(tau2.data, f, k);

            if (t2 >= 0.0)
            {
                DIRECT_A2D_ELEM(tauRel.data, f, k) = sqrt(t2)/tauAvg[k];
            }
            else
            {
                DIRECT_A2D_ELEM(tauRel.data, f, k) = 0.0;
            }
        }
    }

    for (int f = 0; f < ic; f++)
    {
        std::cout << f << ": ";
        BFactorFit bf(tauRel, weight, f, cutoff, Bw, Cw);

        std::vector<double> initial(2);
       // initial[0] = -0.001;
        //initial[1] = -10.0;
        initial[0] = -0.001/Bw;
        initial[1] = -10.0/Cw;

        std::vector<double> opt = NelderMead::optimize(initial, bf, 0.01, 0.000001, 1000000);

        out[f][0] = opt[0]*Bw;
        out[f][1] = opt[1]*Cw;

        std::cout << out[f][0] << ", " << out[f][1]  << "\n";
    }

    return out;
}

std::vector<d2Vector> FscHelper::fitBfactors(const Image<RFLOAT> &table, const Image<RFLOAT> &weight)
{
    const int kc = table.data.xdim;
    const int ic = table.data.ydim;

    std::vector<d2Vector> out(ic);
    std::vector<RFLOAT> avgFsc(kc, 0.0);

    for (int k = 0; k < kc; k++)
    {
        RFLOAT wsum = 0.0;

        for (int i = 0; i < ic; i++)
        {
            if (DIRECT_A2D_ELEM(table.data, i, k) < 0.0) continue;
            avgFsc[k] += DIRECT_A2D_ELEM(weight.data, i, k) * DIRECT_A2D_ELEM(table.data, i, k);
            wsum += DIRECT_A2D_ELEM(weight.data, i, k);
        }

        avgFsc[k] /= wsum;
    }

    const double eps = 1e-3;

    for (int i = 0; i < ic; i++)
    {
        d2Matrix A(0.0, 0.0, 0.0, 0.0);
        d2Vector b(0.0, 0.0);

        for (int k = 1; k < kc; k++)
        {
            RFLOAT fsc = DIRECT_A2D_ELEM(table.data, i, k);
            RFLOAT fscA = avgFsc[k];

            if (fsc < eps || fscA < eps) continue;

            RFLOAT t = sqrt((fsc - fsc * fscA) / (fscA - fsc * fscA));

            if (t < eps) continue;

            RFLOAT w = t * t * DIRECT_A2D_ELEM(weight.data, i, k);
            RFLOAT x = k * k;

            A(0,0) += w * x * x;
            A(1,0) += w * x;
            A(0,1) += w * x;
            A(1,1) += w;

            b[0] += w * x * log(t);
            b[1] += w * log(t);
        }

        A.invert();
        d2Vector mq = A * b;

        out[i][0] = -4.0 * mq[0];
        out[i][1] = exp(mq[1]);
    }

    return out;
}

Image<RFLOAT> FscHelper::tauRatio(const Image<RFLOAT> &table, const Image<RFLOAT> &weight)
{
    const int kc = table.data.xdim;
    const int ic = table.data.ydim;

    Image<RFLOAT> out(kc,ic);
    std::vector<RFLOAT> avgFsc(kc, 0.0);

    for (int k = 0; k < kc; k++)
    {
        RFLOAT wsum = 0.0;

        for (int i = 0; i < ic; i++)
        {
            if (DIRECT_A2D_ELEM(table.data, i, k) < 0.0) continue;

            avgFsc[k] += DIRECT_A2D_ELEM(weight.data, i, k) * DIRECT_A2D_ELEM(table.data, i, k);
            wsum += DIRECT_A2D_ELEM(weight.data, i, k);
        }

        avgFsc[k] /= wsum;
    }

    const double eps = 1e-3;

    for (int i = 0; i < ic; i++)
    {
        for (int k = 0; k < kc; k++)
        {
            RFLOAT fsc = DIRECT_A2D_ELEM(table.data, i, k);
            RFLOAT fscA = avgFsc[k];

            if (fsc < eps || fscA < eps) continue;

            DIRECT_A2D_ELEM(out.data, i, k) = (fsc - fsc * fscA) / (fscA - fsc * fscA);
        }
    }

    return out;
}

void FscHelper::computeBfactors(const std::vector<d2Vector>& bfacs, Image<RFLOAT> &table)
{
    const int kc = table.data.xdim;
    const int ic = table.data.ydim;

    for (int i = 0; i < ic; i++)
    for (int k = 0; k < kc; k++)
    {
        const double Bf = bfacs[i][0];
        const double Cf = bfacs[i][1];

        DIRECT_A2D_ELEM(table.data, i, k) = exp(Bf*k*k/4.0 + Cf);
    }
}

std::vector<double> FscHelper::powerSpectrum3D(const Image<Complex> &img)
{
    const int sh = img.data.xdim;
    const int s = img.data.ydim;

    std::vector<double> sum(sh, 0.0), wgh(sh, 0.0);

    for (int z = 0; z < img.data.zdim; z++)
    for (int y = 0; y < img.data.ydim; y++)
    for (int x = 0; x < img.data.xdim; x++)
    {
        const double xx = x;
        const double yy = y < sh? y : y - s;
        const double zz = z < sh? z : z - s;

        const double r = sqrt(xx*xx + yy*yy + zz*zz);
        const int ri = ROUND(r);

        if (ri < sh)
        {
            sum[ri] += DIRECT_NZYX_ELEM(img.data, 0, z, y, x).norm();
            wgh[ri] += 1.0;
        }
    }

    for (int i = 0; i < sh; i++)
    {
        if (wgh[i] > 0.0)
        {
            sum[i] /= wgh[i];
        }
    }

    return sum;
}

BFactorFit::BFactorFit(const Image<RFLOAT> &tau2, const Image<RFLOAT> &weight, int frame, int cutoff,
                       double Bscale, double Cscale)
:   tau2(tau2),
    weight(weight),
    frame(frame),
    cutoff(cutoff),
    Bscale(Bscale),
    Cscale(Cscale)
{
}

double BFactorFit::f(const std::vector<double>& x, void* tempStorage) const
{
    double Bf = x[0]*Bscale;
    double Cf = x[1]*Cscale;

    double sum = 0.0;

    const int kc = tau2.data.xdim;

    for (int k = cutoff; k < kc; k++)
    {
        double w = DIRECT_A2D_ELEM(weight.data, frame, k);
        RFLOAT t2 = DIRECT_A2D_ELEM(tau2.data, frame, k);
        double pv = exp(Bf*k*k/4.0 + Cf);

        double e = pv - t2;

        sum += w * e*e;
    }

    return sum;
}


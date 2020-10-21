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
#include <src/jaz/single_particle/damage_helper.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/single_particle/image_log.h>

using namespace gravis;

void DamageHelper::fitDamage(const Image<RFLOAT> &frcData, const Image<RFLOAT> &frcWeight,
                          std::vector<double> &amp, std::vector<double> &dec, int t0, double dosePerFrame, bool root)
{
    const int kc = frcData.data.xdim;
    const int tc = frcData.data.ydim;

    amp.resize(kc);
    dec.resize(kc);

    amp[0] = 0;
    dec[0] = 1;

    std::vector<double> snrData(tc), snrWeight(tc);

    const double eps = 1e-3;

    for (int k = 1; k < kc; k++)
    {
        double wSum = 0.0;
        double fSum = 0.0;

        for (int t = 0; t < tc; t++)
        {
            double frc = DIRECT_A2D_ELEM(frcData.data, t, k);
            if (frc > 1.0 - eps) frc = 1.0 - eps;
            if (frc < 0.0) frc = 0.0;

            if (root)
            {
                snrData[t] = sqrt(frc/(1.0 - frc));
            }
            else
            {
                snrData[t] = 2.0 * frc/(1.0 - frc);
            }

            snrWeight[t] = DIRECT_A2D_ELEM(frcWeight.data, t, k);

            wSum += snrWeight[t];
            fSum += snrData[t];
        }

        if (wSum == 0.0)
        {
            amp[k] = 0.0;
            dec[k] = 0.0;
        }
        else
        {
            DamageFit df(snrData, snrWeight, k, t0, 1.0, 1.0);

            std::vector<double> initial(2);
            initial[0] = 1.0;
            initial[1] = 1.0;

            std::vector<double> opt = NelderMead::optimize(initial, df, 0.01, 0.000001, 100000);

            std::cout << k << ": " << opt[0] << ", " << dosePerFrame * opt[1] << "\n";

            amp[k] = opt[0];
            dec[k] = dosePerFrame * opt[1];
        }
    }
}

Image<RFLOAT> DamageHelper::plotDamage(const std::vector<double> &amp, const std::vector<double> &dec, int frames)
{
    Image<RFLOAT> out(amp.size(), frames);

    const int kc = amp.size();
    const int tc = frames;

    for (int t = 0; t < tc; t++)
    for (int k = 0; k < kc; k++)
    {
        const double pred = amp[k] * exp(-t/dec[k]);
        DIRECT_A2D_ELEM(out.data, t, k) = pred;
    }

    return out;
}

Image<RFLOAT> DamageHelper::plotDamage(const std::vector<double> &dec, int frames)
{
    std::vector<double> amp(dec.size(), 1.0);
    return plotDamage(amp, dec, frames);
}

Image<RFLOAT> DamageHelper::plotDamageFrc(const std::vector<double> &amp, const std::vector<double> &dec, int frames)
{
    Image<RFLOAT> out(amp.size(), frames);

    const int kc = amp.size();
    const int tc = frames;

    for (int t = 0; t < tc; t++)
    for (int k = 0; k < kc; k++)
    {
        const double pred = amp[k] * exp(-t/dec[k]);
        DIRECT_A2D_ELEM(out.data, t, k) = pred / (pred + 2.0);
    }

    return out;
}

Image<RFLOAT> DamageHelper::plotDamageFrc(const std::vector<double> &dec, int frames)
{
    std::vector<double> amp(dec.size(), 1.0);
    return plotDamageFrc(amp, dec, frames);
}

void DamageHelper::fitGlobalDamage(const Image<RFLOAT> &frcData,
                                const Image<RFLOAT> &frcWeight,
                                std::vector<double> &amp,
                                double *a, double *b, double *c,
                                int k0, int k1, int t0,
                                double angpix, double dosePerFrame, bool L1)
{
    const int kc = frcData.data.xdim;
    const int tc = frcData.data.ydim;

    amp.resize(kc);
    amp[0] = 0.0;

    Image<RFLOAT> snrData(kc,tc);

    const double eps = 1e-3;

    for (int k = 1; k < kc; k++)
    for (int t = 0; t < tc; t++)
    {
        double frc = DIRECT_A2D_ELEM(frcData.data, t, k);
        if (frc > 1.0 - eps) frc = 1.0 - eps;
        if (frc < 0.0) frc = 0.0;

        DIRECT_A2D_ELEM(snrData.data, t, k) = 2.0 * frc/(1.0 - frc);
        //DIRECT_A2D_ELEM(snrData.data, t, k) = frc;
    }

    GlobalDamageFit gdf(snrData, frcWeight, k0, k1, t0, L1);

    std::vector<double> initial(3);
    initial[0] = 100.0;
    initial[1] = -1.0;
    initial[2] = 0.0;

    std::vector<double> opt = NelderMead::optimize(initial, gdf, 0.001, 0.000001, 1000000);

    const double rho = 1.0 / (2.0 * (kc-1) * angpix);

    *a = dosePerFrame * opt[0] / pow(rho, opt[1]);
    *b = opt[1];
    *c = dosePerFrame * opt[2];

    std::cout << (*a) << ", " << (*b) << ", " << (*c) << "\n";

    const double epsT = 1e-20;

    for (int k = 1; k < kc; k++)
    {
        double tau = opt[0] * pow(k,opt[1]) + opt[2];
        if (tau < epsT) tau = epsT;

        amp[k] = gdf.getScale(k, tau);
    }
}

Image<RFLOAT> DamageHelper::plotGlobalDamage(double a, double b, double c,
                                          const std::vector<double> &amp,
                                          int freqs, int frames,
                                          double angpix, double dosePerFrame, bool frc)
{
    Image<RFLOAT> out(freqs, frames);

    const int kc = freqs;
    const int tc = frames;

    const double eps = 1e-20;
    const double rho = 1.0 / (2.0 * (kc-1) * angpix);

    for (int t = 0; t < tc; t++)
    for (int k = 0; k < kc; k++)
    {
        double tau = a * pow(rho*k,b) + c;
        if (tau < eps) tau = eps;

        const double pred = amp[k] * exp(-t*dosePerFrame/tau);

        if (frc)
        {
            DIRECT_A2D_ELEM(out.data, t, k) = pred / (pred + 2.0);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, t, k) = pred;
        }
    }

    return out;
}

Image<RFLOAT> DamageHelper::plotGlobalDamage(double a, double b, double c, int freqs, int frames,
                                          double angpix, double dosePerFrame, bool frc)
{
    std::vector<double> amp(freqs, 1.0);
    return plotGlobalDamage(a, b, c, amp, freqs, frames, angpix, dosePerFrame, frc);
}

std::vector<Image<RFLOAT> > DamageHelper::damageWeights(
        int s, RFLOAT angpix,
        int f0, int fc, RFLOAT dosePerFrame,
        RFLOAT a, RFLOAT b, RFLOAT c)
{
    std::vector<Image<RFLOAT> > out(fc);

    for (long f = 0; f < fc; f++)
    {
        out[f] = damageWeight(s, angpix, (f+f0)*dosePerFrame, a, b, c);
    }

    return out;
}

Image<RFLOAT> DamageHelper::damageWeight(
        int s, RFLOAT angpix,
        RFLOAT dose,
        RFLOAT a, RFLOAT b, RFLOAT c)
{
    const int kc = s/2 + 1;

    const double rho = 1.0 / (2.0 * (kc-1) * angpix);

    Image<RFLOAT> out(kc, s);

    for (long y = 0; y < s; y++)
    for (long x = 0; x < kc; x++)
    {
        const double xf = x;
        const double yf = y < kc? y : (y - s);

        double k = sqrt(xf*xf + yf*yf);

        double tau = a * pow(rho*k, b) + c;

        DIRECT_A2D_ELEM(out.data, y, x) = exp(-dose/tau);
    }

    return out;
}

RFLOAT DamageHelper::damage(double k, int kc, RFLOAT angpix, RFLOAT dose, RFLOAT a, RFLOAT b, RFLOAT c)
{
    const double rho = 1.0 / (2.0 * (kc-1) * angpix);
    const double tau = a * pow(rho*k, b) + c;

    return exp(-dose/tau);
}

std::vector<double> DamageHelper::fitBFactors(const Image<RFLOAT> &fcc, int k0, int k1, int verb)
{
    const int kc = fcc.data.xdim;
    const int fc = fcc.data.ydim;

    double maxSumf = 0.0;
    int bestF = 0;

    for (int f = 0; f < fc; f++)
    {
        double sumf = 0.0;

        for (int k = 0; k < kc; k++)
        {
            sumf += fcc(f,k);
        }

        if (sumf > maxSumf)
        {
            maxSumf = sumf;
            bestF = f;
        }
    }

    if (verb > 0)
    	std::cout << "best f: " << bestF << "\n";

    std::vector<double> scale(kc, 0.0);

    for (int k = 0; k < kc; k++)
    {
        scale[k] = fcc(bestF,k);
    }

    std::vector<double> sig(fc);

    for (int it = 0; it < 5; it++)
    {
        for (int f = 0; f < fc; f++)
        {
            sig[f] = findSigmaRec(fcc, f, scale, 1, 10*kc, 20, 4, 0.1);
        }

        for (int k = 0; k < kc; k++)
        {
            double num = 0.0, denom = 0.0;

            for (int f = 0; f < fc; f++)
            {
                double p = fcc(f,k);
                double q = exp(-0.5 * k * k / (sig[f] * sig[f]));

                num += q*p;
                denom += q*q;
            }

            const double eps = 1e-20;
            scale[k] = denom > eps? num / denom : num / eps;
        }
    }

    Image<RFLOAT> debug(kc,fc);

    for (int k = 0; k < kc; k++)
    for (int f = 0; f < fc; f++)
    {
        debug(f,k) = scale[k] * exp(-0.5 * k * k / (sig[f] * sig[f]));
    }

    ImageLog::write(debug, "bfacs/debug");

    return sig;
}

std::pair<std::vector<d2Vector>,std::vector<double>> DamageHelper::fitBkFactors(
        const Image<RFLOAT> &fcc, int k0, int k1, int verb)
{
    const int kc = fcc.data.xdim;
    const int fc = fcc.data.ydim;

    double maxSumf = 0.0;
    int bestF = 0;

    for (int f = 0; f < fc; f++)
    {
        double sumf = 0.0;

        for (int k = 0; k < kc; k++)
        {
            sumf += fcc(f,k);
        }

        if (sumf > maxSumf)
        {
            maxSumf = sumf;
            bestF = f;
        }
    }

    if (verb > 0)
    	std::cout << "best f: " << bestF << "\n";

    std::vector<double> scale(kc, 0.0), wgh(kc, 1.0);

    for (int k = 0; k < kc; k++)
    {
        scale[k] = XMIPP_MAX(0.0, fcc(bestF,k));
    }

    std::vector<d2Vector> sig(fc);

    for (int it = 0; it < 5; it++)
    {
        for (int f = 0; f < fc; f++)
        {
            sig[f] = findSigmaKRec(fcc, f, scale, wgh, k0, k1, 1, 10*kc, 20, 4, 0.1);
        }

        for (int k = 0; k < kc; k++)
        {
            double num = 0.0, denom = 0.0;

            for (int f = 0; f < fc; f++)
            {
                double p = fcc(f,k);
                double q = sig[f].y * exp(-0.5 * k * k / (sig[f].x * sig[f].x));

                num += q*p;
                denom += q*q;
            }

            const double eps = 1e-20;
            scale[k] = denom > eps? num / denom : num / eps;
            scale[k] = XMIPP_MAX(0.0, scale[k]);
        }
    }

    return std::make_pair(sig,scale);
}

std::vector<d2Vector> DamageHelper::fitBkFactors(
        const Image<RFLOAT> &fcc,
        const Image<RFLOAT> &env,
        const Image<RFLOAT> &wgh,
        int k0, int k1)
{
    const int kc = fcc.data.xdim;
    const int fc = fcc.data.ydim;

    std::vector<double> envelope(kc), weight(kc);
    std::vector<d2Vector> sig(fc);

    for (int f = 0; f < fc; f++)
    {
        for (int k = 0; k < kc; k++)
        {
            envelope[k] = XMIPP_MAX(0.0, env(f,k));
            weight[k] = XMIPP_MAX(0.0, wgh(f,k));
        }

        sig[f] = findSigmaKRec(fcc, f, envelope, weight, k0, k1, 1, 10*kc, 20, 4, 0.1);
    }

    return sig;
}

Image<RFLOAT> DamageHelper::renderBkFit(
        const std::pair<std::vector<d2Vector>,std::vector<double>>& sigScale,
        int kc, int fc, bool noScale)
{
    Image<RFLOAT> out(kc,fc);

    for (int k = 0; k < kc; k++)
    for (int f = 0; f < fc; f++)
    {
        const double sigma = sigScale.first[f].x;
        const double a = sigScale.first[f].y;
        const double scale = noScale? 1.0 : sigScale.second[k];

        out(f,k) = scale * a * exp(-0.5 * k * k / (sigma * sigma));
    }

    return out;
}

Image<RFLOAT> DamageHelper::renderBkFit(std::vector<d2Vector> sig, int kc, int fc)
{
    Image<RFLOAT> out(kc,fc);

    for (int k = 0; k < kc; k++)
    for (int f = 0; f < fc; f++)
    {
        const double sigma = sig[f].x;
        const double a = sig[f].y;

        out(f,k) = a * exp(-0.5 * k * k / (sigma * sigma));
    }

    return out;
}

double DamageHelper::findSigmaRec(
        const Image<RFLOAT> &fcc, int f,
        const std::vector<double> &scale,
        double sig0, double sig1,
        int steps, int depth, double q)
{
    const int kc = fcc.data.xdim;

    double minErr = std::numeric_limits<double>::max();
    double bestSig = sig0;

    const double eps = 1e-20;

    for (int s = 0; s < steps; s++)
    {
        const double sig = sig0 + s*(sig1 - sig0)/(steps-1);
        const double sig2 = sig*sig;

        if (sig2 < eps) continue;

        double sum = 0.0;

        for (int k = 0; k < kc; k++)
        {
            const double d = fcc(f,k) - scale[k] * exp(-0.5 * k * k / sig2);
            sum += d*d;
        }

        if (sum < minErr)
        {
            minErr = sum;
            bestSig = sig;
        }
    }

    if (depth > 0)
    {
        const double hrange = 0.5 * (sig1 - sig0);
        double snext0 = bestSig - q*hrange;
        double snext1 = bestSig + q*hrange;

        if (snext0 < eps) snext0 = eps;

        return findSigmaRec(fcc, f, scale, snext0, snext1, steps, depth - 1, q);
    }

    return bestSig;
}


d2Vector DamageHelper::findSigmaKRec(
        const Image<RFLOAT> &fcc, int f,
        const std::vector<double> &envelope,
        const std::vector<double> &weight,
        int k0, int k1,
        double sig0, double sig1,
        int steps, int depth, double q)
{
    double minErr = std::numeric_limits<double>::max();
    double bestSig = sig0;
    double bestA = 1.0;

    const double eps = 1e-20;

    for (int s = 0; s < steps; s++)
    {
        const double sig = sig0 + s*(sig1 - sig0)/(steps-1);
        const double sig2 = sig*sig;

        if (sig2 < eps) continue;

        // find a
        double num = 0.0, denom = 0.0;

        for (int k = k0; k < k1; k++)
        {
            double p = fcc(f,k);
            double q = envelope[k] * exp(-0.5 * k * k / sig2);

            num += q*p;
            denom += q*q;
        }

        const double eps = 1e-20;
        double a = denom > eps? num / denom : num / eps;

        double sum = 0.0;

        for (int k = k0; k < k1; k++)
        {
            const double d = fcc(f,k) - envelope[k] * a * exp(-0.5 * k * k / sig2);
            sum += weight[k] * d*d;
        }

        if (sum < minErr)
        {
            minErr = sum;
            bestSig = sig;
            bestA = a;
        }
    }

    if (depth > 0)
    {
        const double hrange = 0.5 * (sig1 - sig0);
        double snext0 = bestSig - q*hrange;
        double snext1 = bestSig + q*hrange;

        if (snext0 < eps) snext0 = eps;

        return findSigmaKRec(fcc, f, envelope, weight, k0, k1, snext0, snext1, steps, depth - 1, q);
    }

    return d2Vector(bestSig, bestA);
}

std::vector<Image<RFLOAT>> DamageHelper::computeWeights(
        const std::vector<double> &bFacs, int kc, bool normalize)
{
    const int fc = bFacs.size();
    const int kc2 = 2 * (kc - 1);

    std::vector<Image<RFLOAT>> out(fc);

    for (int f = 0; f < fc; f++)
    {
        out[f] = Image<RFLOAT>(kc,kc2);

        for (int y = 0; y < kc2; y++)
        for (int x = 0; x < kc; x++)
        {
            double yy = y < kc? y : y - kc2;
            double r2 = x*x + yy*yy;

            out[f](y,x) = exp(-0.5*r2/(bFacs[f]*bFacs[f]));
        }
    }

    if (normalize)
    {
        for (int y = 0; y < kc2; y++)
        for (int x = 0; x < kc; x++)
        {
            double sum = 0.0;

            for (int f = 0; f < fc; f++)
            {
                sum += out[f](y,x);
            }

            for (int f = 0; f < fc; f++)
            {
                out[f](y,x) /= sum;
            }
        }
    }

    return out;
}

std::vector<Image<RFLOAT>> DamageHelper::computeWeights(
        const std::vector<d2Vector> &bkFacs, int kc, bool normalize)
{
    const int fc = bkFacs.size();
    const int kc2 = 2 * (kc - 1);

    std::vector<Image<RFLOAT>> out(fc);

    for (int f = 0; f < fc; f++)
    {
        out[f] = Image<RFLOAT>(kc,kc2);

        for (int y = 0; y < kc2; y++)
        for (int x = 0; x < kc; x++)
        {
            double yy = y < kc? y : y - kc2;
            double r2 = x*x + yy*yy;

            out[f](y,x) = bkFacs[f].y * exp(-0.5*r2/(bkFacs[f].x*bkFacs[f].x));
        }
    }

    if (normalize)
    {
        for (int y = 0; y < kc2; y++)
        for (int x = 0; x < kc; x++)
        {
            double sum = 0.0;

            for (int f = 0; f < fc; f++)
            {
                sum += out[f](y,x);
            }

            for (int f = 0; f < fc; f++)
            {
                out[f](y,x) /= sum;
            }
        }
    }

    return out;
}

std::vector<Image<RFLOAT> > DamageHelper::computeWeights(
        const std::vector<d2Vector> &bkFacs,
        int kc, RFLOAT angpix, RFLOAT totalDose,
        RFLOAT dmga, RFLOAT dmgb, RFLOAT dmgc, bool normalize)
{
    const int fc = bkFacs.size();
    const int kc2 = 2 * (kc - 1);

    std::vector<Image<RFLOAT>> out(fc);

    for (int f = 0; f < fc; f++)
    {
        out[f] = Image<RFLOAT>(kc,kc2);

        const double dose = totalDose*f/fc;

        for (int y = 0; y < kc2; y++)
        for (int x = 0; x < kc; x++)
        {
            double yy = y < kc? y : y - kc2;
            double r2 = x*x + yy*yy;

            out[f](y,x) = bkFacs[f].y * exp(-0.5*r2/(bkFacs[f].x*bkFacs[f].x))
                    * damage(sqrt(r2), kc, angpix, dose, dmga, dmgb, dmgc);
        }
    }

    if (normalize)
    {
        for (int y = 0; y < kc2; y++)
        for (int x = 0; x < kc; x++)
        {
            double sum = 0.0;

            for (int f = 0; f < fc; f++)
            {
                sum += out[f](y,x);
            }

            if (sum > 0.0)
            {
                for (int f = 0; f < fc; f++)
                {
                    out[f](y,x) /= sum;
                }
            }
        }
    }

    return out;
}


std::vector<Image<RFLOAT>> DamageHelper::computeWeights(const Image<RFLOAT>& fcc)
{
    const int kc = fcc.data.xdim;
    const int fc = fcc.data.ydim;
    const int kc2 = 2 * (kc - 1);

    std::vector<Image<RFLOAT>> out(fc);

    for (int f = 0; f < fc; f++)
    {
        out[f] = Image<RFLOAT>(kc,kc2);

        for (int y = 0; y < kc2; y++)
        for (int x = 0; x < kc; x++)
        {
            const double yy = y < kc? y : y - kc2;
            const double r2 = x*x + yy*yy;
            const double r = sqrt(r2);

            int idx = (int)r;

            if (idx >= kc)
            {
                out[f](y,x) = 0.0;
            }
            else if (idx == kc-1)
            {
                out[f](y,x) = XMIPP_MAX(fcc(f, idx), 0.0);
            }
            else
            {
                const double rf = r - idx;
                const double v0 = XMIPP_MAX(fcc(f, idx), 0.0);
                const double v1 = XMIPP_MAX(fcc(f, idx+1), 0.0);

                out[f](y,x) = rf * v1 + (1.0 - rf) * v0;
            }
        }
    }

    for (int y = 0; y < kc2; y++)
    for (int x = 0; x < kc; x++)
    {
        double sum = 0.0;

        for (int f = 0; f < fc; f++)
        {
            sum += out[f](y,x);
        }

        if (sum > 0.0)
        {
            for (int f = 0; f < fc; f++)
            {
                out[f](y,x) /= sum;
            }
        }
    }

    return out;
}

DamageFit::DamageFit(const std::vector<double> &snrData,
                     const std::vector<double> &snrWeight,
                     int k, int t0,
                     double ampScale, double decScale)
:   snrData(snrData),
    snrWeight(snrWeight),
    k(k), t0(t0),
    ampScale(ampScale), decScale(decScale)
{}

double DamageFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int tc = snrData.size();
    const double amp = ampScale * x[0];
    const double dec = decScale * x[1];

    double sum = 0.0;

    for (int t = t0; t < tc; t++)
    {
        const double pred = amp * exp(-t/dec);
        const double obs = snrData[t];
        const double d = pred - obs;

        sum += snrWeight[t] * d*d;
    }

    return sum;
}

GlobalDamageFit::GlobalDamageFit(const Image<RFLOAT> &snrData, const Image<RFLOAT> &snrWeight, int k0, int k1, int t0, bool L1)
: snrData(snrData),
  snrWeight(snrWeight),
  k0(k0),
  k1(k1),
  t0(t0),
  L1(L1)
{}

double GlobalDamageFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int tc = snrData.data.ydim;

    const double a = x[0];
    const double b = x[1];
    const double c = x[2];

    double sum = 0.0;

    const double eps = 1e-20;

    for (int k = k0; k < k1; k++)
    {
        double tau = a * pow(k,b) + c;
        if (tau < eps) tau = eps;

        const double n = getScale(k, tau);

        for (int t = t0; t < tc; t++)
        {
            const double d = exp(-t/tau);
            const double e = DIRECT_A2D_ELEM(snrData.data, t, k) - n * d;

            if (L1)
            {
                sum += DIRECT_A2D_ELEM(snrWeight.data, t, k) * std::abs(e);
            }
            else
            {
                sum += DIRECT_A2D_ELEM(snrWeight.data, t, k) * e * e;
            }
        }
    }

    return sum;
}

double GlobalDamageFit::getScale(int k, double tau) const
{
    double nn = 0.0;
    double nd = 0.0;

    const int tc = snrData.data.ydim;

    for (int t = t0; t < tc; t++)
    {
        const double d = exp(-t/tau);
        const double y = DIRECT_A2D_ELEM(snrData.data, t, k);
        const double w = DIRECT_A2D_ELEM(snrWeight.data, t, k);

        nn += w * y * d;
        nd += w * d * d;
    }

    return nd > 0.0? nn/nd : 0.0;
}

PerFrameBFactorFit::PerFrameBFactorFit(const Image<RFLOAT> &fcc, int k0, int k1)
:   fcc(fcc),
    kc(fcc.data.xdim),
    fc(fcc.data.ydim),
    k0(k0), k1(k1)
{
}

double PerFrameBFactorFit::f(const std::vector<double> &x) const
{
    std::vector<double> scale(k1-k0+1), bfac(fc);

    for (int f = 0; f < fc; f++)
    {
        bfac[f] = x[f];
    }

    for (int k = 0; k < k1-k0+1; k++)
    {
        scale[k] = x[k+fc];
    }

    double sum = 0.0;

    for (int f = 0; f < fc; f++)
    for (int k = k0; k <= k1; k++)
    {
        double d = fcc(f,k) - scale[k-k0] * exp(-0.5*f*f/(bfac[f]*bfac[f]));

        sum += d*d;
    }

    return sum;
}

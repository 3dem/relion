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

#ifndef DAMAGE_FIT_H
#define DAMAGE_FIT_H

#include <src/ctf.h>
#include <src/image.h>
#include <vector>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/optimization/optimization.h>

class DamageFit : public Optimization
{
    public:

        DamageFit(const std::vector<double>& snrData,
                  const std::vector<double>& snrWeight,
                  int k, int t0,
                  double ampScale, double decScale);

        double f(const std::vector<double>& x, void* tempStorage) const;

    private:

        const std::vector<double>& snrData;
        const std::vector<double>& snrWeight;
        int k, t0;
        double ampScale, decScale;
};

class GlobalDamageFit : public Optimization
{
    public:

        GlobalDamageFit(const Image<RFLOAT>& snrData,
                  const Image<RFLOAT>& snrWeight,
                  int k0, int k1, int t0, bool L1);

        double f(const std::vector<double>& x, void* tempStorage) const;

        double getScale(int k, double tau) const;

    private:

        const Image<RFLOAT>& snrData;
        const Image<RFLOAT>& snrWeight;
        int k0, k1, t0;
        bool L1;
};

class PerFrameBFactorFit : public Optimization
{
    public:

        PerFrameBFactorFit(const Image<RFLOAT>& fcc, int k0, int k1);

        double f(const std::vector<double>& x) const;

    private:

        const Image<RFLOAT> &fcc;
        int kc, fc, k0, k1;
};

class DamageHelper
{
    public:

        static void fitDamage(const Image<RFLOAT>& frcData,
                const Image<RFLOAT>& frcWeight,
                std::vector<double>& amp,
                std::vector<double>& dec,
                int t0, double dosePerFrame, bool root = false);

        static Image<RFLOAT> plotDamage(
                const std::vector<double>& amp,
                const std::vector<double>& dec,
                int frames);

        static Image<RFLOAT> plotDamage(
                const std::vector<double>& dec,
                int frames);

        static Image<RFLOAT> plotDamageFrc(
                const std::vector<double>& amp,
                const std::vector<double>& dec,
                int frames);

        static Image<RFLOAT> plotDamageFrc(
                const std::vector<double>& dec,
                int frames);


        static void fitGlobalDamage(
                const Image<RFLOAT>& frcData,
                const Image<RFLOAT>& frcWeight,
                std::vector<double>& amp,
                double* a, double* b, double* c,
                int k0, int k1, int t0, double angpix,
                double dosePerFrame, bool L1 = false);

        static Image<RFLOAT> plotGlobalDamage(
                double a, double b, double c,
                const std::vector<double>& amp,
                int freqs, int frames, double angpix, double dosePerFrame,
                bool frc = false);

        static Image<RFLOAT> plotGlobalDamage(
                double a, double b, double c,
                int freqs, int frames, double angpix, double dosePerFrame,
                bool frc = false);

        static std::vector<Image<RFLOAT>> damageWeights(
                int s, RFLOAT angpix,
                int f0, int fc, RFLOAT dosePerFrame,
                RFLOAT a, RFLOAT b, RFLOAT c);

        static Image<RFLOAT> damageWeight(
                int s, RFLOAT angpix, RFLOAT dose,
                RFLOAT a, RFLOAT b, RFLOAT c);

        static RFLOAT damage(
                double k, int kc, RFLOAT angpix, RFLOAT dose,
                RFLOAT a, RFLOAT b, RFLOAT c);

        static std::vector<double> fitBFactors(const Image<RFLOAT>& fcc, int k0, int k1, int verb = 0);

        static std::pair<std::vector<gravis::d2Vector>,std::vector<double>>
            fitBkFactors(const Image<RFLOAT>& fcc, int k0, int k1, int verb = 0);

        static std::vector<gravis::d2Vector>
            fitBkFactors(const Image<RFLOAT>& fcc,
                         const Image<RFLOAT>& env,
                         const Image<RFLOAT>& wgh,
                         int k0, int k1);

        static Image<RFLOAT> renderBkFit(
                const std::pair<std::vector<gravis::d2Vector>,std::vector<double>>& sigScale,
                int kc, int fc, bool noScale = false);

        static Image<RFLOAT> renderBkFit(
                std::vector<gravis::d2Vector> sig,
                int kc, int fc);

        static double findSigmaRec(const Image<RFLOAT>& fcc, int f,
                                   const std::vector<double>& scale,
                                   double sig0, double sig1,
                                   int steps, int depth, double q);

        static gravis::d2Vector findSigmaKRec(
                const Image<RFLOAT>& fcc, int f,
                const std::vector<double>& envelope,
                const std::vector<double>& weight,
                int k0, int k1,
                double sig0, double sig1,
                int steps, int depth, double q);

        static std::vector<Image<RFLOAT>> computeWeights(
                const std::vector<double>& bFacs, int kc, bool normalize = true);

        static std::vector<Image<RFLOAT>> computeWeights(
                const std::vector<gravis::d2Vector>& bkFacs, int kc, bool normalize = true);

        static std::vector<Image<RFLOAT>> computeWeights(
                const std::vector<gravis::d2Vector>& bkFacs,
                int kc, RFLOAT angpix, RFLOAT totalDose, RFLOAT dmga, RFLOAT dmgb, RFLOAT dmgc,
                bool normalize = true);

        static std::vector<Image<RFLOAT>> computeWeights(
                const Image<RFLOAT>& fcc);
};

#endif

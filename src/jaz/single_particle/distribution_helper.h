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

#ifndef DISTRIBUTION_HELPER_H
#define DISTRIBUTION_HELPER_H

#include <vector>
#include <src/image.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t4Vector.h>
#include <src/jaz/optimization/optimization.h>

class CorruptedGaussFit : public Optimization
{
    public:

        CorruptedGaussFit(const std::vector<double>& values, double area);

            const std::vector<double>& values;
            const double area;

        double f(const std::vector<double>& x, void* tempStorage) const;
};

class CorruptedExponentialFit : public Optimization
{
    public:

        CorruptedExponentialFit(const std::vector<double>& velocities, double area);

            const std::vector<double>& velocities;
            const double area;

        double f(const std::vector<double>& x, void* tempStorage) const;
};

class BivariateStudentFit : public Optimization
{
    public:

        BivariateStudentFit(const std::vector<double>& velocities, double fixedNu);

            const std::vector<double>& values;
            double fixedNu;

        double f(const std::vector<double>& x, void* tempStorage) const;
};

class GaussStudentFit : public Optimization
{
    public:

        GaussStudentFit(const std::vector<double>& values, double fixedNu);

            const std::vector<double>& values;
            double fixedNu;

        double f(const std::vector<double>& x, void* tempStorage) const;
};

class DoubleStudentFit : public Optimization
{
    public:

        DoubleStudentFit(const std::vector<double>& values,
                         const std::vector<int>& counts,
                         double nu0,
                         bool fixedNu0,
                         double nuS,
                         bool fixedNuS,
                         double area,
                         bool relTheta);

            const std::vector<double>& values;
            const std::vector<int>& counts;
            double myNu0, myNuS, area;
            bool fixedNu0, fixedNuS, relTheta;

        double f(const std::vector<double>& x, void* tempStorage) const;
};

class DoubleGaussFit : public Optimization
{
    public:

        DoubleGaussFit(const std::vector<double>& values,
                       const std::vector<int>& counts,
                       const std::vector<int>& totalCounts,
                       double area,
                       bool centered);

            const std::vector<double>& values;
            const std::vector<int>& counts, &totalCounts;
            double area;
            bool centered;

        double f(const std::vector<double>& x, void* tempStorage) const;
};

class MultiStudentFit : public Optimization
{
    public:

        MultiStudentFit( const std::vector<double>& values,
                         const std::vector<int>& counts,
                         int minPN, int maxPN,
                         double area);

            const std::vector<double>& values;
            const std::vector<int>& counts;
            int minPN, maxPN;
            double area;

        double f(const std::vector<double>& x, void* tempStorage) const;
};



class CollectiveGaussStudentFit : public Optimization
{
    public:

        CollectiveGaussStudentFit(const std::vector<std::vector<double>>& values, double fixedNu,
                        double mGauss, double f0Gauss, double betaStudent);

            const std::vector<std::vector<double>>& values;
            double fixedNu;
            double mGauss, f0Gauss, betaStudent;

        double f(const std::vector<double>& x, void* tempStorage) const;
};

class CorruptedBivariateStudentFit : public Optimization
{
    public:

        CorruptedBivariateStudentFit(const std::vector<double>& values, double area, double fixedNu);

            const std::vector<double>& values;
            const double area;
            const double fixedNu;

        double f(const std::vector<double>& x, void* tempStorage) const;
};

class DistributionHelper
{
    public:

        static double sampleGauss(double mu, double sigma);
        static Complex sampleGauss(Complex mu, double sigma);

        static std::vector<std::vector<double> > ringHistograms(
                const std::vector<std::vector<gravis::d2Vector> >& coords,
                int bins,
                double step);

        static Image<double> drawHistogramStack(
                const std::vector<std::vector<double> >& data,
                int height);

        static std::vector<double> normalizeRing(
                const std::vector<double>& histogram,
                double step);

        static std::vector<double> toEnergy(
                const std::vector<double>& probability);

        static std::vector<double> toProbability(
                const std::vector<double>& energy);

        static void writeRing(
                const std::vector<double>& values,
                double step,
                std::string filename,
                bool skipZero);

        static std::vector<double> corruptedExponentialProbability(
                double beta,
                double theta,
                double step,
                int bins,
                double area);

        static std::vector<double> corruptedGaussProbability(
                double sigma,
                double theta,
                double area,
                double step,
                int bins);

        static std::vector<double> bivariateStudentProbability(
                double sigma,
                double nu,
                double step,
                int bins);

        static std::vector<double> corruptedBivariateStudentProbability(
                double sigma,
                double nu,
                double theta,
                double rho,
                double step,
                int bins);

        static std::vector<double> bivariateGaussStudentProbability(
                double sigmaS,
                double sigmaG,
                double theta,
                double nu,
                double step,
                int bins);

        static std::vector<double> doubleStudentProbability(
                double sigma0,
                double sigmaS,
                double thetaS,
                double nu0,
                double nuS,
                double N,
                double area,
                double step,
                int bins);

        static std::vector<double> doubleGaussProbability(
                double sigma0,
                double sigmaS,
                double thetaS,
                double theta0,
                double thetaMin,
                bool relTheta,
                bool erfTheta,
                double N,
                double area,
                double step,
                int bins);

        static std::vector<double> multiStudentProbability(
                std::vector<double> params,
                int N, int minPC, int maxPC,
                double area,
                double step,
                int bins);

        static gravis::d2Vector findCorrupted2DExponential(
                const std::vector<std::vector<gravis::d2Vector> >& coords,
                int f,
                double betaMax,
                double thetaMax,
                int betaSteps,
                int thetaSteps,
                double area,
                std::string outFilename);

        static gravis::d2Vector fitCorrupted2DExponential(
                const std::vector<std::vector<gravis::d2Vector> >& coords,
                int f,
                double area);

        static gravis::d2Vector fitBivariateStudent(
                const std::vector<std::vector<gravis::d2Vector> >& coords,
                int f, double fixedNu = -1);

        static gravis::d3Vector fitCorruptedBivariateStudent(
                const std::vector<std::vector<gravis::d2Vector> >& coords,
                int f, double area, double fixedNu = -1);

        static std::vector<double> fitDoubleStudent(
                const std::vector<std::vector<gravis::d2Vector> >& coords,
                const std::vector<int>& counts,
                int f, double area,
                bool fixedNu0, double nu0,
                bool fixedNuS, double nuS, bool relTheta);

        static std::vector<double> fitDoubleGauss(const std::vector<std::vector<gravis::d2Vector> >& coords,
                const std::vector<int>& counts, const std::vector<int> &totalCounts,
                int f, double area, bool centered, std::vector<double> initial);

        static std::vector<double> fitCorruptedGauss(
                const std::vector<std::vector<gravis::d2Vector> >& coords,
                int f, double area);

        static std::vector<double> fitMultiStudent(
                const std::vector<std::vector<gravis::d2Vector> >& coords,
                const std::vector<int>& counts,
                int f, double area,
                int minPC, int maxPC);

        static gravis::d4Vector fitGaussStudent(const std::vector<std::vector<gravis::d2Vector> >& coords,
                int f, double fixedNu = -1);


};

#endif

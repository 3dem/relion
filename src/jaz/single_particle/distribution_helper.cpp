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

#include <src/jaz/single_particle/distribution_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/optimization/nelder_mead.h>

using namespace gravis;

double DistributionHelper::sampleGauss(double mu, double sigma)
{
    double u1 = rand()/(double)RAND_MAX;
    double u2 = rand()/(double)RAND_MAX;

    double z = sqrt(-2.0*log(u1)) * cos(2.0*PI*u2);

    return mu + sigma * z;
}

Complex DistributionHelper::sampleGauss(Complex mu, double sigma)
{
    double u1r = rand()/(double)RAND_MAX;
    double u2r = rand()/(double)RAND_MAX;
    double u1i = rand()/(double)RAND_MAX;
    double u2i = rand()/(double)RAND_MAX;

    double xr = sqrt(-2.0*log(u1r)) * cos(2.0*PI*u2r);
    double xi = sqrt(-2.0*log(u1i)) * cos(2.0*PI*u2i);

    return mu + sigma * Complex(xr,xi);

}

std::vector<std::vector<double> > DistributionHelper::ringHistograms(
        const std::vector<std::vector<gravis::d2Vector> >& coords, int bins, double step)
{
    int pc = coords.size();
    int fc = coords[0].size();

    std::vector<std::vector<double> > out(fc);

    for (int f = 0; f < fc; f++)
    {
        out[f] = std::vector<double>(bins,0.0);
    }

    for (int f = 0; f < fc-1; f++)
    {
        for (int p = 0; p < pc; p++)
        {
            d2Vector p0 = coords[p][f];
            d2Vector p1 = coords[p][f+1];

            d2Vector v = p1 - p0;

            int vbin = (int)(v.length()/step);

            if (vbin < bins)
            {
                double ringArea = 2*vbin + 1;
                out[f][vbin] += 1.0/ringArea;
            }
        }
    }

    return out;
}

Image<double> DistributionHelper::drawHistogramStack(const std::vector<std::vector<double> >& data, int height)
{
    const int fc = data.size();
    const int bins = data[0].size();

    Image<double> hist(bins,height,fc);
    hist.data.initZeros();
    hist.data.xinit = 0;
    hist.data.yinit = 0;
    hist.data.zinit = 0;

    double bmax = 0.0;

    for (int f = 0; f < fc; f++)
    {
        for (int b = 0; b < bins; b++)
        {
            if (data[f][b] > bmax) bmax = data[f][b];
        }
    }

    if (bmax == 0.0) return hist;

    double scale = (double)(height - 1)/bmax;

    for (int f = 0; f < fc; f++)
    {
        for (int b = 0; b < bins; b++)
        {
            double val = scale * data[f][b];

            int vi = (int)val;
            if (vi > height) vi = height;

            for (int y = 0; y < vi; y++)
            {
                hist(f,y,b) = 1.0;
            }

            if (vi < height-1)
            {
                hist(f,vi,b) = val - vi;
            }
        }
    }

    return hist;
}

std::vector<double> DistributionHelper::normalizeRing(const std::vector<double>& histogram, double step)
{
    const int s = histogram.size();
    std::vector<double> out(s,0.0);

    double sum = 0.0;

    for (int i = 0; i < s; i++)
    {
        double x0 = step*i;
        double x1 = step*(i+1);

        double a0 = 4.0*PI*x0*x0;
        double a1 = 4.0*PI*x1*x1;

        double ringArea = a1 - a0;

        sum += histogram[i]*ringArea;
    }

    if (sum > 0.0)
    {
        for (int i = 0; i < s; i++)
        {
            out[i] = histogram[i] / sum;
        }
    }

    return out;
}

std::vector<double> DistributionHelper::toEnergy(const std::vector<double>& probability)
{
    const int s = probability.size();
    std::vector<double> out(s,0.0);

    for (int i = 0; i < s; i++)
    {
        out[i] = -log(probability[i]);
    }

    return out;
}

std::vector<double> DistributionHelper::toProbability(const std::vector<double> &energy)
{
    const int s = energy.size();
    std::vector<double> out(s,0.0);

    for (int i = 0; i < s; i++)
    {
        out[i] = exp(-energy[i]);
    }

    return out;
}

void DistributionHelper::writeRing(const std::vector<double>& values, double step,
                                   std::string filename, bool skipZero)
{
    std::ofstream ofs(filename);

    const int s = values.size();

    for (int i = 0; i < s; i++)
    {
        if (!(values[i] == values[i])
                || std::isinf(values[i])
                || skipZero && values[i] == 0.0) continue;

        ofs << (i + 0.5)*step << " " << values[i] << "\n";
    }
}

std::vector<double> DistributionHelper::corruptedExponentialProbability(
        double beta, double theta, double step, int bins, double area)
{
    std::vector<double> out(bins);
    const double twoPiBetaSq = 2.0*PI*beta*beta;

    for (int i = 0; i < bins; i++)
    {
        double x = (i+0.5)*step;
        out[i] = theta/area + (1.0 - theta)*exp(-x/beta)/twoPiBetaSq;
    }

    return out;
}

std::vector<double> DistributionHelper::corruptedGaussProbability(
        double sigma, double theta, double area, double step, int bins)
{
    std::vector<double> out(bins);

    for (int i = 0; i < bins; i++)
    {
        double x = (i+0.5)*step;
        out[i] = theta/area + (1.0 - theta)*exp(-0.5*x*x/(sigma*sigma))/(2.0*PI*sigma*sigma);
    }

    return out;
}

std::vector<double> DistributionHelper::bivariateStudentProbability(double sigma, double nu, double step, int bins)
{
    std::vector<double> out(bins);
    const double twoPiSigmaSq = 2.0*PI*sigma*sigma;
    const double nuSigmaSq = nu*sigma*sigma;

    for (int i = 0; i < bins; i++)
    {
        double x = (i+0.5)*step;
        out[i] = (1.0/twoPiSigmaSq) * pow(1.0 + x*x/nuSigmaSq, -(1.0 + nu/2.0));
    }

    return out;
}

std::vector<double> DistributionHelper::corruptedBivariateStudentProbability(
        double sigma, double nu, double theta, double rho, double step, int bins)
{
    std::vector<double> out(bins);
    const double twoPiSigmaSq = 2.0*PI*sigma*sigma;
    const double nuSigmaSq = nu*sigma*sigma;

    for (int i = 0; i < bins; i++)
    {
        double x = (i+0.5)*step;
        out[i] = theta * rho + (1.0 - theta) * (1.0/twoPiSigmaSq) * pow(1.0 + x*x/nuSigmaSq, -(1.0 + nu/2.0));
    }

    return out;
}

std::vector<double> DistributionHelper::bivariateGaussStudentProbability(
        double sigmaS, double sigmaG, double theta, double nu, double step, int bins)
{
    std::vector<double> out(bins);

    const double twoPiSigmaSSq = 2.0*PI*sigmaS*sigmaS;
    const double nuSigmaSSq = nu*sigmaS*sigmaS;
    const double twoPiSigmaGSq = 2.0*PI*sigmaG*sigmaG;


    for (int i = 0; i < bins; i++)
    {
        double x = (i+0.5)*step;
        out[i] = theta * (1.0/(twoPiSigmaGSq)) * exp(-0.5*x*x/(sigmaG*sigmaG))
                  + (1.0 - theta) * (1.0/twoPiSigmaSSq) * pow(1.0 + x*x/nuSigmaSSq, -(1.0 + nu/2.0));
    }

    return out;
}

std::vector<double> DistributionHelper::doubleStudentProbability(
        double sigma0, double sigmaS, double thetaS, double nu0, double nuS, double N, double area,
        double step, int bins)
{
    std::vector<double> out(bins);

    double sig2, nu, theta;

    if (N > 0)
    {
        const double Ninv = 1.0 / N;
        sig2 = sigma0 + Ninv * sigmaS;
        nu = nu0 + Ninv * nuS;
        theta = thetaS*Ninv;
    }
    else
    {
        sig2 = sigma0;
        nu = nu0;
        theta = 0.0;
    }

    const double twoPiSigmaSq = 2.0*PI*sig2;
    const double nuSigmaSq = nu*sig2;

    const double rho = 1.0 / area;

    for (int i = 0; i < bins; i++)
    {
        double x = (i+0.5)*step;
        out[i] = theta * rho
                + (1.0 - theta) * (1.0/twoPiSigmaSq) * pow(1.0 + x*x/nuSigmaSq, -(1.0 + nu/2.0));
    }

    return out;
}

std::vector<double> DistributionHelper::doubleGaussProbability(
        double sigma0, double sigmaS, double thetaS, double theta0, double thetaMin,
        bool relTheta, bool erfTheta, double N, double area, double step, int bins)
{
    std::vector<double> out(bins);

    double sig2, theta;

    if (N > 0)
    {
        const double Ninv = 1.0 / N;
        sig2 = sigma0 + Ninv * sigmaS;

        if (relTheta)
        {
            if (erfTheta)
            {
                theta = 1.0 - exp(log(0.5*(1.0 + sqrt(1.0 - exp(-N*thetaS)))) * area);
            }
            else
            {
                theta = theta0 + thetaS*Ninv;

                if (theta < thetaMin)
                {
                    theta = thetaMin;
                }
                if (theta > 1.0)
                {
                    theta = 1.0;
                }
            }

        }
        else
        {
            theta = thetaS;
        }


        if (theta > 1.0) theta = 1.0;
    }
    else
    {
        sig2 = sigma0;
        theta = thetaMin;
    }

    const double twoPiSigmaSq = 2.0*PI*sig2;

    const double rho = 1.0 / area;

    for (int i = 0; i < bins; i++)
    {
        double x = (i+0.5)*step;
        out[i] = theta * rho
                + (1.0 - theta) * (1.0/twoPiSigmaSq) * exp(-0.5*x*x/sig2);
    }

    return out;

}

std::vector<double> DistributionHelper::multiStudentProbability(
        std::vector<double> params, int N, int minPC, int maxPC, double area, double step, int bins)
{
    std::vector<double> out(bins);

    const double sig0 = params[0];
    const double sigS = params[1];

    const int pbc = maxPC - minPC + 1;

    int pb;
    double sig2;

    if (N == 0)
    {
        pb = maxPC - minPC;
        sig2 = sig0;
    }
    else
    {
        pb = N - minPC;
        sig2 = sig0 + sigS/(double)N;
    }

    if (pb >= pbc) pb = pbc;
    if (pb < 0) pb = 0;

    const double nu = params[2 + 2*pb];
    const double theta = params[2 + 2*pb + 1];

    const double twoPiSigmaSq = 2.0*PI*sig2*nu;
    const double nuSigmaSq = nu*sig2;

    const double rho = 1.0 / area;

    for (int i = 0; i < bins; i++)
    {
        double x = (i+0.5)*step;
        out[i] = theta * rho
                + (1.0 - theta) * (1.0/twoPiSigmaSq) * pow(1.0 + x*x/nuSigmaSq, -(1.0 + nu/2.0));
    }

    return out;
}

gravis::d2Vector DistributionHelper::findCorrupted2DExponential(
        const std::vector<std::vector<gravis::d2Vector> > &coords,
        int f, double betaMax, double thetaMax, int betaSteps, int thetaSteps, double area,
        std::string outFilename)
{
    int pc = coords.size();

    double eMin = std::numeric_limits<double>::max();
    double bestBeta = 0.0, bestTheta = 0.0;

    Image<RFLOAT> eImg;

    bool drawPlot = outFilename != "";

    if (drawPlot)
    {
        eImg = Image<RFLOAT>(thetaSteps, betaSteps);
        eImg.data.xinit = 0;
        eImg.data.yinit = 0;
        eImg.data.zinit = 0;
    }

    for (int i = 0; i < betaSteps; i++)
    for (int j = 0; j < thetaSteps; j++)
    {
        double beta = betaMax*(i+1)/(double)betaSteps;
        double beta2 = beta*beta;
        double theta = thetaMax*(j+1)/(double)thetaSteps;

        double e = 0.0;

        for (int p = 0; p < pc; p++)
        {
            d2Vector p0 = coords[p][f];
            d2Vector p1 = coords[p][f+1];

            d2Vector v = p1 - p0;
            double l1 = v.length();

            e += -log(theta/area + (1.0 - theta)*exp(-l1/beta)/(2.0*PI*beta2));
        }

        if (e < eMin)
        {
            eMin = e;
            bestBeta = beta;
            bestTheta = theta;
        }

        if (drawPlot)
        {
            eImg(i,j) = e;
        }
    }

    if (drawPlot)
    {
        VtkHelper::writeVTK(eImg, outFilename,
                        0, 0, 0, 100*thetaMax/(double)thetaSteps,
                        betaMax/(double)betaSteps, 1.0);
    }

    return d2Vector(bestBeta, bestTheta);
}

gravis::d2Vector DistributionHelper::fitCorrupted2DExponential(
        const std::vector<std::vector<gravis::d2Vector> > &coords,
        int f, double area)
{
    int pc = coords.size();

    std::vector<double> velocities(pc);

    for (int p = 0; p < pc; p++)
    {
        if (coords[p].size() < f+2)
        {
            std::stringstream sts0;
            sts0 << p;
            std::stringstream sts1;
            sts1 << (f+1);
            std::stringstream sts2;
            sts2 << (coords[p].size());

            REPORT_ERROR("not enough data points for particle "+sts0.str()+"; required: "+sts1.str()+", present: "+sts2.str()+"\n");
        }
        velocities[p] = (coords[p][f+1] - coords[p][f]).length();
    }

    CorruptedExponentialFit cef(velocities, area);

    std::vector<double> initial = {1.0, 0.0};

    std::vector<double> params = NelderMead::optimize(initial, cef, 0.1, 0.001, 10000);

    return d2Vector(params[0], params[1]);
}

d2Vector DistributionHelper::fitBivariateStudent(
        const std::vector<std::vector<d2Vector> > &coords, int f, double fixedNu)
{
    int pc = coords.size();

    std::vector<double> velocities(pc);

    for (int p = 0; p < pc; p++)
    {
        velocities[p] = (coords[p][f+1] - coords[p][f]).length();
    }


    if (fixedNu > 0.0)
    {
        BivariateStudentFit cef(velocities, fixedNu);

        std::vector<double> initial = {1.0};
        std::vector<double> params = NelderMead::optimize(initial, cef, 0.1, 0.0001, 100000);
        return d2Vector(params[0], fixedNu);
    }
    else
    {
        BivariateStudentFit cef(velocities, -1);

        std::vector<double> initial = {1.0, 2.0};
        std::vector<double> params = NelderMead::optimize(initial, cef, 0.1, 0.0001, 100000);
        return d2Vector(params[0], params[1]);
    }

}

d3Vector DistributionHelper::fitCorruptedBivariateStudent(
        const std::vector<std::vector<d2Vector> > &coords, int f, double area, double fixedNu)
{
    int pc = coords.size();

    std::vector<double> velocities(pc);

    for (int p = 0; p < pc; p++)
    {
        velocities[p] = (coords[p][f+1] - coords[p][f]).length();
    }

    CorruptedBivariateStudentFit csf(velocities, area, fixedNu);

    if (fixedNu > 0.0)
    {
        std::vector<double> initial = {1.0, 0.01};

        std::vector<double> params = NelderMead::optimize(initial, csf, 0.1, 0.001, 10000);

        return d3Vector(params[0], fixedNu, params[1]);
    }
    else
    {
        std::vector<double> initial = {1.0, 1.0, 0.01};

        std::vector<double> params = NelderMead::optimize(initial, csf, 0.1, 0.001, 10000);

        return d3Vector(params[0], params[1], params[2]);
    }

}

std::vector<double> DistributionHelper::fitDoubleStudent(
        const std::vector<std::vector<d2Vector> >& coords,
        const std::vector<int>& counts,
        int f, double area,
        bool fixedNu0, double nu0,
        bool fixedNuS, double nuS,
        bool relTheta)
{
    int pc = coords.size();

    std::vector<double> velocities(pc);

    for (int p = 0; p < pc; p++)
    {
        velocities[p] = (coords[p][f+1] - coords[p][f]).length();
    }

    DoubleStudentFit dsf(velocities, counts, nu0, fixedNu0, nuS, fixedNuS, area, relTheta);

    if (fixedNu0 && fixedNuS)
    {
        std::vector<double> initial = {1.0, 0.0, 0.01};
        std::vector<double> params = NelderMead::optimize(initial, dsf, 0.1, 0.001, 10000);
        return std::vector<double>{params[0], params[1], params[2], nu0, nuS};
    }
    else if (fixedNu0 && !fixedNuS) // pathological
    {
        std::vector<double> initial = {1.0, 0.0, 0.01, 0.0};
        std::vector<double> params = NelderMead::optimize(initial, dsf, 0.1, 0.001, 10000);
        return std::vector<double>{params[0], params[1], params[2], nu0, params[3]};
    }
    else if (!fixedNu0 && fixedNuS)
    {
        std::vector<double> initial = {1.0, 0.0, 0.01, 2.0};
        std::vector<double> params = NelderMead::optimize(initial, dsf, 0.1, 0.001, 10000);
        return std::vector<double>{params[0], params[1], params[2], params[3], nuS};
    }
    else // (!fixedNu0 && !fixedNuS)
    {
        std::vector<double> initial = {1.0, 0.0, 0.01, 2.0, 0.0};
        std::vector<double> params = NelderMead::optimize(initial, dsf, 0.1, 0.001, 10000);
        return std::vector<double>{params[0], params[1], params[2], params[3], params[4]};
    }

}

std::vector<double> DistributionHelper::fitDoubleGauss(
        const std::vector<std::vector<d2Vector> > &coords,
        const std::vector<int> &counts,
        const std::vector<int> &totalCounts,
        int f, double area, bool centered, std::vector<double> initial)
{
    int pc = coords.size();

    std::vector<double> velocities(pc);

    for (int p = 0; p < pc; p++)
    {
        velocities[p] = (coords[p][f+1] - coords[p][f]).length();
    }

    DoubleGaussFit dgf(velocities, counts, totalCounts, area, centered);

    std::vector<double> params = NelderMead::optimize(initial, dgf, 0.02, 0.00001, 100000);
    return std::vector<double>{params[0], params[1], params[2], params[3], params[4]};
}

std::vector<double> DistributionHelper::fitCorruptedGauss(
        const std::vector<std::vector<d2Vector> > &coords, int f, double area)
{
    int pc = coords.size();

    std::vector<double> velocities(pc);

    for (int p = 0; p < pc; p++)
    {
        velocities[p] = (coords[p][f+1] - coords[p][f]).length();
    }

    CorruptedGaussFit cgf(velocities, area);

    std::vector<double> initial = {1.0, 0.01};
    std::vector<double> params = NelderMead::optimize(initial, cgf, 0.1, 0.001, 10000);

    return params;
}

std::vector<double> DistributionHelper::fitMultiStudent(
        const std::vector<std::vector<d2Vector> > &coords,
        const std::vector<int> &counts,
        int f, double area, int minPC, int maxPC)
{
    int pc = coords.size();

    std::vector<double> velocities(pc);

    for (int p = 0; p < pc; p++)
    {
        velocities[p] = (coords[p][f+1] - coords[p][f]).length();
    }

    MultiStudentFit msf(velocities, counts, minPC, maxPC, area);
    const int pbc = maxPC - minPC + 1;

    std::vector<double> initial(2 + 2*pbc);

    initial[0] = 2.0;
    initial[1] = 0.0;

    for (int i = 0; i < pbc; i++)
    {
        initial[2 + 2*i] = 3.0; // nu
        initial[2 + 2*i + 1] = 0.2; // theta
    }

    std::vector<double> params = NelderMead::optimize(initial, msf, 0.1, 0.001, 100000);

    return params;
}

d4Vector DistributionHelper::fitGaussStudent(
        const std::vector<std::vector<d2Vector> > &coords, int f, double fixedNu)
{
    int pc = coords.size();

    std::vector<double> velocities(pc);

    for (int p = 0; p < pc; p++)
    {
        velocities[p] = (coords[p][f+1] - coords[p][f]).length();
    }

    GaussStudentFit gsf(velocities, fixedNu);

    if (fixedNu > 0.0)
    {
        std::vector<double> initial = {0.5, 2.0, 0.1};

        std::vector<double> params = NelderMead::optimize(initial, gsf, 0.1, 0.001, 10000);

        return d4Vector(params[0], params[1], params[2], fixedNu);
    }
    else
    {
        std::vector<double> initial = {1.0, 1.0, 0.5, 2.0};

        std::vector<double> params = NelderMead::optimize(initial, gsf, 0.1, 0.001, 10000);

        return d4Vector(params[0], params[1], params[2], params[3]);
    }

}


CorruptedExponentialFit::CorruptedExponentialFit(const std::vector<double> &velocities, double area)
    :   velocities(velocities),
        area(area)
{}

double CorruptedExponentialFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int pc = velocities.size();
    double e = 0.0;

    const double beta = x[0];
    const double theta = x[1];
    const double twoPiBetaSq = 2.0*PI*beta*beta;

    if (beta <= 0.0 || theta <= 0.0 || theta > 1.0) return std::numeric_limits<double>::max();

    for (int p = 0; p < pc; p++)
    {
        e += -log(theta/area + (1.0 - theta)*exp(-velocities[p]/beta)/twoPiBetaSq);
    }

    return e;
}


BivariateStudentFit::BivariateStudentFit(const std::vector<double> &values, double fixedNu)
:   values(values),
    fixedNu(fixedNu)
{
}

double BivariateStudentFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int pc = values.size();
    double e = 0.0;

    const double sigma = x[0];
    const double nu = fixedNu > 0.0? fixedNu : x[1];
    const double logTwoPiSigmaSq = log(2.0*PI*sigma*sigma*nu);
    const double nuSigmaSq = nu*sigma*sigma;

    if (sigma <= 0.0 || nu <= 0.0) return std::numeric_limits<double>::max();

    for (int p = 0; p < pc; p++)
    {
        const double v = values[p];
        e += logTwoPiSigmaSq + (1.0 + nu/2.0) * log(1.0 + v*v/nuSigmaSq);
    }

    return e;
}

CorruptedBivariateStudentFit::CorruptedBivariateStudentFit(const std::vector<double> &values, double area, double fixedNu)
    :   values(values),
        area(area),
        fixedNu(fixedNu)

{
}

double CorruptedBivariateStudentFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int pc = values.size();
    double e = 0.0;

    const double sigma = x[0];
    const double nu = fixedNu > 0.0? fixedNu : x[1];
    const double theta = fixedNu > 0.0? x[1] : x[2];
    const double rho = 1.0/area;

    const double twoPiSigmaSq = 2.0*PI*sigma*sigma;
    const double nuSigmaSq = nu*sigma*sigma;

    if (sigma <= 0.0 || nu <= 0.0 || theta <= 0.0 || theta >= 1.0) return std::numeric_limits<double>::max();

    for (int p = 0; p < pc; p++)
    {
        const double y = values[p];
        e += -log(theta * rho + (1.0 - theta) * (1.0/twoPiSigmaSq) * pow(1.0 + y*y/nuSigmaSq, -(1.0 + nu/2.0)));
    }

    return e;

}

GaussStudentFit::GaussStudentFit(const std::vector<double> &values, double fixedNu)
:   values(values),
    fixedNu(fixedNu)
{
}

double GaussStudentFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int pc = values.size();
    double e = 0.0;

    const double sigmaS = x[0];
    const double sigmaG = x[1];
    const double theta = x[2];
    const double nu = fixedNu > 0.0? fixedNu : x[3];

    const double twoPiSigmaSSq = 2.0*PI*sigmaS*sigmaS;
    const double nuSigmaSSq = nu*sigmaS*sigmaS;
    const double twoPiSigmaGSq = 2.0*PI*sigmaG*sigmaG;

    if (sigmaG <= 0.0 || sigmaS <= 0.0 || nu <= 0.0 || theta <= 0.0 || theta >= 1.0)
    {
        return std::numeric_limits<double>::max();
    }

    for (int p = 0; p < pc; p++)
    {
        const double y = values[p];
        e += -log(theta * (1.0/twoPiSigmaGSq) * exp(-0.5*y*y/(sigmaG*sigmaG))
             + (1.0 - theta) * (1.0/twoPiSigmaSSq) * pow(1.0 + y*y/nuSigmaSSq, -(1.0 + nu/2.0)));
    }

    return e;
}

DoubleStudentFit::DoubleStudentFit(const std::vector<double> &values,
        const std::vector<int> &counts,
        double nu0, bool fixedNu0,
        double nuS, bool fixedNuS,
        double area, bool relTheta)
:   values(values),
    counts(counts),
    myNu0(nu0),
    myNuS(nuS),
    area(area),
    fixedNu0(fixedNu0),
    fixedNuS(fixedNuS),
    relTheta(relTheta)
{
}

double DoubleStudentFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int pc = values.size();
    double e = 0.0;

    const double sigma0 = x[0];
    const double sigmaS = x[1];
    const double thetaS = x[2];

    double nu0, nuS;

    if (fixedNu0 && fixedNuS)
    {
        nu0 = myNu0;
        nuS = myNuS;
    }
    else if (!fixedNu0 && fixedNuS)
    {
        nu0 = x[3];
        nuS = myNuS;
    }
    else if (fixedNu0 && !fixedNuS) // pathological
    {
        nu0 = myNu0;
        nuS = x[3];
    }
    else if (!fixedNu0 && !fixedNuS)
    {
        nu0 = x[3];
        nuS = x[4];
    }

    const double rho = 1.0 / area;

    if (sigma0 <= 0.0 || nu0 <= 0.0 || thetaS <= 0.0)
    {
        return std::numeric_limits<double>::max();
    }

    for (int p = 0; p < pc; p++)
    {
        if (counts[p] <= 0) continue;

        const double Ninv = 1.0 / (double)counts[p];
        const double sig2 = sigma0 + Ninv * sigmaS;
        const double nu = nu0 + Ninv * nuS;

        if (sig2 <= 0.0 || nu <= 0.0)
        {
            return std::numeric_limits<double>::max();
        }

        double theta = relTheta? thetaS*Ninv : thetaS;

        if (theta < 0.0) theta = 0.0;
        if (theta > 1.0) theta = 1.0;

        const double twoPiSigmaSq = 2.0*PI*sig2;
        const double nuSigmaSq = nu*sig2;

        const double y = values[p];

        e += -log(theta * rho + (1.0 - theta) * (1.0/twoPiSigmaSq) * pow(1.0 + y*y/nuSigmaSq, -(1.0 + nu/2.0)));
    }

    return e;
}

MultiStudentFit::MultiStudentFit(const std::vector<double> &values,
        const std::vector<int> &counts, int minPN,
        int maxPN, double area)
:   values(values),
    counts(counts),
    minPN(minPN),
    maxPN(maxPN),
    area(area)
{
}

double MultiStudentFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int pc = values.size();
    double e = 0.0;

    const double sigma0 = x[0];
    const double sigmaS = x[1];

    if (sigma0 <= 0.0)
    {
        return std::numeric_limits<double>::max();
    }

    const double rho = 1.0 / area;
    const int pbn = maxPN - minPN + 1;

    std::vector<double> nu(pbn);
    std::vector<double> theta(pbn);

    for (int i = 0; i < pbn; i++)
    {
        nu[i] = x[2 + 2*i];
        theta[i] = x[2 + 2*i + 1];

        if (nu[i] <= 0.0 || theta[i] < 0.0 || theta[i] > 1.0)
        {
            return std::numeric_limits<double>::max();
        }
    }

    for (int p = 0; p < pc; p++)
    {
        const int pb = counts[p] - minPN;

        if (pb < 0 || pb >= pbn) continue;

        const double Ninv = 1.0 / (double)counts[p];
        const double sig2 = sigma0 + Ninv * sigmaS;

        if (sig2 <= 0.0)
        {
            return std::numeric_limits<double>::max();
        }

        const double thetaP = theta[pb];
        const double nuP = nu[pb];

        const double twoPiSigmaSq = 2.0*PI*sig2;
        const double nuSigmaSq = nuP*sig2;

        const double y = values[p];

        e += -log(thetaP * rho + (1.0 - thetaP) * (1.0/twoPiSigmaSq) * pow(1.0 + y*y/nuSigmaSq, -(1.0 + nuP/2.0)));
    }

    return e;
}

DoubleGaussFit::DoubleGaussFit(const std::vector<double> &values,
        const std::vector<int> &counts, const std::vector<int> &totalCounts,
        double area,
        bool centered)
:   values(values),
    counts(counts),
    totalCounts(totalCounts),
    area(area),
    centered(centered)
{
}

double DoubleGaussFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int pc = values.size();
    double e = 0.0;

    const double sigma0 = x[0];
    const double sigmaS = x[1];
    const double theta0 = x[2];
    const double thetaS = x[3];
    const double thetaMin = x[4];

    const double rho = 1.0 / area;

    if (sigma0 <= 0.0 || sigmaS <= 0.0 || thetaS < 0.0 || thetaMin < 0.0)
    {
        return std::numeric_limits<double>::max();
    }

    for (int p = 0; p < pc; p++)
    {
        if (counts[p] <= 0) continue;

        const double Ninv = centered? 1.0 / (double)counts[p] + 1.0 / (double)totalCounts[p] : 1.0 / (double)counts[p];
        const double correction = centered? (1.0 - (double)counts[p] / (double)totalCounts[p]) : 1.0;

        if (correction == 0) continue;

        const double sig2 = correction * (sigma0 + Ninv * sigmaS);

        if (sig2 <= 0.0)
        {
            return std::numeric_limits<double>::max();
        }

        double theta = theta0 + Ninv * thetaS;

        if (theta < thetaMin) theta = thetaMin;
        if (theta > 1.0) theta = 1.0;

        const double twoPiSigmaSq = 2.0*PI*sig2;

        const double y = values[p];

        e += -log(theta * rho + (1.0 - theta) * (1.0/twoPiSigmaSq) * exp(-0.5*y*y/sig2));
    }

    return e;
}

CorruptedGaussFit::CorruptedGaussFit(const std::vector<double> &values, double area)
    : values(values),
      area(area)
{
}

double CorruptedGaussFit::f(const std::vector<double> &x, void* tempStorage) const
{
    const int pc = values.size();
    double e = 0.0;

    const double sigma = x[0];
    const double theta = x[1];

    const double rho = 1.0 / area;

    if (sigma <= 0.0 || theta < 0.0 || theta > 1.0)
    {
        return std::numeric_limits<double>::max();
    }

    const double sig2 = sigma * sigma;
    const double twoPiSigmaSq = 2.0*PI*sig2;

    for (int p = 0; p < pc; p++)
    {
        const double y = values[p];

        e += -log(theta * rho + (1.0 - theta) * (1.0/twoPiSigmaSq) * exp(-0.5*y*y/sig2));
    }

    return e;
}

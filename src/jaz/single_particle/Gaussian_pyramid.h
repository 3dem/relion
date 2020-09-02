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

#ifndef GAUSSIAN_PYRAMID_H
#define GAUSSIAN_PYRAMID_H

#include <src/image.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/resampling_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <vector>
#include <omp.h>

template <typename T>
class GaussianPyramid
{
    public:

        GaussianPyramid(const Image<T>& img, double sigma = 1.0,
                        double cutoff = 2.5, int maxLev = -1);

            double sigma, cutoff;
            std::vector<Image<T>> levels;
            std::vector<double> scales;

        T value(double x, double y, double s);

        std::vector<Image<T>> getUpsampled();
        Image<T> getUpsampledLevel(int i);
        Image<T> getInterpolated(double sig);

        static void test(const Image<T>& img);
        static void interpolationTest(const Image<T>& img);
        static void timeTest(const Image<T>& img0);

};

template <typename T>
GaussianPyramid<T>::GaussianPyramid(const Image<T>& img, double sigma, double cutoff, int maxLev)
:   sigma(sigma), cutoff(cutoff), levels(0), scales(0)
{
    int w0 = img.data.xdim;
    int h0 = img.data.ydim;

    int w = w0;
    int h = h0;

    Image<T> curr = img;
    double currScale = 1.0;
    double currL0Var = 0.0;
    double wLast = w;

    if (PI * sigma < cutoff)
    {
        REPORT_ERROR("GaussianPyramid::GaussianPyramid: unable to shrink image levels - cutoff too large.\n");
    }

    levels.push_back(curr);
    scales.push_back(currScale);

    double dScale = PI * sigma / cutoff;
    int i = 0;

    if (maxLev < 0) maxLev = floor(log((double)w0)/log(dScale));

    while (i < maxLev)
    {
        double tgtSigma_l0 = sigma * pow(dScale, i);

        double dSigma = sqrt(tgtSigma_l0 * tgtSigma_l0 - currL0Var) / currScale;
        currL0Var += currScale * currScale * dSigma * dSigma;

        int wNext = (int)(w / currScale + 0.5);
        double dScaleAct = wLast / wNext;
        wLast = wNext;

        currScale *= dScaleAct;

        Image<T> filt, next;
        FilterHelper::separableGaussianXY(curr, filt, dSigma, 3*dSigma, true);
        ResamplingHelper::subsample2D_cubic(filt, dScaleAct, next);

        levels.push_back(next);
        scales.push_back(currScale);
        curr = next;

        i++;
    }
}

template<typename T>
std::vector<Image<T>> GaussianPyramid<T>::getUpsampled()
{
    std::vector<Image<T>> out(levels.size());

    out[0] = levels[0];

    for (int i = 1; i < levels.size(); i++)
    {
        ResamplingHelper::upsample2D_cubic(levels[i], scales[i], out[i], true);
    }

    return out;
}

template<typename T>
Image<T> GaussianPyramid<T>::getUpsampledLevel(int i)
{
    Image<T> out;
    ResamplingHelper::upsample2D_cubic(levels[i], scales[i], out, true);
    return out;
}

template<typename T>
Image<T> GaussianPyramid<T>::getInterpolated(double sig)
{
    Image<T> out(levels[0].data.xdim, levels[0].data.ydim);

    if (sig <= 0) return levels[0];
    else if (sig < sigma)
    {
        int k = 3*sig < 1? 1 : 3*sig;
        FilterHelper::separableGaussianXY(levels[0], out, sig, k, true);
    }
    else
    {
        double lambda = PI * sigma / cutoff;
        double l = log(sig/sigma) / log(lambda) + 1.0;

        int l0 = (int)l;
        double f = l - l0;

        if (l0 >= levels.size()-1) return levels[levels.size()-1];


        Image<T> img0 = getUpsampledLevel(l0);
        Image<T> img1 = getUpsampledLevel(l0+1);

        ImageOp::linearCombination(img0, img1, 1.0 - f, f, out);
    }

    return out;
}

template <typename T>
void GaussianPyramid<T>::test(const Image<T>& img)
{
    int lev = 10;
    int mp = 5;
    double k = 2.5;
    double sg = 1.0;

    GaussianPyramid<T> gp(img, 1.0, 2.5, lev);

    std::vector<Image<T>> test(lev*mp), baseline(lev*mp);

    double lambda = PI*sg/k;

    for (int i = 0; i < lev; i++)
    for (int j = 0; j < mp; j++)
    {
        double t = (mp*i + j)/(double)(mp);
        double sig = t < 1.0? sg * t : sg * pow(lambda,t-1);

        test[mp*i + j] = gp.getInterpolated(sig);
        FilterHelper::separableGaussianXY(img, baseline[mp*i + j], sig, 3*sig, true);
    }

    ImageLog::write(test, "debug/pyr_interpolated");
    ImageLog::write(baseline, "debug/pyr_baseline");
}


template<typename T>
void GaussianPyramid<T>::interpolationTest(const Image<T> &img)
{
    const double ds = 0.25;

    for (double s0 = 0.0; s0 < 2.0; s0 += ds)
    {
        double s1 = s0 + ds;

        Image<T> img0, img1, imgP, imgT;

        FilterHelper::separableGaussianXY(img, img0, s0, 6, true);
        FilterHelper::separableGaussianXY(img, img1, s1, 6, true);

        imgT.data.resize(img.data);

        std::stringstream stss0;
        stss0 << s0;
        std::stringstream stss1;
        stss1 << s1;

        std::ofstream ofs("debug/interp_"+stss0.str()+"-"+stss1.str()+".dat");

        for (double p = s0; p <= s1; p += ds/50.0)
        {
            FilterHelper::separableGaussianXY(img, imgP, p, 6, true);

            double c_min = 1e20;
            double t_opt = 0.0;

            for (double t = 0.0; t < 1.0; t += 0.005)
            {
                ImageOp::linearCombination(img0, img1, 1.0 - t, t, imgT);

                double c = FilterHelper::L2distance(imgT, imgP);

                if (c < c_min)
                {
                    c_min = c;
                    t_opt = t;
                }
            }

            ofs << p << " " << s0 + t_opt*(s1 - s0) << "\n";
            std::cout << p << " between " << s0 << " and " << s1 << ": t = " << t_opt << "\n";
        }

        std::cout << "\n";
    }
}

template <typename T>
void GaussianPyramid<T>::timeTest(const Image<T>& img0)
{
    int lev = 10;
    int pc = 20;
    double sig = 5.0;

    const int s = img0.data.xdim;
    const int sh = img0.data.xdim/2 + 1;

    Image<RFLOAT> img = img0;
    Image<RFLOAT> sum(s,s);
    sum.data.initZeros();

    double t0 = omp_get_wtime();

    GaussianPyramid<T> gp(img, 1.0, 2.5, lev);

    for (int p = 0; p < pc-1; p++)
    {
        gp = GaussianPyramid<T>(img, 1.0, 2.5, lev);
    }

    double t1 = omp_get_wtime();

    std::cout << "creating " << pc << " pyramids: " << (t1 - t0) << "s\n";

    Image<T> bl;

    for (int p = 0; p < pc; p++)
    for (int q = 0; q < pc-1; q++)
    {
        bl = gp.getInterpolated(5.0);
        ImageOp::linearCombination(sum, bl, 1, 1, sum);
    }

    ImageLog::write(sum, "debug/sum_pyr");

    double t2 = omp_get_wtime();

    std::cout << "reading " << pc*(pc-1) << " images: " << (t2 - t1) << "s\n";
    std::cout << "both: " << (t2 - t0) << "s\n";


    FourierTransformer ft;
    Image<Complex> fs(sh,s);
    Image<RFLOAT> rs(s,s);

    double sig_hat = s/(2.0 * PI * sig);
    double sh2 = sig_hat * sig_hat;
    sum.data.initZeros();

    for (int p = 0; p < pc; p++)
    for (int q = 0; q < pc-1; q++)
    {
        ft.FourierTransform(img(), fs());

        for (int y = 0; y < s; y++)
        for (int x = 0; x < sh; x++)
        {
            double yy = y < sh? y : y - s;
            double xx = x;

            fs(y,x) *= exp(-0.5*(xx*xx + yy*yy)/sh2);
        }

        ft.inverseFourierTransform(fs(), rs());

        ImageOp::linearCombination(sum, rs, 1, 1, sum);
    }

    ImageLog::write(sum, "debug/sum_fftw");

    double t3 = omp_get_wtime();

    std::cout << "performing " << pc*(pc-1) << " Fourier convolutions: " << (t3 - t2) << "s\n";

}

#endif

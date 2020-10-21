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

#ifndef FILTER_HELPER_H
#define FILTER_HELPER_H

#include <src/image.h>
#include <src/fftw.h>
#include <src/strings.h>
#include <src/ctf.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/gravis/t3Matrix.h>
#include <src/jaz/single_particle/volume.h>

class FilterHelper
{
    public:

        template<typename T>
        static void binomial3x3_2D(const Image<T>& src, Image<T>& dest, bool wrap = false);

        template<typename T>
        static void separableGaussian(const Volume<T>& src, Volume<T>& dest, double sigma, int k = -1);

        template<typename T>
        static void separableGaussian(const MultidimArray<T>& src, MultidimArray<T>& dest, double sigma, int k = -1);

        template<typename T>
        static void separableGaussianWrap(const MultidimArray<T>& src, MultidimArray<T>& dest, double sigma, int k = -1);

        static void separableGaussianFreq(
                const MultidimArray<Complex>& src, MultidimArray<Complex>& dest,
                double sigma, int k = -1);

        static void separableGaussianFreqXY(
                const MultidimArray<Complex>& src, MultidimArray<Complex>& dest,
                double sigma, int k = -1);

        static void drawTestPattern(Image<RFLOAT>& img, int squareSize);
        static void drawTestPattern(Volume<RFLOAT>& volume, int squareSize);

        static Image<RFLOAT> expImg(Image<RFLOAT>& img, double scale = 1.0);
        static Image<RFLOAT> logImg(Image<RFLOAT>& img, double thresh = 1e-20, double scale = 1.0);
        static Image<RFLOAT> padCorner2D(Image<RFLOAT>& img, double factor);
        static Image<Complex> padCorner2D(Image<Complex>& img, double factor);
        static Image<RFLOAT> padCorner2D(const Image<RFLOAT> &img, int w, int h);
        static Image<RFLOAT> cropCorner2D(const Image<RFLOAT> &img, int w, int h);
        static Image<Complex> cropCorner2D(const Image<Complex> &img, int w, int h);
        static Image<RFLOAT> zeroOutsideCorner2D(Image<RFLOAT>& img, double radius);
        static void GaussianEnvelopeCorner2D(Image<RFLOAT>& img, double sigma);
		static Image<RFLOAT> raisedCosEnvCorner2D(Image<RFLOAT>& img, double radIn, double radOut);
		static Image<Complex> raisedCosEnvCorner2DFull(Image<Complex>& img, double radIn, double radOut);
		static Image<RFLOAT> raisedCosEnvCorner3D(Image<RFLOAT>& img, double radIn, double radOut);
		static Image<RFLOAT> raisedCosEnvFreq2D(const Image<RFLOAT>& img, double radIn, double radOut);
		static Image<RFLOAT> raisedCosEnvRingFreq2D(const Image<RFLOAT>& img, double rad0, double rad1, double stepWidth);


        static void lowPassFilter(Image<RFLOAT>& img, double maxFreq0, double maxFreq1, Image<RFLOAT>& dest);
        static void lowPassFilterSpectrum(MultidimArray<Complex>& spectrum, double maxFreq0, double maxFreq1);

        static RFLOAT averageValue(Image<RFLOAT>& img);
        static RFLOAT maxValue(Image<RFLOAT>& img);

        static void phaseFlip(Image<RFLOAT>& img, CTF& ctf, RFLOAT angpix, Image<RFLOAT>& dest);
        static void applyBeamTilt(Image<RFLOAT>& img, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
                                  RFLOAT lambda, RFLOAT Cs, RFLOAT angpix, int s, Image<RFLOAT>& dest);
        static void modulate(Image<RFLOAT>& img, CTF& ctf, RFLOAT angpix, Image<RFLOAT>& dest);
        static void modulate(Image<Complex>& imgFreq, CTF& ctf, RFLOAT angpix, Image<RFLOAT>& dest);
        static void modulate(MultidimArray<Complex>& imgFreq, CTF& ctf, RFLOAT angpix);
        static void drawCtf(CTF& ctf, RFLOAT angpix, Image<Complex>& dest);
        static void wienerFilter(Image<RFLOAT>& img, CTF& ctf, RFLOAT angpix, RFLOAT eps, RFLOAT Bfac, Image<RFLOAT>& dest);
        static void richardsonLucy(Image<RFLOAT>& img, CTF& ctf, RFLOAT angpix, RFLOAT eps, int iterations, Image<RFLOAT>& dest);
        static void rampFilter(Image<RFLOAT>& img, RFLOAT s0, RFLOAT t1, double ux, double uy, Image<RFLOAT>& dest);
        static void rampFilter3D(Image<Complex>& img, RFLOAT s0, RFLOAT t1, double tx, double ty, double tz);
        static void doubleRampFilter3D(Image<Complex>& img, RFLOAT s0, RFLOAT t1, double tx, double ty, double tz);
        static void getPhase(const Image<Complex>& img, Image<RFLOAT>& dest);
        static void getAbs(const Image<Complex>& img, Image<RFLOAT>& dest);
        static void getReal(const Image<Complex>& img, Image<RFLOAT>& dest);
        static void getImag(const Image<Complex>& img, Image<RFLOAT>& dest);

        static void powerSpectrum2D(Image<RFLOAT>& img, Volume<RFLOAT>& spectrum);
        static void equiphaseAverage2D(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest);

        static void threshold(Image<RFLOAT>& src, RFLOAT t, Image<RFLOAT>& dest);
        static void fill(Image<RFLOAT>& dest, RFLOAT v);
        static void linearTransform(Image<RFLOAT>& src, RFLOAT m, RFLOAT q, Image<RFLOAT>& dest);
        static void linearCombination(Image<RFLOAT>& src0, Image<RFLOAT>& src1, RFLOAT a0, RFLOAT a1, Image<RFLOAT>& dest);
        static void linearCombination(const Volume<RFLOAT>& src0, const Volume<RFLOAT>& src1, RFLOAT a0, RFLOAT a1, Volume<RFLOAT>& dest);
        static void sumUp(const std::vector<Image<RFLOAT> >& src, Image<RFLOAT>& dest);

        static double L1distance(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0 = 0, int y0 = 0, int w = -1, int h = -1);
        static double L2distance(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0 = 0, int y0 = 0, int w = -1, int h = -1);
        static double NCC(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0 = 0, int y0 = 0, int w = -1, int h = -1);

        static void multiply(Image<RFLOAT>& i0, Image<RFLOAT>& i1, Image<RFLOAT>& dest);
        static void multiply(Image<Complex>& i0, Image<Complex>& i1, Image<Complex>& dest);
        static void wienerDivide(Image<RFLOAT>& num, Image<RFLOAT>& denom, RFLOAT eps, Image<RFLOAT>& dest);
        static void divide(Image<Complex>& num, Volume<RFLOAT>& denom, RFLOAT eps, Image<Complex>& dest);
        static void divide(Image<Complex>& num, Image<RFLOAT>& denom, RFLOAT eps, Image<Complex>& dest);
        static void divideExcessive(Image<Complex>& num, Volume<RFLOAT>& denom, RFLOAT theta, Image<Complex>& dest);
        static void wienerDeconvolve(Image<Complex>& num, Image<Complex>& denom, RFLOAT theta, Image<Complex>& dest);

        static void extract2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                            long int x0, long int y0,
                            long int w, long int h);

        static void extract(const Volume<RFLOAT>& src,
                            Volume<RFLOAT>& dest,
                            long int x0, long int y0, long int z0,
                            long int w, long int h, long int d);

        static void signedDist(const Image<RFLOAT>& src, Image<RFLOAT>& dest);
        static void erode3x3(Image<RFLOAT>& src, Image<RFLOAT>& dest);
        static void localMinima(Image<RFLOAT>& src, Image<RFLOAT>& dest, RFLOAT thresh);
        static std::vector<gravis::d3Vector> localMinima(Image<RFLOAT>& src, RFLOAT thresh);

        static void uniqueInfluenceMask(std::vector<gravis::d3Vector> pts, Image<RFLOAT>& dest, Image<RFLOAT>& indexDest, RFLOAT thresh);
        static void polarRemap(gravis::d2Vector pos, const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                               const Image<RFLOAT>& mask, Image<RFLOAT>& maskDest, int phiRes, int rRes, double rMax);
        static void polarRemap(gravis::d2Vector pos, const Image<RFLOAT>& distTransf, const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                               const Image<RFLOAT>& mask, Image<RFLOAT>& maskDest, int phiRes, int rRes, double rMax);

        static Image<RFLOAT> cartToPolar(const Image<RFLOAT>& img);
        static Image<RFLOAT> polarToCart(const Image<RFLOAT>& img);
        static Image<RFLOAT> polarBlur(const Image<RFLOAT>& img, double sigma);
        static Image<RFLOAT> sectorBlend(const Image<RFLOAT>& img0, const Image<RFLOAT>& img1, int sectors);


        static void diffuseAlongIsocontours2D(const Image<RFLOAT>& src, const Image<RFLOAT>& guide,
                                              Image<RFLOAT>& dest, int iters, RFLOAT sigma, RFLOAT lambda, RFLOAT delta);

        static void EED_2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest, int iters, double sigma, double delta, double tau);

        static void descendTV(const Image<RFLOAT>& src, Image<RFLOAT>& dest, double delta);
        static void descendTV2(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                               Volume<gravis::d3Vector>& xi, Volume<RFLOAT>& uBar,
                               int iters, double sigma, double tau);

        static void segmentTV(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                               Volume<gravis::d3Vector>& xi, Volume<RFLOAT>& uBar,
                               int iters, double sigma, double tau, double nu);

        static void segmentTVAniso2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                               Volume<gravis::d2Vector>& xi, Volume<RFLOAT>& uBar,
                               int iters, double sigma, double tau, double nu,
                               double rho, double theta, double alpha);

        static void fwdGrad(const Image<RFLOAT>& u, Volume<gravis::d3Vector>& dest);
        static void fwdGrad2D(const Image<RFLOAT>& u, Volume<gravis::d2Vector>& dest);
        static void centralGrad2D(const Image<RFLOAT>& u, Volume<gravis::d2Vector>& dest);
        static void centralGrad2D(const Image<Complex>& u, Volume<gravis::d2Vector>& destRe, Volume<gravis::d2Vector>& destIm);

        static void blendSoft(const Image<Complex>& src0, const Image<Complex>& src1,
                              const Volume<RFLOAT>& mask, Image<Complex>& dest, RFLOAT bias1 = 1.0);

        static double totalVariation(const Image<RFLOAT>& src);
        static double totalLogVariation(const Image<RFLOAT>& src, double delta = 1.0);

        static void separableGaussianXYZ(const Image<RFLOAT>& src, Image<RFLOAT>& dest, RFLOAT sigma, int k = -1);
        static void separableGaussianXY(const Image<RFLOAT>& src, Image<RFLOAT>& dest, RFLOAT sigma, int k = -1, bool wrap = false);
        static void separableGaussianX_wrap(const Image<RFLOAT>& src, const Image<RFLOAT>& mask, Image<RFLOAT>& dest, RFLOAT sigma, int k = -1);
        static void separableGaussianX_wrap(const Image<RFLOAT>& src, Image<RFLOAT>& dest, RFLOAT sigma, int k = -1);
        static void averageX(const Image<RFLOAT>& src, const Image<RFLOAT>& mask, Image<RFLOAT>& dest);

        static void centralGradient(const Volume<RFLOAT>& src, Volume<gravis::t3Vector<RFLOAT> >& dest);
        static gravis::t3Vector<RFLOAT> centralGradient(const Volume<RFLOAT>& src, size_t x, size_t y, size_t z);
		
		static MultidimArray<Complex> FriedelExpand(const MultidimArray<Complex>& half);
		
		static Image<RFLOAT> normaliseToUnitInterval(const Image<RFLOAT>& img);
		static Image<RFLOAT> normaliseToUnitIntervalSigned(const Image<RFLOAT>& img);


};


template<typename T>
void FilterHelper::binomial3x3_2D(const Image<T>& src, Image<T>& dest, bool wrap)
{
    const size_t w = src.data.xdim;
    const size_t h = src.data.ydim;
    const size_t d = src.data.zdim;

    dest.data.reshape(d,h,w);

    Image<T> temp(w,h);
    std::vector<double> kernel = {0.25, 0.5, 0.25};

    for (size_t z = 0; z < d; z++)
    for (size_t y = 0; y < h; y++)
    for (size_t x = 0; x < w; x++)
    {
        T v = 0;
        double m = 0;

        for (int i = -1; i <= 1; i++)
        {
            int xx = x + i;
            if (wrap) xx = (xx + w) % w;
            else if (xx < 0 || xx >= w) continue;

            v += kernel[i+1] * DIRECT_NZYX_ELEM(src(), 0, z, y, xx);
            m += kernel[i+1];
        }

        DIRECT_NZYX_ELEM(temp(), 0, z, y, x) = v/m;
    }

    for (size_t z = 0; z < d; z++)
    for (size_t y = 0; y < h; y++)
    for (size_t x = 0; x < w; x++)
    {
        T v = 0;
        double m = 0;

        for (int i = -1; i <= 1; i++)
        {
            int yy = y + i;
            if (wrap) yy = (yy + h) % h;
            else if (yy < 0 || yy >= h) continue;

            v += kernel[i+1] * DIRECT_NZYX_ELEM(temp(), 0, z, yy, x);
            m += kernel[i+1];
        }

        DIRECT_NZYX_ELEM(dest(), 0, z, y, x) = v/m;
    }
}

template<typename T>
void FilterHelper::separableGaussian(const Volume<T>& src, Volume<T>& dest, double sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.resize(src);

    std::vector<double> kernel(2*k+1);
    const double s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-i*i/s2);
    }

    Volume<T> temp(src.dimx, src.dimy, src.dimz);

    for (size_t z = 0; z < src.dimz; z++)
    for (size_t y = 0; y < src.dimy; y++)
    for (size_t x = 0; x < src.dimx; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int xx = x + i;
            if (xx < 0 || xx >= src.dimx) continue;

            v += kernel[i+k] * src(xx,y,z);
            m += kernel[i+k];
        }

        dest(x,y,z) = v/m;
    }

    for (size_t z = 0; z < src.dimz; z++)
    for (size_t y = 0; y < src.dimy; y++)
    for (size_t x = 0; x < src.dimx; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = y + i;
            if (yy < 0 || yy >= src.dimy) continue;

            v += kernel[i+k] * dest(x,yy,z);
            m += kernel[i+k];
        }

        temp(x,y,z) = v/m;
    }

    for (size_t z = 0; z < src.dimz; z++)
    for (size_t y = 0; y < src.dimy; y++)
    for (size_t x = 0; x < src.dimx; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int zz = z + i;
            if (zz < 0 || zz >= src.dimz) continue;

            v += kernel[i+k] * temp(x,y,zz);
            m += kernel[i+k];
        }

        dest(x,y,z) = v/m;
    }
}

template<typename T>
void FilterHelper::separableGaussian(const MultidimArray<T>& src, MultidimArray<T>& dest, double sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.reshape(src);

    std::vector<double> kernel(2*k+1);
    const double s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-0.5*i*i/s2);
    }

    MultidimArray<T> temp(src.zdim, src.ydim, src.xdim);

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int xx = x + i;
            if (xx < 0 || xx >= src.xdim) continue;

            v += kernel[i+k] * DIRECT_NZYX_ELEM(src, 0, z, y, xx);
            m += kernel[i+k];
        }

        DIRECT_NZYX_ELEM(dest, 0, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = y + i;
            if (yy < 0 || yy >= src.ydim) continue;

            v += kernel[i+k] * DIRECT_NZYX_ELEM(dest, 0, z, yy, x);
            m += kernel[i+k];
        }

        DIRECT_NZYX_ELEM(temp, 0, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int zz = z + i;
            if (zz < 0 || zz >= src.zdim) continue;

            v += kernel[i+k] * DIRECT_NZYX_ELEM(temp, 0, zz, y, x);
            m += kernel[i+k];
        }

        DIRECT_NZYX_ELEM(dest, 0, z, y, x) = v/m;
    }
}

template<typename T>
void FilterHelper::separableGaussianWrap(const MultidimArray<T>& src, MultidimArray<T>& dest, double sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.reshape(src);

    std::vector<double> kernel(2*k+1);
    const double s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-0.5*i*i/s2);
    }

    MultidimArray<T> temp(src.zdim, src.ydim, src.xdim);

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int xx = (src.xdim + x + i) % src.xdim;

            v += kernel[i+k] * DIRECT_NZYX_ELEM(src, 0, z, y, xx);
            m += kernel[i+k];
        }

        DIRECT_NZYX_ELEM(dest, 0, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = (src.ydim + y + i) % src.ydim;

            v += kernel[i+k] * DIRECT_NZYX_ELEM(dest, 0, z, yy, x);
            m += kernel[i+k];
        }

        DIRECT_NZYX_ELEM(temp, 0, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int zz = (src.zdim + z + i) % src.zdim;

            v += kernel[i+k] * DIRECT_NZYX_ELEM(temp, 0, zz, y, x);
            m += kernel[i+k];
        }

        DIRECT_NZYX_ELEM(dest, 0, z, y, x) = v/m;
    }
}

#endif

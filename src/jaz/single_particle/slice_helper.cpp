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

#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/interpolation.h>

using namespace gravis;

void SliceHelper::affineTransform(const Image<RFLOAT>& img, d4Matrix A, Image<RFLOAT>& dest)
{
    d4Matrix Ai = A;
    Ai.invert();

    dest.data.resize(1, 1, img.data.ydim, img.data.xdim);

    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        d4Vector s0(x,y,0,1);
        d4Vector s1 = Ai * s0;

        DIRECT_A2D_ELEM(dest.data, y, x) = Interpolation::linearXY(img, s1.x, s1.y, 0);
    }
}

void SliceHelper::downsample(Image<RFLOAT>& img, Image<RFLOAT>& dest)
{
    double q = dest.data.xdim / (double) img.data.xdim;

    Image<RFLOAT> slice0(img.data.xdim, img.data.ydim, 1);

    FilterHelper::lowPassFilter(img, 0.9*q, q, slice0);
    subsample(slice0, dest);
}

void SliceHelper::downsampleSlices(const Image<RFLOAT>& img, Image<RFLOAT>& dest)
{
    double q = dest.data.xdim / (double) img.data.xdim;

    Image<RFLOAT> slice0(img.data.xdim, img.data.ydim, 1);
    Image<RFLOAT> slice1(dest.data.xdim, dest.data.ydim, 1);

    for (long int n = 0; n < img.data.ndim; n++)
    {
        std::cout << n << "/" << img.data.ndim << "\n";

        extractStackSlice(img, slice0, n);
        FilterHelper::lowPassFilter(slice0, 0.9*q, q, slice0);
        subsample(slice0, slice1);
        insertStackSlice(slice1, dest, n);
    }
}

void SliceHelper::downsampleSlicesReal(const Image<RFLOAT>& img, Image<RFLOAT>& dest)
{
    double q = dest.data.xdim / (double) img.data.xdim;

    Image<RFLOAT> slice0(img.data.xdim, img.data.ydim, 1);
    Image<RFLOAT> sliceT(img.data.xdim, img.data.ydim, 1);
    Image<RFLOAT> slice1(dest.data.xdim, dest.data.ydim, 1);

    for (long int n = 0; n < img.data.ndim; n++)
    {
        extractStackSlice(img, slice0, n);
        FilterHelper::separableGaussianXYZ(slice0, sliceT, 1.5/q);
        subsample(sliceT, slice1);
        insertStackSlice(slice1, dest, n);
    }
}

void SliceHelper::lowPassFilterSlicewise(Image<RFLOAT>& img, double maxFreq0, double maxFreq1)
{
    Image<RFLOAT> slice0(img.data.xdim, img.data.ydim, 1);

    for (long int n = 0; n < img.data.ndim; n++)
    {
        extractStackSlice(img, slice0, n);
        FilterHelper::lowPassFilter(slice0, maxFreq0, maxFreq1, slice0);
        insertStackSlice(slice0, img, n);
    }
}

void SliceHelper::lowPassFilterSlice(Image<RFLOAT>& img, long int n, double maxFreq0, double maxFreq1)
{
    Image<RFLOAT> slice0(img.data.xdim, img.data.ydim, 1);

    extractStackSlice(img, slice0, n);
    FilterHelper::lowPassFilter(slice0, maxFreq0, maxFreq1, slice0);
    insertStackSlice(slice0, img, n);
}

void SliceHelper::subsample(const Image<RFLOAT>& img, Image<RFLOAT>& dest)
{
    double q = img.data.xdim / (double) dest.data.xdim;

    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_A2D_ELEM(dest.data, y, x) = DIRECT_A2D_ELEM(img.data, (long int)(q*y + 0.5), (long int)(q*x + 0.5));
    }
}

void SliceHelper::avgPad(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest, double ratio)
{
    int padX = (int)(ratio * src.dimx);
    int padY = (int)(ratio * src.dimy);
    int padZ = (int)(ratio * src.dimz);

    double avg = 0.0;

    for (long int z = 0; z < src.dimz; z++)
    for (long int y = 0; y < src.dimy; y++)
    for (long int x = 0; x < src.dimx; x++)
    {
        avg += src(x,y,z);
    }

    avg /= src.dimx * src.dimy * src.dimz;

    dest.resize(src.dimx + 2*padX, src.dimy + 2*padY, src.dimz + 2*padZ);
    dest.fill(avg);

    for (long int z = 0; z < src.dimz; z++)
    for (long int y = 0; y < src.dimy; y++)
    for (long int x = 0; x < src.dimx; x++)
    {
        dest(x+padX, y+padY, z+padZ) = src(x,y,z);
    }
}

void SliceHelper::avgPad2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest, double ratio)
{
    int padX = (int)(ratio * src.data.xdim);
    int padY = (int)(ratio * src.data.ydim);

    double avg = 0.0;

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        avg += DIRECT_A2D_ELEM(src.data, y, x);
    }

    avg /= src.data.xdim * src.data.ydim;

    dest = Image<RFLOAT>(src.data.xdim + 2*padX, src.data.ydim + 2*padY);

    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_A2D_ELEM(dest.data, y, x) = avg;
    }

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        DIRECT_A2D_ELEM(dest.data, y+padY, x+padX) = DIRECT_A2D_ELEM(src.data, y, x);
    }
}

void SliceHelper::halveSpectrum2D(Image<Complex>& src, Image<Complex>& dest)
{
    dest = Image<Complex>(src.data.xdim/2 + 1, src.data.ydim);

    const int xo = src.data.xdim/2 + 1;
    const int yo = src.data.ydim/2 + 1;

    const int wd = dest.data.xdim;
    const int hd = dest.data.ydim;

    for (long int y = 0; y < hd; y++)
    for (long int x = 0; x < wd; x++)
    {
        /*if (x == 0)
        {
            DIRECT_A2D_ELEM(dest.data, y, 0) = DIRECT_A2D_ELEM(src.data, y, xo);
        }
        else if (xo + x < src.data.xdim)
        {
            DIRECT_A2D_ELEM(dest.data, y, x) = 0.5 * (DIRECT_A2D_ELEM(src.data, y, xo + x)
                                                      + DIRECT_A2D_ELEM(src.data, yo - (y - yo), xo - x));
        }
        else
        {
            DIRECT_A2D_ELEM(dest.data, y, x) = DIRECT_A2D_ELEM(src.data, yo - (y - yo), xo - x);
        }*/

        const int yr = (int)y;
        const int yw = (yr+yo)%hd;

        DIRECT_A2D_ELEM(dest.data, y, x) = DIRECT_A2D_ELEM(src.data, yw, xo+x);
    }
}

void SliceHelper::extractSpectralSlice(Image<Complex>& src, Image<RFLOAT>& dest, d3Matrix proj,
                                       d2Vector volCentImg, double oversample)
{
    const int wi = (double)dest.data.xdim;
    const int hi = (double)dest.data.ydim;

    const double wv = (double)src.data.xdim;
    const double hv = (double)src.data.ydim;
    const double dv = (double)src.data.zdim;

    const double wios = oversample*wi;
    const double hios = oversample*hi;

    const int wiosI = ((int)wios)/2 + 1;
    const int hiosI = ((int)hios);

    const int ciosX = ((int)wios)/2;
    const int ciosY = ((int)hios)/2;

    Image<Complex> dest2(wiosI,hiosI);
    Image<Complex> weight(wiosI,hiosI);

    d2Vector shift(volCentImg.x - ciosX, volCentImg.y - ciosY);

    for (long int y = 0; y < dest2.data.ydim; y++)
    for (long int x = 0; x < dest2.data.xdim; x++)
    {
        d3Vector pi((double)x/(double)(wiosI-1), 2.0*(double)y/(double)hiosI, 0.0);

        if (pi.y >= 1.0) pi.y = pi.y - 2.0;

        if (pi.norm2() > 1.0)
        {
            DIRECT_A2D_ELEM(dest2.data, y, x) = Complex(0,0);
            continue;
        }

        d3Vector pv = proj * pi;

        bool conj = false;

        if (pv.x < 0.0)
        {
            pv = -pv;
            conj = true;
        }

        if (pv.norm2() > 1.0)
        {
            DIRECT_A2D_ELEM(dest2.data, y, x) = Complex(0,0);
            continue;
        }

        double xxd = (wv-1) * pv.x;
        double yyd = hv * pv.y / 2.0;
        double zzd = dv * pv.z / 2.0;

        double ax = std::abs(xxd);
        double ay = std::abs(yyd);
        double az = std::abs(zzd);

        double phi = - PI * (pi.x * shift.x + pi.y * shift.y);
        Complex z0(cos(phi), sin(phi));

        if (ax < 1.0 && ay < 1.0 && az < 1.0)
        {
            DIRECT_A2D_ELEM(weight.data, y, x) = z0 * Complex((1.0 - ax) * (1.0 - ay) * (1.0 - az), 0.0);
        }
        else
        {
            DIRECT_A2D_ELEM(weight.data, y, x) = Complex(0.0, 0.0);
        }

        if (yyd < 0.0) yyd += hv;
        if (zzd < 0.0) zzd += dv;

        if (conj)
        {
            DIRECT_A2D_ELEM(dest2.data, y, x) = z0 * Interpolation::linearFFTW3D(src, xxd, yyd, zzd).conj();
        }
        else
        {
            DIRECT_A2D_ELEM(dest2.data, y, x) = z0 * Interpolation::linearFFTW3D(src, xxd, yyd, zzd);
        }
     }

    Image<RFLOAT> dest2r = Image<RFLOAT>(2 * (dest2.data.xdim - 1), dest2.data.ydim);
    Image<RFLOAT> weightr = Image<RFLOAT>(2 * (dest2.data.xdim - 1), dest2.data.ydim);

    FourierTransformer ft;
    ft.inverseFourierTransform(dest2.data, dest2r());
    CenterFFT(dest2r.data, true);

    FourierTransformer ftw;
    ftw.inverseFourierTransform(weight.data, weightr());
    CenterFFT(weightr.data, true);

    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_A2D_ELEM(dest.data, y, x) = DIRECT_A2D_ELEM(dest2r.data, y, x)
                                         / DIRECT_A2D_ELEM(weightr.data, y, x);
    }
}

void SliceHelper::insertSpectralSlices(
        std::vector<Image<RFLOAT> >& src,
        std::vector<gravis::d3Matrix> proj,
        std::vector<gravis::d2Vector> volCentImg,
        Image<Complex>& dest, double thickness, double thicknessSlope, double imgPad)
{
    const double wv = dest.data.xdim;
    const double hv = dest.data.ydim;
    const double dv = dest.data.zdim;

    const double wir = src[0].data.xdim;
    const double hir = src[0].data.ydim;

    const int ic = src.size();

    std::vector<Image<Complex> > srcSpectra(ic);
    std::vector<d2Vector> shifts(ic);
    std::vector<double> thz(ic);

    Image<RFLOAT> img;

    for (int i = 0; i < ic; i++)
    {
        avgPad2D(src[i], img, imgPad);

        FourierTransformer ft;
        CenterFFT(img.data, false);
        ft.FourierTransform(img(), srcSpectra[i].data, true);

        shifts[i] = d2Vector(volCentImg[i].x - wir/2.0, volCentImg[i].y - hir/2.0);

        thz[i] = 0.5*d3Vector(wv*proj[i](2,0), hv*proj[i](2,1), dv*proj[i](2,2)).length();
    }

    const double wif = srcSpectra[0].data.xdim;
    const double hif = srcSpectra[0].data.ydim;

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        d3Vector pv((double)x/(wv-1), 2.0*(double)y/hv, 2.0*(double)z/dv);

        if (pv.y > 1.0) pv.y = pv.y - 2.0;
        if (pv.z > 1.0) pv.z = pv.z - 2.0;

        if (pv.norm2() >= 1.0)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = Complex(0,0);
            continue;
        }

        const double r = sqrt(pv.x*pv.x + pv.z*pv.z);

        Complex zs(0.0, 0.0);
        double wgh = 0.0;

        for (int i = 0; i < ic; i++)
        {
            d3Vector pi3 = proj[i] * pv;

            if (pi3.x*pi3.x + pi3.y*pi3.y >= 1.0)
            {
                continue;
            }

            const double za = thz[i] * std::abs(pi3.z);

            const double th_r = thickness + r * thz[i] * thicknessSlope;
            if (za > th_r) continue;

            bool conj = false;

            if (pi3.x < 0.0)
            {
                pi3 = -pi3;
                conj = true;
            }

            double xi = (wif-1) * pi3.x;
            double yi = hif * pi3.y / 2.0;



            if (yi < 0.0) yi += hif;

            double phi = PI * (pi3.x * shifts[i].x + pi3.y * shifts[i].y);
            Complex z0(cos(phi), sin(phi));

            //double wgi = (1.0 - za/th_r) * (thickness / th_r);
            double wgi = 1.0 - za/th_r;
            Complex zz = z0 * Interpolation::linearFFTW2D(srcSpectra[i], xi, yi);

            if (conj)
            {
                zz = zz.conj();
            }

            zs += wgi * zz;
            wgh += wgi;
        }

        if (wgh > 1.0) zs /= wgh;

        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = zs;
    }
}


void SliceHelper::insertWeightedSpectralSlices(
        std::vector<Image<RFLOAT> >& src,
        std::vector<gravis::d3Matrix> proj,
        std::vector<gravis::d2Vector> volCentImg,
        std::vector<double> imgWeights,
        Image<Complex>& dest, double thickness, double imgPad)
{
    const double wv = dest.data.xdim;
    const double hv = dest.data.ydim;
    const double dv = dest.data.zdim;

    const double wir = src[0].data.xdim;
    const double hir = src[0].data.ydim;

    const int ic = src.size();

    std::vector<Image<Complex> > srcSpectra(ic);
    std::vector<d2Vector> shifts(ic);
    std::vector<double> thz(ic);

    Image<RFLOAT> img;

    for (int i = 0; i < ic; i++)
    {
        avgPad2D(src[i], img, imgPad);

        FourierTransformer ft;
        CenterFFT(img.data, false);
        ft.FourierTransform(img(), srcSpectra[i].data, true);

        shifts[i] = d2Vector(volCentImg[i].x - wir/2.0, volCentImg[i].y - hir/2.0);

        thz[i] = 0.5*d3Vector(wv*proj[i](2,0), hv*proj[i](2,1), dv*proj[i](2,2)).length();
    }

    const double wif = srcSpectra[0].data.xdim;
    const double hif = srcSpectra[0].data.ydim;

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        d3Vector pv((double)x/(wv-1), 2.0*(double)y/hv, 2.0*(double)z/dv);

        if (pv.y > 1.0) pv.y = pv.y - 2.0;
        if (pv.z > 1.0) pv.z = pv.z - 2.0;

        if (pv.norm2() >= 1.0)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = Complex(0,0);
            continue;
        }

        Complex zs(0.0, 0.0);
        double wgh = 0.0;

        for (int i = 0; i < ic; i++)
        {
            d3Vector pi3 = proj[i] * pv;

            if (pi3.x*pi3.x + pi3.y*pi3.y >= 1.0)
            {
                continue;
            }

            const double za = thz[i] * std::abs(pi3.z);

            const double th_r = thickness;
            if (za > th_r) continue;

            bool conj = false;

            if (pi3.x < 0.0)
            {
                pi3 = -pi3;
                conj = true;
            }

            double xi = (wif-1) * pi3.x;
            double yi = hif * pi3.y / 2.0;



            if (yi < 0.0) yi += hif;

            double phi = PI * (pi3.x * shifts[i].x + pi3.y * shifts[i].y);
            Complex z0(cos(phi), sin(phi));

            double wgi = imgWeights[i] * (1.0 - za/th_r);
            Complex zz = z0 * Interpolation::linearFFTW2D(srcSpectra[i], xi, yi);

            if (conj)
            {
                zz = zz.conj();
            }

            zs += wgi * zz;
            wgh += wgi;
        }

        if (wgh > 1.0) zs /= wgh;

        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = zs;
    }
}

void SliceHelper::extractStackSlice(const Image<RFLOAT>& src, Image<RFLOAT>& dest, long int s)
{
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim)
    {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = DIRECT_NZYX_ELEM(src.data, s, 0, y, x);
    }
}

void SliceHelper::extractStackSlices(const Image<double>& src, Image<RFLOAT>& dest, long int s)
{
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim)
    {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int n = 0; n < dest.data.ndim; n++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, 0, y, x) = DIRECT_NZYX_ELEM(src.data, n+s, 0, y, x);
    }
}

void SliceHelper::extractStackSlices(const Image<float>& src, Image<RFLOAT>& dest, long int s)
{
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim)
    {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int n = 0; n < dest.data.ndim; n++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, 0, y, x) = (RFLOAT) DIRECT_NZYX_ELEM(src.data, n+s, 0, y, x);
    }
}

Image<RFLOAT> SliceHelper::getStackSlice(const Image<RFLOAT> &src, long n)
{
    const long int w = src().xdim;
    const long int h = src().ydim;

    Image<RFLOAT> out(w,h);

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        DIRECT_NZYX_ELEM(out.data, 0, 0, y, x) = DIRECT_NZYX_ELEM(src.data, n, 0, y, x);
    }

    return out;
}

void SliceHelper::insertStackSlice(const Image<double>& src, Image<double>& dest, long int s)
{
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim)
    {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, s, 0, y, x) = DIRECT_NZYX_ELEM(src.data, 0, 0, y, x);
    }
}

void SliceHelper::insertStackSlice(const Image<float>& src, Image<float>& dest, long int s)
{
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim)
    {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, s, 0, y, x) = DIRECT_NZYX_ELEM(src.data, 0, 0, y, x);
    }
}

void SliceHelper::insertZSlice(const Image<double>& src, Image<double>& dest, long int s)
{
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim)
    {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, s, y, x) = DIRECT_NZYX_ELEM(src.data, 0, 0, y, x);
    }
}

void SliceHelper::insertZSlice(const Image<float>& src, Image<float>& dest, long int s)
{
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim)
    {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, s, y, x) = DIRECT_NZYX_ELEM(src.data, 0, 0, y, x);
    }
}

Image<double> SliceHelper::consolidate(const std::vector<Image<double> >& src, bool toN)
{
    const int w = src[0].data.xdim;
    const int h = src[0].data.ydim;
    const int ic = src.size();

    const int zc = toN? 1 : ic;
    const int nc = toN? ic : 1;

    Image<double> out(w,h,zc,nc);

    for (int i = 0; i < ic; i++)
    {
        if (src[i].data.xdim != w || src[i].data.ydim != h)
        {
            REPORT_ERROR("SliceHelper::consolidate(): images are of unequal size.\n");
        }

        if (toN) insertStackSlice(src[i], out, i);
        else insertZSlice(src[i], out, i);
    }

    return out;
}

Image<float> SliceHelper::consolidate(const std::vector<Image<float> >& src, bool toN)
{
    const int w = src[0].data.xdim;
    const int h = src[0].data.ydim;
    const int ic = src.size();

    Image<float> out(w,h,1,ic);

    for (int i = 0; i < ic; i++)
    {
        if (src[i].data.xdim != w || src[i].data.ydim != h)
        {
            REPORT_ERROR("SliceHelper::consolidate(): images are of unequal size.\n");
        }

        if (toN) insertStackSlice(src[i], out, i);
        else insertZSlice(src[i], out, i);
    }

    return out;
}

void SliceHelper::stat(const Image<RFLOAT>& img)
{
    std::cout << "xdim: " << img.data.xdim << "\n";
    std::cout << "ydim: " << img.data.ydim << "\n";
    std::cout << "zdim: " << img.data.zdim << "\n";
    std::cout << "ndim: " << img.data.ndim << "\n";
}

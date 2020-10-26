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

#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/index_sort.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/math/tensor2x2.h>
#include <src/jaz/single_particle/interpolation.h>

#include <limits>

extern "C"
{
        #include <src/jaz/single_particle/d3x3/dsyev2.h>
}

using namespace gravis;

void FilterHelper::separableGaussianFreq(
        const MultidimArray<Complex> &src,
        MultidimArray<Complex> &dest, double sigma, int k)
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

    MultidimArray<Complex> temp(src.zdim, src.ydim, src.xdim);

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        Complex v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            long xx = x + i;

            bool conj = false;

            if (xx < 0)
            {
                xx = -xx;
                conj = true;
            }
            else if (xx >= src.xdim)
            {
                xx = 2*src.xdim - 1 - xx;
                conj = true;
            }

            long yy = conj? (src.ydim - y) % src.ydim : y;
            long zz = conj? (src.zdim - z) % src.zdim : z;

            v += kernel[i+k] * DIRECT_NZYX_ELEM(src, 0, zz, yy, xx);
            m += kernel[i+k];

        }

        DIRECT_NZYX_ELEM(dest, 0, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        Complex v = 0;
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
        Complex v = 0;
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

void FilterHelper::separableGaussianFreqXY(
        const MultidimArray<Complex> &src,
        MultidimArray<Complex> &dest, double sigma, int k)
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

    MultidimArray<Complex> temp(src.zdim, src.ydim, src.xdim);

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        Complex v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            long xx = x + i;

            bool conj = false;

            if (xx < 0)
            {
                xx = -xx;
                conj = true;
            }
            else if (xx >= src.xdim)
            {
                xx = 2*src.xdim - 2 - xx;
                conj = true;
            }

            long yy = conj? (src.ydim - y) % src.ydim : y;
            long zz = conj? (src.zdim - z) % src.zdim : z;

            Complex vv = DIRECT_NZYX_ELEM(src, 0, zz, yy, xx);

            if (conj) vv = vv.conj();

            v += kernel[i+k] * vv;
            m += kernel[i+k];

        }

        DIRECT_NZYX_ELEM(temp, 0, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        Complex v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = (src.ydim + y + i) % src.ydim;

            v += kernel[i+k] * DIRECT_NZYX_ELEM(temp, 0, z, yy, x);
            m += kernel[i+k];
        }

        DIRECT_NZYX_ELEM(dest, 0, z, y, x) = v/m;
    }
}

void FilterHelper::drawTestPattern(Image<RFLOAT>& img, int squareSize)
{
    for (long int z = 0; z < img.data.zdim; z++)
    for (long int y = 0; y < img.data.ydim; y++)
    for (long int x = 0; x < img.data.xdim; x++)
    {
        int xi = (int)(x/squareSize) % 2;
        int yi = (int)(y/squareSize) % 2;
        int zi = (int)(z/squareSize) % 2;
        int v = (xi + yi + zi) % 2;

        DIRECT_A3D_ELEM(img.data, z, y, x) = (RFLOAT) v;
    }
}

void FilterHelper::drawTestPattern(Volume<RFLOAT>& volume, int squareSize)
{
    for (size_t z = 0; z < volume.dimz; z++)
    for (size_t y = 0; y < volume.dimy; y++)
    for (size_t x = 0; x < volume.dimx; x++)
    {
        int xi = (int)(x/squareSize) % 2;
        int yi = (int)(y/squareSize) % 2;
        int zi = (int)(z/squareSize) % 2;
        int v = (xi + yi + zi) % 2;

        volume(x,y,z) = (RFLOAT) v;
    }
}

Image<RFLOAT> FilterHelper::expImg(Image<RFLOAT> &img, double scale)
{
    Image<RFLOAT> out = img;

    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
    {
        DIRECT_NZYX_ELEM(out.data, l, k, i, j) = exp(scale*DIRECT_NZYX_ELEM(img.data, l, k, i, j));
    }

    return out;
}

Image<RFLOAT> FilterHelper::logImg(Image<RFLOAT> &img, double thresh, double scale)
{
    Image<RFLOAT> out = img;

    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
    {
        double v = DIRECT_NZYX_ELEM(img.data, l, k, i, j);
        if (v < thresh) v = thresh;
        DIRECT_NZYX_ELEM(out.data, l, k, i, j) = log(scale*v);
    }

    return out;
}

Image<RFLOAT> FilterHelper::padCorner2D(Image<RFLOAT>& img, double factor)
{
    const int w0 = img.data.xdim;
    const int h0 = img.data.ydim;

    const int w1 = factor * w0;
    const int h1 = factor * h0;

    Image<RFLOAT> out(w1,h1);

    for (int y = 0; y < h1; y++)
    for (int x = 0; x < w1; x++)
    {
        int x1 = x < w1/2? x : x - w1;
        int y1 = y < h1/2? y : y - h1;

        if (x1 < w0/2 && y1 < h0/2 && x1 >= -w0/2 && y1 >= -h0/2)
        {
            int x0 = x1 < 0? x1 + w0 : x1;
            int y0 = y1 < 0? y1 + h0 : y1;

            DIRECT_A2D_ELEM(out.data, y, x) = DIRECT_A2D_ELEM(img.data, y0, x0);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, y, x) = 0.0;
        }
    }

    return out;
}

Image<Complex> FilterHelper::padCorner2D(Image<Complex>& img, double factor)
{
    const int w0 = img.data.xdim;
    const int h0 = img.data.ydim;

    const int w1 = factor * w0;
    const int h1 = factor * h0;

    Image<Complex> out(w1,h1);

    for (int y = 0; y < h1; y++)
    for (int x = 0; x < w1; x++)
    {
        int x1 = x;
        int y1 = y < h1/2? y : y - h1;

        if (x1 < w0/2 && y1 < h0/2 && x1 >= -w0/2 && y1 >= -h0/2)
        {
            int x0 = x1;
            int y0 = y1 < 0? y1 + h0 : y1;

            DIRECT_A2D_ELEM(out.data, y, x) = DIRECT_A2D_ELEM(img.data, y0, x0);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, y, x) = 0.0;
        }
    }

    return out;
}

Image<RFLOAT> FilterHelper::padCorner2D(const Image<RFLOAT>& img, int w, int h)
{
    const int w0 = img.data.xdim;
    const int h0 = img.data.ydim;

    Image<RFLOAT> out(w,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        int x1 = x < w/2? x : x - w;
        int y1 = y < h/2? y : y - h;

        if (x1 < w0/2 && y1 < h0/2 && x1 >= -w0/2 && y1 >= -h0/2)
        {
            int x0 = x1 < 0? x1 + w0 : x1;
            int y0 = y1 < 0? y1 + h0 : y1;

            DIRECT_A2D_ELEM(out.data, y, x) = DIRECT_A2D_ELEM(img.data, y0, x0);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, y, x) = 0.0;
        }
    }

    return out;
}

Image<RFLOAT> FilterHelper::cropCorner2D(const Image<RFLOAT>& img, int w, int h)
{
    const int w1 = img.data.xdim;
    const int h1 = img.data.ydim;

    if (w > w1 || h > h1) return img;

    Image<RFLOAT> out(w,h);

    for (int y = 0; y < h1; y++)
    for (int x = 0; x < w1; x++)
    {
        int x1 = x < w1/2? x : x - w1;
        int y1 = y < h1/2? y : y - h1;

        if (x1 < w/2 && y1 < h/2 && x1 >= -w/2 && y1 >= -h/2)
        {
            int x0 = x1 < 0? x1 + w : x1;
            int y0 = y1 < 0? y1 + h : y1;

            DIRECT_A2D_ELEM(out.data, y0, x0) = DIRECT_A2D_ELEM(img.data, y, x);
        }
    }

    return out;
}

Image<Complex> FilterHelper::cropCorner2D(const Image<Complex>& img, int w, int h)
{
    const int w1 = img.data.xdim;
    const int h1 = img.data.ydim;

    Image<Complex> out(w,h);

    for (int y = 0; y < h1; y++)
    for (int x = 0; x < w1; x++)
    {
        int x1 = x;
        int y1 = y < h1/2? y : y - h1;

        if (x1 < w && y1 < h/2 && y1 >= -h/2)
        {
            int x0 = x1;
            int y0 = y1 < 0? y1 + h : y1;

            DIRECT_A2D_ELEM(out.data, y0, x0) = DIRECT_A2D_ELEM(img.data, y, x);
        }
    }

    return out;
}

Image<RFLOAT> FilterHelper::zeroOutsideCorner2D(Image<RFLOAT> &img, double radius)
{
    const int w = img.data.xdim;
    const int h = img.data.ydim;

    const double rad2 = radius*radius;

    Image<RFLOAT> out(w,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        int xx = x < w/2? x : x - w;
        int yy = y < h/2? y : y - h;

        int r2 = xx*xx + yy*yy;

        if (r2 <= rad2)
        {
            DIRECT_A2D_ELEM(out.data, y, x) = DIRECT_A2D_ELEM(img.data, y, x);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, y, x) = 0.0;
        }
    }

    return out;
}

void FilterHelper::GaussianEnvelopeCorner2D(Image<RFLOAT> &img, double sigma)
{
    const int w = img.data.xdim;
    const int h = img.data.ydim;

    const double s2 = sigma * sigma;

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        double xx = x < w/2? x : x - w;
        double yy = y < h/2? y : y - h;

        double r2 = xx*xx + yy*yy;

        DIRECT_A2D_ELEM(img.data, y, x) *= exp(-0.5*r2/s2);
    }
}

Image<RFLOAT> FilterHelper::raisedCosEnvCorner2D(Image<RFLOAT> &img, double radIn, double radOut)
{
    const int w = img.data.xdim;
    const int h = img.data.ydim;

    Image<RFLOAT> out(w,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        double xx = x < w/2? x : x - w;
        double yy = y < h/2? y : y - h;

        double r = sqrt(xx*xx + yy*yy);

        if (r < radIn)
        {
            DIRECT_A2D_ELEM(out.data, y, x) = DIRECT_A2D_ELEM(img.data, y, x);
        }
        else if (r < radOut)
        {
            double t = (r - radIn)/(radOut - radIn);
            double a = 0.5 * (1.0 + cos(PI * t));

            DIRECT_A2D_ELEM(out.data, y, x) = a * DIRECT_A2D_ELEM(img.data, y, x);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, y, x) = 0.0;
        }
    }

    return out;
}

Image<Complex> FilterHelper::raisedCosEnvCorner2DFull(Image<Complex> &img, double radIn, double radOut)
{
	const int w = img.data.xdim;
    const int h = img.data.ydim;

    Image<Complex> out(w,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        double xx = x < w/2? x : x - w;
        double yy = y < h/2? y : y - h;

        double r = sqrt(xx*xx + yy*yy);

        if (r < radIn)
        {
            DIRECT_A2D_ELEM(out.data, y, x) = DIRECT_A2D_ELEM(img.data, y, x);
        }
        else if (r < radOut)
        {
            double t = (r - radIn)/(radOut - radIn);
            double a = 0.5 * (1.0 + cos(PI * t));

            DIRECT_A2D_ELEM(out.data, y, x) = a * DIRECT_A2D_ELEM(img.data, y, x);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, y, x) = 0.0;
        }
    }

    return out;
}

Image<RFLOAT> FilterHelper::raisedCosEnvCorner3D(Image<RFLOAT> &img, double radIn, double radOut)
{
    const int w = img.data.xdim;
	const int h = img.data.ydim;
	const int d = img.data.zdim;

    Image<RFLOAT> out(w,h,d);
	
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        double xx = x < w/2? x : x - w;
		double yy = y < h/2? y : y - h;
		double zz = z < d/2? z : z - d;

        double r = sqrt(xx*xx + yy*yy + zz*zz);

        if (r < radIn)
        {
            DIRECT_A3D_ELEM(out.data, z, y, x) = DIRECT_A3D_ELEM(img.data, z, y, x);
        }
        else if (r < radOut)
        {
            double t = (r - radIn)/(radOut - radIn);
            double a = 0.5 * (1.0 + cos(PI * t));

            DIRECT_A3D_ELEM(out.data, z, y, x) = a * DIRECT_A3D_ELEM(img.data, z, y, x);
        }
        else
        {
            DIRECT_A3D_ELEM(out.data, z, y, x) = 0.0;
        }
    }

    return out;

}

Image<RFLOAT> FilterHelper::raisedCosEnvFreq2D(const Image<RFLOAT>& img, double radIn, double radOut)
{
    const int w = img.data.xdim;
    const int h = img.data.ydim;

    Image<RFLOAT> out(w,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        double xx = x;
        double yy = y <= h/2? y : y - h;

        double r = sqrt(xx*xx + yy*yy);

        if (r < radIn)
        {
            DIRECT_A2D_ELEM(out.data, y, x) = DIRECT_A2D_ELEM(img.data, y, x);
        }
        else if (r < radOut)
        {
            double t = (r - radIn)/(radOut - radIn);
            double a = 0.5 * (1.0 + cos(PI * t));

            DIRECT_A2D_ELEM(out.data, y, x) = a * DIRECT_A2D_ELEM(img.data, y, x);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, y, x) = 0.0;
        }
    }

	return out;
}

Image<RFLOAT> FilterHelper::raisedCosEnvRingFreq2D(
		const Image<RFLOAT> &img, 
		double rad0, double rad1, double stepWidth)
{
	const int w = img.data.xdim;
    const int h = img.data.ydim;

    Image<RFLOAT> out(w,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        double xx = x;
        double yy = y <= h/2? y : y - h;

        double r = sqrt(xx*xx + yy*yy);
		double r0 = rad0 > 0.0? r - rad0 : stepWidth/2;
		double r1 = rad1 - r;
		
		double re = 2.0 * XMIPP_MIN(r0, r1) / stepWidth;
		
        if (re > 1.0)
        {
            DIRECT_A2D_ELEM(out.data, y, x) = DIRECT_A2D_ELEM(img.data, y, x);
        }
        else if (re > -1.0)
        {
            double t = (re + 1.0)/2.0;
            double a = 0.5 * (1.0 - cos(PI * t));

            DIRECT_A2D_ELEM(out.data, y, x) = a * DIRECT_A2D_ELEM(img.data, y, x);
        }
        else
        {
            DIRECT_A2D_ELEM(out.data, y, x) = 0.0;
        }
    }

	return out;	
}

void FilterHelper::lowPassFilter(Image<RFLOAT>& img, double maxFreq0, double maxFreq1, Image<RFLOAT>& dest)
{
    MultidimArray<Complex> imgFreq;
    FourierTransformer ft;
    ft.FourierTransform(img(), imgFreq, false);

    lowPassFilterSpectrum(imgFreq, maxFreq0, maxFreq1);

    FourierTransformer ft2;
    ft2.inverseFourierTransform(imgFreq, dest());
}

void FilterHelper::lowPassFilterSpectrum(MultidimArray<Complex>& spectrum, double maxFreq0, double maxFreq1)
{
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(spectrum)
    {
        double xi = j/(double)spectrum.xdim;
        double yi = 2.0*i/(double)spectrum.ydim;
        double zi = 2.0*k/(double)spectrum.zdim;

        if (yi > 1.0) yi = 2.0 - yi;
        if (zi > 1.0) zi = 2.0 - zi;

        double r = sqrt(xi*xi + yi*yi + zi*zi);

        if (r > maxFreq1)
        {
            DIRECT_A3D_ELEM(spectrum, k, i, j) = Complex(0.0, 0.0);
        }
        else if (r > maxFreq0)
        {
            const double t = (r - maxFreq0)/(maxFreq1 - maxFreq0);
            const double q = 0.5 * (cos(PI*t) + 1.0);
            DIRECT_A3D_ELEM(spectrum, k, i, j) *= q;
        }
    }
}

RFLOAT FilterHelper::averageValue(Image<RFLOAT>& img)
{
    RFLOAT sum;

    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
    {
        sum += DIRECT_NZYX_ELEM(img.data, l, k, i, j);
    }

    return sum / (double)(img.data.xdim * img.data.ydim * img.data.zdim * img.data.ndim);
}

RFLOAT FilterHelper::maxValue(Image<RFLOAT> &img)
{
    RFLOAT vMax = -std::numeric_limits<double>::max();

    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data)
    {
        RFLOAT v = DIRECT_NZYX_ELEM(img.data, l, k, i, j);
        if (v > vMax) vMax = v;
    }

    return vMax;
}

void FilterHelper::phaseFlip(Image<RFLOAT>& img, CTF& ctf, RFLOAT angpix, Image<RFLOAT>& dest)
{
    MultidimArray<Complex> imgFreq;
    FourierTransformer ft;
    ft.FourierTransform(img(), imgFreq, false);

    RFLOAT xs = (RFLOAT)img.data.xdim * angpix;
    RFLOAT ys = (RFLOAT)img.data.ydim * angpix;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imgFreq)
    {
        const int x = j;
        const int y = i < imgFreq.ydim/2? i : i - imgFreq.ydim;

        RFLOAT c = ctf.getCTF(x/xs, y/ys);

        if (c < 0)
        {
            DIRECT_A2D_ELEM(imgFreq, i, j) *= -1;
        }
    }

    if (dest.data.xdim != img.data.xdim || dest.data.ydim != img.data.ydim)
    {
        dest.data.resize(img.data);
    }

    FourierTransformer ft2;
    ft2.inverseFourierTransform(imgFreq, dest());
}

void FilterHelper::applyBeamTilt(Image<RFLOAT> &img, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
                                 RFLOAT lambda, RFLOAT Cs, RFLOAT angpix, int s, Image<RFLOAT>& dest)
{
    MultidimArray<Complex> imgFreq;
    FourierTransformer ft;
    ft.FourierTransform(img(), imgFreq, false);

    selfApplyBeamTilt(imgFreq, beamtilt_x, beamtilt_y, lambda, Cs, angpix, s);

    FourierTransformer ft2;
    ft2.inverseFourierTransform(imgFreq, dest());
}

void FilterHelper::modulate(Image<RFLOAT>& img, CTF& ctf, RFLOAT angpix, Image<RFLOAT>& dest)
{
    Image<Complex> imgFreq;
    FourierTransformer ft;
    ft.FourierTransform(img(), imgFreq(), false);

    modulate(imgFreq, ctf, angpix, dest);
}

void FilterHelper::modulate(Image<Complex>& imgFreq, CTF& ctf, RFLOAT angpix, Image<RFLOAT>& dest)
{
    const int w = imgFreq.data.xdim;
	const int h = imgFreq.data.ydim;

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix);
	
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imgFreq())
    {
        DIRECT_A2D_ELEM(imgFreq(), i, j) *= DIRECT_A2D_ELEM(ctfImg(), i, j);
    }

    if (dest.data.xdim != 2*(imgFreq.data.xdim-1) || dest.data.ydim != imgFreq.data.ydim)
    {
        dest.data.resize(imgFreq.data.ydim, 2*(imgFreq.data.xdim-1));
    }

    FourierTransformer ft2;
    ft2.inverseFourierTransform(imgFreq(), dest());
}

void FilterHelper::modulate(MultidimArray<Complex>& imgFreq, CTF& ctf, RFLOAT angpix)
{
	const int w = imgFreq.xdim;
	const int h = imgFreq.ydim;

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imgFreq)
    {
        DIRECT_A2D_ELEM(imgFreq, i, j) *= DIRECT_A2D_ELEM(ctfImg(), i, j);
    }
}

void FilterHelper::drawCtf(CTF &ctf, RFLOAT angpix, Image<Complex> &dest)
{
	const int w = dest.data.xdim;
	const int h = dest.data.ydim;

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(dest())
    {
        DIRECT_A2D_ELEM(dest(), i, j) = DIRECT_A2D_ELEM(ctfImg(), i, j);
    }
}

void FilterHelper::wienerFilter(Image<RFLOAT>& img, CTF& ctf, RFLOAT angpix, RFLOAT eps, RFLOAT Bfac, Image<RFLOAT>& dest)
{
    MultidimArray<Complex> imgFreq;
    FourierTransformer ft;
    ft.FourierTransform(img(), imgFreq, false);

    RFLOAT xs = (RFLOAT)img.data.xdim * angpix;
    RFLOAT ys = (RFLOAT)img.data.ydim * angpix;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imgFreq)
    {
        const int x = j;
        const int y = i < imgFreq.ydim/2? i : i - imgFreq.ydim;

        RFLOAT c;
        if (Bfac > 0.0) c = ctf.getCTF(x/xs, y/ys) * exp(-Bfac*(x*x + y*y)/4.0);
        else c = ctf.getCTF(x/xs, y/ys);

        DIRECT_A2D_ELEM(imgFreq, i, j) = (c * DIRECT_A2D_ELEM(imgFreq, i, j))/(c*c + eps);
    }

    if (dest.data.xdim != img.data.xdim || dest.data.ydim != img.data.ydim)
    {
        dest.data.resize(img.data);
    }

    FourierTransformer ft2;
    ft2.inverseFourierTransform(imgFreq, dest());
}

void FilterHelper::richardsonLucy(Image<RFLOAT>& img, CTF& ctf, RFLOAT angpix, RFLOAT eps, int iterations, Image<RFLOAT>& dest)
{
    const int w = img.data.xdim;
    const int h = img.data.ydim;

    Image<RFLOAT> img0(w,h,1,1), img1(w,h,1,1), img1M(w,h,1,1), imgR(w,h,1,1), imgRM(w,h,1,1);

    double vmin = 0;
    double Bfac = (double)w/4.0;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img.data)
    {
        double v0 = DIRECT_A2D_ELEM(img.data, i, j);
        if (v0 < vmin) vmin = v0;
    }

    vmin -= 10;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(img.data)
    {
        DIRECT_A2D_ELEM(img0.data, i, j) = DIRECT_A2D_ELEM(img.data, i, j) + vmin;
    }

    wienerFilter(img0, ctf, angpix, eps, Bfac, img1);

    VtkHelper::writeVTK(img1, "rl_it0.vtk");

    for (int it = 0; it < iterations; it++)
    {
        // img1 = img1 * conv(psf, img / conv(psf, img1) )
        //      = img1 * IFT( ctf * FT(img / IFT( ctf * FT(img1) ) ) )
        //      = img1 * ctf_mod( img / ctf_mod(img1) )

        modulate(img1, ctf, angpix, img1M);
        wienerDivide(img0, img1M, eps, imgR);
        modulate(imgR, ctf, angpix, imgRM);
        multiply(imgRM, img1, img1);

        std::stringstream sts;
        sts << (it+1);
        std::string fn;
        sts >> fn;

        VtkHelper::writeVTK(img1, "rl_it"+fn+".vtk");
    }
}

void FilterHelper::rampFilter(Image<RFLOAT>& img, RFLOAT s0, RFLOAT t1, double ux, double uy, Image<RFLOAT>& dest)
{
    MultidimArray<Complex> imgFreq;
    FourierTransformer ft;
    ft.FourierTransform(img(), imgFreq, false);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imgFreq)
    {
        const int x = j;
        const int y = i < imgFreq.ydim/2? i : i - imgFreq.ydim;

        RFLOAT t = std::abs(x*ux + y*uy);
        RFLOAT s = t < t1? (s0 + (1-s0)*t/t1) : 1.0;

        DIRECT_A2D_ELEM(imgFreq, i, j) = s * DIRECT_A2D_ELEM(imgFreq, i, j);
    }

    if (dest.data.xdim != img.data.xdim || dest.data.ydim != img.data.ydim)
    {
        dest.data.resize(img.data);
    }

    FourierTransformer ft2;
    ft2.inverseFourierTransform(imgFreq, dest());
}

void FilterHelper::rampFilter3D(Image<Complex>& img, RFLOAT s0, RFLOAT t1, double tx, double ty, double tz)
{
    d3Vector ta(tx,ty,tz);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img.data)
    {
        const int x = j;
        const int y = i < img.data.ydim/2? i : i - img.data.ydim;
        const int z = k < img.data.zdim/2? k : k - img.data.zdim;

        d3Vector p(x,y,z);
        d3Vector q = p - p.dot(ta)*ta;
        double t = q.length();

        RFLOAT s = t < t1? (s0 + (1-s0)*t/t1) : 1.0;

        DIRECT_A3D_ELEM(img.data, k, i, j) = s * DIRECT_A3D_ELEM(img.data, k, i, j);
    }
}

void FilterHelper::doubleRampFilter3D(Image<Complex>& img, RFLOAT s0, RFLOAT t1, double tx, double ty, double tz)
{
    d3Vector ta(tx,ty,tz);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img.data)
    {
        const int x = j;
        const int y = i < img.data.ydim/2? i : i - img.data.ydim;
        const int z = k < img.data.zdim/2? k : k - img.data.zdim;

        d3Vector p(x,y,z);
        d3Vector q = p - p.dot(ta)*ta;
        double t = q.length();

        RFLOAT s = t < t1? (s0 + (1-s0)*t/t1) : 1.0 + t1 - t;
        if (s < 0) s = 0;

        DIRECT_A3D_ELEM(img.data, k, i, j) = s * DIRECT_A3D_ELEM(img.data, k, i, j);
    }
}

void FilterHelper::getPhase(const Image<Complex> &img, Image<RFLOAT> &dest)
{
    const long w = img.data.xdim;
    const long h = img.data.ydim;
    const long d = img.data.zdim;

    dest = Image<RFLOAT>(w,h,d);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        if (DIRECT_NZYX_ELEM(img.data, 0, z, y, x).norm() > 0)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(img.data, 0, z, y, x).arg();
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = 0;
        }
    }
}

void FilterHelper::getAbs(const Image<Complex> &img, Image<RFLOAT> &dest)
{
    const long w = img.data.xdim;
    const long h = img.data.ydim;
    const long d = img.data.zdim;

    dest = Image<RFLOAT>(w,h,d);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(img.data, 0, z, y, x).abs();
    }
}

void FilterHelper::getReal(const Image<Complex> &img, Image<RFLOAT> &dest)
{
    const long w = img.data.xdim;
    const long h = img.data.ydim;
    const long d = img.data.zdim;

    dest = Image<RFLOAT>(w,h,d);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(img.data, 0, z, y, x).real;
    }
}

void FilterHelper::getImag(const Image<Complex> &img, Image<RFLOAT> &dest)
{
    const long w = img.data.xdim;
    const long h = img.data.ydim;
    const long d = img.data.zdim;

    dest = Image<RFLOAT>(w,h,d);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(img.data, 0, z, y, x).imag;
    }
}

void FilterHelper::powerSpectrum2D(Image<RFLOAT>& img, Volume<RFLOAT>& spectrum)
{
    MultidimArray<Complex> imgFreq;
    FourierTransformer ft;
    ft.FourierTransform(img(), imgFreq, false);

    spectrum.resize(imgFreq.xdim, imgFreq.ydim, 1);

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(imgFreq)
    {
        Complex z = DIRECT_A2D_ELEM(imgFreq, i, j);

        spectrum(j,i,0) = z.abs();
    }
}

void FilterHelper::equiphaseAverage2D(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest)
{
    int n = src.dimx;
    std::vector<RFLOAT> val(n), wgh(n);

    for (long int i = 0; i < n; i++)
    {
        val[i] = 0.0;
        wgh[i] = 0.0;
    }

    for (long int y = 0; y < src.dimy; y++)
    for (long int x = 0; x < src.dimx; x++)
    {
        double id;

        if (y < src.dimy/2)
        {
            id = sqrt(x*x + y*y);
        }
        else
        {
            id = sqrt(x*x + (src.dimy - y)*(src.dimy - y));
        }

        int i = (int)id;
        double f = id - i;

        if (i >= 0 && i < n)
        {
            val[i] += (1.0 - f) * src(x,y,0);
            wgh[i] += (1.0 - f);
        }

        if (i >= -1 && i < n-1)
        {
            val[i+1] += f * src(x,y,0);
            wgh[i+1] += f;
        }
    }

    for (long int i = 0; i < n; i++)
    {
        if (wgh[i] > 0.0)
        {
            val[i] /= wgh[i];
        }
    }

    dest.resize(src);

    for (long int y = 0; y < src.dimy; y++)
    for (long int x = 0; x < src.dimx; x++)
    {
        double id;

        if (y < src.dimy/2)
        {
            id = sqrt(x*x + y*y);
        }
        else
        {
            id = sqrt(x*x + (src.dimy - y)*(src.dimy - y));
        }

        int i = (int)id;
        double f = id - i;

        if (i >= 0 && i < n-1)
        {
            dest(x,y,0) = (1.0 - f) * val[i] + f * val[i+1];
        }
    }
}

void FilterHelper::threshold(Image<RFLOAT>& src, RFLOAT t, Image<RFLOAT>& dest)
{
    for (long int n = 0; n < src.data.ndim; n++)
    for (long int z = 0; z < src.data.zdim; z++)
    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        if (DIRECT_NZYX_ELEM(src.data, n, z, y, x) > t)
        {
            DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = 1.0;
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = 0.0;
        }
    }
}

void FilterHelper::fill(Image<RFLOAT>& dest, RFLOAT v)
{
    for (long int n = 0; n < dest.data.ndim; n++)
    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = v;
    }
}

void FilterHelper::linearTransform(Image<RFLOAT>& src, RFLOAT m, RFLOAT q, Image<RFLOAT>& dest)
{
    for (long int n = 0; n < src.data.ndim; n++)
    for (long int z = 0; z < src.data.zdim; z++)
    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = m * DIRECT_NZYX_ELEM(src.data, n, z, y, x) + q;
    }
}

void FilterHelper::linearCombination(Image<RFLOAT>& src0, Image<RFLOAT>& src1, RFLOAT a0, RFLOAT a1, Image<RFLOAT>& dest)
{
    for (long int n = 0; n < src0.data.ndim; n++)
    for (long int z = 0; z < src0.data.zdim; z++)
    for (long int y = 0; y < src0.data.ydim; y++)
    for (long int x = 0; x < src0.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = a0 * DIRECT_NZYX_ELEM(src0.data, n, z, y, x) + a1 * DIRECT_NZYX_ELEM(src1.data, n, z, y, x);
    }
}

void FilterHelper::linearCombination(const Volume<RFLOAT>& src0, const Volume<RFLOAT>& src1, RFLOAT a0, RFLOAT a1, Volume<RFLOAT>& dest)
{
    for (long int z = 0; z < src0.dimz; z++)
    for (long int y = 0; y < src0.dimy; y++)
    for (long int x = 0; x < src0.dimx; x++)
    {
        dest(x,y,z) = a0 * src0(x,y,z) + a1 * src1(x,y,z);
    }
}

void FilterHelper::sumUp(const std::vector<Image<RFLOAT> > & src, Image<RFLOAT> &dest)
{
    const int w = src[0].data.xdim;
    const int h = src[0].data.ydim;
    const int d = src[0].data.zdim;
    const int m = src[0].data.ndim;
    const int ic = src.size();

    dest = Image<RFLOAT>(w,h,d,m);
    dest.data.initZeros();

    for (long int i = 0; i < ic; i++)
    {
        if (   src[i].data.xdim != w || src[i].data.ydim != h
            || src[i].data.zdim != d || src[i].data.ndim != m)
        {
            REPORT_ERROR("FilterHelper::sumUp(): image dimension mismatch.\n");
        }

        for (long int n = 0; n < m; n++)
        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++)
        {
            DIRECT_NZYX_ELEM(dest.data, n, z, y, x) += DIRECT_NZYX_ELEM(src[i].data, n, z, y, x);
        }
    }
}

double FilterHelper::L1distance(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0, int y0, int w, int h)
{
    double d = 0.0;

    if (w < 0) w = i0.data.xdim;
    if (h < 0) h = i0.data.ydim;

    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++)
    {
        RFLOAT v0 = DIRECT_NZYX_ELEM(i0.data, n, z, y, x);
        RFLOAT v1 = DIRECT_NZYX_ELEM(i1.data, n, z, y, x);

        double di = v1 - v0;
        d += std::abs(di);
    }

    return d;
}

double FilterHelper::L2distance(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0, int y0, int w, int h)
{
    double d = 0.0;

    if (w < 0) w = i0.data.xdim;
    if (h < 0) h = i0.data.ydim;

    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++)
    {
        RFLOAT v0 = DIRECT_NZYX_ELEM(i0.data, n, z, y, x);
        RFLOAT v1 = DIRECT_NZYX_ELEM(i1.data, n, z, y, x);

        double di = v1 - v0;
        d += di*di;
    }

    return d;
}

double FilterHelper::NCC(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0, int y0, int w, int h)
{
    double d = 0.0;

    if (w < 0) w = i0.data.xdim;
    if (h < 0) h = i0.data.ydim;

    double mu0 = 0.0, mu1 = 0.0, cnt = 0.0;

    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++)
    {
        RFLOAT v0 = DIRECT_NZYX_ELEM(i0.data, n, z, y, x);
        RFLOAT v1 = DIRECT_NZYX_ELEM(i1.data, n, z, y, x);

        mu0 += v0;
        mu1 += v1;
        cnt += 1.0;
    }

    mu0 /= cnt;
    mu1 /= cnt;

    double sig0 = 0.0, sig1 = 0.0;

    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++)
    {
        RFLOAT v0 = DIRECT_NZYX_ELEM(i0.data, n, z, y, x) - mu0;
        RFLOAT v1 = DIRECT_NZYX_ELEM(i1.data, n, z, y, x) - mu1;

        sig0 += v0*v0;
        sig1 += v1*v1;
    }

    sig0 = sqrt(sig0/(cnt - 1.0));
    sig1 = sqrt(sig1/(cnt - 1.0));

    double ncc = 0.0;

    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++)
    {
        RFLOAT v0 = (DIRECT_NZYX_ELEM(i0.data, n, z, y, x) - mu0);
        RFLOAT v1 = (DIRECT_NZYX_ELEM(i1.data, n, z, y, x) - mu1);

        ncc += v0*v1;
    }

    ncc /= sig0*sig1*cnt;

    return ncc;
}

void FilterHelper::multiply(Image<RFLOAT>& i0, Image<RFLOAT>& i1, Image<RFLOAT>& dest)
{
    dest = Image<RFLOAT>(i0.data.xdim, i0.data.ydim, i0.data.zdim, i0.data.ndim);

    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = 0; y < i0.data.ydim; y++)
    for (long int x = 0; x < i0.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = DIRECT_NZYX_ELEM(i0.data, n, z, y, x) * DIRECT_NZYX_ELEM(i1.data, n, z, y, x);
    }
}

void FilterHelper::multiply(Image<Complex>& i0, Image<Complex>& i1, Image<Complex>& dest)
{
    dest = Image<Complex>(i0.data.xdim, i0.data.ydim, i0.data.zdim, i0.data.ndim);

    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = 0; y < i0.data.ydim; y++)
    for (long int x = 0; x < i0.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = DIRECT_NZYX_ELEM(i0.data, n, z, y, x) * DIRECT_NZYX_ELEM(i1.data, n, z, y, x);
    }
}

void FilterHelper::wienerDivide(Image<RFLOAT>& num, Image<RFLOAT>& denom, RFLOAT eps, Image<RFLOAT>& dest)
{
    for (long int n = 0; n < num.data.ndim; n++)
    for (long int z = 0; z < num.data.zdim; z++)
    for (long int y = 0; y < num.data.ydim; y++)
    for (long int x = 0; x < num.data.xdim; x++)
    {
        RFLOAT d = DIRECT_NZYX_ELEM(denom.data, n, z, y, x);
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = d * DIRECT_NZYX_ELEM(num.data, n, z, y, x) / (d*d + eps);
    }
}

void FilterHelper::divide(Image<Complex>& num, Volume<RFLOAT>& denom, RFLOAT eps, Image<Complex>& dest)
{
    for (long int n = 0; n < num.data.ndim; n++)
    for (long int z = 0; z < num.data.zdim; z++)
    for (long int y = 0; y < num.data.ydim; y++)
    for (long int x = 0; x < num.data.xdim; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = DIRECT_NZYX_ELEM(num.data, n, z, y, x) / (denom(x,y,z) + eps);
    }
}

void FilterHelper::divide(Image<Complex>& num, Image<RFLOAT>& denom, RFLOAT eps, Image<Complex>& dest)
{
    FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(num.data)
    {
        DIRECT_NZYX_ELEM(dest.data, l, k, i, j) = DIRECT_NZYX_ELEM(num.data, l, k, i, j)
                                                  / (DIRECT_NZYX_ELEM(denom.data, l, k, i, j) + eps);
    }
}

void FilterHelper::divideExcessive(Image<Complex>& num, Volume<RFLOAT>& denom, RFLOAT theta, Image<Complex>& dest)
{
    for (long int n = 0; n < num.data.ndim; n++)
    for (long int z = 0; z < num.data.zdim; z++)
    for (long int y = 0; y < num.data.ydim; y++)
    for (long int x = 0; x < num.data.xdim; x++)
    {
        RFLOAT t = denom(x,y,z)/theta;

        if (t > 1)
        {
            DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = DIRECT_NZYX_ELEM(num.data, n, z, y, x) / t;
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, n, z, y, x) = DIRECT_NZYX_ELEM(num.data, n, z, y, x);
        }

    }
}

void FilterHelper::wienerDeconvolve(Image<Complex>& num, Image<Complex>& denom, RFLOAT theta, Image<Complex>& dest)
{
    for (long int z = 0; z < num.data.zdim; z++)
    for (long int y = 0; y < num.data.ydim; y++)
    for (long int x = 0; x < num.data.xdim; x++)
    {
        Complex zz = DIRECT_NZYX_ELEM(denom.data, 0, z, y, x);
        Complex z0 = DIRECT_NZYX_ELEM(num.data, 0, z, y, x);

        /*std::cout << "z0 = " << z0.real << " + " << z0.imag << " * i\n";
        std::cout << "zz = " << zz.real << " + " << zz.imag << " * i\n";
        std::cout << "zzB * z0 = " << (zz.conj() * z0).real << " + " << (zz.conj() * z0).imag << " * i\n";
        std::cout << "((zz.conj() * zz).real + theta) = " << ((zz.conj() * zz).real + theta) << " * i\n";*/

        //DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = (zz.conj() * z0) / ((zz.conj() * zz).real + theta);

        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = (zz.real * z0) / (zz.real * zz.real + theta);

        /*RFLOAT t = zz.abs()/theta;

        if (t > 1)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(num.data, 0, z, y, x) / t;
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(num.data, 0, z, y, x);
        }*/
    }
}

void FilterHelper::extract2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                            long int x0, long int y0,
                            long int w, long int h)
{
    dest = Image<RFLOAT>(w,h);

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        long int xx = x0 + x;
        long int yy = y0 + y;

        if (   xx >= 0 && xx < src.data.xdim
            && yy >= 0 && yy < src.data.ydim)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = DIRECT_NZYX_ELEM(src.data, 0, 0, yy, xx);
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = 0;
        }
    }
}

void FilterHelper::extract(
            const Volume<RFLOAT>& src,
            Volume<RFLOAT>& dest,
            long int x0, long int y0, long int z0,
            long int w, long int h, long int d)
{
    dest.resize(w,h,d);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        long int xx = x0 + x;
        long int yy = y0 + y;
        long int zz = z0 + z;

        if (   xx >= 0 && xx < src.dimx
            && yy >= 0 && yy < src.dimy
            && zz >= 0 && zz < src.dimz)
        {
            dest(x, y, z) = src(xx, yy, zz);
        }
    }
}

void FilterHelper::signedDist(const Image<RFLOAT>& src, Image<RFLOAT>& dest)
{
    dest = Image<RFLOAT>(src.data.xdim, src.data.ydim, src.data.zdim);

    Image<RFLOAT>
            ggp(src.data.xdim, src.data.ydim, src.data.zdim),
            ggn(src.data.xdim, src.data.ydim, src.data.zdim),
            gp(src.data.xdim, src.data.ydim, src.data.zdim),
            gn(src.data.xdim, src.data.ydim, src.data.zdim),
            hp(src.data.xdim, src.data.ydim, src.data.zdim),
            hn(src.data.xdim, src.data.ydim, src.data.zdim),
            s(src.data.xdim, src.data.ydim, src.data.zdim);

    double rmax2 = 4.0 * (src.data.xdim*src.data.xdim + src.data.ydim*src.data.ydim + src.data.zdim*src.data.zdim);

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    {
        DIRECT_A3D_ELEM(ggp.data, z, y, 0) = rmax2;
        DIRECT_A3D_ELEM(ggn.data, z, y, 0) = rmax2;

        for (long int x = 1; x < dest.data.xdim; x++)
        {
            if (DIRECT_A3D_ELEM(src.data, z, y, x) < 0.0)
            {
                DIRECT_A3D_ELEM(ggp.data, z, y, x) = 0;

                double d = sqrt(DIRECT_A3D_ELEM(ggn.data, z, y, x-1)) + 1.0;
                DIRECT_A3D_ELEM(ggn.data, z, y, x) = d*d;
            }
            else
            {
                DIRECT_A3D_ELEM(ggn.data, z, y, x) = 0;

                double d = sqrt(DIRECT_A3D_ELEM(ggp.data, z, y, x-1)) + 1.0;
                DIRECT_A3D_ELEM(ggp.data, z, y, x) = d*d;
            }
        }

        DIRECT_A3D_ELEM(gp.data, z, y, dest.data.xdim-1) = DIRECT_A3D_ELEM(ggp.data, z, y, dest.data.xdim-1);
        DIRECT_A3D_ELEM(gn.data, z, y, dest.data.xdim-1) = DIRECT_A3D_ELEM(ggn.data, z, y, dest.data.xdim-1);

        for (long int x = dest.data.xdim-2; x >= 0; x--)
        {
            double dp = sqrt(DIRECT_A3D_ELEM(gp.data, z, y, x+1)) + 1.0;
            double ddp = dp*dp;

            double dn = sqrt(DIRECT_A3D_ELEM(gn.data, z, y, x+1)) + 1.0;
            double ddn = dn*dn;

            if (ddp < DIRECT_A3D_ELEM(ggp.data, z, y, x))
            {
                DIRECT_A3D_ELEM(gp.data, z, y, x) = ddp;
            }
            else
            {
                DIRECT_A3D_ELEM(gp.data, z, y, x) = DIRECT_A3D_ELEM(ggp.data, z, y, x);
            }

            if (ddn < DIRECT_A3D_ELEM(ggn.data, z, y, x))
            {
                DIRECT_A3D_ELEM(gn.data, z, y, x) = ddn;
            }
            else
            {
                DIRECT_A3D_ELEM(gn.data, z, y, x) = DIRECT_A3D_ELEM(ggn.data, z, y, x);
            }
        }
    }

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        long int rp = (long int) sqrt(DIRECT_A3D_ELEM(gp.data, z, y, x));
        long int rn = (long int) sqrt(DIRECT_A3D_ELEM(gn.data, z, y, x));

        double minValP = rmax2;
        double minValN = rmax2;

        for (long int yy = y-rp; yy <= y+rp; yy++)
        {
            if (yy < 0 || yy >= dest.data.ydim) continue;

            double dy = yy - y;
            double vgp = DIRECT_A3D_ELEM(gp.data, z, yy, x) + dy*dy;

            if (vgp < minValP) minValP = vgp;
        }

        for (long int yy = y-rn; yy <= y+rn; yy++)
        {
            if (yy < 0 || yy >= dest.data.ydim) continue;

            double dy = yy - y;
            double vgn = DIRECT_A3D_ELEM(gn.data, z, yy, x) + dy*dy;

            if (vgn < minValN) minValN = vgn;
        }

        DIRECT_A3D_ELEM(hp.data, z, y, x) = minValP;
        DIRECT_A3D_ELEM(hn.data, z, y, x) = minValN;
    }

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++)
    {
        if (DIRECT_A3D_ELEM(src.data, z, y, x) < 0.0)
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = -sqrt(DIRECT_A3D_ELEM(hn.data, z, y, x));
        }
        else
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = sqrt(DIRECT_A3D_ELEM(hp.data, z, y, x));
        }
    }
}

void FilterHelper::erode3x3(Image<RFLOAT>& src, Image<RFLOAT>& dest)
{
    for (long int z = 0; z < src.data.zdim; z++)
    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        double v = std::numeric_limits<double>::max();

        for (long int zz = z-1; zz <= z+1; zz++)
        for (long int yy = y-1; yy <= y+1; yy++)
        for (long int xx = x-1; xx <= x+1; xx++)
        {
            if (   xx >= 0 && xx < src.data.xdim
                && yy >= 0 && yy < src.data.ydim
                && zz >= 0 && zz < src.data.zdim
                && DIRECT_A3D_ELEM(src.data, zz, yy, xx) < v)
            {
                v = DIRECT_A3D_ELEM(src.data, zz, yy, xx);
            }
        }

        DIRECT_A3D_ELEM(dest.data, z, y, x) = v;
    }
}

void FilterHelper::localMinima(Image<RFLOAT>& src, Image<RFLOAT>& dest, RFLOAT thresh)
{
    dest = Image<RFLOAT>(src.data.xdim, src.data.ydim, src.data.zdim);

    for (long int z = 0; z < src.data.zdim; z++)
    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        if (DIRECT_A3D_ELEM(src.data, z, y, x) > thresh)
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = 0.f;
            continue;
        }

        double v = std::numeric_limits<double>::max();

        for (long int zz = z-1; zz <= z+1; zz++)
        for (long int yy = y-1; yy <= y+1; yy++)
        for (long int xx = x-1; xx <= x+1; xx++)
        {
            if (   xx >= 0 && xx < src.data.xdim
                && yy >= 0 && yy < src.data.ydim
                && zz >= 0 && zz < src.data.zdim
                && DIRECT_A3D_ELEM(src.data, zz, yy, xx) < v)
            {
                v = DIRECT_A3D_ELEM(src.data, zz, yy, xx);
            }
        }

        if (v == DIRECT_A3D_ELEM(src.data, z, y, x))
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = 1.f;
        }
        else
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = 0.f;
        }
    }
}

std::vector<gravis::d3Vector> FilterHelper::localMinima(Image<RFLOAT>& src, RFLOAT thresh)
{
    std::vector<gravis::d3Vector> out(0);

    for (long int z = 0; z < src.data.zdim; z++)
    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++)
    {
        if (DIRECT_A3D_ELEM(src.data, z, y, x) > thresh)
        {
            continue;
        }

        double v = std::numeric_limits<double>::max();

        for (long int zz = z-1; zz <= z+1; zz++)
        for (long int yy = y-1; yy <= y+1; yy++)
        for (long int xx = x-1; xx <= x+1; xx++)
        {
            if (   xx >= 0 && xx < src.data.xdim
                && yy >= 0 && yy < src.data.ydim
                && zz >= 0 && zz < src.data.zdim
                && DIRECT_A3D_ELEM(src.data, zz, yy, xx) < v)
            {
                v = DIRECT_A3D_ELEM(src.data, zz, yy, xx);
            }
        }

        if (v == DIRECT_A3D_ELEM(src.data, z, y, x))
        {
            out.push_back(d3Vector(x,y,z));
        }
    }

    return out;
}

void FilterHelper::centralGradient(const Volume<RFLOAT>& src, Volume<t3Vector<RFLOAT> >& dest)
{
    const size_t dimx = src.dimx;
    const size_t dimy = src.dimy;
    const size_t dimz = src.dimz;

    dest.resize(dimx, dimy, dimz);

    FOR_ALL_VOXELS(src)
    {
        if (dimx == 0)
        {
            dest(x,y,z).x = 0;
        }
        else if (x == 0)
        {
            dest(x,y,z).x = src(x+1,y,z) - src(x,y,z);
        }
        else if (x < dimx - 1)
        {
            dest(x,y,z).x = 0.5 * (src(x+1,y,z) - src(x-1,y,z));
        }
        else
        {
            dest(x,y,z).x = src(x,y,z) - src(x-1,y,z);
        }

        if (dimy == 0)
        {
            dest(x,y,z).y = 0;
        }
        else if (y == 0)
        {
            dest(x,y,z).y = src(x,y+1,z) - src(x,y,z);
        }
        else if (y < dimy - 1)
        {
            dest(x,y,z).y = 0.5 * (src(x,y+1,z) - src(x,y-1,z));
        }
        else
        {
            dest(x,y,z).y = src(x,y,z) - src(x,y-1,z);
        }

        if (dimz == 0)
        {
            dest(x,y,z).z = 0;
        }
        else if (z == 0)
        {
            dest(x,y,z).z = src(x,y,z+1) - src(x,y,z);
        }
        else if (z < dimz - 1)
        {
            dest(x,y,z).z = 0.5 * (src(x,y,z+1) - src(x,y,z-1));
        }
        else
        {
            dest(x,y,z).z = src(x,y,z) - src(x,y,z-1);
        }
    }
}

t3Vector<RFLOAT> FilterHelper::centralGradient(const Volume<RFLOAT>& src, size_t x, size_t y, size_t z)
{
    t3Vector<RFLOAT> out;

    if (src.dimx == 0)
    {
        out.x = 0;
    }
    else if (x == 0)
    {
        out.x = src(x+1,y,z) - src(x,y,z);
    }
    else if (x < src.dimx - 1)
    {
        out.x = 0.5 * (src(x+1,y,z) - src(x-1,y,z));
    }
    else
    {
        out.x = src(x,y,z) - src(x-1,y,z);
    }

    if (src.dimy == 0)
    {
        out.y = 0;
    }
    else if (y == 0)
    {
        out.y = src(x,y+1,z) - src(x,y,z);
    }
    else if (y < src.dimy - 1)
    {
        out.y = 0.5 * (src(x,y+1,z) - src(x,y-1,z));
    }
    else
    {
        out.y = src(x,y,z) - src(x,y-1,z);
    }

    if (src.dimz == 0)
    {
        out.z = 0;
    }
    else if (z == 0)
    {
        out.z = src(x,y,z+1) - src(x,y,z);
    }
    else if (z < src.dimz - 1)
    {
        out.z = 0.5 * (src(x,y,z+1) - src(x,y,z-1));
    }
    else
    {
        out.z = src(x,y,z) - src(x,y,z-1);
    }

	return out;
}

MultidimArray<Complex> FilterHelper::FriedelExpand(const MultidimArray<Complex> &half)
{
	const int wh = half.xdim;
	const int h = half.ydim;
	const int d = half.zdim;
	const int c = half.ndim;
	
	const int w = 2*(wh-1);
	
	MultidimArray<Complex> out(d,h,w);
	
	for (int n = 0; n < c; n++)
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	{
		const int zz = (d - z) % d;
		const int yy = (h - y) % h;
		
		for (int x = 0; x < wh; x++)
		{
			DIRECT_NZYX_ELEM(out, n, z, y, x) = DIRECT_NZYX_ELEM(half, n, z, y, x);
		}
		
		for (int x = wh; x < w; x++)
		{
			DIRECT_NZYX_ELEM(out, n, z, y, x) = DIRECT_NZYX_ELEM(half, n, zz, yy, w-x).conj();
		}
	}	
	
	return out;
}

Image<RFLOAT> FilterHelper::normaliseToUnitInterval(const Image<RFLOAT> &img)
{
	const int w = img.data.xdim;
	const int h = img.data.ydim;
	const int d = img.data.zdim;
	const int c = img.data.ndim;
	
	RFLOAT minVal = std::numeric_limits<RFLOAT>::max();
	RFLOAT maxVal = -std::numeric_limits<RFLOAT>::max();
			
	for (int n = 0; n < c; n++)
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		RFLOAT v = DIRECT_NZYX_ELEM(img.data, n, z, y, x);
		
		if (v > maxVal) maxVal = v;
		if (v < minVal) minVal = v;
	}
	
	Image<RFLOAT> out(w,h,d,c);
	
	for (int n = 0; n < c; n++)
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		RFLOAT v = DIRECT_NZYX_ELEM(img.data, n, z, y, x);
		DIRECT_NZYX_ELEM(out.data, n, z, y, x) = (v - minVal)/(maxVal - minVal);
	}
	
	return out;
}

Image<RFLOAT> FilterHelper::normaliseToUnitIntervalSigned(const Image<RFLOAT> &img)
{
	const int w = img.data.xdim;
	const int h = img.data.ydim;
	const int d = img.data.zdim;
	const int c = img.data.ndim;
	
	RFLOAT maxAbs = 0;
			
	for (int n = 0; n < c; n++)
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		RFLOAT v = std::abs(DIRECT_NZYX_ELEM(img.data, n, z, y, x));
		
		if (v > maxAbs) maxAbs = v;
	}
	
	Image<RFLOAT> out(w,h,d,c);
	
	for (int n = 0; n < c; n++)
	for (int z = 0; z < d; z++)
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		RFLOAT v = DIRECT_NZYX_ELEM(img.data, n, z, y, x);
		DIRECT_NZYX_ELEM(out.data, n, z, y, x) = v / maxAbs;
	}
	
	return out;
}

void FilterHelper::uniqueInfluenceMask(std::vector<gravis::d3Vector> pts, Image<RFLOAT>& dest, Image<RFLOAT>& indexDest, RFLOAT thresh)
{
    const long int w = dest.data.xdim;
    const long int h = dest.data.ydim;
    const long int pc = pts.size();

    indexDest = Image<RFLOAT>(w,h);

    const double t2 = thresh * thresh;

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        int closer = 0;
        int lastIndex = -1;

        for (long int p = 0; p < pc; p++)
        {
            d2Vector d(x - pts[p].x, y - pts[p].y);

            if (d.norm2() < t2)
            {
                closer++;
                lastIndex = p;
            }
        }

        if (closer == 1)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = 1.0;
            DIRECT_NZYX_ELEM(indexDest.data, 0, 0, y, x) = lastIndex;
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = 0.0;
            DIRECT_NZYX_ELEM(indexDest.data, 0, 0, y, x) = -1.0;
        }
    }
}

void FilterHelper::polarRemap(d2Vector pos, const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                   const Image<RFLOAT>& mask, Image<RFLOAT>& maskDest, int phiRes, int rRes, double rMax)
{
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;

    dest = Image<RFLOAT>(phiRes, rRes, 1);
    maskDest = Image<RFLOAT>(phiRes, rRes, 1);

    for (long int ri = 0; ri < rRes; ri++)
    for (long int p = 0; p < phiRes; p++)
    {
        const double r = rMax * ri / (double)rRes;
        const double phi = 2.0 * PI * p / (double)phiRes;

        d2Vector pp = pos + r * d2Vector(cos(phi),sin(phi));

        int ppnnx = (int)(pp.x + 0.5);
        int ppnny = (int)(pp.y + 0.5);

        if (   ppnnx > 0 && ppnnx < w - 1
            && ppnny > 0 && ppnny < h - 1
            && DIRECT_NZYX_ELEM(mask.data, 0, 0, ppnny, ppnnx) > 0.5)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, 0, ri, p) = Interpolation::linearXY(src, pp.x, pp.y, 0);
            DIRECT_NZYX_ELEM(maskDest.data, 0, 0, ri, p) = 1.0;
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, 0, 0, ri, p) = 0.0;
            DIRECT_NZYX_ELEM(maskDest.data, 0, 0, ri, p) = 0.0;
        }
    }
}



void FilterHelper::polarRemap(d2Vector pos, const Image<RFLOAT>& distTransf, const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                   const Image<RFLOAT>& mask, Image<RFLOAT>& maskDest, int phiRes, int rRes, double rMax)
{
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;

    dest = Image<RFLOAT>(phiRes, rRes, 1);
    maskDest = Image<RFLOAT>(phiRes, rRes, 1);

    for (long int r = 0; r < rRes; r++)
    for (long int p = 0; p < phiRes; p++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, 0, r, p) = 0.0;
        DIRECT_NZYX_ELEM(maskDest.data, 0, 0, r, p) = 0.0;
    }

    const int x0 = (int)(pos.x - rMax + 0.5);
    const int x1 = (int)(pos.x + rMax + 0.5);
    const int y0 = (int)(pos.y - rMax + 0.5);
    const int y1 = (int)(pos.y + rMax + 0.5);

    for (int y = y0; y <= y1; y++)
    for (int x = x0; x <= x1; x++)
    {
        const double dx = x - pos.x;
        const double dy = y - pos.y;

        if (x < 1 || x >= w-1 || y < 1 || y >= h-1 || (dx == 0.0 && dy == 0.0))
        {
            continue;
        }

        double phiR = std::atan2(dy,dx);

        if (phiR < 0.0) phiR += 2.0*PI;

        const double phiD = phiRes * phiR / (2.0*PI);
        const int phi0 = ((int)(phiD)) % phiRes;
        const int phi1 = ((int)(phiD)+1) % phiRes;
        const double phiF = phiD - (double)phi0;

        const double rD = rRes * DIRECT_NZYX_ELEM(distTransf.data, 0, 0, y, x) / rMax;
        const int r0 = (int)rD;
        const int r1 = (int)rD + 1;
        const double rF = rD - r0;

        const double v = DIRECT_NZYX_ELEM(src.data, 0, 0, y, x);

        if (r0 >= 0 && r0 < rRes)
        {
            DIRECT_NZYX_ELEM(dest.data,     0, 0, r0, phi0) += (1.0 - rF) * (1.0 - phiF) * v;
            DIRECT_NZYX_ELEM(maskDest.data, 0, 0, r0, phi0) += (1.0 - rF) * (1.0 - phiF);

            DIRECT_NZYX_ELEM(dest.data,     0, 0, r0, phi1) += (1.0 - rF) * phiF * v;
            DIRECT_NZYX_ELEM(maskDest.data, 0, 0, r0, phi1) += (1.0 - rF) * phiF;
        }

        if (r1 >= 0 && r1 < rRes)
        {
            DIRECT_NZYX_ELEM(dest.data,     0, 0, r1, phi0) += rF * (1.0 - phiF) * v;
            DIRECT_NZYX_ELEM(maskDest.data, 0, 0, r1, phi0) += rF * (1.0 - phiF);

            DIRECT_NZYX_ELEM(dest.data,     0, 0, r1, phi1) += rF * phiF * v;
            DIRECT_NZYX_ELEM(maskDest.data, 0, 0, r1, phi1) += rF * phiF;
        }

    }

    for (long int r = 0; r < rRes; r++)
    for (long int p = 0; p < phiRes; p++)
    {
        if (DIRECT_NZYX_ELEM(maskDest.data, 0, 0, r, p) > 0.0)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, 0, r, p) /= DIRECT_NZYX_ELEM(maskDest.data, 0, 0, r, p);
        }
    }
}

Image<RFLOAT> FilterHelper::cartToPolar(const Image<RFLOAT> &img)
{
    const int w0 = img.data.xdim;
    const int h0 = img.data.ydim;

    const double w0h = w0/2.0;
    const double h0h = h0/2.0;

    const double cx = w0h + 1;
    const double cy = h0h + 0.5;

    const int w = (int)(2.0*PI*w0h + 1);
    const int h = w0h;

    Image<RFLOAT> out(w,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        const double phi = 2.0 * PI * x / (double)w;
        const double r = w0h * y / (double)h;

        double xx = cx + r * cos(phi);
        double yy = cy + r * sin(phi);

        out(y,x) = Interpolation::cubicXY(img, xx, yy, 0, 0);
    }

    return out;
}

Image<RFLOAT> FilterHelper::polarToCart(const Image<RFLOAT> &img)
{
    const int wp = img.data.xdim;
    const int hp = img.data.ydim;

    const double w0h = hp;

    const double cx = w0h + 1;
    const double cy = w0h + 0.5;

    const int w = 2.0*w0h;
    const int h = 2.0*w0h;

    Image<RFLOAT> out(w,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        const double xd = x - cx;
        const double yd = y - cy;

        const double r = sqrt(xd*xd + yd*yd);
        double phi = (xd == 0 && yd == 0)? 0.0 : atan2(yd,xd);
        if (phi < 0.0) phi += 2.0*PI;

        out(y,x) = Interpolation::cubicXY(img, wp*phi/(2.0*PI), r, 0, 0);
    }

    return out;
}

Image<RFLOAT> FilterHelper::polarBlur(const Image<RFLOAT> &img, double sigma)
{
    Image<RFLOAT> img1 = FilterHelper::cartToPolar(img);
    Image<RFLOAT> img2 = img1;
    separableGaussianX_wrap(img1, img2, sigma);

    return FilterHelper::polarToCart(img2);
}

Image<RFLOAT> FilterHelper::sectorBlend(const Image<RFLOAT>& img0, const Image<RFLOAT>& img1, int sectors)
{
    const int w = img0.data.xdim;
    const int h = img0.data.ydim;

    if (img1.data.xdim != w || img1.data.ydim != h)
    {
        std::cerr << "FilterHelper::sectorBlend: unequal image size: " << w << "x" << h
                  << " vs. " << img1.data.xdim << "x" << img1.data.ydim << "\n";

        REPORT_ERROR("FilterHelper::sectorBlend: unequal image size.");
    }

    Image<RFLOAT> out(w,h);

    const double cx = w/2.0;
    const double cy = h/2.0;

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    {
        const double xd = x - cx;
        const double yd = y - cy;

        double phi = (xd == 0 && yd == 0)? 0.0 : atan2(yd,xd) + PI;
        double a = sectors*phi/(2.0*PI);

        out(y,x) = a - (int)a < 0.5? img0(y,x) : img1(y,x);
    }

    return out;
}

void FilterHelper::diffuseAlongIsocontours2D(const Image<RFLOAT>& src, const Image<RFLOAT>& guide,
                                              Image<RFLOAT>& dest, int iters, RFLOAT sigma, RFLOAT lambda, RFLOAT delta)
{
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;

    const bool sobel = true;

    dest = Image<RFLOAT>(w, h, 1);

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = DIRECT_NZYX_ELEM(src.data, 0, 0, y, x);
    }

    Volume<d2Vector> flux(w,h,1);
    flux.fill(d2Vector(0,0));

    Volume<Tensor2x2<RFLOAT> > D0(w,h,1), D(w,h,1), J(w,h,1);
    D0.fill(Tensor2x2<RFLOAT>(0.0));
    D.fill(Tensor2x2<RFLOAT>(0.0));
    J.fill(Tensor2x2<RFLOAT>(0.0));

    for (long int y = 1; y < h-1; y++)
    for (long int x = 1; x < w-1; x++)
    {
        d2Vector g;

        if (sobel)
        {
            double gxp =   0.25 * DIRECT_NZYX_ELEM(guide.data, 0, 0, y-1, x+1)
                         + 0.5  * DIRECT_NZYX_ELEM(guide.data, 0, 0, y,   x+1)
                         + 0.25 * DIRECT_NZYX_ELEM(guide.data, 0, 0, y+1, x+1);

            double gxn =   0.25 * DIRECT_NZYX_ELEM(guide.data, 0, 0, y-1, x-1)
                         + 0.5  * DIRECT_NZYX_ELEM(guide.data, 0, 0, y,   x-1)
                         + 0.25 * DIRECT_NZYX_ELEM(guide.data, 0, 0, y+1, x-1);

            double gyp =   0.25 * DIRECT_NZYX_ELEM(guide.data, 0, 0, y+1, x-1)
                         + 0.5  * DIRECT_NZYX_ELEM(guide.data, 0, 0, y+1, x)
                         + 0.25 * DIRECT_NZYX_ELEM(guide.data, 0, 0, y+1, x+1);

            double gyn =   0.25 * DIRECT_NZYX_ELEM(guide.data, 0, 0, y-1, x-1)
                         + 0.5  * DIRECT_NZYX_ELEM(guide.data, 0, 0, y-1, x)
                         + 0.25 * DIRECT_NZYX_ELEM(guide.data, 0, 0, y-1, x+1);

            g.x = 0.5 * (gxp - gxn);
            g.y = 0.5 * (gyp - gyn);
        }
        else
        {
            g.x = 0.5 * (DIRECT_NZYX_ELEM(guide.data, 0, 0, y, x+1) - DIRECT_NZYX_ELEM(guide.data, 0, 0, y, x-1));
            g.y = 0.5 * (DIRECT_NZYX_ELEM(guide.data, 0, 0, y+1, x) - DIRECT_NZYX_ELEM(guide.data, 0, 0, y-1, x));
        }

        D0(x,y,0) = Tensor2x2<RFLOAT>::autoDyadicProduct(t2Vector<RFLOAT>(g.x, g.y));
    }

    separableGaussian(D0, D, sigma);

    //Volume<RFLOAT> dbg0(w,h,1), dbg1(w,h,1);

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        t2Matrix<RFLOAT> DxyR = D(x,y,0).toMatrix();
        d2Matrix Dxy(DxyR(0,0), DxyR(0,1), DxyR(1,0), DxyR(1,1));

        double qx, qy, l0, l1;
        dsyev2(Dxy(0,0), Dxy(0,1), Dxy(1,1), &l0, &l1, &qx, &qy);


        double dl = l0 - l1;
        RFLOAT ani = 1.0 - exp(-0.5*dl*dl/(lambda*lambda));

        d2Vector f(-qy, qx);

        //dbg0(x,y,0) = f.length();

        J(x,y,0) = ani * Tensor2x2<RFLOAT>::autoDyadicProduct(t2Vector<RFLOAT>(f.x, f.y));
    }

    //VtkHelper::writeVTK(dbg0, "f_len.vtk");

    for (int it = 0; it < iters; it++)
    {

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 1; y < h-1; y++)
        for (long int x = 1; x < w-1; x++)
        {
            d2Vector g;

            if (sobel)
            {
                double gxp =   0.25 * DIRECT_NZYX_ELEM(dest.data, 0, 0, y-1, x+1)
                             + 0.5  * DIRECT_NZYX_ELEM(dest.data, 0, 0, y,   x+1)
                             + 0.25 * DIRECT_NZYX_ELEM(dest.data, 0, 0, y+1, x+1);

                double gxn =   0.25 * DIRECT_NZYX_ELEM(dest.data, 0, 0, y-1, x-1)
                             + 0.5  * DIRECT_NZYX_ELEM(dest.data, 0, 0, y,   x-1)
                             + 0.25 * DIRECT_NZYX_ELEM(dest.data, 0, 0, y+1, x-1);

                double gyp =   0.25 * DIRECT_NZYX_ELEM(dest.data, 0, 0, y+1, x-1)
                             + 0.5  * DIRECT_NZYX_ELEM(dest.data, 0, 0, y+1,   x)
                             + 0.25 * DIRECT_NZYX_ELEM(dest.data, 0, 0, y+1, x+1);

                double gyn =   0.25 * DIRECT_NZYX_ELEM(dest.data, 0, 0, y-1, x-1)
                             + 0.5  * DIRECT_NZYX_ELEM(dest.data, 0, 0, y-1,   x)
                             + 0.25 * DIRECT_NZYX_ELEM(dest.data, 0, 0, y-1, x+1);

                g.x = 0.5 * (gxp - gxn);
                g.y = 0.5 * (gyp - gyn);
            }
            else
            {
                g.x = 0.5 * (DIRECT_NZYX_ELEM(guide.data, 0, 0, y, x+1) - DIRECT_NZYX_ELEM(guide.data, 0, 0, y, x-1));
                g.y = 0.5 * (DIRECT_NZYX_ELEM(guide.data, 0, 0, y+1, x) - DIRECT_NZYX_ELEM(guide.data, 0, 0, y-1, x));
            }

            t2Vector<RFLOAT> fR = J(x,y,0).toMatrix() * t2Vector<RFLOAT>(g.x, g.y);
            flux(x,y,0) = d2Vector(fR.x, fR.y);
        }

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 1; y < h-1; y++)
        for (long int x = 1; x < w-1; x++)
        {
            double div = 0.0;

            div += flux(x+1,y,0).x - flux(x-1,y,0).x;
            div += flux(x,y+1,0).y - flux(x,y-1,0).y;

            DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) +=  delta * div;
        }
    }
}

void FilterHelper::EED_2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest, int iters, double sigma, double delta, double tau)
{
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = DIRECT_NZYX_ELEM(src.data, 0, 0, y, x);
    }

    Image<RFLOAT> smooth;
    separableGaussianXY(dest, smooth, sigma);

    Volume<d2Vector> flux(w,h,1);
    flux.fill(d2Vector(0,0));

    double tt = tau*tau;

    for (int it = 0; it < iters; it++)
    {

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 1; y < h-1; y++)
        for (long int x = 1; x < w-1; x++)
        {
            d2Vector g, gs;

            g.x = DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x+1) - DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x);
            g.y = DIRECT_NZYX_ELEM(dest.data, 0, 0, y+1, x) - DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x);
            gs.x = DIRECT_NZYX_ELEM(smooth.data, 0, 0, y, x+1) - DIRECT_NZYX_ELEM(smooth.data, 0, 0, y, x);
            gs.y = DIRECT_NZYX_ELEM(smooth.data, 0, 0, y+1, x) - DIRECT_NZYX_ELEM(smooth.data, 0, 0, y, x);

            double iso = exp(-0.5*gs.norm2()/tt);

            double gsl = gs.length();
            if (gsl > 0.0) gs /= gsl;

            d2Vector gn = g.dot(gs) * gs;
            d2Vector gp = g - gn;

            flux(x,y,0) = iso * g + (1.0 - iso) * gp;
        }


        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 1; y < h-1; y++)
        for (long int x = 1; x < w-1; x++)
        {
            double div = 0.0;

            div += flux(x,y,0).x - flux(x-1,y,0).x;
            div += flux(x,y,0).y - flux(x,y-1,0).y;

            DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) +=  delta * div;
        }
    }
}

void FilterHelper::descendTV(const Image<RFLOAT>& src, Image<RFLOAT>& dest, double delta)
{
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;
    const long int d = src.data.zdim;

    std::vector<RFLOAT> vals;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        const double v0 = DIRECT_NZYX_ELEM(src.data, 0, z, y, x);

        vals.clear();
        vals.reserve(6);

        if (x > 0)   vals.push_back(DIRECT_NZYX_ELEM(src.data, 0, z, y, x-1));
        if (x < w-1) vals.push_back(DIRECT_NZYX_ELEM(src.data, 0, z, y, x+1));
        if (y > 0)   vals.push_back(DIRECT_NZYX_ELEM(src.data, 0, z, y-1, x));
        if (y < h-1) vals.push_back(DIRECT_NZYX_ELEM(src.data, 0, z, y+1, x));
        if (z > 0)   vals.push_back(DIRECT_NZYX_ELEM(src.data, 0, z-1, y, x));
        if (z < d-1) vals.push_back(DIRECT_NZYX_ELEM(src.data, 0, z+1, y, x));

        std::vector<int> order = IndexSort<RFLOAT>::sortIndices(vals);
        const int c = vals.size();

        double vm;

        if (vals.size() % 2 == 0)
        {
            vm = 0.5 * (vals[order[c/2]] + vals[order[c/2 - 1]]);
        }
        else
        {
            vm = vals[order[c/2]];
        }

        if (std::abs(v0 - vm) < delta)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = vm;
        }
        else if (v0 < vm)
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = v0 + delta;
        }
        else
        {
            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = v0 - delta;
        }
    }
}

void FilterHelper::descendTV2(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                              Volume<gravis::d3Vector>& xi, Volume<RFLOAT>& uBar,
                              int iters, double sigma, double tau)
{
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;
    const long int d = src.data.zdim;

    fwdGrad(src,xi);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        uBar(x,y,z) = DIRECT_NZYX_ELEM(src.data, 0, z, y, x);
        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = DIRECT_NZYX_ELEM(src.data, 0, z, y, x);
    }

    for (int it = 0; it < iters; it++)
    {
        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++)
        {
            d3Vector gradUBar;

            if (w == 1)
            {
                gradUBar.x = 0;
            }
            else if (x < w - 1)
            {
                gradUBar.x = uBar(x+1,y,z) - uBar(x,y,z);
            }
            else
            {
                gradUBar.x = uBar(x,y,z) - uBar(x-1,y,z);
            }

            if (h == 1)
            {
                gradUBar.y = 0;
            }
            else if (y < h - 1)
            {
                gradUBar.y = uBar(x,y+1,z) - uBar(x,y,z);
            }
            else
            {
                gradUBar.y = uBar(x,y,z) - uBar(x,y-1,z);
            }

            if (d == 1)
            {
                gradUBar.z = 0;
            }
            else if (z < d - 1)
            {
                gradUBar.z = uBar(x,y,z+1) - uBar(x,y,z);
            }
            else
            {
                gradUBar.z = uBar(x,y,z) - uBar(x,y,z-1);
            }

            d3Vector nextXi = xi(x,y,z) + sigma * gradUBar;

            double nxl = nextXi.length();

            xi(x,y,z) = nxl > 0.0? nextXi/nxl : d3Vector(0,0,0);
        }

        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++)
        {
            double divXi = 0.0;

            if (x > 0)
            {
                divXi += xi(x,y,z).x - xi(x-1,y,z).x;
            }

            if (y > 0)
            {
                divXi += xi(x,y,z).y - xi(x,y-1,z).y;
            }

            if (z > 0)
            {
                divXi += xi(x,y,z).z - xi(x,y,z-1).z;
            }

            double du = tau * divXi;

            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) += du;
            uBar(x,y,z) = DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) + 0.5*du;
        }
    }
}

void FilterHelper::segmentTV(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                               Volume<gravis::d3Vector>& xi, Volume<RFLOAT>& uBar,
                               int iters, double sigma, double tau, double nu)
{
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;
    const long int d = src.data.zdim;

    xi.fill(d3Vector(0.0,0.0,0.0));
    uBar.fill(0.0);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        uBar(x,y,z) = DIRECT_NZYX_ELEM(src.data, 0, z, y, x);
        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = 0.0;
    }

    for (int it = 0; it < iters; it++)
    {
        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++)
        {
            d3Vector gradUBar;

            if (w == 1)
            {
                gradUBar.x = 0;
            }
            else if (x < w - 1)
            {
                gradUBar.x = uBar(x+1,y,z) - uBar(x,y,z);
            }
            else
            {
                gradUBar.x = uBar(x,y,z) - uBar(x-1,y,z);
            }

            if (h == 1)
            {
                gradUBar.y = 0;
            }
            else if (y < h - 1)
            {
                gradUBar.y = uBar(x,y+1,z) - uBar(x,y,z);
            }
            else
            {
                gradUBar.y = uBar(x,y,z) - uBar(x,y-1,z);
            }

            if (d == 1)
            {
                gradUBar.z = 0;
            }
            else if (z < d - 1)
            {
                gradUBar.z = uBar(x,y,z+1) - uBar(x,y,z);
            }
            else
            {
                gradUBar.z = uBar(x,y,z) - uBar(x,y,z-1);
            }

            d3Vector nextXi = xi(x,y,z) + sigma * gradUBar;

            double nxl = nextXi.length();

            xi(x,y,z) = nxl > 0.0? nextXi/nxl : d3Vector(0,0,0);
        }

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++)
        {
            double divXi = 0.0;

            if (x > 0)
            {
                divXi += xi(x,y,z).x - xi(x-1,y,z).x;
            }

            if (y > 0)
            {
                divXi += xi(x,y,z).y - xi(x,y-1,z).y;
            }

            if (z > 0)
            {
                divXi += xi(x,y,z).z - xi(x,y,z-1).z;
            }

            double u = DIRECT_NZYX_ELEM(dest.data, 0, z, y, x);
            double du = tau * (nu * divXi + DIRECT_NZYX_ELEM(src.data, 0, z, y, x));

            double nextU = u + du;

            if (nextU > 1.0) nextU = 1.0;
            else if (nextU < 0) nextU = 0;


            DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = nextU;
            uBar(x,y,z) = 2.0 * nextU - u;
        }
    }
}

void FilterHelper::segmentTVAniso2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                               Volume<gravis::d2Vector>& xi, Volume<RFLOAT>& uBar,
                               int iters, double sigma, double tau, double nu,
                               double rho, double theta, double alpha)
{
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;


    Image<RFLOAT> smooth;
    separableGaussianXY(src, smooth, rho);

    Volume<d2Vector> smoothGrad(w,h,1);
    fwdGrad2D(smooth, smoothGrad);

    xi.fill(d2Vector(0.0,0.0));
    uBar.fill(0.0);

    Volume<d2Matrix> D(w,h,1);

    const double tt = theta * theta;

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        uBar(x,y,0) = DIRECT_NZYX_ELEM(src.data, 0, 0, y, x);
        DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = 0.0;

        d2Vector gs = smoothGrad(x,y,0);
        double iso = exp(-0.5*gs.norm2()/tt);

        double gsl = gs.length();
        if (gsl > 0.0) gs /= gsl;

        d2Matrix I;

        //d2Matrix G = E - d2Matrix(gs.x*gs.x, gs.y*gs.x, gs.x*gs.y, gs.y*gs.y);
        // G x = x - (x dot gs) gs

        d2Matrix F = d2Matrix(gs.x*gs.x, gs.y*gs.x, gs.x*gs.y, gs.y*gs.y);

        d2Matrix G = sqrt(alpha) * F + sqrt((3.0 - alpha)/2.0)*(I - F);

        D(x,y,0) = sqrt(nu) * (iso * I + (1.0 - iso) * G);
    }


    for (int it = 0; it < iters; it++)
    {
        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++)
        {
            d2Vector gradUBar;

            if (w == 1)
            {
                gradUBar.x = 0;
            }
            else if (x < w - 1)
            {
                gradUBar.x = uBar(x+1,y,0) - uBar(x,y,0);
            }
            else
            {
                gradUBar.x = uBar(x,y,0) - uBar(x-1,y,0);
            }

            if (h == 1)
            {
                gradUBar.y = 0;
            }
            else if (y < h - 1)
            {
                gradUBar.y = uBar(x,y+1,0) - uBar(x,y,0);
            }
            else
            {
                gradUBar.y = uBar(x,y,0) - uBar(x,y-1,0);
            }

            d2Vector nextXi = D(x,y,0) * (xi(x,y,0) + sigma * D(x,y,0) * gradUBar);

            double nxl = nextXi.length();

            xi(x,y,0) = nxl > 0.0? (nextXi/nxl) : d2Vector(0,0);
        }

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++)
        {
            double divXi = 0.0;

            if (x > 0)
            {
                divXi += xi(x,y,0).x - xi(x-1,y,0).x;
            }

            if (y > 0)
            {
                divXi += xi(x,y,0).y - xi(x,y-1,0).y;
            }

            double u = DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x);
            double du = tau * (divXi + DIRECT_NZYX_ELEM(src.data, 0, 0, y, x));

            double nextU = u + du;

            if (nextU > 1.0) nextU = 1.0;
            else if (nextU < 0) nextU = 0;


            DIRECT_NZYX_ELEM(dest.data, 0, 0, y, x) = nextU;
            uBar(x,y,0) = 2.0 * nextU - u;
        }
    }
}


void FilterHelper::fwdGrad(const Image<RFLOAT>& u, Volume<gravis::d3Vector>& dest)
{
    const long int w = u.data.xdim;
    const long int h = u.data.ydim;
    const long int d = u.data.zdim;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        if (w == 1)
        {
            dest(x,y,z).x = 0;
        }
        else if (x < w - 1)
        {
            dest(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x+1) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x);
        }
        else
        {
            dest(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x-1);
        }

        if (h == 1)
        {
            dest(x,y,z).y = 0;
        }
        else if (y < h - 1)
        {
            dest(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y+1, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x);
        }
        else
        {
            dest(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y-1, x);
        }

        if (d == 1)
        {
            dest(x,y,z).z = 0;
        }
        else if (z < d - 1)
        {
            dest(x,y,z).z = DIRECT_NZYX_ELEM(u.data, 0, z+1, y, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x);
        }
        else
        {
            dest(x,y,z).z = DIRECT_NZYX_ELEM(u.data, 0, z, y, x) - DIRECT_NZYX_ELEM(u.data, 0, z-1, y, x);
        }
    }
}


void FilterHelper::fwdGrad2D(const Image<RFLOAT>& u, Volume<d2Vector>& dest)
{
    const long int w = u.data.xdim;
    const long int h = u.data.ydim;
    const long int d = 1;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        if (w == 1)
        {
            dest(x,y,z).x = 0;
        }
        else if (x < w - 1)
        {
            dest(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x+1) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x);
        }
        else
        {
            dest(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x-1);
        }

        if (h == 1)
        {
            dest(x,y,z).y = 0;
        }
        else if (y < h - 1)
        {
            dest(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y+1, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x);
        }
        else
        {
            dest(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y-1, x);
        }
    }
}

void FilterHelper::centralGrad2D(const Image<RFLOAT> &u, Volume<d2Vector> &dest)
{
    const long int w = u.data.xdim;
    const long int h = u.data.ydim;
    const long int d = 1;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        if (w == 1)
        {
            dest(x,y,z).x = 0;
        }
        else if (x < w-1 && x > 0)
        {
            dest(x,y,z).x = (DIRECT_NZYX_ELEM(u.data, 0, z, y, x+1) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x-1))/2.0;
        }
        else if (x == 0)
        {
            dest(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x+1) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x);
        }
        else if (x == w-1)
        {
            dest(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x-1);
        }

        if (h == 1)
        {
            dest(x,y,z).y = 0;
        }
        else if (y < h-1 && y > 0)
        {
            dest(x,y,z).y = (DIRECT_NZYX_ELEM(u.data, 0, z, y+1, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y-1, x))/2;
        }
        else if (y == 0)
        {
            dest(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y+1, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y, x);
        }
        else if (y == h-1)
        {
            dest(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y, x) - DIRECT_NZYX_ELEM(u.data, 0, z, y-1, x);
        }
    }

}

void FilterHelper::centralGrad2D(const Image<Complex> &u, Volume<d2Vector> &destRe, Volume<d2Vector> &destIm)
{
    const long int w = u.data.xdim;
    const long int h = u.data.ydim;
    const long int d = 1;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        if (w == 1)
        {
            destRe(x,y,z).x = 0;
            destIm(x,y,z).x = 0;
        }
        else if (x < w-1 && x > 0)
        {
            destRe(x,y,z).x = (DIRECT_NZYX_ELEM(u.data, 0, z, y, x+1).real - DIRECT_NZYX_ELEM(u.data, 0, z, y, x-1).real)/2.0;
            destIm(x,y,z).x = (DIRECT_NZYX_ELEM(u.data, 0, z, y, x+1).imag - DIRECT_NZYX_ELEM(u.data, 0, z, y, x-1).imag)/2.0;
        }
        else if (x == 0)
        {
            destRe(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x+1).real - DIRECT_NZYX_ELEM(u.data, 0, z, y, x).real;
            destIm(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x+1).imag - DIRECT_NZYX_ELEM(u.data, 0, z, y, x).imag;
        }
        else if (x == w-1)
        {
            destRe(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x).real - DIRECT_NZYX_ELEM(u.data, 0, z, y, x-1).real;
            destIm(x,y,z).x = DIRECT_NZYX_ELEM(u.data, 0, z, y, x).imag - DIRECT_NZYX_ELEM(u.data, 0, z, y, x-1).imag;
        }

        if (h == 1)
        {
            destRe(x,y,z).y = 0;
            destIm(x,y,z).y = 0;
        }
        else if (y < h-1 && y > 0)
        {
            destRe(x,y,z).y = (DIRECT_NZYX_ELEM(u.data, 0, z, y+1, x).real - DIRECT_NZYX_ELEM(u.data, 0, z, y-1, x).real)/2;
            destIm(x,y,z).y = (DIRECT_NZYX_ELEM(u.data, 0, z, y+1, x).imag - DIRECT_NZYX_ELEM(u.data, 0, z, y-1, x).imag)/2;
        }
        else if (y == 0)
        {
            destRe(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y+1, x).real - DIRECT_NZYX_ELEM(u.data, 0, z, y, x).real;
            destIm(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y+1, x).imag - DIRECT_NZYX_ELEM(u.data, 0, z, y, x).imag;
        }
        else if (y == h-1)
        {
            destRe(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y, x).real - DIRECT_NZYX_ELEM(u.data, 0, z, y-1, x).real;
            destIm(x,y,z).y = DIRECT_NZYX_ELEM(u.data, 0, z, y, x).imag - DIRECT_NZYX_ELEM(u.data, 0, z, y-1, x).imag;
        }
    }

}

void FilterHelper::blendSoft(const Image<Complex>& src0, const Image<Complex>& src1,
                             const Volume<RFLOAT>& mask, Image<Complex>& dest, RFLOAT bias1)
{
    const long int w = src0.data.xdim;
    const long int h = src0.data.ydim;
    const long int d = src0.data.zdim;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++)
    {
        const Complex v0 = DIRECT_NZYX_ELEM(src0.data, 0, z, y, x);
        const Complex v1 = DIRECT_NZYX_ELEM(src1.data, 0, z, y, x);
        const RFLOAT m = mask(x,y,z);

        DIRECT_NZYX_ELEM(dest.data, 0, z, y, x) = (v0 + bias1*m*v1)/(1.0 + bias1*m);
    }
}

double FilterHelper::totalVariation(const Image<RFLOAT>& src)
{
    double sum = 0.0;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(src.data)
    {
        if (i == src.data.ydim - 1 || j == src.data.xdim - 1) continue;

        double v0 = DIRECT_A2D_ELEM(src.data, i, j);
        double vx = DIRECT_A2D_ELEM(src.data, i, j+1);
        double vy = DIRECT_A2D_ELEM(src.data, i+1, j);

        double dx = vx - v0;
        double dy = vy - v0;

        double dtv = sqrt(dx*dx + dy*dy);

        sum += dtv;
    }

    return sum;
}

double FilterHelper::totalLogVariation(const Image<RFLOAT>& src, double delta)
{
    double sum = 0.0;

    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(src.data)
    {
        if (i == src.data.ydim - 1 || j == src.data.xdim - 1) continue;

        double v0 = DIRECT_A2D_ELEM(src.data, i, j);
        double vx = DIRECT_A2D_ELEM(src.data, i, j+1);
        double vy = DIRECT_A2D_ELEM(src.data, i+1, j);

        double dx = vx - v0;
        double dy = vy - v0;

        double dtv = log(delta + sqrt(dx*dx + dy*dy));

        sum += dtv;
    }

    return sum;
}

void FilterHelper::separableGaussianXYZ(const Image<RFLOAT>& src, Image<RFLOAT>& dest, RFLOAT sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.data.resize(src.data);

    std::vector<RFLOAT> kernel(2*k+1);
    const RFLOAT s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-0.5*i*i/s2);
    }

    Image<RFLOAT> temp(src.data.xdim, src.data.ydim, src.data.zdim);

    for (size_t z = 0; z < src.data.zdim; z++)
    for (size_t y = 0; y < src.data.ydim; y++)
    for (size_t x = 0; x < src.data.xdim; x++)
    {
        RFLOAT v = 0;
        RFLOAT m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int xx = x + i;
            if (xx < 0 || xx >= src.data.xdim) continue;

            v += kernel[i+k] * DIRECT_A3D_ELEM(src.data, z, y, xx);
            m += kernel[i+k];
        }

        DIRECT_A3D_ELEM(dest.data, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.data.zdim; z++)
    for (size_t y = 0; y < src.data.ydim; y++)
    for (size_t x = 0; x < src.data.xdim; x++)
    {
        RFLOAT v = 0;
        RFLOAT m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = y + i;
            if (yy < 0 || yy >= src.data.ydim) continue;

            v += kernel[i+k] * DIRECT_A3D_ELEM(dest.data, z, yy, x);
            m += kernel[i+k];
        }

        DIRECT_A3D_ELEM(temp.data, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.data.zdim; z++)
    for (size_t y = 0; y < src.data.ydim; y++)
    for (size_t x = 0; x < src.data.xdim; x++)
    {
        RFLOAT v = 0;
        RFLOAT m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int zz = z + i;
            if (zz < 0 || zz >= src.data.zdim) continue;

            v += kernel[i+k] * DIRECT_A3D_ELEM(temp.data, zz, y, x);
            m += kernel[i+k];
        }

        DIRECT_A3D_ELEM(dest.data, z, y, x) = v/m;
    }
}


void FilterHelper::separableGaussianXY(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                                       RFLOAT sigma, int k, bool wrap)
{
    if (!dest.data.sameShape(src.data))
    {
        dest.data.resize(src.data);
    }

    if (sigma <= 0.0)
    {
        for (size_t z = 0; z < src.data.zdim; z++)
        for (size_t y = 0; y < src.data.ydim; y++)
        for (size_t x = 0; x < src.data.xdim; x++)
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = DIRECT_A3D_ELEM(src.data, z, y, x);
        }

        return;
    }

    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    std::vector<RFLOAT> kernel(2*k+1);
    const RFLOAT s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-0.5*i*i/s2);
    }

    Image<RFLOAT> temp(src.data.xdim, src.data.ydim, src.data.zdim);

    for (size_t z = 0; z < src.data.zdim; z++)
    for (size_t y = 0; y < src.data.ydim; y++)
    for (size_t x = 0; x < src.data.xdim; x++)
    {
        RFLOAT v = 0;
        RFLOAT m = 0;

        for (long int i = -k; i <= k; i++)
        {
            long int xx = x + i;

            if (wrap) xx = (xx + src.data.xdim) % src.data.xdim;
            else if (xx < 0 || xx >= src.data.xdim) continue;

            v += kernel[i+k] * DIRECT_A3D_ELEM(src.data, z, y, xx);
            m += kernel[i+k];
        }

        DIRECT_A3D_ELEM(temp.data, z, y, x) = v/m;
    }

    for (size_t z = 0; z < src.data.zdim; z++)
    for (size_t y = 0; y < src.data.ydim; y++)
    for (size_t x = 0; x < src.data.xdim; x++)
    {
        RFLOAT v = 0;
        RFLOAT m = 0;

        for (long int i = -k; i <= k; i++)
        {
            long int yy = y + i;

            if (wrap) yy = (yy + src.data.ydim) % src.data.ydim;
            else if (yy < 0 || yy >= temp.data.ydim) continue;

            v += kernel[i+k] * DIRECT_A3D_ELEM(temp.data, z, yy, x);
            m += kernel[i+k];
        }

        DIRECT_A3D_ELEM(dest.data, z, y, x) = v/m;
    }
}

void FilterHelper::separableGaussianX_wrap(const Image<RFLOAT>& src, const Image<RFLOAT>& mask, Image<RFLOAT>& dest, RFLOAT sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.data.resize(src.data);

    std::vector<RFLOAT> kernel(2*k+1);
    const RFLOAT s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-0.5*i*i/s2);
    }

    for (size_t z = 0; z < src.data.zdim; z++)
    for (size_t y = 0; y < src.data.ydim; y++)
    for (size_t x = 0; x < src.data.xdim; x++)
    {
        RFLOAT v = 0;
        RFLOAT m = 0;

        for (long int i = -k; i <= k; i++)
        {
            long int xx = (x + i + src.data.xdim) % src.data.xdim;
            if (xx < 0 || xx >= src.data.xdim) continue;

            v += kernel[i+k] * DIRECT_A3D_ELEM(mask.data, z, y, xx) * DIRECT_A3D_ELEM(src.data, z, y, xx);
            m += kernel[i+k] * DIRECT_A3D_ELEM(mask.data, z, y, xx);
        }

        if (m > 0.0)
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = v/m;
        }
    }
}

void FilterHelper::separableGaussianX_wrap(const Image<RFLOAT>& src, Image<RFLOAT>& dest, RFLOAT sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.data.resize(src.data);

    std::vector<RFLOAT> kernel(2*k+1);
    const RFLOAT s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-0.5*i*i/s2);
    }

    for (size_t z = 0; z < src.data.zdim; z++)
    for (size_t y = 0; y < src.data.ydim; y++)
    for (size_t x = 0; x < src.data.xdim; x++)
    {
        RFLOAT v = 0;
        RFLOAT m = 0;

        for (long int i = -k; i <= k; i++)
        {
            long int xx = (x + i + src.data.xdim) % src.data.xdim;
            if (xx < 0 || xx >= src.data.xdim) continue;

            v += kernel[i+k] * DIRECT_A3D_ELEM(src.data, z, y, xx);
            m += kernel[i+k];
        }

        if (m > 0.0)
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = v/m;
        }
    }
}

void FilterHelper::averageX(const Image<RFLOAT>& src, const Image<RFLOAT>& mask, Image<RFLOAT>& dest)
{
    dest.data.resize(src.data);

    for (size_t z = 0; z < src.data.zdim; z++)
    for (size_t y = 0; y < src.data.ydim; y++)
    {
        RFLOAT v = 0;
        RFLOAT m = 0;

        for (size_t x = 0; x < src.data.xdim; x++)
        {
            v += DIRECT_A3D_ELEM(mask.data, z, y, x) * DIRECT_A3D_ELEM(src.data, z, y, x);
            m += DIRECT_A3D_ELEM(mask.data, z, y, x);
        }

        if (m > 0.0)
        {
            v /= m;
        }

        for (size_t x = 0; x < src.data.xdim; x++)
        {
            DIRECT_A3D_ELEM(dest.data, z, y, x) = v;
        }
    }
}

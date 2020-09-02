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

#include <src/jaz/single_particle/Fourier_helper.h>

void FourierHelper::FourierShift2D(MultidimArray<Complex>& img, RFLOAT xshift, RFLOAT yshift)
{
    const long w = img.xdim;
    const long h = img.ydim;

    xshift /= h;
    yshift /= h;

    if (ABS(xshift) < XMIPP_EQUAL_ACCURACY && ABS(yshift) < XMIPP_EQUAL_ACCURACY)
    {
        return;
    }

    for (long int yy = 0; yy < h; yy++)
    for (long int xx = 0; xx < w; xx++)
    {
        RFLOAT x = xx;
        RFLOAT y = yy < w? yy : yy - h;

        RFLOAT dotp = -2.0 * PI * (x * xshift + y * yshift);

        RFLOAT a, b;

        #ifdef RELION_SINGLE_PRECISION
            SINCOSF(dotp, &b, &a);
        #else
            SINCOS(dotp, &b, &a);
        #endif

        RFLOAT c = DIRECT_A2D_ELEM(img, yy, xx).real;
        RFLOAT d = DIRECT_A2D_ELEM(img, yy, xx).imag;
        RFLOAT ac = a * c;
        RFLOAT bd = b * d;
        RFLOAT ab_cd = (a + b) * (c + d);

        DIRECT_A2D_ELEM(img, yy, xx) = Complex(ac - bd, ab_cd - ac - bd);
    }
}

void FourierHelper::FourierShift2D(MultidimArray<RFLOAT> &img, RFLOAT xshift, RFLOAT yshift)
{
    FourierTransformer ft;
    MultidimArray<Complex> imgC;

    ft.FourierTransform(img, imgC);
    //FourierShift2D(imgC, xshift, yshift);
    shiftImageInFourierTransform(imgC, imgC, img.ydim, xshift, yshift);

    ft.inverseFourierTransform(imgC, img);
}

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

#ifndef CONVOLUTION_HELPER_H
#define CONVOLUTION_HELPER_H

#include <src/image.h>
#include <src/fftw.h>

class ConvolutionHelper
{
    public:

        static Image<RFLOAT> convolve2D(Image<RFLOAT>& img0, Image<RFLOAT>& img1, FourierTransformer &ft);

        static Image<RFLOAT> convolve2D(Image<RFLOAT>& img0, Image<RFLOAT>& img1);

        static Image<RFLOAT> gaussianKernel2D(double sigma, int w, int h,
                                              bool normalize = true,
                                              bool centered = false,
                                              bool half = false);

        inline static double sigmaFreq(double sigmaReal, int h)
        {
            return h/(2.0*PI*sigmaReal);
        }

        inline static double sigmaReal(double sigmaFreq, int h)
        {
            return h/(2.0*PI*sigmaFreq);
        }

};

#endif

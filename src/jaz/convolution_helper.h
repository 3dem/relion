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

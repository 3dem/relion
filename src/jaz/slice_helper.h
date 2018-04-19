#ifndef SLICE_HELPER_H
#define SLICE_HELPER_H

#include <src/image.h>
#include <src/fftw.h>
#include <src/strings.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>


class SliceHelper
{
    public:


        static void affineTransform(const Image<RFLOAT>& img, gravis::d4Matrix A, Image<RFLOAT>& dest);

        static void downsample(Image<RFLOAT>& img, Image<RFLOAT>& dest);
        static void downsampleSlices(const Image<RFLOAT>& img, Image<RFLOAT>& dest);
        static void downsampleSlicesReal(const Image<RFLOAT>& img, Image<RFLOAT>& dest);

        static void lowPassFilterSlicewise(Image<RFLOAT>& img, double maxFreq0, double maxFreq1);
        static void lowPassFilterSlice(Image<RFLOAT>& img, long int n, double maxFreq0, double maxFreq1);

        static void subsample(const Image<RFLOAT>& img, Image<RFLOAT>& dest);

        static void avgPad(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest, double ratio);
        static void avgPad2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest, double ratio);

        static void halveSpectrum2D(Image<Complex>& src, Image<Complex>& dest);

        static void extractSpectralSlice(Image<Complex>& src, Image<RFLOAT>& dest,
                                         gravis::d3Matrix proj, gravis::d2Vector volCentImg, double oversample = 4.0);

        static void insertSpectralSlices(std::vector<Image<RFLOAT> >& src,
                                         std::vector<gravis::d3Matrix> proj,
                                         std::vector<gravis::d2Vector> volCentImg,
                                         Image<Complex>& dest, double thickness = 1.0, double thicknessSlope = 0.0, double imgPad = 0.5);

        static void insertWeightedSpectralSlices(std::vector<Image<RFLOAT> >& src,
                                         std::vector<gravis::d3Matrix> proj,
                                         std::vector<gravis::d2Vector> volCentImg,
                                         std::vector<double> imgWeights,
                                         Image<Complex>& dest, double thickness = 1.0, double imgPad = 0.5);

        static void extractStackSlice(const Image<RFLOAT>& src, Image<RFLOAT>& dest, long int s);
        static void extractStackSlices(const Image<double>& src, Image<RFLOAT>& dest, long int s);
        static void extractStackSlices(const Image<float>& src, Image<RFLOAT>& dest, long int s);
        static Image<RFLOAT> getStackSlice(const Image<RFLOAT>& src, long int n);

        static void insertStackSlice(const Image<double>& src, Image<double>& dest, long int s);
        static void insertStackSlice(const Image<float>& src, Image<float>& dest, long int s);
        static void insertZSlice(const Image<double>& src, Image<double>& dest, long int s);
        static void insertZSlice(const Image<float>& src, Image<float>& dest, long int s);

        static Image<double> consolidate(const std::vector<Image<double> >& src, bool toN = false);
        static Image<float> consolidate(const std::vector<Image<float> >& src, bool toN = false);

        static void stat(const Image<RFLOAT>& img);
};

#endif

#ifndef TOMO_STACK_H
#define TOMO_STACK_H

#include <string>
#include <src/image.h>
#include <src/ctf.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>

class TomoStack
{
    public:

        TomoStack(){}
        TomoStack(std::string imagesFn, int imgCount, std::string angles, std::string affineTransforms, std::string ctfPath,
                  double angpix, double scaleFactor = 1.0, bool loadImgs = true);

        TomoStack extractSubStack(gravis::d3Vector center, int w, int h);

        std::vector<Image<RFLOAT> > images;
        std::vector<gravis::d4Matrix> affineXforms;
        std::vector<gravis::d4Matrix> tiltProjs;
        std::vector<gravis::d4Matrix> worldToImage;
        std::vector<CTF> ctfs;
        double angpix, scaleFactor;

        void downsample(int factor, int f0 = 0, int fc = -1);
        void phaseFlip(int f0 = 0, int fc = -1);
        void ctfModulate(int f0 = 0, int fc = -1);
        void wienerFilter(RFLOAT eps, RFLOAT Bfac, int f0 = 0, int fc = -1);
        void richardsonLucy(int iterations, RFLOAT eps, int f0 = 0, int fc = -1);
        void rampFilter(RFLOAT s0, RFLOAT t1, int f0 = 0, int fc = -1);
        void safeLog(RFLOAT eps, int f0 = 0, int fc = -1);
        void scaledExp(RFLOAT scale, int f0 = 0, int fc = -1);

        void defocusStack(int f, double dz0, double dz1, double eps, double Bfac, std::vector<Image<RFLOAT> >& dest, int x0 = 0, int y0 = 0, int w = -1, int h = -1);

        void saveImages(std::string path, int f0 = 0, int fc = -1);

        std::vector<std::vector<gravis::d2Vector> > loadFiducials(std::string file, double scale);
};

#endif

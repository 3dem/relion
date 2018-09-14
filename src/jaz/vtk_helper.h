#ifndef VTK_HELPER_H
#define VTK_HELPER_H

#include <src/image.h>
#include <src/strings.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/volume.h>
#include <src/jaz/tensor3x3.h>

class VtkHelper
{
    public:

        static Image<RFLOAT> allToZ(const Image<RFLOAT>& img);

        static void writeVTK(Image<double>& img, std::string fn,
                             double originX = 0.0, double originY = 0.0, double originZ = 0.0,
                             double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0,
                             bool binary = false);

        static void writeVTK(Image<float>& img, std::string fn,
                             double originX = 0.0, double originY = 0.0, double originZ = 0.0,
                             double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0,
                             bool binary = false);

        static void writeVTK(Image<Complex>& img, std::string fn,
                         double originX = 0.0, double originY = 0.0, double originZ = 0.0,
                         double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0,
                         bool binary = false);

        static void writeVTK(MultidimArray<RFLOAT>& img, std::string fn,
                             double originX = 0.0, double originY = 0.0, double originZ = 0.0,
                             double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0,
                             bool binary = false);

        static void writeVTK_Complex(const MultidimArray<Complex>& img, std::string fn, bool binary = false);
        static void writeVTK_d3(MultidimArray<gravis::t3Vector<RFLOAT> >& img, std::string fn, bool binary = false);
        static void writeTomoVTK(Image<RFLOAT>& img, std::string fn, bool binary = false, 
								 double pixelSize = 1.0, 
								 gravis::d3Vector origin = gravis::d3Vector(0.0,0.0,0.0));

        static void write(std::vector<Image<double> >& img, std::string fn,
                          double originX = 0.0, double originY = 0.0,
                          double spacingX = 1.0, double spacingY = 1.0,
                          bool binary = false);

        static void writeCentered(std::vector<Image<RFLOAT> >& img, std::string fn,
                          double originX = 0.0, double originY = 0.0,
                          double spacingX = 1.0, double spacingY = 1.0,
                          bool binary = false);

        static void write(std::vector<Image<float> >& img, std::string fn,
                          double originX = 0.0, double originY = 0.0,
                          double spacingX = 1.0, double spacingY = 1.0,
                          bool binary = false);


        static void readVTK(std::string fn, Volume<RFLOAT>& vol, gravis::d3Vector& origin, gravis::d3Vector& spacing);

        static void writeVTK(Volume<RFLOAT>& vol, std::string fn,
                             double originX = 0.0, double originY = 0.0, double originZ = 0.0,
                             double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0,
                             bool binary = false);

        static void writeVTK(Volume<gravis::t3Vector<RFLOAT> >& vol, std::string fn,
                             double originX = 0.0, double originY = 0.0, double originZ = 0.0,
                             double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0,
                             bool binary = false);

        static void writeVTK(Volume<Tensor3x3<RFLOAT> >& vol, std::string fn,
                             double originX = 0.0, double originY = 0.0, double originZ = 0.0,
                             double spacingX = 1.0, double spacingY = 1.0, double spacingZ = 1.0,
                             bool binary = false);
};

#endif

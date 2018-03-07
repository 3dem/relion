#ifndef STACK_HELPER_H
#define STACK_HELPER_H

#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/parallel_ft.h>
#include <vector>

class Projector;

class StackHelper
{
    public:

        typedef enum {BoxBin, GaussBin, FourierCrop} BinningType;

        static std::vector<MetaDataTable> splitByStack(const MetaDataTable* mdt);
        static std::vector<Image<RFLOAT> > loadStack(const MetaDataTable* mdt, std::string path = "", int threads = 1);
        static std::vector<Image<Complex> > loadStackFS(const MetaDataTable* mdt, std::string path = "",
                                                       int threads = 1,
                                                       std::vector<ParFourierTransformer>* fts = 0);
        static void saveStack(std::vector<Image<RFLOAT> >& stack, std::string fn);

        static std::vector<std::vector<Image<RFLOAT> > > loadMovieStack(const MetaDataTable* mdt, std::string moviePath);

        static std::vector<std::vector<Image<Complex> > > loadMovieStackFS(
                const MetaDataTable* mdt, std::string moviePath, bool center = false, int threads = 1,
                std::vector<ParFourierTransformer>* fts = 0);

        static std::vector<std::vector<Image<Complex> > > extractMovieStackFS(
                const MetaDataTable* mdt,
                std::string metaPath, std::string moviePath, std::string movie_ending,
                double outPs, double coordsPs, double moviePs,
                int squareSize, int threads,
                bool loadData = true, RFLOAT hot = -1.0, bool verbose = false);

        static Image<Complex> projectView(Projector* projector, const MetaDataTable* mdt, int index);
        //static Image<Complex> projectView2(Projector* projector, const MetaDataTable* mdt, int index);
        static std::vector<Image<Complex> > projectStack(Projector* projector, const MetaDataTable* mdt);
        static std::vector<Image<Complex> > projectStackPar(Projector* projector, const MetaDataTable* mdt, int numThreads);
        static std::vector<Image<Complex> > FourierTransform(std::vector<Image<RFLOAT> >& stack);
        static std::vector<Image<RFLOAT> > inverseFourierTransform(std::vector<Image<Complex> >& stack);

        static std::vector<Image<Complex> > applyBeamTilt(std::vector<Image<Complex> >& stack,
                                  RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
                                  RFLOAT tilt_x, RFLOAT tilt_y,
                                  RFLOAT tilt_xx, RFLOAT tilt_xy, RFLOAT tilt_yy);

        static std::vector<Image<Complex> > applyBeamTiltPar(std::vector<Image<Complex> >& stack,
                                  RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
                                  RFLOAT tilt_x, RFLOAT tilt_y,
                                  RFLOAT tilt_xx, RFLOAT tilt_xy, RFLOAT tilt_yy,
                                  int numThreads);

        static void varianceNormalize(std::vector<Image<Complex>>& movie, bool circleCropped = false);

        static std::vector<double> powerSpectrum(const std::vector<std::vector<Image<Complex> > > &stack);
        static std::vector<double> varSpectrum(const std::vector<std::vector<Image<Complex> > > &stack);
        static std::vector<double> powerSpectrum(const std::vector<std::vector<Image<Complex>>>& obs,
                                                 const std::vector<Image<Complex> >& signal);
};

#endif

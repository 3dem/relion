#ifndef REFINEMENT_PROGRAM_H
#define REFINEMENT_PROGRAM_H

#include <string>

#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/backprojector.h>

#include <src/jaz/obs_model.h>
#include <src/jaz/stack_helper.h>

#include <omp.h>

class RefinementProgram
{
    public:

        RefinementProgram(bool singleReference = false, bool doesMovies = false);

            // options:

            bool singleReference, doesMovies, debug, applyTilt, anisoTilt, useFsc,
                optStar, noStar, optReference, noReference, noTilt,
                preextracted, coordsAtMgRes;

            long maxMG, minMG;

            RFLOAT angpix, paddingFactor,
                beamtilt_x, beamtilt_y,
                beamtilt_xx, beamtilt_xy, beamtilt_yy,
                hotCutoff;

            int nr_omp_threads;

            double movie_scale, movie_angpix, coords_angpix;

            std::string
                starFn, reconFn0, reconFn1, maskFn,
                outPath, imgPath, fscFn,
                meta_path, movie_ending,
                movie_toReplace, movie_replaceBy;

            // data:

            Image<RFLOAT> maps[2];
            Image<RFLOAT> powSpec[2];
            Image<RFLOAT> freqWeight;
            Projector projectors[2];

            MetaDataTable mdt0;
            std::vector<MetaDataTable> mdts;
            RFLOAT Cs, lambda, kV;
            ObservationModel obsModel;

            int s, sh, fc;
            long g0, gc;


        int init(int argc, char *argv[]);
        int run();

        virtual int readMoreOptions(IOParser& parser, int argc, char *argv[]) {}
        virtual int _init(){return 0;}
        virtual int _run() = 0;

        double angstToPixFreq(double a);
        double pixToAngstFreq(double p);
};

#endif

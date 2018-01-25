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

#include <omp.h>

class RefinementProgram
{
    public:

        RefinementProgram(bool singleReference = false);

            // options:

            bool singleReference, debug, applyTilt, useFsc;

            long maxMG, minMG;
            RFLOAT angpix, paddingFactor, beamtilt_x, beamtilt_y;
            int nr_omp_threads;

            std::string starFn, reconFn0, reconFn1, maskFn, outPath, imgPath, fscFn;

            // data:

            Image<RFLOAT> maps[2];
            Image<RFLOAT> powSpec[2];
            Image<RFLOAT> freqWeight;
            Projector projectors[2];

            MetaDataTable mdt0;
            std::vector<MetaDataTable> mdts;
            RFLOAT Cs, lambda, kV;
            ObservationModel obsModel;

            int s, sh;
            long g0, gc;


        int init(int argc, char *argv[]);
        int run();

        virtual void readMoreOptions(IOParser& parser, int argc, char *argv[]) {}
        virtual int _init(){return 0;}
        virtual int _run() = 0;
};

#endif

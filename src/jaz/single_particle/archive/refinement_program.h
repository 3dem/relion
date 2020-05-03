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

#include <src/jaz/single_particle/legacy_obs_model.h>
#include <src/jaz/single_particle/stack_helper.h>

#include "src/micrograph_model.h"

#include <omp.h>

class RefinementProgram
{
    public:

        RefinementProgram(bool singleReference = false, bool doesMovies = false);

            // options:

            bool singleReference, doesMovies, debug, debugMov, applyTilt, anisoTilt, useFsc,
                optStar, noStar, optReference, noReference, noTilt,
                preextracted, coordsAtMgRes, hasCorrMic, saveMem;

            long maxMG, minMG;
            int firstFrame, lastFrame;

            RFLOAT angpix, paddingFactor,
                beamtilt_x, beamtilt_y,
                beamtilt_xx, beamtilt_xy, beamtilt_yy,
                hotCutoff;

            int nr_omp_threads;

            double movie_angpix, coords_angpix;

            std::string
                starFn, reconFn0, reconFn1, maskFn,
                outPath, imgPath, fscFn,
                meta_path, movie_ending,
                movie_toReplace, movie_replaceBy,
                corrMicFn, gain_path, last_gainFn;

            std::map<std::string, std::string> mic2meta;

            // data:

            Image<RFLOAT> maps[2];
            Image<RFLOAT> powSpec[2];
            Image<RFLOAT> freqWeight, lastGainRef;
            Projector projectors[2];

            Micrograph micrograph;

            MetaDataTable mdt0;
            std::vector<MetaDataTable> mdts;
            RFLOAT Cs, lambda, kV;
            LegacyObservationModel obsModel;

            int s, sh, fc;
            long g0, gc;

            std::vector<double> freqWeight1D;


        int init(int argc, char *argv[]);
        int run();

        virtual int readMoreOptions(IOParser& parser, int argc, char *argv[]) {return 0;}
        virtual int _init(){return 0;}
        virtual int _run() = 0;

        double angstToPixFreq(double a);
        double pixToAngstFreq(double p);

        void loadInitialMovieValues();

        std::vector<std::vector<Image<Complex>>> loadMovie(
                int g, int pc, std::vector<ParFourierTransformer>& fts);

        void setForAll(EMDLabel label, RFLOAT value);

        std::string getMicrographTag(int m);
};

#endif

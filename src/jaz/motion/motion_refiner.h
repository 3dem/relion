/***************************************************************************
 *
 * Author: "Jasenko Zivanov & Sjors H.W. Scheres"
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

#ifndef MOTION_REFINER_H_
#define MOTION_REFINER_H_


#include <src/ctf.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/backprojector.h>
#include <src/micrograph_model.h>
#include <src/jaz/obs_model.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/parallel_ft.h>

#include "motion_param_estimator.h"
#include "motion_estimator.h"

#include <omp.h>


class MotionRefiner
{
    public:

        MotionRefiner();


        // Verbosity
        int verb;

        // Allow continuation of crashed jobs
        bool only_do_unfinished;

        // Write out debugging information
        bool debug, debugMov;

        bool preextracted, hasCorrMic, saveMem;

        double movie_angpix, coords_angpix;

        long maxMG, minMG;

        RFLOAT angpix, paddingFactor, hotCutoff;

        int nr_omp_threads, firstFrame, lastFrame;

        std::string
            starFn, reconFn0, reconFn1, maskFn,
            outPath, fscFn,
            meta_path,
            movie_path, movie_ending, movie_toReplace, movie_replaceBy,
            corrMicFn, gain_path, last_gainFn;

        std::map<std::string, std::string> mic2meta;
        Micrograph micrograph;

        // For recombining frames
        bool doCombineFrames, hasBfacs, bfac_debug;
        int k0, k1;
        double k0a, k1a;
        std::string trackFn, bfacFn;

        // data:
        Image<RFLOAT> maps[2];
        Image<RFLOAT> powSpec[2];
        Image<RFLOAT> freqWeight, lastGainRef;
        std::vector<double> freqWeight1D;
        Projector projectors[2];

        MetaDataTable mdt0;
        std::vector<MetaDataTable> mdts;

        ObservationModel obsModel;

        MotionParamEstimator motionParamEstimator;
        MotionEstimator motionEstimator;

        // Q: Jasenko, can we have more informative names for these important variables?
        // A: They are so important and common that their names should be short!
        // (s: full image size, sh: half-size + 1, fc: frame count
        //  - these are consistent throughout the codebase.)
        int s, sh, fc;

        int micrograph_xsize, micrograph_ysize;



        // Read command line arguments
        void read(int argc, char **argv);

        // Initialise some general stuff after reading
        void init();

        // Re-initialise vector mdts to allow only_do_unfinished for combine_frames
        void initialiseCombineFrames();

        // General Running
        void run();


        double angToPix(double a);
        double pixToAng(double p);


        // Combine frames on a subset of the micrographa
        void combineFramesSubsetMicrographs(long g_start, long g_end);

        // combine all EPS files into one logfile.pdf
        void combineEPSAndSTARfiles();


        // For original particle-polishing-like Bfactors
        void calculateSingleFrameReconstruction(int iframe);


        // Get output STAR file name for the gth entry in the mdts
        FileName getOutputFileNameRoot(long int g);

        std::string getMicrographTag(long g);

        // load a movie and extract all particles
        // returns a per-particle vector of per-frame images of size sh x s
        std::vector<std::vector<Image<Complex>>> loadMovie(
                long g, int pc, std::vector<ParFourierTransformer>& fts);

        // load the header of the first movie only to learn the frame number and micrograph size
        // (also the dose and movie pixel size if a corrected_micrographs.star is available)
        // - this should be given by the global part of run_data.star, once it exists
        std::vector<std::vector<Image<Complex>>> loadInitialMovie();



        // recombine_frames
        std::vector<Image<RFLOAT>> weightsFromFCC();

        std::vector<Image<RFLOAT>> weightsFromBfacs();


    protected:

        // apply changes to micrograph-filenames implied by
        // movie_path, movie_ending, movie_toReplace, movie_replaceBy
        void adaptMovieNames();
};



#endif

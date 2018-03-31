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

#ifndef MOTION_FITTER_H_
#define MOTION_FITTER_H_


#include "src/ctf.h"
#include "src/image.h"
#include "src/fftw.h"
#include "src/backprojector.h"
#include "src/micrograph_model.h"
#include "src/jaz/obs_model.h"
#include "src/jaz/gravis/t2Vector.h"
#include "src/jaz/parallel_ft.h"
#include "src/jaz/motion_param_estimator.h"
// Includes not required in the header moved to the .cpp file.
// This prevents recompilation cascades.

#include <omp.h>

//using namespace gravis;

/*
   Do not use 'using' in headers!

   Headers are contagious (they include each other).
   This can lead to unexpected symbols from another namespace
   overwriting the expected symbols.
   Names in namespaces are chosen without regard for collisions.

   harmless example: abs (returns ints) vs. std::abs (returns doubles)

    --JZ
*/


class MotionFitter
{
    public:

        MotionFitter();


        // I/O Parser
        IOParser parser;

        // Verbosity
        int verb;

        // Allow continuation of crashed jobs
        bool only_do_unfinished;

        // Write out debugging information
        bool debug, debugMov;

        bool unregGlob, noGlobOff,
            debugOpt, diag, expKer, global_init,
            preextracted, hasCorrMic, saveMem;

        int maxEDs, maxIters;

        double dmga, dmgb, dmgc, dosePerFrame,
            sig_vel, sig_div, sig_acc,
            k_cutoff, k_cutoff_Angst, optEps,
            movie_angpix, coords_angpix;

        long maxMG, minMG;

        RFLOAT angpix, paddingFactor, hotCutoff;

        int nr_omp_threads, firstFrame, lastFrame;

        std::string
            starFn, reconFn0, reconFn1, maskFn,
            outPath, imgPath, fscFn,
            meta_path,
            movie_ending, movie_toReplace, movie_replaceBy,
            corrMicFn, gain_path, last_gainFn;

        std::map<std::string, std::string> mic2meta;
        std::vector<Image<RFLOAT>> dmgWeight, dmgWeightEval;
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

        // Q: Jasenko, can we have more informative names for these important variables?
        // A: They are so important and common that their names should be short!
        // (s: full image size, sh: half-size + 1, fc: frame count
        //  - these are consistent throughout the codebase.)
        int s, sh, fc;
        int micrograph_xsize, micrograph_ysize;


    public:

        // Read command line arguments
        void read(int argc, char **argv);

        // Print usage instructions
        void usage();

        // Initialise some general stuff after reading
        void initialise();

        // Re-initialise vector mdts to allow only_do_unfinished for combine_frames
        void initialiseCombineFrames();

        // General Running
        void run();

        // Fit CTF parameters for all particles on a subset of the micrographs
        void processSubsetMicrographs(long g_start, long g_end);

        // Combine frames on a subset of the micrographa
        void combineFramesSubsetMicrographs(long g_start, long g_end);

        // combine all EPS files into one logfile.pdf
        void combineEPSAndSTARfiles();

        // For original particle-polishing-like Bfactors
        void calculateSingleFrameReconstruction(int iframe);

        // load micrograph g and compute all data required for the optimization;
        // also used by MotionParamEstimator
        void prepMicrograph(
            // in:
            long g, std::vector<ParFourierTransformer>& fts,
            const std::vector<Image<RFLOAT>>& dmgWeight,
            // out:
            std::vector<std::vector<Image<Complex>>>& movie,
            std::vector<std::vector<Image<RFLOAT>>>& movieCC,
            std::vector<gravis::d2Vector>& positions,
            std::vector<std::vector<gravis::d2Vector>>& initialTracks,
            std::vector<gravis::d2Vector>& globComp);

        // Actual optimization method; also used by MotionParamEstimator
        std::vector<std::vector<gravis::d2Vector>> optimize(
            const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
            const std::vector<std::vector<gravis::d2Vector>>& inTracks,
            double sig_vel_px, double sig_acc_px, double sig_div_px,
            const std::vector<gravis::d2Vector>& positions,
            const std::vector<gravis::d2Vector>& globComp);


    private:

        // Get output STAR file name for the gth entry in the mdts
        FileName getOutputFileNameRoot(long int g);

        std::string getMicrographTag(long g);

        std::vector<std::vector<Image<Complex>>> loadMovie(
                long g, int pc, std::vector<ParFourierTransformer>& fts, int only_this_frame = -1);

        void updateFCC(
                const std::vector<std::vector<Image<Complex>>>& movie,
                const std::vector<std::vector<gravis::d2Vector>>& tracks,
                const MetaDataTable& mdt,
                std::vector<Image<RFLOAT>>& tables,
                std::vector<Image<RFLOAT>>& weights0,
                std::vector<Image<RFLOAT>>& weights1);

        void writeOutput(
                const std::vector<std::vector<gravis::d2Vector>>& tracks,
                const std::vector<Image<RFLOAT>>& fccData,
                const std::vector<Image<RFLOAT>>& fccWeight0,
                const std::vector<Image<RFLOAT>>& fccWeight1,
                const std::vector<gravis::d2Vector>& positions,
                std::string outPath, long g,
                double visScale);

        // recombine_frames
        std::vector<Image<RFLOAT>> weightsFromFCC();

        std::vector<Image<RFLOAT>> weightsFromBfacs();
};



#endif /* MOTIONFITTER_H_ */

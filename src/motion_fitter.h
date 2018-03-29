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
#include <src/jaz/image_log.h>
#include <src/jaz/slice_helper.h>
#include <src/jaz/spectral_helper.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/backprojection_helper.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/fsc_helper.h>
#include <src/jaz/gp_motion_fit.h>
#include <src/jaz/gradient_descent.h>
#include <src/jaz/motion_refinement.h>
#include <src/jaz/image_op.h>
#include <src/jaz/parallel_ft.h>

#include <omp.h>

using namespace gravis;


class MotionFitter
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Allow continuation of crashed jobs
	bool only_do_unfinished;

	// Write out debugging information
	bool debug, debugMov;

	// FOR NOW: copied all params from Jasenko's refinement_program class
    bool unregGlob, noGlobOff,
        paramEstim2, paramEstim3,
        debugOpt, diag, expKer, global_init,
    	preextracted, coordsAtMgRes, hasCorrMic, saveMem;

    int maxIters, paramEstimIters, paramEstimSteps;
    double dmga, dmgb, dmgc, dosePerFrame,
        sig_vel, sig_div, sig_acc,
        k_cutoff, maxStep, maxDistDiv,
        param_rV, param_rD, param_rA,
		movie_angpix, coords_angpix;

    long maxMG, minMG;

    RFLOAT angpix, paddingFactor, hotCutoff;

    int nr_omp_threads, bin, firstFrame, lastFrame, coords_bin, movie_bin;

    std::string
        starFn, reconFn0, reconFn1, maskFn,
        outPath, imgPath, fscFn,
        meta_path, bin_type_str,
		movie_ending, movie_toReplace, movie_replaceBy,
		corrMicFn, gain_path, last_gainFn;

    std::map<std::string, std::string> mic2meta;
    std::vector<Image<RFLOAT>> dmgWeight;
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

    // Jasenko, can we have more informative names for these important variables?
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


	// Helper functions
private:

	// Get output STAR file name for the gth entry in the mdts
	FileName getOutputFileNameRoot(long int g);

    std::string getMicrographTag(long g);

    std::vector<std::vector<Image<Complex>>> loadMovie(
	        long g, int pc, std::vector<ParFourierTransformer>& fts, int only_this_frame = -1);

	void prepMicrograph(
	        long g, std::vector<ParFourierTransformer>& fts,
	        const std::vector<Image<RFLOAT>>& dmgWeight,
	        std::vector<std::vector<Image<Complex>>>& movie,
	        std::vector<std::vector<Image<RFLOAT>>>& movieCC,
	        std::vector<d2Vector>& positions,
	        std::vector<std::vector<d2Vector>>& initialTracks,
	        std::vector<d2Vector>& globComp);

    std::vector<std::vector<d2Vector>> optimize(
            const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
            const std::vector<std::vector<d2Vector>>& inTracks,
            double sig_vel_px, double sig_acc_px, double sig_div_px,
            const std::vector<d2Vector>& positions,
            const std::vector<d2Vector>& globComp,
            double step, double minStep, double minDiff,
            long maxIters, double inertia);

    void updateFCC(
            const std::vector<std::vector<Image<Complex>>>& movie,
            const std::vector<std::vector<d2Vector>>& tracks,
            const MetaDataTable& mdt,
            std::vector<Image<RFLOAT>>& tables,
            std::vector<Image<RFLOAT>>& weights0,
            std::vector<Image<RFLOAT>>& weights1);

    void writeOutput(
            const std::vector<std::vector<d2Vector>>& tracks,
            const std::vector<Image<RFLOAT>>& fccData,
            const std::vector<Image<RFLOAT>>& fccWeight0,
            const std::vector<Image<RFLOAT>>& fccWeight1,
            const std::vector<d2Vector>& positions,
            std::string outPath, long g,
            double visScale);

    // recombine_frames
    std::vector<Image<RFLOAT>> weightsFromFCC();

    std::vector<Image<RFLOAT>> weightsFromBfacs();



};



#endif /* MOTIONFITTER_H_ */

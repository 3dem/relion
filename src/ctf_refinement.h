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

#ifndef CTF_REFINEMENT_H_
#define CTF_REFINEMENT_H_


#include "src/ctf.h"
#include "src/image.h"
#include "src/fftw.h"
#include <src/backprojector.h>
#include <src/jaz/image_log.h>
#include <src/jaz/slice_helper.h>
#include <src/jaz/spectral_helper.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/backprojection_helper.h>
#include <src/jaz/volume_converter.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/resampling_helper.h>
#include <src/jaz/ctf_helper.h>
#include <src/jaz/defocus_refinement.h>
#include <src/jaz/magnification_refinement.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/tilt_refinement.h>
#include <src/jaz/motion_refinement.h>
#include <src/jaz/image_op.h>
#include <src/jaz/refinement_program.h>
#include <src/jaz/parallel_ft.h>

#include <omp.h>

using namespace gravis;


class CtfRefiner
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Input & Output rootname
	FileName fn_opt, fn_fsc, fn_out;

	// Allow continuation of crashed jobs
	bool only_do_unfinished;

	// Perform per-particle defocus estimation?
	bool do_defocus_fit;

	// Perform beamtilt estimation?
	bool do_tilt_fit;

    // FOR NOW: copied all params from Jasenko's refinement_program class
	// TODO: throw away not-needed ones
	bool singleReference, doesMovies, debug, applyTilt, anisoTilt, useFsc,
        optStar, noStar, optReference, noReference, noTilt,
        preextracted, nogain, ctfTilt;

    long maxMG, minMG;

    RFLOAT angpix, paddingFactor,
        beamtilt_x, beamtilt_y,
        beamtilt_xx, beamtilt_xy, beamtilt_yy,
        hotCutoff;

    int nr_omp_threads, bin, coords_bin, movie_bin;

    std::string
        starFn, reconFn0, reconFn1, maskFn,
        outPath, imgPath, fscFn,
        meta_path, bin_type_str;

    StackHelper::BinningType binType;

    // data:

    Image<RFLOAT> maps[2];
    Image<RFLOAT> powSpec[2];
    Image<RFLOAT> freqWeight;
    std::vector<double> freqWeight1D;
    Projector projectors[2];

    MetaDataTable mdt0;
    std::vector<MetaDataTable> mdts;
    RFLOAT Cs, lambda, kV;
    ObservationModel obsModel;

    // Jasenko, can we have more informative names for these important variables?
    int s, sh, fc;
    long g0, gc;

    // Defocus_fit options
    RFLOAT defocusRange;
    bool fitAstigmatism, noGlobAstig, diag;


    // Tilt fit options
    RFLOAT kmin,
        testtilt_x, testtilt_y,
        testtilt_xx, testtilt_xy, testtilt_yy;

    bool precomputed, aniso;
    std::string precomp;

    Image<Complex> lastXY;
    Image<RFLOAT> lastW;


public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	// Fir defocus for all particles on one micrograph
	void fitDefocusOneMicrograph(long g, const std::vector<Image<Complex> > &obsF, const std::vector<Image<Complex> > &preds);

	// Perform beamtilt calculations for one micrograph
	void fitBeamtiltOneMicrograph(long g, const std::vector<Image<Complex> > &obsF, const std::vector<Image<Complex> > &pred,
			 std::vector<Image<Complex> > &xyAcc, std::vector<Image<RFLOAT> > &wAcc);

	// After sums of phase shifts have been accumulated over all micrographs: fit the actual beamtilt
	void fitBeamTiltFromSumsAllMicrographs(Image<Complex> &xyAccSum, Image<RFLOAT> &wAccSum);

	// Fit CTF parameters for all particles on a subset of the micrographs micrograph
	void processSubsetMicrographs(long g_start, long g_end, Image<Complex> &xyAccSum, Image<RFLOAT> &wAccSum);

	// General Running
	void run();


};



#endif /* CTF_REFINEMENT_H_ */

/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
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

#ifndef MOTIONCORR_RUNNER_H_
#define MOTIONCORR_RUNNER_H_

#include <glob.h>
#include <vector>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <algorithm>
#include <src/time.h>
#include "src/metadata_table.h"
#include "src/image.h"
#include "src/micrograph_model.h"
#include <src/jaz/single_particle/obs_model.h>

class MotioncorrRunner
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Number of threads per process
	int n_threads;
	int max_io_threads;

	// Output rootname
	FileName fn_in, fn_out, fn_movie;

	// Filenames of all the micrographs to run Motioncorr on
	std::vector<FileName> fn_micrographs, fn_ori_micrographs;

	// Optics group number for all original micrographs
	std::vector<int> optics_group_micrographs, optics_group_ori_micrographs;

	// Information about the optics groups
	ObservationModel obsModel;

	// Skip generation of logfile
	bool do_skip_logfile;

	// Use our own implementation
	bool do_own;
	bool interpolate_shifts;

	// Write in float16 (MRC mode 12)?
	bool write_float16;

	// Maximum number of iterations
	int max_iter;

	// Save aligned but non-dose weighted micrograph.
	// With MOTIONCOR2, this flag is always assumed to be true
	bool save_noDW;

	// Use MOTIONCOR2 instead of UNBLUR?
	bool do_motioncor2;

	// First and last movie frames to use in alignment and written-out corrected average and movie (default: do all)
	int first_frame_ali, last_frame_ali, first_frame_sum, last_frame_sum;

	// Group this number of frames and write summed power spectrum. -1 == do not write
	int grouping_for_ps;
	int ps_size;

	// Binning factor for binning inside MOTIONCORR/MOTIONCOR2
	double bin_factor;

	// Do binning before processing
	bool early_binning;

	// B-factor for MOTIONCOR2
	double bfactor;

	// Downsampling rate of CCF
	double ccf_downsample;

	// Dose at which to distinguish between early/late global motion in output statistics
	double dose_motionstats_cutoff;

	// Additional arguments that need to be passed to MOTIONCORR
	FileName fn_other_motioncor2_args;

	// MOTIONCOR2 executable
	FileName fn_motioncor2_exe;

	// Voltage and dose per frame for MOTIONCOR2/UNBLUR dose-weighting
	bool do_dose_weighting;
	double voltage;
	double dose_per_frame;
	double pre_exposure;

	// Gain reference file
	FileName fn_gain_reference;
	int gain_rotation, gain_flip;

	// Defect file
	FileName fn_defect;

	// Skip hot pixel detection in own motioncorr
	bool skip_defect;

	// Archive directory
	FileName fn_archive;

	// Number of patches in X, Y direction for MOTIONCOR2
	int patch_x, patch_y;

	// How many frames to group in MOTIONCOR2
	int group;

	// Pixel size for UNBLUR
	double angpix;

	// Continue an old run: only estimate CTF if logfile WITH Final Values line does not yet exist, otherwise skip the micrograph
	bool continue_old;

	// Process at most this number of (unprocessed) micrographs
	long do_at_most;

	// EER parameters
	int eer_upsampling, eer_grouping;

	// Output STAR file
	MetaDataTable MDavg, MDmov;

	// Which GPU devices to use?
	int devCount;
	std::string gpu_ids;
	std::vector < std::vector < std::string > > allThreadIDs;

	// Read command line arguments
	void read(int argc, char **argv, int rank = 0);

	// Print usage instructions
	void usage();

	// Initialise some stuff after reading
	void initialise();

	void prepareGainReference(bool write_gain);

	// Execute all MOTIONCORR jobs
	void run();

	// Given an input fn_mic filename, this function will determine the names of the output corrected image (fn_avg) and the corrected movie (fn_mov).
	FileName getOutputFileNames(FileName fn_mic);

	// Execute MOTIONCOR2 for a single micrograph
	bool executeMotioncor2(Micrograph &mic, int rank = 0);

	// Get the shifts from MOTIONCOR2
	void getShiftsMotioncor2(FileName fn_log, Micrograph &mic);

	// Execute our own implementation for a single micrograph
	bool executeOwnMotionCorrection(Micrograph &mic);

	// Plot the shifts
	void plotShifts(FileName fn_mic, Micrograph &mic);

	// Save micrograph model
	void saveModel(Micrograph &mic);

	// Make a PDF file with all the shifts and write output STAR files
	void generateLogFilePDFAndWriteStarFiles();

	// Write out final STAR file
	void writeSTAR();

	// Read fn_defect (defect map, where 1 is bad, or defect text in the UCSF MotionCor2 format, x y w h) and fill bBad.
	static void fillDefectMask(MultidimArray<bool> &bBad, FileName fn_defect, int n_threads=1);

	// Check if fn_defect is Serial EM's defect file
	static bool detectSerialEMDefectText(FileName fn_defect);

private:
	// shiftx, shifty is relative to the (real space) image size
	void shiftNonSquareImageInFourierTransform(MultidimArray<fComplex> &frame, RFLOAT shiftx, RFLOAT shifty);

	bool alignPatch(std::vector<MultidimArray<fComplex> > &Fframes, const int pnx, const int pny, const RFLOAT scaled_B, std::vector<RFLOAT> &xshifts, std::vector<RFLOAT> &yshifts, std::ostream &logfile);

	void binNonSquareImage(Image<float> &Iwork, RFLOAT bin_factor);

	int findGoodSize(int request);

	void doseWeighting(std::vector<MultidimArray<fComplex> > &Fframes, std::vector<RFLOAT> doses, RFLOAT apix);

	void realSpaceInterpolation(Image <float> &Isum, std::vector<Image<float> > &Iframes, MotionModel *model, std::ostream &logfile);

	void realSpaceInterpolation_ThirdOrderPolynomial(Image <float> &Isum, std::vector<Image<float> > &Iframes, ThirdOrderPolynomialModel &model, std::ostream &logfile);

	void interpolateShifts(std::vector<int> &group_start, std::vector<int> &group_size,
	                       std::vector<RFLOAT> &xshifts, std::vector<RFLOAT> &yshifts,
	                       int n_frames,
	                       std::vector<RFLOAT> &interpolated_xshifts, std::vector<RFLOAT> &interpolated_yshifts);
};


#endif /* MOTIONCORR_RUNNER_H_ */

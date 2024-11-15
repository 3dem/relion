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
#ifndef RELION_ALIGN_TILTSERIES_RUNNER_H
#define RELION_ALIGN_TILTSERIES_RUNNER_H


#include  <glob.h>
#include  <vector>
#include  <string>
#include  <stdlib.h>
#include  <unistd.h>
#include  <stdio.h>
#include <src/metadata_table.h>
#include <src/image.h>
#include <src/time.h>
#include <src/jaz/single_particle/obs_model.h>
#include "src/jaz/tomography/tomogram_set.h"
const std::string fiducial_directive =
        "setupset.copyarg.userawtlt = 1 \n"
        "setupset.copyarg.stackext = mrc \n"
//        "setupset.copyarg.rotation = 85 \n"
//        "setupset.copyarg.pixel = 0.0675 \n"
//        "setupset.copyarg.gold = 10 \n"
        "setupset.copyarg.dual = 0 \n"
        "runtime.Positioning.any.wholeTomogram = 1 \n"
        "runtime.Fiducials.any.trackingMethod = 0 \n"
        "runtime.Fiducials.any.seedingMethod = 3 \n"
        "runtime.Excludeviews.any.deleteOldFiles = 0 \n"
        "runtime.AlignedStack.any.binByFactor = 8 \n"
        "comparam.prenewst.newstack.BinByFactor = 8 \n"
        "comparam.prenewst.newstack.AntialiasFilter = -1 \n"
        "comparam.autofidseed.autofidseed.TargetNumberOfBeads = 50 \n"
        "comparam.track.beadtrack.SobelFilterCentering = 1 \n"
        "comparam.track.beadtrack.ScalableSigmaForSobel = 0.12 \n"
        "comparam.newst.newstack.TaperAtFill = 1,1 \n"
        "comparam.newst.newstack.AntialiasFilter = -1 \n"
        "comparam.golderaser.ccderaser.ExpandCircleIterations = 3 \n"
        "comparam.eraser.ccderaser.PeakCriterion = 8.0 \n"
        "comparam.eraser.ccderaser.DiffCriterion = 6.0 \n"
        "comparam.align.tiltalign.TiltOption = 0 \n"
        "comparam.align.tiltalign.SurfacesToAnalyze = 1 \n"
        "comparam.align.tiltalign.RotOption = -1 \n"
        "comparam.align.tiltalign.MagOption = 0 \n"
        "comparam.align.tiltalign.BeamTiltOption = 2 \n";

const std::string patchtrack_directive =
        "setupset.copyarg.userawtlt = 1 \n"
        "setupset.copyarg.stackext = mrc \n"
//        "setupset.copyarg.rotation = 85 \n"
//        "setupset.copyarg.pixel = 0.0675 \n"
        "setupset.copyarg.gold = 10 \n"
        "setupset.copyarg.dual = 0 \n"
        "runtime.Positioning.any.wholeTomogram = 1 \n"
        "runtime.PatchTracking.any.adjustTiltAngles = 0 \n"
        "runtime.Fiducials.any.trackingMethod = 1 \n"
        "runtime.Fiducials.any.seedingMethod = 3 \n"
        "runtime.Excludeviews.any.deleteOldFiles = 0 \n"
        "runtime.AlignedStack.any.binByFactor = 8 \n"
//        "comparam.xcorr_pt.tiltxcorr.SizeOfPatchesXandY = 180,180 \n"
//        "comparam.xcorr_pt.tiltxcorr.OverlapOfPatchesXandY = 0.33 \n"
        "comparam.xcorr_pt.tiltxcorr.ShiftLimitsXandY = 200,200 \n"
        "comparam.xcorr_pt.tiltxcorr.LengthAndOverlap = 16,4 \n"
        "comparam.xcorr_pt.tiltxcorr.IterateCorrelations = 20 \n"
        "comparam.xcorr_pt.tiltxcorr.FilterSigma2 = 0.03 \n"
        "comparam.xcorr_pt.tiltxcorr.FilterRadius2 = 0.125 \n"
        "comparam.track.beadtrack.SobelFilterCentering = 1 \n"
        "comparam.track.beadtrack.ScalableSigmaForSobel = 0.12 \n"
//        "comparam.tilt.tilt.THICKNESS = 3000 \n"
//        "comparam.prenewst.newstack.BinByFactor = 8 \n"
        "comparam.prenewst.newstack.AntialiasFilter = -1 \n"
        "comparam.newst.newstack.TaperAtFill = 1,1 \n"
        "comparam.newst.newstack.AntialiasFilter = -1 \n"
        "comparam.golderaser.ccderaser.ExpandCircleIterations = 3 \n"
        "comparam.eraser.ccderaser.PeakCriterion = 8.0 \n"
        "comparam.eraser.ccderaser.DiffCriterion = 6.0 \n"
        "comparam.align.tiltalign.TiltOption = 0 \n"
        "comparam.align.tiltalign.SurfacesToAnalyze = 1 \n"
        "comparam.align.tiltalign.RotOption = -1 \n"
        "comparam.align.tiltalign.MagOption = 0 \n"
        "comparam.align.tiltalign.BeamTiltOption = 2 \n";


class AlignTiltseriesRunner
{
public:

    // I/O Parser
    IOParser parser;

    // Verbosity
    int verb;

    // Output rootname
    FileName fn_in, fn_out;

    // Filenames of all the tomograms to estimate the CTF from
    std::vector<long> idx_tomograms, idx_tomograms_all;

    // Information about tomography experiment
    TomogramSet tomogramSet;

    // Imod executable
    FileName fn_batchtomo_exe;

    // AreTomo executable
    FileName fn_aretomo_exe;

    // Filename for IMOD directives file
    FileName fn_adoc_template;

    std::string my_adoc_template;

    // Use IMOD:fiducials
    bool do_imod_fiducials;

    // Nominal fiducial diameter (nm)
    RFLOAT fiducial_diam;

    // Use IMOD:patch-tracking
    bool do_imod_patchtrack;

    // Unbinned patch size (pixels)
    RFLOAT patch_size;

    // Patch overlap (percentage 0-100)
    RFLOAT patch_overlap;

    // Use AreTomo
    bool do_aretomo;

    // Perform tilt angle correction in AreTomo
    bool do_aretomo_tiltcorrect;

    // User-specified value for tilt angle correction
    RFLOAT aretomo_tilcorrect_angle;

    // Do CTF estimation in aretomo?
    bool do_aretomo_ctf;

    // Do phase shift estimation in AreTomo?
    bool do_aretomo_phaseshift;

    // estimated tomogram thickness (for -AlignZ)
    RFLOAT tomogram_thickness;

    // Which GPU devices to use?
    int devCount;
    std::string gpu_ids;
    std::vector < std::vector < std::string > > allThreadIDs;

    // Additional gctf command line options
    std::string other_wrapper_args;

    // Continue an old run: only estimate CTF if logfile WITH Final Values line does not yet exist, otherwise skip the tomogram
    bool continue_old;

    // Process at most this number of unprocessed tomograms
    long do_at_most;

public:
    // Read command line arguments
    void read(int argc, char **argv, int rank = 0);

    // Print usage instructions
    void usage();

    // Initialise some stuff after reading
    void initialise(bool is_leader = true);

    // Execute all CTFFIND jobs to get CTF parameters
    void run();

    // Check STAR file for tomogram exists and has the correct labels
    bool checkResults(long idx_tomo);

    // Generate MRC stack and raw tilt file
    void generateMRCStackAndRawTiltFile(long idx_tomo, bool is_aretomo);

    // Execute IMOD for a single tomogram
    void executeIMOD(long idx_tomo, int rank = 0);

    // Execute AreTomo for a single tomogram
    void executeAreTomo(long idx_tomo, int rank = 0);

    // Make per-tiltseries EPS files for the logfile
    void makePerTiltSeriesEPSFiles(long idx_tomo, bool do_ctf = false);

    // Read IMOD results files and insert data into relion's MetaDataTable
    bool readIMODResults(long idx_tomo, std::string &error_message);

    // Read AreTomo results files and insert data into relion's MetaDataTable
    bool readAreTomoResults(long idx_tomo, std::string &error_message);

    // Harvest all AreTomo results into the single tomograms set starfile, and write it out
    void joinResults();


};


#endif //RELION_ALIGN_TILTSERIES_RUNNER_H

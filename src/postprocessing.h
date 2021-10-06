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

#ifndef POSTPROCESSING_H_
#define POSTPROCESSING_H_

#include "src/image.h"
#include "src/multidim_array.h"
#include "src/metadata_table.h"
#include "src/healpix_sampling.h"
#include "src/fftw.h"
#include "src/time.h"
#include "src/mask.h"
#include "src/funcs.h"
#include "src/CPlot2D.h"
#include "src/mpi.h"


class Postprocessing
{
public:

	// I/O Parser
	IOParser parser;

	// Verbosity
	int verb;

	// Input & Output rootname
	FileName fn_out, fn_I1, fn_I2, fn_OS;

	// Optimiser set table for tomography
	MetaDataTable optimisationSet;

	// Images for the two half-reconstructions and the mask
	Image<RFLOAT> I1, I2, Im;

	// Pixel size in Angstroms
	RFLOAT angpix, mtf_angpix;

	/////// Masking

	// Perform automated masking (based on a density threshold)
	bool do_auto_mask;

	// Write half maps as well
	bool write_halfmaps;

	// Density threshold below which to calculate initial mask seed
	RFLOAT ini_mask_density_threshold;

	// Number of pixels to extend the mask beyond the initial mask seed
	RFLOAT extend_ini_mask;

	// Width (in pixels) for soft mask edge
	RFLOAT width_soft_mask_edge;

	// From the resolution where the FSC drops below this value, randomize the phases in the two maps
	RFLOAT randomize_fsc_at;

	// Or user-specified resolution (in A) for randomisation of phases
	RFLOAT randomize_at_A;

	// Filename for a user-provided mask
	FileName fn_mask;

	// Fraction of pixels corresponding to protein (i.e. white) in the solvent mask or calculated from moweight
	RFLOAT frac_solvent_mask, frac_molweight;

	// The user-provided Molecular Weight (in kDa) or the molecular weight contained in the solvent mask
	RFLOAT molweight, molweight_frommask;

	// Use the mask even when the resolution becomes worse
	bool force_mask;

	/////// Local-resolution

	// Perform automated masking (based on a density threshold)
	bool do_locres;

	// Sampling rate for local resolution estimates (in A)
	RFLOAT locres_sampling;

	// Radius of mask for local resolution estimates (in A)
	RFLOAT locres_maskrad;

	// Width of soft edge on local resolution masks (in A)
	RFLOAT locres_edgwidth;

	// Randomize phases from this resolution for local resolution calculation
	RFLOAT locres_randomize_fsc;

	// Lowest resolution allowed in the locres map
	RFLOAT locres_minres;

	//////// Sharpening

	// Filename for the STAR-file with the MTF of the detector
	FileName fn_mtf;

	// Flag to indicate whether to perform Rosenthal&Henderson-2003 like B-factor estimation
	bool do_auto_bfac;

	// Flag to indicate whether we'll do masking
	bool do_mask;

	// Flag to indicate whether we're also doing amplitude correlations and DPR, in that case re-normalize amplitudes for n=uneven orientational distributions
	bool do_ampl_corr;

	// Minimum and maximum resolution to use in the fit
	RFLOAT fit_minres, fit_maxres;

	// User-provided (ad hoc) B-factor
	RFLOAT adhoc_bfac;

	///////// Filtering

	// Flag to indicate whether to use FSC-weighting before B-factor sharpening
	bool do_fsc_weighting;

	// Frequency at which to low-pass filter the final map
	RFLOAT low_pass_freq;

	// Width of raised cosine edge on low-pass filter
	int filter_edge_width;

	// Shell for randomisation of phases
	int randomize_at;

	// Allow the user to interpolate the FSC to obtain a fractional-shell resolution estimate (needed to evaluate new methods)
	bool do_interpolate;

	// Arrays to store FSC, Guinier curves etc
	MultidimArray<RFLOAT> fsc_unmasked, acorr_unmasked, acorr_masked, dpr_unmasked, dpr_masked;
	MultidimArray<RFLOAT> fsc_masked, fsc_random_masked, fsc_true, fsc_part_molweight, fsc_part_fracmask;
	RFLOAT global_intercept, global_slope, global_corr_coeff, global_bfactor, global_resol, fract_resol;
	// The Guinier plots
	std::vector<fit_point2D>  guinierin, guinierinvmtf, guinierweighted, guiniersharpen;

public:
	// Read command line arguments
	void read(int argc, char **argv);

	// Print usage instructions
	void usage();

	// Set parameters to some useful defaults
	void clear();

	// Initialise some stuff after reading
	void initialise();

	// Generate the mask (or read it from file)
	// Returns true if masking needs to be done, false otherwise
	bool getMask();

	// Divide by MTF and perform FSC-weighted B-factor sharpening, as in Rosenthal and Henderson, 2003
	// Returns the applied filter resolution
	RFLOAT sharpenMap();

	// Map nay 3D FFTW pixel onto the surface of a sphere with radius myradius_count
	bool findSurfacePixel(int idx, int kp, int ip, int jp,
			int &best_kpp, int &best_ipp, int &best_jpp,
			int myradius_count, int search=2);

	// Correct amplitudes inhomogeneity due to non-uniform orientational distributions
	void correctRadialAmplitudeDistribution(MultidimArray<RFLOAT > &I);

	// This divides the input FT by the mtf (if fn_mtf !="")
	void divideByMtf(MultidimArray<Complex > &FT);

	// Make a Guinier plot from the Fourier transform of an image
	void makeGuinierPlot(MultidimArray<Complex > &FT, std::vector<fit_point2D> &guinier);

	// Use Richard's formula to calculate FSC_true
	void calculateFSCtrue(MultidimArray<RFLOAT> &fsc_true, MultidimArray<RFLOAT> &fsc_unmasked,
			MultidimArray<RFLOAT> &fsc_masked, MultidimArray<RFLOAT> &fsc_random_masked, int randomize_at );

	// cisTEM-like FSC corrected for fraction of solvent mask
	void calculateFSCpart(const MultidimArray<RFLOAT> fsc_unmasked, RFLOAT fraction, MultidimArray<RFLOAT> &fsc_part);

	// Apply sqrt(2FSC/(FSC=1)) weighting prior to B-factor sharpening
	void applyFscWeighting(MultidimArray<Complex > &FT, MultidimArray<RFLOAT> my_fsc);

	// Output map and masked map
	void writeMaps(FileName fn_root);

	// Output map and STAR files with metadata, also write final resolution to screen
	void writeOutput();

	// Write XML file for EMDB submission
	void writeFscXml(MetaDataTable &MDfsc);

	// Write DAT file for easier plotting in xmgrace
	void writeFscDat(MetaDataTable &MDfsc);

	// Local-resolution running
	void run_locres(int rank = 0, int size = 1);

	// General Running
	void run();

};

#endif /* POSTPROCESSING_H_ */

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
/***************************************************************************
 *
 * Authors:	J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307 USA
 *
 *	All comments concerning this program package may be sent to the
 *	e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef METADATA_LABEL_H
#define METADATA_LABEL_H

#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include "src/funcs.h"

class EMDLabelData;
class StaticInitialization;

enum EMDLabel
{
	EMDL_UNDEFINED = -1, // Keep the order the same as in StaticInitialization below!!
	EMDL_FIRST_LABEL,
	EMDL_OBJID = EMDL_FIRST_LABEL, ///< object id (int), NOTE: This label is special and shouldn't be used

	EMDL_AREA_ID, ///< ID for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
	EMDL_AREA_NAME, ///< Name for the area (or field of view). If one does not use (tilt) series, area would be the same as micrograph...
	EMDL_COMMENT, // The EMDL_COMMENT is handled specially as well

	EMDL_BODY_MASK_NAME, ///< For multi-body refinements
	EMDL_BODY_KEEP_FIXED, ///< For multi-body refinements
	EMDL_BODY_REFERENCE_NAME,
	EMDL_BODY_ROTATE_DIRECTION_X,
	EMDL_BODY_ROTATE_DIRECTION_Y,
	EMDL_BODY_ROTATE_DIRECTION_Z,
	EMDL_BODY_ROTATE_RELATIVE_TO,
	EMDL_BODY_SIGMA_ANG,
	EMDL_BODY_SIGMA_OFFSET, // deprecated
	EMDL_BODY_SIGMA_OFFSET_ANGSTROM,
	EMDL_BODY_SIGMA_ROT,
	EMDL_BODY_SIGMA_TILT,
	EMDL_BODY_SIGMA_PSI,
	EMDL_BODY_STAR_FILE,

	EMDL_CTF_ASTIGMATISM,
	EMDL_CTF_BFACTOR, ///< B-factor
	EMDL_CTF_MAXRES, ///< Maximum resolution with Thon rings
	EMDL_CTF_VALIDATIONSCORE, ///< Gctf-based validation score for CTF fit
	EMDL_CTF_SCALEFACTOR, ///< linear scale-factor
	EMDL_CTF_SAMPLING_RATE, ///< Sampling rate
	EMDL_CTF_VOLTAGE, ///< Microscope voltage (kV)
	EMDL_CTF_DEFOCUSU, ///< Defocus U (Angstroms)
	EMDL_CTF_DEFOCUSV, ///< Defocus V (Angstroms)
	EMDL_CTF_DEFOCUS_ANGLE, ///< Defocus angle (degrees)
	EMDL_CTF_CS, ///< Spherical aberration
	EMDL_CTF_CA, ///< Chromatic aberration
	EMDL_CTF_DETECTOR_PIXEL_SIZE, ///< Pixel size for detector as used in CTF-determination (deprecated)
	EMDL_CTF_POWER_SPECTRUM,
	EMDL_CTF_ENERGY_LOSS, ///< Energy loss
	EMDL_CTF_FOM, ///< ctffind FOM (CC) for quality of CTF-fit
	EMDL_CTF_IMAGE, ///< name of an image describing the CTF model
	EMDL_CTF_LENS_STABILITY, ///< Lens stability
	EMDL_CTF_MAGNIFICATION, ///< Magnification used for CTF-determination (deprecated)
	EMDL_CTF_PHASESHIFT, ///< Phase-shift from a phase plate
	EMDL_CTF_CONVERGENCE_CONE, ///< Convergence cone
	EMDL_CTF_LONGITUDINAL_DISPLACEMENT, ///< Longitudinal displacement
	EMDL_CTF_TRANSVERSAL_DISPLACEMENT, ///< Transversal displacemente
	EMDL_CTF_Q0, ///< Amplitude contrast
	EMDL_CTF_K, ///< CTF gain
	EMDL_CTF_VALUE, ///< CTF value

	EMDL_IMAGE_NAME,
	EMDL_IMAGE_ORI_NAME,
	EMDL_IMAGE_RECONSTRUCT_NAME,
	EMDL_IMAGE_ID,
	EMDL_IMAGE_ENABLED,
	EMDL_IMAGE_DATATYPE,
	EMDL_IMAGE_DIMENSIONALITY,
	EMDL_IMAGE_BEAMTILT_X,
	EMDL_IMAGE_BEAMTILT_Y,
	EMDL_IMAGE_MTF_FILENAME,
	EMDL_IMAGE_OPTICS_GROUP,
	EMDL_IMAGE_OPTICS_GROUP_NAME,
	EMDL_IMAGE_ODD_ZERNIKE_COEFFS,
	EMDL_IMAGE_EVEN_ZERNIKE_COEFFS,
	EMDL_IMAGE_PIXEL_SIZE,
	EMDL_IMAGE_MAG_MATRIX_00,
	EMDL_IMAGE_MAG_MATRIX_01,
	EMDL_IMAGE_MAG_MATRIX_10,
	EMDL_IMAGE_MAG_MATRIX_11,

	EMDL_IMAGE_COORD_X,
	EMDL_IMAGE_COORD_Y,
	EMDL_IMAGE_COORD_Z,
	EMDL_IMAGE_FRAME_NR,
	EMDL_IMAGE_MAGNIFICATION_CORRECTION,
	EMDL_IMAGE_NORM_CORRECTION,
	EMDL_IMAGE_SAMPLINGRATE,
	EMDL_IMAGE_SAMPLINGRATE_X,
	EMDL_IMAGE_SAMPLINGRATE_Y,
	EMDL_IMAGE_SAMPLINGRATE_Z,
	EMDL_IMAGE_SIZE,
	EMDL_IMAGE_SIZE_X,
	EMDL_IMAGE_SIZE_Y,
	EMDL_IMAGE_SIZE_Z,
	EMDL_IMAGE_STATS_MIN,
	EMDL_IMAGE_STATS_MAX,
	EMDL_IMAGE_STATS_AVG,
	EMDL_IMAGE_STATS_STDDEV,
	EMDL_IMAGE_STATS_SKEW,
	EMDL_IMAGE_STATS_KURT,
	EMDL_IMAGE_WEIGHT,

	EMDL_JOB_IS_CONTINUE,
	EMDL_JOB_TYPE,
	EMDL_JOB_TYPE_NAME,

	EMDL_JOBOPTION_TYPE,
	EMDL_JOBOPTION_VARIABLE,
	EMDL_JOBOPTION_VALUE,
	EMDL_JOBOPTION_LABEL,
	EMDL_JOBOPTION_DEFAULT_VALUE,
	EMDL_JOBOPTION_MINVAL,
	EMDL_JOBOPTION_MAXVAL,
	EMDL_JOBOPTION_STEPVAL,
	EMDL_JOBOPTION_HELPTEXT,
	EMDL_JOBOPTION_PATTERN,
	EMDL_JOBOPTION_DIRECTORY,
	EMDL_JOBOPTION_MENUOPTIONS,

	EMDL_MATRIX_1_1,
	EMDL_MATRIX_1_2,
	EMDL_MATRIX_1_3,
	EMDL_MATRIX_2_1,
	EMDL_MATRIX_2_2,
	EMDL_MATRIX_2_3,
	EMDL_MATRIX_3_1,
	EMDL_MATRIX_3_2,
	EMDL_MATRIX_3_3,

	EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL,
	EMDL_MICROGRAPH_ACCUM_MOTION_EARLY,
	EMDL_MICROGRAPH_ACCUM_MOTION_LATE,
	EMDL_MICROGRAPH_ID,
	EMDL_MICROGRAPH_NAME,
	EMDL_MICROGRAPH_GAIN_NAME,
	EMDL_MICROGRAPH_DEFECT_FILE,
	EMDL_MICROGRAPH_NAME_WODOSE,
	EMDL_MICROGRAPH_MOVIE_NAME,
	EMDL_MICROGRAPH_METADATA_NAME,
	EMDL_MICROGRAPH_TILT_ANGLE,
	EMDL_MICROGRAPH_TILT_AXIS_DIRECTION,
	EMDL_MICROGRAPH_TILT_AXIS_OUTOFPLANE,
	EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE,
	EMDL_MICROGRAPH_PIXEL_SIZE,
	EMDL_MICROGRAPH_PRE_EXPOSURE,
	EMDL_MICROGRAPH_DOSE_RATE,
	EMDL_MICROGRAPH_BINNING,
	EMDL_MICROGRAPH_FRAME_NUMBER,
	EMDL_MICROGRAPH_MOTION_MODEL_VERSION,
	EMDL_MICROGRAPH_START_FRAME,
	EMDL_MICROGRAPH_END_FRAME,
	EMDL_MICROGRAPH_SHIFT_X,
	EMDL_MICROGRAPH_SHIFT_Y,
	EMDL_MICROGRAPH_MOTION_COEFFS_IDX,
	EMDL_MICROGRAPH_MOTION_COEFF,

	EMDL_MASK_NAME,

	EMDL_MLMODEL_ACCURACY_ROT,
	EMDL_MLMODEL_ACCURACY_TRANS, // deprecated
	EMDL_MLMODEL_ACCURACY_TRANS_ANGSTROM,
	EMDL_MLMODEL_AVE_PMAX,
	EMDL_MLMODEL_CURRENT_RESOLUTION,
	EMDL_MLMODEL_CURRENT_SIZE,
	EMDL_MLMODEL_DATA_VS_PRIOR_REF,
	EMDL_MLMODEL_DIMENSIONALITY,
	EMDL_MLMODEL_DIMENSIONALITY_DATA,
	EMDL_MLMODEL_DIFF2_HALVES_REF,
	EMDL_MLMODEL_ESTIM_RESOL_REF,
	EMDL_MLMODEL_FOURIER_COVERAGE_REF,
	EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF,
	EMDL_MLMODEL_FSC_HALVES_REF,
	EMDL_MLMODEL_GROUP_NAME,
	EMDL_MLMODEL_GROUP_NO,
	EMDL_MLMODEL_GROUP_NR_PARTICLES,
	EMDL_MLMODEL_GROUP_SCALE_CORRECTION,
	EMDL_MLMODEL_HELICAL_NR_ASU,
	EMDL_MLMODEL_HELICAL_TWIST,
	EMDL_MLMODEL_HELICAL_TWIST_MIN,
	EMDL_MLMODEL_HELICAL_TWIST_MAX,
	EMDL_MLMODEL_HELICAL_TWIST_INITIAL_STEP,
	EMDL_MLMODEL_HELICAL_RISE,
	EMDL_MLMODEL_HELICAL_RISE_MIN,
	EMDL_MLMODEL_HELICAL_RISE_MAX,
	EMDL_MLMODEL_HELICAL_RISE_INITIAL_STEP,
	EMDL_MLMODEL_IS_HELIX,
	EMDL_MLMODEL_INTERPOLATOR,
	EMDL_MLMODEL_LL,
	EMDL_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION,
	EMDL_MLMODEL_NORM_CORRECTION_AVG,
	EMDL_MLMODEL_NR_BODIES,
	EMDL_MLMODEL_NR_CLASSES,
	EMDL_MLMODEL_NR_GROUPS,
	EMDL_MLMODEL_ORIGINAL_SIZE,
	EMDL_MLMODEL_ORIENTABILITY_CONTRIBUTION,
	EMDL_MLMODEL_PADDING_FACTOR,
	EMDL_MLMODEL_PDF_CLASS,
	EMDL_MLMODEL_PRIOR_OFFX_CLASS,
	EMDL_MLMODEL_PRIOR_OFFY_CLASS,
	EMDL_MLMODEL_PDF_ORIENT,
	EMDL_MLMODEL_PIXEL_SIZE,
	EMDL_MLMODEL_POWER_REF,
	EMDL_MLMODEL_PRIOR_MODE,
	EMDL_MLMODEL_SIGMA_OFFSET, // deprecated
	EMDL_MLMODEL_SIGMA_OFFSET_ANGSTROM,
	EMDL_MLMODEL_SIGMA_ROT,
	EMDL_MLMODEL_SIGMA_TILT,
	EMDL_MLMODEL_SIGMA_PSI,
	EMDL_MLMODEL_REF_IMAGE,
	EMDL_MLMODEL_GRADIENT_MOMENT1_IMAGE,
	EMDL_MLMODEL_GRADIENT_MOMENT2_IMAGE,
	EMDL_MLMODEL_SIGMA2_NOISE,
	EMDL_MLMODEL_SIGMA2_REF,
	EMDL_MLMODEL_SSNR_REF,
	EMDL_MLMODEL_TAU2_FUDGE_FACTOR,
	EMDL_MLMODEL_TAU2_REF,

	EMDL_OPTIMISER_ACCURACY_ROT,
	EMDL_OPTIMISER_ACCURACY_TRANS, // deprecated
	EMDL_OPTIMISER_ACCURACY_TRANS_ANGSTROM,
	EMDL_OPTIMISER_ADAPTIVE_FRACTION,
	EMDL_OPTIMISER_ADAPTIVE_OVERSAMPLING,
	EMDL_OPTIMISER_AUTO_LOCAL_HP_ORDER,
	EMDL_OPTIMISER_AVAILABLE_MEMORY,
	EMDL_OPTIMISER_BEST_RESOL_THUS_FAR,
	EMDL_OPTIMISER_CHANGES_OPTIMAL_OFFSETS,
	EMDL_OPTIMISER_CHANGES_OPTIMAL_ORIENTS,
	EMDL_OPTIMISER_CHANGES_OPTIMAL_CLASSES,
	EMDL_OPTIMISER_COARSE_SIZE,
	EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED,
	EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED,
	EMDL_OPTIMISER_DATA_STARFILE,
	EMDL_OPTIMISER_DO_AUTO_REFINE,
	EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES,
	EMDL_OPTIMISER_DO_CORRECT_CTF,
	EMDL_OPTIMISER_DO_CORRECT_MAGNIFICATION,
	EMDL_OPTIMISER_DO_CORRECT_NORM,
	EMDL_OPTIMISER_DO_CORRECT_SCALE,
	EMDL_OPTIMISER_DO_EXTERNAL_RECONSTRUCT,
	EMDL_OPTIMISER_DO_REALIGN_MOVIES,
	EMDL_OPTIMISER_DO_MAP,
	EMDL_OPTIMISER_DO_VMGD,
	EMDL_OPTIMISER_DO_STOCHASTIC_EM,
	EMDL_OPTIMISER_EXTERNAL_RECONS_DATA_REAL,
	EMDL_OPTIMISER_EXTERNAL_RECONS_DATA_IMAG,
	EMDL_OPTIMISER_EXTERNAL_RECONS_WEIGHT,
	EMDL_OPTIMISER_EXTERNAL_RECONS_RESULT,
	EMDL_OPTIMISER_EXTERNAL_RECONS_NEWSTAR,
	EMDL_OPTIMISER_FAST_SUBSETS,
	EMDL_OPTIMISER_SGD_INI_ITER,
	EMDL_OPTIMISER_SGD_FIN_ITER,
	EMDL_OPTIMISER_SGD_INBETWEEN_ITER,
	EMDL_OPTIMISER_SGD_INI_RESOL,
	EMDL_OPTIMISER_SGD_FIN_RESOL,
	EMDL_OPTIMISER_SGD_INI_SUBSET_SIZE,
	EMDL_OPTIMISER_SGD_FIN_SUBSET_SIZE,
	EMDL_OPTIMISER_SGD_MU,
	EMDL_OPTIMISER_SGD_SIGMA2FUDGE_INI,
	EMDL_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE,
	EMDL_OPTIMISER_SGD_SKIP_ANNNEAL,
	EMDL_OPTIMISER_SGD_SUBSET_SIZE,
	EMDL_OPTIMISER_SGD_WRITE_EVERY_SUBSET,
	EMDL_OPTIMISER_SGD_MAX_SUBSETS,
	EMDL_OPTIMISER_SGD_STEPSIZE,
	EMDL_OPTIMISER_SGD_DO_MOM1,
	EMDL_OPTIMISER_SGD_DO_MOM2,
	EMDL_OPTIMISER_DO_SOLVENT_FLATTEN,
	EMDL_OPTIMISER_DO_SOLVENT_FSC,
	EMDL_OPTIMISER_DO_SKIP_ALIGN,
	EMDL_OPTIMISER_DO_SKIP_ROTATE,
	EMDL_OPTIMISER_DO_SPLIT_RANDOM_HALVES,
	EMDL_OPTIMISER_DO_ZERO_MASK,
	EMDL_OPTIMISER_FIX_SIGMA_NOISE,
	EMDL_OPTIMISER_FIX_SIGMA_OFFSET,
	EMDL_OPTIMISER_FIX_TAU,
	EMDL_OPTIMISER_HAS_CONVERGED,
	EMDL_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT,
	EMDL_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO,
	EMDL_OPTIMISER_DO_HELICAL_REFINE,
	EMDL_OPTIMISER_IGNORE_HELICAL_SYMMETRY,
	EMDL_OPTIMISER_FOURIER_MASK,
	EMDL_OPTIMISER_HELICAL_TWIST_INITIAL,
	EMDL_OPTIMISER_HELICAL_RISE_INITIAL,
	EMDL_OPTIMISER_HELICAL_Z_PERCENTAGE,
	EMDL_OPTIMISER_HELICAL_NSTART,
	EMDL_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER,
	EMDL_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER,
	EMDL_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT,
	EMDL_OPTIMISER_HELICAL_SIGMA_DISTANCE,
	EMDL_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED,
	EMDL_OPTIMISER_LOWRES_LIMIT_EXP,
	EMDL_OPTIMISER_HIGHRES_LIMIT_EXP,
	EMDL_OPTIMISER_HIGHRES_LIMIT_SGD,
	EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK,
	EMDL_OPTIMISER_INCR_SIZE,
	EMDL_OPTIMISER_ITERATION_NO,
	EMDL_OPTIMISER_LOCAL_SYMMETRY_FILENAME,
	EMDL_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES,
	EMDL_OPTIMISER_MAGNIFICATION_RANGE,
	EMDL_OPTIMISER_MAGNIFICATION_STEP,
	EMDL_OPTIMISER_MAX_COARSE_SIZE,
	EMDL_OPTIMISER_MAX_NR_POOL,
	EMDL_OPTIMISER_MODEL_STARFILE,
	EMDL_OPTIMISER_MODEL_STARFILE2,
	EMDL_OPTIMISER_NR_ITERATIONS,
	EMDL_OPTIMISER_NR_ITER_WO_RESOL_GAIN,
	EMDL_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES,
	EMDL_OPTIMISER_OPTICS_STARFILE,
	EMDL_OPTIMISER_OUTPUT_ROOTNAME,
	EMDL_OPTIMISER_PARTICLE_DIAMETER,
	EMDL_OPTIMISER_RADIUS_MASK_3D_MAP,
	EMDL_OPTIMISER_RADIUS_MASK_EXP_PARTICLES,
	EMDL_OPTIMISER_RANDOM_SEED,
	EMDL_OPTIMISER_REFS_ARE_CTF_CORRECTED,
	EMDL_OPTIMISER_SAMPLING_STARFILE,
	EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES,
	EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS,
	EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS,
	EMDL_OPTIMISER_SOLVENT_MASK_NAME,
	EMDL_OPTIMISER_SOLVENT_MASK2_NAME,
	EMDL_OPTIMISER_TAU_SPECTRUM_NAME,
	EMDL_OPTIMISER_USE_TOO_COARSE_SAMPLING,
	EMDL_OPTIMISER_WIDTH_MASK_EDGE,

	EMDL_ORIENT_FLIP,
	EMDL_ORIENT_ID,
	EMDL_ORIENT_ORIGIN_X, // (deprecated)
	EMDL_ORIENT_ORIGIN_Y, // (deprecated)
	EMDL_ORIENT_ORIGIN_Z, // (deprecated)
	EMDL_ORIENT_ORIGIN_X_PRIOR, // (deprecated)
	EMDL_ORIENT_ORIGIN_Y_PRIOR, // (deprecated)
	EMDL_ORIENT_ORIGIN_Z_PRIOR, // (deprecated)
	EMDL_ORIENT_ORIGIN_X_ANGSTROM,
	EMDL_ORIENT_ORIGIN_Y_ANGSTROM,
	EMDL_ORIENT_ORIGIN_Z_ANGSTROM,
	EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM,
	EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM,
	EMDL_ORIENT_ORIGIN_Z_PRIOR_ANGSTROM,
	EMDL_ORIENT_ROT,
	EMDL_ORIENT_ROT_PRIOR,
	EMDL_ORIENT_ROT_PRIOR_FLIP_RATIO,	// KThurber
	EMDL_ORIENT_TILT,
	EMDL_ORIENT_TILT_PRIOR,
	EMDL_ORIENT_PSI,
	EMDL_ORIENT_PSI_PRIOR,
	EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO,
	EMDL_ORIENT_PSI_PRIOR_FLIP,  // KThurber

	EMDL_PARTICLE_AUTOPICK_FOM,
	EMDL_PARTICLE_HELICAL_TUBE_ID,
	EMDL_PARTICLE_HELICAL_TUBE_PITCH,
	EMDL_PARTICLE_HELICAL_TRACK_LENGTH, //deprecated
	EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM,
	EMDL_PARTICLE_CLASS,
	EMDL_PARTICLE_DLL,
	EMDL_PARTICLE_ID,
	EMDL_PARTICLE_FOM,
	EMDL_PARTICLE_KL_DIVERGENCE,
	EMDL_PARTICLE_RANDOM_SUBSET,
	EMDL_PARTICLE_BEAM_TILT_CLASS,
	EMDL_PARTICLE_NAME,
	EMDL_PARTICLE_ORI_NAME,
	EMDL_PARTICLE_NR_SIGNIFICANT_SAMPLES,
	EMDL_PARTICLE_NR_FRAMES,
	EMDL_PARTICLE_NR_FRAMES_AVG,
	EMDL_PARTICLE_MOVIE_RUNNING_AVG,
	EMDL_PARTICLE_PMAX,
	EMDL_PARTICLE_NUMBER,

	EMDL_PIPELINE_JOB_COUNTER,
	EMDL_PIPELINE_NODE_NAME,
	EMDL_PIPELINE_NODE_TYPE,
	EMDL_PIPELINE_PROCESS_ALIAS,
	EMDL_PIPELINE_PROCESS_NAME,
	EMDL_PIPELINE_PROCESS_TYPE,
	EMDL_PIPELINE_PROCESS_STATUS,
	EMDL_PIPELINE_EDGE_FROM,
	EMDL_PIPELINE_EDGE_TO,
	EMDL_PIPELINE_EDGE_PROCESS,

	EMDL_POSTPROCESS_BFACTOR,
	EMDL_POSTPROCESS_FINAL_RESOLUTION,
	EMDL_POSTPROCESS_FRACTION_MOLWEIGHT,
	EMDL_POSTPROCESS_FRACTION_SOLVENT_MASK,
	EMDL_POSTPROCESS_FSC_GENERAL,
	EMDL_POSTPROCESS_FSC_TRUE,
	EMDL_POSTPROCESS_FSC_PART_MOLWEIGHT,
	EMDL_POSTPROCESS_FSC_PART_FRACMASK,
	EMDL_POSTPROCESS_FSC_MASKED,
	EMDL_POSTPROCESS_FSC_UNMASKED,
	EMDL_POSTPROCESS_FSC_RANDOM_MASKED,
	EMDL_POSTPROCESS_AMPLCORR_MASKED,
	EMDL_POSTPROCESS_AMPLCORR_UNMASKED,
	EMDL_POSTPROCESS_DPR_MASKED,
	EMDL_POSTPROCESS_DPR_UNMASKED,
	EMDL_POSTPROCESS_GUINIER_FIT_CORRELATION,
	EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT,
	EMDL_POSTPROCESS_GUINIER_FIT_SLOPE,
	EMDL_POSTPROCESS_GUINIER_VALUE_IN,
	EMDL_POSTPROCESS_GUINIER_VALUE_INVMTF,
	EMDL_POSTPROCESS_GUINIER_VALUE_WEIGHTED,
	EMDL_POSTPROCESS_GUINIER_VALUE_SHARPENED,
	EMDL_POSTPROCESS_GUINIER_VALUE_INTERCEPT,
	EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED,
	EMDL_POSTPROCESS_MOLWEIGHT,
	EMDL_POSTPROCESS_MTF_VALUE, ///< Detector MTF value
	EMDL_POSTPROCESS_RANDOMISE_FROM,
	EMDL_POSTPROCESS_UNFIL_HALFMAP1,
	EMDL_POSTPROCESS_UNFIL_HALFMAP2,

	EMDL_SAMPLING_IS_3D,
	EMDL_SAMPLING_IS_3D_TRANS,
	EMDL_SAMPLING_HEALPIX_ORDER,
	EMDL_SAMPLING_HEALPIX_ORDER_ORI,
	EMDL_SAMPLING_LIMIT_TILT,
	EMDL_SAMPLING_OFFSET_RANGE,
	EMDL_SAMPLING_OFFSET_STEP,
	EMDL_SAMPLING_OFFSET_RANGE_ORI,
	EMDL_SAMPLING_OFFSET_STEP_ORI,
	EMDL_SAMPLING_HELICAL_OFFSET_STEP,
	EMDL_SAMPLING_PERTURB,
	EMDL_SAMPLING_PERTURBATION_FACTOR,
	EMDL_SAMPLING_PRIOR_MODE,
	EMDL_SAMPLING_PSI_STEP,
	EMDL_SAMPLING_PSI_STEP_ORI,
	EMDL_SAMPLING_SIGMA_ROT,
	EMDL_SAMPLING_SIGMA_TILT,
	EMDL_SAMPLING_SIGMA_PSI,
	EMDL_SAMPLING_SYMMETRY,

	EMDL_SCHEDULE_EDGE_NUMBER,
	EMDL_SCHEDULE_EDGE_INPUT,
	EMDL_SCHEDULE_EDGE_OUTPUT,
	EMDL_SCHEDULE_EDGE_IS_FORK,
	EMDL_SCHEDULE_EDGE_OUTPUT_TRUE,
	EMDL_SCHEDULE_EDGE_BOOLEAN,
	EMDL_SCHEDULE_GENERAL_CURRENT_NODE,
	EMDL_SCHEDULE_GENERAL_ORIGINAL_START_NODE,
	EMDL_SCHEDULE_GENERAL_EMAIL,
	EMDL_SCHEDULE_GENERAL_NAME,
	EMDL_SCHEDULE_JOB_NAME,
	EMDL_SCHEDULE_JOB_ORI_NAME,
	EMDL_SCHEDULE_JOB_MODE,
	EMDL_SCHEDULE_JOB_HAS_STARTED,
	EMDL_SCHEDULE_OPERATOR_NAME,
	EMDL_SCHEDULE_OPERATOR_TYPE,
	EMDL_SCHEDULE_OPERATOR_INPUT1,
	EMDL_SCHEDULE_OPERATOR_INPUT2,
	EMDL_SCHEDULE_OPERATOR_OUTPUT,
	EMDL_SCHEDULE_VAR_BOOL_NAME,
	EMDL_SCHEDULE_VAR_BOOL_VALUE,
	EMDL_SCHEDULE_VAR_BOOL_ORI_VALUE,
	EMDL_SCHEDULE_VAR_FLOAT_NAME,
	EMDL_SCHEDULE_VAR_FLOAT_VALUE,
	EMDL_SCHEDULE_VAR_FLOAT_ORI_VALUE,
	EMDL_SCHEDULE_VAR_STRING_NAME,
	EMDL_SCHEDULE_VAR_STRING_VALUE,
	EMDL_SCHEDULE_VAR_STRING_ORI_VALUE,


	EMDL_SELECTED,
	EMDL_SELECT_PARTICLES_ZSCORE,
	EMDL_SORTED_IDX,
	EMDL_STARFILE_MOVIE_PARTICLES,
	EMDL_PERFRAME_CUMULATIVE_WEIGHT,
	EMDL_PERFRAME_RELATIVE_WEIGHT,

	EMDL_RESOLUTION,
	EMDL_RESOLUTION_ANGSTROM,
	EMDL_RESOLUTION_INVPIXEL,
	EMDL_SPECTRAL_IDX,

	EMDL_UNKNOWN_LABEL,

	EMDL_LAST_LABEL // **** NOTE ****: Do keep this label always at the end
	// it is here for looping purposes
};//close enum Label

enum EMDLabelType
{
	EMDL_INT, EMDL_BOOL, EMDL_DOUBLE, EMDL_STRING, EMDL_DOUBLE_VECTOR, EMDL_UNKNOWN
};

class EMDL
{
public:
	// This enum defines what MetaDataLabels this class can manage, if
	// you need a new one add it here and modify affected methods:
	//
	// - static EMDLabel codifyLabel( std::string strLabel ); EMDL::addLabel(EMDL_OPTIMISER_RANDOM_SEED, EMDL_INT, "randomSeed");

	// - static std::string EMDL::label2Str( EMDLabel inputLabel );
	// - void writeValuesToFile( std::ofstream &outfile, EMDLabel inputLabel );
	// - void addValue( std::string name, std::string value );
	//
	// Keep this special structure (using EMDL_FIRSTLABEL and EMDL_LAST_LABEL) so the
	// programmer can iterate through it like this:
	//
	//	for( EMDLabel mdl = EMDL_FIRST_LABEL ; mdl < EMDL_LAST_LABEL ; EMDLabel( mdl+1 ) )
	//

	static EMDLabel str2Label(const std::string &labelName);
	static std::string label2Str(const EMDLabel &label);

	static bool isInt(const EMDLabel &label);
	static bool isBool(const EMDLabel &label);
	static bool isString(const EMDLabel &label);
	static bool isDouble(const EMDLabel &label);
	static bool isNumber(const EMDLabel &label);
	static bool isDoubleVector(const EMDLabel &label);
	static bool isVector(const EMDLabel &label);
	static bool isUnknown(const EMDLabel &label);

	static bool isValidLabel(const EMDLabel &label);
	static bool isValidLabel(const std::string &labelName);

	static void printDefinitions(std::ostream& out);

private:
	static std::map<EMDLabel, EMDLabelData> data;
	static std::map<std::string, EMDLabel> names;
	static std::map<std::string, std::string> definitions;
	static StaticInitialization initialization; //Just for initialization

	static void addLabel(EMDLabel label, EMDLabelType type, std::string name, std::string definition = "undocumented");
	static void addAltLabel(EMDLabel label, std::string name);

	friend class StaticInitialization;
}
;//close class MLD definition

//Just an struct to store type and string alias
class EMDLabelData
{
public:
	EMDLabelType type;
	std::string str;
	//Default constructor
	EMDLabelData()
	{
	}
	EMDLabelData(EMDLabelType t, std::string s)
	{
		type = t;
		str = s;
	}
};//close class EMDLabelData c

//Just a class for static initialization
class StaticInitialization
{
private:
	StaticInitialization()
	{
		///==== Add labels entries from here in the SAME ORDER as declared in ENUM ==========
		EMDL::addLabel(EMDL_COMMENT, EMDL_STRING, "rlnComment", "A metadata comment (This is treated in a special way)");

		EMDL::addLabel(EMDL_AREA_ID, EMDL_INT, "rlnAreaId", "ID (i.e. a unique number) of an area (i.e. field-of-view)");
		EMDL::addLabel(EMDL_AREA_NAME, EMDL_STRING, "rlnAreaName", "Name of an area (i.e. field-of-view)");

		EMDL::addLabel(EMDL_BODY_MASK_NAME, EMDL_STRING, "rlnBodyMaskName", "Name of an image that contains a [0,1] body mask for multi-body refinement");
		EMDL::addLabel(EMDL_BODY_KEEP_FIXED, EMDL_INT, "rlnBodyKeepFixed", "Flag to indicate whether to keep a body fixed (value 1) or keep on refining it (0)");
		EMDL::addLabel(EMDL_BODY_REFERENCE_NAME, EMDL_STRING, "rlnBodyReferenceName", "Name of an image that contains the initial reference for one body of a multi-body refinement");
		EMDL::addLabel(EMDL_BODY_ROTATE_DIRECTION_X, EMDL_DOUBLE, "rlnBodyRotateDirectionX", "X-component of axis around which to rotate this body");
		EMDL::addLabel(EMDL_BODY_ROTATE_DIRECTION_Y, EMDL_DOUBLE, "rlnBodyRotateDirectionY", "Y-component of axis around which to rotate this body");
		EMDL::addLabel(EMDL_BODY_ROTATE_DIRECTION_Z, EMDL_DOUBLE, "rlnBodyRotateDirectionZ", "Z-component of axis around which to rotate this body");
		EMDL::addLabel(EMDL_BODY_ROTATE_RELATIVE_TO, EMDL_INT, "rlnBodyRotateRelativeTo", "Number of the body relative to which this body rotates (if negative, use rlnBodyRotateDirectionXYZ)");
		EMDL::addLabel(EMDL_BODY_SIGMA_ANG, EMDL_DOUBLE, "rlnBodySigmaAngles", "Width of prior on all three Euler angles of a body in multibody refinement (in degrees)");
		EMDL::addLabel(EMDL_BODY_SIGMA_OFFSET, EMDL_DOUBLE, "rlnBodySigmaOffset", "Width of prior on origin offsets of a body in multibody refinement (in pixels)");
		EMDL::addLabel(EMDL_BODY_SIGMA_OFFSET_ANGSTROM, EMDL_DOUBLE, "rlnBodySigmaOffsetAngst", "Width of prior on origin offsets of a body in multibody refinement (in Angstroms)");
		EMDL::addLabel(EMDL_BODY_SIGMA_ROT, EMDL_DOUBLE, "rlnBodySigmaRot", "Width of prior on rot angles of a body in multibody refinement (in degrees)");
		EMDL::addLabel(EMDL_BODY_SIGMA_TILT, EMDL_DOUBLE, "rlnBodySigmaTilt", "Width of prior on tilt angles of a body in multibody refinement (in degrees)");
		EMDL::addLabel(EMDL_BODY_SIGMA_PSI, EMDL_DOUBLE, "rlnBodySigmaPsi", "Width of prior on psi angles of a body in multibody refinement (in degrees)");
		EMDL::addLabel(EMDL_BODY_STAR_FILE, EMDL_STRING, "rlnBodyStarFile", "Name of STAR file with body masks and metadata");

		EMDL::addLabel(EMDL_CTF_ASTIGMATISM, EMDL_DOUBLE, "rlnCtfAstigmatism", "Absolute value of the difference between defocus in U- and V-direction (in A)");
		EMDL::addLabel(EMDL_CTF_BFACTOR, EMDL_DOUBLE, "rlnCtfBfactor", "B-factor (in A^2) that describes CTF power spectrum fall-off");
		EMDL::addLabel(EMDL_CTF_MAXRES, EMDL_DOUBLE, "rlnCtfMaxResolution", "Estimated maximum resolution (in A) of significant CTF Thon rings");
		EMDL::addLabel(EMDL_CTF_VALIDATIONSCORE, EMDL_DOUBLE, "rlnCtfValidationScore", "Gctf-based validation score for the quality of the CTF fit");
		EMDL::addLabel(EMDL_CTF_SCALEFACTOR, EMDL_DOUBLE, "rlnCtfScalefactor", "Linear scale-factor on the CTF (values between 0 and 1)");
		EMDL::addLabel(EMDL_CTF_VOLTAGE, EMDL_DOUBLE, "rlnVoltage", "Voltage of the microscope (in kV)");
		EMDL::addLabel(EMDL_CTF_DEFOCUSU, EMDL_DOUBLE, "rlnDefocusU", "Defocus in U-direction (in Angstroms, positive values for underfocus)");
		EMDL::addLabel(EMDL_CTF_DEFOCUSV, EMDL_DOUBLE, "rlnDefocusV", "Defocus in V-direction (in Angstroms, positive values for underfocus)");
		EMDL::addLabel(EMDL_CTF_DEFOCUS_ANGLE, EMDL_DOUBLE, "rlnDefocusAngle", "Angle between X and defocus U direction (in degrees)");
		EMDL::addLabel(EMDL_CTF_CS, EMDL_DOUBLE, "rlnSphericalAberration", "Spherical aberration (in millimeters)");
		EMDL::addLabel(EMDL_CTF_CA, EMDL_DOUBLE, "rlnChromaticAberration", "Chromatic aberration (in millimeters)");
		EMDL::addLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE, EMDL_DOUBLE, "rlnDetectorPixelSize", "Pixel size of the detector (in micrometers)");
		EMDL::addLabel(EMDL_CTF_POWER_SPECTRUM, EMDL_STRING, "rlnCtfPowerSpectrum", "Power spectrum for CTF estimation");
		EMDL::addLabel(EMDL_CTF_ENERGY_LOSS, EMDL_DOUBLE, "rlnEnergyLoss", "Energy loss (in eV)");
		EMDL::addLabel(EMDL_CTF_FOM, EMDL_DOUBLE, "rlnCtfFigureOfMerit", "Figure of merit for the fit of the CTF (not used inside relion_refine)");
		EMDL::addLabel(EMDL_CTF_IMAGE, EMDL_STRING, "rlnCtfImage", "Name of an image with all CTF values");
		EMDL::addLabel(EMDL_CTF_LENS_STABILITY, EMDL_DOUBLE, "rlnLensStability", "Lens stability (in ppm)");
		EMDL::addLabel(EMDL_CTF_MAGNIFICATION, EMDL_DOUBLE, "rlnMagnification", "Magnification at the detector (in times)");
		EMDL::addLabel(EMDL_CTF_PHASESHIFT, EMDL_DOUBLE, "rlnPhaseShift", "Phase-shift from a phase-plate (in degrees)");
		EMDL::addLabel(EMDL_CTF_CONVERGENCE_CONE, EMDL_DOUBLE, "rlnConvergenceCone", "Convergence cone (in mrad)");
		EMDL::addLabel(EMDL_CTF_LONGITUDINAL_DISPLACEMENT, EMDL_DOUBLE, "rlnLongitudinalDisplacement", "Longitudinal displacement (in Angstroms)");
		EMDL::addLabel(EMDL_CTF_TRANSVERSAL_DISPLACEMENT, EMDL_DOUBLE, "rlnTransversalDisplacement", "Transversal displacement (in Angstroms)");
		EMDL::addLabel(EMDL_CTF_Q0, EMDL_DOUBLE, "rlnAmplitudeContrast", "Amplitude contrast (as a fraction, i.e. 10% = 0.1)");
		EMDL::addLabel(EMDL_CTF_VALUE, EMDL_DOUBLE, "rlnCtfValue", "Value of the Contrast Transfer Function");


		EMDL::addLabel(EMDL_IMAGE_NAME, EMDL_STRING, "rlnImageName", "Name of an image");
		EMDL::addLabel(EMDL_IMAGE_ORI_NAME, EMDL_STRING, "rlnImageOriginalName", "Original name of an image");
		EMDL::addLabel(EMDL_IMAGE_RECONSTRUCT_NAME, EMDL_STRING, "rlnReconstructImageName", "Name of an image to be used for reconstruction only");
		EMDL::addLabel(EMDL_IMAGE_ID, EMDL_INT, "rlnImageId", "ID (i.e. a unique number) of an image");
		EMDL::addLabel(EMDL_IMAGE_ENABLED, EMDL_BOOL, "rlnEnabled", "Not used in RELION, only included for backward compatibility with XMIPP selfiles");
		EMDL::addLabel(EMDL_IMAGE_DATATYPE, EMDL_INT, "rlnDataType", "Type of data stored in an image (e.g. int, RFLOAT etc)");
		EMDL::addLabel(EMDL_IMAGE_DIMENSIONALITY, EMDL_INT, "rlnImageDimensionality", "Dimensionality of data stored in an image (i.e. 2 or 3)");
		EMDL::addLabel(EMDL_IMAGE_BEAMTILT_X, EMDL_DOUBLE, "rlnBeamTiltX", "Beam tilt in the X-direction (in mrad)");
		EMDL::addLabel(EMDL_IMAGE_BEAMTILT_Y, EMDL_DOUBLE, "rlnBeamTiltY", "Beam tilt in the Y-direction (in mrad)");
		EMDL::addLabel(EMDL_IMAGE_MTF_FILENAME, EMDL_STRING, "rlnMtfFileName", "The filename of a STAR file with the MTF for this optics group or image");
		EMDL::addLabel(EMDL_IMAGE_OPTICS_GROUP, EMDL_INT, "rlnOpticsGroup", "Group of particles with identical optical properties");
		EMDL::addLabel(EMDL_IMAGE_OPTICS_GROUP_NAME, EMDL_STRING, "rlnOpticsGroupName", "The name of a group of particles with identical optical properties");
		EMDL::addLabel(EMDL_IMAGE_ODD_ZERNIKE_COEFFS, EMDL_DOUBLE_VECTOR, "rlnOddZernike", "Coefficients for the antisymmetrical Zernike polynomials");
		EMDL::addLabel(EMDL_IMAGE_EVEN_ZERNIKE_COEFFS, EMDL_DOUBLE_VECTOR, "rlnEvenZernike", "Coefficients for the symmetrical Zernike polynomials");
		EMDL::addLabel(EMDL_IMAGE_PIXEL_SIZE, EMDL_DOUBLE, "rlnImagePixelSize", "Pixel size (in Angstrom)");
		EMDL::addLabel(EMDL_IMAGE_MAG_MATRIX_00, EMDL_DOUBLE, "rlnMagMat00", "Anisotropic magnification matrix, element 1,1");
		EMDL::addLabel(EMDL_IMAGE_MAG_MATRIX_01, EMDL_DOUBLE, "rlnMagMat01", "Anisotropic magnification matrix, element 1,2");
		EMDL::addLabel(EMDL_IMAGE_MAG_MATRIX_10, EMDL_DOUBLE, "rlnMagMat10", "Anisotropic magnification matrix, element 2,1");
		EMDL::addLabel(EMDL_IMAGE_MAG_MATRIX_11, EMDL_DOUBLE, "rlnMagMat11", "Anisotropic magnification matrix, element 2,2");

		EMDL::addLabel(EMDL_IMAGE_COORD_X, EMDL_DOUBLE, "rlnCoordinateX", "X-Position of an image in a micrograph (in pixels)");
		EMDL::addLabel(EMDL_IMAGE_COORD_Y, EMDL_DOUBLE, "rlnCoordinateY", "Y-Position of an image in a micrograph (in pixels)");
		EMDL::addLabel(EMDL_IMAGE_COORD_Z, EMDL_DOUBLE, "rlnCoordinateZ", "Z-Position of an image in a 3D micrograph, i.e. tomogram (in pixels)");
		EMDL::addLabel(EMDL_IMAGE_FRAME_NR, EMDL_INT, "rlnMovieFrameNumber", "Number of a movie frame");
		EMDL::addLabel(EMDL_IMAGE_NORM_CORRECTION, EMDL_DOUBLE, "rlnNormCorrection", "Normalisation correction value for an image");
		EMDL::addLabel(EMDL_IMAGE_MAGNIFICATION_CORRECTION, EMDL_DOUBLE, "rlnMagnificationCorrection", "Magnification correction value for an image");
		EMDL::addLabel(EMDL_IMAGE_SAMPLINGRATE, EMDL_DOUBLE, "rlnSamplingRate", "Sampling rate of an image (in Angstrom/pixel)");
		EMDL::addLabel(EMDL_IMAGE_SAMPLINGRATE_X, EMDL_DOUBLE, "rlnSamplingRateX", "Sampling rate in X-direction of an image (in Angstrom/pixel)");
		EMDL::addLabel(EMDL_IMAGE_SAMPLINGRATE_Y, EMDL_DOUBLE, "rlnSamplingRateY", "Sampling rate in Y-direction of an image (in Angstrom/pixel)");
		EMDL::addLabel(EMDL_IMAGE_SAMPLINGRATE_Z, EMDL_DOUBLE, "rlnSamplingRateZ", "Sampling rate in Z-direction of an image (in Angstrom/pixel)");
		EMDL::addLabel(EMDL_IMAGE_SIZE, EMDL_INT, "rlnImageSize", "Size of an image (in pixels)");
		EMDL::addLabel(EMDL_IMAGE_SIZE_X, EMDL_INT, "rlnImageSizeX", "Size of an image in the X-direction (in pixels)");
		EMDL::addLabel(EMDL_IMAGE_SIZE_Y, EMDL_INT, "rlnImageSizeY", "Size of an image in the Y-direction (in pixels)");
		EMDL::addLabel(EMDL_IMAGE_SIZE_Z, EMDL_INT, "rlnImageSizeZ", "Size of an image in the Z-direction (in pixels)");
		EMDL::addLabel(EMDL_IMAGE_STATS_MIN, EMDL_DOUBLE, "rlnMinimumValue", "Minimum value for the pixels in an image");
		EMDL::addLabel(EMDL_IMAGE_STATS_MAX, EMDL_DOUBLE, "rlnMaximumValue", "Maximum value for the pixels in an image");
		EMDL::addLabel(EMDL_IMAGE_STATS_AVG, EMDL_DOUBLE, "rlnAverageValue", "Average value for the pixels in an image");
		EMDL::addLabel(EMDL_IMAGE_STATS_STDDEV, EMDL_DOUBLE, "rlnStandardDeviationValue", "Standard deviation for the pixel values in an image");
		EMDL::addLabel(EMDL_IMAGE_STATS_SKEW, EMDL_DOUBLE, "rlnSkewnessValue", "Skewness (3rd moment) for the pixel values in an image");
		EMDL::addLabel(EMDL_IMAGE_STATS_KURT, EMDL_DOUBLE, "rlnKurtosisExcessValue", "Kurtosis excess (4th moment - 3) for the pixel values in an image");
		EMDL::addLabel(EMDL_IMAGE_WEIGHT, EMDL_DOUBLE, "rlnImageWeight", "Relative weight of an image");

		EMDL::addLabel(EMDL_MASK_NAME, EMDL_STRING, "rlnMaskName", "Name of an image that contains a [0,1] mask");

		EMDL::addLabel(EMDL_JOB_IS_CONTINUE, EMDL_BOOL, "rlnJobIsContinue", "Is tthis a continuation job?");
		EMDL::addLabel(EMDL_JOB_TYPE, EMDL_INT, "rlnJobType", "Which type of job is this?");
		EMDL::addLabel(EMDL_JOB_TYPE_NAME, EMDL_STRING, "rlnJobTypeName", "The name for this type of job (also name of main directory for output jobs)");

		EMDL::addLabel(EMDL_JOBOPTION_TYPE, EMDL_INT, "rlnJoboptionType", "Which type of joboption is this?");
		EMDL::addLabel(EMDL_JOBOPTION_VARIABLE, EMDL_STRING, "rlnJobOptionVariable", "Name of the joboption variable");
		EMDL::addLabel(EMDL_JOBOPTION_VALUE, EMDL_STRING, "rlnJobOptionValue", "Value of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_LABEL, EMDL_STRING, "rlnJobOptionGUILabel", "GUI label of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_DEFAULT_VALUE, EMDL_STRING, "rlnJobOptionDefaultValue", "Default value of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_MINVAL, EMDL_DOUBLE, "rlnJobOptionSliderMin", "Minimum value for slider of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_MAXVAL, EMDL_DOUBLE, "rlnJobOptionSliderMax", "Maximum value for slider of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_STEPVAL, EMDL_DOUBLE, "rlnJobOptionSliderStep", "Step value for slider of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_HELPTEXT, EMDL_STRING, "rlnJobOptionHelpText", "Extra helptext of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_PATTERN, EMDL_STRING, "rlnJobOptionFilePattern", "Pattern for file browser of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_DIRECTORY, EMDL_STRING, "rlnJobOptionDirectoryDefault", "Default directory for file browser of a joboption");
		EMDL::addLabel(EMDL_JOBOPTION_MENUOPTIONS, EMDL_STRING, "rlnJobOptionMenuOptions", "Options for pull-down menu");

		EMDL::addLabel(EMDL_MATRIX_1_1, EMDL_DOUBLE, "rlnMatrix_1_1", "Matrix element (1,1) of a 3x3 matrix");
		EMDL::addLabel(EMDL_MATRIX_1_2, EMDL_DOUBLE, "rlnMatrix_1_2", "Matrix element (1,2) of a 3x3 matrix");
		EMDL::addLabel(EMDL_MATRIX_1_3, EMDL_DOUBLE, "rlnMatrix_1_3", "Matrix element (1,3) of a 3x3 matrix");
		EMDL::addLabel(EMDL_MATRIX_2_1, EMDL_DOUBLE, "rlnMatrix_2_1", "Matrix element (2,1) of a 3x3 matrix");
		EMDL::addLabel(EMDL_MATRIX_2_2, EMDL_DOUBLE, "rlnMatrix_2_2", "Matrix element (2,1) of a 3x3 matrix");
		EMDL::addLabel(EMDL_MATRIX_2_3, EMDL_DOUBLE, "rlnMatrix_2_3", "Matrix element (2,1) of a 3x3 matrix");
		EMDL::addLabel(EMDL_MATRIX_3_1, EMDL_DOUBLE, "rlnMatrix_3_1", "Matrix element (3,1) of a 3x3 matrix");
		EMDL::addLabel(EMDL_MATRIX_3_2, EMDL_DOUBLE, "rlnMatrix_3_2", "Matrix element (3,1) of a 3x3 matrix");
		EMDL::addLabel(EMDL_MATRIX_3_3, EMDL_DOUBLE, "rlnMatrix_3_3", "Matrix element (3,1) of a 3x3 matrix");

		EMDL::addLabel(EMDL_MICROGRAPH_ACCUM_MOTION_TOTAL, EMDL_DOUBLE, "rlnAccumMotionTotal","Accumulated global motion during the entire movie (in A)");
		EMDL::addLabel(EMDL_MICROGRAPH_ACCUM_MOTION_EARLY, EMDL_DOUBLE, "rlnAccumMotionEarly","Accumulated global motion during the first frames of the movie (in A)");
		EMDL::addLabel(EMDL_MICROGRAPH_ACCUM_MOTION_LATE, EMDL_DOUBLE, "rlnAccumMotionLate","Accumulated global motion during the last frames of the movie (in A)");
		EMDL::addLabel(EMDL_MICROGRAPH_ID, EMDL_INT, "rlnMicrographId", "ID (i.e. a unique number) of a micrograph");
		EMDL::addLabel(EMDL_MICROGRAPH_NAME, EMDL_STRING, "rlnMicrographName", "Name of a micrograph");
		EMDL::addLabel(EMDL_MICROGRAPH_GAIN_NAME, EMDL_STRING, "rlnMicrographGainName", "Name of a gain reference");
		EMDL::addLabel(EMDL_MICROGRAPH_DEFECT_FILE, EMDL_STRING, "rlnMicrographDefectFile", "Name of a defect list file");
		EMDL::addLabel(EMDL_MICROGRAPH_NAME_WODOSE, EMDL_STRING, "rlnMicrographNameNoDW", "Name of a micrograph without dose weighting");
		EMDL::addLabel(EMDL_MICROGRAPH_MOVIE_NAME, EMDL_STRING, "rlnMicrographMovieName", "Name of a micrograph movie stack");
		EMDL::addLabel(EMDL_MICROGRAPH_METADATA_NAME, EMDL_STRING, "rlnMicrographMetadata", "Name of a micrograph metadata file");
		EMDL::addLabel(EMDL_MICROGRAPH_TILT_ANGLE, EMDL_DOUBLE, "rlnMicrographTiltAngle", "Tilt angle (in degrees) used to collect a micrograph");
		EMDL::addLabel(EMDL_MICROGRAPH_TILT_AXIS_DIRECTION, EMDL_DOUBLE, "rlnMicrographTiltAxisDirection", "Direction of the tilt-axis (in degrees) used to collect a micrograph");
		EMDL::addLabel(EMDL_MICROGRAPH_TILT_AXIS_OUTOFPLANE, EMDL_DOUBLE, "rlnMicrographTiltAxisOutOfPlane", "Out-of-plane angle (in degrees) of the tilt-axis used to collect a micrograph (90=in-plane)");
		EMDL::addLabel(EMDL_MICROGRAPH_ORIGINAL_PIXEL_SIZE, EMDL_DOUBLE, "rlnMicrographOriginalPixelSize", "Pixel size of original movie before binning in Angstrom/pixel.");
		EMDL::addLabel(EMDL_MICROGRAPH_PIXEL_SIZE, EMDL_DOUBLE, "rlnMicrographPixelSize", "Pixel size of (averaged) micrographs after binning in Angstrom/pixel.");
		EMDL::addLabel(EMDL_MICROGRAPH_PRE_EXPOSURE, EMDL_DOUBLE, "rlnMicrographPreExposure", "Pre-exposure dose in electrons per square Angstrom");
		EMDL::addLabel(EMDL_MICROGRAPH_DOSE_RATE, EMDL_DOUBLE, "rlnMicrographDoseRate", "Dose rate in electrons per square Angstrom per frame");
		EMDL::addLabel(EMDL_MICROGRAPH_BINNING, EMDL_DOUBLE, "rlnMicrographBinning", "Micrograph binning factor");
		EMDL::addLabel(EMDL_MICROGRAPH_FRAME_NUMBER, EMDL_INT, "rlnMicrographFrameNumber", "Micrograph frame number");
		EMDL::addLabel(EMDL_MICROGRAPH_MOTION_MODEL_VERSION, EMDL_INT, "rlnMotionModelVersion", "Version of micrograph motion model");
		EMDL::addLabel(EMDL_MICROGRAPH_START_FRAME, EMDL_INT, "rlnMicrographStartFrame", "Start frame of a motion model");
		EMDL::addLabel(EMDL_MICROGRAPH_END_FRAME, EMDL_INT, "rlnMicrographEndFrame", "End frame of a motion model");
		EMDL::addLabel(EMDL_MICROGRAPH_SHIFT_X, EMDL_DOUBLE, "rlnMicrographShiftX", "X shift of a (patch of) micrograph");
		EMDL::addLabel(EMDL_MICROGRAPH_SHIFT_Y, EMDL_DOUBLE, "rlnMicrographShiftY", "Y shift of a (patch of) micrograph");
		EMDL::addLabel(EMDL_MICROGRAPH_MOTION_COEFFS_IDX, EMDL_INT, "rlnMotionModelCoeffsIdx", "Index of a coefficient of a motion model");
		EMDL::addLabel(EMDL_MICROGRAPH_MOTION_COEFF, EMDL_DOUBLE, "rlnMotionModelCoeff", "A coefficient of a motion model");

		EMDL::addLabel(EMDL_MLMODEL_ACCURACY_ROT, EMDL_DOUBLE, "rlnAccuracyRotations", "Estimated accuracy (in degrees) with which rotations can be assigned");
		EMDL::addLabel(EMDL_MLMODEL_ACCURACY_TRANS, EMDL_DOUBLE, "rlnAccuracyTranslations", "Estimated accuracy (in pixels) with which translations can be assigned");
		EMDL::addLabel(EMDL_MLMODEL_ACCURACY_TRANS_ANGSTROM, EMDL_DOUBLE, "rlnAccuracyTranslationsAngst", "Estimated accuracy (in Angstroms) with which translations can be assigned");
		EMDL::addLabel(EMDL_MLMODEL_AVE_PMAX, EMDL_DOUBLE, "rlnAveragePmax", "Average value (over all images) of the maxima of the probability distributions");
		EMDL::addLabel(EMDL_MLMODEL_CURRENT_RESOLUTION, EMDL_DOUBLE, "rlnCurrentResolution", "Current resolution where SSNR^MAP drops below 1 (in 1/Angstroms)");
		EMDL::addLabel(EMDL_MLMODEL_CURRENT_SIZE, EMDL_INT, "rlnCurrentImageSize", "Current size of the images used in the refinement");
		EMDL::addLabel(EMDL_MLMODEL_DATA_VS_PRIOR_REF, EMDL_DOUBLE, "rlnSsnrMap", "Spectral signal-to-noise ratio as defined for MAP estimation (SSNR^MAP)");
		EMDL::addLabel(EMDL_MLMODEL_DIMENSIONALITY, EMDL_INT, "rlnReferenceDimensionality", "Dimensionality of the references (2D/3D)");
		EMDL::addLabel(EMDL_MLMODEL_DIMENSIONALITY_DATA, EMDL_INT, "rlnDataDimensionality", "Dimensionality of the data (2D/3D)");
		EMDL::addLabel(EMDL_MLMODEL_DIFF2_HALVES_REF, EMDL_DOUBLE, "rlnDiff2RandomHalves", "Power of the differences between two independent reconstructions from random halves of the data");
		EMDL::addLabel(EMDL_MLMODEL_ESTIM_RESOL_REF, EMDL_DOUBLE, "rlnEstimatedResolution", "Estimated resolution (in A) for a reference");
		EMDL::addLabel(EMDL_MLMODEL_FOURIER_COVERAGE_REF, EMDL_DOUBLE, "rlnFourierCompleteness", "Fraction of Fourier components (per resolution shell) with SNR>1");
		EMDL::addLabel(EMDL_MLMODEL_FOURIER_COVERAGE_TOTAL_REF, EMDL_DOUBLE, "rlnOverallFourierCompleteness", "Fraction of all Fourier components up to the current resolution with SNR>1");
		EMDL::addLabel(EMDL_MLMODEL_FSC_HALVES_REF, EMDL_DOUBLE, "rlnGoldStandardFsc", "Fourier shell correlation between two independent reconstructions from random halves of the data");
		EMDL::addLabel(EMDL_MLMODEL_GROUP_NAME, EMDL_STRING, "rlnGroupName", "The name of a group of images (e.g. all images from a micrograph)");
		EMDL::addLabel(EMDL_MLMODEL_GROUP_NO, EMDL_INT, "rlnGroupNumber", "The number of a group of images");
		EMDL::addLabel(EMDL_MLMODEL_GROUP_NR_PARTICLES, EMDL_INT, "rlnGroupNrParticles", "Number particles in a group of images");
		EMDL::addLabel(EMDL_MLMODEL_GROUP_SCALE_CORRECTION, EMDL_DOUBLE, "rlnGroupScaleCorrection", "Intensity-scale correction for a group of images");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_NR_ASU, EMDL_INT, "rlnNrHelicalAsymUnits", "How many new helical asymmetric units are there in each box");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_TWIST, EMDL_DOUBLE, "rlnHelicalTwist", "The helical twist (rotation per subunit) in degrees");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_TWIST_MIN, EMDL_DOUBLE, "rlnHelicalTwistMin", "Minimum helical twist (in degrees, + for right-handedness)");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_TWIST_MAX, EMDL_DOUBLE, "rlnHelicalTwistMax", "Maximum helical twist (in degrees, + for right-handedness)");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_TWIST_INITIAL_STEP, EMDL_DOUBLE, "rlnHelicalTwistInitialStep", "Initial step of helical twist search (in degrees)");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_RISE, EMDL_DOUBLE, "rlnHelicalRise", "The helical rise (translation per subunit) in Angstroms");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_RISE_MIN, EMDL_DOUBLE, "rlnHelicalRiseMin", "Minimum helical rise (in Angstroms)");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_RISE_MAX, EMDL_DOUBLE, "rlnHelicalRiseMax", "Maximum helical rise (in Angstroms)");
		EMDL::addLabel(EMDL_MLMODEL_HELICAL_RISE_INITIAL_STEP, EMDL_DOUBLE, "rlnHelicalRiseInitialStep", "Initial step of helical rise search (in Angstroms)");
		EMDL::addLabel(EMDL_MLMODEL_IS_HELIX, EMDL_BOOL, "rlnIsHelix", "Flag to indicate that helical refinement should be performed");
		EMDL::addLabel(EMDL_MLMODEL_INTERPOLATOR, EMDL_INT, "rlnFourierSpaceInterpolator", "The kernel used for Fourier-space interpolation (NN=0, linear=1)");
		EMDL::addLabel(EMDL_MLMODEL_LL, EMDL_DOUBLE, "rlnLogLikelihood", "Value of the log-likelihood target function");
		EMDL::addLabel(EMDL_MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION,  EMDL_INT, "rlnMinRadiusNnInterpolation","Minimum radius for NN-interpolation (in Fourier pixels), for smaller radii linear int. is used");
		EMDL::addLabel(EMDL_MLMODEL_NORM_CORRECTION_AVG, EMDL_DOUBLE, "rlnNormCorrectionAverage", "Average value (over all images) of the normalisation correction values");
		EMDL::addLabel(EMDL_MLMODEL_NR_CLASSES, EMDL_INT, "rlnNrClasses", "The number of references (i.e. classes) to be used in refinement");
		EMDL::addLabel(EMDL_MLMODEL_NR_BODIES, EMDL_INT, "rlnNrBodies", "The number of independent rigid bodies to be refined in multi-body refinement");
		EMDL::addLabel(EMDL_MLMODEL_NR_GROUPS, EMDL_INT, "rlnNrGroups", "The number of different groups of images (each group has its own noise spectrum, and intensity-scale correction)");
		EMDL::addLabel(EMDL_MLMODEL_ORIENTABILITY_CONTRIBUTION, EMDL_DOUBLE, "rlnSpectralOrientabilityContribution", "Spectral SNR contribution to the orientability of individual particles");
		EMDL::addLabel(EMDL_MLMODEL_ORIGINAL_SIZE, EMDL_INT, "rlnOriginalImageSize", "Original size of the images (in pixels)");
		EMDL::addLabel(EMDL_MLMODEL_PADDING_FACTOR, EMDL_DOUBLE, "rlnPaddingFactor", "Oversampling factor for Fourier transforms of the references");
		EMDL::addLabel(EMDL_MLMODEL_PDF_CLASS, EMDL_DOUBLE, "rlnClassDistribution", "Probability Density Function of the different classes (i.e. fraction of images assigned to each class)");
		EMDL::addLabel(EMDL_MLMODEL_PRIOR_OFFX_CLASS, EMDL_DOUBLE, "rlnClassPriorOffsetX", "Prior in the X-offset for a class (in pixels)");
		EMDL::addLabel(EMDL_MLMODEL_PRIOR_OFFY_CLASS, EMDL_DOUBLE, "rlnClassPriorOffsetY", "Prior in the Y-offset for a class (in pixels)");
		EMDL::addLabel(EMDL_MLMODEL_PDF_ORIENT, EMDL_DOUBLE, "rlnOrientationDistribution", "Probability Density Function of the orientations  (i.e. fraction of images assigned to each orient)");
		EMDL::addLabel(EMDL_MLMODEL_PIXEL_SIZE, EMDL_DOUBLE, "rlnPixelSize", "Size of the pixels in the references and images (in Angstroms)");
		EMDL::addLabel(EMDL_MLMODEL_POWER_REF, EMDL_DOUBLE, "rlnReferenceSpectralPower", "Spherical average of the power of the reference");
		EMDL::addLabel(EMDL_MLMODEL_PRIOR_MODE, EMDL_INT, "rlnOrientationalPriorMode", "Mode for prior distributions on the orientations (0=no prior; 1=(rot,tilt,psi); 2=(rot,tilt); 3=rot; 4=tilt; 5=psi) ");
		EMDL::addLabel(EMDL_MLMODEL_REF_IMAGE, EMDL_STRING, "rlnReferenceImage", "Name of a reference image");
		EMDL::addLabel(EMDL_MLMODEL_GRADIENT_MOMENT1_IMAGE, EMDL_STRING, "rlnGradMoment1", "Name of image containing the first moment of the gradient");
		EMDL::addLabel(EMDL_MLMODEL_GRADIENT_MOMENT2_IMAGE, EMDL_STRING, "rlnGradMoment2", "Name of image containing the second moment of the gradient");
		EMDL::addLabel(EMDL_MLMODEL_SIGMA_OFFSET, EMDL_DOUBLE, "rlnSigmaOffsets","Standard deviation in the origin offsets (in pixels)");
		EMDL::addLabel(EMDL_MLMODEL_SIGMA_OFFSET_ANGSTROM, EMDL_DOUBLE, "rlnSigmaOffsetsAngst","Standard deviation in the origin offsets (in Angstroms)");
		EMDL::addLabel(EMDL_MLMODEL_SIGMA2_NOISE, EMDL_DOUBLE, "rlnSigma2Noise", "Spherical average of the standard deviation in the noise (sigma)");
		EMDL::addLabel(EMDL_MLMODEL_SIGMA2_REF, EMDL_DOUBLE, "rlnReferenceSigma2", "Spherical average of the estimated power in the noise of a reference");
		EMDL::addLabel(EMDL_MLMODEL_SIGMA_ROT, EMDL_DOUBLE, "rlnSigmaPriorRotAngle", "Standard deviation of the prior on the rot (i.e. first Euler) angle");
		EMDL::addLabel(EMDL_MLMODEL_SIGMA_TILT, EMDL_DOUBLE, "rlnSigmaPriorTiltAngle", "Standard deviation of the prior on the tilt (i.e. second Euler) angle");
		EMDL::addLabel(EMDL_MLMODEL_SIGMA_PSI, EMDL_DOUBLE, "rlnSigmaPriorPsiAngle", "Standard deviation of the prior on the psi (i.e. third Euler) angle");
		EMDL::addLabel(EMDL_MLMODEL_SSNR_REF, EMDL_DOUBLE, "rlnSignalToNoiseRatio", "Spectral signal-to-noise ratio for a reference");
		EMDL::addLabel(EMDL_MLMODEL_TAU2_FUDGE_FACTOR, EMDL_DOUBLE, "rlnTau2FudgeFactor", "Regularisation parameter with which estimates for the power in the references will be multiplied (T in original paper)");
		EMDL::addLabel(EMDL_MLMODEL_TAU2_REF, EMDL_DOUBLE, "rlnReferenceTau2", "Spherical average of the estimated power in the signal of a reference");

		EMDL::addLabel(EMDL_OPTIMISER_ACCURACY_ROT, EMDL_DOUBLE, "rlnOverallAccuracyRotations", "Overall accuracy of the rotational assignments (in degrees)");
		EMDL::addLabel(EMDL_OPTIMISER_ACCURACY_TRANS, EMDL_DOUBLE, "rlnOverallAccuracyTranslations", "Overall accuracy of the translational assignments (in pixels)");
		EMDL::addLabel(EMDL_OPTIMISER_ACCURACY_TRANS_ANGSTROM, EMDL_DOUBLE, "rlnOverallAccuracyTranslationsAngst", "Overall accuracy of the translational assignments (in Angstroms)");
		EMDL::addLabel(EMDL_OPTIMISER_ADAPTIVE_FRACTION, EMDL_DOUBLE, "rlnAdaptiveOversampleFraction", "Fraction of the weights that will be oversampled in a second pass of the adaptive oversampling strategy");
		EMDL::addLabel(EMDL_OPTIMISER_ADAPTIVE_OVERSAMPLING, EMDL_INT, "rlnAdaptiveOversampleOrder", "Order of the adaptive oversampling (0=no oversampling, 1= 2x oversampling; 2= 4x oversampling, etc)");
		EMDL::addLabel(EMDL_OPTIMISER_AUTO_LOCAL_HP_ORDER, EMDL_INT, "rlnAutoLocalSearchesHealpixOrder", "Healpix order (before oversampling) from which autosampling procedure will use local angular searches");
		EMDL::addLabel(EMDL_OPTIMISER_AVAILABLE_MEMORY, EMDL_DOUBLE, "rlnAvailableMemory", "Available memory per computing node (i.e. per MPI-process)");
		EMDL::addLabel(EMDL_OPTIMISER_BEST_RESOL_THUS_FAR, EMDL_DOUBLE, "rlnBestResolutionThusFar", "The highest resolution that has been obtained in this optimization thus far");
		EMDL::addLabel(EMDL_OPTIMISER_COARSE_SIZE, EMDL_INT, "rlnCoarseImageSize", "Current size of the images to be used in the first pass of the adaptive oversampling strategy (may be smaller than the original image size)");
		EMDL::addLabel(EMDL_OPTIMISER_CHANGES_OPTIMAL_OFFSETS, EMDL_DOUBLE, "rlnChangesOptimalOffsets", "The average change in optimal translation in the last iteration (in pixels) ");
		EMDL::addLabel(EMDL_OPTIMISER_CHANGES_OPTIMAL_ORIENTS, EMDL_DOUBLE, "rlnChangesOptimalOrientations", "The average change in optimal orientation in the last iteration (in degrees) ");
		EMDL::addLabel(EMDL_OPTIMISER_CHANGES_OPTIMAL_CLASSES, EMDL_DOUBLE, "rlnChangesOptimalClasses", "The number of particles that changed their optimal clsas assignment in the last iteration");
		EMDL::addLabel(EMDL_OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, EMDL_BOOL, "rlnCtfDataArePhaseFlipped", "Flag to indicate that the input images have been phase-flipped");
		EMDL::addLabel(EMDL_OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, EMDL_BOOL, "rlnCtfDataAreCtfPremultiplied", "Flag to indicate that the input images have been premultiplied with their CTF");
		EMDL::addLabel(EMDL_OPTIMISER_DATA_STARFILE, EMDL_STRING, "rlnExperimentalDataStarFile", "STAR file with metadata for the experimental images");
		EMDL::addLabel(EMDL_OPTIMISER_DO_CORRECT_CTF, EMDL_BOOL, "rlnDoCorrectCtf", "Flag to indicate that CTF-correction should be performed");
		EMDL::addLabel(EMDL_OPTIMISER_DO_CORRECT_MAGNIFICATION, EMDL_BOOL, "rlnDoCorrectMagnification", "Flag to indicate that (per-group) magnification correction should be performed");
		EMDL::addLabel(EMDL_OPTIMISER_DO_CORRECT_NORM, EMDL_BOOL, "rlnDoCorrectNorm", "Flag to indicate that (per-image) normalisation-error correction should be performed");
		EMDL::addLabel(EMDL_OPTIMISER_DO_CORRECT_SCALE, EMDL_BOOL, "rlnDoCorrectScale", "Flag to indicate that internal (per-group) intensity-scale correction should be performed");
		EMDL::addLabel(EMDL_OPTIMISER_DO_EXTERNAL_RECONSTRUCT, EMDL_BOOL, "rlnDoExternalReconstruct", "Flag to indicate that the reconstruction will be performed outside relion_refine, e.g. for learned priors");
		EMDL::addLabel(EMDL_OPTIMISER_DO_REALIGN_MOVIES, EMDL_BOOL, "rlnDoRealignMovies", "Flag to indicate that individual frames of movies are being re-aligned");
		EMDL::addLabel(EMDL_OPTIMISER_DO_MAP, EMDL_BOOL, "rlnDoMapEstimation", "Flag to indicate that MAP estimation should be performed (otherwise ML estimation)");
		EMDL::addLabel(EMDL_OPTIMISER_DO_VMGD, EMDL_BOOL, "rlnDoStochasticGradientDescent", "Flag to indicate that SGD-optimisation should be performed (otherwise expectation maximisation)");
		EMDL::addLabel(EMDL_OPTIMISER_DO_STOCHASTIC_EM,EMDL_BOOL, "rlnDoStochasticEM", "Flag to indicate that stochastic EM-optimisation should be performed (an alternative to SGD)");
		EMDL::addLabel(EMDL_OPTIMISER_EXTERNAL_RECONS_DATA_REAL, EMDL_STRING, "rlnExtReconsDataReal", "Name of the map with the real components of the input data array for the external reconstruction program");
		EMDL::addLabel(EMDL_OPTIMISER_EXTERNAL_RECONS_DATA_IMAG, EMDL_STRING, "rlnExtReconsDataImag", "Name of the map with the imaginary components of the input data array for the external reconstruction program");
		EMDL::addLabel(EMDL_OPTIMISER_EXTERNAL_RECONS_WEIGHT, EMDL_STRING, "rlnExtReconsWeight", "Name of the map with the input weight array for the external reconstruction program");
		EMDL::addLabel(EMDL_OPTIMISER_EXTERNAL_RECONS_RESULT, EMDL_STRING, "rlnExtReconsResult", "Name of the output reconstruction from the external reconstruction program");
		EMDL::addLabel(EMDL_OPTIMISER_EXTERNAL_RECONS_NEWSTAR, EMDL_STRING, "rlnExtReconsResultStarfile", "Name of the output STAR file with updated FSC or tau curves");
		EMDL::addLabel(EMDL_OPTIMISER_FAST_SUBSETS, EMDL_BOOL, "rlnDoFastSubsetOptimisation", "Use subsets of the data in the earlier iterations to speed up convergence");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_INI_ITER, EMDL_INT, "rlnSgdInitialIterations", "Number of initial SGD iterations (at rlnSgdInitialResolution and with rlnSgdInitialSubsetSize)");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_FIN_ITER, EMDL_INT, "rlnSgdFinalIterations", "Number of final SGD iterations (at rlnSgdFinalResolution and with rlnSgdFinalSubsetSize)");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_INBETWEEN_ITER, EMDL_INT, "rlnSgdInBetweenIterations", "Number of SGD iteration in between the initial ones to the final ones (with linear interpolation of resolution and subset size)");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_INI_RESOL, EMDL_DOUBLE, "rlnSgdInitialResolution", "Resolution (in A) to use during the initial SGD iterations");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_FIN_RESOL, EMDL_DOUBLE, "rlnSgdFinalResolution", "Resolution (in A) to use during the final SGD iterations");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_INI_SUBSET_SIZE, EMDL_INT, "rlnSgdInitialSubsetSize", "Number of particles in a mini-batch (subset) during the initial SGD iterations");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_FIN_SUBSET_SIZE, EMDL_INT, "rlnSgdFinalSubsetSize", "Number of particles in a mini-batch (subset) during the final SGD iteration");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_MU, EMDL_DOUBLE, "rlnSgdMuFactor", "The mu-parameter that controls the momentum of the SGD gradients");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_SIGMA2FUDGE_INI, EMDL_DOUBLE, "rlnSgdSigma2FudgeInitial", "The variance of the noise will initially be multiplied with this value (larger than 1)");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE, EMDL_INT, "rlnSgdSigma2FudgeHalflife", "After processing this many particles the multiplicative factor for the noise variance will have halved");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_SKIP_ANNNEAL, EMDL_BOOL, "rlnSgdSkipAnneal", "Option to switch off annealing of multiple references in SGD");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_SUBSET_SIZE, EMDL_INT, "rlnSgdSubsetSize", "The number of particles in the random subsets for SGD");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_WRITE_EVERY_SUBSET, EMDL_INT, "rlnSgdWriteEverySubset", "Every this many iterations the model is written to disk in SGD");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_MAX_SUBSETS, EMDL_INT, "rlnSgdMaxSubsets", "Stop SGD after doing this many subsets (possibly spanning more than 1 iteration)");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_STEPSIZE, EMDL_DOUBLE, "rlnSgdStepsize", "Stepsize in SGD updates)");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_DO_MOM1, EMDL_BOOL, "rlnSgdDoMoment1", "Do first moment in SGD updates)");
		EMDL::addLabel(EMDL_OPTIMISER_SGD_DO_MOM2, EMDL_BOOL, "rlnSgdDoMoment2", "Do second moment in SGD updates)");
		EMDL::addLabel(EMDL_OPTIMISER_DO_AUTO_REFINE, EMDL_BOOL, "rlnDoAutoRefine", "Flag to indicate that 3D auto-refine procedure is being used");
		EMDL::addLabel(EMDL_OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, EMDL_BOOL, "rlnDoOnlyFlipCtfPhases", "Flag to indicate that CTF-correction should only comprise phase-flipping");
		EMDL::addLabel(EMDL_OPTIMISER_DO_SOLVENT_FLATTEN, EMDL_BOOL, "rlnDoSolventFlattening", "Flag to indicate that the references should be masked to set their solvent areas to a constant density");
		EMDL::addLabel(EMDL_OPTIMISER_DO_SOLVENT_FSC, EMDL_BOOL, "rlnDoSolventFscCorrection", "Flag to indicate that the FSCs should be solvent-corrected during refinement");
		EMDL::addLabel(EMDL_OPTIMISER_DO_SKIP_ALIGN, EMDL_BOOL, "rlnDoSkipAlign", "Flag to indicate that orientational (i.e. rotational and translational) searches will be omitted from the refinement, only marginalisation over classes will take place");
		EMDL::addLabel(EMDL_OPTIMISER_DO_SKIP_ROTATE, EMDL_BOOL, "rlnDoSkipRotate", "Flag to indicate that rotational searches will be omitted from the refinement, only marginalisation over classes and translations will take place");
		EMDL::addLabel(EMDL_OPTIMISER_DO_SPLIT_RANDOM_HALVES, EMDL_BOOL, "rlnDoSplitRandomHalves", "Flag to indicate that the data should be split into two completely separate, random halves");
		EMDL::addLabel(EMDL_OPTIMISER_DO_ZERO_MASK, EMDL_BOOL, "rlnDoZeroMask", "Flag to indicate that the surrounding solvent area in the experimental particles will be masked to zeros (by default random noise will be used");
		EMDL::addLabel(EMDL_OPTIMISER_FIX_SIGMA_NOISE, EMDL_BOOL, "rlnFixSigmaNoiseEstimates", "Flag to indicate that the estimates for the power spectra of the noise should be kept constant");
		EMDL::addLabel(EMDL_OPTIMISER_FIX_SIGMA_OFFSET ,EMDL_BOOL, "rlnFixSigmaOffsetEstimates", "Flag to indicate that the estimates for the stddev in the origin offsets should be kept constant");
		EMDL::addLabel(EMDL_OPTIMISER_FIX_TAU, EMDL_BOOL, "rlnFixTauEstimates", "Flag to indicate that the estimates for the power spectra of the signal (i.e. the references) should be kept constant");
		EMDL::addLabel(EMDL_OPTIMISER_HAS_CONVERGED, EMDL_BOOL, "rlnHasConverged", "Flag to indicate that the optimization has converged");
		EMDL::addLabel(EMDL_OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, EMDL_BOOL, "rlnHasHighFscAtResolLimit", "Flag to indicate that the FSC at the resolution limit is significant");
		EMDL::addLabel(EMDL_OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, EMDL_INT, "rlnHasLargeSizeIncreaseIterationsAgo", "How many iterations have passed since the last large increase in image size");
		EMDL::addLabel(EMDL_OPTIMISER_DO_HELICAL_REFINE, EMDL_BOOL, "rlnDoHelicalRefine", "Flag to indicate that helical refinement should be performed");
		EMDL::addLabel(EMDL_OPTIMISER_IGNORE_HELICAL_SYMMETRY, EMDL_BOOL, "rlnIgnoreHelicalSymmetry", "Flag to indicate that helical symmetry is ignored in 3D reconstruction");
		EMDL::addLabel(EMDL_OPTIMISER_FOURIER_MASK, EMDL_STRING, "rlnFourierMask", "Name of an FFTW-centred Fourier mask to be applied to the Projector for refinement.");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_TWIST_INITIAL, EMDL_DOUBLE, "rlnHelicalTwistInitial", "The intial helical twist (rotation per subunit) in degrees before refinement");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_RISE_INITIAL, EMDL_DOUBLE, "rlnHelicalRiseInitial", "The initial helical rise (translation per subunit) in Angstroms before refinement");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_Z_PERCENTAGE, EMDL_DOUBLE, "rlnHelicalCentralProportion", "Only expand this central fraction of the Z axis when imposing real-space helical symmetry");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_NSTART, EMDL_INT, "rlnNrHelicalNStart", "The N-number for an N-start helix");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_TUBE_INNER_DIAMETER, EMDL_DOUBLE, "rlnHelicalMaskTubeInnerDiameter", "Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER, EMDL_DOUBLE, "rlnHelicalMaskTubeOuterDiameter", "Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT, EMDL_BOOL, "rlnHelicalSymmetryLocalRefinement", "Flag to indicate that local refinement of helical parameters should be performed");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_SIGMA_DISTANCE, EMDL_DOUBLE, "rlnHelicalSigmaDistance", "Sigma of distance along the helical tracks");
		EMDL::addLabel(EMDL_OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED, EMDL_BOOL, "rlnHelicalKeepTiltPriorFixed", "Flag to indicate that helical tilt priors are kept fixed (at 90 degrees) in global angular searches");
		EMDL::addLabel(EMDL_OPTIMISER_LOWRES_LIMIT_EXP, EMDL_DOUBLE, "rlnLowresLimitExpectation", "Low-resolution-limit (in Angstrom) for the expectation step");
		EMDL::addLabel(EMDL_OPTIMISER_HIGHRES_LIMIT_EXP, EMDL_DOUBLE, "rlnHighresLimitExpectation", "High-resolution-limit (in Angstrom) for the expectation step");
		EMDL::addLabel(EMDL_OPTIMISER_HIGHRES_LIMIT_SGD, EMDL_DOUBLE, "rlnHighresLimitSGD", "High-resolution-limit (in Angstrom) for Stochastic Gradient Descent");
		EMDL::addLabel(EMDL_OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, EMDL_BOOL, "rlnDoIgnoreCtfUntilFirstPeak", "Flag to indicate that the CTFs should be ignored until their first peak");
		EMDL::addLabel(EMDL_OPTIMISER_INCR_SIZE, EMDL_INT, "rlnIncrementImageSize", "Number of Fourier shells to be included beyond the resolution where SSNR^MAP drops below 1");
		EMDL::addLabel(EMDL_OPTIMISER_ITERATION_NO, EMDL_INT, "rlnCurrentIteration", "The number of the current iteration");
		EMDL::addLabel(EMDL_OPTIMISER_LOCAL_SYMMETRY_FILENAME, EMDL_STRING, "rlnLocalSymmetryFile", "Local symmetry description file containing list of masks and their operators");
		EMDL::addLabel(EMDL_OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, EMDL_DOUBLE, "rlnJoinHalvesUntilThisResolution", "Resolution (in Angstrom) to join the two random half-reconstructions to prevent their diverging orientations (for C-symmetries)");
		EMDL::addLabel(EMDL_OPTIMISER_MAGNIFICATION_RANGE, EMDL_DOUBLE, "rlnMagnificationSearchRange", "Search range for magnification correction");
		EMDL::addLabel(EMDL_OPTIMISER_MAGNIFICATION_STEP, EMDL_DOUBLE, "rlnMagnificationSearchStep", "Step size  for magnification correction");
		EMDL::addLabel(EMDL_OPTIMISER_MAX_COARSE_SIZE, EMDL_INT, "rlnMaximumCoarseImageSize", "Maximum size of the images to be used in the first pass of the adaptive oversampling strategy (may be smaller than the original image size)");
		EMDL::addLabel(EMDL_OPTIMISER_MAX_NR_POOL, EMDL_INT, "rlnMaxNumberOfPooledParticles", "Maximum number particles that are processed together to speed up calculations");
		EMDL::addLabel(EMDL_OPTIMISER_MODEL_STARFILE, EMDL_STRING, "rlnModelStarFile", "STAR file with metadata for the model that is being refined");
		EMDL::addLabel(EMDL_OPTIMISER_MODEL_STARFILE2, EMDL_STRING, "rlnModelStarFile2", "STAR file with metadata for the second model that is being refined (from random halves of the data)");
		EMDL::addLabel(EMDL_OPTIMISER_NR_ITERATIONS, EMDL_INT, "rlnNumberOfIterations", "Maximum number of iterations to be performed");
		EMDL::addLabel(EMDL_OPTIMISER_NR_ITER_WO_RESOL_GAIN, EMDL_INT, "rlnNumberOfIterWithoutResolutionGain", "Number of iterations that have passed without a gain in resolution");
		EMDL::addLabel(EMDL_OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, EMDL_INT, "rlnNumberOfIterWithoutChangingAssignments", "Number of iterations that have passed without large changes in orientation and class assignments");
		EMDL::addLabel(EMDL_OPTIMISER_OPTICS_STARFILE, EMDL_STRING, "rlnOpticsStarFile", "STAR file with metadata for the optical groups (new as of version 3.1)");
		EMDL::addLabel(EMDL_OPTIMISER_OUTPUT_ROOTNAME, EMDL_STRING, "rlnOutputRootName", "Rootname for all output files (this may include a directory structure, which should then exist)");
		EMDL::addLabel(EMDL_OPTIMISER_PARTICLE_DIAMETER, EMDL_DOUBLE, "rlnParticleDiameter", "Diameter of the circular mask to be applied to all experimental images (in Angstroms)");
		EMDL::addLabel(EMDL_OPTIMISER_RADIUS_MASK_3D_MAP, EMDL_INT, "rlnRadiusMaskMap", "Radius of the spherical mask to be applied to all references (in Angstroms)");
		EMDL::addLabel(EMDL_OPTIMISER_RADIUS_MASK_EXP_PARTICLES, EMDL_INT, "rlnRadiusMaskExpImages", "Radius of the circular mask to be applied to all experimental images (in Angstroms)");
		EMDL::addLabel(EMDL_OPTIMISER_RANDOM_SEED, EMDL_INT, "rlnRandomSeed", "Seed (i.e. a number) for the random number generator");
		EMDL::addLabel(EMDL_OPTIMISER_REFS_ARE_CTF_CORRECTED, EMDL_BOOL, "rlnRefsAreCtfCorrected", "Flag to indicate that the input references have been CTF-amplitude corrected");
		EMDL::addLabel(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, EMDL_INT, "rlnSmallestChangesClasses", "Smallest changes thus far in the optimal class assignments (in numer of particles).");
		EMDL::addLabel(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, EMDL_DOUBLE, "rlnSmallestChangesOffsets", "Smallest changes thus far in the optimal offset assignments (in pixels).");
		EMDL::addLabel(EMDL_OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, EMDL_DOUBLE, "rlnSmallestChangesOrientations", "Smallest changes thus far in the optimal orientation assignments (in degrees).");
		EMDL::addLabel(EMDL_OPTIMISER_SAMPLING_STARFILE, EMDL_STRING, "rlnOrientSamplingStarFile", "STAR file with metadata for the orientational sampling");
		EMDL::addLabel(EMDL_OPTIMISER_SOLVENT_MASK_NAME, EMDL_STRING, "rlnSolventMaskName", "Name of an image that contains a (possibly soft) mask for the solvent area (values=0 for solvent, values =1 for protein)");
		EMDL::addLabel(EMDL_OPTIMISER_SOLVENT_MASK2_NAME, EMDL_STRING, "rlnSolventMask2Name", "Name of a secondary solvent mask (e.g. to flatten density inside an icosahedral virus)");
		EMDL::addLabel(EMDL_OPTIMISER_TAU_SPECTRUM_NAME, EMDL_STRING, "rlnTauSpectrumName", "Name of a STAR file that holds a tau2-spectrum");
		EMDL::addLabel(EMDL_OPTIMISER_USE_TOO_COARSE_SAMPLING, EMDL_BOOL, "rlnUseTooCoarseSampling", "Flag to indicate that the angular sampling on the sphere will be one step coarser than needed to speed up calculations");
		EMDL::addLabel(EMDL_OPTIMISER_WIDTH_MASK_EDGE, EMDL_INT, "rlnWidthMaskEdge", "Width (in pixels) of the soft edge for spherical/circular masks to be used for solvent flattening");

		EMDL::addLabel(EMDL_ORIENT_FLIP, EMDL_BOOL, "rlnIsFlip", "Flag to indicate that an image should be mirrored");
		EMDL::addLabel(EMDL_ORIENT_ID, EMDL_INT, "rlnOrientationsID", "ID (i.e. a unique number) for an orientation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_X, EMDL_DOUBLE, "rlnOriginX", "X-coordinate (in pixels) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_Y, EMDL_DOUBLE, "rlnOriginY", "Y-coordinate (in pixels) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_Z, EMDL_DOUBLE, "rlnOriginZ", "Z-coordinate (in pixels) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_X_PRIOR, EMDL_DOUBLE, "rlnOriginXPrior", "Center of the prior on the X-coordinate (in pixels) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR, EMDL_DOUBLE, "rlnOriginYPrior", "Center of the prior on the Y-coordinate (in pixels) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_Z_PRIOR, EMDL_DOUBLE, "rlnOriginZPrior", "Center of the prior on the Z-coordinate (in pixels) for the origin of rotation");

		EMDL::addLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM, EMDL_DOUBLE, "rlnOriginXAngst", "X-coordinate (in Angstrom) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, EMDL_DOUBLE, "rlnOriginYAngst", "Y-coordinate (in Angstrom) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, EMDL_DOUBLE, "rlnOriginZAngst", "Z-coordinate (in Angstrom) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_X_PRIOR_ANGSTROM, EMDL_DOUBLE, "rlnOriginXPriorAngst", "Center of the prior on the X-coordinate (in Angstrom) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_Y_PRIOR_ANGSTROM, EMDL_DOUBLE, "rlnOriginYPriorAngst", "Center of the prior on the Y-coordinate (in Angstrom) for the origin of rotation");
		EMDL::addLabel(EMDL_ORIENT_ORIGIN_Z_PRIOR_ANGSTROM, EMDL_DOUBLE, "rlnOriginZPriorAngst", "Center of the prior on the Z-coordinate (in Angstrom) for the origin of rotation");

		EMDL::addLabel(EMDL_ORIENT_ROT, EMDL_DOUBLE, "rlnAngleRot", "First Euler angle (rot, in degrees)");
		EMDL::addLabel(EMDL_ORIENT_ROT_PRIOR, EMDL_DOUBLE, "rlnAngleRotPrior", "Center of the prior (in degrees) on the first Euler angle (rot)");
		EMDL::addLabel(EMDL_ORIENT_ROT_PRIOR_FLIP_RATIO, EMDL_DOUBLE, "rlnAngleRotFlipRatio", "Flip ratio of bimodal rot prior (0~0.5, 0 means an ordinary prior, 0.5 means a perfect bimodal prior)");   // KThurber
		EMDL::addLabel(EMDL_ORIENT_TILT, EMDL_DOUBLE, "rlnAngleTilt", "Second Euler angle (tilt, in degrees)");
		EMDL::addLabel(EMDL_ORIENT_TILT_PRIOR, EMDL_DOUBLE, "rlnAngleTiltPrior", "Center of the prior (in degrees) on the second Euler angle (tilt)");
		EMDL::addLabel(EMDL_ORIENT_PSI, EMDL_DOUBLE, "rlnAnglePsi", "Third Euler, or in-plane angle (psi, in degrees)");
		EMDL::addLabel(EMDL_ORIENT_PSI_PRIOR, EMDL_DOUBLE, "rlnAnglePsiPrior", "Center of the prior (in degrees) on the third Euler angle (psi)");
		EMDL::addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, EMDL_DOUBLE, "rlnAnglePsiFlipRatio", "Flip ratio of bimodal psi prior (0~0.5, 0 means an ordinary prior, 0.5 means a perfect bimodal prior)");
		EMDL::addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP, EMDL_BOOL, "rlnAnglePsiFlip", "Flag to indicate that psi prior angle has been flipped");  // KThurber

		EMDL::addLabel(EMDL_PARTICLE_AUTOPICK_FOM, EMDL_DOUBLE, "rlnAutopickFigureOfMerit", "Autopicking FOM for a particle");
		EMDL::addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID, EMDL_INT, "rlnHelicalTubeID", "Helical tube ID for a helical segment");
		EMDL::addLabel(EMDL_PARTICLE_HELICAL_TUBE_PITCH, EMDL_DOUBLE, "rlnHelicalTubePitch", "Cross-over distance for a helical segment (A)");
		EMDL::addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH, EMDL_DOUBLE, "rlnHelicalTrackLength", "Distance (in pix) from the position of this helical segment to the starting point of the tube");
		EMDL::addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, EMDL_DOUBLE, "rlnHelicalTrackLengthAngst", "Distance (in A) from the position of this helical segment to the starting point of the tube");
		EMDL::addLabel(EMDL_PARTICLE_CLASS, EMDL_INT, "rlnClassNumber", "Class number for which a particle has its highest probability");
		EMDL::addLabel(EMDL_PARTICLE_DLL, EMDL_DOUBLE, "rlnLogLikeliContribution", "Contribution of a particle to the log-likelihood target function");
		EMDL::addLabel(EMDL_PARTICLE_ID, EMDL_INT, "rlnParticleId", "ID (i.e. a unique number) for a particle");
		EMDL::addLabel(EMDL_PARTICLE_FOM, EMDL_DOUBLE, "rlnParticleFigureOfMerit", "Developmental FOM for a particle");
		EMDL::addLabel(EMDL_PARTICLE_KL_DIVERGENCE, EMDL_DOUBLE, "rlnKullbackLeiblerDivergence", "Kullback-Leibler divergence for a particle");
		EMDL::addAltLabel(EMDL_PARTICLE_KL_DIVERGENCE,			 "rlnKullbackLeibnerDivergence"); // wrong spelling for backwards compatibility
		EMDL::addLabel(EMDL_PARTICLE_RANDOM_SUBSET, EMDL_INT, "rlnRandomSubset", "Random subset to which this particle belongs");
		EMDL::addLabel(EMDL_PARTICLE_BEAM_TILT_CLASS, EMDL_INT, "rlnBeamTiltClass", "Beam-tilt class of a particle");
		EMDL::addLabel(EMDL_PARTICLE_NAME, EMDL_STRING, "rlnParticleName", "Name for a particle");
		EMDL::addLabel(EMDL_PARTICLE_ORI_NAME, EMDL_STRING, "rlnOriginalParticleName", "Original name for a particles");
		EMDL::addLabel(EMDL_PARTICLE_NR_SIGNIFICANT_SAMPLES, EMDL_INT, "rlnNrOfSignificantSamples", "Number of orientational/class assignments (for a particle) with sign.probabilities in the 1st pass of adaptive oversampling"); /**< particle, Number of orientations contributing to weights*/
		EMDL::addLabel(EMDL_PARTICLE_NR_FRAMES, EMDL_INT, "rlnNrOfFrames", "Number of movie frames that were collected for this particle");
		EMDL::addLabel(EMDL_PARTICLE_NR_FRAMES_AVG, EMDL_INT, "rlnAverageNrOfFrames", "Number of movie frames that one averages over upon extraction of movie-particles");
		EMDL::addLabel(EMDL_PARTICLE_MOVIE_RUNNING_AVG, EMDL_INT, "rlnMovieFramesRunningAverage", "Number of movie frames inside the running average that will be used for movie-refinement");
		EMDL::addLabel(EMDL_PARTICLE_PMAX, EMDL_DOUBLE, "rlnMaxValueProbDistribution", "Maximum value of the (normalised) probability function for a particle"); /**< particle, Maximum value of probability distribution */
		EMDL::addLabel(EMDL_PARTICLE_NUMBER, EMDL_INT, "rlnParticleNumber", "Number of particles");

		EMDL::addLabel(EMDL_PIPELINE_JOB_COUNTER, EMDL_INT, "rlnPipeLineJobCounter", "Number of the last job in the pipeline");
		EMDL::addLabel(EMDL_PIPELINE_NODE_NAME, EMDL_STRING , "rlnPipeLineNodeName", "Name of a Node in the pipeline");
		EMDL::addLabel(EMDL_PIPELINE_NODE_TYPE, EMDL_INT, "rlnPipeLineNodeType", "Type of a Node in the pipeline");
		EMDL::addLabel(EMDL_PIPELINE_PROCESS_ALIAS, EMDL_STRING , "rlnPipeLineProcessAlias", "Alias of a Process in the pipeline");
		EMDL::addLabel(EMDL_PIPELINE_PROCESS_NAME, EMDL_STRING , "rlnPipeLineProcessName", "Name of a Process in the pipeline");
		EMDL::addLabel(EMDL_PIPELINE_PROCESS_TYPE, EMDL_INT, "rlnPipeLineProcessType", "Type of a Process in the pipeline");
		EMDL::addLabel(EMDL_PIPELINE_PROCESS_STATUS, EMDL_INT, "rlnPipeLineProcessStatus", "Status of a Process in the pipeline (running, scheduled, finished or cancelled)");
		EMDL::addLabel(EMDL_PIPELINE_EDGE_FROM, EMDL_STRING , "rlnPipeLineEdgeFromNode", "Name of the origin of an edge");
		EMDL::addLabel(EMDL_PIPELINE_EDGE_TO, EMDL_STRING ,"rlnPipeLineEdgeToNode", "Name of the to-Node in an edge");
		EMDL::addLabel(EMDL_PIPELINE_EDGE_PROCESS, EMDL_STRING ,"rlnPipeLineEdgeProcess", "Name of the destination of an edge");

		EMDL::addLabel(EMDL_POSTPROCESS_FINAL_RESOLUTION, EMDL_DOUBLE, "rlnFinalResolution", "Final estimated resolution after postprocessing (in Angstroms)");
		EMDL::addLabel(EMDL_POSTPROCESS_BFACTOR, EMDL_DOUBLE, "rlnBfactorUsedForSharpening", "Applied B-factor in the sharpening of the map");
		EMDL::addLabel(EMDL_POSTPROCESS_FRACTION_MOLWEIGHT, EMDL_DOUBLE, "rlnParticleBoxFractionMolecularWeight", "Fraction of protein voxels in the box, based on ordered molecular weight estimate, for calculating cisTEM-like part_FSC");
		EMDL::addLabel(EMDL_POSTPROCESS_FRACTION_SOLVENT_MASK, EMDL_DOUBLE, "rlnParticleBoxFractionSolventMask", "Fraction of protein voxels in the box, based on the solvent mask, for calculating cisTEM-like part_FSC");
		EMDL::addLabel(EMDL_POSTPROCESS_FSC_GENERAL, EMDL_DOUBLE, "rlnFourierShellCorrelation", "FSC value (of unspecified type, e.g. masked or unmasked)");
		EMDL::addLabel(EMDL_POSTPROCESS_FSC_TRUE, EMDL_DOUBLE, "rlnFourierShellCorrelationCorrected", "Final FSC value: i.e. after correction based on masking of randomized-phases maps");
		EMDL::addLabel(EMDL_POSTPROCESS_FSC_PART_MOLWEIGHT, EMDL_DOUBLE, "rlnFourierShellCorrelationParticleMolWeight", "CisTEM-like correction of unmasked FSCs, based on ordered molecular weight estimate");
		EMDL::addLabel(EMDL_POSTPROCESS_FSC_PART_FRACMASK, EMDL_DOUBLE, "rlnFourierShellCorrelationParticleMaskFraction", "CisTEM-like correction of unmasked FSCs, based on fraction of white pixels in solvent mask");
		EMDL::addLabel(EMDL_POSTPROCESS_FSC_MASKED, EMDL_DOUBLE, "rlnFourierShellCorrelationMaskedMaps", "FSC value after masking of the original maps");
		EMDL::addLabel(EMDL_POSTPROCESS_FSC_UNMASKED, EMDL_DOUBLE, "rlnFourierShellCorrelationUnmaskedMaps", "FSC value before masking of the original maps");
		EMDL::addLabel(EMDL_POSTPROCESS_FSC_RANDOM_MASKED, EMDL_DOUBLE, "rlnCorrectedFourierShellCorrelationPhaseRandomizedMaskedMaps", "FSC value after masking of the randomized-phases maps");
		EMDL::addLabel(EMDL_POSTPROCESS_AMPLCORR_MASKED, EMDL_DOUBLE, "rlnAmplitudeCorrelationMaskedMaps", "Correlation coefficient between amplitudes in Fourier shells of masked maps");
		EMDL::addLabel(EMDL_POSTPROCESS_AMPLCORR_UNMASKED, EMDL_DOUBLE, "rlnAmplitudeCorrelationUnmaskedMaps", "Correlation coefficient between amplitudes in Fourier shells of unmasked maps");
		EMDL::addLabel(EMDL_POSTPROCESS_DPR_MASKED, EMDL_DOUBLE, "rlnDifferentialPhaseResidualMaskedMaps", "Differential Phase Residual in Fourier shells of masked maps");
		EMDL::addLabel(EMDL_POSTPROCESS_DPR_UNMASKED,  EMDL_DOUBLE, "rlnDifferentialPhaseResidualUnmaskedMaps", "Differential Phase Residual in Fourier shells of unmasked maps");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, EMDL_DOUBLE, "rlnFittedInterceptGuinierPlot", "The fitted intercept of the Guinier-plot");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_FIT_SLOPE, EMDL_DOUBLE, "rlnFittedSlopeGuinierPlot", "The fitted slope of the Guinier-plot");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_FIT_CORRELATION, EMDL_DOUBLE, "rlnCorrelationFitGuinierPlot", "The correlation coefficient of the fitted line through the Guinier-plot");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_VALUE_IN, EMDL_DOUBLE, "rlnLogAmplitudesOriginal", "Y-value for Guinier plot: the logarithm of the radially averaged amplitudes of the input map");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_VALUE_INVMTF, EMDL_DOUBLE, "rlnLogAmplitudesMTFCorrected", "Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after MTF correction");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_VALUE_WEIGHTED, EMDL_DOUBLE, "rlnLogAmplitudesWeighted", "Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after FSC-weighting");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_VALUE_SHARPENED, EMDL_DOUBLE, "rlnLogAmplitudesSharpened", "Y-value for Guinier plot: the logarithm of the radially averaged amplitudes after sharpening");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_VALUE_INTERCEPT, EMDL_DOUBLE, "rlnLogAmplitudesIntercept", "Y-value for Guinier plot: the fitted plateau of the logarithm of the radially averaged amplitudes");
		EMDL::addLabel(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, EMDL_DOUBLE, "rlnResolutionSquared", "X-value for Guinier plot: squared resolution in 1/Angstrom^2");
		EMDL::addLabel(EMDL_POSTPROCESS_MOLWEIGHT, EMDL_DOUBLE, "rlnMolecularWeight", "Molecular weight of the ordered mass inside the box for calculating cisTEM-like part.FSC (in kDa)");
		EMDL::addLabel(EMDL_POSTPROCESS_MTF_VALUE, EMDL_DOUBLE, "rlnMtfValue", "Value of the detectors modulation transfer function (between 0 and 1)");
		EMDL::addLabel(EMDL_POSTPROCESS_RANDOMISE_FROM, EMDL_DOUBLE, "rlnRandomiseFrom", "Resolution (in A) from which the phases are randomised in the postprocessing step");
		EMDL::addLabel(EMDL_POSTPROCESS_UNFIL_HALFMAP1, EMDL_STRING, "rlnUnfilteredMapHalf1", "Name of the unfiltered map from halfset 1");
		EMDL::addLabel(EMDL_POSTPROCESS_UNFIL_HALFMAP2, EMDL_STRING, "rlnUnfilteredMapHalf2", "Name of the unfiltered map from halfset 2");

		EMDL::addLabel(EMDL_SAMPLING_IS_3D, EMDL_BOOL, "rlnIs3DSampling", "Flag to indicate this concerns a 3D sampling ");
		EMDL::addLabel(EMDL_SAMPLING_IS_3D_TRANS, EMDL_BOOL, "rlnIs3DTranslationalSampling", "Flag to indicate this concerns a x,y,z-translational sampling ");
		EMDL::addLabel(EMDL_SAMPLING_HEALPIX_ORDER, EMDL_INT, "rlnHealpixOrder", "Healpix order for the sampling of the first two Euler angles (rot, tilt) on the 3D sphere");
		EMDL::addLabel(EMDL_SAMPLING_HEALPIX_ORDER_ORI, EMDL_INT, "rlnHealpixOrderOriginal", "Original healpix order for the sampling of the first two Euler angles (rot, tilt) on the 3D sphere");
		EMDL::addLabel(EMDL_SAMPLING_LIMIT_TILT, EMDL_DOUBLE, "rlnTiltAngleLimit", "Values to which to limit the tilt angles (positive for keeping side views, negative for keeping top views)");
		EMDL::addLabel(EMDL_SAMPLING_OFFSET_RANGE, EMDL_DOUBLE, "rlnOffsetRange", "Search range for the origin offsets (in Angstroms)");
		EMDL::addLabel(EMDL_SAMPLING_OFFSET_STEP, EMDL_DOUBLE, "rlnOffsetStep", "Step size for the searches in the origin offsets (in Angstroms)");
		EMDL::addLabel(EMDL_SAMPLING_OFFSET_RANGE_ORI, EMDL_DOUBLE, "rlnOffsetRangeOriginal", "Original search range for the origin offsets (in Angstroms)");
		EMDL::addLabel(EMDL_SAMPLING_OFFSET_STEP_ORI, EMDL_DOUBLE, "rlnOffsetStepOriginal", "Original step size for the searches in the origin offsets (in Angstroms)");
		EMDL::addLabel(EMDL_SAMPLING_HELICAL_OFFSET_STEP, EMDL_DOUBLE, "rlnHelicalOffsetStep", "Step size for the searches of offsets along helical axis (in Angstroms)");
		EMDL::addLabel(EMDL_SAMPLING_PERTURB, EMDL_DOUBLE, "rlnSamplingPerturbInstance", "Random instance of the random perturbation on the orientational sampling");
		EMDL::addLabel(EMDL_SAMPLING_PERTURBATION_FACTOR, EMDL_DOUBLE, "rlnSamplingPerturbFactor", "Factor for random perturbation on the orientational sampling (between 0 no perturbation and 1 very strong perturbation)");
		EMDL::addLabel(EMDL_SAMPLING_PSI_STEP, EMDL_DOUBLE, "rlnPsiStep", "Step size (in degrees) for the sampling of the in-plane rotation angle (psi)");
		EMDL::addLabel(EMDL_SAMPLING_PSI_STEP_ORI, EMDL_DOUBLE, "rlnPsiStepOriginal", "Original step size (in degrees) for the sampling of the in-plane rotation angle (psi)");
		EMDL::addLabel(EMDL_SAMPLING_SYMMETRY, EMDL_STRING, "rlnSymmetryGroup", "Symmetry group (e.g., C1, D7, I2, I5, etc.)");

		EMDL::addLabel(EMDL_SCHEDULE_EDGE_NUMBER, EMDL_INT, "rlnScheduleEdgeNumber", "Numbered index of an edge inside a Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_EDGE_INPUT, EMDL_STRING, "rlnScheduleEdgeInputNodeName" , "Name of the input Node for a schedule Edge");
		EMDL::addLabel(EMDL_SCHEDULE_EDGE_OUTPUT, EMDL_STRING, "rlnScheduleEdgeOutputNodeName", "Name of the output Node for a schedule Edge");
		EMDL::addLabel(EMDL_SCHEDULE_EDGE_IS_FORK, EMDL_BOOL, "rlnScheduleEdgeIsFork", "Flag to indicate that this Edge is a Fork, dependent on a Boolean Schedule variable");
		EMDL::addLabel(EMDL_SCHEDULE_EDGE_OUTPUT_TRUE, EMDL_STRING, "rlnScheduleEdgeOutputNodeNameIfTrue", "Name of the output Node for a schedule Fork if the associated Boolean is True");
		EMDL::addLabel(EMDL_SCHEDULE_EDGE_BOOLEAN, EMDL_STRING, "rlnScheduleEdgeBooleanVariable", "Name of the associated Boolean variable if this Edge is a Fork");
		EMDL::addLabel(EMDL_SCHEDULE_GENERAL_CURRENT_NODE, EMDL_STRING, "rlnScheduleCurrentNodeName", "Name of the current Node for this Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_GENERAL_ORIGINAL_START_NODE, EMDL_STRING, "rlnScheduleOriginalStartNodeName", "Name of the original starting Node for this Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_GENERAL_EMAIL, EMDL_STRING, "rlnScheduleEmailAddress", "Email address to send message when Schedule finishes");
		EMDL::addLabel(EMDL_SCHEDULE_GENERAL_NAME, EMDL_STRING, "rlnScheduleName", "Name for this Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_JOB_NAME, EMDL_STRING, "rlnScheduleJobName", "Name of a Job in a Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_JOB_ORI_NAME, EMDL_STRING, "rlnScheduleJobNameOriginal", "Original name of a Job in a Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_JOB_MODE, EMDL_STRING, "rlnScheduleJobMode", "Mode on how to execute a Job");
		EMDL::addLabel(EMDL_SCHEDULE_JOB_HAS_STARTED, EMDL_BOOL, "rlnScheduleJobHasStarted", "Flag to indicate whether a Job has started already in the execution of the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_OPERATOR_NAME, EMDL_STRING, "rlnScheduleOperatorName", "Name of a Boolean operator in the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_OPERATOR_TYPE, EMDL_STRING, "rlnScheduleOperatorType", "Type of an operator in the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_OPERATOR_INPUT1, EMDL_STRING, "rlnScheduleOperatorInput1", "Name of the 1st input to the operator");
		EMDL::addLabel(EMDL_SCHEDULE_OPERATOR_INPUT2, EMDL_STRING, "rlnScheduleOperatorInput2", "Name of the 2nd input to the operator");
		EMDL::addLabel(EMDL_SCHEDULE_OPERATOR_OUTPUT, EMDL_STRING, "rlnScheduleOperatorOutput", "Name of the output variable on which this operator acts");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_BOOL_NAME, EMDL_STRING, "rlnScheduleBooleanVariableName", "Name of a Boolean variable in the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_BOOL_VALUE, EMDL_BOOL, "rlnScheduleBooleanVariableValue", "Value of a Boolean variable in the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_BOOL_ORI_VALUE, EMDL_BOOL, "rlnScheduleBooleanVariableResetValue", "Value which a Boolean variable will take upon a reset");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_FLOAT_NAME, EMDL_STRING, "rlnScheduleFloatVariableName", "Name of a Float variable in the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_FLOAT_VALUE, EMDL_DOUBLE, "rlnScheduleFloatVariableValue", "Value of a Float variable in the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_FLOAT_ORI_VALUE, EMDL_DOUBLE, "rlnScheduleFloatVariableResetValue", "Value which a Float variable will take upon a reset");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_STRING_NAME, EMDL_STRING, "rlnScheduleStringVariableName", "Name of a String variable in the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_STRING_VALUE, EMDL_STRING, "rlnScheduleStringVariableValue", "Value of a String variable in the Schedule");
		EMDL::addLabel(EMDL_SCHEDULE_VAR_STRING_ORI_VALUE, EMDL_STRING, "rlnScheduleStringVariableResetValue", "Value which a String variable will take upon a reset");

		EMDL::addLabel(EMDL_SELECTED, EMDL_INT, "rlnSelected", "Flag whether an entry in a metadatatable is selected (1) in the viewer or not (0)");
		EMDL::addLabel(EMDL_SELECT_PARTICLES_ZSCORE, EMDL_DOUBLE, "rlnParticleSelectZScore", "Sum of Z-scores from particle_select. High Z-scores are likely to be outliers.");
		EMDL::addLabel(EMDL_SORTED_IDX, EMDL_INT, "rlnSortedIndex", "Index of a metadata entry after sorting (first sorted index is 0).");
		EMDL::addLabel(EMDL_STARFILE_MOVIE_PARTICLES, EMDL_STRING, "rlnStarFileMovieParticles", "Filename of a STAR file with movie-particles in it");
		EMDL::addLabel(EMDL_PERFRAME_CUMULATIVE_WEIGHT, EMDL_DOUBLE, "rlnPerFrameCumulativeWeight", "Sum of the resolution-dependent relative weights from the first frame until the given frame");
		EMDL::addLabel(EMDL_PERFRAME_RELATIVE_WEIGHT, EMDL_DOUBLE, "rlnPerFrameRelativeWeight", "The resolution-dependent relative weights for a given frame");

		EMDL::addLabel(EMDL_RESOLUTION, EMDL_DOUBLE, "rlnResolution", "Resolution (in 1/Angstroms)");
		EMDL::addLabel(EMDL_RESOLUTION_ANGSTROM, EMDL_DOUBLE, "rlnAngstromResolution", "Resolution (in Angstroms)");
		EMDL::addLabel(EMDL_RESOLUTION_INVPIXEL, EMDL_DOUBLE, "rlnResolutionInversePixel", "Resolution (in 1/pixel, Nyquist = 0.5)");
		EMDL::addLabel(EMDL_SPECTRAL_IDX, EMDL_INT, "rlnSpectralIndex", "Spectral index (i.e. distance in pixels to the origin in Fourier space) ");

		EMDL::addLabel(EMDL_UNKNOWN_LABEL, EMDL_UNKNOWN, "rlnUnknownLabel", "NON-RELION label: values will be ignored, yet maintained in the STAR file.");
	 }

	~StaticInitialization()
	{
	}
	friend class EMDL;
};

/**Just an utility function */
bool vectorContainsLabel(const std::vector<EMDLabel>& labelsVector, const EMDLabel label);

#endif

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

#ifndef LOCAL_SYMMETRY_H_
#define LOCAL_SYMMETRY_H_

#include "src/args.h"
#include "src/euler.h"
#include "src/funcs.h"
#include "src/macros.h"
#include "src/matrix1d.h"
#include "src/matrix2d.h"
#include "src/image.h"
#include "src/fftw.h"
#include "src/transformations.h"
#include "src/healpix_sampling.h"
#include "src/time.h"
#include <queue>

// DM (ccp4) operator types
// http://www.ccp4.ac.uk/html/rotationmatrices.html
#define ROTA_EULER_TYPE 1
#define ROTA_POLAR_TYPE 2
#define ROTA_MATRIX_TYPE 4
#define OMAT_TYPE 8

// Positions of parameters in Matrix1D operators
#define AA_POS 0
#define BB_POS 1
#define GG_POS 2
#define DX_POS 3
#define DY_POS 4
#define DZ_POS 5
#define CC_POS 6
#define NR_LOCALSYM_PARAMETERS 7

#define LOCALSYM_OP_DO_INVERT (true)
#define LOCALSYM_OP_DONT_INVERT (false)

template <typename T>
bool isMultidimArray3DCubic(const MultidimArray<T>& v)
{
	if ( (NSIZE(v) != 1) || (ZSIZE(v) <= 1) || (YSIZE(v) <= 1) || (XSIZE(v) <= 1)
			|| (ZSIZE(v) != YSIZE(v)) || (ZSIZE(v) != XSIZE(v))
			|| (ZSIZE(v) % 2) )
		return false;
	return true;
}

void sum3DCubicMask(
		const MultidimArray<RFLOAT> v,
		RFLOAT& val_sum,
		RFLOAT& val_ctr);

bool similar3DCubicMasks(
		RFLOAT mask1_sum,
		RFLOAT mask1_ctr,
		RFLOAT mask2_sum,
		RFLOAT mask2_ctr);

void truncateMultidimArray(
		MultidimArray<RFLOAT>& v,
		RFLOAT minval = 0.,
		RFLOAT maxval = 0.);

void Localsym_outputOperator(
		const Matrix1D<RFLOAT>& op,
		std::ostream* o_ptr,
		RFLOAT scale_angpix = 1.);

void Localsym_composeOperator(
		Matrix1D<RFLOAT>& op,
		RFLOAT aa = 0., RFLOAT bb = 0., RFLOAT gg = 0.,
		RFLOAT dx = 0., RFLOAT dy = 0., RFLOAT dz = 0.,
		RFLOAT cc = (1e10));

void Localsym_decomposeOperator(
		const Matrix1D<RFLOAT>& op,
		RFLOAT& aa, RFLOAT& bb, RFLOAT& gg,
		RFLOAT& dx, RFLOAT& dy, RFLOAT& dz,
		RFLOAT& cc);

void Localsym_scaleTranslations(
		Matrix1D<RFLOAT>& op,
		RFLOAT factor = 1.);

void Localsym_shiftTranslations(
		Matrix1D<RFLOAT>& op,
		const Matrix1D<RFLOAT>& voffset);

void Localsym_translations2vector(
		const Matrix1D<RFLOAT>& vec,
		Matrix1D<RFLOAT>& trans_vec,
		bool invert = LOCALSYM_OP_DONT_INVERT);

void Localsym_angles2matrix(
		const Matrix1D<RFLOAT>& vec,
		Matrix2D<RFLOAT>& mat,
		bool invert = LOCALSYM_OP_DONT_INVERT);

void Localsym_operator2matrix(
		const Matrix1D<RFLOAT>& vec,
		Matrix2D<RFLOAT>& mat,
		bool invert = LOCALSYM_OP_DONT_INVERT);

void standardiseEulerAngles(
		RFLOAT aa_old, RFLOAT bb_old, RFLOAT gg_old,
		RFLOAT& aa_new, RFLOAT& bb_new, RFLOAT& gg_new);

bool sameLocalsymOperators(
		const Matrix1D<RFLOAT>& lhs,
		const Matrix1D<RFLOAT>& rhs);

void parseDMFormatMasksAndOperators(
		FileName fn_in,
		FileName fn_out);

void readRelionFormatMasksAndOperators(
		FileName fn_info,
		std::vector<FileName>& fn_mask_list,
		std::vector<std::vector<Matrix1D<RFLOAT> > >& ops,
		RFLOAT angpix = 1.,
		bool verb = false);

void readRelionFormatMasksWithoutOperators(
		FileName fn_info,
		std::vector<FileName>& fn_mask_list,
		std::vector<std::vector<Matrix1D<RFLOAT> > >& ops,
		std::vector<std::vector<FileName> >& op_masks,
		bool all_angular_search_ranges_are_global = true,
		bool verb = false);

void writeRelionFormatMasksAndOperators(
		FileName fn_info,
		const std::vector<FileName>& fn_mask_list,
		const std::vector<std::vector<Matrix1D<RFLOAT> > >& ops,
		RFLOAT angpix = 1.);

void writeRelionFormatLocalSearchOperatorResults(
		FileName fn_out,
		const std::vector<Matrix1D<RFLOAT> >& op_samplings,
		RFLOAT angpix = 1.);

void readDMFormatMasksAndOperators(
		FileName fn_info,
		std::vector<FileName>& fn_mask_list,
		std::vector<std::vector<Matrix1D<RFLOAT> > >& op_list,
		RFLOAT angpix = 1.,
		bool verb = false);

void writeDMFormatMasksAndOperators(
		FileName fn_info,
		const std::vector<FileName>& fn_mask_list,
		const std::vector<std::vector<Matrix1D<RFLOAT> > >& ops,
		RFLOAT angpix = 1.);

void duplicateLocalSymmetry(
		MultidimArray<RFLOAT>& out_map,
		const MultidimArray<RFLOAT>& ori_map,
		const std::vector<FileName> fn_masks,
		const std::vector<std::vector<Matrix1D<RFLOAT> > > ops,
		bool duplicate_masks_only = false);

void applyLocalSymmetry(
		MultidimArray<RFLOAT>& sym_map,
		const MultidimArray<RFLOAT>& ori_map,
		const std::vector<FileName> fn_masks,
		const std::vector<std::vector<Matrix1D<RFLOAT> > > ops,
		RFLOAT radius = -1.,
		RFLOAT cosine_width_pix = 5.);

void applyLocalSymmetry(
		MultidimArray<RFLOAT>& map,
		const std::vector<FileName> fn_masks,
		const std::vector<std::vector<Matrix1D<RFLOAT> > > ops,
		RFLOAT radius = -1.,
		RFLOAT cosine_width_pix = 5.);

void getMinCropSize(
		MultidimArray<RFLOAT>& vol,
		Matrix1D<RFLOAT>& center,
		long int& mindim,
		RFLOAT edge = 0.);

bool compareOperatorsByCC(
		const Matrix1D<RFLOAT>& lhs,
		const Matrix1D<RFLOAT>& rhs);

void getLocalSearchOperatorSamplings(
		const Matrix1D<RFLOAT>& op_old,
		const Matrix1D<RFLOAT>& op_search_ranges,
		std::vector<Matrix1D<RFLOAT> >& op_samplings,
		RFLOAT ang_search_step = 1.,
		RFLOAT trans_search_step = 1.,
		bool use_healpix = false,
		bool verb = true);

void calculateOperatorCC(
		const MultidimArray<RFLOAT>& src,
		const MultidimArray<RFLOAT>& dest,
		const MultidimArray<RFLOAT>& mask,
		std::vector<Matrix1D<RFLOAT> >& op_samplings,
		bool do_sort = true,
		bool verb = true);

void separateMasksBFS(
		const FileName& fn_in,
		const int K = 2,
		RFLOAT val_thres = XMIPP_EQUAL_ACCURACY);

/*
void separateMasksKMeans(
		const FileName& fn_in,
		const int K = 2,
		int random_seed = -1);
*/

class local_symmetry_parameters
{
public:
	IOParser parser;

	// Available options
	// PLEASE MAKE SURE THAT ALL THESE OPTIONS ARE INITIALISED IN THE PARSING STEP!
	// ----------------------------------------
	bool show_usage_for_an_option;

	bool do_apply_local_symmetry;
	bool do_duplicate_local_symmetry;
	bool do_local_search_local_symmetry_ops;
	bool do_txt2rln;
	bool do_transform;
	bool do_debug;

	FileName fn_unsym, fn_sym, fn_mask;

	// Input file with mask filenames and rotational / translational operators
	FileName fn_info_in, fn_op_mask_info_in, fn_info_out, fn_info_in_parsed_ext;

	// Manually reset pixel size (in Angstroms)
	RFLOAT angpix_image;

	// Local searches of local symmetry operators
	RFLOAT ang_rot_range, ang_tilt_range, ang_psi_range, ang_range, ang_step;

	RFLOAT offset_x_range, offset_y_range, offset_z_range, offset_range, offset_step;

	RFLOAT rot, tilt, psi, xoff, yoff, zoff;

	RFLOAT binning_factor;

	// Width of soft edge
	RFLOAT width_edge_pix;

	// % of box size as the 2D / 3D spherical mask
	RFLOAT sphere_percentage;

	int nr_masks;

	RFLOAT ini_threshold;

	bool use_healpix_sampling;

	// Verbose output?
	bool verb;

	void initBoolOptions();

	void clear();

	void displayEmptyLine();

	void usage();

	void read(int argc, char **argv);

	void run();

	void writeCommand(FileName fn_cmd, std::string str_executable_name);

	local_symmetry_parameters() { clear(); };

	~local_symmetry_parameters() { clear(); };
};

#endif /* LOCAL_SYMMETRY_H_ */

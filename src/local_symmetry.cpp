/***************************************************************************
 *
 * Author: "Shaoda He"
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

#include "src/local_symmetry.h"

//#define DEBUG
#define NEW_APPLY_SYMMETRY_METHOD

static std::string str_new_mask = "NEW_MASK_AND_OPERATORS";
static std::string str_mask_filename = "MASKFILENAME";

void sum3DCubicMask(
		const MultidimArray<RFLOAT> v,
		RFLOAT& val_sum,
		RFLOAT& val_ctr)
{
	RFLOAT val = 0.;
	val_sum = val_ctr = 0.;

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(v)
	{
		val = DIRECT_A3D_ELEM(v, k, i, j);
		if ( (val < -(XMIPP_EQUAL_ACCURACY)) || ((val - 1.) > (XMIPP_EQUAL_ACCURACY)))
			REPORT_ERROR("ERROR: mask - values are not in range [0,1]!");
		if (val > XMIPP_EQUAL_ACCURACY)
		{
			val_sum += val;
			val_ctr += 1.;
		}
	}

	if ( (val_ctr < 0.9) || (val_sum < 0.01) )
		REPORT_ERROR("ERROR: mask is empty!");
}

bool similar3DCubicMasks(
		RFLOAT mask1_sum,
		RFLOAT mask1_ctr,
		RFLOAT mask2_sum,
		RFLOAT mask2_ctr)
{
	RFLOAT q_sum = 1., q_ctr = 1.;

	if ( (mask1_ctr < 0.9) || (mask1_sum < 0.01) || (mask2_ctr < 0.9) || (mask2_sum < 0.01) )
		REPORT_ERROR("ERROR: mask1 and/or mask2 are empty!");

	q_sum = (mask1_sum > mask2_sum) ? (mask1_sum / mask2_sum) : (mask2_sum / mask1_sum);
	q_ctr = (mask1_ctr > mask2_ctr) ? (mask1_ctr / mask2_ctr) : (mask2_ctr / mask1_ctr);

	if ( (q_sum > 1.1) || (q_ctr > 1.1) )
		return false;
	return true;
}

void truncateMultidimArray(
		MultidimArray<RFLOAT>& v,
		RFLOAT minval,
		RFLOAT maxval)
{
	RFLOAT val = 0.;

	if (minval > maxval)
		REPORT_ERROR("ERROR: minval should be smaller than maxval!");
	/*
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(v)
	{
		val = DIRECT_A3D_ELEM(v, k, i, j);
		if (val < minval)
			DIRECT_A3D_ELEM(v, k, i, j) = minval;
		if (val > maxval)
			DIRECT_A3D_ELEM(v, k, i, j) = maxval;
	}
	*/
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(v)
	{
		val = DIRECT_MULTIDIM_ELEM(v, n);
		if (val < minval)
			DIRECT_MULTIDIM_ELEM(v, n) = minval;
		if (val > maxval)
			DIRECT_MULTIDIM_ELEM(v, n) = maxval;
	}
}

void Localsym_outputOperator(
		const Matrix1D<RFLOAT>& op,
		std::ostream* o_ptr,
		RFLOAT scale_angpix)
{
	if (VEC_XSIZE(op) != NR_LOCALSYM_PARAMETERS)
		REPORT_ERROR("ERROR: op is not a local symmetry operator!");
	if (o_ptr == NULL)
		REPORT_ERROR("ERROR: std::ostream* o_ptr == NULL !");
	if (scale_angpix < 0.001)
		REPORT_ERROR("ERROR: Invalid scale of pixel size!");

	// Enable bold fonts in Unix OS
#ifdef __unix__
	(*o_ptr) << "Angles (rot, tilt, psi) = (" << "\e[1m" << VEC_ELEM(op, AA_POS) << ", " << VEC_ELEM(op, BB_POS) << ", " << VEC_ELEM(op, GG_POS)
			<< "\e[0m" << ") degree(s). Translations (dx, dy, dz) = (" << "\e[1m" << scale_angpix * VEC_ELEM(op, DX_POS) << ", "
			<< scale_angpix * VEC_ELEM(op, DY_POS) << ", " << scale_angpix * VEC_ELEM(op, DZ_POS) << "\e[0m" << ") Angstrom(s)." << std::flush;
#else
	(*o_ptr) << "Angles (rot, tilt, psi) = (" << VEC_ELEM(op, AA_POS) << ", " << VEC_ELEM(op, BB_POS) << ", " << VEC_ELEM(op, GG_POS)
			<< ") degree(s). Translations (dx, dy, dz) = (" << scale_angpix * VEC_ELEM(op, DX_POS) << ", "
			<< scale_angpix * VEC_ELEM(op, DY_POS) << ", " << scale_angpix * VEC_ELEM(op, DZ_POS) << ") Angstrom(s)." << std::flush;
#endif
}

void Localsym_composeOperator(
		Matrix1D<RFLOAT>& op,
		RFLOAT aa, RFLOAT bb, RFLOAT gg,
		RFLOAT dx, RFLOAT dy, RFLOAT dz,
		RFLOAT cc)
{
	op.initZeros(NR_LOCALSYM_PARAMETERS);

	VEC_ELEM(op, AA_POS) = aa; VEC_ELEM(op, BB_POS) = bb; VEC_ELEM(op, GG_POS) = gg;
	VEC_ELEM(op, DX_POS) = dx; VEC_ELEM(op, DY_POS) = dy; VEC_ELEM(op, DZ_POS) = dz;
	VEC_ELEM(op, CC_POS) = cc;
}

void Localsym_decomposeOperator(
		const Matrix1D<RFLOAT>& op,
		RFLOAT& aa, RFLOAT& bb, RFLOAT& gg,
		RFLOAT& dx, RFLOAT& dy, RFLOAT& dz,
		RFLOAT& cc)
{
	aa = bb = gg = dx = dy = dz = 0.;
	cc = (1e10);

	if (VEC_XSIZE(op) != NR_LOCALSYM_PARAMETERS)
		REPORT_ERROR("ERROR: op is not a local symmetry operator!");

	aa = VEC_ELEM(op, AA_POS); bb = VEC_ELEM(op, BB_POS); gg = VEC_ELEM(op, GG_POS);
	dx = VEC_ELEM(op, DX_POS); dy = VEC_ELEM(op, DY_POS); dz = VEC_ELEM(op, DZ_POS);
	cc = VEC_ELEM(op, CC_POS);
}

void Localsym_scaleTranslations(
		Matrix1D<RFLOAT>& op,
		RFLOAT factor)
{
	if (VEC_XSIZE(op) != NR_LOCALSYM_PARAMETERS)
		REPORT_ERROR("ERROR: op is not a local symmetry operator!");

	VEC_ELEM(op, DX_POS) *= factor;
	VEC_ELEM(op, DY_POS) *= factor;
	VEC_ELEM(op, DZ_POS) *= factor;
}

void Localsym_shiftTranslations(
		Matrix1D<RFLOAT>& op,
		const Matrix1D<RFLOAT>& voffset)
{
	if (VEC_XSIZE(op) != NR_LOCALSYM_PARAMETERS)
		REPORT_ERROR("ERROR: op is not a local symmetry operator!");

	if (VEC_XSIZE(voffset) != 3)
		REPORT_ERROR("ERROR: voffset is not a vectorR3!");

	VEC_ELEM(op, DX_POS) += XX(voffset);
	VEC_ELEM(op, DY_POS) += YY(voffset);
	VEC_ELEM(op, DZ_POS) += ZZ(voffset);
}

void Localsym_translations2vector(
		const Matrix1D<RFLOAT>& vec,
		Matrix1D<RFLOAT>& trans_vec,
		bool invert)
{
	trans_vec.clear();

	if (vec.size() != NR_LOCALSYM_PARAMETERS)
		REPORT_ERROR("ERROR: Syntax error in input vector!");

	trans_vec.initZeros(3);
	XX(trans_vec) = VEC_ELEM(vec, DX_POS);
	YY(trans_vec) = VEC_ELEM(vec, DY_POS);
	ZZ(trans_vec) = VEC_ELEM(vec, DZ_POS);

	if (invert == LOCALSYM_OP_DO_INVERT)
	{
		XX(trans_vec) *= -1.;
		YY(trans_vec) *= -1.;
		ZZ(trans_vec) *= -1.;
	}
}

void Localsym_angles2matrix(
		const Matrix1D<RFLOAT>& vec,
		Matrix2D<RFLOAT>& mat,
		bool invert)
{
	RFLOAT aa = 0., bb = 0., gg = 0.;

	mat.clear();

	if (vec.size() != NR_LOCALSYM_PARAMETERS)
		REPORT_ERROR("ERROR: Syntax error in input vector!");

	aa = VEC_ELEM(vec, AA_POS);
	bb = VEC_ELEM(vec, BB_POS);
	gg = VEC_ELEM(vec, GG_POS);
	Euler_angles2matrix(aa, bb, gg, mat);
	if (invert == LOCALSYM_OP_DO_INVERT)
		mat = mat.transpose();
	mat.resize(4, 4);
	MAT_ELEM(mat, 3, 3) = 1.;
}

void Localsym_operator2matrix(
		const Matrix1D<RFLOAT>& vec,
		Matrix2D<RFLOAT>& mat,
		bool invert)
{
	RFLOAT aa = 0., bb = 0., gg = 0.;
	Matrix1D<RFLOAT> trans_vec;

	mat.clear();

	if (vec.size() != NR_LOCALSYM_PARAMETERS)
		REPORT_ERROR("ERROR: Syntax error in input vector!");

	aa = VEC_ELEM(vec, AA_POS);
	bb = VEC_ELEM(vec, BB_POS);
	gg = VEC_ELEM(vec, GG_POS);
	Euler_angles2matrix(aa, bb, gg, mat);

	if (invert == LOCALSYM_OP_DO_INVERT)
	{
		mat = mat.transpose();

		trans_vec.initZeros(3);
		XX(trans_vec)= (-1.) * VEC_ELEM(vec, DX_POS);
		YY(trans_vec)= (-1.) * VEC_ELEM(vec, DY_POS);
		ZZ(trans_vec)= (-1.) * VEC_ELEM(vec, DZ_POS);

		trans_vec = mat * trans_vec;

		mat.resize(4, 4);
		MAT_ELEM(mat, 0, 3) = XX(trans_vec);
		MAT_ELEM(mat, 1, 3) = YY(trans_vec);
		MAT_ELEM(mat, 2, 3) = ZZ(trans_vec);
	}
	else
	{
		mat.resize(4, 4);
		MAT_ELEM(mat, 0, 3) = VEC_ELEM(vec, DX_POS);
		MAT_ELEM(mat, 1, 3) = VEC_ELEM(vec, DY_POS);
		MAT_ELEM(mat, 2, 3) = VEC_ELEM(vec, DZ_POS);
	}

	MAT_ELEM(mat, 3, 3) = 1.;
}

void standardiseEulerAngles(
		RFLOAT aa_old, RFLOAT bb_old, RFLOAT gg_old,
		RFLOAT& aa_new, RFLOAT& bb_new, RFLOAT& gg_new)
{
	Matrix2D<RFLOAT> rot_mat;
	rot_mat.clear();

	// Re-calculate angles so that they follow the conventions in RELION!
	if ( (ABS(aa_old) > 179.) || (bb_old < 1.) || (bb_old > 179.) || (ABS(gg_old) > 179.) )
	{
		Euler_angles2matrix(aa_old, bb_old, gg_old, rot_mat);
		Euler_matrix2angles(rot_mat, aa_new, bb_new, gg_new);
		return;
	}
	aa_new = aa_old; bb_new = bb_old; gg_new = gg_old;
}

bool sameLocalsymOperators(
		const Matrix1D<RFLOAT>& lhs,
		const Matrix1D<RFLOAT>& rhs)
{
	RFLOAT aa1 = 0., bb1 = 0., gg1 = 0., dx1 = 0., dy1 = 0., dz1 = 0., cc1 = 0.;
	RFLOAT aa2 = 0., bb2 = 0., gg2 = 0., dx2 = 0., dy2 = 0., dz2 = 0., cc2 = 0.;
	const RFLOAT eps = (XMIPP_EQUAL_ACCURACY);

	Localsym_decomposeOperator(lhs, aa1, bb1, gg1, dx1, dy1, dz1, cc1);
	Localsym_decomposeOperator(rhs, aa2, bb2, gg2, dx2, dy2, dz2, cc2);

	standardiseEulerAngles(aa1, bb1, gg1, aa1, bb1, gg1);
	standardiseEulerAngles(aa2, bb2, gg2, aa2, bb2, gg2);

	if ( (ABS(aa1 - aa2) < eps) && (ABS(bb1 - bb2) < eps) && (ABS(gg1 - gg2) < eps)
			&& (ABS(dx1 - dx2) < eps) && (ABS(dy1 - dy2) < eps) && (ABS(dz1 - dz2) < eps) )
	{
		return true;
	}
	return false;
}

// Parsing only. Don't validate data here.
void parseDMFormatMasksAndOperators(
		FileName fn_in,
		FileName fn_out)
{
	std::ifstream fin;
	std::ofstream fout;
	std::string line;
	std::vector<std::string> words;

	fin.open(fn_in.c_str(), std::ios_base::in);
	if (fin.fail())
		REPORT_ERROR("ERROR: Cannot open file: " + (std::string)(fn_in));
	//if (exists(fn_out))
	//	REPORT_ERROR("ERROR: output file: " + (std::string)(fn_out) + " exists! Please use another file name!");
	fout.open(fn_out.c_str(), std::ios_base::out);
	if (fout.fail())
		REPORT_ERROR("ERROR: Cannot open file: " + (std::string)(fn_out));

	while (getline(fin, line, '\n'))
	{
		tokenize(line, words);

		// Empty line
		if (words.size() < 1)
			continue;

		// Commented line
		if (words[0][0] == '#')
			continue;

		// Line with mask filename
		if (words[0] == str_mask_filename)
			fout << str_new_mask << std::endl;

		for (int i = 0; i < words.size(); i++)
			fout << words[i] << " " << std::flush;
		fout << std::endl;
	}
	fout.close();
	fin.close();
}

void readRelionFormatMasksAndOperators(
		FileName fn_info,
		std::vector<FileName>& fn_mask_list,
		std::vector<std::vector<Matrix1D<RFLOAT> > >& ops,
		RFLOAT angpix,
		bool verb)
{
	MetaDataTable MD;
	FileName fn_mask;
	std::vector<Matrix1D<RFLOAT> > dummy;
	Matrix1D<RFLOAT> op, op_i;
	bool is_maskname_found = false;
	RFLOAT aa = 0., bb = 0., gg = 0., dx = 0., dy = 0., dz = 0.;

	// Initialisation
	fn_mask_list.clear();
	ops.clear();
	MD.clear();
	dummy.clear();
	op.clear();
	op_i.clear();

	if (angpix < 0.001)
		REPORT_ERROR("ERROR: Pixel size is invalid!");
	if (fn_info.getExtension() != "star")
		REPORT_ERROR("ERROR: " + (std::string)(fn_info) + " is not a STAR file!");
	if (verb)
		std::cout << " Reading list of masks from " << fn_info << "..." << std::endl;

	MD.read(fn_info);
	if (MD.numberOfObjects() < 1)
		REPORT_ERROR("ERROR: STAR file " + (std::string)(fn_info) + " is empty!");
	if ( (!MD.containsLabel(EMDL_MASK_NAME))
			|| (!MD.containsLabel(EMDL_ORIENT_ROT))
			|| (!MD.containsLabel(EMDL_ORIENT_TILT))
			|| (!MD.containsLabel(EMDL_ORIENT_PSI))
			|| (!MD.containsLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM))
			|| (!MD.containsLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM))
			|| (!MD.containsLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM)) )
		REPORT_ERROR("ERROR: You need rlnMaskName, rlnAngleRot, rlnAngleTilt, rlnAnglePsi, rlnOriginXAngst, rlnOriginYAngst and rlnOriginZAngst columns. Some of them are missing in your STAR file " + (std::string)(fn_info) + ". Note that rlnOriginX/Y/Z were changed to rlnOriginX/Y/ZAngst in RELION 3.1. Since the values in the symmetry definition file were in Angstrom from the beginning, please only edit the column names, not values.");

	// Load mask names
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		is_maskname_found = false;
		MD.getValue(EMDL_MASK_NAME, fn_mask);

		for (int id_mask = 0; id_mask < fn_mask_list.size(); id_mask++)
		{
			if (fn_mask_list[id_mask] == fn_mask)
			{
				is_maskname_found = true;
				break;
			}
		}
		if (!is_maskname_found)
			fn_mask_list.push_back(fn_mask);
	}
	if (fn_mask_list.size() < 1)
		REPORT_ERROR("ERROR: No mask filenames in " + (std::string)(fn_info) + " !");

	// Load all operators
	op.initZeros(NR_LOCALSYM_PARAMETERS);
	op_i.initZeros(NR_LOCALSYM_PARAMETERS);
	for (int id_mask = 0; id_mask < fn_mask_list.size(); id_mask++)
	{
		dummy.clear();
		if (verb)
		{
			std::cout << " * Mask #" << (id_mask + 1) << " = " << fn_mask_list[id_mask] << std::endl;
			std::cout << "   --> Operator #" << int(0) << " = " << std::flush;
			Localsym_outputOperator(op_i, &std::cout);
			std::cout << " (the original)" << std::endl;
		}

		// Find all operators for this mask
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			MD.getValue(EMDL_MASK_NAME, fn_mask);
			if (fn_mask != fn_mask_list[id_mask])
				continue;

			// Get this operator
			MD.getValue(EMDL_ORIENT_ROT, aa);
			MD.getValue(EMDL_ORIENT_TILT, bb);
			MD.getValue(EMDL_ORIENT_PSI, gg);
			MD.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, dx);
			MD.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, dy);
			MD.getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, dz);

			// Re-calculate angles so that they follow the conventions in RELION!
			standardiseEulerAngles(aa, bb, gg, aa, bb, gg);
			Localsym_composeOperator(op, aa, bb, gg, dx, dy, dz);

			// Do nothing if it is an identical operator
			if (sameLocalsymOperators(op, op_i))
				continue;

			if (verb)
			{
				std::cout << "   --> Operator #" << (dummy.size() + 1) << " = " << std::flush;
				Localsym_outputOperator(op, &std::cout);
				std::cout << std::endl;
			}

			Localsym_scaleTranslations(op, 1. / angpix);

			// Push back the operator
			dummy.push_back(op);
		}

		if (dummy.size() < 1)
			REPORT_ERROR("ERROR: Please provide at least one non-identical operator for mask file " + fn_mask_list[id_mask] + " !");
		ops.push_back(dummy);
	}

	// Verify mask filenames and operators (detect duplication)
	if ((fn_mask_list.size() < 1) || (ops.size() < 1))
		REPORT_ERROR("ERROR: number of masks and/or operator lists are zero!");
	if (fn_mask_list.size() != ops.size())
		REPORT_ERROR("ERROR: number of masks and operator lists do not match!");
	// Check mask filenames
	for (int imask = 0; imask < fn_mask_list.size() - 1; imask++)
	{
		for (int jmask = imask + 1; jmask < fn_mask_list.size(); jmask++)
		{
			if ( (fn_mask_list[imask].afterLastOf("/").length() > 0)
					&& (fn_mask_list[jmask].afterLastOf("/").length() > 0)
					&& (fn_mask_list[imask].afterLastOf("/") == fn_mask_list[jmask].afterLastOf("/")) )
				REPORT_ERROR("ERROR: Ambiguous mask filenames: " + fn_mask_list[imask] + " and " + fn_mask_list[jmask] + " !");
		}
	}
	// Detect duplicated operators
	// Identical operators have already been removed
	for (int imask = 0; imask < fn_mask_list.size(); imask++)
	{
		for (int iop = 0; iop < ops[imask].size() - 1; iop++)
		{
			for (int jop = iop + 1; jop < ops[imask].size(); jop++)
			{
				if (sameLocalsymOperators(ops[imask][iop], ops[imask][jop]))
					REPORT_ERROR("ERROR: mask filename: " + fn_mask_list[imask] + " contain duplicated operators!");
			}
		}
	}
}

void readRelionFormatMasksWithoutOperators(
		FileName fn_info,
		std::vector<FileName>& fn_mask_list,
		std::vector<std::vector<Matrix1D<RFLOAT> > >& ops,
		std::vector<std::vector<FileName> >& op_masks,
		bool all_angular_search_ranges_are_global,
		bool verb)
{
	FileName fn;
	MetaDataTable MD;
	std::vector<FileName> fns, fns_empty;
	long int id = 0, ide = 0;
	std::vector<long int> ids;
	std::vector<std::vector<FileName> > op_masks_tmp;
	Matrix1D<RFLOAT> op_empty;
	std::vector<Matrix1D<RFLOAT> > ops_empty;

	// Initialisation
	fn_mask_list.clear();
	ops.clear();
	op_masks.clear();

	if (fn_info.getExtension() != "star")
		REPORT_ERROR("ERROR: " + (std::string)(fn_info) + " is not a STAR file!");
	if (verb)
		std::cout << " Reading list of masks from " << fn_info << "..." << std::endl;

	MD.clear();
	MD.read(fn_info);
	if ( (MD.numberOfObjects() < 2) || (MD.numberOfObjects() > 999) )
		REPORT_ERROR("ERROR: STAR file " + (std::string)(fn_info) + " should have 2~999 entries!");
	if ( (!MD.containsLabel(EMDL_MASK_NAME)) || (!MD.containsLabel(EMDL_AREA_ID)) )
		REPORT_ERROR("ERROR: Label EMDL_MASK_NAME and/or EMDL_AREA_ID are missing in STAR file " + (std::string)(fn_info) + " !");
	if (verb)
		std::cout << " Reading list of masks for all operators from " << fn_info << "..." << std::endl;

	// Collect all entries
	fns.clear();
	ids.clear();
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		MD.getValue(EMDL_MASK_NAME, fn);
		MD.getValue(EMDL_AREA_ID, id);
		fns.push_back(fn);
		ids.push_back(id);
		if (id <= 0)
			REPORT_ERROR("ERROR: EMDL_AREA_ID is not a positive integer: " + (std::string)(fn));
	}

	// Check whether there exist duplicated mask filenames
	// fns.size() = id_max
	for (int ii = 0; ii < fns.size() - 1; ii++)
	{
		for (int jj = ii + 1; jj < fns.size(); jj++)
		{
			if (fns[ii] == fns[jj])
				REPORT_ERROR("ERROR: Duplicated mask filenames have been detected: " + fns[ii]);
		}
	}

	// Initialise empty op_masks_tmp[0 ... id_max]
	op_masks_tmp.clear();
	fns_empty.clear();
	for (int ii = 0; ii <= fns.size(); ii++)
		op_masks_tmp.push_back(fns_empty);

	// Collect mask filenames according to their IDs
	for (int ii = 0; ii < fns.size(); ii++)
	{
		// 1 <= ids[ii] <= id_max, op_masks_tmp.size() = id_max + 1
		if (ids[ii] >= op_masks_tmp.size())
			REPORT_ERROR("ERROR: Mask filename contains invalid ID: " + fns[ii]);
		op_masks_tmp[ids[ii]].push_back(fns[ii]);
	}

	// Find the largest area ID (id_e) with mask filenames
	ide = 0;
	// op_masks_tmp.size() = id_max + 1
	for (int ii = op_masks_tmp.size() - 1; ii >= 1; ii--)
	{
		if (op_masks_tmp[ii].size() > 0)
		{
			ide = ii;
			break;
		}
	}
	if (ide <= 0)
		REPORT_ERROR("ERROR: No masks (this should not happen)!");

	// All area IDs 1 ... id_e should be assigned with >= 2 mask filenames
	for (int ii = 1; ii <= ide; ii++)
	{
		if (op_masks_tmp[ii].size() < 2)
			REPORT_ERROR("ERROR: There should be multiple (>= 2) masks for each set of regions!");
	}

	// Input files are valid. Now output arrays for the program.
	fn_mask_list.clear();
	ops.clear();
	op_masks.clear();
	Localsym_composeOperator(op_empty, 0., (all_angular_search_ranges_are_global) ? (90.) : (0.), 0.);
	fns_empty.clear();
	ops_empty.clear();
	for (int ii = 1; ii <= ide; ii++)
	{
		// For each set of N regions:
		// There is a same mask filename in the local symmetry description file
		fn_mask_list.push_back(op_masks_tmp[ii][0]);

		// There is a list of N-1 mask filenames used for global searches
		op_masks.push_back(fns_empty);

		// There is a list of N-1 operators in the local symmetry description file
		ops.push_back(ops_empty);

		// Fill in N-1 mask filenames and N-1 operators into the arrays
		for (int jj = 1; jj < op_masks_tmp[ii].size(); jj++)
		{
			op_masks[op_masks.size() - 1].push_back(op_masks_tmp[ii][jj]);
			ops[ops.size() - 1].push_back(op_empty);
		}
	}

	// Screen output
	if (verb)
	{
		op_empty.initZeros(NR_LOCALSYM_PARAMETERS);

		for (int imask = 0; imask < fn_mask_list.size(); imask++)
		{
			std::cout << " * Mask #" << (imask + 1) << " = " << fn_mask_list[imask] << std::endl;
			std::cout << "   --> Operator #" << int(0) << " = " << std::flush;
			Localsym_outputOperator(op_empty, &std::cout);
			std::cout << " (the original)" << std::endl;

			for (int iop = 0; iop < ops[imask].size(); iop++)
			{
				std::cout << "   --> Operator #" << (iop + 1) << " = " << std::flush;
				Localsym_outputOperator(ops[imask][iop], &std::cout);
				std::cout << " (undefined) - from mask " << op_masks[imask][iop] << std::endl;
			}
		}
	}

	// TODO: this function needs thorough tests!!!
}

void writeRelionFormatMasksAndOperators(
		FileName fn_info,
		const std::vector<FileName>& fn_mask_list,
		const std::vector<std::vector<Matrix1D<RFLOAT> > >& ops,
		RFLOAT angpix)
{
	MetaDataTable MD;

	if (fn_info.getExtension() != "star")
		REPORT_ERROR("ERROR: Output file should have .star extension!");

	if (fn_mask_list.size() != ops.size())
		REPORT_ERROR("ERROR: number of masks and operator lists do not match!");
	if (fn_mask_list.size() < 1)
		REPORT_ERROR("No masks!");

	if (angpix < 0.001)
		REPORT_ERROR("ERROR: Invalid pixel size!");

	MD.clear();
	MD.addLabel(EMDL_MASK_NAME);
	MD.addLabel(EMDL_ORIENT_ROT);
	MD.addLabel(EMDL_ORIENT_TILT);
	MD.addLabel(EMDL_ORIENT_PSI);
	MD.addLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM);
	for (int imask = 0; imask < fn_mask_list.size(); imask++)
	{
		if (ops[imask].size() < 1)
			REPORT_ERROR("ERROR: no operators for mask: " + fn_mask_list[imask]);
		for (int iop = 0; iop < ops[imask].size(); iop++)
		{
			MD.addObject();
			MD.setValue(EMDL_MASK_NAME, fn_mask_list[imask]);
			MD.setValue(EMDL_ORIENT_ROT, VEC_ELEM(ops[imask][iop], AA_POS));
			MD.setValue(EMDL_ORIENT_TILT, VEC_ELEM(ops[imask][iop], BB_POS));
			MD.setValue(EMDL_ORIENT_PSI, VEC_ELEM(ops[imask][iop], GG_POS));
			MD.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, angpix * VEC_ELEM(ops[imask][iop], DX_POS));
			MD.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, angpix * VEC_ELEM(ops[imask][iop], DY_POS));
			MD.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, angpix * VEC_ELEM(ops[imask][iop], DZ_POS));
		}
	}
	MD.write(fn_info);
}

void writeRelionFormatLocalSearchOperatorResults(
		FileName fn_out,
		const std::vector<Matrix1D<RFLOAT> >& op_samplings,
		RFLOAT angpix)
{
	MetaDataTable MD;

	if (angpix < 0.001)
		REPORT_ERROR("ERROR: Invalid pixel size!");
	if (fn_out.getExtension() != "star")
		REPORT_ERROR("ERROR: Output file should have .star extension!");
	if (op_samplings.size() < 1)
		REPORT_ERROR("ERROR: No results!");

	MD.clear();
	MD.addLabel(EMDL_ORIENT_ROT);
	MD.addLabel(EMDL_ORIENT_TILT);
	MD.addLabel(EMDL_ORIENT_PSI);
	MD.addLabel(EMDL_ORIENT_ORIGIN_X_ANGSTROM);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Z_ANGSTROM);
	MD.addLabel(EMDL_IMAGE_WEIGHT);

	for (int iop = 0; iop < op_samplings.size(); iop++)
	{
		if (VEC_XSIZE(op_samplings[iop]) != NR_LOCALSYM_PARAMETERS)
			REPORT_ERROR("ERROR: syntax errors in results!");

		MD.addObject();
		MD.setValue(EMDL_ORIENT_ROT, VEC_ELEM(op_samplings[iop], AA_POS));
		MD.setValue(EMDL_ORIENT_TILT, VEC_ELEM(op_samplings[iop], BB_POS));
		MD.setValue(EMDL_ORIENT_PSI, VEC_ELEM(op_samplings[iop], GG_POS));
		MD.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, angpix * VEC_ELEM(op_samplings[iop], DX_POS));
		MD.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, angpix * VEC_ELEM(op_samplings[iop], DY_POS));
		MD.setValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, angpix * VEC_ELEM(op_samplings[iop], DZ_POS));
		MD.setValue(EMDL_IMAGE_WEIGHT, VEC_ELEM(op_samplings[iop], CC_POS));
	}
	MD.write(fn_out);
}

void readDMFormatMasksAndOperators(FileName fn_info,
		std::vector<FileName>& fn_mask_list,
		std::vector<std::vector<Matrix1D<RFLOAT> > >& op_list,
		RFLOAT angpix,
		bool verb)
{
	// http://www.ccp4.ac.uk/html/rotationmatrices.html

	std::ifstream fin;
	std::string line;
	FileName fn_mask;
	std::vector<std::string> words;

	const std::string str_mask = str_mask_filename;
	const std::string str_rota = "ROTA ";
	const std::string str_euler = "ROTA EULER";
	const std::string str_polar = "ROTA POLAR";
	const std::string str_matrix = "ROTA MATRIX";
	const std::string str_omat = "OMAT";
	const std::string str_trans = "TRAN";
	Matrix1D<RFLOAT> op, op_i;
	std::vector<Matrix1D<RFLOAT> > ops;

	int id_matrix_type = 0;
	RFLOAT a11 = 0., a12 = 0., a13 = 0., a21 = 0., a22 = 0., a23 = 0., a31 = 0., a32 = 0., a33 = 0.;
	RFLOAT dx = 0., dy = 0., dz = 0., aa = 0., bb = 0., gg = 0.;

	// Initialisation
	fn_mask_list.clear();
	op_list.clear();
	op.clear();
	op_i.clear();
	ops.clear();

	op.initZeros(NR_LOCALSYM_PARAMETERS);
	op_i.initZeros(NR_LOCALSYM_PARAMETERS);

	// Open info file
	fin.open(fn_info.c_str(), std::ios_base::in);
	if (fin.fail())
		REPORT_ERROR("ERROR: Cannot open file: " + (std::string)(fn_info));
	if (verb)
		std::cout << " Reading list of masks from " << fn_info << "..." << std::endl;

	while (getline(fin, line, '\n'))
	{
		if (line.find(str_new_mask) != std::string::npos)
			continue;

		// Mask filename is found
		if (line.find(str_mask) != std::string::npos)
		{
			tokenize(line.substr(str_mask.length() + 1), words);
			fn_mask = words[0];
			if (!exists(fn_mask))
				REPORT_ERROR("ERROR: Mask file " + fn_mask + " does not exist!");
			fn_mask_list.push_back(fn_mask);
			if (verb)
			{
				std::cout << " * Mask #" << fn_mask_list.size() << " = " << fn_mask << std::endl;
				std::cout << "   --> Operator #" << int(0) << " = " << std::flush;
				Localsym_outputOperator(op_i, &std::cout);
				std::cout << " (the original)" << std::endl;
			}

			// Get all the operators for this mask
			ops.clear();
			while (getline(fin, line, '\n'))
			{
				if (line.find(str_new_mask) != std::string::npos)
					break;

				id_matrix_type = 0;

				if (line.find(str_rota) == std::string::npos)
				{
					if (line.find(str_omat) == std::string::npos)
						REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
					else
						id_matrix_type += OMAT_TYPE;
				}
				if (line.find(str_euler) != std::string::npos)
					id_matrix_type += ROTA_EULER_TYPE;
				if (line.find(str_polar) != std::string::npos)
					id_matrix_type += ROTA_POLAR_TYPE;
				if (line.find(str_matrix) != std::string::npos)
					id_matrix_type += ROTA_MATRIX_TYPE;

				if ((id_matrix_type != ROTA_EULER_TYPE)
						&& (id_matrix_type != ROTA_POLAR_TYPE)
						&& (id_matrix_type != ROTA_MATRIX_TYPE)
						&& (id_matrix_type != OMAT_TYPE))
					REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);

				if (id_matrix_type == ROTA_EULER_TYPE)
					tokenize(line.substr(str_euler.length() + 1), words);
				else if (id_matrix_type == ROTA_POLAR_TYPE)
					tokenize(line.substr(str_polar.length() + 1), words);
				else if (id_matrix_type == ROTA_MATRIX_TYPE)
					tokenize(line.substr(str_matrix.length() + 1), words);
				else if (id_matrix_type == OMAT_TYPE)
				{
					if (!getline(fin, line, '\n')) // Read a new line
						REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
					tokenize(line, words);
				}
				if (words.size() < 3)
					REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);

				a11 = textToFloat(words[0]);
				a12 = textToFloat(words[1]);
				a13 = textToFloat(words[2]);
				if ((id_matrix_type == ROTA_MATRIX_TYPE) || (id_matrix_type == OMAT_TYPE))
				{
					if (getline(fin, line, '\n'))
					{
						tokenize(line, words);
						if (words.size() < 3)
							REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
						a21 = textToFloat(words[0]);
						a22 = textToFloat(words[1]);
						a23 = textToFloat(words[2]);
					}
					else
						REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);

					if (getline(fin, line, '\n'))
					{
						tokenize(line, words);
						if (words.size() < 3)
							REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
						a31 = textToFloat(words[0]);
						a32 = textToFloat(words[1]);
						a33 = textToFloat(words[2]);
					}
					else
						REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
				}
				if (id_matrix_type == OMAT_TYPE)
				{
					if (getline(fin, line, '\n'))
					{
						tokenize(line, words);
						if (words.size() < 3)
							REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
						dx = textToFloat(words[0]);
						dy = textToFloat(words[1]);
						dz = textToFloat(words[2]);
					}
					else
						REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
				}

				if (id_matrix_type == ROTA_EULER_TYPE)
				{
					standardiseEulerAngles(a11, a12, a13, aa, bb, gg);
				}
				else if (id_matrix_type == ROTA_POLAR_TYPE)
				{
					// omega, phi, kappa
					RFLOAT omega = a11, phi = a12, kappa = a13;
					RFLOAT ll = sin(DEG2RAD(omega)) * cos(DEG2RAD(phi));
					RFLOAT mm = sin(DEG2RAD(omega)) * sin(DEG2RAD(phi));
					RFLOAT nn = cos(DEG2RAD(omega));
					RFLOAT ck = cos(DEG2RAD(kappa));
					RFLOAT sk = sin(DEG2RAD(kappa));

					a11 = ll * ll + (mm * mm + nn * nn) * ck;
					a12 = ll * mm * (1. - ck) - nn * sk;
					a13 = nn * ll * (1. - ck) + mm * sk;
					a21 = ll * mm * (1. - ck) + nn * sk;
					a22 = mm * mm + (ll * ll + nn * nn) * ck;
					a23 = mm * nn * (1. - ck) - ll * sk;
					a31 = nn * ll * (1. - ck) - mm * sk;
					a32 = mm * nn * (1. - ck) + ll * sk;
					a33 = nn * nn + (ll * ll + mm * mm) * ck;

				}
				// These three type of operators contain angular matrices
				if ((id_matrix_type == ROTA_POLAR_TYPE)
						|| (id_matrix_type == ROTA_MATRIX_TYPE)
						|| (id_matrix_type == OMAT_TYPE))
				{
					Matrix2D<RFLOAT> A;
					A.resize(3, 3);

					MAT_ELEM(A, 0, 0) = a11; MAT_ELEM(A, 0, 1) = a12; MAT_ELEM(A, 0, 2) = a13;
					MAT_ELEM(A, 1, 0) = a21; MAT_ELEM(A, 1, 1) = a22; MAT_ELEM(A, 1, 2) = a23;
					MAT_ELEM(A, 2, 0) = a31; MAT_ELEM(A, 2, 1) = a32; MAT_ELEM(A, 2, 2) = a33;

					Euler_matrix2angles(A.transpose(), aa, bb, gg); // TODO: do we need transpose here?
				}

				// Read TRANS
				if ((id_matrix_type == ROTA_EULER_TYPE)
						|| (id_matrix_type == ROTA_POLAR_TYPE)
						|| (id_matrix_type == ROTA_MATRIX_TYPE))
				{
					if (getline(fin, line, '\n') && (line.find(str_trans) != std::string::npos))
					{
						tokenize(line.substr(str_trans.length() + 1), words);
						if (words.size() < 3)
							REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
						dx = textToFloat(words[0]);
						dy = textToFloat(words[1]);
						dz = textToFloat(words[2]);
					}
					else
						REPORT_ERROR("ERROR: Syntax error: Operators of mask file " + fn_mask);
				}

				// New matrix has been processed
				Localsym_composeOperator(op, aa, bb, gg, dx, dy, dz);

				// Check whether it is an identical operator, then push back
				if (sameLocalsymOperators(op, op_i))
					continue;

				if (verb)
				{
					std::cout << "   --> Operator #" << (ops.size() + 1) << " = " << std::flush;
					Localsym_outputOperator(op, &std::cout);
					std::cout << std::endl;
				}
				Localsym_scaleTranslations(op, 1. / angpix);
				ops.push_back(op);

			}
			// All the operators for this mask are read
			if (ops.size() < 1)
				REPORT_ERROR("ERROR: Please provide at least one non-identical operator for mask file " + fn_mask + " !");
			op_list.push_back(ops);
		}
		else
		{
			// Mask filename is not found
			REPORT_ERROR("ERROR: Syntax error: mask filename is not found! " + fn_info);
		}
	}

	// Verify mask filenames and operators (detect duplication)
	if ((fn_mask_list.size() < 1) || (op_list.size() < 1))
		REPORT_ERROR("ERROR: number of masks and/or operator lists are zero!");
	if (fn_mask_list.size() != op_list.size())
		REPORT_ERROR("ERROR: number of masks and operator lists do not match!");
	// Check mask filenames
	for (int imask = 0; imask < fn_mask_list.size() - 1; imask++)
	{
		for (int jmask = imask + 1; jmask < fn_mask_list.size(); jmask++)
		{
			if ( (fn_mask_list[imask].afterLastOf("/").length() > 0)
					&& (fn_mask_list[jmask].afterLastOf("/").length() > 0)
					&& (fn_mask_list[imask].afterLastOf("/") == fn_mask_list[jmask].afterLastOf("/")) )
				REPORT_ERROR("ERROR: Ambiguous mask filenames: " + fn_mask_list[imask] + " and " + fn_mask_list[jmask] + " !");
		}
	}
	// Detect duplicated operators
	// Identical operators have already been removed
	for (int imask = 0; imask < fn_mask_list.size(); imask++)
	{
		for (int iop = 0; iop < op_list[imask].size() - 1; iop++)
		{
			for (int jop = iop + 1; jop < op_list[imask].size(); jop++)
			{
				if (sameLocalsymOperators(op_list[imask][iop], op_list[imask][jop]))
					REPORT_ERROR("ERROR: mask filename: " + fn_mask_list[imask] + " contain duplicated operators!");
			}
		}
	}

#ifdef DEBUG
	for (int imask = 0; imask < fn_mask_list.size(); imask++)
	{
		std::cout << " * Mask #" << (imask + 1) << " = " << fn_mask_list[imask] << std::endl;
		for (int iop = 0; iop < op_list[imask].size(); iop++)
		{
			std::cout << "   --> Operator #" << (iop + 1) << " = "
					<< VEC_ELEM(op_list[imask][iop], AA_POS) << ", " << VEC_ELEM(op_list[imask][iop], BB_POS) << ", " << VEC_ELEM(op_list[imask][iop], GG_POS) << "; "
					<< VEC_ELEM(op_list[imask][iop], DX_POS) << ", " << VEC_ELEM(op_list[imask][iop], DY_POS) << ", " << VEC_ELEM(op_list[imask][iop], DZ_POS) << std::endl;
		}
	}
#endif
}

void writeDMFormatMasksAndOperators(
		FileName fn_info,
		const std::vector<FileName>& fn_mask_list,
		const std::vector<std::vector<Matrix1D<RFLOAT> > >& ops,
		RFLOAT angpix)
{
	if (fn_info.getExtension() == "star")
		REPORT_ERROR("ERROR: Output file should not have .star extension!");

	if (fn_mask_list.size() != ops.size())
		REPORT_ERROR("ERROR: number of masks and operator lists do not match!");
	if (fn_mask_list.size() < 1)
		REPORT_ERROR("No masks!");

	if (angpix < 0.001)
		REPORT_ERROR("ERROR: Invalid pixel size!");

	for (int imask = 0; imask < fn_mask_list.size(); imask++)
	{
		if (ops[imask].size() < 1)
			REPORT_ERROR("ERROR: no operators for mask: " + fn_mask_list[imask]);
	}

	int str_w = 15;
	std::ofstream fout;
	fout.open(fn_info.c_str(), std::ios::out);
	if (!fout)
		REPORT_ERROR("ERROR: Cannot write to file: " + fn_info);
	for (int imask = 0; imask < fn_mask_list.size(); imask++)
	{
		fout << std::endl << str_mask_filename << " " << fn_mask_list[imask] << std::endl;
		for (int iop = 0; iop < ops[imask].size(); iop++)
		{
			fout << " ROTA EULER " << std::setiosflags(std::ios::fixed)
					<< std::setw(str_w) << VEC_ELEM(ops[imask][iop], AA_POS) << " "
					<< std::setw(str_w) << VEC_ELEM(ops[imask][iop], BB_POS) << " "
					<< std::setw(str_w) << VEC_ELEM(ops[imask][iop], GG_POS)
					<< std::resetiosflags(std::ios::fixed) << std::endl;
			fout << "       TRAN " << std::setiosflags(std::ios::fixed)
					<< std::setw(str_w) << angpix * VEC_ELEM(ops[imask][iop], DX_POS) << " "
					<< std::setw(str_w) << angpix * VEC_ELEM(ops[imask][iop], DY_POS) << " "
					<< std::setw(str_w) << angpix * VEC_ELEM(ops[imask][iop], DZ_POS)
					<< std::resetiosflags(std::ios::fixed) << std::endl;
		}
	}
	fout.close();
}

void duplicateLocalSymmetry(
		MultidimArray<RFLOAT>& out_map,
		const MultidimArray<RFLOAT>& ori_map,
		const std::vector<FileName> fn_masks,
		const std::vector<std::vector<Matrix1D<RFLOAT> > > ops,
		bool duplicate_masks_only)
{
	Image<RFLOAT> mask;
	MultidimArray<RFLOAT> vol1, ori_map_masked;
	Matrix1D<RFLOAT> trans_vec;
	Matrix2D<RFLOAT> op_mat;

	out_map.clear();

	if ((fn_masks.size() < 1) || (ops.size() < 1))
		REPORT_ERROR("ERROR: number of masks and/or operator lists are zero!");
	if (fn_masks.size() != ops.size())
		REPORT_ERROR("ERROR: number of masks and operator lists do not match!");

	// Open the first mask header, or copy original map for initialisation of output map
	if (duplicate_masks_only)
	{
		if (!exists(fn_masks[0]))
			REPORT_ERROR("ERROR: mask " + std::string(fn_masks[0]) + " does not exist!");
		mask.read(fn_masks[0], false);
		if ((NSIZE(mask()) != 1) || (ZSIZE(mask()) <= 1) || (YSIZE(mask()) <= 1) || (XSIZE(mask()) <= 1))
			REPORT_ERROR("ERROR: input mask is not 3D!");

		out_map.initZeros(mask());
	}
	else
		out_map.initZeros(ori_map);
	vol1.clear();
	ori_map_masked.clear();

	// Loop over all masks
	for (int imask = 0; imask < fn_masks.size(); imask++)
	{
		// Load this mask
		if (!exists(fn_masks[imask]))
			REPORT_ERROR("ERROR: mask " + std::string(fn_masks[imask]) + " does not exist!");
		mask.clear();
		mask.read(fn_masks[imask]);
		if ((NSIZE(out_map) != NSIZE(mask())) || (ZSIZE(out_map) != ZSIZE(mask())) || (YSIZE(out_map) != YSIZE(mask())) || (XSIZE(out_map) != XSIZE(mask())))
			REPORT_ERROR("ERROR: All masks (and input map) should have the same sizes!");
		// Masks and the original map may not have the same origin!
		mask().copyShape(out_map); // VERY IMPORTANT!

		// Add this mask (or masked original map) to final result
		if (duplicate_masks_only)
			out_map += mask();
		else
		{
			ori_map_masked = ori_map * mask();
			out_map += ori_map_masked;
		}

		// Loop over all operators for this mask
		if (ops[imask].size() < 1)
			REPORT_ERROR("ERROR: number of operators for mask " + std::string(fn_masks[imask]) + " is less than 1!");
		for (int iop = 0; iop < ops[imask].size(); iop++)
		{
#ifdef NEW_APPLY_SYMMETRY_METHOD
			Localsym_operator2matrix(ops[imask][iop], op_mat);

			if (duplicate_masks_only)
				applyGeometry(mask(), vol1, op_mat, IS_NOT_INV, DONT_WRAP);
			else
				applyGeometry(ori_map_masked, vol1, op_mat, IS_NOT_INV, DONT_WRAP);
#else
			Localsym_angles2matrix(ops[imask][iop], op_mat);
			Localsym_translations2vector(ops[imask][iop], trans_vec);

			if (duplicate_masks_only)
				applyGeometry(mask(), vol1, op_mat, IS_NOT_INV, DONT_WRAP);
			else
				applyGeometry(ori_map_masked, vol1, op_mat, IS_NOT_INV, DONT_WRAP);
			selfTranslate(vol1, trans_vec, DONT_WRAP);
#endif
			out_map += vol1;
		}
	}
}

void applyLocalSymmetry(MultidimArray<RFLOAT>& sym_map,
		const MultidimArray<RFLOAT>& ori_map,
		const std::vector<FileName> fn_masks,
		const std::vector<std::vector<Matrix1D<RFLOAT> > > ops,
		RFLOAT radius,
		RFLOAT cosine_width_pix)
{
	MultidimArray<RFLOAT> w, vol1, vol2;
	Image<RFLOAT> mask;
	Matrix1D<RFLOAT> trans_vec;
	Matrix2D<RFLOAT> op_mat;
	RFLOAT mask_val = 0., sym_val = 0., radius2 = 0., radiusw2 = 0., dist2 = 0., xinit = 0., yinit = 0., zinit = 0.;

	// Initialise the result
	sym_map.clear();

	if ((NSIZE(ori_map) != 1) || (ZSIZE(ori_map) <= 1) || (YSIZE(ori_map) <= 1) || (XSIZE(ori_map) <= 1))
		REPORT_ERROR("ERROR: input unsymmetrised map is not 3D!");
	// Support 3D maps which are not cubic
	// Support 3D maps and masks which do not share the same origins

	if ( (radius > 0.) && (cosine_width_pix < (XMIPP_EQUAL_ACCURACY)) )
		REPORT_ERROR("ERROR: Cosine width should be larger than 0!");

	if ((fn_masks.size() < 1) || (ops.size() < 1))
		REPORT_ERROR("ERROR: number of masks and/or operator lists are zero!");

	if (fn_masks.size() != ops.size())
		REPORT_ERROR("ERROR: number of masks and operator lists do not match!");

	sym_map.initZeros(ori_map);
	w.initZeros(ori_map);
	vol1.clear();
	vol2.clear(); // Use vol2 only as the output from 'applyGeometry()'

	// Loop over all the masks
	for (int imask = 0; imask < fn_masks.size(); imask++)
	{
		vol1 = ori_map;

		// Loop over all operators for this mask
		RFLOAT nr_ops = RFLOAT(ops[imask].size());
		if (nr_ops < 0.9)
			REPORT_ERROR("ERROR: number of operators for mask " + std::string(fn_masks[imask]) + " is less than 1!");
		for (int iop = 0; iop < ops[imask].size(); iop++)
		{
			// A0 op(--rot--> --trans-->) A1
			// Now we want A1 op'(?) A0
			// Transform A1 to A0 (original) position, then superimpose with the original
#ifdef NEW_APPLY_SYMMETRY_METHOD
			Localsym_operator2matrix(ops[imask][iop], op_mat, LOCALSYM_OP_DO_INVERT);

			applyGeometry(ori_map, vol2, op_mat, IS_NOT_INV, DONT_WRAP);
#else
			Localsym_translations2vector(ops[imask][iop], trans_vec, LOCALSYM_OP_DO_INVERT);
			Localsym_angles2matrix(ops[imask][iop], op_mat, LOCALSYM_OP_DO_INVERT);

			translate(ori_map, vol2, trans_vec, DONT_WRAP);
			selfApplyGeometry(vol2, op_mat, IS_NOT_INV, DONT_WRAP);
#endif
			vol1 += vol2;
		}

		// Load this mask
		if (!exists(fn_masks[imask]))
			REPORT_ERROR("ERROR: mask " + std::string(fn_masks[imask]) + " does not exist!");
		mask.clear();
		mask.read(fn_masks[imask]);
		if ((NSIZE(ori_map) != NSIZE(mask())) || (ZSIZE(ori_map) != ZSIZE(mask())) || (YSIZE(ori_map) != YSIZE(mask())) || (XSIZE(ori_map) != XSIZE(mask())))
			REPORT_ERROR("ERROR: sizes of input and masks do not match!");
		// Masks and the original map may not have the same origin!
		mask().copyShape(ori_map); // VERY IMPORTANT!

		// 'Vol1' contains one symmetrised subunit, make it into perfect "mask-weighted sum"
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(vol1)
		{
			// Get values of mask in every voxel
			mask_val = DIRECT_A3D_ELEM(mask(), k, i, j); // "weights from mask" - w
			if ((mask_val < -(XMIPP_EQUAL_ACCURACY)) || ((mask_val - 1.) > (XMIPP_EQUAL_ACCURACY)))
				REPORT_ERROR("ERROR: mask " + std::string(fn_masks[imask]) + " - values are not in range [0,1]!");

			// This voxel is inside the mask
			if (mask_val > (XMIPP_EQUAL_ACCURACY))
			{
				DIRECT_A3D_ELEM(vol1, k, i, j) *= mask_val / (nr_ops + 1.); // "mask-weighted sum" - wsum
			}
			else
			{
				// This voxel is not inside the mask
				DIRECT_A3D_ELEM(vol1, k, i, j) = 0.;
			}
		}

		// Make various copies of vol1 and mask to wsum and w
		sym_map += vol1;
		w += mask();
		for (int iop = 0; iop < ops[imask].size(); iop++)
		{
#ifdef NEW_APPLY_SYMMETRY_METHOD
			Localsym_operator2matrix(ops[imask][iop], op_mat);

			applyGeometry(vol1, vol2, op_mat, IS_NOT_INV, DONT_WRAP);
			sym_map += vol2;
			applyGeometry(mask(), vol2, op_mat, IS_NOT_INV, DONT_WRAP);
			w += vol2;
#else
			Localsym_angles2matrix(ops[imask][iop], op_mat);
			Localsym_translations2vector(ops[imask][iop], trans_vec);

			applyGeometry(vol1, vol2, op_mat, IS_NOT_INV, DONT_WRAP);
			selfTranslate(vol2, trans_vec, DONT_WRAP);
			sym_map += vol2;
			applyGeometry(mask(), vol2, op_mat, IS_NOT_INV, DONT_WRAP);
			selfTranslate(vol2, trans_vec, DONT_WRAP);
			w += vol2;
#endif
		}

		// Unload this mask
		mask.clear();
	}
	vol1.clear();
	vol2.clear();
	mask.clear();

	// TODO: check! please always ensure - free memory space in time!

	// sym_map and w contain all symmetised subunits (wsum) and mask coefficients (w) needed
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(sym_map)
	{
		mask_val = DIRECT_A3D_ELEM(w, k, i, j); // get weights

		// TODO: check radius2 here!

		// This voxel is inside one of the masks
		if (mask_val > (XMIPP_EQUAL_ACCURACY)) // weight > 0
		{
			if ((mask_val - 1.) > (XMIPP_EQUAL_ACCURACY)) // weight > 1
			{
				// ncs = wsum / w
				DIRECT_A3D_ELEM(sym_map, k, i, j) /= mask_val;
			}
			else if ((mask_val - 1.) < (-(XMIPP_EQUAL_ACCURACY))) // 0 < weight < 1
			{
				// ncs = w * (wsum / w) + (1 - w) * ori_val
				sym_val = DIRECT_A3D_ELEM(sym_map, k, i, j);
				DIRECT_A3D_ELEM(sym_map, k, i, j) = sym_val + (1. - mask_val) * DIRECT_A3D_ELEM(ori_map, k, i, j);
			}
			// weight = 1, ncs = wsum / w, nothing to do...
		}
		else
		{
			// weight <= 0, ncs = ori_val
			DIRECT_A3D_ELEM(sym_map, k, i, j) = DIRECT_A3D_ELEM(ori_map, k, i, j);
		}
	}

	if (radius > 0.)
	{
		radius2 = radius * radius;
		radiusw2 = (radius + cosine_width_pix) * (radius + cosine_width_pix);
		xinit = FIRST_XMIPP_INDEX(XSIZE(sym_map));
		yinit = FIRST_XMIPP_INDEX(YSIZE(sym_map));
		zinit = FIRST_XMIPP_INDEX(ZSIZE(sym_map));

		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(sym_map)
		{
			dist2 = (k + zinit) * (k + zinit) + (i + yinit) * (i + yinit) + (j + xinit) * (j + xinit);
			if (dist2 > radiusw2)
				DIRECT_A3D_ELEM(sym_map, k, i, j) = 0.;
			else if (dist2 > radius2)
				DIRECT_A3D_ELEM(sym_map, k, i, j) *= 0.5 + 0.5 * cos(PI * (radius + cosine_width_pix - sqrt(dist2)) / cosine_width_pix);
		}
	}

	// Done!
}

void applyLocalSymmetry(
		MultidimArray<RFLOAT>& map,
		const std::vector<FileName> fn_masks,
		const std::vector<std::vector<Matrix1D<RFLOAT> > > ops,
		RFLOAT radius,
		RFLOAT cosine_width_pix)
{
	MultidimArray<RFLOAT> vol;
	applyLocalSymmetry(vol, map, fn_masks, ops, radius, cosine_width_pix);
	map = vol;
}

void getMinCropSize(
		MultidimArray<RFLOAT>& vol,
		Matrix1D<RFLOAT>& center,
		long int& mindim,
		RFLOAT edge)
{
	RFLOAT val = 0., dist2 = 0., dist2_max = 0.;
	RFLOAT xori = 0., yori = 0., zori = 0.;
	Matrix1D<RFLOAT> new_center;

	mindim = -1;
	center.initZeros(3);
	new_center.initZeros(3);

	if ((NSIZE(vol) != 1) || (ZSIZE(vol) <= 1) || (YSIZE(vol) <= 1) || (XSIZE(vol) <= 1))
		REPORT_ERROR("ERROR: input mask is not 3D!");

	vol.setXmippOrigin();
	vol.centerOfMass(center);
	xori = XX(center); yori = YY(center); zori = ZZ(center);

	dist2_max = -999.;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(vol)
	{
		val = A3D_ELEM(vol, k, i, j);

		if (val < -(XMIPP_EQUAL_ACCURACY))
			REPORT_ERROR("ERROR: all voxels in the input map should have positive values!");

		if (val > (XMIPP_EQUAL_ACCURACY))
		{
			dist2 = (RFLOAT(k) - zori) * (RFLOAT(k) - zori) + (RFLOAT(i) - yori) * (RFLOAT(i) - yori) + (RFLOAT(j) - xori) * (RFLOAT(j) - xori);
			if (dist2 > dist2_max)
				dist2_max = dist2;
		}
	}

	if (dist2_max < 0.)
		REPORT_ERROR("ERROR: the input map is empty!");
	if (dist2_max > 99999999. * 99999999.)
		REPORT_ERROR("ERROR: size of the input map is too large (> 99999999)!");

	dist2_max = sqrt(dist2_max);
	if (edge > 0.)
		dist2_max += edge;
	mindim = 2 * (long int)(ceil(dist2_max)); // bestdim % 2 = 0
}

bool compareOperatorsByCC(
		const Matrix1D<RFLOAT>& lhs,
		const Matrix1D<RFLOAT>& rhs)
{
	return (VEC_ELEM(lhs, CC_POS) < VEC_ELEM(rhs, CC_POS));
}

void getLocalSearchOperatorSamplings(
		const Matrix1D<RFLOAT>& op_old,
		const Matrix1D<RFLOAT>& op_search_ranges,
		std::vector<Matrix1D<RFLOAT> >& op_samplings,
		RFLOAT ang_search_step,
		RFLOAT trans_search_step,
		bool use_healpix,
		bool verb)
{
	RFLOAT aa = 0., bb = 0., gg = 0., dx = 0., dy = 0., dz = 0., cc = 0.;
	RFLOAT aa_range = 0., bb_range = 0., gg_range = 0., dx_range = 0., dy_range = 0., dz_range = 0.;
	RFLOAT aa_residue = 0., bb_residue = 0., gg_residue = 0., dx_residue = 0., dy_residue = 0., dz_residue = 0.;
	RFLOAT aa_init = 0., bb_init = 0., gg_init = 0., dx_init = 0., dy_init = 0., dz_init = 0.;;
	RFLOAT val = 0., r2 = 0.;
	long int nr_dir = 0, nr_all_samplings = 0;

	std::vector<RFLOAT> aas, bbs, ggs, dxs, dys, dzs;
	Matrix1D<RFLOAT> op_tmp;
	Matrix2D<RFLOAT> op_mat;

	HealpixSampling sampling;
	std::vector<int> pointer_dir_nonzeroprior, pointer_psi_nonzeroprior;

	op_samplings.clear();

	op_tmp.clear();
	op_mat.clear();
	sampling.clear();
	pointer_dir_nonzeroprior.clear(); pointer_psi_nonzeroprior.clear();

	if ( (VEC_XSIZE(op_old) != NR_LOCALSYM_PARAMETERS) || (VEC_XSIZE(op_search_ranges) != NR_LOCALSYM_PARAMETERS) )
		REPORT_ERROR("ERROR: Input operator contains syntax error!");

	if ( (ang_search_step < 0.0001) || (ang_search_step > 30.) )
		REPORT_ERROR("ERROR: Angular searching step should be within range (+0.0001, +30.0000) degrees!");
	if ( (trans_search_step < 0.0001) || (trans_search_step > 5.) )
		REPORT_ERROR("ERROR: Translational searching step should be within range (+0.0001, +5.0000) rescaled / binned pixels!");

	Localsym_decomposeOperator(op_old, aa_init, bb_init, gg_init, dx_init, dy_init, dz_init, cc);
	Localsym_decomposeOperator(op_search_ranges, aa_range, bb_range, gg_range, dx_range, dy_range, dz_range, cc);

	// Angular searching ranges
	if (!use_healpix)
	{
		//aa_range = ( (aa_range > 180.) || (aa_range < 0.) ) ? (180.) : aa_range;
		//bb_range = ( (bb_range >  90.) || (bb_range < 0.) ) ? ( 90.) : bb_range;
		//gg_range = ( (gg_range > 180.) || (gg_range < 0.) ) ? (180.) : gg_range;

		aa_range = (aa_range > 180.) ? (180.) : aa_range;
		bb_range = (bb_range >  90.) ? ( 90.) : bb_range;
		gg_range = (gg_range > 180.) ? (180.) : gg_range;
		aa_range = (aa_range > 0.) ? aa_range : 0.;
		bb_range = (bb_range > 0.) ? bb_range : 0.;
		gg_range = (gg_range > 0.) ? gg_range : 0.;
	}
	if ( ( (aa_range < ang_search_step) && (aa_range > XMIPP_EQUAL_ACCURACY) )
			|| ( (bb_range < ang_search_step) && (bb_range > XMIPP_EQUAL_ACCURACY) )
			|| ( (gg_range < ang_search_step) && (gg_range > XMIPP_EQUAL_ACCURACY) ) )
		REPORT_ERROR("ERROR: Angular searching step should be smaller than its searching range!");
	if (!use_healpix)
	{
		// aa, bb, gg ranges >= 0, ang_search_step > 0.01
		aa_residue = aa_range - ang_search_step * floor(aa_range / ang_search_step);
		bb_residue = bb_range - ang_search_step * floor(bb_range / ang_search_step);
		gg_residue = gg_range - ang_search_step * floor(gg_range / ang_search_step);
	}

	// Translational searching ranges
	dx_range = (dx_range > 0.) ? dx_range : 0.;
	dy_range = (dy_range > 0.) ? dy_range : 0.;
	dz_range = (dz_range > 0.) ? dz_range : 0.;
	if ( ( (dx_range < trans_search_step) && (dx_range > XMIPP_EQUAL_ACCURACY) )
			|| ( (dy_range < trans_search_step) && (dy_range > XMIPP_EQUAL_ACCURACY) )
			|| ( (dz_range < trans_search_step) && (dz_range > XMIPP_EQUAL_ACCURACY) ) )
		REPORT_ERROR("ERROR: Translational searching step should be smaller than its searching range!");
	// dx, dy, dz ranges >= 0, ang_search_step > 0.01
	dx_residue = dx_range - trans_search_step * floor(dx_range / trans_search_step);
	dy_residue = dy_range - trans_search_step * floor(dy_range / trans_search_step);
	dz_residue = dz_range - trans_search_step * floor(dz_range / trans_search_step);

	if (verb)
	{
		//std::cout << " + Local searches of local symmetry operator: Angles (rot, tilt, psi) = ("
		//		<< aa_init << ", " << bb_init << ", " << gg_init << ") degree(s), center of mass (x, y, z; cropped, rescaled, binned) = ("
		//		<< dx_init << ", " << dy_init << ", " << dz_init << ") pixel(s)..." << std::endl;
		std::cout << " + Generating sampling points with ranges: Angles (rot, tilt, psi) = +/- ("
				<< aa_range << ", " << bb_range << ", " << gg_range << ") degree(s), center of mass (x, y, z; cropped, rescaled, binned) = +/- ("
				<< dx_range << ", " << dy_range << ", " << dz_range << ") pixel(s)." << std::endl;
		std::cout << " + Generating sampling points with step sizes: " << ang_search_step << " degree(s), "
				<< trans_search_step << " rescaled (binned) pixel(s)." << std::endl;
	}

	// Angular samplings
	if (use_healpix)
	{
		std::vector<RFLOAT> dummy1, dummy2;

		pointer_dir_nonzeroprior.clear(); pointer_psi_nonzeroprior.clear();
		dummy1.clear(); dummy2.clear();

		// Get healpix order and mode
		int healpix_order = 0, prior_mode = 0;
		for (healpix_order = 0; healpix_order <= 100; healpix_order++)
		{
			if (ang_search_step > (360. / (6. * ROUND(std::pow(2., healpix_order)))) )
				break;
		}
		if (healpix_order >= 100)
			REPORT_ERROR("ERROR: healpix_order is larger than 100!");
		prior_mode = ((aa_range < 0.) && (bb_range < 0) && (gg_range < 0.)) ? (NOPRIOR) : (PRIOR_ROTTILT_PSI);

		// Initialise healpix sampling
		sampling.clear();
		sampling.healpix_order = healpix_order;
		sampling.is_3D = sampling.is_3d_trans = true;
		sampling.limit_tilt = -91.; // Don't limit tilts
		sampling.psi_step = 360. / (6. * ROUND(std::pow(2., healpix_order)));
		sampling.offset_range = sampling.offset_step = 1.; // I don't use Healpix translational samplings
		sampling.random_perturbation = sampling.perturbation_factor = 0.;

		// Get all orientations
		sampling.initialise(3, true, false, false, (prior_mode == NOPRIOR) ? (false) : (true));
		sampling.setOrientations();

		// Select orientations
		if (prior_mode == PRIOR_ROTTILT_PSI)
		{
			sampling.selectOrientationsWithNonZeroPriorProbability(
					aa_init, bb_init, gg_init, aa_range, bb_range, gg_range,
		    		pointer_dir_nonzeroprior, dummy1, pointer_psi_nonzeroprior, dummy2, false, 1.);
		}
		else
		{
			// Just push all directions
			for (int idir = 0; idir < sampling.rot_angles.size(); idir++)
				pointer_dir_nonzeroprior.push_back(idir);
			for (int ipsi = 0; ipsi < sampling.psi_angles.size(); ipsi++)
				pointer_psi_nonzeroprior.push_back(ipsi);
		}

		if ( (sampling.rot_angles.size() < 1) || (sampling.tilt_angles.size() < 1) || (sampling.psi_angles.size() < 1)
				|| (sampling.rot_angles.size() != sampling.tilt_angles.size()) )
			REPORT_ERROR("ERROR: sampling.rot, tilt, psi_angles.size() are invalid!");
		if ( (pointer_dir_nonzeroprior.size() < 1) || (pointer_dir_nonzeroprior.size() > sampling.rot_angles.size())
				|| (pointer_psi_nonzeroprior.size() < 1) || (pointer_psi_nonzeroprior.size() > sampling.psi_angles.size()) )
			REPORT_ERROR("ERROR: pointer_dir_nonzeroprior.size() and/or pointer_psi_nonzeroprior.size() are invalid!");

		nr_dir = pointer_dir_nonzeroprior.size() * pointer_psi_nonzeroprior.size();
	}
	else
	{
		aas.clear(); bbs.clear(); ggs.clear();
		if (aa_range > XMIPP_EQUAL_ACCURACY)
		{
			for (val = aa_init + aa_residue - aa_range; val < aa_init + aa_range + XMIPP_EQUAL_ACCURACY; val += ang_search_step)
				aas.push_back(val);
		}
		else
			aas.push_back(aa_init);
		if (bb_range > XMIPP_EQUAL_ACCURACY)
		{
			for (val = bb_init + bb_residue - bb_range; val < bb_init + bb_range + XMIPP_EQUAL_ACCURACY; val += ang_search_step)
				bbs.push_back(val);
		}
		else
			bbs.push_back(bb_init);
		if (gg_range > XMIPP_EQUAL_ACCURACY)
		{
			for (val = gg_init + gg_residue - gg_range; val < gg_init + gg_range + XMIPP_EQUAL_ACCURACY; val += ang_search_step)
				ggs.push_back(val);
		}
		else
			ggs.push_back(gg_init);

		nr_dir = aas.size() * bbs.size() * ggs.size();
	}

	// Translational samplings
	dxs.clear(); dys.clear(); dzs.clear();
	if (dx_range > XMIPP_EQUAL_ACCURACY)
	{
		for (val = dx_residue - dx_range; val < dx_range + XMIPP_EQUAL_ACCURACY; val += trans_search_step)
			dxs.push_back(val);
	}
	else
		dxs.push_back(0.);

	if (dy_range > XMIPP_EQUAL_ACCURACY)
	{
		for (val = dy_residue - dy_range; val < dy_range + XMIPP_EQUAL_ACCURACY; val += trans_search_step)
			dys.push_back(val);
	}
	else
		dys.push_back(0.);

	if (dz_range > XMIPP_EQUAL_ACCURACY)
	{
		for (val = dz_residue - dz_range; val < dz_range + XMIPP_EQUAL_ACCURACY; val += trans_search_step)
			dzs.push_back(val);
	}
	else
		dzs.push_back(0.);

#ifdef DEBUG
	if (verb)
	{
		if (use_healpix)
		{
			std::cout << "  PSI = " << std::flush;
			for (int ii = 0; ii < pointer_psi_nonzeroprior.size(); ii++)
				std::cout << sampling.psi_angles[pointer_psi_nonzeroprior[ii]] << ", " << std::flush;
		}
		else
		{
			std::cout << "  ROT = " << std::flush;
			for (int ii = 0; ii < aas.size(); ii++)
				std::cout << aas[ii] << ", " << std::flush;
			std::cout << std::endl << " TILT = " << std::flush;
			for (int ii = 0; ii < bbs.size(); ii++)
				std::cout << bbs[ii] << ", " << std::flush;
			std::cout << std::endl << "  PSI = " << std::flush;
			for (int ii = 0; ii < ggs.size(); ii++)
				std::cout << ggs[ii] << ", " << std::flush;
		}
		std::cout << std::endl << "   DX = " << std::flush;
		for (int ii = 0; ii < dxs.size(); ii++)
			std::cout << dxs[ii] << ", " << std::flush;
		std::cout << std::endl << "   DY = " << std::flush;
		for (int ii = 0; ii < dys.size(); ii++)
			std::cout << dys[ii] << ", " << std::flush;
		std::cout << std::endl << "   DZ = " << std::flush;
		for (int ii = 0; ii < dzs.size(); ii++)
			std::cout << dzs[ii] << ", " << std::flush;
		std::cout << std::endl << "  NR_TOTAL_DIR = " << nr_dir << ", NR_TOTAL_TRANS <= " << dxs.size() * dys.size() * dzs.size() << std::endl;
	}
#endif

	// Get all sampling points
	op_samplings.clear();
	op_tmp.initZeros(NR_LOCALSYM_PARAMETERS);
	nr_all_samplings = 0;
	// For translations: op_ori = op_int + op_res
	if (dx_range < XMIPP_EQUAL_ACCURACY)
		dx_range = (1e+10);
	if (dy_range < XMIPP_EQUAL_ACCURACY)
		dy_range = (1e+10);
	if (dz_range < XMIPP_EQUAL_ACCURACY)
		dz_range = (1e+10);
	for (int idz = 0; idz < dzs.size(); idz++)
	{
		for (int idy = 0; idy < dys.size(); idy++)
		{
			for (int idx = 0; idx < dxs.size(); idx++)
			{
				dz = dzs[idz]; dy = dys[idy]; dx = dxs[idx];
				r2 = (dz * dz) / (dz_range * dz_range) + (dy * dy) / (dy_range * dy_range) + (dx * dx) / (dx_range * dx_range);
				if ( (r2 - XMIPP_EQUAL_ACCURACY) > 1.)
					continue;

				if (use_healpix)
				{
					for (int idir = 0; idir < pointer_dir_nonzeroprior.size(); idir++)
					{
						aa = sampling.rot_angles[pointer_dir_nonzeroprior[idir]];
						bb = sampling.tilt_angles[pointer_dir_nonzeroprior[idir]];
						for (int ipsi = 0; ipsi < pointer_psi_nonzeroprior.size(); ipsi++)
						{
							gg = sampling.psi_angles[pointer_psi_nonzeroprior[ipsi]];

							// Re-calculate op_old so that they follow the conventions in RELION!
							standardiseEulerAngles(aa, bb, gg, aa, bb, gg);

							Localsym_composeOperator(op_tmp, aa, bb, gg, dx + dx_init, dy + dy_init, dz + dz_init, (1e10));

							op_samplings.push_back(op_tmp);
							nr_all_samplings++;
						}
					}
				}
				else
				{
					for (int iaa = 0; iaa < aas.size(); iaa++)
					{
						for (int ibb = 0; ibb < bbs.size(); ibb++)
						{
							for (int igg = 0; igg < ggs.size(); igg++)
							{
								// Re-calculate op_old so that they follow the conventions in RELION!
								standardiseEulerAngles(aas[iaa], bbs[ibb], ggs[igg], aa, bb, gg);

								Localsym_composeOperator(op_tmp, aa, bb, gg, dx + dx_init, dy + dy_init, dz + dz_init, (1e10));

								op_samplings.push_back(op_tmp);
								nr_all_samplings++;
							}
						}
					}
				}
			}
		}
	}

	if (verb)
	{
#ifdef __unix__
		std::cout << " + Total sampling points = " << "\e[1m" << op_samplings.size() << "\e[0m" << std::flush;
#else
		std::cout << " + Total sampling points = " << op_samplings.size() << std::flush;
#endif
		std::cout << ", calculating cross-correlation (CC) values ..." << std::endl;
	}

	if (op_samplings.size() < 1)
		REPORT_ERROR("ERROR: No sampling points!");
}

void calculateOperatorCC(
		const MultidimArray<RFLOAT>& src,
		const MultidimArray<RFLOAT>& dest,
		const MultidimArray<RFLOAT>& mask,
		std::vector<Matrix1D<RFLOAT> >& op_samplings,
		bool do_sort,
		bool verb)
{
	RFLOAT val = 0., mask_val = 0., mask_val_sum = 0., mask_val_ctr = 0., cc = 0.;
	int barstep = 0, updatebar = 0, totalbar = 0;
	Matrix2D<RFLOAT> op_mat;
	MultidimArray<RFLOAT> vol;

	if (op_samplings.size() < 1)
		REPORT_ERROR("ERROR: No sampling points!");

	if ( (!isMultidimArray3DCubic(src)) || (!isMultidimArray3DCubic(dest)) || (!isMultidimArray3DCubic(mask)) )
		REPORT_ERROR("ERROR: MultidimArray src, dest, mask should all be 3D cubic!");
	if ( (!src.sameShape(dest)) || (!src.sameShape(mask)) )
		REPORT_ERROR("ERROR: MultidimArray src, dest, mask should have the same sizes!");

	// Check the mask, calculate the sum of mask values
	sum3DCubicMask(mask, mask_val_sum, mask_val_ctr);
	if (mask_val_sum < 1.)
		std::cout << " + WARNING: sum of mask values is smaller than 1! Please check whether it is a correct mask!" << std::endl;

	// Calculate all CCs
	if (verb)
	{
		//std::cout << " + Calculate CCs for all sampling points ..." << std::endl;
		init_progress_bar(op_samplings.size());
		barstep = op_samplings.size() / 100;
		updatebar = totalbar = 0;
	}
	for (int iop = 0; iop < op_samplings.size(); iop++)
	{
		Localsym_operator2matrix(op_samplings[iop], op_mat, LOCALSYM_OP_DO_INVERT);
		applyGeometry(dest, vol, op_mat, IS_NOT_INV, DONT_WRAP);

		cc = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(vol)
		{
			mask_val = DIRECT_A3D_ELEM(mask, k, i, j);
			if (mask_val < XMIPP_EQUAL_ACCURACY)
				continue;

			val = DIRECT_A3D_ELEM(vol, k, i, j) - DIRECT_A3D_ELEM(src, k, i, j);
			//cc += val * val;
			cc += mask_val * val * val; // weighted by mask value ?
		}
		VEC_ELEM(op_samplings[iop], CC_POS) = sqrt(cc / mask_val_sum);

		if (verb)
		{
			if (updatebar > barstep)
			{
				updatebar = 0;
				progress_bar(totalbar);
			}
			updatebar++;
			totalbar++;
		}
	}
	if (verb)
		progress_bar(op_samplings.size());

	// Sort cc, in descending order
	if (do_sort)
		std::stable_sort(op_samplings.begin(), op_samplings.end(), compareOperatorsByCC);
}

void separateMasksBFS(
		const FileName& fn_in,
		const int K,
		RFLOAT val_thres)
{
	MetaDataTable MD;
	MultidimArray<int> vol_rec;
	Image<RFLOAT> img, img_out;
	FileName fn_out;
	RFLOAT x_angpix = 0., y_angpix = 0., z_angpix = 0., float_val = 0.;
	long int pos_val_ctr = 0, xx = 0, yy = 0, zz = 0;
	int id = 0, int_val = 0;
	std::queue<Matrix1D<int> > q;
	Matrix1D<int> vec1;
	const int K_max = 999;

	// Check K
	if ( (K < 2) || (K > K_max) )
		REPORT_ERROR("ERROR: number of sub-masks should be at 2~999 !");
	if (K > 20)
		std::cerr << " WARNING: K = " << K << " seems too large!" << std::endl;
#ifdef DEBUG
	std::cout << " K = " << K << std::endl;
#endif

    // Read the header of input map
	img.read(fn_in);
	//img().setXmippOrigin();
	if ((NSIZE(img()) != 1) || (ZSIZE(img()) <= 10) || (YSIZE(img()) <= 10) || (XSIZE(img()) <= 10))
		REPORT_ERROR("ERROR: Image file " + fn_in + " is an invalid 3D map! (< 10 X 10 X 10 pixels)");
	if ( (XSIZE(img()) != YSIZE(img())) || (XSIZE(img()) != ZSIZE(img())) )
		REPORT_ERROR("ERROR: Image file " + fn_in + " is not a 3D cubic map!");
	img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, x_angpix);
	img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Y, y_angpix);
	img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Z, z_angpix);

	// Initialise vol_rec
	vol_rec.initZeros(img());
	//vol_rec.setXmippOrigin();

	// Count voxels with positive values
	pos_val_ctr = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
	{
		float_val = DIRECT_A3D_ELEM(img(), k, i, j);
		//if (val < -(XMIPP_EQUAL_ACCURACY))
		//	REPORT_ERROR("ERROR: Image file " + fn_in + " contains negative values!");
		if (float_val > val_thres)
			pos_val_ctr++;
		else
			DIRECT_A3D_ELEM(vol_rec, k, i, j) = -1; // Mark as invalid!
	}
	if (pos_val_ctr <= K)
		REPORT_ERROR("ERROR: Image file " + fn_in + " has nearly no voxels with positive values!");
#ifdef DEBUG
	std::cout << " pos_val_ctr = " << pos_val_ctr << std::endl;
#endif

	id = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(vol_rec)
	{
		int_val = DIRECT_A3D_ELEM(vol_rec, k, i, j);
		if (int_val != 0)
			continue;

		id++;
#ifdef DEBUG
		std::cout << " id= " << id << ", kij= " << k << ", " << i << ", " << j << ", " << std::endl;
#endif
		q.push(vectorR3(int(j), int(i), int(k)));
		DIRECT_A3D_ELEM(vol_rec, k, i, j) = id;
		while(!q.empty())
		{
			vec1 = q.front();
			q.pop();
			DIRECT_A3D_ELEM(vol_rec, ZZ(vec1), YY(vec1), XX(vec1)) = id;
			for (int dz = -1; dz <= 1; dz++)
			{
				for (int dy = -1; dy <= 1; dy++)
				{
					for (int dx = -1; dx <= 1; dx++)
					{
						if ( (dx * dy * dz) != 0)
							continue;
						zz = ZZ(vec1) + dz; yy = YY(vec1) + dy; xx = XX(vec1) + dx;
						if ( (zz < 0) || (zz >= ZSIZE(vol_rec)) || (yy < 0) || (yy >= YSIZE(vol_rec)) || (xx < 0) || (xx >= XSIZE(vol_rec)) )
							continue;

						if (DIRECT_A3D_ELEM(vol_rec, zz, yy, xx) == 0)
						{
							q.push(vectorR3(int(xx), int(yy), int(zz)));
							DIRECT_A3D_ELEM(vol_rec, zz, yy, xx) = id;
						}
					}
				}
			}
		}
	}

	std::cout << " " << id << " region(s) detected on map " << fn_in << "." << std::endl;
	if (K != id)
	{
		std::cout << " But the number of regions specified is " << K << "! Please check your input map. Exit now with no file output..." << std::endl;
#ifndef DEBUG
		return;
#endif
	}

	// Write output maps and STAR file
	MD.clear();
	MD.addLabel(EMDL_MASK_NAME);
	for (int icen = 0; (icen < K) && (icen < id); icen++)
	{
		fn_out = fn_in.withoutExtension() + "_sub" + integerToString(icen + 1, 3, '0') + ".mrc";
		img_out().initZeros(img());
		//img_out().setXmippOrigin();

		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(vol_rec)
		{
			if (DIRECT_A3D_ELEM(vol_rec, k, i, j) == (icen + 1) )
				DIRECT_A3D_ELEM(img_out(), k, i, j) = 1.;
		}

		img_out.setStatisticsInHeader();
		img_out.setSamplingRateInHeader(x_angpix, y_angpix, z_angpix);
		img_out.write(fn_out);

		MD.addObject();
		MD.setValue(EMDL_MASK_NAME, fn_out);
	}
	fn_out = fn_in.withoutExtension() + "_masklist.star";
	MD.write(fn_out);
}

/*
void separateMasksKMeans(
		const FileName& fn_in,
		const int K,
		int random_seed)
{
	Image<RFLOAT> img, img_out;
	std::vector<Matrix1D<RFLOAT> > ocen, ncen;
	std::vector<RFLOAT> wcen;
	std::vector<int> vec_rec;
	Matrix1D<int> vec;
	Matrix2D<RFLOAT> mat;
	RFLOAT a = 0., b = 0., g = 0., x = 0., y = 0., z = 0., val = 0., dist2 = 0, dist2_min = 0.;
	RFLOAT x_angpix = 0., y_angpix = 0., z_angpix = 0.;
	int best_cen = -1, pos_val_ctr = 0;
	long int cen_ptr = 0;
	FileName fn_out;
	MultidimArray<int> vol_rec;
	const int K_max = 999;
	int vec_len_max = 1024000, q = 0;

	// Check K
	if ( (K < 2) || (K > K_max) )
		REPORT_ERROR("ERROR: number of sub-masks should be at 2~999 !");
	if (K > 20)
		std::cerr << " WARNING: K = " << K << " seems too large!" << std::endl;
#ifdef DEBUG
	std::cout << " K = " << K << std::endl;
#endif

	// Initialise arrays for centroids
	ocen.clear(); ncen.clear(); wcen.clear();
	for (int ii = 0; ii < K; ii++)
	{
		ocen.push_back(vectorR3(0., 0., 0.));
		ncen.push_back(vectorR3(0., 0., 0.));
		wcen.push_back(0.);
	}

	// Initialise random number generator
	if (random_seed < 0)
		random_seed = time(NULL);
    init_random_generator(random_seed);

    // Read the header of input map
	img.read(fn_in);
	img().setXmippOrigin();
	if ((NSIZE(img()) != 1) || (ZSIZE(img()) <= 10) || (YSIZE(img()) <= 10) || (XSIZE(img()) <= 10))
		REPORT_ERROR("ERROR: Image file " + fn_in + " is an invalid 3D map! (< 10 X 10 X 10 pixels)");
	if ( (XSIZE(img()) != YSIZE(img())) || (XSIZE(img()) != ZSIZE(img())) )
		REPORT_ERROR("ERROR: Image file " + fn_in + " is not a 3D cubic map!");
	img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, x_angpix);
	img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Y, y_angpix);
	img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_Z, z_angpix);

	// Count voxels with positive values
	pos_val_ctr = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img())
	{
		val = DIRECT_A3D_ELEM(img(), k, i, j);
		//if (val < -(XMIPP_EQUAL_ACCURACY))
		//	REPORT_ERROR("ERROR: Image file " + fn_in + " contains negative values!");
		if (val > XMIPP_EQUAL_ACCURACY)
			pos_val_ctr++;
	}
	if (pos_val_ctr <= K)
		REPORT_ERROR("ERROR: Image file " + fn_in + " has nearly no voxels with positive values!");
#ifdef DEBUG
	std::cout << " pos_val_ctr = " << pos_val_ctr << std::endl;
#endif

	// ????
	vol_rec.initZeros(img());
	vol_rec.setXmippOrigin();

	// Randomly select K centroids
	vec_rec.clear();
	vec_len_max = (vec_len_max >= pos_val_ctr) ? (pos_val_ctr) : (vec_len_max);
	q = pos_val_ctr / vec_len_max;
	q = (q <= 1) ? (1) : (q);
	vec.initZeros(vec_len_max + 1);
	for (int ii = 1; ii <= K; ii++) // Knuth shuffle
	{
		cen_ptr = (long int)(rnd_unif(RFLOAT(ii), RFLOAT(vec_len_max)));
		cen_ptr = (cen_ptr < ii) ? (ii) : (cen_ptr);
		cen_ptr = (cen_ptr > vec_len_max) ? (vec_len_max) : (cen_ptr);
		if (VEC_ELEM(vec, cen_ptr) != 0)
			vec_rec.push_back(q * VEC_ELEM(vec, cen_ptr));
		else
			vec_rec.push_back(q * cen_ptr);
		VEC_ELEM(vec, cen_ptr) = ii;
	}
#ifdef DEBUG
	std::cout << " " << vec_rec.size() << " voxel IDs in total: " << std::flush;
	for (int ii = 0; ii < vec_rec.size(); ii++)
		std::cout << vec_rec[ii] << ", " << std::flush;
	std::cout << std::endl;
#endif
	best_cen = pos_val_ctr = 0;
	FOR_ALL_ELEMENTS_IN_ARRAY3D(img())
	{
		if (best_cen >= K)
			break;
		if (A3D_ELEM(img(), k, i, j) > XMIPP_EQUAL_ACCURACY)
		{
			pos_val_ctr++;
			if (vec_rec[best_cen] == pos_val_ctr)
			{
				ocen[best_cen] = vectorR3(RFLOAT(k), RFLOAT(i), RFLOAT(j));
				best_cen++;
			}
		}
	}
	for (int ii = 0; ii < K; ii++)
	{
#ifdef DEBUG
		std::cout << " Centroid #" << ii + 1 << " : XYZ= " << XX(ocen[ii]) << ", " << YY(ocen[ii]) << ", " << ZZ(ocen[ii]) << std::endl;
#endif
	}

	//a = rnd_unif(-179., 179.);
	//b = rnd_unif(1., 179.);
	//g = rnd_unif(-179., 179.);
	//Euler_angles2matrix(a, b, g, mat);
	//for (int ii = 0; ii < K; ii++)
	//{
	//	z = RFLOAT(ii) * RFLOAT(ZSIZE(img())) / RFLOAT(K) + RFLOAT(STARTINGZ(img()));
	//	y = rnd_unif(0., YSIZE(img())) + RFLOAT(STARTINGY(img()));
	//	x = rnd_unif(0., XSIZE(img())) + RFLOAT(STARTINGX(img()));
	//	z /= sqrt(3.); y /= sqrt(3.); x /= sqrt(3.);

	//	ocen[ii] = mat * vectorR3(x, y, z);
//#ifdef DEBUG
	//	std::cout << " Centroid #" << ii + 1 << " : XYZ= " << XX(ocen[ii]) << ", " << YY(ocen[ii]) << ", " << ZZ(ocen[ii]) << std::endl;
//#endif
	//}

	// K-means
	for (int iter = 1; iter <= 100; iter++)
	{
#ifdef DEBUG
		std::cout << std::endl;
#endif
		FOR_ALL_ELEMENTS_IN_ARRAY3D(img())
		{
			// For voxels with positive values
			val = A3D_ELEM(img(), k, i, j);
			if (val < XMIPP_EQUAL_ACCURACY)
				continue;

			// Find the smallest distance to one of the centroids
			dist2_min = 1e+30;
			best_cen = -1;
			for (int icen = 0; icen < K; icen++)
			{
				z = ZZ(ocen[icen]); y = YY(ocen[icen]); x = XX(ocen[icen]);
				dist2 = (RFLOAT(k) - z) * (RFLOAT(k) - z) + (RFLOAT(i) - y) * (RFLOAT(i) - y) + (RFLOAT(j) - x) * (RFLOAT(j) - x);
				if (dist2 < dist2_min)
				{
					dist2_min = dist2;
					best_cen = icen;
				}
			}
			if (best_cen < 0)
				REPORT_ERROR("ERROR: best_cen < 0 !");

			ZZ(ncen[best_cen]) += k * val;
			YY(ncen[best_cen]) += i * val;
			XX(ncen[best_cen]) += j * val;
			wcen[best_cen] += val;

			A3D_ELEM(vol_rec, k, i, j) = best_cen + 1;
		}

		// Update centroids
		for (int ii = 0; ii < K; ii++)
		{
			if (wcen[ii] < XMIPP_EQUAL_ACCURACY)
				REPORT_ERROR("ERROR: wcen[ii] <= 0 !");
			ocen[ii] = ncen[ii] / wcen[ii];

			ncen[ii] = vectorR3(0., 0., 0.);
			wcen[ii] = 0.;
#ifdef DEBUG
			std::cout << " Centroid #" << ii + 1 << " : XYZ= " << XX(ocen[ii]) << ", " << YY(ocen[ii]) << ", " << ZZ(ocen[ii]) << std::endl;
#endif
		}
	}

	// Write output maps
	for (int icen = 0; icen < K; icen++)
	{
		fn_out = fn_in.withoutExtension() + "_sub" + integerToString(icen + 1, 3, '0') + ".mrc";
		img_out().initZeros(img());
		img_out().setXmippOrigin();

		FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_rec)
		{
			if (A3D_ELEM(vol_rec, k, i, j) == (icen + 1) )
				A3D_ELEM(img_out(), k, i, j) = A3D_ELEM(img(), k, i, j);
		}

		img_out.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, x_angpix);
		img_out.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, y_angpix);
		img_out.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, z_angpix);
		img_out.write(fn_out);
	}
}
*/

void local_symmetry_parameters::initBoolOptions()
{
	show_usage_for_an_option = false;

	do_apply_local_symmetry = false;
	do_duplicate_local_symmetry = false;
	do_local_search_local_symmetry_ops = false;
	do_txt2rln = false;
	do_transform = false;
	do_debug = false;
}

void local_symmetry_parameters::clear()
{
	parser.clear();
	initBoolOptions();
}

void local_symmetry_parameters::displayEmptyLine()
{
	std::cout << "=========================================================================" << std::endl;
}

void local_symmetry_parameters::usage()
{
	parser.writeUsage(std::cerr);
}

void local_symmetry_parameters::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int init_section = parser.addSection("Show usage");
	show_usage_for_an_option = parser.checkOption("--function_help", "Show usage for the selected function (JUN 30, 2017)");

	int options_section = parser.addSection("Options");
	do_apply_local_symmetry = parser.checkOption("--apply", "Apply local symmetry to a 3D cryo-EM density map");
	do_duplicate_local_symmetry = parser.checkOption("--duplicate", "Duplicate subunits/masks according to local symmetry operators");
	do_local_search_local_symmetry_ops = parser.checkOption("--search", "Local searches of local symmetry operators");
	do_transform = parser.checkOption("--transform", "Transform a map according to three Euler angles and XYZ translations");
	do_txt2rln = parser.checkOption("--txt2rln", "Convert operators from DM to RELION STAR format");
	do_debug = parser.checkOption("--debug", "(DEBUG ONLY)");

	int params_section = parser.addSection("Parameters (alphabetically ordered)");
	angpix_image = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms) of input image", "1."));
	ang_range = textToFloat(parser.getOption("--ang_range", "Angular search range of operators (in degrees), overwrite rot-tilt-psi ranges if set to positive", "0."));
	ang_rot_range = textToFloat(parser.getOption("--ang_rot_range", "Angular (rot) search range of operators (in degrees)", "0."));
	ang_tilt_range = textToFloat(parser.getOption("--ang_tilt_range", "Angular (tilt) search range of operators (in degrees)", "0."));
	ang_psi_range = textToFloat(parser.getOption("--ang_psi_range", "Angular (psi) search range of operators (in degrees)", "0."));
	ang_step = textToFloat(parser.getOption("--ang_step", "Angular search step of operators (in degrees)", "1."));
	binning_factor = textToFloat(parser.getOption("--bin", "Binning factor (<= 1 means no binning)", "-1."));
	ini_threshold = textToFloat(parser.getOption("--ini_threshold", "Initial threshold for binarization", "0.01"));
	fn_unsym = parser.getOption("--i_map", "Input 3D unsymmetrised map", "");
	fn_info_in = parser.getOption("--i_mask_info", "Input file with mask filenames and rotational / translational operators (for local searches)", "maskinfo.txt");
	fn_op_mask_info_in = parser.getOption("--i_op_mask_info", "Input file with mask filenames for all operators (for global searches)", "None");
	nr_masks = textToInteger(parser.getOption("--n", "Create this number of masks according to the input density map", "2"));
	offset_range = textToFloat(parser.getOption("--offset_range", "Translational search range of operators (in Angstroms), overwrite x-y-z ranges if set to positive", "0."));
	offset_x_range = textToFloat(parser.getOption("--offset_x_range", "Translational (x) search range of operators (in Angstroms)", "0."));
	offset_y_range = textToFloat(parser.getOption("--offset_y_range", "Translational (y) search range of operators (in Angstroms)", "0."));
	offset_z_range = textToFloat(parser.getOption("--offset_z_range", "Translational (z) search range of operators (in Angstroms)", "0."));
	offset_step = textToFloat(parser.getOption("--offset_step", "Translational search step of operators (in Angstroms)", "1."));
	fn_sym = parser.getOption("--o_map", "Output 3D symmetrised map", "");
	fn_info_out = parser.getOption("--o_mask_info", "Output file with mask filenames and rotational / translational operators", "maskinfo_refined.txt");
	psi = textToFloat(parser.getOption("--psi", "Third Euler angle (psi, in degrees)", "0."));
	rot = textToFloat(parser.getOption("--rot", "First Euler angle (rot, in degrees)", "0."));
	sphere_percentage = textToFloat(parser.getOption("--sphere_percentage", "Diameter of spherical mask divided by the box size (< 0.99)", "-1."));
	tilt = textToFloat(parser.getOption("--tilt", "Second Euler angle (tilt, in degrees)", "0."));
	xoff = textToFloat(parser.getOption("--xoff", "X-offset (in Angstroms)", "0."));
	yoff = textToFloat(parser.getOption("--yoff", "Y-offset (in Angstroms)", "0."));
	zoff = textToFloat(parser.getOption("--zoff", "Z-offset (in Angstroms)", "0."));
	verb = parser.checkOption("--verb", "Verbose output?");

	int expert_section = parser.addSection("Parameters (expert options - alphabetically ordered)");
	fn_mask = parser.getOption("--i_mask", "(DEBUG) Input mask", "mask.mrc");
	fn_info_in_parsed_ext = parser.getOption("--i_mask_info_parsed_ext", "Extension of parsed input file with mask filenames and rotational / translational operators", "parsed");
	use_healpix_sampling = parser.checkOption("--use_healpix", "Use Healpix for angular samplings?");
	width_edge_pix = textToFloat(parser.getOption("--width", "Width of cosine soft edge (in pixels)", "5."));

   	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
};

void local_symmetry_parameters::run()
{
	bool do_sort = true, do_verb = true;
	FileName fn_parsed, fn_tmp;
	std::vector<FileName> fn_mask_list;
	std::vector<std::vector<Matrix1D<RFLOAT> > > op_list;

	fn_mask_list.clear();
	op_list.clear();

	// Check options
	int valid_options = 0;
	valid_options += (do_apply_local_symmetry) ? (1) : (0);
	valid_options += (do_duplicate_local_symmetry) ? (1) : (0);
	valid_options += (do_local_search_local_symmetry_ops) ? (1) : (0);
	valid_options += (do_txt2rln) ? (1) : (0);
	valid_options += (do_transform) ? (1) : (0);
	valid_options += (do_debug) ? (1) : (0);

	if (valid_options <= 0)
		REPORT_ERROR("Please specify one option!");
	if (valid_options > 1)
		REPORT_ERROR("Only one option can be specified at one time! valid_options = " + integerToString(valid_options));

	if (do_apply_local_symmetry)
	{
		if (show_usage_for_an_option)
		{
			displayEmptyLine();
			std::cout << " Apply local symmetry to a 3D cryo-EM density map" << std::endl;
			std::cout << "  USAGE: --apply --angpix 1.34 --i_map unsym.mrc --i_mask_info maskinfo.star --o_map sym.mrc (--sphere_percentage 0.9)" << std::endl;
			displayEmptyLine();
			return;
		}

		Image<RFLOAT> unsym_map, sym_map;

		if (angpix_image < 0.001)
			REPORT_ERROR("Invalid pixel size!");
		if (sphere_percentage > 0.991)
			REPORT_ERROR("Diameter of spherical mask divided by the box size should be smaller than 0.99!");

		// Parse mask info file
		if (fn_info_in.getExtension() == "star")
		{
			readRelionFormatMasksAndOperators(fn_info_in, fn_mask_list, op_list, angpix_image, do_verb);
		}
		else
		{
			fn_parsed = fn_info_in + std::string(".") + fn_info_in_parsed_ext;
			parseDMFormatMasksAndOperators(fn_info_in, fn_parsed);
			readDMFormatMasksAndOperators(fn_parsed, fn_mask_list, op_list, angpix_image, do_verb);
		}

		unsym_map.clear();
		unsym_map.read(fn_unsym);
		//sym_map.clear();

		int box_size = ((XSIZE(unsym_map())) < (YSIZE(unsym_map()))) ? (XSIZE(unsym_map())) : (YSIZE(unsym_map()));
		box_size = (box_size < (ZSIZE(unsym_map()))) ? box_size : (ZSIZE(unsym_map()));

		applyLocalSymmetry(sym_map(), unsym_map(), fn_mask_list, op_list, (RFLOAT(box_size) * sphere_percentage) / 2., width_edge_pix);
		sym_map().setXmippOrigin();

		sym_map.setSamplingRateInHeader(angpix_image, angpix_image, angpix_image);
		sym_map.setStatisticsInHeader();

		sym_map.write(fn_sym);
	}
	else if (do_duplicate_local_symmetry)
	{
		if (show_usage_for_an_option)
		{
			displayEmptyLine();
			std::cout << " Duplicate subunits/masks according to local symmetry operators" << std::endl;
			std::cout << "  USAGE: --duplicate (--i_map unsym.mrc) --angpix 1.34 --i_mask_info maskinfo.txt --o_map duplicated.mrc" << std::endl;
			std::cout << "  Leave '--i_map' empty if you want to duplicate masks only." << std::endl;
			displayEmptyLine();
			return;
		}

		Image<RFLOAT> map_in, map_out;
		bool duplicate_masks_only = true;

		if (angpix_image < 0.001)
			REPORT_ERROR("Invalid pixel size!");

		// Parse mask info file
		if (fn_info_in.getExtension() == "star")
		{
			readRelionFormatMasksAndOperators(fn_info_in, fn_mask_list, op_list, angpix_image, do_verb);
		}
		else
		{
			fn_parsed = fn_info_in + std::string(".") + fn_info_in_parsed_ext;
			parseDMFormatMasksAndOperators(fn_info_in, fn_parsed);
			readDMFormatMasksAndOperators(fn_parsed, fn_mask_list, op_list, angpix_image, do_verb);
		}

		map_in.clear();
		//map_out.clear();
		if (exists(fn_unsym))
		{
			duplicate_masks_only = false;
			map_in.read(fn_unsym);
		}
		else
			duplicate_masks_only = true;
		duplicateLocalSymmetry(map_out(), map_in(), fn_mask_list, op_list, duplicate_masks_only);
		map_out().setXmippOrigin();

		map_out.setSamplingRateInHeader(angpix_image, angpix_image, angpix_image);
		map_out.setStatisticsInHeader();

		map_out.write(fn_sym);
	}
	else if (do_local_search_local_symmetry_ops)
	{
		if (show_usage_for_an_option)
		{
			displayEmptyLine();
			std::cout << " Searches of local symmetry operators" << std::endl;
			std::cout << "    MPI: mpirun -n 23 relion_localsym_mpi ..." << std::endl;
			std::cout << "  USAGE FOR GLOBAL SEARCHES:" << std::endl;
			std::cout << "         --search --i_map unsym.mrc --i_op_mask_info mask_list.star --o_mask_info maskinfo_iter000.star --angpix 1.34 (--bin 2)" << std::endl;
			std::cout << "         --ang_step 5 (--offset_range 2 --offset_step 1)" << std::endl;
			std::cout << "  USAGE FOR LOCAL SEARCHES:" << std::endl;
			std::cout << "         --search --i_map unsym.mrc --i_mask_info maskinfo_iter001.star --o_mask_info maskinfo_iter002.star --angpix 1.34 (--bin 2)" << std::endl;
			std::cout << "         --ang_range 2 (--ang_rot_range 2 --ang_tilt_range 2 --ang_psi_range 2) --ang_step 0.5" << std::endl;
			std::cout << "         --offset_range 2 (--offset_x_range 2 --offset_y_range 2 --offset_z_range 2) --offset_step 1" << std::endl;
			std::cout << "  Ranges/steps of angular and translational searches are in degrees and Angstroms respectively." << std::endl;
			displayEmptyLine();
			return;
		}

		// Adapted from local_symmetry_mpi.cpp

		long int newdim = 0, cropdim = 0, z0 = 0, y0 = 0, x0 = 0, zf = 0, yf = 0, xf = 0;
		RFLOAT aa = 0, bb = 0, gg = 0., dx = 0., dy = 0., dz = 0., cc = 0., tmp_binning_factor = 1.;
		RFLOAT mask_sum = 0., mask_ctr = 0., mask2_sum = 0., mask2_ctr = 0.;

		Image<RFLOAT> map, mask, mask2;
		Matrix1D<RFLOAT> op_search_ranges, op, com0_int, com1_int, com1_float, com1_diff, vecR3;
		//std::vector<FileName> fn_mask_list;
		//std::vector<std::vector<Matrix1D<RFLOAT> > > op_list;
		std::vector<std::vector<FileName> > op_mask_list;
		std::vector<Matrix1D<RFLOAT> > op_samplings;
		MultidimArray<RFLOAT> src_cropped, dest_cropped, mask_cropped;
		Matrix2D<RFLOAT> mat1;
		FileName fn_searched_op_samplings;
		//FileName fn_parsed, fn_tmp;

		map.clear(); mask.clear(); mask2.clear();
		op_search_ranges.clear(); op.clear(); com0_int.clear(); com1_int.clear(); com1_float.clear(); com1_diff.clear(); vecR3.clear();
		fn_mask_list.clear();
		op_list.clear();
		op_mask_list.clear();
		op_samplings.clear();
		src_cropped.clear(); dest_cropped.clear(); mask_cropped.clear();
		mat1.clear();
		fn_parsed.clear(); fn_tmp.clear(); fn_searched_op_samplings.clear();

		displayEmptyLine();

		// Master gets search ranges (in degrees and pixels), sets offset_step (in pixels).
		if (angpix_image < 0.001)
			REPORT_ERROR("Invalid pixel size!");
		if (fn_op_mask_info_in != "None")
		{
			if ( (ang_range < (XMIPP_EQUAL_ACCURACY) )
					&& (ang_rot_range < (XMIPP_EQUAL_ACCURACY) )
					&& (ang_tilt_range < (XMIPP_EQUAL_ACCURACY) )
					&& (ang_psi_range < (XMIPP_EQUAL_ACCURACY) ) )
			{
				ang_range = 180.;
				std::cout << " Initial searches: reset searching ranges of all 3 Euler angles to +/-180 degrees." << std::endl;
			}
			else
			{
				if (ang_range > (XMIPP_EQUAL_ACCURACY) )
					std::cout << " User-defined initial searches: searching ranges of all 3 Euler angles are set to +/-" << ang_range << " degree(s)." << std::endl;
				else
					std::cout << " User-defined initial searches: (rot, tilt, psi) ranges are +/- (" << ang_rot_range << ", " << ang_tilt_range << ", " << ang_psi_range << ") degree(s)." << std::endl;
			}
		}
		Localsym_composeOperator(
				op_search_ranges,
				(ang_range > (XMIPP_EQUAL_ACCURACY)) ? (ang_range) : (ang_rot_range),
				(ang_range > (XMIPP_EQUAL_ACCURACY)) ? (ang_range) : (ang_tilt_range),
				(ang_range > (XMIPP_EQUAL_ACCURACY)) ? (ang_range) : (ang_psi_range),
				(offset_range > (XMIPP_EQUAL_ACCURACY)) ? (offset_range) : (offset_x_range),
				(offset_range > (XMIPP_EQUAL_ACCURACY)) ? (offset_range) : (offset_y_range),
				(offset_range > (XMIPP_EQUAL_ACCURACY)) ? (offset_range) : (offset_z_range) );
		Localsym_scaleTranslations(op_search_ranges, 1. / angpix_image);
		offset_step /= angpix_image;

		// Master parses and reads mask info file
		// Local searches
		if (fn_op_mask_info_in == "None")
		{
			if (fn_info_in.getExtension() == "star")
			{
				readRelionFormatMasksAndOperators(fn_info_in, fn_mask_list, op_list, angpix_image, true);
			}
			else
			{
				fn_parsed = fn_info_in + std::string(".") + fn_info_in_parsed_ext;
				parseDMFormatMasksAndOperators(fn_info_in, fn_parsed);
				readDMFormatMasksAndOperators(fn_parsed, fn_mask_list, op_list, angpix_image, true);
			}
		}
		else
		{
			// Global searches
			std::cout << " Global searches: option --i_mask_info " << fn_info_in << " is ignored." << std::endl;
			readRelionFormatMasksWithoutOperators(fn_op_mask_info_in, fn_mask_list, op_list, op_mask_list, (ang_range > 179.99), true);
		}

		// Master reads input map
		std::cout << std::endl << " Pixel size = " << angpix_image << " Angstrom(s)" << std::endl;
		std::cout << " Read input map " << fn_unsym << " ..." << std::endl;
		map.read(fn_unsym);
		map().setXmippOrigin();
		if (!isMultidimArray3DCubic(map()))
			REPORT_ERROR("ERROR: Input map " + fn_unsym + " is not 3D cube!");

		// All nodes loop over all masks
		for (int imask = 0; imask < fn_mask_list.size(); imask++)
		{
			displayEmptyLine();

			// Master reads and checks the mask
			std::cout << " Read mask #" << imask + 1 << ": " << fn_mask_list[imask] << " ..." << std::endl;
			mask.read(fn_mask_list[imask]);
			mask().setXmippOrigin();
			if (!isMultidimArray3DCubic(mask()))
				REPORT_ERROR("ERROR: Input mask " + fn_mask_list[imask] + " is not 3D cube!");
			if (!map().sameShape(mask()))
				REPORT_ERROR("ERROR: Input map " + fn_unsym + " and mask " + fn_mask_list[imask] + " should have the same size!");
			sum3DCubicMask(mask(), mask_sum, mask_ctr);

			// Get com0 of this mask. Assume that com0 has all integer values!
			getMinCropSize(mask(), com0_int, cropdim, offset_range / angpix_image);
			if (cropdim < 2)
				REPORT_ERROR("ERROR: Mask " + fn_mask_list[imask] + " is too small!");
			XX(com0_int) = round(XX(com0_int));
			YY(com0_int) = round(YY(com0_int));
			ZZ(com0_int) = round(ZZ(com0_int));
			std::cout << " Mask #" << imask + 1 << " : center of mass XYZ = (" << XX(com0_int) << ", " << YY(com0_int) << ", " << ZZ(com0_int) << ") pixel(s)."<< std::endl;

			// Crop the mask and the corresponding region of the map
			z0 = ROUND(ZZ(com0_int)) + FIRST_XMIPP_INDEX(cropdim);
			zf = ROUND(ZZ(com0_int)) + LAST_XMIPP_INDEX(cropdim);
			y0 = ROUND(YY(com0_int)) + FIRST_XMIPP_INDEX(cropdim);
			yf = ROUND(YY(com0_int)) + LAST_XMIPP_INDEX(cropdim);
			x0 = ROUND(XX(com0_int)) + FIRST_XMIPP_INDEX(cropdim);
			xf = ROUND(XX(com0_int)) + LAST_XMIPP_INDEX(cropdim);

			std::cout << " Mask #" << imask + 1 << " : cropped box size = " << cropdim << " pixels." << std::endl;
#ifdef DEBUG
			std::cout << " Window: x0, y0, z0 = " << x0 << ", " << y0 << ", " << z0 << "; xf, yf, zf = " << xf << ", " << yf << ", " << zf << std::endl;
#endif
			mask().window(mask_cropped, z0, y0, x0, zf, yf, xf);
			mask_cropped.setXmippOrigin();
			map().window(src_cropped, z0, y0, x0, zf, yf, xf);
			src_cropped.setXmippOrigin();

			// Rescale the map and the mask (if binning_factor > 1), set 'newdim'.
			tmp_binning_factor = 1.;
			newdim = cropdim;
			if ((binning_factor - 1.) > XMIPP_EQUAL_ACCURACY)
			{
				newdim = (long int)(ceil(RFLOAT(cropdim) / binning_factor));
				if (newdim < 2)
					REPORT_ERROR("ERROR: Binning factor is too large / Mask is too small!");
				if ((newdim + 1) < cropdim) // Need rescaling
				{
					// Dimension should always be even
					if (newdim % 2)
						newdim++;
					resizeMap(mask_cropped, newdim);
					mask_cropped.setXmippOrigin();
					resizeMap(src_cropped, newdim);
					src_cropped.setXmippOrigin();
					tmp_binning_factor = RFLOAT(cropdim) / RFLOAT(newdim);
					std::cout << " + Rescale cropped box size from " << cropdim << " to " << newdim << " pixels. Binning factor = " << tmp_binning_factor << std::endl;

					// Mask values might go out of range after rescaling. Fix it if it happens
					truncateMultidimArray(mask_cropped, 0., 1.);
				}
				else
					newdim = cropdim;
			}
#ifdef DEBUG
			std::cout << " newdim= " << newdim << ", cropdim= " << cropdim << std::endl;
#endif

			// All nodes loop over all operators of this mask
			for (int iop = 0; iop < op_list[imask].size(); iop++)
			{
				std::cout << std::endl;

				// Master gets sampling points
				com1_float.initZeros(3);
				com1_int.initZeros(3);
				com1_diff.initZeros(3);

				Localsym_decomposeOperator(op_list[imask][iop], aa, bb, gg, dx, dy, dz, cc);

				if (fn_op_mask_info_in == "None")
				{
					// Local searches
					// Get com1_float. (floating point numbers)
					// Com1f = R * Com0 + v
					Euler_angles2matrix(aa, bb, gg, mat1);
					com1_float = mat1 * com0_int;
					com1_float += vectorR3(dx, dy, dz);
				}
				else
				{
					// Global searches
					// Master reads and checks the mask
					std::cout << " Read mask #" << imask + 1 << " operator #" << iop + 1 << " : " << op_mask_list[imask][iop] << " ..." << std::endl;
					mask2.read(op_mask_list[imask][iop]);
					mask2().setXmippOrigin();
					if (!isMultidimArray3DCubic(mask2()))
						REPORT_ERROR("ERROR: Input mask " + op_mask_list[imask][iop] + " is not 3D cube!");
					if (!map().sameShape(mask2()))
						REPORT_ERROR("ERROR: Input map " + fn_unsym + " and mask " + op_mask_list[imask][iop] + " should have the same size!");
					sum3DCubicMask(mask2(), mask2_sum, mask2_ctr);
					if (!similar3DCubicMasks(mask_sum, mask_ctr, mask2_sum, mask2_ctr))
						std::cerr << " WARNING: masks " << fn_mask_list[imask] << " and " << op_mask_list[imask][iop] << " seem different! Please check whether they are covering regions from the same set!" << std::endl;

					// Calculate Com1f of this mask
					mask2().centerOfMass(com1_float);
					std::cout << " Mask #" << imask + 1 << " operator #" << iop + 1 << " : center of mass XYZ = (" << XX(com1_float) << ", " << YY(com1_float) << ", " << ZZ(com1_float) << ") pixel(s)."<< std::endl;
				}

				// Get com1_int and com1_diff
				// diff = Com1f - Com1i
				XX(com1_int) = round(XX(com1_float));
				YY(com1_int) = round(YY(com1_float));
				ZZ(com1_int) = round(ZZ(com1_float));
				XX(com1_diff) = XX(com1_float) - XX(com1_int);
				YY(com1_diff) = YY(com1_float) - YY(com1_int);
				ZZ(com1_diff) = ZZ(com1_float) - ZZ(com1_int);

				// Crop this region
				z0 = ROUND(ZZ(com1_int)) + FIRST_XMIPP_INDEX(cropdim);
				zf = ROUND(ZZ(com1_int)) + LAST_XMIPP_INDEX(cropdim);
				y0 = ROUND(YY(com1_int)) + FIRST_XMIPP_INDEX(cropdim);
				yf = ROUND(YY(com1_int)) + LAST_XMIPP_INDEX(cropdim);
				x0 = ROUND(XX(com1_int)) + FIRST_XMIPP_INDEX(cropdim);
				xf = ROUND(XX(com1_int)) + LAST_XMIPP_INDEX(cropdim);
#ifdef DEBUG
				std::cout << " Window: x0, y0, z0 = " << x0 << ", " << y0 << ", " << z0 << "; xf, yf, zf = " << xf << ", " << yf << ", " << zf << std::endl;
#endif
				map().window(dest_cropped, z0, y0, x0, zf, yf, xf);
				dest_cropped.setXmippOrigin();

				// Do the same rescaling
				if (newdim != cropdim)
				{
					resizeMap(dest_cropped, newdim);
					dest_cropped.setXmippOrigin();
				}

				// Master gets sampling points
				// Get sampling points - Rescale translational search ranges and steps
				Localsym_composeOperator(op, aa, bb, gg, XX(com1_diff), YY(com1_diff), ZZ(com1_diff), cc);
				if (newdim != cropdim)
				{
					Localsym_scaleTranslations(op_search_ranges, 1. / tmp_binning_factor);
					offset_step *= 1. / tmp_binning_factor;
					Localsym_scaleTranslations(op, 1. / tmp_binning_factor);
				}
#ifdef __unix__
				std::cout << " + Refining " << "\e[1m" << "Mask #" << imask + 1 << " Operator #" << iop + 1 << "\e[0m" << ": " << std::flush;
#else
				std::cout << " + Refining Mask #" << imask + 1 << " Operator #" << iop + 1 << ": " << std::flush;
#endif
				Localsym_outputOperator(op_list[imask][iop], &std::cout, angpix_image);
				std::cout << std::endl;
				getLocalSearchOperatorSamplings(
						op,
						op_search_ranges,
						op_samplings,
						ang_step,
						offset_step,
						use_healpix_sampling,
						true);
				if (newdim != cropdim)
				{
					Localsym_scaleTranslations(op_search_ranges, tmp_binning_factor);
					offset_step *= tmp_binning_factor;
					Localsym_scaleTranslations(op, tmp_binning_factor);
				}

				// TODO: test this!!!
				if (op_samplings.size() <= 0)
					REPORT_ERROR("ERROR: No sampling points!");

				// Calculate all CCs for the sampling points
				calculateOperatorCC(src_cropped, dest_cropped, mask_cropped, op_samplings, false, do_verb);

				// TODO: For rescaled maps
				if (newdim != cropdim)
				{
					for (int isamp = 0; isamp < op_samplings.size(); isamp++)
						Localsym_scaleTranslations(op_samplings[isamp], tmp_binning_factor);
				}
				// Now translations are all unscaled.

				// TODO: add vectors together!!!
				// Update com1_float
				for (int isamp = 0; isamp < op_samplings.size(); isamp++)
				{
					// Get new_com1
					// newCom1f = Com1f + best_trans_samp - diff
					Localsym_shiftTranslations(op_samplings[isamp], com1_float - com1_diff); // equivalently, com1_int

					// Update v = newCom1f + ( - newR * com0)
					Localsym_decomposeOperator(op_samplings[isamp], aa, bb, gg, dx, dy, dz, cc);
					Euler_angles2matrix(aa, bb, gg, mat1);
					vecR3 = vectorR3(dx, dy, dz) - mat1 * com0_int;
					Localsym_composeOperator(op_samplings[isamp], aa, bb, gg, XX(vecR3), YY(vecR3), ZZ(vecR3), cc);
				}

				// Master sorts the results
				std::stable_sort(op_samplings.begin(), op_samplings.end(), compareOperatorsByCC);

				// Master outputs the local searches results
				fn_tmp.compose(fn_info_out.withoutExtension() + "_cc_mask", imask + 1, "tmp", 3);  // "*_cc_mask001.tmp"
				fn_tmp = fn_tmp.withoutExtension(); // "*_cc_mask001"
				fn_searched_op_samplings.compose(fn_tmp + "_op", iop + 1, "star", 3); // "*_cc_mask001_op001.star"
				writeRelionFormatLocalSearchOperatorResults(fn_searched_op_samplings, op_samplings, angpix_image);
				std::cout << " + List of sampling points for this local symmetry operator: " << fn_searched_op_samplings << std::endl;

				// Master updates this operator and do screen output
				op_list[imask][iop] = op_samplings[0];
				std::cout << " + Done! Refined operator: " << std::flush;
				Localsym_outputOperator(op_samplings[0], &std::cout, angpix_image);
				std::cout << std::endl;
			}
		}

		// Master writes out new mask info file
		if (fn_info_out.getExtension() == "star")
			writeRelionFormatMasksAndOperators(fn_info_out, fn_mask_list, op_list, angpix_image);
		else
			writeDMFormatMasksAndOperators(fn_info_out, fn_mask_list, op_list, angpix_image);

		displayEmptyLine();
#ifdef __unix__
		std::cout << " Done! New local symmetry description file: " << "\e[1m" << fn_info_out << "\e[0m" << std::endl;
#else
		std::cout << " Done! New local symmetry description file: " << fn_info_out << std::endl;
#endif
	}
	else if (do_txt2rln)
	{
		if (show_usage_for_an_option)
		{
			displayEmptyLine();
			std::cout << " Convert operators from DM to RELION STAR format" << std::endl;
			std::cout << "  USAGE: --txt2rln --i_mask_info in.txt --o_mask_info out.star" << std::endl;
			displayEmptyLine();
			return;
		}

		if ( (fn_info_in.getExtension() == "star") || (fn_info_out.getExtension() != "star") )
			REPORT_ERROR("ERROR: input and output text files should be in plain-text (not .star) and .star formats respectively!");

		fn_parsed = fn_info_in + std::string(".") + fn_info_in_parsed_ext;
		parseDMFormatMasksAndOperators(fn_info_in, fn_parsed);
		readDMFormatMasksAndOperators(fn_parsed, fn_mask_list, op_list, 1., do_verb);
		writeRelionFormatMasksAndOperators(fn_info_out, fn_mask_list, op_list, 1.);
	}
	else if (do_transform)
	{
		if (show_usage_for_an_option)
		{
			displayEmptyLine();
			std::cout << " Transform a map according to three Euler angles and XYZ translations" << std::endl;
			std::cout << "  USAGE: --transform --angpix 1.34 --i_map in.mrc --o_map out.mrc --rot 5 --tilt 5 --psi 5 --xoff 5 --yoff 5 --zoff 5" << std::endl;
			displayEmptyLine();
			return;
		}

		Image<RFLOAT> img;
		Matrix2D<RFLOAT> op_mat;
		Matrix1D<RFLOAT> op;

		img.read(fn_unsym);
		standardiseEulerAngles(rot, tilt, psi, rot, tilt, psi);
		Localsym_composeOperator(op, rot, tilt, psi, xoff / angpix_image, yoff / angpix_image, zoff / angpix_image);

		std::cout << " Pixel size = " << angpix_image << " Angstrom(s)" << std::endl;
		std::cout << " Transform input map " << fn_unsym << " : " << std::flush;
		Localsym_outputOperator(op, &std::cout, angpix_image);
		std::cout << std::endl;

		Localsym_operator2matrix(op, op_mat, LOCALSYM_OP_DONT_INVERT);
		selfApplyGeometry(img(), op_mat, IS_NOT_INV, DONT_WRAP);
		img.write(fn_sym);

		std::cout << " Done writing " << fn_sym << std::endl;
	}
	else if (do_debug)
	{
		//separateMasksKMeans(fn_unsym, 4);
		separateMasksBFS(fn_unsym, nr_masks, ini_threshold);

		/*
		Image<RFLOAT> img1, img2;
		std::vector<Matrix1D<RFLOAT> > op_samplings;
		RFLOAT aa = 0., bb = 0., gg = 0., dx = 0., dy = 0., dz = 0.;
		Matrix2D<RFLOAT> op_mat1, op_mat2;
		Matrix1D<RFLOAT> op_old, op_search_ranges, op_new, trans_vec1, trans_vec2;

		img1.read(fn_unsym);
		img2.read(fn_mask);

		Localsym_composeOperator(op_old, 36., 130., -110., -6., 4., -5.);
		Localsym_composeOperator(op_search_ranges, 2., 2., 2., 2., 2., 2.);

		localRefineOneOperator(
				img1(),
				img2(),
				op_old,
				op_search_ranges,
				op_samplings,
				0.5,
				1.);

		return;

		aa = 37.6; bb = 129.3; gg = -111.9; dx = -4.87; dy = 5.22; dz = -3.8;
		Localsym_composeOperator(op_old, aa, bb, gg, dx, dy, dz);

		Localsym_operator2matrix(op_old, op_mat1);
		Localsym_operator2matrix(op_old, op_mat2, LOCALSYM_OP_DO_INVERT);

		img1.read(fn_unsym);
		img2 = img1;
		applyGeometry(img1(), img2(), op_mat1, IS_NOT_INV, DONT_WRAP);
		img2.write(fn_sym);

		return;
		*/
	}
	else
	{
		REPORT_ERROR("Please specify an option!");
	}
	if ( (!show_usage_for_an_option) && (!do_debug) )
	{
		writeCommand("relion_localsym.log", "`which relion_localsym`");
	}
}

void local_symmetry_parameters::writeCommand(FileName fn_cmd, std::string str_executable_name)
{
	std::ofstream ofs;
	ofs.open(fn_cmd.c_str(), std::ofstream::out | std::ofstream::app);

	time_t now = time(0);
    char nodename[64] = "undefined";
    gethostname(nodename,sizeof(nodename));
    std::string hostname(nodename);
	ofs << std::endl << " ++++ Executed the following command at host " << hostname << " on " << ctime(&now);
	ofs << "  " << str_executable_name << " " << std::flush;
	parser.writeCommandLine(ofs);
	ofs.close();
};

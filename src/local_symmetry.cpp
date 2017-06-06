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

#define DEBUG
//#define NEW_APPLY_SYMMETRY_METHOD

static std::string str_new_mask = "NEW_MASK_AND_OPERATORS";
static std::string str_mask_filename = "MASKFILENAME";

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
	VEC_ELEM(op, DX_POS) *= factor;
	VEC_ELEM(op, DY_POS) *= factor;
	VEC_ELEM(op, DZ_POS) *= factor;
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
			|| (!MD.containsLabel(EMDL_ORIENT_ORIGIN_X))
			|| (!MD.containsLabel(EMDL_ORIENT_ORIGIN_Y))
			|| (!MD.containsLabel(EMDL_ORIENT_ORIGIN_Z)) )
		REPORT_ERROR("ERROR: Some of the MetaDataLabels are missing in STAR file " + (std::string)(fn_info) + " !");

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
			std::cout << "   --> Operator #" << int(0) << " = " << RFLOAT(0.) << ", " << RFLOAT(0.)
					<< ", " << RFLOAT(0.) << "; " << RFLOAT(0.) << ", " << RFLOAT(0.) << ", " << RFLOAT(0.) << " (the original)"<< std::endl;
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
			MD.getValue(EMDL_ORIENT_ORIGIN_X, dx);
			MD.getValue(EMDL_ORIENT_ORIGIN_Y, dy);
			MD.getValue(EMDL_ORIENT_ORIGIN_Z, dz);

			// Re-calculate angles so that they follow the conventions in RELION!
			standardiseEulerAngles(aa, bb, gg, aa, bb, gg);
			Localsym_composeOperator(op, aa, bb, gg, dx, dy, dz);

			// Do nothing if it is an identical operator
			if (sameLocalsymOperators(op, op_i))
				continue;

			if (verb)
			{
				std::cout << "   --> Operator #" << (dummy.size() + 1) << " = " << aa << ", " << bb << ", " << gg << "; " << dx << ", " << dy << ", " << dz << std::endl;
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
	MD.addLabel(EMDL_ORIENT_ORIGIN_X);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Y);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Z);
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
			MD.setValue(EMDL_ORIENT_ORIGIN_X, angpix * VEC_ELEM(ops[imask][iop], DX_POS));
			MD.setValue(EMDL_ORIENT_ORIGIN_Y, angpix * VEC_ELEM(ops[imask][iop], DY_POS));
			MD.setValue(EMDL_ORIENT_ORIGIN_Z, angpix * VEC_ELEM(ops[imask][iop], DZ_POS));
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
	MD.addLabel(EMDL_ORIENT_ORIGIN_X);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Y);
	MD.addLabel(EMDL_ORIENT_ORIGIN_Z);
	MD.addLabel(EMDL_IMAGE_WEIGHT);

	for (int iop = 0; iop < op_samplings.size(); iop++)
	{
		if (VEC_XSIZE(op_samplings[iop]) != NR_LOCALSYM_PARAMETERS)
			REPORT_ERROR("ERROR: syntax errors in results!");

		MD.addObject();
		MD.setValue(EMDL_ORIENT_ROT, VEC_ELEM(op_samplings[iop], AA_POS));
		MD.setValue(EMDL_ORIENT_TILT, VEC_ELEM(op_samplings[iop], BB_POS));
		MD.setValue(EMDL_ORIENT_PSI, VEC_ELEM(op_samplings[iop], GG_POS));
		MD.setValue(EMDL_ORIENT_ORIGIN_X, angpix * VEC_ELEM(op_samplings[iop], DX_POS));
		MD.setValue(EMDL_ORIENT_ORIGIN_Y, angpix * VEC_ELEM(op_samplings[iop], DY_POS));
		MD.setValue(EMDL_ORIENT_ORIGIN_Z, angpix * VEC_ELEM(op_samplings[iop], DZ_POS));
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
				std::cout << "   --> Operator #" << int(0) << " = " << RFLOAT(0.) << ", " << RFLOAT(0.)
						<< ", " << RFLOAT(0.) << "; " << RFLOAT(0.) << ", " << RFLOAT(0.) << ", " << RFLOAT(0.) << " (the original)"<< std::endl;
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
					std::cout << "   --> Operator #" << (ops.size() + 1) << " = " << aa << ", " << bb << ", " << gg << "; " << dx << ", " << dy << ", " << dz << std::endl;

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

/*
void getMaxSizeOfMask(
		const MultidimArray<RFLOAT>& mask,
		long int& xcen,
		long int& ycen,
		long int& zcen,
		long int& newdim)
{
	RFLOAT val = 0.;
	long int xmin = 0, xmax = 0, ymin = 0, ymax = 0, zmin = 0, zmax = 0;
	const long int my_long_int_max = 99999999;

	xmin = ymin = zmin =  my_long_int_max;
	xmax = ymax = zmax = -my_long_int_max;
	xcen = ycen = zcen = newdim = -1;

	if ((NSIZE(mask) != 1) || (ZSIZE(mask) <= 1) || (YSIZE(mask) <= 1) || (XSIZE(mask) <= 1))
		REPORT_ERROR("ERROR: input mask is not 3D!");

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mask)
	{
		val = DIRECT_A3D_ELEM(mask, k, i, j);
		if ((val < -(XMIPP_EQUAL_ACCURACY)) || ((val - 1.) > (XMIPP_EQUAL_ACCURACY)))
			REPORT_ERROR("ERROR: mask values are not in range [0,1]!");

		if (val > (XMIPP_EQUAL_ACCURACY))
		{
			if (k < zmin)
				zmin = k;
			if (k > zmax)
				zmax = k;
			if (j < ymin)
				ymin = j;
			if (j > ymax)
				ymax = j;
			if (i < xmin)
				xmin = i;
			if (i > xmax)
				xmax = i;
		}
	}

	if ( (xmin == my_long_int_max) || (ymin == my_long_int_max) || (zmin == my_long_int_max)
			|| (xmax == (-my_long_int_max)) || (ymax == (-my_long_int_max)) || (zmax == (-my_long_int_max)) )
		REPORT_ERROR("ERROR: mask is empty!");

	// TODO: check this re-centering is doing the right thing!!!
	if ((xmax - xmin) % 2 == 0)
		xmax++;
	xcen = FIRST_XMIPP_INDEX(xmax - xmin + 1) + xmax;
	if ((ymax - ymin) % 2 == 0)
		ymax++;
	ycen = FIRST_XMIPP_INDEX(ymax - ymin + 1) + ymax;
	if ((zmax - zmin) % 2 == 0)
		zmax++;
	zcen = FIRST_XMIPP_INDEX(zmax - zmin + 1) + zmax;

	newdim = CEIL(sqrt(RFLOAT((zmax - zmin) * (zmax - zmin) + (ymax - ymin) * (ymax - ymin) + (xmax - xmin) * (xmax - xmin))));
	if (newdim % 2)
		newdim++;
}
*/

void getSmallerSizeOfMask(
		const MultidimArray<RFLOAT>& mask,
		const Matrix1D<RFLOAT>& op_search_ranges,
		long int& newdim)
{
	RFLOAT xx = 0., yy = 0., zz = 0., xinit = 0., yinit = 0., zinit = 0., val = 0., rr = 0., rrmax = 0.;
	long int calcdim = 0;

	if ((NSIZE(mask) != 1) || (ZSIZE(mask) <= 1) || (YSIZE(mask) <= 1) || (XSIZE(mask) <= 1)
			|| (ZSIZE(mask) != YSIZE(mask)) || (ZSIZE(mask) != XSIZE(mask))
			|| (ZSIZE(mask) % 2) )
		REPORT_ERROR("ERROR: Input mask is not 3D cube!");

	newdim = ZSIZE(mask);

    zinit = RFLOAT(FIRST_XMIPP_INDEX(ZSIZE(mask)));
    yinit = RFLOAT(FIRST_XMIPP_INDEX(YSIZE(mask)));
    xinit = RFLOAT(FIRST_XMIPP_INDEX(XSIZE(mask)));
    rrmax = -1.;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mask)
	{
		val = DIRECT_A3D_ELEM(mask, k, i, j);
		if ((val < -(XMIPP_EQUAL_ACCURACY)) || ((val - 1.) > (XMIPP_EQUAL_ACCURACY)))
			REPORT_ERROR("ERROR: mask values are not in range [0,1]!");

		if (val > (XMIPP_EQUAL_ACCURACY))
		{
			zz = RFLOAT(k) + zinit;
			yy = RFLOAT(i) + yinit;
			xx = RFLOAT(j) + xinit;

			rr = xx * xx + yy * yy + zz * zz;
			if (rr > rrmax)
				rrmax = rr;
		}
	}
	rrmax = sqrt(rrmax);

	Localsym_decomposeOperator(op_search_ranges, val, val, val, xx, yy, zz, val);
	rrmax += sqrt(xx * xx + yy * yy + zz * zz);
	calcdim = (long int)(ceil(rrmax));
	if (calcdim % 2)
		calcdim++;
	if (calcdim < newdim)
		newdim = calcdim;
}

bool compareOperatorsByCC(
		const Matrix1D<RFLOAT>& lhs,
		const Matrix1D<RFLOAT>& rhs)
{
	return (VEC_ELEM(lhs, CC_POS) < VEC_ELEM(rhs, CC_POS));
}

bool useHealpixAngularSamplings(
		const Matrix1D<RFLOAT>& op,
		const Matrix1D<RFLOAT>& op_search_ranges,
		RFLOAT use_healpix_tilt_min)
{
	RFLOAT aa = 0., bb = 0., gg = 0., aa_range = 0., bb_range = 0., gg_range = 0., dummy = 0.;

	if (use_healpix_tilt_min < (XMIPP_EQUAL_ACCURACY))
		return true;
	if (use_healpix_tilt_min > (90. - (XMIPP_EQUAL_ACCURACY)))
		return false;

	Localsym_decomposeOperator(op_search_ranges, aa_range, bb_range, gg_range, dummy, dummy, dummy, dummy);

	// Always use Healpix if all rot and tilt angles are to be searched
	aa_range = (aa_range > 180.) ? (-1.) : (aa_range);
	bb_range = (bb_range >  90.) ? (-1.) : (bb_range);
	if ( (aa_range < 0.) && (bb_range < 0.) )
		return true;

	Localsym_decomposeOperator(op, aa, bb, gg, dummy, dummy, dummy, dummy);
	standardiseEulerAngles(aa, bb, gg, aa, bb, gg);
	if ( (bb > use_healpix_tilt_min) && (bb < (180. - use_healpix_tilt_min)) )
		return true;

	return false;
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

	if ( (ang_search_step < 0.01) || (ang_search_step > 30.) )
		REPORT_ERROR("ERROR: Angular searching step should be within range (+0.01, +30.00) degrees!");
	if ( (trans_search_step < 0.01) || (trans_search_step > 5.) )
		REPORT_ERROR("ERROR: Translational searching step should be within range (+0.01, +5.00) rescaled / binned pixels!");

	Localsym_decomposeOperator(op_old, aa_init, bb_init, gg_init, dx_init, dy_init, dz_init, cc);
	Localsym_decomposeOperator(op_search_ranges, aa_range, bb_range, gg_range, dx_range, dy_range, dz_range, cc);

	// Angular searching ranges
	if (!use_healpix)
	{
		aa_range = ( (aa_range > 180.) || (aa_range < 0.) ) ? (180.) : aa_range;
		bb_range = ( (bb_range >  90.) || (bb_range < 0.) ) ? ( 90.) : bb_range;
		gg_range = ( (gg_range > 180.) || (gg_range < 0.) ) ? (180.) : gg_range;
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
		std::cout << " + Local searches of local symmetry operator: Angles = ("
				<< aa_init << ", " << bb_init << ", " << gg_init << "), Translations (rescaled, binned) = ("
				<< dx_init << ", " << dy_init << ", " << dz_init << ") ..." << std::endl;
		std::cout << " + Generating sampling points with ranges: Angles = +/- ("
				<< aa_range << ", " << bb_range << ", " << gg_range << "), Translations (rescaled, binned) = +/- ("
				<< dx_range << ", " << dy_range << ", " << dz_range << ")." << std::endl;
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
		sampling.orientational_prior_mode = prior_mode;
		sampling.is_3D = sampling.is_3d_trans = true;
		sampling.limit_tilt = -91.; // Don't limit tilts
		sampling.psi_step = 360. / (6. * ROUND(std::pow(2., healpix_order)));
		sampling.offset_range = sampling.offset_step = 1.; // I don't use Healpix translational samplings
		sampling.random_perturbation = sampling.perturbation_factor = 0.;

		// Get all orientations
		sampling.initialise(prior_mode, 3, true, false, false, (prior_mode == NOPRIOR) ? (false) : (true));
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
			std::cout << std::endl << "  PSI = " << std::flush;
			for (int ii = 0; ii < pointer_psi_nonzeroprior.size(); ii++)
				std::cout << sampling.psi_angles[pointer_psi_nonzeroprior[ii]] << ", " << std::flush;
		}
		else
		{
			std::cout << std::endl << "  ROT = " << std::flush;
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
		std::cout << std::endl;
		std::cout << "  NR_TOTAL_DIR = " << nr_dir << ", NR_TOTAL_TRANS <= " << dxs.size() * dys.size() * dzs.size() << std::endl;
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
		std::cout << " + Total sampling points = " << op_samplings.size() << std::endl;

	if (op_samplings.size() < 1)
		REPORT_ERROR("ERROR: No sampling points!");
}

void checkSamplingRatesForMask(
		const MultidimArray<RFLOAT>& mask,
		RFLOAT ang_search_step,
		RFLOAT trans_search_step)
{
	RFLOAT mask_val = 0., rr = 0., rrmax = -1., xinit = 0., yinit = 0., zinit = 0., xx = 0., yy = 0., zz = 0.;

	if ((NSIZE(mask) != 1) || (ZSIZE(mask) <= 1) || (YSIZE(mask) <= 1) || (XSIZE(mask) <= 1))
		REPORT_ERROR("ERROR: Input mask is not 3D!");
	if ( (ang_search_step < 0.01) || (trans_search_step < 0.01) )
		REPORT_ERROR("ERROR: Angular / translational samplings are smaller than 0.01!");

    zinit = RFLOAT(FIRST_XMIPP_INDEX(ZSIZE(mask)));
    yinit = RFLOAT(FIRST_XMIPP_INDEX(YSIZE(mask)));
    xinit = RFLOAT(FIRST_XMIPP_INDEX(XSIZE(mask)));

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mask)
	{
		mask_val = DIRECT_A3D_ELEM(mask, k, i, j);
		if ((mask_val < -(XMIPP_EQUAL_ACCURACY)) || ((mask_val - 1.) > (XMIPP_EQUAL_ACCURACY)))
			REPORT_ERROR("ERROR: mask values are not in range [0,1]!");

		if (mask_val > (XMIPP_EQUAL_ACCURACY))
		{
			zz = RFLOAT(k) + zinit;
			yy = RFLOAT(i) + yinit;
			xx = RFLOAT(j) + xinit;

			rr = xx * xx + yy * yy + zz * zz;
			if (rr > rrmax)
				rrmax = rr;
		}
	}
	rrmax = sqrt(rrmax);

	std::cout << " + Angular search step = " << ang_search_step << " degree(s) is equivalent to translational search step = "
			<< rrmax * tan(DEG2RAD(ang_search_step)) << " pixel(s) with the current rescaled (binned) mask." << std::endl;
	std::cout << " + You need an angular step of " << RAD2DEG(atan2(trans_search_step, rrmax)) << " degree(s) to match translational search step = "
			<< trans_search_step << " pixel(s) with the current rescaled (binned) mask." << std::endl;
}

void calculateOperatorCC(
		const MultidimArray<RFLOAT>& map,
		const MultidimArray<RFLOAT>& mask,
		std::vector<Matrix1D<RFLOAT> >& op_samplings,
		bool do_sort,
		bool verb)
{
	RFLOAT val = 0., mask_val = 0., mask_val_sum = 0., cc = 0.;
	int barstep = 0, updatebar = 0, totalbar = 0;

	Matrix1D<RFLOAT> op_tmp;
	Matrix2D<RFLOAT> op_mat;
	MultidimArray<RFLOAT> vol;

	if (op_samplings.size() < 1)
		REPORT_ERROR("ERROR: No sampling points!");

	if ((NSIZE(map) != 1) || (ZSIZE(map) <= 1) || (YSIZE(map) <= 1) || (XSIZE(map) <= 1))
		REPORT_ERROR("ERROR: Input map is not 3D!");
	if ((NSIZE(map) != NSIZE(mask)) || (ZSIZE(map) != ZSIZE(mask)) || (YSIZE(map) != YSIZE(mask)) || (XSIZE(map) != XSIZE(mask)))
		REPORT_ERROR("ERROR: Input map and mask should have the same sizes!");

	// Check the mask, calculate the sum of mask values
	mask_val_sum = 0.;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mask)
	{
		mask_val = DIRECT_A3D_ELEM(mask, k, i, j);
		if ((mask_val < -(XMIPP_EQUAL_ACCURACY)) || ((mask_val - 1.) > (XMIPP_EQUAL_ACCURACY)))
			REPORT_ERROR("ERROR: mask values are not in range [0,1]!");

		if (mask_val > (XMIPP_EQUAL_ACCURACY))
			mask_val_sum += mask_val;
	}
	if (mask_val_sum < 1.)
		std::cout << " + WARNING: sum of mask values is smaller than 1! Please check whether it is a correct mask!" << std::endl;

	// Calculate all CCs
	if (verb)
	{
		std::cout << " + Calculate CCs for all sampling points ..." << std::endl;
		init_progress_bar(op_samplings.size());
		barstep = op_samplings.size() / 100;
		updatebar = totalbar = 0;
	}
	for (int iop = 0; iop < op_samplings.size(); iop++)
	{
		Localsym_operator2matrix(op_samplings[iop], op_mat, LOCALSYM_OP_DO_INVERT);

		applyGeometry(map, vol, op_mat, IS_NOT_INV, DONT_WRAP);

		cc = 0.;
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(vol)
		{
			mask_val = DIRECT_A3D_ELEM(mask, k, i, j);
			if (mask_val < XMIPP_EQUAL_ACCURACY)
				continue;

			val = DIRECT_A3D_ELEM(vol, k, i, j) - DIRECT_A3D_ELEM(map, k, i, j);
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

void local_symmetry_parameters::initBoolOptions()
{
	show_usage_for_an_option = false;

	do_apply_local_symmetry = false;
	do_duplicate_local_symmetry = false;
	do_local_search_local_symmetry_ops = false;
	do_txt2rln = false;
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
	show_usage_for_an_option = parser.checkOption("--help", "Show usage for the selected function (MAY 05, 2017)");

	int options_section = parser.addSection("Options");
	do_apply_local_symmetry = parser.checkOption("--apply", "Apply local symmetry to a 3D cryo-EM density map");
	do_duplicate_local_symmetry = parser.checkOption("--duplicate", "Duplicate subunits/masks according to local symmetry operators");
	do_local_search_local_symmetry_ops = parser.checkOption("--search", "Local searches of local symmetry operators");
	do_txt2rln = parser.checkOption("--txt2rln", "Convert operators from DM to RELION STAR format");
	do_debug = parser.checkOption("--debug", "(DEBUG ONLY)");

	int params_section = parser.addSection("Parameters (alphabetically ordered)");
	//angpix_ops = textToFloat(parser.getOption("--angpix_ops", "Pixel size (in Angstroms) of local symmetry operators", "-1."));
	angpix_image = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms) of input image", "-1."));
	ang_range = textToFloat(parser.getOption("--ang_range", "Angular search range of operators (in degrees), overwrite rot-tilt-psi ranges if set to positive", "2."));
	ang_rot_range = textToFloat(parser.getOption("--ang_rot_range", "Angular (rot) search range of operators (in degrees)", "2."));
	ang_tilt_range = textToFloat(parser.getOption("--ang_tilt_range", "Angular (tilt) search range of operators (in degrees)", "2."));
	ang_psi_range = textToFloat(parser.getOption("--ang_psi_range", "Angular (psi) search range of operators (in degrees)", "2."));
	ang_step = textToFloat(parser.getOption("--ang_step", "Angular search step of operators (in degrees)", "1."));
	binning_factor = textToFloat(parser.getOption("--bin", "Binning factor (< 1 means no binning)", "-1."));
	fn_unsym = parser.getOption("--i_map", "Input 3D unsymmetised map", "");
	fn_info_in = parser.getOption("--i_mask_info", "Input file with mask filenames and rotational / translational operators", "maskinfo.txt");
	offset_range = textToFloat(parser.getOption("--offset_range", "Translational search range of operators (in Angstroms), overwrite x-y-z ranges if set to positive", "2."));
	offset_x_range = textToFloat(parser.getOption("--offset_x_range", "Translational (x) search range of operators (in Angstroms)", "2."));
	offset_y_range = textToFloat(parser.getOption("--offset_y_range", "Translational (y) search range of operators (in Angstroms)", "2."));
	offset_z_range = textToFloat(parser.getOption("--offset_z_range", "Translational (z) search range of operators (in Angstroms)", "2."));
	offset_step = textToFloat(parser.getOption("--offset_step", "Translational search step of operators (in Angstroms)", "1."));
	fn_sym = parser.getOption("--o_map", "Output 3D symmetised map", "");
	fn_info_out = parser.getOption("--o_mask_info", "Output file with mask filenames and rotational / translational operators", "maskinfo_refined.txt");
	sphere_percentage = textToFloat(parser.getOption("--sphere_percentage", "Diameter of spherical mask divided by the box size (< 0.99)", "-1."));
	verb = parser.checkOption("--verb", "Verbose output?");

	int expert_section = parser.addSection("Parameters (expert options - alphabetically ordered)");
	fn_mask = parser.getOption("--i_mask", "(DEBUG) Input mask", "mask.mrc");
	fn_info_in_parsed_ext = parser.getOption("--i_mask_info_parsed_ext", "Extension of parsed input file with mask filenames and rotational / translational operators", "parsed");
	use_healpix_tilt_min = textToFloat(parser.getOption("--use_healpix_tilt_min", "Don't use Healpix for angular samplings if the original tilt value is <X or >(180-X), always use Healpix if set to negative", "30."));
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
			std::cout << "  USAGE: " << std::endl;
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
			std::cout << " Local searches of local symmetry operators" << std::endl;
			std::cout << "  USAGE: --search --i_map unsym.mrc --i_mask_info maskinfo.txt --o_mask_info maskinfo_new.txt --ang_range 2 --ang_step 0.5 --offset_range 2 --offset_step 1" << std::endl;
			displayEmptyLine();
			return;
		}

		FileName fn_searched_op_samplings;
		Matrix1D<RFLOAT> op_search_ranges;
		std::vector<Matrix1D<RFLOAT> > op_samplings;
		Image<RFLOAT> map, mask;
		long int olddim = 0, newdim = 0;
		RFLOAT mask_val = 0.;
		bool use_healpix = true;

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

		// Read the unsymmetrised map
		map.read(fn_unsym);
		if ((NSIZE(map()) != 1) || (ZSIZE(map()) <= 1) || (YSIZE(map()) <= 1) || (XSIZE(map()) <= 1)
				|| (ZSIZE(map()) != YSIZE(map())) || (ZSIZE(map()) != XSIZE(map()))
				|| (ZSIZE(map()) % 2) )
			REPORT_ERROR("ERROR: Input map is not 3D cube!");
		olddim = newdim = ZSIZE(map());
		if ((binning_factor - 1.) > XMIPP_EQUAL_ACCURACY)
		{
			newdim = (long int)(ceil(RFLOAT(olddim) / binning_factor));
			if (newdim < 2)
				REPORT_ERROR("ERROR: Binning factor is too large!");
			if ((newdim + 1) < olddim) // Need rescaling
			{
				// Dimension should always be even
				if (newdim % 2)
					newdim++;
				resizeMap(map(), newdim);
			}
			else
				newdim = olddim;
			binning_factor = RFLOAT(olddim) / RFLOAT(newdim);
			if (newdim != olddim)
				std::cout << " + Rescale box size from " << olddim << " to " << newdim << ". Binning factor = " << binning_factor << std::endl;
		}
		else
			binning_factor = 1.;

		// Get search ranges
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

		// Rescale translational search ranges and steps
		if (newdim != olddim)
		{
			Localsym_scaleTranslations(op_search_ranges, 1. / binning_factor);
			offset_step /= binning_factor;

			for (int imask = 0; imask < op_list.size(); imask++)
			{
				for (int iop = 0; iop < op_list[imask].size(); iop++)
					Localsym_scaleTranslations(op_list[imask][iop], 1. / binning_factor);
			}
		}

		// Loop over all masks
		for (int imask = 0; imask < fn_mask_list.size(); imask++)
		{
			// Read the mask
			mask.read(fn_mask_list[imask]);
			if ((NSIZE(mask()) != 1) || (ZSIZE(mask()) <= 1) || (YSIZE(mask()) <= 1) || (XSIZE(mask()) <= 1)
					|| (ZSIZE(mask()) != YSIZE(mask())) || (ZSIZE(mask()) != XSIZE(mask()))
					|| (ZSIZE(mask()) % 2) )
				REPORT_ERROR("ERROR: Input mask is not 3D cube!");
			if (olddim != ZSIZE(mask()))
				REPORT_ERROR("ERROR: Input map and masks should have the same size!");
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mask())
			{
				mask_val = DIRECT_A3D_ELEM(mask(), k, i, j);
				if ((mask_val < -(XMIPP_EQUAL_ACCURACY)) || ((mask_val - 1.) > (XMIPP_EQUAL_ACCURACY)))
					REPORT_ERROR("ERROR: mask " + std::string(fn_mask_list[imask]) + " - values are not in range [0,1]!");
			}
			if (newdim != olddim)
			{
				resizeMap(mask(), newdim);

				// Mask values might go out of range after rescaling. Fix it if it happens
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(mask())
				{
					mask_val = DIRECT_A3D_ELEM(mask(), k, i, j);
					if (mask_val < 0.)
						DIRECT_A3D_ELEM(mask(), k, i, j) = 0.;
					if (mask_val > 1.)
						DIRECT_A3D_ELEM(mask(), k, i, j) = 1.;
				}
			}

			checkSamplingRatesForMask(mask(), ang_step, offset_step);

			// Loop over all operators of this mask
			for (int iop = 0; iop < op_list[imask].size(); iop++)
			{
				// Check whether use Healpix
				use_healpix = useHealpixAngularSamplings(op_list[imask][iop], op_search_ranges, use_healpix_tilt_min);

				// Get sampling points
				getLocalSearchOperatorSamplings(
						op_list[imask][iop],
						op_search_ranges,
						op_samplings,
						ang_step,
						offset_step,
						use_healpix,
						true);

				// Calculate all CCs for the sampling points
				calculateOperatorCC(map(), mask(), op_samplings, do_sort, do_verb);

				// For rescaled maps
				if (newdim != olddim)
				{
					for (int isamp = 0; isamp < op_samplings.size(); isamp++)
						Localsym_scaleTranslations(op_samplings[isamp], binning_factor);
				}

				// Output the local searches results
				fn_tmp.compose("cc_mask", imask + 1, "tmp", 3);  // "cc_mask001.tmp"
				fn_tmp = fn_tmp.withoutExtension(); // "cc_mask001"
				fn_searched_op_samplings.compose(fn_tmp + "_op", iop + 1, "star", 3); // "cc_mask001_op001.star"
				writeRelionFormatLocalSearchOperatorResults(fn_searched_op_samplings, op_samplings, angpix_image);

				// Update this operator and do screen output
				op_list[imask][iop] = op_samplings[0];
				std::cout << " + Done! Local refined local symmetry operator: Angles = ("
						<< VEC_ELEM(op_samplings[0], AA_POS) << ", " << VEC_ELEM(op_samplings[0], BB_POS) << ", " << VEC_ELEM(op_samplings[0], GG_POS) << "), Translations = ("
						<< angpix_image * VEC_ELEM(op_samplings[0], DX_POS) << ", " << angpix_image * VEC_ELEM(op_samplings[0], DY_POS) << ", " << angpix_image * VEC_ELEM(op_samplings[0], DZ_POS) << ")." << std::endl;
			}
		}

		// Write out new mask info file
		if (fn_info_out.getExtension() == "star")
			writeRelionFormatMasksAndOperators(fn_info_out, fn_mask_list, op_list, angpix_image);
		else
			writeDMFormatMasksAndOperators(fn_info_out, fn_mask_list, op_list, angpix_image);
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
	else if (do_debug)
	{
		Image<RFLOAT> img1;
		RFLOAT mask_val = 0.;
		img1.read(fn_mask);
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(img1())
		{
			mask_val = DIRECT_A3D_ELEM(img1(), k, i, j);
			if (mask_val > 1.)
				DIRECT_A3D_ELEM(img1(), k, i, j) = 1.;
			else if (mask_val < 0.)
				DIRECT_A3D_ELEM(img1(), k, i, j) = 0.;
		}
		img1.write(fn_sym);
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
}

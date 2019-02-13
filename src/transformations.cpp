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
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *              Sjors H.W. Scheres
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/transformations.h"

/* Rotation 2D ------------------------------------------------------------- */
void rotation2DMatrix(RFLOAT ang, Matrix2D< RFLOAT > &result, bool homogeneous)
{
    RFLOAT cosine, sine;

    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    if (homogeneous)
    {
        if (MAT_XSIZE(result)!=3 || MAT_YSIZE(result)!=3)
            result.resize(3,3);
        MAT_ELEM(result,0, 2) = 0;
        MAT_ELEM(result,1, 2) = 0;
        MAT_ELEM(result,2, 0) = 0;
        MAT_ELEM(result,2, 1) = 0;
        MAT_ELEM(result,2, 2) = 1;
    }
    else
        if (MAT_XSIZE(result)!=2 || MAT_YSIZE(result)!=2)
            result.resize(2,2);

    MAT_ELEM(result,0, 0) = cosine;
    MAT_ELEM(result,0, 1) = -sine;

    MAT_ELEM(result,1, 0) = sine;
    MAT_ELEM(result,1, 1) = cosine;
}

/* Translation 2D ---------------------------------------------------------- */
void translation2DMatrix(const Matrix1D<RFLOAT> &v,
                         Matrix2D< RFLOAT > &result)
{
    // if (VEC_XSIZE(v) != 2)
    //    REPORT_ERROR("Translation2D_matrix: vector is not in R2");

    result.initIdentity(3);
    MAT_ELEM(result,0, 2) = XX(v);
    MAT_ELEM(result,1, 2) = YY(v);
}

/* Rotation 3D around the system axes -------------------------------------- */
void rotation3DMatrix(RFLOAT ang, char axis, Matrix2D< RFLOAT > &result,
                      bool homogeneous)
{
    if (homogeneous)
    {
        result.initZeros(4,4);
        MAT_ELEM(result,3, 3) = 1;
    }
    else
        result.initZeros(3,3);

    RFLOAT cosine, sine;
    ang = DEG2RAD(ang);
    cosine = cos(ang);
    sine = sin(ang);

    switch (axis)
    {
    case 'Z':
        MAT_ELEM(result,0, 0) = cosine;
        MAT_ELEM(result,0, 1) = -sine;
        MAT_ELEM(result,1, 0) = sine;
        MAT_ELEM(result,1, 1) = cosine;
        MAT_ELEM(result,2, 2) = 1;
        break;
    case 'Y':
        MAT_ELEM(result,0, 0) = cosine;
        MAT_ELEM(result,0, 2) = -sine;
        MAT_ELEM(result,2, 0) = sine;
        MAT_ELEM(result,2, 2) = cosine;
        MAT_ELEM(result,1, 1) = 1;
        break;
    case 'X':
        MAT_ELEM(result,1, 1) = cosine;
        MAT_ELEM(result,1, 2) = -sine;
        MAT_ELEM(result,2, 1) = sine;
        MAT_ELEM(result,2, 2) = cosine;
        MAT_ELEM(result,0, 0) = 1;
        break;
    default:
        REPORT_ERROR("rotation3DMatrix: Unknown axis");
    }
}

/* Align a vector with Z axis */
void alignWithZ(const Matrix1D<RFLOAT> &axis, Matrix2D<RFLOAT>& result,
                bool homogeneous)
{
    if (axis.size() != 3)
        REPORT_ERROR("alignWithZ: Axis is not in R3");
    if (homogeneous)
    {
        result.initZeros(4,4);
        MAT_ELEM(result,3, 3) = 1;
    }
    else
        result.initZeros(3,3);
    Matrix1D<RFLOAT>  Axis(axis);
    Axis.selfNormalize();

    // Compute length of the projection on YZ plane
    RFLOAT proj_mod = sqrt(YY(Axis) * YY(Axis) + ZZ(Axis) * ZZ(Axis));
    if (proj_mod > XMIPP_EQUAL_ACCURACY)
    {   // proj_mod!=0
        // Build Matrix result, which makes the turning axis coincident with Z
        MAT_ELEM(result,0, 0) = proj_mod;
        MAT_ELEM(result,0, 1) = -XX(Axis) * YY(Axis) / proj_mod;
        MAT_ELEM(result,0, 2) = -XX(Axis) * ZZ(Axis) / proj_mod;
        MAT_ELEM(result,1, 0) = 0;
        MAT_ELEM(result,1, 1) = ZZ(Axis) / proj_mod;
        MAT_ELEM(result,1, 2) = -YY(Axis) / proj_mod;
        MAT_ELEM(result,2, 0) = XX(Axis);
        MAT_ELEM(result,2, 1) = YY(Axis);
        MAT_ELEM(result,2, 2) = ZZ(Axis);
    }
    else
    {
        // I know that the Axis is the X axis, EITHER POSITIVE OR NEGATIVE!!
        MAT_ELEM(result,0, 0) = 0;
        MAT_ELEM(result,0, 1) = 0;
        MAT_ELEM(result,0, 2) = (XX(Axis) > 0)? -1 : 1;
        MAT_ELEM(result,1, 0) = 0;
        MAT_ELEM(result,1, 1) = 1;
        MAT_ELEM(result,1, 2) = 0;
        MAT_ELEM(result,2, 0) = (XX(Axis) > 0)? 1 : -1;
        MAT_ELEM(result,2, 1) = 0;
        MAT_ELEM(result,2, 2) = 0;
    }
}

/* Rotation 3D around any axis -------------------------------------------- */
void rotation3DMatrix(RFLOAT ang, const Matrix1D<RFLOAT> &axis,
                      Matrix2D<RFLOAT> &result, bool homogeneous)
{
    // Compute a matrix which makes the turning axis coincident with Z
    // And turn around this axis
    Matrix2D<RFLOAT> A,R;
    alignWithZ(axis,A,homogeneous);
    rotation3DMatrix(ang, 'Z', R, homogeneous);
    result=A.transpose() * R * A;
}

/* Translation 3D ---------------------------------------------------------- */
void translation3DMatrix(const Matrix1D<RFLOAT> &v, Matrix2D<RFLOAT> &result)
{
    if (VEC_XSIZE(v) != 3)
        REPORT_ERROR("Translation3D_matrix: vector is not in R3");

    result.initIdentity(4);
    MAT_ELEM(result,0, 3) = XX(v);
    MAT_ELEM(result,1, 3) = YY(v);
    MAT_ELEM(result,2, 3) = ZZ(v);
}

/* Scale 3D ---------------------------------------------------------------- */
void scale3DMatrix(const Matrix1D<RFLOAT> &sc, Matrix2D<RFLOAT>& result,
                   bool homogeneous)
{
    if (VEC_XSIZE(sc) != 3)
        REPORT_ERROR("Scale3D_matrix: vector is not in R3");

    if (homogeneous)
    {
        result.initZeros(4,4);
        MAT_ELEM(result,3, 3) = 1;
    }
    else
        result.initZeros(3,3);
    MAT_ELEM(result,0, 0) = XX(sc);
    MAT_ELEM(result,1, 1) = YY(sc);
    MAT_ELEM(result,2, 2) = ZZ(sc);
}

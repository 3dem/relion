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

#include "src/matrix1d.h"

Matrix1D<RFLOAT> vectorR2(RFLOAT x, RFLOAT y)
{
    Matrix1D<RFLOAT> result(2);
    result( 0) = x;
    result( 1) = y;
    return result;
}

Matrix1D<RFLOAT> vectorR3(RFLOAT x, RFLOAT y, RFLOAT z)
{
    Matrix1D<RFLOAT> result(3);
    result( 0) = x;
    result( 1) = y;
    result( 2) = z;
    return result;
}

// This function only makes sense after all code has been modified with 'sed' to allow single-precision runs
#ifdef RELION_SINGLE_PRECISION
Matrix1D<float> vectorR3(double xx, double yy, double zz)
{
	return vectorR3((float)xx, (float)yy, (float)zz);
}
#endif

Matrix1D<int> vectorR3(int x, int y, int z)
{
    Matrix1D<int> result(3);
    result( 0) = x;
    result( 1) = y;
    result( 2) = z;
    return result;
}

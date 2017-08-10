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

#include <iostream>
#include <math.h>

#include "src/euler.h"
#include "src/funcs.h"

/* Euler angles --> matrix ------------------------------------------------- */
void Euler_angles2matrix(RFLOAT alpha, RFLOAT beta, RFLOAT gamma,
                         Matrix2D<RFLOAT> &A, bool homogeneous)
{
    RFLOAT ca, sa, cb, sb, cg, sg;
    RFLOAT cc, cs, sc, ss;

    if (homogeneous)
    {
        A.initZeros(4,4);
        MAT_ELEM(A,3,3)=1;
    }
    else
        if (MAT_XSIZE(A) != 3 || MAT_YSIZE(A) != 3)
            A.resize(3, 3);

    alpha = DEG2RAD(alpha);
    beta  = DEG2RAD(beta);
    gamma = DEG2RAD(gamma);

    ca = cos(alpha);
    cb = cos(beta);
    cg = cos(gamma);
    sa = sin(alpha);
    sb = sin(beta);
    sg = sin(gamma);
    cc = cb * ca;
    cs = cb * sa;
    sc = sb * ca;
    ss = sb * sa;

    A(0, 0) =  cg * cc - sg * sa;
    A(0, 1) =  cg * cs + sg * ca;
    A(0, 2) = -cg * sb;
    A(1, 0) = -sg * cc - cg * sa;
    A(1, 1) = -sg * cs + cg * ca;
    A(1, 2) = sg * sb;
    A(2, 0) =  sc;
    A(2, 1) =  ss;
    A(2, 2) = cb;
}

/* Euler direction --------------------------------------------------------- */
void Euler_angles2direction(RFLOAT alpha, RFLOAT beta,
						    Matrix1D<RFLOAT> &v)
{
    RFLOAT ca, sa, cb, sb;
    RFLOAT sc, ss;

    v.resize(3);
    alpha = DEG2RAD(alpha);
    beta  = DEG2RAD(beta);

    ca = cos(alpha);
    cb = cos(beta);
    sa = sin(alpha);
    sb = sin(beta);
    sc = sb * ca;
    ss = sb * sa;

    v(0) = sc;
    v(1) = ss;
    v(2) = cb;
}

/* Euler direction2angles ------------------------------- */
//gamma is useless but I keep it for simmetry
//with Euler_direction
void Euler_direction2angles(Matrix1D<RFLOAT> &v0,
                            RFLOAT &alpha, RFLOAT &beta)
{
	// Aug25,2015 - Shaoda
	// This function can recover tilt (b) as small as 0.0001 degrees
	// It replaces a more complicated version in the code before Aug2015
    Matrix1D<RFLOAT> v;

    // Make sure the vector is normalised
    v.resize(3);
    v = v0;
    v.selfNormalize();

    // Tilt (b) should be [0, +180] degrees. Rot (a) should be [-180, +180] degrees
    alpha = RAD2DEG(atan2(v(1), v(0))); // 'atan2' returns an angle within [-pi, +pi] radians for rot
    beta = RAD2DEG(acos(v(2))); // 'acos' returns an angle within [0, +pi] radians for tilt

    // The following is done to keep in line with the results from old codes
    // If tilt (b) = 0 or 180 degrees, sin(b) = 0, rot (a) cannot be calculated from the direction
    if ( (fabs(beta) < 0.001) || (fabs(beta - 180.) < 0.001) )
    	alpha = 0.;

    return;

}

/* Matrix --> Euler angles ------------------------------------------------- */
#define CHECK
//#define DEBUG_EULER
void Euler_matrix2angles(const Matrix2D<RFLOAT> &A, RFLOAT &alpha,
                         RFLOAT &beta, RFLOAT &gamma)
{
    RFLOAT abs_sb, sign_sb;

    if (MAT_XSIZE(A) != 3 || MAT_YSIZE(A) != 3)
        REPORT_ERROR( "Euler_matrix2angles: The Euler matrix is not 3x3");

    abs_sb = sqrt(A(0, 2) * A(0, 2) + A(1, 2) * A(1, 2));
    if (abs_sb > 16*FLT_EPSILON)
    {
        gamma = atan2(A(1, 2), -A(0, 2));
        alpha = atan2(A(2, 1), A(2, 0));
        if (ABS(sin(gamma)) < FLT_EPSILON)
            sign_sb = SGN(-A(0, 2) / cos(gamma));
        // if (sin(alpha)<FLT_EPSILON) sign_sb=SGN(-A(0,2)/cos(gamma));
        // else sign_sb=(sin(alpha)>0) ? SGN(A(2,1)):-SGN(A(2,1));
        else
            sign_sb = (sin(gamma) > 0) ? SGN(A(1, 2)) : -SGN(A(1, 2));
        beta  = atan2(sign_sb * abs_sb, A(2, 2));
    }
    else
    {
        if (SGN(A(2, 2)) > 0)
        {
            // Let's consider the matrix as a rotation around Z
            alpha = 0;
            beta  = 0;
            gamma = atan2(-A(1, 0), A(0, 0));
        }
        else
        {
            alpha = 0;
            beta  = PI;
            gamma = atan2(A(1, 0), -A(0, 0));
        }
    }

    gamma = RAD2DEG(gamma);
    beta  = RAD2DEG(beta);
    alpha = RAD2DEG(alpha);

#ifdef DEBUG_EULER
    std::cout << "abs_sb " << abs_sb << std::endl;
    std::cout << "A(1,2) " << A(1, 2) << " A(0,2) " << A(0, 2) << " gamma "
    << gamma << std::endl;
    std::cout << "A(2,1) " << A(2, 1) << " A(2,0) " << A(2, 0) << " alpha "
    << alpha << std::endl;
    std::cout << "sign sb " << sign_sb << " A(2,2) " << A(2, 2)
    << " beta " << beta << std::endl;
#endif
}
#undef CHECK

#ifdef NEVERDEFINED
// Michael's method
void Euler_matrix2angles(Matrix2D<RFLOAT> A, RFLOAT *alpha, RFLOAT *beta,
                         RFLOAT *gamma)
{
    RFLOAT abs_sb;

    if (ABS(A(1, 1)) > FLT_EPSILON)
    {
        abs_sb = sqrt((-A(2, 2) * A(1, 2) * A(2, 1) - A(0, 2) * A(2, 0)) / A(1, 1));
    }
    else if (ABS(A(0, 1)) > FLT_EPSILON)
    {
        abs_sb = sqrt((-A(2, 1) * A(2, 2) * A(0, 2) + A(2, 0) * A(1, 2)) / A(0, 1));
    }
    else if (ABS(A(0, 0)) > FLT_EPSILON)
    {
        abs_sb = sqrt((-A(2, 0) * A(2, 2) * A(0, 2) - A(2, 1) * A(1, 2)) / A(0, 0));
    }
    else
        EXIT_ERROR(1, "Don't know how to extract angles");

    if (abs_sb > FLT_EPSILON)
    {
        *beta  = atan2(abs_sb, A(2, 2));
        *alpha = atan2(A(2, 1) / abs_sb, A(2, 0) / abs_sb);
        *gamma = atan2(A(1, 2) / abs_sb, -A(0, 2) / abs_sb);
    }
    else
    {
        *alpha = 0;
        *beta  = 0;
        *gamma = atan2(A(1, 0), A(0, 0));
    }

    *gamma = rad2deg(*gamma);
    *beta  = rad2deg(*beta);
    *alpha = rad2deg(*alpha);
}
#endif
/* Euler up-down correction ------------------------------------------------ */
void Euler_up_down(RFLOAT rot, RFLOAT tilt, RFLOAT psi,
                   RFLOAT &newrot, RFLOAT &newtilt, RFLOAT &newpsi)
{
    newrot  = rot;
    newtilt = tilt + 180;
    newpsi  = -(180 + psi);
}

/* Same view, differently expressed ---------------------------------------- */
void Euler_another_set(RFLOAT rot, RFLOAT tilt, RFLOAT psi,
                       RFLOAT &newrot, RFLOAT &newtilt, RFLOAT &newpsi)
{
    newrot  = rot + 180;
    newtilt = -tilt;
    newpsi  = -180 + psi;
}

/* Euler mirror Y ---------------------------------------------------------- */
void Euler_mirrorY(RFLOAT rot, RFLOAT tilt, RFLOAT psi,
                   RFLOAT &newrot, RFLOAT &newtilt, RFLOAT &newpsi)
{
    newrot  = rot;
    newtilt = tilt + 180;
    newpsi  = -psi;
}

/* Euler mirror X ---------------------------------------------------------- */
void Euler_mirrorX(RFLOAT rot, RFLOAT tilt, RFLOAT psi,
                   RFLOAT &newrot, RFLOAT &newtilt, RFLOAT &newpsi)
{
    newrot  = rot;
    newtilt = tilt + 180;
    newpsi  = 180 - psi;
}

/* Euler mirror XY --------------------------------------------------------- */
void Euler_mirrorXY(RFLOAT rot, RFLOAT tilt, RFLOAT psi,
                    RFLOAT &newrot, RFLOAT &newtilt, RFLOAT &newpsi)
{
    newrot  = rot;
    newtilt = tilt;
    newpsi  = 180 + psi;
}

/* Apply a transformation matrix to Euler angles --------------------------- */
void Euler_apply_transf(const Matrix2D<RFLOAT> &L,
                        const Matrix2D<RFLOAT> &R,
                        RFLOAT rot,
                        RFLOAT tilt,
                        RFLOAT psi,
                        RFLOAT &newrot,
                        RFLOAT &newtilt,
                        RFLOAT &newpsi)
{

    Matrix2D<RFLOAT> euler(3, 3), temp;
    Euler_angles2matrix(rot, tilt, psi, euler);
    temp = L * euler * R;
    Euler_matrix2angles(temp, newrot, newtilt, newpsi);
}

/* Rotate (3D) MultidimArray with 3 Euler angles ------------------------------------- */
void Euler_rotation3DMatrix(RFLOAT rot, RFLOAT tilt, RFLOAT psi, Matrix2D<RFLOAT> &result)
{
    Euler_angles2matrix(rot, tilt, psi, result, true);
}



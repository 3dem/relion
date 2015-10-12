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
    RFLOAT abs_ca, sb, cb;
    RFLOAT aux_alpha;
    RFLOAT aux_beta;
    RFLOAT error, newerror;
    Matrix1D<RFLOAT> v_aux;
    Matrix1D<RFLOAT> v;

    //if not normalized do it so
    v.resize(3);
    v = v0;
    v.selfNormalize();

    v_aux.resize(3);
    cb = v(2);

    if (fabs((cb)) > 0.999847695)/*one degree */
    {
        std::cerr << "\nWARNING: Routine Euler_direction2angles is not reliable\n"
        "for small tilt angles. Up to 0.001 deg it should be OK\n"
        "for most applications but you never know";
    }

    if (fabs((cb - 1.)) < FLT_EPSILON)
    {
        alpha = 0.;
        beta = 0.;
    }
    else
    {/*1*/

        aux_beta = acos(cb); /* beta between 0 and PI */


        sb = sin(aux_beta);

        abs_ca = fabs(v(0)) / sb;
        if (fabs((abs_ca - 1.)) < FLT_EPSILON)
            aux_alpha = 0.;
        else
            aux_alpha = acos(abs_ca);

        v_aux(0) = sin(aux_beta) * cos(aux_alpha);
        v_aux(1) = sin(aux_beta) * sin(aux_alpha);
        v_aux(2) = cos(aux_beta);

        error = fabs(dotProduct(v, v_aux) - 1.);
        alpha = aux_alpha;
        beta = aux_beta;

        v_aux(0) = sin(aux_beta) * cos(-1. * aux_alpha);
        v_aux(1) = sin(aux_beta) * sin(-1. * aux_alpha);
        v_aux(2) = cos(aux_beta);
        newerror = fabs(dotProduct(v, v_aux) - 1.);
        if (error > newerror)
        {
            alpha = -1. * aux_alpha;
            beta  = aux_beta;
            error = newerror;
        }

        v_aux(0) = sin(-aux_beta) * cos(-1. * aux_alpha);
        v_aux(1) = sin(-aux_beta) * sin(-1. * aux_alpha);
        v_aux(2) = cos(-aux_beta);
        newerror = fabs(dotProduct(v, v_aux) - 1.);
        if (error > newerror)
        {
            alpha = -1. * aux_alpha;
            beta  = -1. * aux_beta;
            error = newerror;
        }

        v_aux(0) = sin(-aux_beta) * cos(aux_alpha);
        v_aux(1) = sin(-aux_beta) * sin(aux_alpha);
        v_aux(2) = cos(-aux_beta);
        newerror = fabs(dotProduct(v, v_aux) - 1.);

        if (error > newerror)
        {
            alpha = aux_alpha;
            beta  = -1. * aux_beta;
            error = newerror;
        }
    }/*else 1 end*/
    beta  = RAD2DEG(beta);
    alpha = RAD2DEG(alpha);
}/*Eulerdirection2angles end*/

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

#ifdef RFLOAT

    Matrix2D<RFLOAT> Ap;
    Euler_angles2matrix(alpha, beta, gamma, Ap);
    if (A != Ap)
    {
        std::cout << "---\n";
        std::cout << "Euler_matrix2angles: I have computed angles "
        " which doesn't match with the original matrix\n";
        std::cout << "Original matrix\n" << A;
        std::cout << "Computed angles alpha=" << alpha << " beta=" << beta
        << " gamma=" << gamma << std::endl;
        std::cout << "New matrix\n" << Ap;
        std::cout << "---\n";
    }
#endif

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
#undef DEBUG

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



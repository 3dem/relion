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

#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "src/multidim_array.h"
#include "src/transformations.h"

#ifndef FLT_EPSILON
#define FLT_EPSILON 1.19209e-07
#endif

/// @name Euler operations
/// @{

/** Euler angles --> "Euler" matrix
 *
 * This function returns the transformation matrix associated to the 3 given
 * Euler angles (in degrees).
 *
 * As an implementation note you might like to know that this function calls
 * always to Matrix2D::resize
 *
 * See http://xmipp.cnb.csic.es/twiki/bin/view/Xmipp/EulerAngles for a
 * description of the Euler angles.
 */
void Euler_angles2matrix(RFLOAT a, RFLOAT b, RFLOAT g, Matrix2D< RFLOAT >& A,
                         bool homogeneous=false);

/** Euler angles2direction
 *
 * This function returns  a vector parallel to the  projection direction.
 * Resizes v if needed
 */
void Euler_angles2direction(RFLOAT alpha, RFLOAT beta, Matrix1D< RFLOAT >& v);

/** Euler direction2angles
 *
 * This function returns the 2 Euler angles (rot&tilt) associated to the direction given by
 * the vector v.
 */
void Euler_direction2angles(Matrix1D< RFLOAT >& v,
                            RFLOAT& alpha,
                            RFLOAT& beta);

/** "Euler" matrix --> angles
 *
 * This function compute a set of Euler angles which result in an "Euler" matrix
 * as the one given. See \ref Euler_angles2matrix to know more about how this
 * matrix is computed and what each row means. The result angles are in degrees.
 * Alpha, beta and gamma are respectively the first, second and third rotation
 * angles. If the input matrix is not 3x3 then an exception is thrown, the
 * function doesn't check that the Euler matrix is truly representing a
 * coordinate system.
 *
 * @code
 * Euler_matrix2angles(Euler, alpha, beta, gamma);
 * @endcode
 */
void Euler_matrix2angles(const Matrix2D< RFLOAT >& A,
                         RFLOAT& alpha,
                         RFLOAT& beta,
                         RFLOAT& gamma);

/** Up-Down projection equivalence
 *
 * As you know a projection view from a point has got its homologous from its
 * diametrized point in the projection sphere. This function takes a projection
 * defined by its 3 Euler angles and computes an equivalent set of Euler angles
 * from which the view is exactly the same but in the other part of the sphere
 * (if the projection is taken from the bottom then the new projection from the
 * top, and viceversa). The defined projections are exactly the same except for
 * a flip over X axis, ie, an up-down inversion. Exactly the correction
 * performed is:
 *
 * @code
 * newrot = rot;
 * newtilt = tilt + 180;
 * newpsi = -(180 + psi);
 * @endcode
 *
 * @code
 * Euler_up_down(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_up_down(RFLOAT rot,
                   RFLOAT tilt,
                   RFLOAT psi,
                   RFLOAT& newrot,
                   RFLOAT& newtilt,
                   RFLOAT& newpsi);

/** The same view but differently expressed
 *
 * As you know a projection view from a point can be expressed with different
 * sets of Euler angles. This function gives you another expression of the Euler
 * angles for this point of view. Exactly the operation performed is:
 *
 * @code
 * newrot = rot + 180;
 * newtilt = -tilt;
 * newpsi = -180 + psi;
 * @endcode
 *
 * @code
 * Euler_another_set(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_another_set(RFLOAT rot,
                       RFLOAT tilt,
                       RFLOAT psi,
                       RFLOAT& newrot,
                       RFLOAT& newtilt,
                       RFLOAT& newpsi);

/** Mirror over Y axis
 *
 * Given a set of Euler angles this function returns a new set which define a
 * mirrored (over Y axis) version of the former projection.
 *
 * @code
 *  -----> X               X<------
 *  |                              |
 *  |                              |
 *  |               ======>        |
 *  v                              v
 *  Y                             Y
 * @endcode
 *
 * The operation performed is
 *
 * @code
 * newrot = rot;
 * newtilt = tilt + 180;
 * newpsi = -psi;
 * @endcode
 *
 * @code
 * Euler_mirrorY(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_mirrorY(RFLOAT rot,
                   RFLOAT tilt,
                   RFLOAT psi,
                   RFLOAT& newrot,
                   RFLOAT& newtilt,
                   RFLOAT& newpsi);

/** Mirror over X axis
 *
 * Given a set of Euler angles this function returns a new set which define a
 * mirrored (over X axis) version of the former projection.
 *
 * @code
 *  -----> X               Y
 *  |                       ^
 *  |                       |
 *  |               ======> |
 *  v                       |
 *  Y                        -----> X
 * @endcode
 *
 * The operation performed is
 *
 * @code
 * newrot = rot;
 * newtilt = tilt + 180;
 * newpsi = 180 - psi;
 * @endcode
 *
 * @code
 * Euler_mirrorX(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_mirrorX(RFLOAT rot,
                   RFLOAT tilt,
                   RFLOAT psi,
                   RFLOAT& newrot,
                   RFLOAT& newtilt,
                   RFLOAT& newpsi);

/** Mirror over X and Y axes
 *
 * Given a set of Euler angles this function returns a new set which define a
 * mirrored (over X and Y axes at the same time) version of the former
 * projection.
 *
 * @code
 *  -----> X                       Y
 *  |                               ^
 *  |                               |
 *  |               ======>         |
 *  v                               |
 *  Y                        X<-----
 * @endcode
 *
 * The operation performed is
 *
 * @code
 * newrot = rot;
 * newtilt = tilt;
 * newpsi = 180 + psi;
 * @endcode
 *
 * @code
 * Euler_mirrorX(rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_mirrorXY(RFLOAT rot,
                    RFLOAT tilt,
                    RFLOAT psi,
                    RFLOAT& newrot,
                    RFLOAT& newtilt,
                    RFLOAT& newpsi);

/** Apply a geometrical transformation
 *
 * The idea behind this function is the following. 3 Euler angles define a point
 * of view for a projection, but also a coordinate system. You might apply a
 * geometrical transformation to this system, and then compute back what the
 * Euler angles for the new system are. This could be used to "mirror" points of
 * view, rotate them and all the stuff. The transformation matrix must be 3x3
 * but it must transform R3 vectors into R3 vectors (that is a normal 3D
 * transformation matrix when vector coordinates are not homogeneous) and it
 * will be applied in the sense:
 *
 * @code
 * New Euler matrix = L * Old Euler matrix * R
 * @endcode
 *
 * where you know that the Euler matrix rows represent the different system
 * axes. See Euler_angles2matrix for more information about the Euler coordinate
 * system.
 *
 * @code
 * Matrix2D< RFLOAT > R60 = rotation3DMatrix(60, 'Z');
 * R60.resize(3, 3); // Get rid of homogeneous part
 * Matrix2D< RFLOAT > I(3, 3);
 * I.initIdentity();
 * Euler_apply_transf(I, R60, rot, tilt, psi, newrot, newtilt, newpsi);
 * @endcode
 */
void Euler_apply_transf(const Matrix2D< RFLOAT >& L,
                        const Matrix2D< RFLOAT >& R,
                        RFLOAT rot,
                        RFLOAT tilt,
                        RFLOAT psi,
                        RFLOAT& newrot,
                        RFLOAT& newtilt,
                        RFLOAT& newpsi);

/** 3D Rotation matrix after 3 Euler angles
 *
 * Creates a rotational matrix (4x4) for volumes around the combination of the 3
 * rotations around ZYZ. All angles are in degrees. You must use it with
 * IS_NOT_INV in applyGeometry.
 *
 * @code
 * Matrix2D< float > euler = Euler_rotation3DMatrix(60, 30, 60);
 * @endcode
 */
void Euler_rotation3DMatrix(RFLOAT rot, RFLOAT tilt, RFLOAT psi,
                            Matrix2D<RFLOAT> &result);

//@}

#endif

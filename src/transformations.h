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

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include "src/multidim_array.h"
#include "src/euler.h"


#define IS_INV true
#define IS_NOT_INV false
#define DONT_WRAP false
#define WRAP true
#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-50
#endif

/// @defgroup GeometricalTransformations Geometrical transformations
/// @ingroup DataLibrary
//@{
/** Creates a rotational matrix (3x3) for images
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees.
 * m must have been already resized to 3x3
 *
 * @code
 *  rotation2DMatrix(60,m);
 * @endcode
 */
void rotation2DMatrix(RFLOAT ang, Matrix2D< RFLOAT > &m, bool homogeneous=true);

/** Creates a translational matrix (3x3) for images
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R2 vector (shift_X, shift_Y). An exception is thrown
 * if the displacement is not a R2 vector.
 *
 * @code
 * // Displacement of 1 pixel to the right
 * m = translation2DMatrix(vectorR2(1, 0));
 * @endcode
 */
void translation2DMatrix(const Matrix1D< RFLOAT > &v, Matrix2D< RFLOAT > &m);

/** Creates a rotational matrix (4x4) for volumes around system axis
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is either 'X', 'Y'
 * or 'Z'. An exception is thrown if the axis given is not one of these.
 *
 * The returned matrices are respectively alpha degrees around Z
 *
 * @code
 * [ cos(A) -sin(A)     0   ]
 * [ sin(A)  cos(A)     0   ]
 * [   0       0        1   ]
 * @endcode
 *
 * alpha degrees around Y
 * @code
 * [ cos(A)    0    -sin(A) ]
 * [   0       1       0    ]
 * [ sin(A)    0     cos(A) ]
 * @endcode
 *
 * alpha degrees around X
 * @code
 * [   1       0       0    ]
 * [   0     cos(A) -sin(A) ]
 * [   0     sin(A)  cos(A) ]
 * @endcode
 *
 * @code
 * m = rotation3DMatrix(60, 'X');
 * @endcode
 */
void rotation3DMatrix(RFLOAT ang, char axis, Matrix2D< RFLOAT > &m,
		bool homogeneous=true);

/** Creates a rotational matrix (4x4) for volumes around any axis
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is given as a R3
 * vector. An exception is thrown if the axis is not a R3 vector. The axis needs
 * not to be unitary.
 *
 * @code
 * m = rotation3DMatrix(60, vectorR3(1, 1, 1));
 * @endcode
 */
void rotation3DMatrix(RFLOAT ang, const Matrix1D< RFLOAT >& axis, Matrix2D< RFLOAT > &m,
		bool homogeneous=true);

/** Matrix which transforms the given axis into Z
 * @ingroup GeometricalTransformations
 *
 * A geometrical transformation matrix (4x4) is returned such that the given
 * axis is rotated until it is aligned with the Z axis. This is very useful in
 * order to produce rotational matrices, for instance, around any axis.
 *
 * @code
 * Matrix2D< RFLOAT > A = alignWithZ(axis);
 * return A.transpose() * rotation3DMatrix(ang, 'Z') * A;
 * @endcode
 *
 * The returned matrix is such that A*axis=Z, where Z and axis are column
 * vectors.
 */
void alignWithZ(const Matrix1D< RFLOAT >& axis, Matrix2D< RFLOAT > &m, bool homogeneous=true);

/** Creates a translational matrix (4x4) for volumes
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R3 vector (shift_X, shift_Y, shift_Z). An exception
 * is thrown if the displacement is not a R3 vector.
 *
 * @code
 * // Displacement of 2 pixels down
 * m = translation3DMatrix(vectorR3(0, 0, 2));
 * @endcode
 */
void translation3DMatrix(const Matrix1D< RFLOAT >& v, Matrix2D< RFLOAT > &m);

/** Creates a scaling matrix (4x4) for volumes
 * @ingroup GeometricalTransformations
 *
 * The scaling factors for the different axis must be given as a vector. So
 * that, XX(sc)=scale for X axis, YY(sc)=...
 */
void scale3DMatrix(const Matrix1D< RFLOAT >& sc, Matrix2D< RFLOAT > &m,
		bool homogeneous=true);

/** Applies a geometrical transformation.
 * @ingroup GeometricalTransformations
 *
 * Any geometrical transformation defined by the matrix A (RFLOAT (4x4)!!
 * ie, in homogeneous R3 coordinates) is applied to the volume V1.
 * The result is stored in V2 (it cannot be the same as the input volume).
 * An exception is thrown if the transformation matrix is not 4x4.
 *
 * Structure of the transformation matrix: It should have the following
 * components
 *
 * r11 r12 r13 x
 * r21 r22 r23 y
 * r31 r32 r33 z
 * 0   0   0   1
 *
 * where (x,y,z) is the translation desired, and Rij are the components of
 * the rotation matrix R. If you want to apply a scaling factor to the
 * transformation, then multiply r11, r22 and r33 by it.
 *
 * The result volume (with ndim=1) is resized to the same
 * dimensions as V1 if V2 is empty (0x0) at the beginning, if it
 * is not, ie, if V2 has got some size then only those values in
 * the volume are filled, this is very useful for resizing the
 * volume, then you manually resize the output volume to the
 * desired size and then call this routine.
 *
 * The relationship between the output coordinates and the input ones are
 *
 * @code
 * out = A * in
 * (x, y, z) = A * (x', y', z')
 * @endcode
 *
 * This function works independently from the logical indexing of each
 * matrix, it sets the logical center and the physical center of the image
 * and work with these 2 coordinate spaces. At the end the original logical
 * indexing of each matrix is kept.
 *
 * The procedure followed goes from coordinates in the output volume
 * to the ones in the input one, so the inverse of the A matrix is
 * needed. There is a flag telling if the given matrix is already
 * the inverse one or the normal one. If it is the normal one internally
 * the matrix is inversed. If you are to do many "rotations" then
 * some time is spent in inverting the matrix. Normally the matrix is the
 * normal one.
 *
 * There is something else to tell about the geometrical tranformation.
 * The value of the voxel in the output volume is computed via
 * bilinear interpolation in the input volume. If any of the voxels
 * participating in the interpolation falls outside the input volume,
 * then automatically the corresponding output voxel is set to 0, unless
 * that the wrap flag has been set to 1. In this case if the voxel
 * falls out by the right hand then it is "wrapped" and the corresponding
 * voxel in the left hand is used. The same is appliable to top-bottom.
 * Usually wrap mode is off. Wrap mode is interesting for translations
 * but not for rotations, for example.
 *
 * The inverse mode and wrapping mode should be taken by default by the
 * routine, g++ seems to have problems with template functions outside
 * a class with default parameters. So, I'm sorry, you will have to
 * put them always. The usual combination is
 *
 * applyGeometry(..., IS_NOT_INV, DONT_WRAP).
 *
 * Although you can also use the constants IS_INV, or WRAP.
 *
 * @code
 * Matrix2D< RFLOAT > A(4,4);
 * A.initIdentity;
 * applyGeometry(V2, A, V1);
 * @endcode
 */
template<typename T>
void applyGeometry(const MultidimArray<T>& V1,
                   MultidimArray<T>& V2,
                   const Matrix2D< RFLOAT > A,
                   bool inv,
                   bool wrap,
                   T outside = 0)
{

    if (&V1 == &V2)
        REPORT_ERROR("ApplyGeometry: Input array cannot be the same as output array");

    if ( V1.getDim()==2 && ((MAT_XSIZE(A) != 3) || (MAT_YSIZE(A) != 3)) )
        REPORT_ERROR("ApplyGeometry: 2D transformation matrix is not 3x3");

    if ( V1.getDim()==3 && ((MAT_XSIZE(A) != 4) || (MAT_YSIZE(A) != 4)) )
        REPORT_ERROR("ApplyGeometry: 3D transformation matrix is not 4x4");

    if (A.isIdentity())
    {
        V2=V1;
        return;
    }

    if (XSIZE(V1) == 0)
    {
        V2.clear();
        return;
    }

    Matrix2D<RFLOAT> Ainv;
    const Matrix2D<RFLOAT> * Aptr=&A;
    if (!inv)
    {
        Ainv = A.inv();
        Aptr=&Ainv;
    }
    const Matrix2D<RFLOAT> &Aref=*Aptr;

    // For scalings the output matrix is resized outside to the final
    // size instead of being resized inside the routine with the
    // same size as the input matrix
    if (XSIZE(V2) == 0)
        V2.resize(V1);

    if (V1.getDim() == 2)
    {
        // 2D transformation

        int m1, n1, m2, n2;
        RFLOAT x, y, xp, yp;
        RFLOAT minxp, minyp, maxxp, maxyp;
        int cen_x, cen_y, cen_xp, cen_yp;
        RFLOAT wx, wy;
        int Xdim, Ydim;

        // Find center and limits of image
        cen_y  = (int)(YSIZE(V2) / 2);
        cen_x  = (int)(XSIZE(V2) / 2);
        cen_yp = (int)(YSIZE(V1) / 2);
        cen_xp = (int)(XSIZE(V1) / 2);
        minxp  = -cen_xp;
        minyp  = -cen_yp;
        maxxp  = XSIZE(V1) - cen_xp - 1;
        maxyp  = YSIZE(V1) - cen_yp - 1;
        Xdim   = XSIZE(V1);
        Ydim   = YSIZE(V1);

        // Now we go from the output image to the input image, ie, for any pixel
        // in the output image we calculate which are the corresponding ones in
        // the original image, make an interpolation with them and put this value
        // at the output pixel

#ifdef DEBUG_APPLYGEO
        std::cout << "A\n" << Aref << std::endl
        << "(cen_x ,cen_y )=(" << cen_x  << "," << cen_y  << ")\n"
        << "(cen_xp,cen_yp)=(" << cen_xp << "," << cen_yp << ")\n"
        << "(min_xp,min_yp)=(" << minxp  << "," << minyp  << ")\n"
        << "(max_xp,max_yp)=(" << maxxp  << "," << maxyp  << ")\n";
#endif

        for (int i = 0; i < YSIZE(V2); i++)
        {
            // Calculate position of the beginning of the row in the output image
            x = -cen_x;
            y = i - cen_y;

            // Calculate this position in the input image according to the
            // geometrical transformation
            // they are related by
            // coords_output(=x,y) = A * coords_input (=xp,yp)
            xp = x * Aref(0, 0) + y * Aref(0, 1) + Aref(0, 2);
            yp = x * Aref(1, 0) + y * Aref(1, 1) + Aref(1, 2);

            for (int j = 0; j < XSIZE(V2); j++)
            {
                bool interp;
                T tmp;

#ifdef DEBUG_APPLYGEO

                std::cout << "Computing (" << i << "," << j << ")\n";
                std::cout << "   (y, x) =(" << y << "," << x << ")\n"
                << "   before wrapping (y',x')=(" << yp << "," << xp << ") "
                << std::endl;
#endif
                // If the point is outside the image, apply a periodic extension
                // of the image, what exits by one side enters by the other
                interp = true;
                if (wrap)
                {
                    if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                        xp > maxxp + XMIPP_EQUAL_ACCURACY)
                        xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);

                    if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                        yp > maxyp + XMIPP_EQUAL_ACCURACY)

                        yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);
                }
                else
                {
                    if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                        xp > maxxp + XMIPP_EQUAL_ACCURACY)
                        interp = false;

                    if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                        yp > maxyp + XMIPP_EQUAL_ACCURACY)
                        interp = false;
                }

#ifdef DEBUG_APPLYGEO
                std::cout << "   after wrapping (y',x')=(" << yp << "," << xp << ") "
                << std::endl;
                std::cout << "   Interp = " << interp << std::endl;
                // The following line sounds dangerous...
                //x++;
#endif

                if (interp)
                {
                        // Linear interpolation

                        // Calculate the integer position in input image, be careful
                        // that it is not the nearest but the one at the top left corner
                        // of the interpolation square. Ie, (0.7,0.7) would give (0,0)
                        // Calculate also weights for point m1+1,n1+1
                        wx = xp + cen_xp;
                        m1 = (int) wx;
                        wx = wx - m1;
                        m2 = m1 + 1;
                        wy = yp + cen_yp;
                        n1 = (int) wy;
                        wy = wy - n1;
                        n2 = n1 + 1;

                        // m2 and n2 can be out by 1 so wrap must be check here
                        if (wrap)
                        {
                            if (m2 >= Xdim)
                                m2 = 0;
                            if (n2 >= Ydim)
                                n2 = 0;
                        }

#ifdef DEBUG_APPLYGEO
                        std::cout << "   From (" << n1 << "," << m1 << ") and ("
                        << n2 << "," << m2 << ")\n";
                        std::cout << "   wx= " << wx << " wy= " << wy << std::endl;
#endif

                        // Perform interpolation
                        // if wx == 0 means that the rightest point is useless for this
                        // interpolation, and even it might not be defined if m1=xdim-1
                        // The same can be said for wy.
                        tmp  = (T)((1 - wy) * (1 - wx) * DIRECT_A2D_ELEM(V1, n1, m1));

                        if (m2 < V1.xdim)
                            tmp += (T)((1 - wy) * wx * DIRECT_A2D_ELEM(V1, n1, m2));

                        if (n2 < V1.ydim)
                        {
                            tmp += (T)(wy * (1 - wx) * DIRECT_A2D_ELEM(V1, n2, m1));

                            if (m2 < V1.xdim)
                                tmp += (T)(wy * wx * DIRECT_A2D_ELEM(V1, n2, m2));
                        }

                        dAij(V2, i, j) = tmp;
#ifdef DEBUG_APPYGEO
                    std::cout << "   val= " << dAij(V2, i, j) << std::endl;
#endif

                } // if interp
                else
                	dAij(V2, i, j) = outside;

                // Compute new point inside input image
                xp += Aref(0, 0);
                yp += Aref(1, 0);
            }
        }
    }
    else
    {
        // 3D transformation

        int m1, n1, o1, m2, n2, o2;
        RFLOAT x, y, z, xp, yp, zp;
        RFLOAT minxp, minyp, maxxp, maxyp, minzp, maxzp;
        int cen_x, cen_y, cen_z, cen_xp, cen_yp, cen_zp;
        RFLOAT wx, wy, wz;

        // Find center of MultidimArray
        cen_z = (int)(V2.zdim / 2);
        cen_y = (int)(V2.ydim / 2);
        cen_x = (int)(V2.xdim / 2);
        cen_zp = (int)(V1.zdim / 2);
        cen_yp = (int)(V1.ydim / 2);
        cen_xp = (int)(V1.xdim / 2);
        minxp = -cen_xp;
        minyp = -cen_yp;
        minzp = -cen_zp;
        maxxp = V1.xdim - cen_xp - 1;
        maxyp = V1.ydim - cen_yp - 1;
        maxzp = V1.zdim - cen_zp - 1;

#ifdef DEBUG

        std::cout << "Geometry 2 center=("
        << cen_z  << "," << cen_y  << "," << cen_x  << ")\n"
        << "Geometry 1 center=("
        << cen_zp << "," << cen_yp << "," << cen_xp << ")\n"
        << "           min=("
        << minzp  << "," << minyp  << "," << minxp  << ")\n"
        << "           max=("
        << maxzp  << "," << maxyp  << "," << maxxp  << ")\n"
        ;
#endif

        // Now we go from the output MultidimArray to the input MultidimArray, ie, for any
        // voxel in the output MultidimArray we calculate which are the corresponding
        // ones in the original MultidimArray, make an interpolation with them and put
        // this value at the output voxel

        // V2 is not initialised to 0 because all its pixels are rewritten
        for (int k = 0; k < V2.zdim; k++)
            for (int i = 0; i < V2.ydim; i++)
            {
                // Calculate position of the beginning of the row in the output
                // MultidimArray
                x = -cen_x;
                y = i - cen_y;
                z = k - cen_z;

                // Calculate this position in the input image according to the
                // geometrical transformation they are related by
                // coords_output(=x,y) = A * coords_input (=xp,yp)
                xp = x * Aref(0, 0) + y * Aref(0, 1) + z * Aref(0, 2) + Aref(0, 3);
                yp = x * Aref(1, 0) + y * Aref(1, 1) + z * Aref(1, 2) + Aref(1, 3);
                zp = x * Aref(2, 0) + y * Aref(2, 1) + z * Aref(2, 2) + Aref(2, 3);

                for (int j = 0; j < V2.xdim; j++)
                {
                    bool interp;
                    T tmp;

#ifdef DEBUG

                    bool show_debug = false;
                    if ((i == 0 && j == 0 && k == 0) ||
                        (i == V2.ydim - 1 && j == V2.xdim - 1 && k == V2.zdim - 1))
                        show_debug = true;

                    if (show_debug)
                        std::cout << "(x,y,z)-->(xp,yp,zp)= "
                        << "(" << x  << "," << y  << "," << z  << ") "
                        << "(" << xp << "," << yp << "," << zp << ")\n";
#endif

                    // If the point is outside the volume, apply a periodic
                    // extension of the volume, what exits by one side enters by
                    // the other
                    interp  = true;
                    if (wrap)
                    {
                        if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                            xp > maxxp + XMIPP_EQUAL_ACCURACY)
                            xp = realWRAP(xp, minxp - 0.5, maxxp + 0.5);

                        if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                            yp > maxyp + XMIPP_EQUAL_ACCURACY)
                            yp = realWRAP(yp, minyp - 0.5, maxyp + 0.5);

                        if (zp < minzp - XMIPP_EQUAL_ACCURACY ||
                            zp > maxzp + XMIPP_EQUAL_ACCURACY)
                            zp = realWRAP(zp, minzp - 0.5, maxzp + 0.5);
                    }
                    else
                    {
                        if (xp < minxp - XMIPP_EQUAL_ACCURACY ||
                            xp > maxxp + XMIPP_EQUAL_ACCURACY)
                            interp = false;

                        if (yp < minyp - XMIPP_EQUAL_ACCURACY ||
                            yp > maxyp + XMIPP_EQUAL_ACCURACY)
                            interp = false;

                        if (zp < minzp - XMIPP_EQUAL_ACCURACY ||
                            zp > maxzp + XMIPP_EQUAL_ACCURACY)
                            interp = false;
                    }

                    if (interp)
                    {

                            // Linear interpolation

                            // Calculate the integer position in input volume, be
                            // careful that it is not the nearest but the one at the
                            // top left corner of the interpolation square. Ie,
                            // (0.7,0.7) would give (0,0)
                            // Calculate also weights for point m1+1,n1+1
                            wx = xp + cen_xp;
                            m1 = (int) wx;
                            wx = wx - m1;
                            m2 = m1 + 1;
                            wy = yp + cen_yp;
                            n1 = (int) wy;
                            wy = wy - n1;
                            n2 = n1 + 1;
                            wz = zp + cen_zp;
                            o1 = (int) wz;
                            wz = wz - o1;
                            o2 = o1 + 1;

#ifdef DEBUG

                            if (show_debug)
                            {
                                std::cout << "After wrapping(xp,yp,zp)= "
                                << "(" << xp << "," << yp << "," << zp << ")\n";
                                std::cout << "(m1,n1,o1)-->(m2,n2,o2)="
                                << "(" << m1 << "," << n1 << "," << o1 << ") "
                                << "(" << m2 << "," << n2 << "," << o2 << ")\n";
                                std::cout << "(wx,wy,wz)="
                                << "(" << wx << "," << wy << "," << wz << ")\n";
                            }
#endif

                            // Perform interpolation
                            // if wx == 0 means that the rightest point is useless for
                            // this interpolation, and even it might not be defined if
                            // m1=xdim-1
                            // The same can be said for wy.
                            tmp  = (T)((1 - wz) * (1 - wy) * (1 - wx) * DIRECT_A3D_ELEM(V1, o1, n1, m1));

                            if (m2 < V1.xdim)
                                tmp += (T)((1 - wz) * (1 - wy) * wx * DIRECT_A3D_ELEM(V1, o1, n1, m2));

                            if (n2 < V1.ydim)
                            {
                                tmp += (T)((1 - wz) * wy * (1 - wx) * DIRECT_A3D_ELEM(V1, o1, n2, m1));
                                if (m2 < V1.xdim)
                                    tmp += (T)((1 - wz) * wy * wx * DIRECT_A3D_ELEM(V1, o1, n2, m2));
                            }

                            if (o2 < V1.zdim)
                            {
                                tmp += (T)(wz * (1 - wy) * (1 - wx) * DIRECT_A3D_ELEM(V1, o2, n1, m1));
                                if (m2 < V1.xdim)
                                    tmp += (T)(wz * (1 - wy) * wx * DIRECT_A3D_ELEM(V1, o2, n1, m2));
                                if (n2 < V1.ydim)
                                {
                                    tmp += (T)(wz * wy * (1 - wx) * DIRECT_A3D_ELEM(V1, o2, n2, m1));
                                    if (m2 < V1.xdim)
                                        tmp += (T)(wz * wy * wx * DIRECT_A3D_ELEM(V1, o2, n2, m2));
                                }
                            }

#ifdef DEBUG
                            if (show_debug)
                                std::cout <<
                                "tmp1=" << DIRECT_A3D_ELEM(V1, o1, n1, m1) << " "
                                << (T)((1 - wz) *(1 - wy) *(1 - wx) * DIRECT_A3D_ELEM(V1, o1, n1, m1))
                                << std::endl <<
                                "tmp2=" << DIRECT_A3D_ELEM(V1, o1, n1, m2) << " "
                                << (T)((1 - wz) *(1 - wy) * wx * DIRECT_A3D_ELEM(V1, o1, n1, m2))
                                << std::endl <<
                                "tmp3=" << DIRECT_A3D_ELEM(V1, o1, n2, m1) << " "
                                << (T)((1 - wz) * wy *(1 - wx) * DIRECT_A3D_ELEM(V1, o1, n2, m1))
                                << std::endl <<
                                "tmp4=" << DIRECT_A3D_ELEM(V1, o1, n2, m2) << " "
                                << (T)((1 - wz) * wy * wx * DIRECT_A3D_ELEM(V1, o2, n1, m1))
                                << std::endl <<
                                "tmp6=" << DIRECT_A3D_ELEM(V1, o2, n1, m2) << " "
                                << (T)(wz * (1 - wy) * wx * DIRECT_A3D_ELEM(V1, o2, n1, m2))
                                << std::endl <<
                                "tmp7=" << DIRECT_A3D_ELEM(V1, o2, n2, m1) << " "
                                << (T)(wz * wy *(1 - wx) * DIRECT_A3D_ELEM(V1, o2, n2, m1))
                                << std::endl <<
                                "tmp8=" << DIRECT_A3D_ELEM(V1, o2, n2, m2) << " "
                                << (T)(wz * wy * wx * DIRECT_A3D_ELEM(V1, o2, n2, m2))
                                << std::endl <<
                                "tmp= " << tmp << std::endl;
#endif

                            dAkij(V2 , k, i, j) = tmp;
                    }
                    else
                        dAkij(V2, k, i, j) = outside;

                    // Compute new point inside input image
                    xp += Aref(0, 0);
                    yp += Aref(1, 0);
                    zp += Aref(2, 0);
                }
            }
    }
}

/** Applies a geometrical transformation and overwrites the input matrix.
 * @ingroup GeometricalTransformations
 *
 * The same as the previous function, but input array is overwritten
 */
template<typename T>
void selfApplyGeometry(MultidimArray<T>& V1,
                       const Matrix2D< RFLOAT > A, bool inv,
                       bool wrap, T outside = 0)
{
    MultidimArray<T> aux = V1;
    V1.initZeros();
    applyGeometry(aux, V1, A, inv, wrap, outside);
}

/** Rotate an array around a given system axis.
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is either
 * 'X', 'Y' or 'Z' for 3D arrays and only 'Z' for 2D arrays. An
 * exception is thrown if the axis given is not one of these.
 *
 * @code
 * rotate(Vin, Vout, 60);
 * @endcode
 */
template<typename T>
void rotate(const MultidimArray<T>& V1,
            MultidimArray<T>& V2,
            RFLOAT ang, char axis = 'Z',
            bool wrap = DONT_WRAP, T outside = 0)
{
    Matrix2D< RFLOAT > tmp;
    if (V1.getDim()==2)
    {
        rotation2DMatrix(ang,tmp);
    }
    else if (V1.getDim()==3)
    {
        rotation3DMatrix(ang, axis, tmp);
    }
    else
        REPORT_ERROR("rotate ERROR: rotate only valid for 2D or 3D arrays");

    applyGeometry(V1, V2, tmp, IS_NOT_INV, wrap, outside);
}

/** Rotate an array around a given system axis.
 * @ingroup GeometricalTransformations
 *
 * The same as the previous function, but input array is overwritten
 */
template<typename T>
void selfRotate(MultidimArray<T>& V1,
                RFLOAT ang, char axis = 'Z',
                bool wrap = DONT_WRAP, T outside = 0)
{
    MultidimArray<T> aux = V1;
    rotate(aux, V1, ang, axis, wrap, outside);
}

/** Translate a array.
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R2 or R3 vector (shift_X, shift_Y, shift_Z) for 2D and 3D arrays, respectively.
 * An exception is thrown if the displacement is not a R3 vector.
 *
 * @code
 * // Displacement of 2 pixels down
 * V2 = V1.translate(vectorR3(0, 0, 2));
 * @endcode
 */
template<typename T>
void translate(const MultidimArray<T> &V1,
               MultidimArray<T> &V2,
               const Matrix1D< RFLOAT >& v,
               bool wrap = WRAP, T outside = 0)
{
    Matrix2D< RFLOAT > tmp;
    if (V1.getDim()==2)
        translation2DMatrix(v, tmp);
    else if (V1.getDim()==3)
        translation3DMatrix(v, tmp);
    else
        REPORT_ERROR("translate ERROR: translate only valid for 2D or 3D arrays");

    applyGeometry(V1, V2, tmp, IS_NOT_INV, wrap, outside);
}

/** Translate an array.
 * @ingroup GeometricalTransformations
 *
 * The same as the previous function, but input array is overwritten
 */
template<typename T>
void selfTranslate(MultidimArray<T>& V1,
                   const Matrix1D< RFLOAT >& v,
                   bool wrap = WRAP, T outside = 0)
{
    MultidimArray<T> aux = V1;
    translate(aux, V1, v, wrap, outside);
}

/** Translate center of mass to center
 * @ingroup GeometricalTransformations
 *
 * If the input has very high values, it is better to rescale it to be
 * between 0 and 1.
 */
template<typename T>
void translateCenterOfMassToCenter(const MultidimArray<T> &V1,
                                   MultidimArray<T> &V2,
                                   bool wrap = WRAP,
								   bool verb = false)
{
    V2 = V1;
    V2.setXmippOrigin();
    Matrix1D< RFLOAT > center;
    V2.centerOfMass(center);
    if (verb)
    {
    	std::cout << " Center of mass: x= " << XX(center) << " y= " << YY(center) << " z= " << ZZ(center) << std::endl;
    }
    center *= -1;
    translate(V1, V2, center, wrap, (RFLOAT)0.);
}

/** Translate center of mass to center
 * @ingroup GeometricalTransformations
 *
 * The same as the previous function, but input array is overwritten
 */
template<typename T>
void selfTranslateCenterOfMassToCenter(MultidimArray<T> &V1,
                                       bool wrap = WRAP,
									   bool verb = false)
{
    MultidimArray<T> aux = V1;
    translateCenterOfMassToCenter(aux, V1, wrap, verb);
}

/** Scales to a new size.
 * @ingroup GeometricalTransformations
 *
 * The volume is scaled (resampled) to fill a new size. It is not the
 * same as "window" in this same class. The size can be larger or smaller
 * than the actual one.
 *
 * @code
 * scaleToSize(Vin, Vout, 128, 128, 128);
 * @endcode
 */
template<typename T>
void scaleToSize(const MultidimArray<T> &V1,
                 MultidimArray<T> &V2,
                 int Xdim, int Ydim, int Zdim = 1)
{

    Matrix2D< RFLOAT > tmp;
    if (V1.getDim()==2)
    {
        tmp.initIdentity(3);
        tmp(0, 0) = (RFLOAT) Xdim / (RFLOAT) XSIZE(V1);
        tmp(1, 1) = (RFLOAT) Ydim / (RFLOAT) YSIZE(V1);
        V2.resize(1, 1, Ydim, Xdim);
    }
    else if (V1.getDim()==3)
    {
        tmp.initIdentity(4);
        tmp(0, 0) = (RFLOAT) Xdim / (RFLOAT) XSIZE(V1);
        tmp(1, 1) = (RFLOAT) Ydim / (RFLOAT) YSIZE(V1);
        tmp(2, 2) = (RFLOAT) Zdim / (RFLOAT) ZSIZE(V1);
        V2.resize(1, Zdim, Ydim, Xdim);
    }
    else
        REPORT_ERROR("scaleToSize ERROR: scaleToSize only valid for 2D or 3D arrays");

    applyGeometry(V1, V2, tmp, IS_NOT_INV, WRAP, (T)0);
}

/** Scales to a new size.
 * @ingroup GeometricalTransformations
 *
 * The same as the previous function, but input array is overwritten
 */
template<typename T>
void selfScaleToSize(MultidimArray<T> &V1,
                     int Xdim, int Ydim, int Zdim = 1)
{
    MultidimArray<T> aux = V1;
    scaleToSize(aux, V1, Xdim, Ydim, Zdim);
}

/** Does a radial average of a 2D/3D image, around the voxel where is the origin.
 * @ingroup GeometricalTransformations
 *
 * A vector radial_mean is returned where:
 * - the first element is the mean of the voxels whose
 *   distance to the origin is (0-1),
 * - the second element is the mean of the voxels
 *   whose distance to the origin is (1-2)
 * - and so on.
 *
 * A second vector radial_count is returned containing the number of voxels
 * over which each radial average was calculated.
 *
 * Sjors nov2003: if rounding=true, element=round(distance);
 * - so the first element is the mean of the voxels whose distance to the
 *   origin is (0.5-1.5),
 * - the second element is the mean of the voxels whose distance to the origin
 *   is (1.5-2.5)
 * - and so on.
 */
template<typename T>
void radialAverage(const MultidimArray< T >& m,
                   Matrix1D< int >& center_of_rot,
                   MultidimArray< T >& radial_mean,
                   MultidimArray< int >& radial_count,
                   const bool& rounding = false)
{
    Matrix1D< RFLOAT > idx(3);

    // If center_of_rot was written for 2D image
    if (center_of_rot.size() < 3)
        center_of_rot.resize(3);

    // First determine the maximum distance that one should expect, to set the
    // dimension of the radial average vector
    MultidimArray< int > distances(8);

    RFLOAT z = STARTINGZ(m) - ZZ(center_of_rot);
    RFLOAT y = STARTINGY(m) - YY(center_of_rot);
    RFLOAT x = STARTINGX(m) - XX(center_of_rot);

    distances(0) = (int) floor(sqrt(x * x + y * y + z * z));
    x = FINISHINGX(m) - XX(center_of_rot);

    distances(1) = (int) floor(sqrt(x * x + y * y + z * z));
    y = FINISHINGY(m) - YY(center_of_rot);

    distances(2) = (int) floor(sqrt(x * x + y * y + z * z));
    x = STARTINGX(m) - XX(center_of_rot);

    distances(3) = (int) floor(sqrt(x * x + y * y + z * z));
    z = FINISHINGZ(m) - ZZ(center_of_rot);

    distances(4) = (int) floor(sqrt(x * x + y * y + z * z));
    x = FINISHINGX(m) - XX(center_of_rot);

    distances(5) = (int) floor(sqrt(x * x + y * y + z * z));
    y = STARTINGY(m) - YY(center_of_rot);

    distances(6) = (int) floor(sqrt(x * x + y * y + z * z));
    x = STARTINGX(m) - XX(center_of_rot);

    distances(7) = (int) floor(sqrt(x * x + y * y + z * z));

    int dim = (int) CEIL(distances.computeMax()) + 1;
    if (rounding)
        dim++;

    // Define the vectors
    radial_mean.resize(dim);
    radial_mean.initZeros();
    radial_count.resize(dim);
    radial_count.initZeros();

    // Perform the radial sum and count pixels that contribute to every
    // distance
    FOR_ALL_ELEMENTS_IN_ARRAY3D(m)
    {
        ZZ(idx) = k - ZZ(center_of_rot);
        YY(idx) = i - YY(center_of_rot);
        XX(idx) = j - XX(center_of_rot);

        // Determine distance to the center
        int distance;
        if (rounding)
            distance = (int) ROUND(idx.module());
        else
            distance = (int) floor(idx.module());

        // Sum te value to the pixels with the same distance
        radial_mean(distance) += A3D_ELEM(m, k, i, j);

        // Count the pixel
        radial_count(distance)++;
    }

    // Perform the mean
    FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_mean)
    radial_mean(i) /= (T) radial_count(i);
}

//@}
#endif

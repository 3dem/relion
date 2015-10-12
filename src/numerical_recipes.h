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
#ifndef _NUMERICAL_HH
#   define _NUMERICAL_HH

#include <math.h>
#include "src/memory.h"
#include "src/macros.h"
#include "src/error.h"

/*****************************************************************************/
/* Variable and prototype definitions for the Numerical Core                 */
/*****************************************************************************/

//@defgroup NumericalRecipes Functions from the Numerical Recipes
//@ingroup DataLibrary
//@{

// Utilities --------------------------------------------------------------
void nrerror(const char error_text[]);

// Bessel functions --------------------------------------------------------
RFLOAT bessj0(RFLOAT x);
RFLOAT bessj3_5(RFLOAT x);
RFLOAT bessj1_5(RFLOAT x);

RFLOAT bessi0(RFLOAT x);
RFLOAT bessi1(RFLOAT x);
RFLOAT bessi0_5(RFLOAT x);
RFLOAT bessi1_5(RFLOAT x);
RFLOAT bessi2(RFLOAT x);
RFLOAT bessi3(RFLOAT x);
RFLOAT bessi2_5(RFLOAT x);
RFLOAT bessi3_5(RFLOAT x);
RFLOAT bessi4(RFLOAT x);

// Special functions -------------------------------------------------------
RFLOAT gammln(RFLOAT xx);
RFLOAT gammp(RFLOAT a, RFLOAT x);
RFLOAT betacf(RFLOAT a, RFLOAT b, RFLOAT x);
RFLOAT betai(RFLOAT a, RFLOAT b, RFLOAT x);

// Singular value descomposition of matrix a (numerical recipes, chapter 2-6 for details)
void svdcmp(RFLOAT *a, int m, int n, RFLOAT *w, RFLOAT *v);
void svbksb(RFLOAT *u, RFLOAT *w, RFLOAT *v, int m, int n, RFLOAT *b, RFLOAT *x);

// Optimization ------------------------------------------------------------
void powell(RFLOAT *p, RFLOAT *xi, int n, RFLOAT ftol, int &iter,
            RFLOAT &fret, RFLOAT(*func)(RFLOAT *, void *), void *prm,
            bool show);

// Working with matrices ---------------------------------------------------
// LU decomposition
#define TINY 1.0e-20;
/* Chapter 2 Section 3: LU DECOMPOSITION */
template <class T>
void ludcmp(T *a, int n, int *indx, T *d)
{
    int i, imax, j, k;
    T big, dum, sum, temp;
    T *vv;

    ask_Tvector(vv, 1, n);
    *d = (T)1.0;
    for (i = 1;i <= n;i++)
    {
        big = (T)0.0;
        for (j = 1;j <= n;j++)
            if ((temp = (T)fabs((RFLOAT)a[i*n+j])) > big)
                big = temp;
        if (big == (T)0.0)
            nrerror("Singular matrix in routine LUDCMP");
        vv[i] = (T)1.0 / big;
    }
    for (j = 1;j <= n;j++)
    {
        for (i = 1;i < j;i++)
        {
            sum = a[i*n+j];
            for (k = 1;k < i;k++)
                sum -= a[i*n+k] * a[k*n+j];
            a[i*n+j] = sum;
        }
        big = (T)0.0;
        for (i = j;i <= n;i++)
        {
            sum = a[i*n+j];
            for (k = 1;k < j;k++)
                sum -= a[i*n+k] * a[k*n+j];
            a[i*n+j] = sum;
            if ((dum = vv[i] * (T)fabs((RFLOAT)sum)) >= big)
            {
                big = dum;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 1;k <= n;k++)
            {
                dum = a[imax*n+k];
                a[imax*n+k] = a[j*n+k];
                a[j*n+k] = dum;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[j*n+j] == 0.0)
            a[j*n+j] = (T) TINY;
        if (j != n)
        {
            dum = (T)1.0 / (a[j*n+j]);
            for (i = j + 1;i <= n;i++)
                a[i*n+j] *= dum;
        }
    }
    free_Tvector(vv, 1, n);
}
#undef TINY

// Solve Ax=b
/* Chapter 2 Section 3: LU BACKWARD-FORWARD SUBSTITUTION */
template <class T>
void lubksb(T *a, int n, int *indx, T b[])
{
    int i, ii = 0, ip, j;
    T sum;

    for (i = 1;i <= n;i++)
    {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii)
            for (j = ii;j <= i - 1;j++)
                sum -= a[i*n+j] * b[j];
        else if (sum)
            ii = i;
        b[i] = sum;
    }
    for (i = n;i >= 1;i--)
    {
        sum = b[i];
        for (j = i + 1;j <= n;j++)
            sum -= a[i*n+j] * b[j];
        b[i] = sum / a[i*n+i];
    }
}

/* Chapter 2, Section 1. Gauss-Jordan equation system resolution ----------- */
// Solve Ax=b (b=matrix)
template <class T>
void gaussj(T *a, int n, T *b, int m)
{
    T temp;
    int *indxc, *indxr, *ipiv;
    int i, icol, irow, j, k, l, ll;
    T big, dum;
    RFLOAT pivinv;

    ask_Tvector(indxc, 1, n);
    ask_Tvector(indxr, 1, n);
    ask_Tvector(ipiv, 1, n);
    for (j = 1;j <= n;j++)
        ipiv[j] = 0;
    for (i = 1;i <= n;i++)
    {
        big = (T)0;
        for (j = 1;j <= n;j++)
            if (ipiv[j] != 1)
                for (k = 1;k <= n;k++)
                {
                    if (ipiv[k] == 0)
                    {
                        if (fabs((RFLOAT)a[j*n+k]) >= (RFLOAT) big)
                        {
                            big = ABS(a[j*n+k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                        nrerror("GAUSSJ: Singular Matrix-1");
                }
        ++(ipiv[icol]);
        if (irow != icol)
        {
            for (l = 1;l <= n;l++)
                SWAP(a[irow*n+l], a[icol*n+l], temp)
                for (l = 1;l <= m;l++)
                    SWAP(b[irow*n+l], b[icol*n+l], temp)
                }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol*n+icol] == 0.0)
            nrerror("GAUSSJ: Singular Matrix-2");
        pivinv = 1.0f / a[icol*n+icol];
        a[icol*n+icol] = (T)1;
        for (l = 1;l <= n;l++)
            a[icol*n+l] = (T)(pivinv * a[icol*n+l]);
        for (l = 1;l <= m;l++)
            b[icol*n+l] = (T)(pivinv * b[icol*n+l]);
        for (ll = 1;ll <= n;ll++)
            if (ll != icol)
            {
                dum = a[ll*n+icol];
                a[ll*n+icol] = (T)0;
                for (l = 1;l <= n;l++)
                    a[ll*n+l] -= a[icol*n+l] * dum;
                for (l = 1;l <= m;l++)
                    b[ll*n+l] -= b[icol*n+l] * dum;
            }
    }
    for (l = n;l >= 1;l--)
    {
        if (indxr[l] != indxc[l])
            for (k = 1;k <= n;k++)
                SWAP(a[k*n+indxr[l]], a[k*n+indxc[l]], temp);
    }
    free_Tvector(ipiv, 1, n);
    free_Tvector(indxr, 1, n);
    free_Tvector(indxc, 1, n);
}


#endif

/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#include "svd_helper.h"
#include <src/jaz/util/index_sort.h>

void SvdHelper::decompose(
        const Matrix2D<double>& A,
        Matrix2D<double>& U,
        Matrix1D<double>& S,
        Matrix2D<double>& Vt)
{
    Matrix2D<double> U0, Vt0;
    Matrix1D<double> S0;

    svdcmp(A, U0, S0, Vt0);

    const int rc = A.mdimy;
    const int cc = A.mdimx;

    std::vector<double> Svec(cc);

    for (int i = 0; i < cc; i++)
    {
        Svec[i] = S0(i);
    }

    std::vector<int> order = IndexSort<double>::sortIndices(Svec);

    U = Matrix2D<double>(rc,cc);
    S = Matrix1D<double>(cc);
    Vt = Matrix2D<double>(cc,cc);

    for (int i = 0; i < cc; i++)
    {
        const int j = order[cc - i - 1];

        for (int c = 0; c < cc; c++)
        {
            Vt(c,i) = Vt0(c,j);
        }

        S(i) = S0(j);

        for (int r = 0; r < rc; r++)
        {
            U(r,i) = U0(r,j);
        }
    }
}


/* SVD code: stolen from XMIPP and slightly altered */
/*
#define SVDMAXITER 1000000

void SvdHelper::svdcmp(const Matrix2D< double >& a,
					   Matrix1D< double >& w,
					   Matrix2D< double >& v)
{
	Matrix2D<double> u = a;
	w.initZeros(u.mdimx);
	v.initZeros(u.mdimx, u.mdimx);
		
	
	double* U = u.mdata;
	const int rows = u.mdimy;
	const int cols = u.mdimx;
	double* W = w.vdata;
	double* V = v.mdata;
		
    double *rv1 = (double *)NULL;
	
    double c, f, h, s;
    double x, y, z;
    long l = 0;
    int     maxIterations = SVDMAXITER;

	std::vector<double> rv1_vec(cols*cols);
	rv1 = rv1_vec.data();
	
    double g = 0.0;
	double scale = 0.0;
	double norm = 0.0;
	
    for (int i = 0; i < cols; i++)
    {
        l = i + 1;
        rv1[i] = scale * g;
		
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		
        if (i < rows)
        {
            for (int k = i; (k < rows); k++)
            {
                scale += std::abs(U[k * cols + i]);
            }
			
            if (scale != 0.0)
            {
				s = 0.0;
				
                for (int k = i; k < rows; k++)
                {
                    U[k * cols + i] /= scale;
                    s += U[k * cols + i] * U[k * cols + i];
                }
				
                f = U[i * cols + i];
                g = (f >= 0.0)? -sqrt(s) : sqrt(s);
                h = f * g - s;
				
                U[i * cols + i] = f - g;
				
                for (int j = l; j < cols; j++)
                {
					s = 0.0;
					
                    for (int k = i; k < rows; k++)
                    {
                        s += U[k * cols + i] * U[k * cols + j];
                    }
					
                    f = s / h;
					
                    for (int k = i; k < rows; k++)
                    {
                        U[k * cols + j] += f * U[k * cols + i];
                    }
                }
				
                for (int k = i; k < rows; k++)
                {
                    U[k * cols + i] *= scale;
                }
            }
        }
		
        W[i] = scale * g;
		
		g = 0.0;
		s = 0.0;
		scale = 0.0;
		
        if ((i < rows) && (i != (cols - 1L)))
        {
            for (int k = l; (k < cols); k++)
            {
                scale += std::abs(U[i * cols + k]);
            }
			
            if (scale != 0.0)
            {
                for (int k = l; (k < cols); k++)
                {
                    U[i * cols + k] /= scale;
                    s += U[i * cols + k] * U[i * cols + k];
                }
				
                f = U[i * cols + l];
                g = (0.0 <= f) ? (-sqrt(s)) : (sqrt(s));
                h = f * g - s;
                U[i * cols + l] = f - g;
				
                for (int k = l; (k < cols); k++)
                {
                    rv1[k] = U[i * cols + k] / h;
                }
				
                for (int j = l; (j < rows); j++)
                {
					s = 0.0;
					
                    for (int k = l; (k < cols); k++)
                    {
                        s += U[j * cols + k] * U[i * cols + k];
                    }
					
                    for (int k = l; (k < cols); k++)
                    {
                        U[j * cols + k] += s * rv1[k];
                    }
                }
				
                for (int k = l; (k < cols); k++)
                {
                    U[i * cols + k] *= scale;
                }
            }
        }
		
		const double nn = std::abs(W[i]) + std::abs(rv1[i]);				
        norm = std::max(nn, norm);
    }
	
    for (int i = cols - 1; i >= 0; i--)
    {
        if (i < cols - 1)
        {
            if (g != 0.0)
            {
                for (int j = l; j < cols; j++)
                {
                    V[j * cols + i] = U[i * cols + j] / (U[i * cols + l] * g);
                }
				
                for (int j = l; j < cols; j++)
                {
					s = 0.0;
					
                    for (int k = l; (k < cols); k++)
                    {
                        s += U[i * cols + k] * V[k * cols + j];
                    }
					
                    for (int k = l; (k < cols); k++)
                    {
                        if (s != 0.0)
                        {
                            V[k * cols + j] += s * V[k * cols + i];
                        }
                    }
                }
            }
			
            for (int j = l; (j < cols); j++)
            {
                V[i * cols + j] = 0.0;
				V[j * cols + i] = 0.0;
            }
        }
		
        V[i * cols + i] = 1.0;
        g = rv1[i];
        l = i;
    }
	
    for (int i = (rows < cols) ? (rows - 1) : (cols - 1); i >= 0; i--)
    {
        l = i + 1L;
        g = W[i];
		
        for (int j = l; (j < cols); j++)
        {
            U[i * cols + j] = 0.0;
        }
		
        if (g != 0.0)
        {
            g = 1.0 / g;
			
            for (int j = l; (j < cols); j++)
            {
				s = 0.0;
				
                for (int k = l; (k < rows); k++)
                {
                    s += U[k * cols + i] * U[k * cols + j];
                }
				
                f = s * g / U[i * cols + i];
				
                for (int k = i; (k < rows); k++)
                {
                    if (f != 0.0)
                    {
                        U[k * cols + j] += f * U[k * cols + i];
                    }
                }
            }
			
            for (int j = i; (j < rows); j++)
            {
                U[j * cols + i] *= g;
            }
        }
        else
        {
            for (int j = i; (j < rows); j++)
            {
                U[j * cols + i] = 0.0;
            }
        }
		
        U[i * cols + i] += 1.0;
    }
	
    for (int k = cols - 1; k >= 0; k--)
    {
        for (int its = 1; its <= maxIterations; its++)
        {
            bool any_zero = true;
			
			int nm = 0;
			
            for (int l = k; l >= 0; l--)
            {
                nm = l - 1;
				
                if ((std::abs(rv1[l]) + norm) == norm)
                {
                    any_zero = false;
                    break;
                }
				
                if ((std::abs(W[nm]) + norm) == norm)
                {
                    break;
                }
            }
			
            if (any_zero)
            {
                c = 0.0;
                s = 1.0;
				
                for (int i = l; (i <= k); i++)
                {
                    f = s * rv1[i];
                    rv1[i] *= c;
					
                    if ((std::abs(f) + norm) == norm)
                    {
                        break;
                    }
					
                    g = W[i];
                    h = f*f + g*g;
                    W[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
					
                    for (int j = 0; j < rows; j++)
                    {
                        y = U[j * cols + nm];
                        z = U[j * cols + i];
                        U[j * cols + nm] = y * c + z * s;
                        U[j * cols + i] = z * c - y * s;
                    }
                }
            }
			
            z = W[k];
			
            if (l == k)
            {
                if (z < 0.0)
                {
                    W[k] = -z;
					
                    for (int j = 0; j < cols; j++)
                    {
                        V[j * cols + k] = -V[j * cols + k];
                    }
                }
				
                break;
            }
			
            if (its == maxIterations)
            {
                return;
            }
			
            x = W[l];
            nm = k - 1;
            y = W[nm];
            g = rv1[nm];
            h = rv1[k];
			
            f = ((y - z)*(y + z) + (g - h)*(g + h)) / (2.0 * h * y);
            g = f*f + 1.0;
            f = ((x - z)*(x + z) + h*((y / (f + ((f >= 0.0) ? (std::abs(g))
                                                : (-std::abs(g))))) - h)) / x;
            c = s = 1.0;
			
            for (int j = l; j <= nm; j++)
            {
                int i = j + 1L;
				
                g = rv1[i];
                y = W[i];
                h = s * g;
                g = c * g;
				
                z = f*f + h*h;
				
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
				
                for (int jj = 0; jj < cols; jj++)
                {
                    x = V[jj * cols + j];
                    z = V[jj * cols + i];
                    V[jj * cols + j] = x * c + z * s;
                    V[jj * cols + i] = z * c - x * s;
                }
				
                z = f*f + h*h;
                W[j] = z;
				
                if (z != 0.0)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
				
                f = c * g + s * y;
                x = c * y - s * g;
				
                for (int jj = 0; jj < rows; jj++)
                {
                    y = U[jj * cols + j];
                    z = U[jj * cols + i];
                    U[jj * cols + j] = y * c + z * s;
                    U[jj * cols + i] = z * c - y * s;
                }
            }
			
            rv1[l] = 0.0;
            rv1[k] = f;
            W[k] = x;
        }
    }
}
*/

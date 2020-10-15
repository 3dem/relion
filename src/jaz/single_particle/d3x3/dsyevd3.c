// ----------------------------------------------------------------------------
// Numerical diagonalization of 3x3 matrcies
// Copyright (C) 2006  Joachim Kopp
// ----------------------------------------------------------------------------
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
// ----------------------------------------------------------------------------
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include "dsyev2.h"
#include "slvsec3.h"
#include "dsytrd3.h"
#include "dsyevd3.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 


// ----------------------------------------------------------------------------
int dsyevd3(double A[3][3], double Q[3][3], double w[3])
// ----------------------------------------------------------------------------
// Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
// matrix A using Cuppen's Divide & Conquer algorithm.
// The function accesses only the diagonal and upper triangular parts of A.
// The access is read-only.
// ----------------------------------------------------------------------------
// Parameters:
//   A: The symmetric input matrix
//   Q: Storage buffer for eigenvectors
//   w: Storage buffer for eigenvalues
// ----------------------------------------------------------------------------
// Return value:
//   0: Success
//  -1: Error
// ----------------------------------------------------------------------------
// Dependencies:
//   dsyev2(), slvsec3(), dsytrd3()
// ----------------------------------------------------------------------------
{
  const int n = 3;
  double R[3][3];                // Householder transformation matrix
  double P[3][3];                // Unitary transformation matrix which diagonalizes D + w w^T
  double e[2];                   // Off-diagonal elements after Householder transformation
  double d[3];                   // Eigenvalues of split matrix in the "divide" step)
  double c, s;                   // Eigenvector of 2x2 block in the "divide" step
  double z[3];                   // Numerators of secular equation / Updating vector
  double t;                      // Miscellaenous temporary stuff

  // Initialize Q
#ifndef EVALS_ONLY
  memset(Q, 0.0, 9*sizeof(double));
#endif
  
  // Transform A to real tridiagonal form by the Householder method
  dsytrd3(A, R, w, e);
 
  
  // "Divide"
  // --------
  
  // Detect matrices that factorize to avoid multiple eigenvalues in the Divide/Conquer algorithm
  for (int i=0; i < n-1; i++)
  {
    t = fabs(w[i]) + fabs(w[i+1]);
    if (fabs(e[i]) <= 8.0*DBL_EPSILON*t)
    {
      if (i == 0)
      {
        dsyev2(w[1], e[1], w[2], &d[1], &d[2], &c, &s);
        w[1] = d[1];
        w[2] = d[2];
#ifndef EVALS_ONLY
        Q[0][0] = 1.0;
        for (int j=1; j < n; j++)
        {
          Q[j][1] = s*R[j][2] + c*R[j][1];
          Q[j][2] = c*R[j][2] - s*R[j][1];
        }
#endif
      }
      else
      {
        dsyev2(w[0], e[0], w[1], &d[0], &d[1], &c, &s);
        w[0] = d[0];
        w[1] = d[1];
#ifndef EVALS_ONLY
        Q[0][0]   = c;
        Q[0][1]   = -s;
        Q[1][0]   = R[1][1]*s;
        Q[1][1]   = R[1][1]*c;
        Q[1][2]   = R[1][2];
        Q[2][0]   = R[2][1]*s;
        Q[2][1]   = R[2][1]*c;
        Q[2][2]   = R[2][2];
#endif
      }

      return 0;
    }
  }
  
  // Calculate eigenvalues and eigenvectors of 2x2 block
  dsyev2(w[1]-e[0], e[1], w[2], &d[1], &d[2], &c, &s);
  d[0] = w[0] - e[0];

  
  // "Conquer"
  // ---------

  // Determine coefficients of secular equation
  z[0] = e[0];
  z[1] = e[0] * SQR(c);
  z[2] = e[0] * SQR(s);

  // Call slvsec3 with d sorted in ascending order. We make
  // use of the fact that dsyev2 guarantees d[1] >= d[2].
  if (d[0] < d[2])
    slvsec3(d, z, w, P, 0, 2, 1);
  else if (d[0] < d[1])
    slvsec3(d, z, w, P, 2, 0, 1);
  else
    slvsec3(d, z, w, P, 2, 1, 0);

#ifndef EVALS_ONLY
  // Calculate eigenvectors of matrix D + beta * z * z^t and store them in the
  // columns of P
  z[0] = sqrt(fabs(e[0]));
  z[1] = c * z[0];
  z[2] = -s * z[0];

  // Detect duplicate elements in d to avoid division by zero
  t = 8.0*DBL_EPSILON*(fabs(d[0]) + fabs(d[1]) + fabs(d[2]));
  if (fabs(d[1] - d[0]) <= t)
  {
    for (int j=0; j < n; j++)
    {
      if (P[0][j] * P[1][j] <= 0.0)
      {
        P[0][j] = z[1];
        P[1][j] = -z[0];
        P[2][j] = 0.0;
      }
      else
        for (int i=0; i < n; i++)
          P[i][j] = z[i]/P[i][j];
    }
  }
  else if (fabs(d[2] - d[0]) <= t)
  {
    for (int j=0; j < n; j++)
    {
      if (P[0][j] * P[2][j] <= 0.0)
      {
        P[0][j] = z[2];
        P[1][j] = 0.0;
        P[2][j] = -z[0];
      }
      else
        for (int i=0; i < n; i++)
          P[i][j] = z[i]/P[i][j];
    }
  }
  else
  {
    for (int j=0; j < n; j++)
      for (int i=0; i < n; i++)
      {
        if (P[i][j] == 0.0)
        {
          P[i][j]       = 1.0;
          P[(i+1)%n][j] = 0.0;
          P[(i+2)%n][j] = 0.0;
          break;
        }
        else
          P[i][j] = z[i]/P[i][j];
      }
  }

  // Normalize eigenvectors of D + beta * z * z^t
  for (int j=0; j < n; j++)
  {
    t = SQR(P[0][j]) + SQR(P[1][j]) + SQR(P[2][j]);
    t = 1.0 / sqrt(t);
    for (int i=0; i < n; i++)
      P[i][j] *= t;
  }
  
  // Undo diagonalization of 2x2 block
  for (int j=0; j < n; j++)
  {
    t       = P[1][j];
    P[1][j] = c*t - s*P[2][j];
    P[2][j] = s*t + c*P[2][j];
  }

  // Undo Householder transformation
  for (int j=0; j < n; j++)
    for (int k=0; k < n; k++)
    {
      t = P[k][j];
      for (int i=0; i < n; i++)
        Q[i][j] += t * R[i][k];
    }
#endif

  return 0;
}




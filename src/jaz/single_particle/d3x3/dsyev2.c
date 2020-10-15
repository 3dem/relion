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
#include <math.h>
#include "dsyev2.h"

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 
#define SQR_ABS(x)  (SQR(creal(x)) + SQR(cimag(x)))  // |x|^2


// ----------------------------------------------------------------------------
inline void dsyev2(double A, double B, double C, double *rt1, double *rt2,
                   double *cs, double *sn)
// ----------------------------------------------------------------------------
// Calculates the eigensystem of a real symmetric 2x2 matrix
//    [ A  B ]
//    [ B  C ]
// in the form
//    [ A  B ]  =  [ cs  -sn ] [ rt1   0  ] [  cs  sn ]
//    [ B  C ]     [ sn   cs ] [  0   rt2 ] [ -sn  cs ]
// where rt1 >= rt2. Note that this convention is different from the one used
// in the LAPACK routine DLAEV2, where |rt1| >= |rt2|.
// ----------------------------------------------------------------------------
{
  double sm = A + C;
  double df = A - C;
  double rt = sqrt(SQR(df) + 4.0*B*B);
  double t;

  if (sm > 0.0)
  {
    *rt1 = 0.5 * (sm + rt);
    t = 1.0/(*rt1);
    *rt2 = (A*t)*C - (B*t)*B;
  }
  else if (sm < 0.0)
  {
    *rt2 = 0.5 * (sm - rt);
    t = 1.0/(*rt2);
    *rt1 = (A*t)*C - (B*t)*B;
  }
  else       // This case needs to be treated separately to avoid div by 0
  {
    *rt1 = 0.5 * rt;
    *rt2 = -0.5 * rt;
  }

  // Calculate eigenvectors
  if (df > 0.0)
    *cs = df + rt;
  else
    *cs = df - rt;

  if (fabs(*cs) > 2.0*fabs(B))
  {
    t   = -2.0 * B / *cs;
    *sn = 1.0 / sqrt(1.0 + SQR(t));
    *cs = t * (*sn);
  }
  else if (fabs(B) == 0.0)
  {
    *cs = 1.0;
    *sn = 0.0;
  }
  else
  {
    t   = -0.5 * (*cs) / B;
    *cs = 1.0 / sqrt(1.0 + SQR(t));
    *sn = t * (*cs);
  }

  if (df > 0.0)
  {
    t   = *cs;
    *cs = -(*sn);
    *sn = t;
  }
}


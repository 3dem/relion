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
#include <float.h>
#include "slvsec3.h"

// Constants
#define M_SQRT3    1.73205080756887729352744634151   // sqrt(3)

// Macros
#define SQR(x)      ((x)*(x))                        // x^2 


// ----------------------------------------------------------------------------
void slvsec3(double d[3], double z[3], double w[3],
                    double R[3][3], int i0, int i1, int i2)
// ----------------------------------------------------------------------------
// Finds the three roots w_j of the secular equation
//   f(w_j) = 1 + Sum[ z_i / (d_i - w_j) ]  ==  0.
// It is assumed that d_0 <= d_1 <= d_2, and that all z_i have the same sign.
// The arrays P_i will contain the information required for the calculation
// of the eigenvectors:
//   P_ij = d_i - w_j.
// These differences can be obtained with better accuracy from intermediate
// results.
// ----------------------------------------------------------------------------
{
  double a[4];            // Bounds of the intervals bracketing the roots
  double delta;           // Shift of the d_i which ensures better accuracy
  double dd[3];           // Shifted coefficients dd_i = d_i - delta
  double xl, xh;          // Interval which straddles the current root. f(xl) < 0, f(xh) > 0
  double x;               // Current estimates for the root
  double x0[3];           // Analytically calculated roots, used as starting values
  double F, dF;           // Function value f(x) and derivative f'(x)
  double dx, dxold;       // Current and last stepsizes
  double error;           // Numerical error estimate, used for termination condition
  double t[3];            // Temporary storage used for evaluating f
  double alpha, beta, gamma;       // Coefficients of polynomial f(x) * Product [ d_i - x ]
  double p, sqrt_p, q, c, s, phi;  // Intermediate results of analytical calculation
  
  // Determine intervals which must contain the roots
  if (z[0] > 0)
  {
    a[0] = d[i0];
    a[1] = d[i1];
    a[2] = d[i2];
    a[3] = fabs(d[0] + 3.0*z[0]) + fabs(d[1] + 3.0*z[1]) + fabs(d[2] + 3.0*z[2]);
  }
  else
  {    
    a[0] = -fabs(d[0] + 3.0*z[0]) - fabs(d[1] + 3.0*z[1]) - fabs(d[2] + 3.0*z[2]);
    a[1] = d[i0];
    a[2] = d[i1];
    a[3] = d[i2];
  }

  // Calculate roots of f(x) = 0 analytically (analogous to ZHEEVC3)
  t[0]  = d[1]*d[2];
  t[1]  = d[0]*d[2];
  t[2]  = d[0]*d[1];
  gamma = t[0]*d[0] + (z[0]*t[0] + z[1]*t[1] + z[2]*t[2]);    // Coefficients
  beta  = (z[0]*(d[1]+d[2]) + z[1]*(d[0]+d[2]) + z[2]*(d[0]+d[1]))
           + (t[0] + t[1] + t[2]);
  alpha = (z[0] + z[1] + z[2]) + (d[0] + d[1] + d[2]);
  
  p = SQR(alpha) - 3.0*beta;    // Transformation that removes the x^2 term
  q = alpha*(p - (3.0/2.0)*beta) + (27.0/2.0)*gamma;
  sqrt_p = sqrt(fabs(p));

  phi = 27.0 * ( 0.25*SQR(beta)*(p - beta) - gamma*(q - 27.0/4.0*gamma));
  phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
  c = sqrt_p*cos(phi);
  s = (1.0/M_SQRT3)*sqrt_p*fabs(sin(phi));

  x0[0] = x0[1] = x0[2] = (1.0/3.0)*(alpha - c);
  if (c > s)             // Make sure the roots are in ascending order
  {
    x0[0] -= s;
    x0[1] += s;
    x0[2] += c;
  }
  else if (c < -s)
  {
    x0[0] += c;
    x0[1] -= s;
    x0[2] += s;
  }
  else
  {
    x0[0] -= s;
    x0[1] += c;
    x0[2] += s;
  }

  // Refine roots with a combined Bisection/Newton-Raphson method
  for (int i=0; i < 3; i++)
  {
    xl = a[i];               // Lower bound of bracketing interval
    xh = a[i+1];             // Upper bound of bracketing interval
    dx = dxold = 0.5 * (xh - xl);

    // Make sure that xl != xh
    if (dx == 0.0)
    {
      w[i] = xl;
      for (int j=0; j < 3; j++)
        R[j][i] = d[j] - xl;
      continue;
    }
    
    // Shift the root close to zero to achieve better accuracy
    if (x0[i] >= xh)
    {
      delta = xh;
      x     = -dx;
      for (int j=0; j < 3; j++)
      {
        dd[j]   = d[j] - delta;
        R[j][i] = dd[j] - x;
      }
    }
    else if (x0[i] <= xl)
    {
      delta = xl;
      x     = dx;
      for (int j=0; j < 3; j++)
      {
        dd[j]   = d[j] - delta;
        R[j][i] = dd[j] - x;
      }
    }
    else
    {
      delta = x0[i];
      x     = 0.0;
      for (int j=0; j < 3; j++)
        R[j][i] = dd[j] = d[j] - delta;
    }
    xl -= delta;
    xh -= delta;
   
    // Make sure that f(xl) < 0 and f(xh) > 0 
    if (z[0] < 0.0)
    {
      double t = xh;
      xh = xl;
      xl = t;
    }

    // Main iteration loop
    for (int nIter=0; nIter < 500; nIter++)
    {
      // Evaluate f and f', and calculate an error estimate
      F     = 1.0;
      dF    = 0.0;
      error = 1.0;
      for (int j=0; j < 3; j++)
      {
        t[0]   = 1.0 / R[j][i];
        t[1]   = z[j] * t[0];
        t[2]   = t[1] * t[0];
        F     += t[1];
        error += fabs(t[1]);
        dF    += t[2];
      }

      // Check for convergence 
      if (fabs(F) <= DBL_EPSILON * (8.0 * error + fabs(x*dF)))
        break;

      // Adjust interval boundaries
      if (F < 0.0)
        xl   = x;
      else
        xh   = x;

      // Check, whether Newton-Raphson would converge fast enough. If so,
      // give it a try. If not, or if it would run out of bounds, use bisection
      if (fabs(2.0 * F) < fabs(dxold * dF))
      {
        dxold = dx;
        dx    = F / dF;
        x     = x - dx;
        if ((x - xh) * (x - xl) >= 0.0)
        {
          dx = 0.5 * (xh - xl);
          x  = xl + dx;
        }
      }
      else
      {
        dx = 0.5 * (xh - xl);
        x  = xl + dx;
      }

      // Prepare next iteration
      for (int j=0; j < 3; j++)
        R[j][i] = dd[j] - x;
    }
     
    // Un-shift result
    w[i] = x + delta;
  }
}

 


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
#include "src/funcs.h"
#include "src/args.h"

#include <stdio.h>
#include <fstream>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>
#include <complex>
#include <fstream>
#include <typeinfo>

void fitStraightLine(const std::vector<fit_point2D> &points, RFLOAT &slope, RFLOAT &intercept, RFLOAT &corr_coeff)
{
	// From: http://mathworld.wolfram.com/LeastSquaresFitting.html
	// ss_xx = Sum_i x_i^2 - n ave_x^2
	// ss_yy = Sum_i y_i^2 - n ave_y^2
	// ss_xy = Sum_i x_i * y_i - n ave_x n_ave_y
	// slope = xx_xy / ss_xx
	// intercept = ave_y - slope * ave_x
	// corr_coeff = ss_xy^2 / (ss_xx * ss_yy)
	RFLOAT ss_xy = 0.;
	RFLOAT ss_xx = 0.;
	RFLOAT ss_yy = 0.;
	RFLOAT ave_x = 0.;
	RFLOAT ave_y = 0.;
	RFLOAT sum_w = 0.;
	for (int i = 0; i < points.size(); i++)
	{
		ave_x += points[i].w * points[i].x;
		ave_y += points[i].w * points[i].y;
		sum_w += points[i].w;
		ss_xx += points[i].w * points[i].x * points[i].x;
		ss_yy += points[i].w * points[i].y * points[i].y;
		ss_xy += points[i].w * points[i].x * points[i].y;
	}
	ave_x /= sum_w;
	ave_y /= sum_w;
	ss_xx -= sum_w * ave_x * ave_x;
	ss_yy -= sum_w * ave_y * ave_y;
	ss_xy -= sum_w * ave_x * ave_y;

	//std::cerr << " ss_xx= " << ss_xx << " ss_yy= " << ss_yy << " ss_xy= " << ss_xy << std::endl;
	//std::cerr << " sum_w= " << sum_w << " ave_x= " << ave_x << " ave_y= " << ave_y << std::endl;
	if (ss_xx > 0.)
	{
		slope = ss_xy / ss_xx;
		intercept = ave_y - slope * ave_x;
		corr_coeff = ss_xy * ss_xy / (ss_xx * ss_yy);
	}
	else
	{
		intercept = slope = corr_coeff = 0.;
	}
}

void fitLeastSquaresPlane(const std::vector<fit_point3D> & points, RFLOAT &plane_a, RFLOAT &plane_b, RFLOAT &plane_c)
{
    RFLOAT  D = 0;
    RFLOAT  E = 0;
    RFLOAT  F = 0;
    RFLOAT  G = 0;
    RFLOAT  H = 0;
    RFLOAT  I = 0;
    RFLOAT  J = 0;
    RFLOAT  K = 0;
    RFLOAT  L = 0;
    RFLOAT  W2 = 0;
    RFLOAT  error = 0;
    RFLOAT  denom = 0;

    for (int i = 0; i < points.size(); i++)
    {
        W2 = points[i].w * points[i].w;
        D += points[i].x * points[i].x * W2 ;
        E += points[i].x * points[i].y * W2 ;
        F += points[i].x * W2 ;
        G += points[i].y * points[i].y * W2 ;
        H += points[i].y * W2 ;
        I += 1 * W2 ;
        J += points[i].x * points[i].z * W2 ;
        K += points[i].y * points[i].z * W2 ;
        L += points[i].z * W2 ;
    }

    denom = F * F * G - 2 * E * F * H + D * H * H + E * E * I - D * G * I;

    // X axis slope
    plane_a = (H * H * J - G * I * J + E * I * K + F * G * L - H * (F * K + E * L)) / denom;
    // Y axis slope
    plane_b = (E * I * J + F * F * K - D * I * K + D * H * L - F * (H * J + E * L)) / denom;
    // Z axis intercept
    plane_c = (F * G * J - E * H * J - E * F * K + D * H * K + E * E * L - D * G * L) / denom;
}



/* Value of a blob --------------------------------------------------------- */
RFLOAT kaiser_value(RFLOAT r, RFLOAT a, RFLOAT alpha, int m)
{
    RFLOAT rda, rdas, arg, w;
    rda = r / a;
    rdas = rda * rda;
    if (rdas <= 1.0)
    {
        arg = alpha * sqrt(1.0 - rdas);
        if (m == 0)
        {
            w = bessi0(arg) / bessi0(alpha);
        }
        else if (m == 1)
        {
            w = sqrt (1.0 - rdas);
            if (alpha != 0.0)
                w *= bessi1(arg) / bessi1(alpha);
        }
        else if (m == 2)
        {
            w = sqrt (1.0 - rdas);
            w = w * w;
            if (alpha != 0.0)
                w *= bessi2(arg) / bessi2(alpha);
        }
        else if (m == 3)
        {
            w = sqrt (1.0 - rdas);
            w = w * w * w;
            if (alpha != 0.0)
                w *= bessi3(arg) / bessi3(alpha);
        }
        else if (m == 4)
        {
            w = sqrt (1.0 - rdas);
            w = w * w * w *w;
            if (alpha != 0.0)
                w *= bessi4(arg) / bessi4(alpha);
        }
        else REPORT_ERROR("m out of range in kaiser_value()");
    }
    else
        w = 0.0;
    return w;
}
/* Line integral through a blob -------------------------------------------- */
/* Value of line integral through Kaiser-Bessel radial function
   (n >=2 dimensions) at distance s from center of function.
   Parameter m = 0, 1, or 2. */
RFLOAT kaiser_proj(RFLOAT s, RFLOAT a, RFLOAT alpha, int m)
{
    RFLOAT sda, sdas, w, arg, p;
    sda = s / a;
    sdas = sda * sda;
    w = 1.0 - sdas;
    if (w > 1.0e-10)
    {
        arg = alpha * sqrt(w);
        if (m == 0)
        {
            if (alpha == 0.0)
                p = 2.0 * a * sqrt(w);
            else
                p = (2.0 * a / alpha) * sinh(arg) / bessi0(alpha);
        }
        else if (m == 1)
        {
            if (alpha == 0.0)
                p = 2.0 * a * w * sqrt(w) * (2.0 / 3.0);
            else
                p = (2.0 * a / alpha) * sqrt(w) * (cosh(arg) - sinh(arg) / arg)
                    / bessi1(alpha);
        }
        else if (m == 2)
        {
            if (alpha == 0.0)
                p = 2.0 * a * w * w * sqrt(w) * (8.0 / 15.0);
            else
                p = (2.0 * a / alpha) * w *
                    ((3.0 / (arg * arg) + 1.0) * sinh(arg) - (3.0 / arg) * cosh(arg)) / bessi2(alpha);
        }
        else REPORT_ERROR("m out of range in kaiser_proj()");
    }
    else
        p = 0.0;
    return p;
}
/* Fourier value of a blob ------------------------------------------------- */
RFLOAT kaiser_Fourier_value(RFLOAT w, RFLOAT a, RFLOAT alpha, int m)
{
    RFLOAT sigma = sqrt(ABS(alpha * alpha - (2. * PI * a * w) * (2. * PI * a * w)));
    if (m == 2)
    {
        if (2.*PI*a*w > alpha)
            return  pow(2.*PI, 3. / 2.)*pow(a, 3.)*pow(alpha, 2.)*bessj3_5(sigma)
                    / (bessi0(alpha)*pow(sigma, 3.5));
        else
            return  pow(2.*PI, 3. / 2.)*pow(a, 3.)*pow(alpha, 2.)*bessi3_5(sigma)
                    / (bessi0(alpha)*pow(sigma, 3.5));
    }
    else if (m == 0)
    {
        if (2*PI*a*w > alpha)
            return  pow(2.*PI, 3. / 2.)*pow(a, 3)*bessj1_5(sigma)
                    / (bessi0(alpha)*pow(sigma, 1.5));
        else
            return  pow(2.*PI, 3. / 2.)*pow(a, 3)*bessi1_5(sigma)
                    / (bessi0(alpha)*pow(sigma, 1.5));
    }
    else
    	REPORT_ERROR("m out of range in kaiser_Fourier_value()");
}
/* Volume integral of a blob ----------------------------------------------- */
RFLOAT  basvolume(RFLOAT a, RFLOAT alpha, int m, int n)
{
    RFLOAT  hn, tpi, v;
    hn = 0.5 * n;
    tpi = 2.0 * PI;
    if (alpha == 0.0)
    {
        if ((n / 2)*2 == n)           /* n even                               */
            v = pow(tpi, hn) * in_zeroarg(n / 2 + m) / in_zeroarg(m);
        else                        /* n odd                                */
            v = pow(tpi, hn) * inph_zeroarg(n / 2 + m) / in_zeroarg(m);
    }
    else
    {                        /* alpha > 0.0                          */
        if ((n / 2)*2 == n)           /* n even                               */
            v = pow(tpi / alpha, hn) * i_n(n / 2 + m, alpha) / i_n(m, alpha);
        else                        /* n odd                                */
            v = pow(tpi / alpha, hn) * i_nph(n / 2 + m, alpha) / i_n(m, alpha);
    }
    return v * pow(a, (RFLOAT)n);
}
/* Bessel function I_n (x),  n = 0, 1, 2, ...
 Use ONLY for small values of n     */
RFLOAT i_n(int n, RFLOAT x)
{
    int i;
    RFLOAT i_ns1, i_n, i_np1;
    if (n == 0)   return bessi0(x);
    if (n == 1)   return bessi1(x);
    if (x == 0.0) return 0.0;
    i_ns1 = bessi0(x);
    i_n   = bessi1(x);
    for (i = 1; i < n; i++)
    {
        i_np1 = i_ns1 - (2 * i) / x * i_n;
        i_ns1 = i_n;
        i_n   = i_np1;
    }
    return i_n;
}
/*.....Bessel function I_(n+1/2) (x),  n = 0, 1, 2, ..........................*/
RFLOAT i_nph(int n, RFLOAT x)
{
    int i;
    RFLOAT r2dpix;
    RFLOAT i_ns1, i_n, i_np1;
    if (x == 0.0) return 0.0;
    r2dpix = sqrt(2.0 / (PI * x));
    i_ns1 = r2dpix * cosh(x);
    i_n   = r2dpix * sinh(x);
    for (i = 1; i <= n; i++)
    {
        i_np1 = i_ns1 - (2 * i - 1) / x * i_n;
        i_ns1 = i_n;
        i_n   = i_np1;
    }
    return i_n;
}
/*....Limit (z->0) of (1/z)^n I_n(z)..........................................*/
RFLOAT in_zeroarg(int n)
{
    int i;
    RFLOAT fact;
    fact = 1.0;
    for (i = 1; i <= n; i++)
    {
        fact *= 0.5 / i;
    }
    return fact;
}
/*.......Limit (z->0) of (1/z)^(n+1/2) I_(n+1/2) (z)..........................*/
RFLOAT inph_zeroarg(int n)
{
    int i;
    RFLOAT fact;
    fact = 1.0;
    for (i = 1; i <= n; i++)
    {
        fact *= 1.0 / (2 * i + 1.0);
    }
    return fact*sqrt(2.0 / PI);
}
/* Zero freq --------------------------------------------------------------- */
RFLOAT blob_freq_zero(struct blobtype b)
{
    return sqrt(b.alpha*b.alpha + 6.9879*6.9879) / (2*PI*b.radius);
}
/* Attenuation ------------------------------------------------------------- */
RFLOAT blob_att(RFLOAT w, struct blobtype b)
{
    return blob_Fourier_val(w, b) / blob_Fourier_val(0, b);
}
/* Number of operations ---------------------------------------------------- */
RFLOAT blob_ops(RFLOAT w, struct blobtype b)
{
    return pow(b.alpha*b.alpha + 6.9879*6.9879, 1.5) / b.radius;
}

/* Gaussian value ---------------------------------------------------------- */
RFLOAT gaussian1D(RFLOAT x, RFLOAT sigma, RFLOAT mu)
{
    x -= mu;
    return 1 / sqrt(2*PI*sigma*sigma)*exp(-0.5*((x / sigma)*(x / sigma)));
}

/* t-student value -------------------------------------------------------- */
RFLOAT tstudent1D(RFLOAT x, RFLOAT df, RFLOAT sigma, RFLOAT mu)
{
    x -= mu;
    RFLOAT norm = exp(gammln((df+1.)/2.)) / exp(gammln(df/2.));
    norm /= sqrt(df*PI*sigma*sigma);
    return norm * pow((1 + (x/sigma)*(x/sigma)/df),-((df+1.)/2.));

}

RFLOAT gaussian2D(RFLOAT x, RFLOAT y, RFLOAT sigmaX, RFLOAT sigmaY,
                  RFLOAT ang, RFLOAT muX, RFLOAT muY)
{
    // Express x,y in the gaussian internal coordinates
    x -= muX;
    y -= muY;
    RFLOAT xp = cos(ang) * x + sin(ang) * y;
    RFLOAT yp = -sin(ang) * x + cos(ang) * y;

    // Now evaluate
    return 1 / sqrt(2*PI*sigmaX*sigmaY)*exp(-0.5*((xp / sigmaX)*(xp / sigmaX) +
                                            (yp / sigmaY)*(yp / sigmaY)));
}

/* ICDF Gaussian ----------------------------------------------------------- */
RFLOAT icdf_gauss(RFLOAT p)
{
    const RFLOAT c[] =
        {
            2.515517, 0.802853, 0.010328
        };
    const RFLOAT d[] =
        {
            1.432788, 0.189269, 0.001308
        };
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        RFLOAT t=sqrt(-2.0*log(p));
        RFLOAT z=t - ((c[2]*t + c[1])*t + c[0]) /
                 (((d[2]*t + d[1])*t + d[0])*t + 1.0);
        return -z;
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        RFLOAT t=sqrt(-2.0*log(1-p));
        RFLOAT z=t - ((c[2]*t + c[1])*t + c[0]) /
                 (((d[2]*t + d[1])*t + d[0])*t + 1.0);
        return z;
    }
}

/* CDF Gaussian ------------------------------------------------------------ */
RFLOAT cdf_gauss(RFLOAT x)
{
    return 0.5 * (1. + erf(x/sqrt(2.)));
}

/*************************************************************************
Student's t distribution

Computes the integral from minus infinity to t of the Student
t distribution with integer k > 0 degrees of freedom:

                                     t
                                     -
                                    | |
             -                      |         2   -(k+1)/2
            | ( (k+1)/2 )           |  (     x   )
      ----------------------        |  ( 1 + --- )        dx
                    -               |  (      k  )
      sqrt( k pi ) | ( k/2 )        |
                                  | |
                                   -
                                  -inf.

Relation to incomplete beta integral:

       1 - stdtr(k,t) = 0.5 * incbet( k/2, 1/2, z )
where
       z = k/(k + t**2).

For t < -2, this is the method of computation.  For higher t,
a direct method is derived from integration by parts.
Since the function is symmetric about t=0, the area under the
right tail of the density is found by calling the function
with -t instead of t.

ACCURACY:

Tested at random 1 <= k <= 25.  The "domain" refers to t.
                     Relative error:
arithmetic   domain     # trials      peak         rms
   IEEE     -100,-2      50000       5.9e-15     1.4e-15
   IEEE     -2,100      500000       2.7e-15     4.9e-17

Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1995, 2000 by Stephen L. Moshier
*************************************************************************/
RFLOAT cdf_tstudent(int k, RFLOAT t)
{
    RFLOAT EPS=5E-16;
    RFLOAT result;
    RFLOAT x;
    RFLOAT rk;
    RFLOAT z;
    RFLOAT f;
    RFLOAT tz;
    RFLOAT p;
    RFLOAT xsqk;
    int j;

    if ( t==0 )
    {
        result = 0.5;
        return result;
    }
    if ( t<-2.0 )
    {
        rk = k;
        z = rk/(rk+t*t);
        result = 0.5*betai(0.5*rk, 0.5, z);
        return result;
    }
    if ( t<0 )
    {
        x = -t;
    }
    else
    {
        x = t;
    }
    rk = k;
    z = 1.0+x*x/rk;
    if ( k%2 != 0 )
    {
        xsqk = x/sqrt(rk);
        p = atan(xsqk);
        if ( k > 1 )
        {
            f = 1.0;
            tz = 1.0;
            j = 3;
            while ( j <= k-2 && tz/f > EPS )
            {
                tz = tz*((j-1)/(z*j));
                f = f+tz;
                j = j+2;
            }
            p = p+f*xsqk/z;
        }
        p = p*2.0/PI;
    }
    else
    {
        f = 1.0;
        tz = 1.0;
        j = 2;
        while ( j<= k-2 && tz/f > EPS)
        {
            tz = tz*((j-1)/(z*j));
            f = f+tz;
            j = j+2;
        }
        p = f*x/sqrt(z*rk);
    }
    if ( t<0 )
    {
        p = -p;
    }
    result = 0.5+0.5*p;
    return result;
}

/* Snedecor's F ------------------------------------------------------------ */
// http://en.wikipedia.org/wiki/F-distribution
RFLOAT cdf_FSnedecor(int d1, int d2, RFLOAT x)
{
    return betai(0.5*d1,0.5*d2,(d1*x)/(d1*x+d2));
}

RFLOAT icdf_FSnedecor(int d1, int d2, RFLOAT p)
{
    RFLOAT xl=0, xr=1e6;
    RFLOAT pl=cdf_FSnedecor(d1,d2,xl);
    RFLOAT pr=cdf_FSnedecor(d1,d2,xr);
    RFLOAT xm, pm;
    do
    {
        xm=(xl+xr)*0.5;
        pm=cdf_FSnedecor(d1,d2,xm);
        if (pm>p)
        {
            xr=xm;
            pr=pm;
        }
        else
        {
            xl=xm;
            pl=pm;
        }
    }
    while (ABS(pm-p)/p>0.001);
    return xm;
}

// Uniform distribution ....................................................
void init_random_generator(int seed)
{
    if (seed < 0)
    	randomize_random_generator();
    else
    	srand(static_cast <unsigned> (seed) );
}

void randomize_random_generator()
{
	srand(static_cast <unsigned> (time(NULL)) );
}

float rnd_unif(float a, float b)
{


	if (a == b)
        return a;
    else
        return a + static_cast <float> (rand()) /( static_cast <float> (RAND_MAX/(b-a)));
}

// Gaussian distribution ...................................................
float rnd_gaus(float mu, float sigma)
{
  float U1, U2, W, mult;
  static float X1, X2;
  static int call = 0;

  if (sigma == 0)
	  return mu;

  if (call == 1)
  {
      call = !call;
      return (mu + sigma * (float) X2);
  }

  do
  {
      U1 = -1 + ((float) rand () / RAND_MAX) * 2;
      U2 = -1 + ((float) rand () / RAND_MAX) * 2;
      W = pow (U1, 2) + pow (U2, 2);
  }
  while (W >= 1 || W == 0);

  mult = sqrt ((-2 * log (W)) / W);
  X1 = U1 * mult;
  X2 = U2 * mult;

  call = !call;

  return (mu + sigma * (float) X1);

}

float rnd_student_t(RFLOAT nu, float mu, float sigma)
{
	REPORT_ERROR("rnd_student_t currently not implemented!");
}


float gaus_within_x0(float x0, float mean, float stddev)
{
    float z0 = (x0 - mean) / stddev;
    return erf(ABS(z0) / sqrt(2.0));
}

float gaus_outside_x0(float x0, float mean, float stddev)
{
    float z0 = (x0 - mean) / stddev;
    return erfc(ABS(z0) / sqrt(2.0));
}

float gaus_up_to_x0(float x0, float mean, float stddev)
{
    if (x0 > mean)
        return 1.0 -gaus_outside_x0(x0, mean, stddev) / 2;
    else if (x0 == mean)
        return 0.5;
    else
        return gaus_outside_x0(x0, mean, stddev) / 2;
}

float gaus_from_x0(float x0, float mean, float stddev)
{
    if (x0 > mean)
        return gaus_outside_x0(x0, mean, stddev) / 2;
    else if (x0 == mean)
        return 0.5;
    else
        return 1.0 -gaus_outside_x0(x0, mean, stddev) / 2;
}

float gaus_outside_probb(float p, float mean, float stddev)
{
    // Make a Bolzano search for the right value
    float p1, p2, pm, x1, x2, xm;
    x1 = mean;
    x2 = mean + 5 * stddev;
    do
    {
        xm = (x1 + x2) / 2;
        p1 = gaus_outside_x0(x1, mean, stddev);
        p2 = gaus_outside_x0(x2, mean, stddev);
        pm = gaus_outside_x0(xm, mean, stddev);
        if (pm > p)
            x1 = xm;
        else
            x2 = xm;
    }
    while (ABS(pm - p) / p > 0.005);
    return xm;
}

// See Numerical Recipes, Chap. 6.3
float student_within_t0(float t0, float degrees_of_freedom)
{
    return 1 -betai(degrees_of_freedom / 2, 0.5,
                    degrees_of_freedom / (degrees_of_freedom + t0*t0));
}

float student_outside_t0(float t0, float degrees_of_freedom)
{
    return 1 -student_within_t0(t0, degrees_of_freedom);
}

float student_up_to_t0(float t0, float degrees_of_freedom)
{
    if (t0 >= 0)
        return 1.0 -student_outside_t0(t0, degrees_of_freedom) / 2;
    else
        return student_outside_t0(t0, degrees_of_freedom) / 2;
}

float student_from_t0(float t0, float degrees_of_freedom)
{
    return 1 -student_up_to_t0(t0, degrees_of_freedom);
}

float student_outside_probb(float p, float degrees_of_freedom)
{
    // Make a Bolzano search for the right value
    float p1, p2, pm, t1, t2, tm;
    t1 = 0;
    t2 = 100;
    do
    {
        tm = (t1 + t2) / 2;
        p1 = student_outside_t0(t1, degrees_of_freedom);
        p2 = student_outside_t0(t2, degrees_of_freedom);
        pm = student_outside_t0(tm, degrees_of_freedom);
        if (pm > p)
            t1 = tm;
        else
            t2 = tm;
    }
    while (ABS(pm - p) / p > 0.005);
    return tm;
}

float chi2_up_to_t0(float t0, float degrees_of_freedom)
{
    return gammp(degrees_of_freedom / 2, t0 / 2);
}

float chi2_from_t0(float t0, float degrees_of_freedom)
{
    return 1 -chi2_up_to_t0(t0, degrees_of_freedom);
}

// Log uniform distribution ................................................
float rnd_log(float a, float b)
{
    if (a == b)
        return a;
    else
        return exp(rnd_unif(log(a), log(b)));
}

// Bsoft function
void swapbytes(char* v, unsigned long n)
{
    char            t;
    for ( int i=0; i<n/2; i++ )
    {
        t = v[i];
        v[i] = v[n-1-i];
        v[n-1-i] = t;
    }
}

void HSL2RGB(RFLOAT H, RFLOAT S, RFLOAT L, RFLOAT &R, RFLOAT &G, RFLOAT &B)
{
	if (S < XMIPP_EQUAL_ACCURACY)
	{
		R = G = B = L;
	}
	else
	{

		RFLOAT temp1 = (L < 0.5) ? L * (1.0 + S) : L + S - L * S;
		RFLOAT temp2 = 2 * L - temp1;
		RFLOAT tR = H + 0.33333;
		RFLOAT tG = H;
		RFLOAT tB = H - 0.33333;
		realWRAP(tR, 0., 1.);
		realWRAP(tG, 0., 1.);
		realWRAP(tB, 0., 1.);

		// Red
		if (6*tR < 1.)
			R = temp2 + (temp1 - temp2) * 6 * tR;
		else if (2*tR < 1.)
			R = temp1;
		else if (3*tR < 2.)
			R = temp2 + (temp1 - temp2) * (0.6666 - tR) * 6;
		else
			R = temp2;

		// Green
		if (6*tG < 1.)
			G = temp2 + (temp1 - temp2) * 6 * tG;
		else if (2*tG < 1.)
			G = temp1;
		else if (3*tG < 2.)
			G = temp2 + (temp1 - temp2) * (0.6666 - tG) * 6;
		else
			G = temp2;

		// Blue
		if (6*tB < 1.)
			B = temp2 + (temp1 - temp2) * 6 * tB;
		else if (2*tB < 1.)
			B = temp1;
		else if (3*tB < 2.)
			B = temp2 + (temp1 - temp2) * (0.6666 - tB) * 6;
		else
			B = temp2;

	}
}



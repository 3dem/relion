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

void fitStraightLine(std::vector<fit_point2D> &points, double &slope, double &intercept, double &corr_coeff)
{
	// From: http://mathworld.wolfram.com/LeastSquaresFitting.html
	// ss_xx = Sum_i x_i^2 - n ave_x^2
	// ss_yy = Sum_i y_i^2 - n ave_y^2
	// ss_xy = Sum_i x_i * y_i - n ave_x n_ave_y
	// slope = xx_xy / ss_xx
	// intercept = ave_y - slope * ave_x
	// corr_coeff = ss_xy^2 / (ss_xx * ss_yy)
	double ss_xy = 0.;
	double ss_xx = 0.;
	double ss_yy = 0.;
	double ave_x = 0.;
	double ave_y = 0.;
	double sum_w = 0.;
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

/* Value of a blob --------------------------------------------------------- */
double kaiser_value(double r, double a, double alpha, int m)
{
    double rda, rdas, arg, w;
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
double kaiser_proj(double s, double a, double alpha, int m)
{
    double sda, sdas, w, arg, p;
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
double kaiser_Fourier_value(double w, double a, double alpha, int m)
{
    double sigma = sqrt(ABS(alpha * alpha - (2. * PI * a * w) * (2. * PI * a * w)));
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
double  basvolume(double a, double alpha, int m, int n)
{
    double  hn, tpi, v;
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
    return v * pow(a, (double)n);
}
/* Bessel function I_n (x),  n = 0, 1, 2, ...
 Use ONLY for small values of n     */
double i_n(int n, double x)
{
    int i;
    double i_ns1, i_n, i_np1;
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
double i_nph(int n, double x)
{
    int i;
    double r2dpix;
    double i_ns1, i_n, i_np1;
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
double in_zeroarg(int n)
{
    int i;
    double fact;
    fact = 1.0;
    for (i = 1; i <= n; i++)
    {
        fact *= 0.5 / i;
    }
    return fact;
}
/*.......Limit (z->0) of (1/z)^(n+1/2) I_(n+1/2) (z)..........................*/
double inph_zeroarg(int n)
{
    int i;
    double fact;
    fact = 1.0;
    for (i = 1; i <= n; i++)
    {
        fact *= 1.0 / (2 * i + 1.0);
    }
    return fact*sqrt(2.0 / PI);
}
/* Zero freq --------------------------------------------------------------- */
double blob_freq_zero(struct blobtype b)
{
    return sqrt(b.alpha*b.alpha + 6.9879*6.9879) / (2*PI*b.radius);
}
/* Attenuation ------------------------------------------------------------- */
double blob_att(double w, struct blobtype b)
{
    return blob_Fourier_val(w, b) / blob_Fourier_val(0, b);
}
/* Number of operations ---------------------------------------------------- */
double blob_ops(double w, struct blobtype b)
{
    return pow(b.alpha*b.alpha + 6.9879*6.9879, 1.5) / b.radius;
}

/* Gaussian value ---------------------------------------------------------- */
double gaussian1D(double x, double sigma, double mu)
{
    x -= mu;
    return 1 / sqrt(2*PI*sigma*sigma)*exp(-0.5*((x / sigma)*(x / sigma)));
}

/* t-student value -------------------------------------------------------- */
double tstudent1D(double x, double df, double sigma, double mu)
{
    x -= mu;
    double norm = exp(gammln((df+1.)/2.)) / exp(gammln(df/2.));
    norm /= sqrt(df*PI*sigma*sigma);
    return norm * pow((1 + (x/sigma)*(x/sigma)/df),-((df+1.)/2.));

}

double gaussian2D(double x, double y, double sigmaX, double sigmaY,
                  double ang, double muX, double muY)
{
    // Express x,y in the gaussian internal coordinates
    x -= muX;
    y -= muY;
    double xp = cos(ang) * x + sin(ang) * y;
    double yp = -sin(ang) * x + cos(ang) * y;

    // Now evaluate
    return 1 / sqrt(2*PI*sigmaX*sigmaY)*exp(-0.5*((xp / sigmaX)*(xp / sigmaX) +
                                            (yp / sigmaY)*(yp / sigmaY)));
}

/* ICDF Gaussian ----------------------------------------------------------- */
double icdf_gauss(double p)
{
    const double c[] =
        {
            2.515517, 0.802853, 0.010328
        };
    const double d[] =
        {
            1.432788, 0.189269, 0.001308
        };
    if (p < 0.5)
    {
        // F^-1(p) = - G^-1(p)
        double t=sqrt(-2.0*log(p));
        double z=t - ((c[2]*t + c[1])*t + c[0]) /
                 (((d[2]*t + d[1])*t + d[0])*t + 1.0);
        return -z;
    }
    else
    {
        // F^-1(p) = G^-1(1-p)
        double t=sqrt(-2.0*log(1-p));
        double z=t - ((c[2]*t + c[1])*t + c[0]) /
                 (((d[2]*t + d[1])*t + d[0])*t + 1.0);
        return z;
    }
}

/* CDF Gaussian ------------------------------------------------------------ */
double cdf_gauss(double x)
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
double cdf_tstudent(int k, double t)
{
    double EPS=5E-16;
    double result;
    double x;
    double rk;
    double z;
    double f;
    double tz;
    double p;
    double xsqk;
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
double cdf_FSnedecor(int d1, int d2, double x)
{
    return betai(0.5*d1,0.5*d2,(d1*x)/(d1*x+d2));
}

double icdf_FSnedecor(int d1, int d2, double p)
{
    double xl=0, xr=1e6;
    double pl=cdf_FSnedecor(d1,d2,xl);
    double pr=cdf_FSnedecor(d1,d2,xr);
    double xm, pm;
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

/* Random functions -------------------------------------------------------- */
int idum;

// Uniform distribution ....................................................
void  init_random_generator(int seed)
{
    idum = -1;
    ran1(&idum);
    if (seed != -1)
    {
        // Prevent seeds larger than 65000
        seed %=0xffff;
        for (int i = 0; i < seed; i++)
            ran1(&idum);
    }
}

void  randomize_random_generator()
{
    static  unsigned int seed;
    int rand_return;

    srand(seed);
    rand_return = rand();

    time_t t;
    time(&t);
    rand_return = abs(rand_return);
    idum = (-(int)(t % 10000)
            - (int)(rand_return % 10000));
    ran1(&idum);
    seed = (unsigned int)rand_return;
}

float rnd_unif()
{
    return ran1(&idum);
}
float rnd_unif(float a, float b)
{
    if (a == b)
        return a;
    else
        return a + (b - a)*ran1(&idum);
}

// t-distribution
float rnd_student_t(double nu)
{
    return tdev(nu, &idum);
}
float rnd_student_t(double nu, float a, float b)
{
    if (b == 0)
        return a;
    else
        return b*tdev(nu, &idum) + a;
}

// Gaussian distribution ...................................................
float rnd_gaus()
{
    return gasdev(&idum);
}
float rnd_gaus(float a, float b)
{
    if (b == 0)
        return a;
    else
        return b*gasdev(&idum) + a;
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

/* Log2 -------------------------------------------------------------------- */
// Does not work with xlc compiler
//#ifndef __xlC__
//double log2(double value)
//{
//    return 3.32192809488736*log10(value);
//    // log10(value)/log10(2)
//}
//#endif

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




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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "src/numerical_recipes.h"

/* NUMERICAL UTILITIES ----------------------------------------------------- */
void nrerror(const char error_text[])
{
    fprintf(stderr, "Numerical Recipes run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}
#define NRSIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

/* RANDOM NUMBERS ---------------------------------------------------------- */
#define M1 259200
#define IA1 7141
#define IC1 54773
#define RM1 (1.0/M1)
#define M2 134456
#define IA2 8121
#define IC2 28411
#define RM2 (1.0/M2)
#define M3 243000
#define IA3 4561
#define IC3 51349

/* Chapter 7 Section 1: UNIFORM RANDOM NUMBERS */
double ran1(int *idum)
{
    static long ix1, ix2, ix3;
    static double r[98];
    double temp;
    static int iff = 0;
    int j;

    if (*idum < 0 || iff == 0)
    {
        iff = 1;
        ix1 = (IC1 - (*idum)) % M1;
        ix1 = (IA1 * ix1 + IC1) % M1;
        ix2 = ix1 % M2;
        ix1 = (IA1 * ix1 + IC1) % M1;
        ix3 = ix1 % M3;
        for (j = 1;j <= 97;j++)
        {
            ix1 = (IA1 * ix1 + IC1) % M1;
            ix2 = (IA2 * ix2 + IC2) % M2;
            r[j] = (ix1 + ix2 * RM2) * RM1;
        }
        *idum = 1;
    }
    ix1 = (IA1 * ix1 + IC1) % M1;
    ix2 = (IA2 * ix2 + IC2) % M2;
    ix3 = (IA3 * ix3 + IC3) % M3;
    j = 1 + ((97 * ix3) / M3);
    if (j > 97 || j < 1)
        nrerror("RAN1: This cannot happen.");
    temp = r[j];
    r[j] = (ix1 + ix2 * RM2) * RM1;
    return temp;
}

#undef M1
#undef IA1
#undef IC1
#undef RM1
#undef M2
#undef IA2
#undef IC2
#undef RM2
#undef M3
#undef IA3
#undef IC3

/* Chapter 7 Section 3: GAUSSIAN RANDOM NUMBERS */
double gasdev(int *idum)
{
    static int iset = 0;
    static double gset;
    double fac, r, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            r = v1 * v1 + v2 * v2;
        }
        while (r >= 1.0);
        fac = sqrt(-2.0 * log(r) / r);
        gset = v1 * fac;
        iset = 1;
        return v2*fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}

// t-distribution (nor Numerical Recipes, but Mathematics of Computation, vol. 62, 779-781.
// I downloaded sem-code from http://ideas.repec.org/c/ega/comcod/200703.html
// Sjors May 2008
double tdev(double nu, int *idum)
{
    static int iset = 0;
    static double gset;
    double fac, r, v1, v2;

    if (iset == 0)
    {
        do
        {
            v1 = 2.0 * ran1(idum) - 1.0;
            v2 = 2.0 * ran1(idum) - 1.0;
            r = v1 * v1 + v2 * v2;
        }
        while (r >= 1.0);
        fac = sqrt(nu*(pow(r,-2.0/nu) -1.0)/r);
        gset = v1 * fac;
        iset = 1;
        return v2*fac;
    }
    else
    {
        iset = 0;
        return gset;
    }
}


/* BESSEL FUNCTIONS -------------------------------------------------------- */
/* CO: They may not come in the numerical recipes but it is not a bad
   idea to put them here, in fact they come from Gabor's group in Feb'84     */
double bessj0(double x)
{
    double ax, z;
    double xx, y, ans, ans1, ans2;

    if ((ax = fabs(x)) < 8.0)
    {
        y = x * x;
        ans1 = 57568490574.0 + y * (-13362590354.0 +
                                    y * (651619640.7
                                         + y * (-11214424.18 +
                                                y * (77392.33017 +
                                                     y * (-184.9052456)))));
        ans2 = 57568490411.0 + y * (1029532985.0 +
                                    y * (9494680.718
                                         + y * (59272.64853 +
                                                y * (267.8532712 +
                                                     y * 1.0))));
        ans = ans1 / ans2;
    }
    else
    {
        z = 8.0 / ax;
        y = z * z;
        xx = ax - 0.785398164;
        ans1 = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4
                          + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
        ans2 = -0.1562499995e-1 + y * (0.1430488765e-3
                                       + y * (-0.6911147651e-5 + y * (0.7621095161e-6
                                                                      - y * 0.934935152e-7)));
        ans = sqrt(0.636619772 / ax) * (cos(xx) * ans1 - z * sin(xx) * ans2);
    }
    return ans;
}

/*............................................................................*/
double bessi0(double x)
{
    double y, ax, ans;
    if ((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        ans = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                                          + y * (0.2659732 + y * (0.360768e-1 + y * 0.45813e-2)))));
    }
    else
    {
        y = 3.75 / ax;
        ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1
                                      + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2
                                                                + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1
                                                                                            + y * 0.392377e-2))))))));
    }
    return ans;
}

/*............................................................................*/
double bessi1(double x)
{
    double ax, ans;
    double y;
    if ((ax = fabs(x)) < 3.75)
    {
        y = x / 3.75;
        y *= y;
        ans = ax * (0.5 + y * (0.87890594 + y * (0.51498869 + y * (0.15084934
                               + y * (0.2658733e-1 + y * (0.301532e-2 + y * 0.32411e-3))))));
    }
    else
    {
        y = 3.75 / ax;
        ans = 0.2282967e-1 + y * (-0.2895312e-1 + y * (0.1787654e-1
                                  - y * 0.420059e-2));
        ans = 0.39894228 + y * (-0.3988024e-1 + y * (-0.362018e-2
                                + y * (0.163801e-2 + y * (-0.1031555e-1 + y * ans))));
        ans *= (exp(ax) / sqrt(ax));
    }
    return x < 0.0 ? -ans : ans;
}

/* General Bessel functions ------------------------------------------------ */
double chebev(double a, double b, double c[], int m, double x)
{
    double d = 0.0, dd = 0.0, sv, y, y2;
    int j;

    if ((x - a)*(x - b) > 0.0)
        nrerror("x not in range in routine chebev");
    y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a));
    for (j = m - 1;j >= 1;j--)
    {
        sv = d;
        d = y2 * d - dd + c[j];
        dd = sv;
    }
    return y*d - dd + 0.5*c[0];
}
#define NUSE1 5
#define NUSE2 5

void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
    double xx;
    static double c1[] =
        {
            -1.142022680371172e0, 6.516511267076e-3,
            3.08709017308e-4, -3.470626964e-6, 6.943764e-9,
            3.6780e-11, -1.36e-13
        };
    static double c2[] =
        {
            1.843740587300906e0, -0.076852840844786e0,
            1.271927136655e-3, -4.971736704e-6, -3.3126120e-8,
            2.42310e-10, -1.70e-13, -1.0e-15
        };

    xx = 8.0 * x * x - 1.0;
    *gam1 = chebev(-1.0, 1.0, c1, NUSE1, xx);
    *gam2 = chebev(-1.0, 1.0, c2, NUSE2, xx);
    *gampl = *gam2 - x * (*gam1);
    *gammi = *gam2 + x * (*gam1);
}

#undef NUSE1
#undef NUSE2

#define EPS 1.0e-16
#define FPMIN 1.0e-30
#define MAXIT 10000
#define XMIN 2.0
void bessjy(double x, double xnu, double *rj, double *ry, double *rjp, double *ryp)
{
    int i, isign, l, nl;
    double a, b, br, bi, c, cr, ci, d, del, del1, den, di, dlr, dli, dr, e, f, fact, fact2,
    fact3, ff, gam, gam1, gam2, gammi, gampl, h, p, pimu, pimu2, q, r, rjl,
    rjl1, rjmu, rjp1, rjpl, rjtemp, ry1, rymu, rymup, rytemp, sum, sum1,
    temp, w, x2, xi, xi2, xmu, xmu2;

    if (x <= 0.0 || xnu < 0.0)
        nrerror("bad arguments in bessjy");
    nl = (x < XMIN ? (int)(xnu + 0.5) : XMIPP_MAX(0, (int)(xnu - x + 1.5)));
    xmu = xnu - nl;
    xmu2 = xmu * xmu;
    xi = 1.0 / x;
    xi2 = 2.0 * xi;
    w = xi2 / PI;
    isign = 1;
    h = xnu * xi;
    if (h < FPMIN)
        h = FPMIN;
    b = xi2 * xnu;
    d = 0.0;
    c = h;
    for (i = 1;i <= MAXIT;i++)
    {
        b += xi2;
        d = b - d;
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = b - 1.0 / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        del = c * d;
        h = del * h;
        if (d < 0.0)
            isign = -isign;
        if (fabs(del - 1.0) < EPS)
            break;
    }
    if (i > MAXIT)
        nrerror("x too large in bessjy; try asymptotic expansion");
    rjl = isign * FPMIN;
    rjpl = h * rjl;
    rjl1 = rjl;
    rjp1 = rjpl;
    fact = xnu * xi;
    for (l = nl;l >= 1;l--)
    {
        rjtemp = fact * rjl + rjpl;
        fact -= xi;
        rjpl = fact * rjtemp - rjl;
        rjl = rjtemp;
    }
    if (rjl == 0.0)
        rjl = EPS;
    f = rjpl / rjl;
    if (x < XMIN)
    {
        x2 = 0.5 * x;
        pimu = PI * xmu;
        fact = (fabs(pimu) < EPS ? 1.0 : pimu / sin(pimu));
        d = -log(x2);
        e = xmu * d;
        fact2 = (fabs(e) < EPS ? 1.0 : sinh(e) / e);
        beschb(xmu, &gam1, &gam2, &gampl, &gammi);
        ff = 2.0 / PI * fact * (gam1 * cosh(e) + gam2 * fact2 * d);
        e = exp(e);
        p = e / (gampl * PI);
        q = 1.0 / (e * PI * gammi);
        pimu2 = 0.5 * pimu;
        fact3 = (fabs(pimu2) < EPS ? 1.0 : sin(pimu2) / pimu2);
        r = PI * pimu2 * fact3 * fact3;
        c = 1.0;
        d = -x2 * x2;
        sum = ff + r * q;
        sum1 = p;
        for (i = 1;i <= MAXIT;i++)
        {
            ff = (i * ff + p + q) / (i * i - xmu2);
            c *= (d / i);
            p /= (i - xmu);
            q /= (i + xmu);
            del = c * (ff + r * q);
            sum += del;
            del1 = c * p - i * del;
            sum1 += del1;
            if (fabs(del) < (1.0 + fabs(sum))*EPS)
                break;
        }
        if (i > MAXIT)
            nrerror("bessy series failed to converge");
        rymu = -sum;
        ry1 = -sum1 * xi2;
        rymup = xmu * xi * rymu - ry1;
        rjmu = w / (rymup - f * rymu);
    }
    else
    {
        a = 0.25 - xmu2;
        p = -0.5 * xi;
        q = 1.0;
        br = 2.0 * x;
        bi = 2.0;
        fact = a * xi / (p * p + q * q);
        cr = br + q * fact;
        ci = bi + p * fact;
        den = br * br + bi * bi;
        dr = br / den;
        di = -bi / den;
        dlr = cr * dr - ci * di;
        dli = cr * di + ci * dr;
        temp = p * dlr - q * dli;
        q = p * dli + q * dlr;
        p = temp;
        for (i = 2;i <= MAXIT;i++)
        {
            a += 2 * (i - 1);
            bi += 2.0;
            dr = a * dr + br;
            di = a * di + bi;
            if (fabs(dr) + fabs(di) < FPMIN)
                dr = FPMIN;
            fact = a / (cr * cr + ci * ci);
            cr = br + cr * fact;
            ci = bi - ci * fact;
            if (fabs(cr) + fabs(ci) < FPMIN)
                cr = FPMIN;
            den = dr * dr + di * di;
            dr /= den;
            di /= -den;
            dlr = cr * dr - ci * di;
            dli = cr * di + ci * dr;
            temp = p * dlr - q * dli;
            q = p * dli + q * dlr;
            p = temp;
            if (fabs(dlr - 1.0) + fabs(dli) < EPS)
                break;
        }
        if (i > MAXIT)
            nrerror("cf2 failed in bessjy");
        gam = (p - f) / q;
        rjmu = sqrt(w / ((p - f) * gam + q));
        rjmu = NRSIGN(rjmu, rjl);
        rymu = rjmu * gam;
        rymup = rymu * (p + q / gam);
        ry1 = xmu * xi * rymu - rymup;
    }
    fact = rjmu / rjl;
    *rj = rjl1 * fact;
    *rjp = rjp1 * fact;
    for (i = 1;i <= nl;i++)
    {
        rytemp = (xmu + i) * xi2 * ry1 - rymu;
        rymu = ry1;
        ry1 = rytemp;
    }
    *ry = rymu;
    *ryp = xnu * xi * rymu - ry1;
}
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN

/*............................................................................*/
double bessi0_5(double x)
{
    return (x == 0) ? 0 : sqrt(2 / (PI*x))*sinh(x);
}
double bessi1_5(double x)
{
    return (x == 0) ? 0 : sqrt(2 / (PI*x))*(cosh(x) - sinh(x) / x);
}
double bessi2(double x)
{
    return (x == 0) ? 0 : bessi0(x) - ((2*1) / x) * bessi1(x);
}
double bessi2_5(double x)
{
    return (x == 0) ? 0 : bessi0_5(x) - ((2*1.5) / x) * bessi1_5(x);
}
double bessi3(double x)
{
    return (x == 0) ? 0 : bessi1(x) - ((2*2) / x) * bessi2(x);
}
double bessi3_5(double x)
{
    return (x == 0) ? 0 : bessi1_5(x) - ((2*2.5) / x) * bessi2_5(x);
}
double bessi4(double x)
{
    return (x == 0) ? 0 : bessi2(x) - ((2*3) / x) * bessi3(x);
}
double bessj1_5(double x)
{
    double rj, ry, rjp, ryp;
    bessjy(x, 1.5, &rj, &ry, &rjp, &ryp);
    return rj;
}
double bessj3_5(double x)
{
    double rj, ry, rjp, ryp;
    bessjy(x, 3.5, &rj, &ry, &rjp, &ryp);
    return rj;
}

/* Special functions ------------------------------------------------------- */
double gammln(double xx)
{
    double x, tmp, ser;
    static double cof[6] =
        {
            76.18009173, -86.50532033, 24.01409822,
            -1.231739516, 0.120858003e-2, -0.536382e-5
        };
    int j;

    x = xx - 1.0;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.0;
    for (j = 0;j <= 5;j++)
    {
        x += 1.0;
        ser += cof[j] / x;
    }
    return -tmp + log(2.50662827465*ser);
}


double betai(double a, double b, double x)
{
    double bt;
    if (x < 0.0 || x > 1.0)
        nrerror("Bad x in routine BETAI");
    if (x == 0.0 || x == 1.0)
        bt = 0.0;
    else
        bt = exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1.0 - x));
    if (x < (a + 1.0) / (a + b + 2.0))
        return bt*betacf(a, b, x) / a;
    else
        return 1.0 -bt*betacf(b, a, 1.0 - x) / b;

}

#define ITMAX 100
#define EPS 3.0e-7
double betacf(double a, double b, double x)
{
    double qap, qam, qab, em, tem, d;
    double bz, bm = 1.0, bp, bpp;
    double az = 1.0, am = 1.0, ap, app, aold;
    int m;

    qab = a + b;
    qap = a + 1.0;
    qam = a - 1.0;
    bz = 1.0 - qab * x / qap;
    for (m = 1;m <= ITMAX;m++)
    {
        em = (double) m;
        tem = em + em;
        d = em * (b - em) * x / ((qam + tem) * (a + tem));
        ap = az + d * am;
        bp = bz + d * bm;
        d = -(a + em) * (qab + em) * x / ((qap + tem) * (a + tem));
        app = ap + d * az;
        bpp = bp + d * bz;
        aold = az;
        am = ap / bpp;
        bm = bp / bpp;
        az = app / bpp;
        bz = 1.0;
        if (fabs(az - aold) < (EPS*fabs(az)))
            return az;
    }
    nrerror("a or b too big, or ITMAX too small in BETACF");
    return 0;
}
#undef ITMAX
#undef EPS

/* Powell optimization ------------------------------------------------------------ */
#undef MAX
#undef SIGN
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define F1DIM(x,f) {\
    for (int j = 1; j<=ncom; j++) \
        xt[j] = pcom[j] + x * xicom[j]; \
    f = (*func)(xt,prm);}

void mnbrak(double *ax, double *bx, double *cx,
            double *fa, double *fb, double *fc, double(*func)(double *, void*),
            void *prm, int ncom, double *pcom, double *xicom)
{
    double ulim, u, r, q, fu, dum;
    double *xt=NULL;
    ask_Tvector(xt, 1, ncom);

    F1DIM(*ax,*fa);
    F1DIM(*bx,*fb);
    if (*fb > *fa)
    {
        SHFT(dum, *ax, *bx, dum)
        SHFT(dum, *fb, *fa, dum)
    }
    *cx = (*bx) + GOLD * (*bx - *ax);
    F1DIM(*cx,*fc);
    while (*fb > *fc)
    {
        r = (*bx - *ax) * (*fb - *fc);
        q = (*bx - *cx) * (*fb - *fa);
        u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
            (2.0 * SIGN(MAX(fabs(q - r), TINY), q - r));
        ulim = (*bx) + GLIMIT * (*cx - *bx);
        if ((*bx - u)*(u - *cx) > 0.0)
        {
            F1DIM(u,fu);
            if (fu < *fc)
            {
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                return;
            }
            else if (fu > *fb)
            {
                *cx = u;
                *fc = fu;
                return;
            }
            u = (*cx) + GOLD * (*cx - *bx);
            F1DIM(u,fu);
        }
        else if ((*cx - u)*(u - ulim) > 0.0)
        {
            F1DIM(u,fu);
            if (fu < *fc)
            {
                SHFT(*bx, *cx, u, *cx + GOLD*(*cx - *bx))
                double aux; F1DIM(u,aux);
                SHFT(*fb, *fc, fu, aux)
            }
        }
        else if ((u - ulim)*(ulim - *cx) >= 0.0)
        {
            u = ulim;
            F1DIM(u,fu);
        }
        else
        {
            u = (*cx) + GOLD * (*cx - *bx);
            F1DIM(u,fu);
        }
        SHFT(*ax, *bx, *cx, u)
        SHFT(*fa, *fb, *fc, fu)
    }
    free_Tvector(xt, 1, ncom);
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
double brent(double ax, double bx, double cx, double(*func)(double *,void*),
             void *prm, double tol, double *xmin,
             int ncom, double *pcom, double *xicom)
{
    int iter;
    double a, b, d, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
    double e = 0.0;
    double *xt=NULL;
    ask_Tvector(xt, 1, ncom);

    a = (ax < cx ? ax : cx);
    b = (ax > cx ? ax : cx);
    x = w = v = bx;
    F1DIM(x,fx);
    fw = fv = fx;
    for (iter = 1;iter <= ITMAX;iter++)
    {
        xm = 0.5 * (a + b);
        tol2 = 2.0 * (tol1 = tol * fabs(x) + ZEPS);
        if (fabs(x - xm) <= (tol2 - 0.5*(b - a)))
        {
            *xmin = x;
            free_Tvector(xt, 1, ncom);
            return fx;
        }
        if (fabs(e) > tol1)
        {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0)
                p = -p;
            q = fabs(q);
            etemp = e;
            e = d;
            if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a - x) || p >= q*(b - x))
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            else
            {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = SIGN(tol1, xm - x);
            }
        }
        else
        {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }
        u = (fabs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
        F1DIM(u,fu);
        if (fu <= fx)
        {
            if (u >= x)
                a = x;
            else
                b = x;
            SHFT(v, w, x, u)
            SHFT(fv, fw, fx, fu)
        }
        else
        {
            if (u < x)
                a = u;
            else
                b = u;
            if (fu <= fw || w == x)
            {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            }
            else if (fu <= fv || v == x || v == w)
            {
                v = u;
                fv = fu;
            }
        }
    }
    nrerror("Too many iterations in brent");
    *xmin = x;
    free_Tvector(xt, 1, ncom);
    return fx;
}
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef F1DIM

#define TOL 2.0e-4
void linmin(double *p, double *xi, int n, double &fret,
            double(*func)(double *, void*), void *prm)
{
    int j;
    double xx, xmin, fx, fb, fa, bx, ax;

    int ncom = n;
    double *pcom=NULL;
    double *xicom=NULL;
    ask_Tvector(pcom, 1, n);
    ask_Tvector(xicom, 1, n);
    for (j = 1;j <= n;j++)
    {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }
    ax = 0.0;
    xx = 1.0;
    bx = 2.0;
    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, func, prm, ncom, pcom, xicom);
    fret = brent(ax, xx, bx, func, prm, TOL, &xmin, ncom, pcom, xicom);
    for (j = 1;j <= n;j++)
    {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_Tvector(xicom, 1, n);
    free_Tvector(pcom, 1, n);
}
#undef TOL

#define ITMAX 200
void powell(double *p, double *xi, int n, double ftol, int &iter,
            double &fret, double(*func)(double *, void *), void *prm,
            bool show)
{
    int i, ibig, j;
    double t, fptt, fp, del;
    double *pt, *ptt, *xit;
    bool   different_from_0;

    ask_Tvector(pt, 1, n);
    ask_Tvector(ptt, 1, n);
    ask_Tvector(xit, 1, n);
    fret = (*func)(p,prm);
    for (j = 1;j <= n;j++)
        pt[j] = p[j];

    for (iter = 1;;(iter)++)
    {
        /* By coss ----- */
        if (show)
        {
            std::cout << iter << " (" << p[1];
            for (int co = 2; co <= n; co++)
                std::cout << "," << p[co];
            std::cout << ")--->" << fret << std::endl;
        }
        /* ------------- */

        fp = fret;
        ibig = 0;
        del = 0.0;
        for (i = 1;i <= n;i++)
        {
            different_from_0 = false; // CO
            for (j = 1;j <= n;j++)
            {
                xit[j] = xi[j*n+i];
                if (xit[j] != 0)
                    different_from_0 = true;
            }
            if (different_from_0)
            {
                fptt = fret;
                linmin(p, xit, n, fret, func, prm);
                if (fabs(fptt - fret) > del)
                {
                    del = fabs(fptt - fret);
                    ibig = i;
                }
                /* By coss ----- */
                if (show)
                {
                    std::cout << "   (";
                    if (i == 1)
                        std::cout << "***";
                    std::cout << p[1];
                    for (int co = 2; co <= n; co++)
                    {
                        std::cout << ",";
                        if (co == i)
                            std::cout << "***";
                        std::cout << p[co];
                    }
                    std::cout << ")--->" << fret << std::endl;
                }
                /* ------------- */
            }
        }
        if (2.0*fabs(fp - fret) <= ftol*(fabs(fp) + fabs(fret)))
        {
            free_Tvector(xit, 1, n);
            free_Tvector(ptt, 1, n);
            free_Tvector(pt, 1, n);
            return;
        }
        if (iter == ITMAX)
            nrerror("Too many iterations in routine POWELL");
        for (j = 1;j <= n;j++)
        {
            ptt[j] = 2.0 * p[j] - pt[j];
            xit[j] = p[j] - pt[j];
            pt[j] = p[j];
        }
        fptt = (*func)(ptt,prm);
        if (fptt < fp)
        {
            #define SQR(a) ((a)*(a))
            t = 2.0 * (fp - 2.0 * fret + fptt) * SQR(fp - fret - del) - del * SQR(fp - fptt);
            if (t < 0.0)
            {
                linmin(p, xit, n, fret, func, prm);
                for (j = 1;j <= n;j++)
                    xi[j*n+ibig] = xit[j];
            }
        }
    }
}
#undef ITMAX
#undef SQR

/* Singular value descomposition ------------------------------------------- */
/* Copied from Bilib library (linearalgebra.h) */
double Pythag(double a, double b)
{
    double absa, absb;
    absa = fabs(a);
    absb = fabs(b);
    if (absb < absa)
        return(absa * sqrt(1.0 + absb * absb / (absa * absa)));
    else
        return((absb == 0.0) ? (0.0) : (absb * sqrt(1.0 + absa * absa / (absb * absb))));
}

#define SVDMAXITER 1000
void svdcmp(double *U, int Lines, int Columns, double *W, double *V)
{
    double *rv1 = (double *)NULL;
    double Norm, Scale;
    double c, f, g, h, s;
    double x, y, z;
    long i, its, j, jj, k, l = 0L, nm = 0L;
    bool    Flag;
    int     MaxIterations = SVDMAXITER;

    ask_Tvector(rv1, 0, Columns*Columns - 1);
    g = Scale = Norm = 0.0;
    for (i = 0L; (i < Columns); i++)
    {
        l = i + 1L;
        rv1[i] = Scale * g;
        g = s = Scale = 0.0;
        if (i < Lines)
        {
            for (k = i; (k < Lines); k++)
            {
                Scale += fabs(U[k * Columns + i]);
            }
            if (Scale != 0.0)
            {
                for (k = i; (k < Lines); k++)
                {
                    U[k * Columns + i] /= Scale;
                    s += U[k * Columns + i] * U[k * Columns + i];
                }
                f = U[i * Columns + i];
                g = (0.0 <= f) ? (-sqrt(s)) : (sqrt(s));
                h = f * g - s;
                U[i * Columns + i] = f - g;
                for (j = l; (j < Columns); j++)
                {
                    for (s = 0.0, k = i; (k < Lines); k++)
                    {
                        s += U[k * Columns + i] * U[k * Columns + j];
                    }
                    f = s / h;
                    for (k = i; (k < Lines); k++)
                    {
                        U[k * Columns + j] += f * U[k * Columns + i];
                    }
                }
                for (k = i; (k < Lines); k++)
                {
                    U[k * Columns + i] *= Scale;
                }
            }
        }
        W[i] = Scale * g;
        g = s = Scale = 0.0;
        if ((i < Lines) && (i != (Columns - 1L)))
        {
            for (k = l; (k < Columns); k++)
            {
                Scale += fabs(U[i * Columns + k]);
            }
            if (Scale != 0.0)
            {
                for (k = l; (k < Columns); k++)
                {
                    U[i * Columns + k] /= Scale;
                    s += U[i * Columns + k] * U[i * Columns + k];
                }
                f = U[i * Columns + l];
                g = (0.0 <= f) ? (-sqrt(s)) : (sqrt(s));
                h = f * g - s;
                U[i * Columns + l] = f - g;
                for (k = l; (k < Columns); k++)
                {
                    rv1[k] = U[i * Columns + k] / h;
                }
                for (j = l; (j < Lines); j++)
                {
                    for (s = 0.0, k = l; (k < Columns); k++)
                    {
                        s += U[j * Columns + k] * U[i * Columns + k];
                    }
                    for (k = l; (k < Columns); k++)
                    {
                        U[j * Columns + k] += s * rv1[k];
                    }
                }
                for (k = l; (k < Columns); k++)
                {
                    U[i * Columns + k] *= Scale;
                }
            }
        }
        Norm = ((fabs(W[i]) + fabs(rv1[i])) < Norm) ? (Norm) : (fabs(W[i]) + fabs(rv1[i]));
    }
    for (i = Columns - 1L; (0L <= i); i--)
    {
        if (i < (Columns - 1L))
        {
            if (g != 0.0)
            {
                for (j = l; (j < Columns); j++)
                {
                    V[j * Columns + i] = U[i * Columns + j] / (U[i * Columns + l] * g);
                }
                for (j = l; (j < Columns); j++)
                {
                    for (s = 0.0, k = l; (k < Columns); k++)
                    {
                        s += U[i * Columns + k] * V[k * Columns + j];
                    }
                    for (k = l; (k < Columns); k++)
                    {
                        if (s != 0.0)
                        {
                            V[k * Columns + j] += s * V[k * Columns + i];
                        }
                    }
                }
            }
            for (j = l; (j < Columns); j++)
            {
                V[i * Columns + j] = V[j * Columns + i] = 0.0;
            }
        }
        V[i * Columns + i] = 1.0;
        g = rv1[i];
        l = i;
    }
    for (i = (Lines < Columns) ? (Lines - 1L) : (Columns - 1L); (0L <= i); i--)
    {
        l = i + 1L;
        g = W[i];
        for (j = l; (j < Columns); j++)
        {
            U[i * Columns + j] = 0.0;
        }
        if (g != 0.0)
        {
            g = 1.0 / g;
            for (j = l; (j < Columns); j++)
            {
                for (s = 0.0, k = l; (k < Lines); k++)
                {
                    s += U[k * Columns + i] * U[k * Columns + j];
                }
                f = s * g / U[i * Columns + i];
                for (k = i; (k < Lines); k++)
                {
                    if (f != 0.0)
                    {
                        U[k * Columns + j] += f * U[k * Columns + i];
                    }
                }
            }
            for (j = i; (j < Lines); j++)
            {
                U[j * Columns + i] *= g;
            }
        }
        else
        {
            for (j = i; (j < Lines); j++)
            {
                U[j * Columns + i] = 0.0;
            }
        }
        U[i * Columns + i] += 1.0;
    }
    for (k = Columns - 1L; (0L <= k); k--)
    {
        for (its = 1L; (its <= MaxIterations); its++)
        {
            Flag = true;
            for (l = k; (0L <= l); l--)
            {
                nm = l - 1L;
                if ((fabs(rv1[l]) + Norm) == Norm)
                {
                    Flag = false;
                    break;
                }
                if ((fabs(W[nm]) + Norm) == Norm)
                {
                    break;
                }
            }
            if (Flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; (i <= k); i++)
                {
                    f = s * rv1[i];
                    rv1[i] *= c;
                    if ((fabs(f) + Norm) == Norm)
                    {
                        break;
                    }
                    g = W[i];
                    h = Pythag(f, g);
                    W[i] = h;
                    h = 1.0 / h;
                    c = g * h;
                    s = -f * h;
                    for (j = 0L; (j < Lines); j++)
                    {
                        y = U[j * Columns + nm];
                        z = U[j * Columns + i];
                        U[j * Columns + nm] = y * c + z * s;
                        U[j * Columns + i] = z * c - y * s;
                    }
                }
            }
            z = W[k];
            if (l == k)
            {
                if (z < 0.0)
                {
                    W[k] = -z;
                    for (j = 0L; (j < Columns); j++)
                    {
                        V[j * Columns + k] = -V[j * Columns + k];
                    }
                }
                break;
            }
            if (its == MaxIterations)
            {
                free_Tvector(rv1, 0, Columns*Columns - 1);
                return;
            }
            x = W[l];
            nm = k - 1L;
            y = W[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = Pythag(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + ((0.0 <= f) ? (fabs(g))
                                                : (-fabs(g))))) - h)) / x;
            c = s = 1.0;
            for (j = l; (j <= nm); j++)
            {
                i = j + 1L;
                g = rv1[i];
                y = W[i];
                h = s * g;
                g = c * g;
                z = Pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (jj = 0L; (jj < Columns); jj++)
                {
                    x = V[jj * Columns + j];
                    z = V[jj * Columns + i];
                    V[jj * Columns + j] = x * c + z * s;
                    V[jj * Columns + i] = z * c - x * s;
                }
                z = Pythag(f, h);
                W[j] = z;
                if (z != 0.0)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (jj = 0L; (jj < Lines); jj++)
                {
                    y = U[jj * Columns + j];
                    z = U[jj * Columns + i];
                    U[jj * Columns + j] = y * c + z * s;
                    U[jj * Columns + i] = z * c - y * s;
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            W[k] = x;
        }
    }
    free_Tvector(rv1, 0, Columns*Columns - 1);
}

void svbksb(double *u, double *w, double *v, int m, int n, double *b, double *x)
{
    int jj, j, i;
    double s, *tmp;

    ask_Tvector(tmp, 1, n);
    for (j = 1;j <= n;j++)
    {
        s = 0.0;
        if (w[j])
        {
            for (i = 1;i <= m;i++)
                s += u[i*n+j] * b[i];
            s /= w[j];
        }
        tmp[j] = s;
    }
    for (j = 1;j <= n;j++)
    {
        s = 0.0;
        for (jj = 1;jj <= n;jj++)
            s += v[j*n+jj] * tmp[jj];
        x[j] = s;
    }
    free_Tvector(tmp, 1, n);
}


/* Gamma function ---------------------------------------------------------- */
#define ITMAX 100
#define EPS 3.0e-7

void gser(double *gamser, double a, double x, double *gln)
{
    int n;
    double sum, del, ap;

    *gln = gammln(a);
    if (x <= 0.0)
    {
        if (x < 0.0)
            nrerror("x less than 0 in routine gser");
        *gamser = 0.0;
        return;
    }
    else
    {
        ap = a;
        del = sum = 1.0 / a;
        for (n = 1;n <= ITMAX;n++)
        {
            ++ap;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum)*EPS)
            {
                *gamser = sum * exp(-x + a * log(x) - (*gln));
                return;
            }
        }
        nrerror("a too large, ITMAX too small in routine gser");
        return;
    }
}
#undef ITMAX
#undef EPS

#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
    int i;
    double an, b, c, d, del, h;

    *gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1;i <= ITMAX;i++)
    {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < FPMIN)
            d = FPMIN;
        c = b + an / c;
        if (fabs(c) < FPMIN)
            c = FPMIN;
        d = 1.0 / d;
        del = d * c;
        h *= del;
        if (fabs(del - 1.0) < EPS)
            break;
    }
    if (i > ITMAX)
        nrerror("a too large, ITMAX too small in gcf");
    *gammcf = exp(-x + a * log(x) - (*gln)) * h;
}
#undef ITMAX
#undef EPS
#undef FPMIN

double gammp(double a, double x)
{
    double gamser, gammcf, gln;

    if (x < 0.0 || a <= 0.0)
        nrerror("Invalid arguments in routine gammp");
    if (x < (a + 1.0))
    {
        gser(&gamser, a, x, &gln);
        return gamser;
    }
    else
    {
        gcf(&gammcf, a, x, &gln);
        return 1.0 -gammcf;
    }
}

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

#ifndef JAZ_INTERPOLATION_H
#define JAZ_INTERPOLATION_H

#include <src/image.h>
#include <src/fftw.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t4Matrix.h>

#define INTERPOL_WRAP(i,n) ((i) >= 0 ? (i) % (n) : -(i) % (n) ? (n) - (-(i) % (n)) : 0)

class Interpolation
{
    public:

        static bool isInSlice(const Image<RFLOAT>& img, double x, double y);
        static double getTaperWeight(const Image<RFLOAT>& img, double x, double y, double rx, double ry);

        static double linearXY(const Image<RFLOAT>& img, double x, double y, int n);
        static Complex linear3D(const Image<Complex>& img, double x, double y, double z);
        static Complex linearFFTW3D(const Image<Complex>& img, double x, double y, double z);
        static Complex linearFFTW2D(const Image<Complex>& img, double x, double y);

        template<typename T>
        static gravis::d2Vector quadraticMaxXY(const Image<T>& img, double eps = 1e-25)
        {
            const int w = img.data.xdim;
            const int h = img.data.ydim;

            int xmax = -1, ymax = -1;
            double vmax = 0.0;

            for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
            {
                T v = DIRECT_A2D_ELEM(img.data, y, x);

                if (xmax < 0 || v > vmax)
                {
                    vmax = v;
                    xmax = x;
                    ymax = y;
                }
            }

            gravis::d2Vector p(xmax, ymax);

            if (xmax > 0 && xmax < w-1)
            {
                const T vp = DIRECT_A2D_ELEM(img.data, ymax, xmax+1);
                const T vn = DIRECT_A2D_ELEM(img.data, ymax, xmax-1);

                if (std::abs(vp + vn - 2.0*vmax) > eps)
                {
                    p.x -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
                }
            }

            if (xmax > 0 && xmax < w-1)
            {
                const T vp = DIRECT_A2D_ELEM(img.data, ymax+1, xmax);
                const T vn = DIRECT_A2D_ELEM(img.data, ymax-1, xmax);

                if (std::abs(vp + vn - 2.0*vmax) > eps)
                {
                    p.y -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
                }
            }

            return p;
        }

        template<typename T>
        static gravis::d2Vector quadraticMaxWrapXY(
				const Image<T>& img, double eps = 1e-25,
				int wMax = -1, int hMax = -1)
        {
            const int w = img.data.xdim;
            const int h = img.data.ydim;
			
			if (wMax < 0) wMax = w;
			if (hMax < 0) hMax = h;

            int xmax = -1, ymax = -1;
            double vmax = 0.0;

            for (int y = 0; y < h; y++)
            for (int x = 0; x < w; x++)
            {
				if ((x > wMax && w - x > wMax)
				  || y > hMax && h - x > hMax)
				{
					continue;
				}
					
                T v = DIRECT_A2D_ELEM(img.data, y, x);

                if (xmax < 0 || v > vmax)
                {
                    vmax = v;
                    xmax = x;
                    ymax = y;
                }
            }

            gravis::d2Vector p(xmax, ymax);

            {
                const T vp = DIRECT_A2D_ELEM(img.data, ymax, (xmax+1)%w);
                const T vn = DIRECT_A2D_ELEM(img.data, ymax, (xmax-1+w)%w);

                if (std::abs(vp + vn - 2.0*vmax) > eps)
                {
                    p.x -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
                }
            }

            {
                const T vp = DIRECT_A2D_ELEM(img.data, (ymax+1)%w, xmax);
                const T vn = DIRECT_A2D_ELEM(img.data, (ymax-1+w)%w, xmax);

                if (std::abs(vp + vn - 2.0*vmax) > eps)
                {
                    p.y -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
                }
            }

            return p;
        }

        static double cubic1D(double y0, double y1, double y2, double y3, double t);

        template<typename T>
        static T cubicXY(const Image<T>& img, double x, double y, int z = 0, int n = 0, bool wrap = false)
        {
            int xi    = (int)std::floor(x);
            int yi    = (int)std::floor(y);
            int xi_n1 = xi-1;
            int yi_n1 = yi-1;
            int xi_p1 = xi+1;
            int yi_p1 = yi+1;
            int xi_p2 = xi+2;
            int yi_p2 = yi+2;

            const double xf = x - xi;
            const double yf = y - yi;

            if (wrap)
            {
				xi_n1 = INTERPOL_WRAP(xi_n1, img.data.xdim);
                yi_n1 = INTERPOL_WRAP(yi_n1, img.data.ydim);
				xi    = INTERPOL_WRAP(xi,    img.data.xdim);
                yi    = INTERPOL_WRAP(yi,    img.data.ydim);
                xi_p1 = INTERPOL_WRAP(xi_p1, img.data.xdim);
                yi_p1 = INTERPOL_WRAP(yi_p1, img.data.ydim);
                xi_p2 = INTERPOL_WRAP(xi_p2, img.data.xdim);
                yi_p2 = INTERPOL_WRAP(yi_p2, img.data.ydim);
            }
            else
            {
                xi    = XMIPP_MAX(0, XMIPP_MIN(img.data.xdim - 1, xi));
                yi    = XMIPP_MAX(0, XMIPP_MIN(img.data.ydim - 1, yi));
                xi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.data.xdim - 1, xi_n1));
                yi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.data.ydim - 1, yi_n1));
                xi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.data.xdim - 1, xi_p1));
                yi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.data.ydim - 1, yi_p1));
                xi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.data.xdim - 1, xi_p2));
                yi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.data.ydim - 1, yi_p2));
            }


            const T f00 = DIRECT_NZYX_ELEM(img.data, n, z, yi_n1,   xi_n1);
            const T f01 = DIRECT_NZYX_ELEM(img.data, n, z, yi_n1,   xi);
            const T f02 = DIRECT_NZYX_ELEM(img.data, n, z, yi_n1,   xi_p1);
            const T f03 = DIRECT_NZYX_ELEM(img.data, n, z, yi_n1,   xi_p2);

            const T f10 = DIRECT_NZYX_ELEM(img.data, n, z, yi,     xi_n1);
            const T f11 = DIRECT_NZYX_ELEM(img.data, n, z, yi,     xi);
            const T f12 = DIRECT_NZYX_ELEM(img.data, n, z, yi,     xi_p1);
            const T f13 = DIRECT_NZYX_ELEM(img.data, n, z, yi,     xi_p2);

            const T f20 = DIRECT_NZYX_ELEM(img.data, n, z, yi_p1,   xi_n1);
            const T f21 = DIRECT_NZYX_ELEM(img.data, n, z, yi_p1,   xi);
            const T f22 = DIRECT_NZYX_ELEM(img.data, n, z, yi_p1,   xi_p1);
            const T f23 = DIRECT_NZYX_ELEM(img.data, n, z, yi_p1,   xi_p2);

            const T f30 = DIRECT_NZYX_ELEM(img.data, n, z, yi_p2,   xi_n1);
            const T f31 = DIRECT_NZYX_ELEM(img.data, n, z, yi_p2,   xi);
            const T f32 = DIRECT_NZYX_ELEM(img.data, n, z, yi_p2,   xi_p1);
            const T f33 = DIRECT_NZYX_ELEM(img.data, n, z, yi_p2,   xi_p2);

            const gravis::d4Matrix A(  -1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
                                1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
                               -1.0/2.0,  0.0,      1.0/2.0,  0.0,
                                0.0,      1.0,      0.0,      0.0);

            const gravis::d4Matrix V(  f00, f10, f20, f30,
                               f01, f11, f21, f31,
                               f02, f12, f22, f32,
                               f03, f13, f23, f33);

            gravis::d4Matrix At = A;
            At.transpose();

            gravis::d4Matrix AVA = A * V * At;

            const gravis::d4Vector xx(xf*xf*xf, xf*xf, xf, 1.0);
            const gravis::d4Vector yy(yf*yf*yf, yf*yf, yf, 1.0);

            return (T)(xx.dot(AVA * yy));
        }

        template<typename T>
        static gravis::t2Vector<T> cubicXYgrad(const Image<T>& img, double x, double y, int z = 0, int n = 0, bool wrap = false)
        {
			int xi    = (int)std::floor(x);
            int yi    = (int)std::floor(y);
            int xi_n1 = xi-1;
            int yi_n1 = yi-1;
            int xi_p1 = xi+1;
            int yi_p1 = yi+1;
            int xi_p2 = xi+2;
            int yi_p2 = yi+2;

            const double xf = x - xi;
            const double yf = y - yi;

            if (wrap)
            {
				xi_n1 = INTERPOL_WRAP(xi_n1, img.data.xdim);
                yi_n1 = INTERPOL_WRAP(yi_n1, img.data.ydim);
				xi    = INTERPOL_WRAP(xi,    img.data.xdim);
                yi    = INTERPOL_WRAP(yi,    img.data.ydim);
                xi_p1 = INTERPOL_WRAP(xi_p1, img.data.xdim);
                yi_p1 = INTERPOL_WRAP(yi_p1, img.data.ydim);
                xi_p2 = INTERPOL_WRAP(xi_p2, img.data.xdim);
                yi_p2 = INTERPOL_WRAP(yi_p2, img.data.ydim);
            }
            else
            {
                xi    = XMIPP_MAX(0, XMIPP_MIN(img.data.xdim - 1, xi));
                yi    = XMIPP_MAX(0, XMIPP_MIN(img.data.ydim - 1, yi));
                xi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.data.xdim - 1, xi_n1));
                yi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.data.ydim - 1, yi_n1));
                xi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.data.xdim - 1, xi_p1));
                yi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.data.ydim - 1, yi_p1));
                xi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.data.xdim - 1, xi_p2));
                yi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.data.ydim - 1, yi_p2));
            }

            const T f00 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_n1,   xi_n1);
            const T f01 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_n1,   xi);
            const T f02 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_n1,   xi_p1);
            const T f03 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_n1,   xi_p2);

            const T f10 = DIRECT_NZYX_ELEM(img.data, n, 0, yi,     xi_n1);
            const T f11 = DIRECT_NZYX_ELEM(img.data, n, 0, yi,     xi);
            const T f12 = DIRECT_NZYX_ELEM(img.data, n, 0, yi,     xi_p1);
            const T f13 = DIRECT_NZYX_ELEM(img.data, n, 0, yi,     xi_p2);

            const T f20 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_p1,   xi_n1);
            const T f21 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_p1,   xi);
            const T f22 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_p1,   xi_p1);
            const T f23 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_p1,   xi_p2);

            const T f30 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_p2,   xi_n1);
            const T f31 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_p2,   xi);
            const T f32 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_p2,   xi_p1);
            const T f33 = DIRECT_NZYX_ELEM(img.data, n, 0, yi_p2,   xi_p2);

            const gravis::d4Matrix A(  -1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
                                1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
                               -1.0/2.0,  0.0,      1.0/2.0,  0.0,
                                0.0,      1.0,      0.0,      0.0);

            const gravis::d4Matrix V(  f00, f10, f20, f30,
                               f01, f11, f21, f31,
                               f02, f12, f22, f32,
                               f03, f13, f23, f33);

            gravis::d4Matrix At = A;
            At.transpose();

            gravis::d4Matrix AVA = A * V * At;

            const gravis::d4Vector xx(xf*xf*xf, xf*xf, xf, 1.0);
            const gravis::d4Vector yy(yf*yf*yf, yf*yf, yf, 1.0);
            const gravis::d4Vector xxd(3.0*xf*xf, 2.0*xf, 1.0, 0.0);
            const gravis::d4Vector yyd(3.0*yf*yf, 2.0*yf, 1.0, 0.0);

            return gravis::t2Vector<T>(xxd.dot(AVA * yy), xx.dot(AVA * yyd));
        }

        static void test2D();
};

#endif

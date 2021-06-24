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

#ifndef JAZ_NEW_INTERPOLATION_H
#define JAZ_NEW_INTERPOLATION_H

#include "raw_image.h"
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t4Matrix.h>


#define INTERPOL_WRAP(i,n) ((i) >= 0 ? (i) % (n) : -(i) % (n) ? (n) - (-(i) % (n)) : 0)

class Interpolation
{
    public:
		
		template<typename T> inline
		static T linearXY_clip(const RawImage<T>& img, double x, double y, int z = 0);

        template<typename T> inline
        static T linearXYZ_clip(const RawImage<T>& img, double x, double y, double z);
		
		template<typename T> inline
		static T cubicXY_clip(const RawImage<T>& img, double x, double y, int z = 0);	

		template<typename T> inline
		static gravis::t2Vector<T> cubicXYGrad_clip(const RawImage<T>& img, double x, double y, int z = 0);

		template<typename T> inline
		static gravis::t3Vector<T> cubicXYGradAndValue_clip(const RawImage<T>& img, double x, double y, int z = 0);

		template<typename T> inline
		static gravis::t3Vector<T> cubicXYGradAndValue_raw(const RawImage<T>& img, double x, double y, int z = 0);
		
		template<typename T> inline
		static T cubicXY_wrap(const RawImage<T>& img, double x, double y, int z = 0);
		
		template<typename T> inline
		static gravis::t2Vector<T> cubicXYGrad_wrap(const RawImage<T>& img, double x, double y, int z = 0);		
		
		static inline double linearXYkernel(double dx, double dy);			
		static inline double cubicXYkernel(double dx, double dy);
		
		
		/*template<typename T> inline
		static tComplex<T> linearXY_FftwHalf(
				const RawImage<tComplex<T>>& img, 
				double xSgn, double ySgn, int z = 0);*/
		
		
		
		template<typename T> inline
		static T linearXY_symmetric_FftwHalf_clip(
				const RawImage<T>& img, 
				double xSgn, double ySgn, int z = 0);
		
		template<typename T> inline
		static tComplex<T> linearXY_complex_FftwHalf_clip(
				const RawImage<tComplex<T>>& img, 
				double xSgn, double ySgn, int z = 0);
		
		
		template<typename T> inline
		static T linearXY_symmetric_FftwHalf_wrap(
				const RawImage<T>& img, 
				double xSgn, double ySgn, int z = 0);
		
		template<typename T> inline
		static tComplex<T> linearXY_complex_FftwHalf_wrap(
				const RawImage<tComplex<T>>& img, 
				double xSgn, double ySgn, int z = 0);
		
		

		template<typename T> inline
		static tComplex<T> linearXYZ_FftwHalf_complex(
				const RawImage<tComplex<T>>& img,
				double xSgn, double ySgn, double zSgn);

		template<typename T> inline
		static T linearXYZ_FftwHalf_generic(
				const RawImage<T>& img,
				double xSgn, double ySgn, double zSgn);
		
		template<typename T> inline
		static T linearXYZ_FftwHalf_real(
				const RawImage<T>& img, 
				double xSgn, double ySgn, double zSgn);

		template<typename T> inline
		static gravis::t3Vector<tComplex<T>> linearXYZGradient_FftwHalf_complex(
				const RawImage<tComplex<T>>& img,
				double xSgn, double ySgn, double zSgn);

		template<typename T> inline
		static gravis::t4Vector<tComplex<T>> linearXYZGradientAndValue_FftwHalf_complex(
				const RawImage<tComplex<T>>& img,
				double xSgn, double ySgn, double zSgn);

		template<typename T1, typename T2> inline
		static void insertLinearXY_FftwHalf(
				T1 val, T2 wgh, 
				double xSgn, double ySgn, 
				RawImage<T1>& dataTgt,
				RawImage<T2>& wghtTgt,
				int z = 0);
		
		
		template<typename T>
        static gravis::d2Vector quadraticMaxXY(const RawImage<T>& img, double eps = 1e-25);
		
		template<typename T>
        static gravis::d3Vector quadraticMaxXYV(const RawImage<T>& img, double eps = 1e-25);
		
		template<typename T> inline
        static gravis::d3Vector localQuadraticMaxXYZ(
				const RawImage<T>& img, int x, int y, int z, double eps = 1e-25);
		
		template<typename T> inline
        static gravis::d3Vector localQuadraticMinXYZ(
				const RawImage<T>& img, int x, int y, int z, double eps = 1e-25);
		
		
		template<typename T>
        static gravis::d2Vector discreteMaxXY(const RawImage<T>& img, double eps = 1e-25);
		
		template<typename T>
        static gravis::d3Vector discreteMaxXYV(const RawImage<T>& img, double eps = 1e-25);
		

        template<typename T>
        static gravis::d2Vector quadraticMaxWrapXY(
				const Image<T>& img, double eps = 1e-25,
				int wMax = -1, int hMax = -1);
		
		template<typename T1, typename T2>
		static inline double cosine(T2 x, T1 v0, T1 v1)
		{
			const T2 a = (T2)(0.5 - cos(PI*x) / 2.0);
			return a * v1 + (1.0 - a) * v0;
		}
		
};

template<typename T> inline
T Interpolation::linearXY_clip(const RawImage<T>& img, double x, double y, int z)
{
	int x0 = FLOOR(x);
	int y0 = FLOOR(y);
	
	int x1 = x0 + 1;
	int y1 = y0 + 1;
	
	const double xf = x - x0;
	const double yf = y - y0;
	
	if (x0 < 0) x0 = 0;
	if (x0 >= img.xdim) x0 = img.xdim - 1;
	if (x1 < 0) x1 = 0;
	if (x1 >= img.xdim) x1 = img.xdim - 1;
	if (y0 < 0) y0 = 0;
	if (y0 >= img.ydim) y0 = img.ydim - 1;
	if (y1 < 0) y1 = 0;
	if (y1 >= img.ydim) y1 = img.ydim - 1;
	
	const T vx0 = (1 - xf) * img(x0,y0,z) + xf * img(x1,y0,z);
	const T vx1 = (1 - xf) * img(x0,y1,z) + xf * img(x1,y1,z);
	
	return (1 - yf) * vx0 + yf * vx1;
}

template<typename T> inline
T Interpolation::linearXYZ_clip(const RawImage<T>& img, double x, double y, double z)
{
    int x0 = FLOOR(x);
    int y0 = FLOOR(y);
    int z0 = FLOOR(z);

    int x1 = x0 + 1;
    int y1 = y0 + 1;
    int z1 = z0 + 1;

    const double xf = x - x0;
    const double yf = y - y0;
    const double zf = z - z0;

    if (x0 < 0) x0 = 0;
    if (x0 >= img.xdim) x0 = img.xdim - 1;
    if (x1 < 0) x1 = 0;
    if (x1 >= img.xdim) x1 = img.xdim - 1;

    if (y0 < 0) y0 = 0;
    if (y0 >= img.ydim) y0 = img.ydim - 1;
    if (y1 < 0) y1 = 0;
    if (y1 >= img.ydim) y1 = img.ydim - 1;

    if (z0 < 0) z0 = 0;
    if (z0 >= img.zdim) z0 = img.zdim - 1;
    if (z1 < 0) z1 = 0;
    if (z1 >= img.zdim) z1 = img.zdim - 1;

    const T v_000 = img(x0,y0,z0);
    const T v_001 = img(x0,y0,z1);
    const T v_010 = img(x0,y1,z0);
    const T v_011 = img(x0,y1,z1);
    const T v_100 = img(x1,y0,z0);
    const T v_101 = img(x1,y0,z1);
    const T v_110 = img(x1,y1,z0);
    const T v_111 = img(x1,y1,z1);

    const T v_00 = (1 - zf) * v_000 + zf * v_001;
    const T v_01 = (1 - zf) * v_010 + zf * v_011;
    const T v_10 = (1 - zf) * v_100 + zf * v_101;
    const T v_11 = (1 - zf) * v_110 + zf * v_111;

    const T v_0 = (1 - yf) * v_00 + yf * v_01;
    const T v_1 = (1 - yf) * v_10 + yf * v_11;

    return (1 - xf) * v_0 + xf * v_1;
}

template<typename T> inline
T Interpolation::cubicXY_clip(const RawImage<T>& img, double x, double y, int z)
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

	xi    = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi));
	yi    = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi));
	xi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_n1));
	yi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_n1));
	xi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_p1));
	yi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_p1));
	xi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_p2));
	yi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_p2));

	const T f00 = img(xi_n1, yi_n1, z);
	const T f01 = img(xi,    yi_n1, z);
	const T f02 = img(xi_p1, yi_n1, z);
	const T f03 = img(xi_p2, yi_n1, z);

	const T f10 = img(xi_n1, yi,    z);
	const T f11 = img(xi,    yi,    z);
	const T f12 = img(xi_p1, yi,    z);
	const T f13 = img(xi_p2, yi,    z);

	const T f20 = img(xi_n1, yi_p1, z);
	const T f21 = img(xi,    yi_p1, z);
	const T f22 = img(xi_p1, yi_p1, z);
	const T f23 = img(xi_p2, yi_p1, z);

	const T f30 = img(xi_n1, yi_p2, z);
	const T f31 = img(xi,    yi_p2, z);
	const T f32 = img(xi_p1, yi_p2, z);
	const T f33 = img(xi_p2, yi_p2, z);

	const gravis::d4Matrix A(  
			-1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
			 1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
			-1.0/2.0,  0.0,      1.0/2.0,  0.0,
			 0.0,      1.0,      0.0,      0.0);

	const gravis::d4Matrix V(  
			f00, f10, f20, f30,
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

template<typename T> inline
gravis::t2Vector<T> Interpolation::cubicXYGrad_clip(const RawImage<T>& img, double x, double y, int z)
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

	xi    = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi));
	yi    = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi));
	xi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_n1));
	yi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_n1));
	xi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_p1));
	yi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_p1));
	xi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_p2));
	yi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_p2));

	const T f00 = img(xi_n1, yi_n1, z);
	const T f01 = img(xi,    yi_n1, z);
	const T f02 = img(xi_p1, yi_n1, z);
	const T f03 = img(xi_p2, yi_n1, z);

	const T f10 = img(xi_n1, yi,    z);
	const T f11 = img(xi,    yi,    z);
	const T f12 = img(xi_p1, yi,    z);
	const T f13 = img(xi_p2, yi,    z);

	const T f20 = img(xi_n1, yi_p1, z);
	const T f21 = img(xi,    yi_p1, z);
	const T f22 = img(xi_p1, yi_p1, z);
	const T f23 = img(xi_p2, yi_p1, z);

	const T f30 = img(xi_n1, yi_p2, z);
	const T f31 = img(xi,    yi_p2, z);
	const T f32 = img(xi_p1, yi_p2, z);
	const T f33 = img(xi_p2, yi_p2, z);

	const gravis::d4Matrix A(
			-1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
			 1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
			-1.0/2.0,  0.0,      1.0/2.0,  0.0,
			 0.0,      1.0,      0.0,      0.0);

	const gravis::d4Matrix V(
			f00, f10, f20, f30,
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


template<typename T> inline
gravis::t3Vector<T> Interpolation::cubicXYGradAndValue_clip(const RawImage<T>& img, double x, double y, int z)
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

	xi    = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi));
	yi    = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi));
	xi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_n1));
	yi_n1 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_n1));
	xi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_p1));
	yi_p1 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_p1));
	xi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.xdim - 1, xi_p2));
	yi_p2 = XMIPP_MAX(0, XMIPP_MIN(img.ydim - 1, yi_p2));

	const T f00 = img(xi_n1, yi_n1, z);
	const T f01 = img(xi,    yi_n1, z);
	const T f02 = img(xi_p1, yi_n1, z);
	const T f03 = img(xi_p2, yi_n1, z);

	const T f10 = img(xi_n1, yi,    z);
	const T f11 = img(xi,    yi,    z);
	const T f12 = img(xi_p1, yi,    z);
	const T f13 = img(xi_p2, yi,    z);

	const T f20 = img(xi_n1, yi_p1, z);
	const T f21 = img(xi,    yi_p1, z);
	const T f22 = img(xi_p1, yi_p1, z);
	const T f23 = img(xi_p2, yi_p1, z);

	const T f30 = img(xi_n1, yi_p2, z);
	const T f31 = img(xi,    yi_p2, z);
	const T f32 = img(xi_p1, yi_p2, z);
	const T f33 = img(xi_p2, yi_p2, z);

	const gravis::d4Matrix A(
			-1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
			 1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
			-1.0/2.0,  0.0,      1.0/2.0,  0.0,
			 0.0,      1.0,      0.0,      0.0);

	const gravis::d4Matrix V(
			f00, f10, f20, f30,
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

	return gravis::t3Vector<T>(xxd.dot(AVA * yy), xx.dot(AVA * yyd), xx.dot(AVA * yy));
}

template<typename T> inline
gravis::t3Vector<T> Interpolation::cubicXYGradAndValue_raw(const RawImage<T>& img, double x, double y, int z)
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

	const T f00 = img(xi_n1, yi_n1, z);
	const T f01 = img(xi,    yi_n1, z);
	const T f02 = img(xi_p1, yi_n1, z);
	const T f03 = img(xi_p2, yi_n1, z);

	const T f10 = img(xi_n1, yi,    z);
	const T f11 = img(xi,    yi,    z);
	const T f12 = img(xi_p1, yi,    z);
	const T f13 = img(xi_p2, yi,    z);

	const T f20 = img(xi_n1, yi_p1, z);
	const T f21 = img(xi,    yi_p1, z);
	const T f22 = img(xi_p1, yi_p1, z);
	const T f23 = img(xi_p2, yi_p1, z);

	const T f30 = img(xi_n1, yi_p2, z);
	const T f31 = img(xi,    yi_p2, z);
	const T f32 = img(xi_p1, yi_p2, z);
	const T f33 = img(xi_p2, yi_p2, z);

	const gravis::d4Matrix A(
			-1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
			 1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
			-1.0/2.0,  0.0,      1.0/2.0,  0.0,
			 0.0,      1.0,      0.0,      0.0);

	const gravis::d4Matrix V(
			f00, f10, f20, f30,
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

	return gravis::t3Vector<T>(xxd.dot(AVA * yy), xx.dot(AVA * yyd), xx.dot(AVA * yy));
}



template<typename T> inline
T Interpolation::cubicXY_wrap(const RawImage<T>& img, double x, double y, int z)
{
	const int w = img.xdim;
	const int h = img.ydim;
	
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
	
	xi    -= std::floor( ((double)xi) / w) * w;
	yi    -= std::floor( ((double)yi) / h) * h;
	xi_n1 -= std::floor( ((double)xi_n1) / w) * w;
	yi_n1 -= std::floor( ((double)yi_n1) / h) * h;
	xi_p1 -= std::floor( ((double)xi_p1) / w) * w;
	yi_p1 -= std::floor( ((double)yi_p1) / h) * h;
	xi_p2 -= std::floor( ((double)xi_p2) / w) * w;
	yi_p2 -= std::floor( ((double)yi_p2) / h) * h;

	const T f00 = img(xi_n1, yi_n1, z);
	const T f01 = img(xi,    yi_n1, z);
	const T f02 = img(xi_p1, yi_n1, z);
	const T f03 = img(xi_p2, yi_n1, z);

	const T f10 = img(xi_n1, yi,    z);
	const T f11 = img(xi,    yi,    z);
	const T f12 = img(xi_p1, yi,    z);
	const T f13 = img(xi_p2, yi,    z);

	const T f20 = img(xi_n1, yi_p1, z);
	const T f21 = img(xi,    yi_p1, z);
	const T f22 = img(xi_p1, yi_p1, z);
	const T f23 = img(xi_p2, yi_p1, z);

	const T f30 = img(xi_n1, yi_p2, z);
	const T f31 = img(xi,    yi_p2, z);
	const T f32 = img(xi_p1, yi_p2, z);
	const T f33 = img(xi_p2, yi_p2, z);

	const gravis::d4Matrix A(  
			-1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
			 1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
			-1.0/2.0,  0.0,      1.0/2.0,  0.0,
			 0.0,      1.0,      0.0,      0.0);

	const gravis::d4Matrix V(  
			f00, f10, f20, f30,
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

template<typename T> inline
gravis::t2Vector<T> Interpolation::cubicXYGrad_wrap(const RawImage<T>& img, double x, double y, int z)
{
	const int w = img.xdim;
	const int h = img.ydim;
	
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
	
	xi    -= std::floor( ((double)xi) / w) * w;
	yi    -= std::floor( ((double)yi) / h) * h;
	xi_n1 -= std::floor( ((double)xi_n1) / w) * w;
	yi_n1 -= std::floor( ((double)yi_n1) / h) * h;
	xi_p1 -= std::floor( ((double)xi_p1) / w) * w;
	yi_p1 -= std::floor( ((double)yi_p1) / h) * h;
	xi_p2 -= std::floor( ((double)xi_p2) / w) * w;
	yi_p2 -= std::floor( ((double)yi_p2) / h) * h;

	const T f00 = img(xi_n1, yi_n1, z);
	const T f01 = img(xi,    yi_n1, z);
	const T f02 = img(xi_p1, yi_n1, z);
	const T f03 = img(xi_p2, yi_n1, z);

	const T f10 = img(xi_n1, yi,    z);
	const T f11 = img(xi,    yi,    z);
	const T f12 = img(xi_p1, yi,    z);
	const T f13 = img(xi_p2, yi,    z);

	const T f20 = img(xi_n1, yi_p1, z);
	const T f21 = img(xi,    yi_p1, z);
	const T f22 = img(xi_p1, yi_p1, z);
	const T f23 = img(xi_p2, yi_p1, z);

	const T f30 = img(xi_n1, yi_p2, z);
	const T f31 = img(xi,    yi_p2, z);
	const T f32 = img(xi_p1, yi_p2, z);
	const T f33 = img(xi_p2, yi_p2, z);

	const gravis::d4Matrix A(  
			-1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
			 1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
			-1.0/2.0,  0.0,      1.0/2.0,  0.0,
			 0.0,      1.0,      0.0,      0.0);

	const gravis::d4Matrix V(  
			f00, f10, f20, f30,
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

inline double Interpolation::linearXYkernel(double dx, double dy)
{
	const double ax = 1 - std::abs(dx);
	const double ay = 1 - std::abs(dy);
	
	if (ax > 0 && ay > 0) return ax*ay;
	else return 0.0;
}

inline double Interpolation::cubicXYkernel(double dx, double dy)
{
	gravis::d4Matrix V(  
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0,
			0, 0, 0, 0);
	
	int xi = (int) std::floor(dx);
	int yi = (int) std::floor(dy);
	
	if (xi < -2 || xi > 1 || yi < -2 || yi > 1) return 0.0;
	
	double xf = dx - xi;
	double yf = dy - yi;
	
	V(1 - xi, 1 - yi) = 1.0;
	
	const gravis::d4Matrix A(  
			-1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
			 1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
			-1.0/2.0,  0.0,      1.0/2.0,  0.0,
			 0.0,      1.0,      0.0,      0.0);
	
	gravis::d4Matrix At = A;
	At.transpose();

	gravis::d4Matrix AVA = A * V * At;

	const gravis::d4Vector xx(xf*xf*xf, xf*xf, xf, 1.0);
	const gravis::d4Vector yy(yf*yf*yf, yf*yf, yf, 1.0);

	return xx.dot(AVA * yy);
}

template<typename T> inline
T Interpolation::linearXY_symmetric_FftwHalf_clip(
		const RawImage<T>& img, 
		double xSgn, double ySgn, int z)
{
	double xd, yd;
	
	if (xSgn >= 0.0)
	{
		xd = xSgn;
		yd = ySgn;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
	}
	
	if (yd < 0.0) yd += img.ydim;
	
	int x0 = (int) xd;
	int y0 = (int) yd;
	
	const double xf = xd - x0;
	const double yf = yd - y0;
	
	if (x0 >= img.xdim) x0 = img.xdim-1;
	
	int x1 = x0 + 1;	
	if (x1 >= img.xdim) x1 = img.xdim-2;
	
	if (y0 < 0) y0 = 0;
	else if (y0 >= img.ydim) y0 = img.ydim-1;
	
	int y1 = (y0 + 1) % img.ydim;
	
	const tComplex<T> vx0 = (1 - xf) * img(x0,y0,z) + xf * img(x1,y0,z);
	const tComplex<T> vx1 = (1 - xf) * img(x0,y1,z) + xf * img(x1,y1,z);
	return (1 - yf) * vx0 + yf * vx1;
}

template<typename T> inline
tComplex<T> Interpolation::linearXY_complex_FftwHalf_clip(
		const RawImage<tComplex<T>>& img, 
		double xSgn, double ySgn, int z)
{
	double xd, yd;
	bool conj;
	
	if (xSgn >= 0.0)
	{
		xd = xSgn;
		yd = ySgn;
		conj = false;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
		conj = true;
	}
	
	if (yd < 0.0) yd += img.ydim;
	
	int x0 = (int) xd;
	int y0 = (int) yd;
	
	const double xf = xd - x0;
	const double yf = yd - y0;
	
	if (x0 >= img.xdim) x0 = img.xdim-1;
	
	int x1 = x0 + 1;	
	if (x1 >= img.xdim) x1 = img.xdim-2;
	
	if (y0 < 0) y0 = 0;
	else if (y0 >= img.ydim) y0 = img.ydim-1;
	
	int y1 = (y0 + 1) % img.ydim;
	
	const tComplex<T> vx0 = (1 - xf) * img(x0,y0,z) + xf * img(x1,y0,z);
	const tComplex<T> vx1 = (1 - xf) * img(x0,y1,z) + xf * img(x1,y1,z);
	const tComplex<T> v = (1 - yf) * vx0 + yf * vx1;
		
	return conj? v.conj() : v;
}



template<typename T> inline
T Interpolation::linearXY_symmetric_FftwHalf_wrap(
		const RawImage<T>& img, 
		double xSgn, double ySgn, int z)
{
	double xd, yd;
	
	if (xSgn >= 0.0)
	{
		xd = xSgn;
		yd = ySgn;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
	}
	
	if (yd < 0.0) yd += img.ydim;
	
	int x0 = (int) std::floor(xd);
	int y0 = (int) std::floor(yd);
	
	const double xf = xd - x0;
	const double yf = yd - y0;
	
	T sum  = 0;
	
	const int wh = img.xdim;
	const int w = 2 * (wh - 1);
	const int h = img.ydim;
	
	for (int dy = 0; dy < 2; dy++)
	for (int dx = 0; dx < 2; dx++)
	{
		int xg = x0 + dx;
		int yg = y0 + dy;
		
		xg = xg % w;
		
		if (xg >= wh)
		{
			xg = 2*wh - xg - 2;
			yg = (h - yg) % h;
		}
		
		if (yg < 0)
		{
			yg += (int)((-1 - yg) / h + 1) * h;
		}
		else if (yg >= img.ydim)
		{
			yg -= (int)(yg / img.ydim);
		}
		
		const double ax = 1 - std::abs(xf - (double)dx);
		const double ay = 1 - std::abs(yf - (double)dy);
		
		sum += ax * ay * img(xg, yg, z);
	}
	
	return sum;
}

template<typename T> inline
tComplex<T> Interpolation::linearXY_complex_FftwHalf_wrap(
		const RawImage<tComplex<T>>& img, 
		double xSgn, double ySgn, int z)
{
	double xd, yd;
	const bool conj0 = xSgn < 0.0;
	
	if (xSgn >= 0.0)
	{
		xd = xSgn;
		yd = ySgn;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
	}
	
	if (yd < 0.0) yd += img.ydim;
	
	int x0 = (int) std::floor(xd);
	int y0 = (int) std::floor(yd);
	
	const double xf = xd - x0;
	const double yf = yd - y0;
	
	tComplex<T> sum  = 0;
	
	const int wh = img.xdim;
	const int w = 2 * (wh - 1);
	const int h = img.ydim;

	for (int dy = 0; dy < 2; dy++)
	for (int dx = 0; dx < 2; dx++)
	{
		int xg = x0 + dx;
		int yg = y0 + dy;

		xg = xg % w;

		bool conj;
		
		if (xg >= wh)
		{
			xg = 2*wh - xg - 2;
			yg = (h - yg) % h;
			
			conj = !conj0;
		}
		else
		{
			conj = conj0;
		}
		
		if (yg < 0)
		{
			yg += (int)((-1 - yg) / h + 1) * h;
		}
		else if (yg >= img.ydim)
		{
			yg -= (int)(yg / img.ydim);
		}
		
		const double ax = 1 - std::abs(xf - (double)dx);
		const double ay = 1 - std::abs(yf - (double)dy);
		
		tComplex<T> val = img(xg, yg, z);
		
		if (conj) val.imag *= -1;
		
		sum += ax * ay * val;
	}
	
	return sum;
}


template<typename T> inline
tComplex<T> Interpolation::linearXYZ_FftwHalf_complex(
		const RawImage<tComplex<T>>& img, 
		double xSgn, double ySgn, double zSgn)
{
	double xd, yd, zd;
	bool conj;
	
	if (xSgn > 0.0)
	{
		xd = xSgn;
		yd = ySgn;
		zd = zSgn;
		conj = false;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
		zd = -zSgn;
		conj = true;
	}
	
	if (yd < 0.0) yd += img.ydim;
	if (zd < 0.0) zd += img.zdim;
	
	int x0 = (int) xd;
	int y0 = (int) yd;
	int z0 = (int) zd;
	
	const double xf = xd - x0;
	const double yf = yd - y0;
	const double zf = zd - z0;
	
	if (x0 >= img.xdim) x0 = img.xdim-1;
	
	int x1 = x0 + 1;	
	if (x1 >= img.xdim) x1 = img.xdim-2;
	
	
	if (y0 < 0) y0 = 0;
	else if (y0 >= img.ydim) y0 = img.ydim-1;
	
	int y1 = (y0 + 1) % img.ydim;
	
	if (z0 < 0) z0 = 0;
	else if (z0 >= img.zdim) z0 = img.zdim-1;
	
	int z1 = (z0 + 1) % img.zdim;
	
	const tComplex<T> vx00 = (1 - xf) * img(x0,y0,z0) + xf * img(x1,y0,z0);
	const tComplex<T> vx10 = (1 - xf) * img(x0,y1,z0) + xf * img(x1,y1,z0);
	const tComplex<T> vx01 = (1 - xf) * img(x0,y0,z1) + xf * img(x1,y0,z1);
	const tComplex<T> vx11 = (1 - xf) * img(x0,y1,z1) + xf * img(x1,y1,z1);
	
	const tComplex<T> vxy0 = (1 - yf) * vx00 + yf * vx10;
	const tComplex<T> vxy1 = (1 - yf) * vx01 + yf * vx11;
	
	const tComplex<T> vxyz = (1 - zf) * vxy0 + zf * vxy1;
		
	return conj? vxyz.conj() : vxyz;
}


template<typename T> inline
gravis::t3Vector<tComplex<T>> Interpolation::linearXYZGradient_FftwHalf_complex(
		const RawImage<tComplex<T>>& img,
		double xSgn, double ySgn, double zSgn)
{
	double xd, yd, zd;
	bool conj;

	if (xSgn > 0.0)
	{
		xd = xSgn;
		yd = ySgn;
		zd = zSgn;
		conj = false;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
		zd = -zSgn;
		conj = true;
	}

	if (yd < 0.0) yd += img.ydim;
	if (zd < 0.0) zd += img.zdim;

	int x0 = (int) xd;
	int y0 = (int) yd;
	int z0 = (int) zd;

	const double xf = xd - x0;
	const double yf = yd - y0;
	const double zf = zd - z0;

	if (x0 >= img.xdim) x0 = img.xdim-1;

	int x1 = x0 + 1;
	if (x1 >= img.xdim) x1 = img.xdim-2;


	if (y0 < 0) y0 = 0;
	else if (y0 >= img.ydim) y0 = img.ydim-1;

	int y1 = (y0 + 1) % img.ydim;

	if (z0 < 0) z0 = 0;
	else if (z0 >= img.zdim) z0 = img.zdim-1;

	int z1 = (z0 + 1) % img.zdim;

	const tComplex<T> vx00 = (1 - xf) * img(x0,y0,z0) + xf * img(x1,y0,z0);
	const tComplex<T> vx10 = (1 - xf) * img(x0,y1,z0) + xf * img(x1,y1,z0);
	const tComplex<T> vx01 = (1 - xf) * img(x0,y0,z1) + xf * img(x1,y0,z1);
	const tComplex<T> vx11 = (1 - xf) * img(x0,y1,z1) + xf * img(x1,y1,z1);

	const tComplex<T> vx00_x = img(x1,y0,z0) - img(x0,y0,z0);
	const tComplex<T> vx10_x = img(x1,y1,z0) - img(x0,y1,z0);
	const tComplex<T> vx01_x = img(x1,y0,z1) - img(x0,y0,z1);
	const tComplex<T> vx11_x = img(x1,y1,z1) - img(x0,y1,z1);

	const tComplex<T> vxy0 = (1 - yf) * vx00 + yf * vx10;
	const tComplex<T> vxy1 = (1 - yf) * vx01 + yf * vx11;

	const tComplex<T> vxy0_y = vx10 - vx00;
	const tComplex<T> vxy1_y = vx11 - vx01;

	const tComplex<T> vxyz_z = vxy1 - vxy0;

	gravis::t3Vector<tComplex<T>> out(
		(1 - zf) * ((1 - yf) * vx00_x + yf * vx10_x) + zf * ((1 - yf) * vx01_x + yf * vx11_x),
		(1 - zf) * (vxy0_y) + zf * (vxy1_y),
		vxyz_z);

	if (conj)
	{
		for (int i = 0; i < 3; i++)
		{
			out[i].real *= -1;
		}
	}

	return out;
}

template<typename T> inline
gravis::t4Vector<tComplex<T>> Interpolation::linearXYZGradientAndValue_FftwHalf_complex(
		const RawImage<tComplex<T>>& img,
		double xSgn, double ySgn, double zSgn)
{
	double xd, yd, zd;
	bool conj;

	if (xSgn > 0.0)
	{
		xd = xSgn;
		yd = ySgn;
		zd = zSgn;
		conj = false;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
		zd = -zSgn;
		conj = true;
	}

	if (yd < 0.0) yd += img.ydim;
	if (zd < 0.0) zd += img.zdim;

	int x0 = (int) xd;
	int y0 = (int) yd;
	int z0 = (int) zd;

	const double xf = xd - x0;
	const double yf = yd - y0;
	const double zf = zd - z0;

	if (x0 >= img.xdim) x0 = img.xdim-1;

	int x1 = x0 + 1;
	if (x1 >= img.xdim) x1 = img.xdim-2;


	if (y0 < 0) y0 = 0;
	else if (y0 >= img.ydim) y0 = img.ydim-1;

	int y1 = (y0 + 1) % img.ydim;

	if (z0 < 0) z0 = 0;
	else if (z0 >= img.zdim) z0 = img.zdim-1;

	int z1 = (z0 + 1) % img.zdim;

	const tComplex<T> vx00 = (1 - xf) * img(x0,y0,z0) + xf * img(x1,y0,z0);
	const tComplex<T> vx10 = (1 - xf) * img(x0,y1,z0) + xf * img(x1,y1,z0);
	const tComplex<T> vx01 = (1 - xf) * img(x0,y0,z1) + xf * img(x1,y0,z1);
	const tComplex<T> vx11 = (1 - xf) * img(x0,y1,z1) + xf * img(x1,y1,z1);

	const tComplex<T> vx00_x = img(x1,y0,z0) - img(x0,y0,z0);
	const tComplex<T> vx10_x = img(x1,y1,z0) - img(x0,y1,z0);
	const tComplex<T> vx01_x = img(x1,y0,z1) - img(x0,y0,z1);
	const tComplex<T> vx11_x = img(x1,y1,z1) - img(x0,y1,z1);

	const tComplex<T> vxy0 = (1 - yf) * vx00 + yf * vx10;
	const tComplex<T> vxy1 = (1 - yf) * vx01 + yf * vx11;

	const tComplex<T> vxy0_y = vx10 - vx00;
	const tComplex<T> vxy1_y = vx11 - vx01;

	const tComplex<T> vxyz_z = vxy1 - vxy0;

	const tComplex<T> vxyz = (1 - zf) * vxy0 + zf * vxy1;

	gravis::t4Vector<tComplex<T>> out(
		(1 - zf) * ((1 - yf) * vx00_x + yf * vx10_x) + zf * ((1 - yf) * vx01_x + yf * vx11_x),
		(1 - zf) * (vxy0_y) + zf * (vxy1_y),
		vxyz_z,
		vxyz);

	if (conj)
	{
		for (int i = 0; i < 3; i++)
		{
			out[i].real *= -1;
		}

		out[3].imag *= -1;
	}

	return out;
}


template<typename Conjugable> inline
Conjugable Interpolation::linearXYZ_FftwHalf_generic(
		const RawImage<Conjugable>& img,
		double xSgn, double ySgn, double zSgn)
{
	double xd, yd, zd;
	bool conj;

	if (xSgn > 0.0)
	{
		xd = xSgn;
		yd = ySgn;
		zd = zSgn;
		conj = false;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
		zd = -zSgn;
		conj = true;
	}

	if (yd < 0.0) yd += img.ydim;
	if (zd < 0.0) zd += img.zdim;

	int x0 = (int) xd;
	int y0 = (int) yd;
	int z0 = (int) zd;

	const double xf = xd - x0;
	const double yf = yd - y0;
	const double zf = zd - z0;

	if (x0 >= img.xdim) x0 = img.xdim-1;

	int x1 = x0 + 1;
	if (x1 >= img.xdim) x1 = img.xdim-2;


	if (y0 < 0) y0 = 0;
	else if (y0 >= img.ydim) y0 = img.ydim-1;

	int y1 = (y0 + 1) % img.ydim;

	if (z0 < 0) z0 = 0;
	else if (z0 >= img.zdim) z0 = img.zdim-1;

	int z1 = (z0 + 1) % img.zdim;

	const Conjugable vx00 = (1 - xf) * img(x0,y0,z0) + xf * img(x1,y0,z0);
	const Conjugable vx10 = (1 - xf) * img(x0,y1,z0) + xf * img(x1,y1,z0);
	const Conjugable vx01 = (1 - xf) * img(x0,y0,z1) + xf * img(x1,y0,z1);
	const Conjugable vx11 = (1 - xf) * img(x0,y1,z1) + xf * img(x1,y1,z1);

	const Conjugable vxy0 = (1 - yf) * vx00 + yf * vx10;
	const Conjugable vxy1 = (1 - yf) * vx01 + yf * vx11;

	const Conjugable vxyz = (1 - zf) * vxy0 + zf * vxy1;

	return conj? vxyz.conj() : vxyz;
}

template<typename T> inline
T Interpolation::linearXYZ_FftwHalf_real(
		const RawImage<T>& img, 
		double xSgn, double ySgn, double zSgn)
{
	double xd, yd, zd;
	
	if (xSgn > 0.0)
	{
		xd = xSgn;
		yd = ySgn;
		zd = zSgn;
	}
	else
	{
		xd = -xSgn;
		yd = -ySgn;
		zd = -zSgn;
	}
	
	if (yd < 0.0) yd += img.ydim;
	if (zd < 0.0) zd += img.zdim;
	
	int x0 = (int) xd;
	int y0 = (int) yd;
	int z0 = (int) zd;
	
	const double xf = xd - x0;
	const double yf = yd - y0;
	const double zf = zd - z0;
	
	if (x0 >= img.xdim) x0 = img.xdim-1;
	
	int x1 = x0 + 1;	
	if (x1 >= img.xdim) x1 = img.xdim-2;
	
	
	if (y0 < 0) y0 = 0;
	else if (y0 >= img.ydim) y0 = img.ydim-1;
	
	int y1 = (y0 + 1) % img.ydim;
	
	if (z0 < 0) z0 = 0;
	else if (z0 >= img.zdim) z0 = img.zdim-1;
	
	int z1 = (z0 + 1) % img.zdim;
	
	const T vx00 = (1 - xf) * img(x0,y0,z0) + xf * img(x1,y0,z0);
	const T vx10 = (1 - xf) * img(x0,y1,z0) + xf * img(x1,y1,z0);
	const T vx01 = (1 - xf) * img(x0,y0,z1) + xf * img(x1,y0,z1);
	const T vx11 = (1 - xf) * img(x0,y1,z1) + xf * img(x1,y1,z1);
	
	const T vxy0 = (1 - yf) * vx00 + yf * vx10;
	const T vxy1 = (1 - yf) * vx01 + yf * vx11;
	
	const T vxyz = (1 - zf) * vxy0 + zf * vxy1;
		
	return vxyz;
}
	
template<typename T1, typename T2> inline
void Interpolation::insertLinearXY_FftwHalf(
		T1 val, T2 wgh, 
		double xSgn, double ySgn, 
		RawImage<T1>& dataTgt,
		RawImage<T2>& wghtTgt,
		int z)
{
	if (xSgn < 0.0)
	{
		xSgn = -xSgn;
		ySgn = -ySgn;
	}
	
	const int s = dataTgt.ydim;
	const int sh = dataTgt.xdim;
	
	const int x0 = (int) xSgn;
	const int x1 = (int) xSgn + 1;
	const int y0 = (int) (ySgn < 0.0? ySgn + s : ySgn);
	const int y1 = (int) (ySgn < -1.0? ySgn + s + 1: ySgn + 1);
	
	const double dx = xSgn - x0;
	const double dy = ySgn - std::floor(ySgn);
	
	if (x0 < sh)
	{
		if (y0 >= 0 && y0 < s)
		{
			dataTgt(x0, y0) += wgh * (1 - dx) * (1 - dy) * val;
			wghtTgt(x0, y0) += wgh * (1 - dx) * (1 - dy);
		}
		
		if (y1 >= 0 && y1 < s)
		{
			dataTgt(x0, y1) += wgh * (1 - dx) * dy * val;
			wghtTgt(x0, y1) += wgh * (1 - dx) * dy;
		}
	}
	
	if (x1 < sh)
	{
		if (y0 >= 0 && y0 < s)
		{
			dataTgt(x1, y0) += wgh * dx * (1 - dy) * val;
			wghtTgt(x1, y0) += wgh * dx * (1 - dy);
		}
		
		if (y1 >= 0 && y1 < s)
		{
			dataTgt(x1, y1) += wgh * dx * dy * val;
			wghtTgt(x1, y1) += wgh * dx * dy;
		}
	}
}

template<typename T>
gravis::d2Vector Interpolation::quadraticMaxXY(const RawImage<T>& img, double eps)
{
	const int w = img.xdim;
	const int h = img.ydim;

	int xmax = -1, ymax = -1;
	double vmax = 0.0;

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		T v = img(x,y);

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
		const T vp = img(xmax + 1, ymax);
		const T vn = img(xmax - 1, ymax);

		if (std::abs(vp + vn - 2.0*vmax) > eps)
		{
			p.x -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
		}
	}

	if (ymax > 0 && ymax < h-1)
	{
		const T vp = img(xmax, ymax + 1);
		const T vn = img(xmax, ymax - 1);

		if (std::abs(vp + vn - 2.0*vmax) > eps)
		{
			p.y -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
		}
	}

	return p;
}

template<typename T>
gravis::d3Vector Interpolation::quadraticMaxXYV(const RawImage<T>& img, double eps)
{
	gravis::d3Vector out;	
	gravis::d2Vector p = quadraticMaxXY(img, eps);
	
	out.x = p.x;
	out.y = p.y;
	out.z = cubicXY_clip(img, p.x, p.y);
			
	return out;
}

template<typename T> inline
gravis::d3Vector Interpolation::localQuadraticMaxXYZ(const RawImage<T>& img, int x, int y, int z, double eps)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	gravis::d3Vector p(x,y,z);
	
	T vmax = img(x,y,z);
	
	if (x > 0 && x < w - 1)
	{
		const T vp = img(x + 1, y, z);
		const T vn = img(x - 1, y, z);

		if (std::abs(vp + vn - 2.0*vmax) > eps)
		{
			p.x -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
		}
	}
	
	if (y > 0 && y < h - 1)
	{
		const T vp = img(x, y + 1, z);
		const T vn = img(x, y - 1, z);

		if (std::abs(vp + vn - 2.0*vmax) > eps)
		{
			p.y -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
		}
	}
	
	if (z > 0 && z < d - 1)
	{
		const T vp = img(x, y, z + 1);
		const T vn = img(x, y, z - 1);

		if (std::abs(vp + vn - 2.0*vmax) > eps)
		{
			p.z -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmax);
		}
	}

	return p;
}



template<typename T> inline
gravis::d3Vector Interpolation::localQuadraticMinXYZ(const RawImage<T>& img, int x, int y, int z, double eps)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	gravis::d3Vector p(x,y,z);
	
	T vmin = img(x,y,z);
	
	if (x > 0 && x < w - 1)
	{
		const T vp = img(x + 1, y, z);
		const T vn = img(x - 1, y, z);

		if (std::abs(vp + vn - 2.0*vmin) > eps)
		{
			p.x -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmin);
		}
	}
	
	if (y > 0 && y < h - 1)
	{
		const T vp = img(x, y + 1, z);
		const T vn = img(x, y - 1, z);

		if (std::abs(vp + vn - 2.0*vmin) > eps)
		{
			p.y -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmin);
		}
	}
	
	if (z > 0 && z < d - 1)
	{
		const T vp = img(x, y, z + 1);
		const T vn = img(x, y, z - 1);

		if (std::abs(vp + vn - 2.0*vmin) > eps)
		{
			p.z -= 0.5 * (vp - vn) / (vp + vn - 2.0*vmin);
		}
	}

	return p;
}


template<typename T>
gravis::d2Vector Interpolation::discreteMaxXY(const RawImage<T>& img, double eps)
{
	const int w = img.xdim;
	const int h = img.ydim;

	int xmax = -1, ymax = -1;
	double vmax = 0.0;

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		T v = img(x,y);

		if (xmax < 0 || v > vmax)
		{
			vmax = v;
			xmax = x;
			ymax = y;
		}
	}

	return gravis::d2Vector(xmax, ymax);
}

template<typename T>
gravis::d3Vector Interpolation::discreteMaxXYV(const RawImage<T>& img, double eps)
{
	gravis::d3Vector out;	
	gravis::d2Vector p = discreteMaxXY(img, eps);
	
	out.x = p.x;
	out.y = p.y;
	out.z = cubicXY_clip(img, p.x, p.y);
			
	return out;
}

template<typename T>
gravis::d2Vector Interpolation::quadraticMaxWrapXY(
		const Image<T>& img, double eps,
		int wMax, int hMax)
{
	const int w = img.xdim;
	const int h = img.ydim;
	
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

#endif

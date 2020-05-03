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

#include <src/jaz/single_particle/interpolation.h>
#include <src/jaz/single_particle/image_log.h>

using namespace gravis;

bool Interpolation::isInSlice(const Image<RFLOAT>& img, double x, double y)
{
    return x >= 0.0 && x < img.data.xdim-1 && y >= 0.0 && y < img.data.ydim-1;
}

double Interpolation::getTaperWeight(const Image<RFLOAT>& img, double x, double y, double rx, double ry)
{
    double wx(1.0), wy(1.0);

    if (x < rx) wx *= (1.0 - cos(PI * (x+1) / rx))/2.0;
    if (x >= img.data.xdim - rx) wx *= (1.0 - cos(PI * (img.data.xdim - x) / rx))/2.0;

    if (y < ry) wy *= (1.0 - cos(PI * (y+1) / ry))/2.0;
    if (y >= img.data.ydim - ry) wy *= (1.0 - cos(PI * (img.data.ydim - y) / ry))/2.0;

    return wx * wy;
}

double Interpolation::linearXY(const Image<RFLOAT>& img, double x, double y, int n)
{
    if (!(x >= 0.0 && x < img.data.xdim-1 && y >= 0.0 && y < img.data.ydim-1))
    {
        return 0.0;
    }

    const int xi = (int)x;
    const int yi = (int)y;

    const double xf = x - xi;
    const double yf = y - yi;

    const double f00 = DIRECT_NZYX_ELEM(img.data, n, 0, yi,   xi);
    const double f01 = DIRECT_NZYX_ELEM(img.data, n, 0, yi,   xi+1);
    const double f10 = DIRECT_NZYX_ELEM(img.data, n, 0, yi+1, xi);
    const double f11 = DIRECT_NZYX_ELEM(img.data, n, 0, yi+1, xi+1);

    const double f0 = xf * f01 + (1.0 - xf) * f00;
    const double f1 = xf * f11 + (1.0 - xf) * f10;

    const double f = yf * f1 + (1.0 - yf) * f0;

    return f;
}

double Interpolation::cubic1D(double y0, double y1, double y2, double y3, double t)
{
    const d4Matrix A(  -1.0/2.0,  3.0/2.0, -3.0/2.0,  1.0/2.0,
                        1.0,     -5.0/2.0,  2.0,     -1.0/2.0,
                       -1.0/2.0,  0.0,      1.0/2.0,  0.0,
                        0.0,      1.0,      0.0,      0.0);

    d4Vector y(y0, y1, y2, y3);
    d4Vector c = A*y;

    d4Vector x(t*t*t, t*t, t, 1.0);

    return x.dot(c);
}

Complex Interpolation::linear3D(const Image<Complex>& img, double x, double y, double z)
{
    if (!(x >= 0.0 && x < img.data.xdim-1 && y >= 0.0 && y < img.data.ydim-1 && z >= 0.0 && z < img.data.zdim-1))
    {
        return 0.0;
    }

    const int xi = (int)x;
    const int yi = (int)y;
    const int zi = (int)z;

    const double xf = x - xi;
    const double yf = y - yi;
    const double zf = z - zi;

    const Complex f000 = DIRECT_NZYX_ELEM(img.data, 0, zi, yi,   xi);
    const Complex f001 = DIRECT_NZYX_ELEM(img.data, 0, zi, yi,   xi+1);

    const Complex f010 = DIRECT_NZYX_ELEM(img.data, 0, zi, yi+1, xi);
    const Complex f011 = DIRECT_NZYX_ELEM(img.data, 0, zi, yi+1, xi+1);

    const Complex f100 = DIRECT_NZYX_ELEM(img.data, 0, zi+1, yi, xi);
    const Complex f101 = DIRECT_NZYX_ELEM(img.data, 0, zi+1, yi, xi+1);

    const Complex f110 = DIRECT_NZYX_ELEM(img.data, 0, zi+1, yi+1, xi);
    const Complex f111 = DIRECT_NZYX_ELEM(img.data, 0, zi+1, yi+1, xi+1);

    const Complex f00 = xf * f001 + (1.0 - xf) * f000;
    const Complex f01 = xf * f011 + (1.0 - xf) * f010;
    const Complex f10 = xf * f101 + (1.0 - xf) * f100;
    const Complex f11 = xf * f111 + (1.0 - xf) * f110;

    const Complex f0 = yf * f01 + (1.0 - yf) * f00;
    const Complex f1 = yf * f11 + (1.0 - yf) * f10;

    const Complex f = zf * f1 + (1.0 - zf) * f0;

    return f;

}

Complex Interpolation::linearFFTW3D(const Image<Complex>& img, double x, double y, double z)
{
    if (x > img.data.xdim-1)
    {
        return 0.0;
    }

    const int xi = (int)x;
    const int yi = (int)y;
    const int zi = (int)z;

    const int xp = xi + 1;
    const int yp = (yi+1)%((int)img.data.ydim);
    const int zp = (zi+1)%((int)img.data.zdim);

    const double xf = x - xi;
    const double yf = y - yi;
    const double zf = z - zi;

    const Complex f000 = DIRECT_NZYX_ELEM(img.data, 0, zi, yi, xi);
    const Complex f001 = DIRECT_NZYX_ELEM(img.data, 0, zi, yi, xp);

    const Complex f010 = DIRECT_NZYX_ELEM(img.data, 0, zi, yp, xi);
    const Complex f011 = DIRECT_NZYX_ELEM(img.data, 0, zi, yp, xp);

    const Complex f100 = DIRECT_NZYX_ELEM(img.data, 0, zp, yi, xi);
    const Complex f101 = DIRECT_NZYX_ELEM(img.data, 0, zp, yi, xp);

    const Complex f110 = DIRECT_NZYX_ELEM(img.data, 0, zp, yp, xi);
    const Complex f111 = DIRECT_NZYX_ELEM(img.data, 0, zp, yp, xp);

    const Complex f00 = xf * f001 + (1.0 - xf) * f000;
    const Complex f01 = xf * f011 + (1.0 - xf) * f010;
    const Complex f10 = xf * f101 + (1.0 - xf) * f100;
    const Complex f11 = xf * f111 + (1.0 - xf) * f110;

    const Complex f0 = yf * f01 + (1.0 - yf) * f00;
    const Complex f1 = yf * f11 + (1.0 - yf) * f10;

    const Complex f = zf * f1 + (1.0 - zf) * f0;

    return f;
}

Complex Interpolation::linearFFTW2D(const Image<Complex>& img, double x, double y)
{
    if (x > img.data.xdim-1)
    {
        return 0.0;
    }

    const int xi = (int)x;
    const int yi = (int)y;

    const int xp = xi + 1;
    const int yp = (yi+1)%((int)img.data.ydim);

    const double xf = x - xi;
    const double yf = y - yi;

    const Complex f00 = DIRECT_NZYX_ELEM(img.data, 0, 0, yi, xi);
    const Complex f01 = DIRECT_NZYX_ELEM(img.data, 0, 0, yi, xp);

    const Complex f10 = DIRECT_NZYX_ELEM(img.data, 0, 0, yp, xi);
    const Complex f11 = DIRECT_NZYX_ELEM(img.data, 0, 0, yp, xp);

    const Complex f0 = xf * f01 + (1.0 - xf) * f00;
    const Complex f1 = xf * f11 + (1.0 - xf) * f10;

    const Complex f = yf * f1 + (1.0 - yf) * f0;

    return f;
}

void Interpolation::test2D()
{
    int w0 = 5, w1 = 1000;
    Image<RFLOAT> 
		img0(w0,w0), img(w1,w1), img1(w1,w1), 
		img2a(w1,w1), img2b(w1,w1), 
		img3a(w1,w1), img3b(w1,w1), img3c(w1,w1);
	
    Image<RFLOAT> gradx(w1,w1), grady(w1,w1);
    Image<RFLOAT> gradnx(w1,w1), gradny(w1,w1);

    for (int y = 0; y < w0; y++)
    for (int x = 0; x < w0; x++)
    {
        DIRECT_NZYX_ELEM(img0.data, 0, 0, y, x) = (x + w0*y)*(1 - 2*((x%2)^(y%2)));
    }

    double eps = 0.001;

    for (int y = 0; y < w1; y++)
    for (int x = 0; x < w1; x++)
    {
        DIRECT_NZYX_ELEM(img.data, 0, 0, y, x)  = DIRECT_NZYX_ELEM(img0.data, 0, 0, w0*y/w1, w0*x/w1);
		
		DIRECT_NZYX_ELEM(img1.data, 0, 0, y, x) = cubicXY(
					img0, w0*x/(double)w1 - 0.5, w0*y/(double)w1 - 0.5, 0);
		
		DIRECT_NZYX_ELEM(img2a.data, 0, 0, y, x) = cubicXY(
					img0, w0*x/(double)w1 - 0.5, w0*y/(double)w1 - 0.5, 0, 0, true);
		
		DIRECT_NZYX_ELEM(img2b.data, 0, 0, y, x) = cubicXY(
					img0, w0*x/(double)w1 - 0.5 - w0/2, w0*y/(double)w1 - 0.5 - w0/2, 0, 0, true);
		
		DIRECT_NZYX_ELEM(img3a.data, 0, 0, y, x) = cubicXY(
					img0, 1e-2*(x-w1/2), 1e-2*(y-w1/2), 0, 0, true);
		
		DIRECT_NZYX_ELEM(img3b.data, 0, 0, y, x) = cubicXY(
					img0, 1e-7*(x-w1/2), 1e-7*(y-w1/2), 0, 0, true);
		
		DIRECT_NZYX_ELEM(img3c.data, 0, 0, y, x) = cubicXY(
					img0, 1e-16*(x-w1/2), 1e-16*(y-w1/2), 0, 0, true);
        
        t2Vector<RFLOAT> g = cubicXYgrad(img0, w0*x/(double)w1 - 0.5, w0*y/(double)w1 - 0.5, 0);

        DIRECT_NZYX_ELEM(gradx.data, 0, 0, y, x) = g.x;
        DIRECT_NZYX_ELEM(grady.data, 0, 0, y, x) = g.y;

        DIRECT_NZYX_ELEM(gradnx.data, 0, 0, y, x) =
               (  cubicXY(img0, w0*x/(double)w1 - 0.5 + eps, w0*y/(double)w1 - 0.5, 0)
                - cubicXY(img0, w0*x/(double)w1 - 0.5 - eps, w0*y/(double)w1 - 0.5, 0))/(2.0*eps);

        DIRECT_NZYX_ELEM(gradny.data, 0, 0, y, x) =
               (  cubicXY(img0, w0*x/(double)w1 - 0.5, w0*y/(double)w1 - 0.5 + eps, 0)
                - cubicXY(img0, w0*x/(double)w1 - 0.5, w0*y/(double)w1 - 0.5 - eps, 0))/(2.0*eps);
    }

    VtkHelper::writeVTK(img,  "debug/interpolationX_0.vtk", 0, 0, 0, 1, 1, 1);
	VtkHelper::writeVTK(img1, "debug/interpolationX_1.vtk", 0, 0, 0, 1, 1, 1);
	VtkHelper::writeVTK(img2a, "debug/interpolationX_1w.vtk", 0, 0, 0, 1, 1, 1);
	VtkHelper::writeVTK(img2b, "debug/interpolationX_2w.vtk", 0, 0, 0, 1, 1, 1);
	
	VtkHelper::writeVTK(img3a, "debug/interpolationX_3w_0.vtk", 0, 0, 0, 1, 1, 1);
	VtkHelper::writeVTK(img3b, "debug/interpolationX_3w_4.vtk", 0, 0, 0, 1, 1, 1);
	VtkHelper::writeVTK(img3c, "debug/interpolationX_3w_16.vtk", 0, 0, 0, 1, 1, 1);
	
    VtkHelper::writeVTK(gradx, "debug/interpolationX_gx.vtk", 0, 0, 0, 1, 1, 1);
    VtkHelper::writeVTK(grady, "debug/interpolationX_gy.vtk", 0, 0, 0, 1, 1, 1);
    VtkHelper::writeVTK(gradnx, "debug/interpolationX_gnx.vtk", 0, 0, 0, 1, 1, 1);
    VtkHelper::writeVTK(gradny, "debug/interpolationX_gny.vtk", 0, 0, 0, 1, 1, 1);
}

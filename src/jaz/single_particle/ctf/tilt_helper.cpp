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

#include "tilt_helper.h"

#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/optimization/nelder_mead.h>

#include <src/jaz/math/Zernike.h>

using namespace gravis;

void TiltHelper::updateTiltShift(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        CTF& ctf, double angpix,
        Image<Complex>& xyDest,
        Image<RFLOAT>& wDest,
		bool do_ctf_padding)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix, false, false, false, true, do_ctf_padding);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);

        RFLOAT c = ctfImg(y,x);

        DIRECT_A2D_ELEM(xyDest.data, y, x) += c * vx.conj() * vy;
        DIRECT_A2D_ELEM(wDest.data, y, x) += c * c * vx.norm();
    }
}

void TiltHelper::updateTiltShiftPar(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        CTF& ctf, double angpix,
        Image<Complex>& xyDest,
        Image<RFLOAT>& wDest,
		bool do_ctf_padding)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix, false, false, false, true, do_ctf_padding);

    #pragma omp parallel for
    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);

        RFLOAT c = ctfImg(y,x);

        DIRECT_A2D_ELEM(xyDest.data, y, x) += c * vx.conj() * vy;
        DIRECT_A2D_ELEM(wDest.data, y, x) += c * c * vx.norm();
    }
}

void TiltHelper::fitTiltShift(
		const Image<RFLOAT>& phase,
		const Image<RFLOAT>& weight,
		double Cs, double lambda,
		double angpix, const Matrix2D<RFLOAT>& mag,
		double* shift_x, double* shift_y,
		double* tilt_x, double* tilt_y,
		Image<RFLOAT>* fit)
{
    const long w = phase.data.xdim;
    const long h = phase.data.ydim;

    double axx = 0.0, axy = 0.0, axz = 0.0,//axw == ayz,
                      ayy = 0.0, ayz = 0.0, ayw = 0.0,
                                 azz = 0.0, azw = 0.0,
                                            aww = 0.0;

    double bx = 0.0, by = 0.0, bz = 0.0, bw = 0.0;

    const RFLOAT as = (RFLOAT)h * angpix;

    for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++)
    {
        const double x0 = xi;
        const double y0 = yi < w? yi : ((yi-h));

		const double x = mag(0,0) * x0 + mag(0,1) * y0;
		const double y = mag(1,0) * x0 + mag(1,1) * y0;

        double q = x*x + y*y;

        double v = DIRECT_A2D_ELEM(phase.data, yi, xi);
        double g = DIRECT_A2D_ELEM(weight.data, yi, xi);

        axx += g     * x * x;
        axy += g     * x * y;
        axz += g * q * x * x;

        ayy += g     * y * y;
        ayz += g * q * x * y;
        ayw += g * q * y * y;

        azz += g * q * q * x * x;
        azw += g * q * q * x * y;

        aww += g * q * q * y * y;

        bx += g * x * v;
        by += g * y * v;
        bz += g * q * x * v;
        bw += g * q * y * v;
    }

    gravis::d4Matrix A;
    gravis::d4Vector b(bx, by, bz, bw);

    A(0,0) = axx;
    A(0,1) = axy;
    A(0,2) = axz;
    A(0,3) = ayz;

    A(1,0) = axy;
    A(1,1) = ayy;
    A(1,2) = ayz;
    A(1,3) = ayw;

    A(2,0) = axz;
    A(2,1) = ayz;
    A(2,2) = azz;
    A(2,3) = azw;

    A(3,0) = ayz;
    A(3,1) = ayw;
    A(3,2) = azw;
    A(3,3) = aww;

    gravis::d4Matrix Ainv = A;
    Ainv.invert();

    gravis::d4Vector opt = Ainv * b;

    //std::cout << opt[0] << ", " << opt[1] << ", " << opt[2] << ", " << opt[3] << "\n";

    *shift_x = opt[0];
    *shift_y = opt[1];

    *tilt_x = -opt[2]*180.0/(0.360 * Cs * 10000000 * lambda * lambda * 3.141592654);
    *tilt_y = -opt[3]*180.0/(0.360 * Cs * 10000000 * lambda * lambda * 3.141592654);

    /*destA = gravis::d2Vector(opt[0], opt[1]).length();
    destPhiA = (180.0/3.1416)*std::atan2(opt[1], opt[0]);

    destB = gravis::d2Vector(opt[2], opt[3]).length();
    destPhiB = (180.0/3.1416)*std::atan2(opt[3], opt[2]);

    std::cout << "linear: " << destA << " @ " << destPhiA << "°\n";
    std::cout << "cubic:  " << destB << " @ " << destPhiB << "°\n";
    std::cout << "    =  -" << destB << " @ " << (destPhiB + 180.0) << "°\n";*/

    //std::cout << "tilt_x = " << *tilt_x << "\n";
    //std::cout << "tilt_y = " << *tilt_y << "\n";

    if (fit != 0)
    {
        *fit = Image<RFLOAT>(w,h);
        drawPhaseShift(opt[0], opt[1], opt[2], opt[3], w, h, as, mag, fit);
    }
}

void TiltHelper::optimizeTilt(
    const Image<Complex> &xy,
    const Image<RFLOAT> &weight,
    double Cs, double lambda,
	double angpix, const Matrix2D<RFLOAT>& mag,
    bool L1,
    double shift0_x, double shift0_y,
    double tilt0_x, double tilt0_y,
    double *shift_x, double *shift_y,
    double *tilt_x, double *tilt_y,
    Image<RFLOAT> *fit)
{
    TiltOptimization prob(xy, weight, angpix, mag, L1, false);

    double scale = 180.0/(0.360 * Cs * 10000000 * lambda * lambda * 3.141592654);

    std::vector<double> initial{shift0_x, shift0_y, -tilt0_x/scale, -tilt0_y/scale};

    std::vector<double> opt = NelderMead::optimize(
                initial, prob, 0.01, 0.000001, 100000,
                1.0, 2.0, 0.5, 0.5, false);

    *shift_x = opt[0];
    *shift_y = opt[1];

    *tilt_x = -opt[2]*scale;
    *tilt_y = -opt[3]*scale;

    //std::cout << opt[0] << ", " << opt[1] << ", " << opt[2] << ", " << opt[3] << "\n";
    //std::cout << "tilt_x = " << *tilt_x << "\n";
    //std::cout << "tilt_y = " << *tilt_y << "\n";

    if (fit != 0)
    {
        const int w = xy.data.xdim;
        const int h = xy.data.ydim;
        const RFLOAT as = (RFLOAT)h * angpix;

        *fit = Image<RFLOAT>(w,h);
        drawPhaseShift(opt[0], opt[1], opt[2], opt[3], w, h, as, mag, fit);
	}
}

std::vector<double> TiltHelper::fitOddZernike(
	const Image<Complex>& xy,
	const Image<RFLOAT>& weight,
	double angpix, const Matrix2D<RFLOAT>& mag,
	int n_max,
	Image<RFLOAT>* fit)
{
	const int w = xy.data.xdim;
	const int h = xy.data.ydim;
	const int cc = Zernike::numberOfOddCoeffs(n_max);

	std::vector<Image<RFLOAT>> basis = TiltHelper::computeOddZernike(h, angpix, mag, n_max);

	std::vector<double> out = fitBasisLin(xy, weight, basis);

	if (fit != 0)
	{
		*fit = Image<RFLOAT>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			for (int c = 0; c < cc; c++)
			{
				(*fit)(y,x) += out[c] * basis[c](y,x);
			}
		}
	}

	return out;
}

std::vector<double> TiltHelper::optimiseOddZernike(
	const Image<Complex> &xy,
	const Image<RFLOAT> &weight,
	double angpix, const Matrix2D<RFLOAT>& mag, int n_max,
	const std::vector<double> &coeffs,
	Image<RFLOAT> *fit)
{
	const int w = xy.data.xdim;
	const int h = xy.data.ydim;
	const int cc = Zernike::numberOfOddCoeffs(n_max);

	std::vector<Image<RFLOAT>> basis = computeOddZernike(h, angpix, mag, n_max);

	std::vector<double> opt = optimiseBasis(xy, weight, basis, coeffs);

	if (fit != 0)
	{
		*fit = Image<RFLOAT>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			for (int c = 0; c < cc; c++)
			{
				(*fit)(y,x) += opt[c] * basis[c](y,x);
			}
		}
	}

	return opt;
}

std::vector<Image<RFLOAT> > TiltHelper::computeOddZernike(
		int s, double angpix, const Matrix2D<RFLOAT>& mag, int n_max)
{
	const int cc = Zernike::numberOfOddCoeffs(n_max);
	const double as = (double)s * angpix;
	const int sh = s/2 + 1;

	std::vector<Image<RFLOAT>> basis(cc);

	for (int c = 0; c < cc; c++)
	{
		basis[c] = Image<RFLOAT>(sh,s);

		int m, n;

		Zernike::oddIndexToMN(c, m, n);

		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			const double xx0 = x/as;
			const double yy0 = y < sh-1? y/as : (y-s)/as;

			const double xx = mag(0,0) * xx0 + mag(0,1) * yy0;
			const double yy = mag(1,0) * xx0 + mag(1,1) * yy0;

			basis[c](y,x) = Zernike::Z_cart(m, n, xx, yy);
		}
	}

	return basis;
}

Image<RFLOAT> TiltHelper::plotOddZernike(
		const std::vector<double>& coeffs,
		int s, double angpix, const Matrix2D<RFLOAT>& mag)
{
	Image<RFLOAT> out(s,s);

	const double as = (double)s * angpix;

	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		const double xx0 = (x - s/2) / as;
		const double yy0 = (y - s/2) / as;

		const double xx = mag(0,0) * xx0 + mag(0,1) * yy0;
		const double yy = mag(1,0) * xx0 + mag(1,1) * yy0;

		for (int c = 0; c < coeffs.size(); c++)
		{
			int m, n;
			Zernike::oddIndexToMN(c, m, n);

			out(y,x) += coeffs[c] * Zernike::Z_cart(m, n, xx, yy);
		}
	}

	return out;
}

Image<RFLOAT> TiltHelper::plotTilt(
		double tx, double ty, int s,
		double angpix, const Matrix2D<RFLOAT>& mag,
		double Cs, double lambda)
{
	Image<RFLOAT> out(s,s);

	/*double boxsize = angpix * s;
    double factor = 0.360 * Cs * 10000000 * lambda * lambda / (boxsize * boxsize * boxsize);*/

	const double scale = Cs * 20000 * lambda * lambda * 3.141592654;
	const double as = (double)s * angpix;

	for (int y = 0; y < s; y++)
	for (int x = 0; x < s; x++)
	{
		/*const double xx = (x - s/2);
		const double yy = (y - s/2);

		out(y,x) = factor * (yy * yy + xx * xx) * (yy * ty + xx * tx);*/

		const double xx0 = (x - s/2) / as;
		const double yy0 = (y - s/2) / as;

		const double xx = mag(0,0) * xx0 + mag(0,1) * yy0;
		const double yy = mag(1,0) * xx0 + mag(1,1) * yy0;

		out(y,x) = scale * (xx*xx + yy*yy) * (xx*tx + yy*ty);
	}

	return out;
}


std::vector<double> TiltHelper::fitEvenZernike(
		const Image<RFLOAT>& phase,
		const Image<RFLOAT>& weight,
		double angpix, const Matrix2D<RFLOAT>& mag, int n_max,
		Image<RFLOAT>* fit)
{
	const int w = phase.data.xdim;
	const int h = phase.data.ydim;
	const int cc = Zernike::numberOfEvenCoeffs(n_max);

	std::vector<Image<RFLOAT>> basis = computeEvenZernike(h, angpix, mag, n_max);

	std::vector<double> out = fitBasisLin(phase, weight, basis);

	if (fit != 0)
	{
		*fit = Image<RFLOAT>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			for (int c = 0; c < cc; c++)
			{
				(*fit)(y,x) += out[c] * basis[c](y,x);
			}
		}
	}

	return out;
}

std::vector<double> TiltHelper::optimiseEvenZernike(
	const Image<Complex>& xy,
	const Image<RFLOAT>& weight,
	double angpix, const Matrix2D<RFLOAT>& mag, int n_max,
	const std::vector<double>& coeffs,
	Image<RFLOAT>* fit)
{
	const int w = xy.data.xdim;
	const int h = xy.data.ydim;
	const int cc = Zernike::numberOfEvenCoeffs(n_max);

	std::vector<Image<RFLOAT>> basis = computeEvenZernike(h, angpix, mag, n_max);

	std::vector<double> opt = optimiseBasis(xy, weight, basis, coeffs);

	if (fit != 0)
	{
		*fit = Image<RFLOAT>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			for (int c = 0; c < cc; c++)
			{
				(*fit)(y,x) += opt[c] * basis[c](y,x);
			}
		}
	}

	return opt;
}

std::vector<double> TiltHelper::optimiseEvenZernike(
	const Image<Complex>& xy,
	const Image<RFLOAT>& weight0,
	const Image<RFLOAT>& Axx,
	const Image<RFLOAT>& Axy,
	const Image<RFLOAT>& Ayy,
	double angpix, const Matrix2D<RFLOAT>& mag, int n_max,
	const std::vector<double>& coeffs,
	Image<RFLOAT>* fit)
{
	const int w = xy.data.xdim;
	const int h = xy.data.ydim;
	const int cc = Zernike::numberOfEvenCoeffs(n_max);

	std::vector<Image<RFLOAT>> basis = computeEvenZernike(h, angpix, mag, n_max);

	std::vector<double> opt = optimiseBasis(xy, weight0, Axx, Axy, Ayy, basis, coeffs);

	if (fit != 0)
	{
		*fit = Image<RFLOAT>(w,h);

		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			for (int c = 0; c < cc; c++)
			{
				(*fit)(y,x) += opt[c] * basis[c](y,x);
			}
		}
	}

	return opt;
}

std::vector<Image<RFLOAT> > TiltHelper::computeEvenZernike(
		int s, double angpix, const Matrix2D<RFLOAT>& mag, int n_max)
{
	const int cc = Zernike::numberOfEvenCoeffs(n_max);
	const double as = (double)s * angpix;
	const int sh = s/2 + 1;

	std::vector<Image<RFLOAT>> basis(cc);

	for (int c = 0; c < cc; c++)
	{
		basis[c] = Image<RFLOAT>(sh,s);

		int m, n;

		Zernike::evenIndexToMN(c, m, n);

		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			const double xx0 = x/as;
			const double yy0 = y < sh-1? y/as : (y-s)/as;

			const double xx = mag(0,0) * xx0 + mag(0,1) * yy0;
			const double yy = mag(1,0) * xx0 + mag(1,1) * yy0;

			basis[c](y,x) = Zernike::Z_cart(m, n, xx, yy);
		}
	}

	return basis;
}

void TiltHelper::extractTilt(
	std::vector<double>& oddZernikeCoeffs,
	double& tilt_x, double& tilt_y,
	double Cs, double lambda)
{
	if (oddZernikeCoeffs.size() <= 5)
		oddZernikeCoeffs.resize(5, 0);

	const double scale = Cs * 20000 * lambda * lambda * 3.141592654;

	const double Z3x = oddZernikeCoeffs[4];
	const double Z3y = oddZernikeCoeffs[3];

	// p = Z1x x + Z3x (3r² - 2) x
	//   = (Z1x - 2 Z3x) x + 3 Z3x r² x
	//   = Z1x' x - tx r² x

	oddZernikeCoeffs[4] = 0.0;
	oddZernikeCoeffs[3] = 0.0;

	oddZernikeCoeffs[1] -= 2.0 * Z3x;
	oddZernikeCoeffs[0] -= 2.0 * Z3y;

	tilt_x = -3.0 * Z3x / scale;
	tilt_y = -3.0 * Z3y / scale;
}

void TiltHelper::insertTilt(
	std::vector<double>& oddZernikeCoeffs,
	double tilt_x, double tilt_y,
	double Cs, double lambda)
{
	if (oddZernikeCoeffs.size() <= 5)
		oddZernikeCoeffs.resize(5, 0);

	const double scale = Cs * 20000 * lambda * lambda * 3.141592654;
	const double Z3x = -scale * tilt_x / 3.0;
	const double Z3y = -scale * tilt_y / 3.0;

	oddZernikeCoeffs[1] += 2.0 * Z3x;
	oddZernikeCoeffs[0] += 2.0 * Z3y;

	oddZernikeCoeffs[4] += Z3x;
	oddZernikeCoeffs[3] += Z3y;
}

std::vector<double> TiltHelper::fitBasisLin(
	const Image<Complex>& xy,
	const Image<RFLOAT>& weight,
	const std::vector<Image<RFLOAT>>& basis)
{
	const int w = xy.data.xdim;
	const int h = xy.data.ydim;

	Image<RFLOAT> phase(w,h);

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		phase(y,x) = xy(y,x).arg();
	}

	return fitBasisLin(phase, weight, basis);
}

std::vector<double> TiltHelper::fitBasisLin(
	const Image<RFLOAT>& phase,
	const Image<RFLOAT>& weight,
	const std::vector<Image<RFLOAT>>& basis)
{
	const int cc = basis.size();
	const int w = phase.data.xdim;
	const int h = phase.data.ydim;

	Matrix2D<RFLOAT> A(cc,cc);
	Matrix1D<RFLOAT> b(cc);

	for (int c1 = 0; c1 < cc; c1++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			b(c1) += weight(y,x) * basis[c1](y,x) * phase(y,x);
		}

		for (int c2 = c1; c2 < cc; c2++)
		{
			for (int y = 0; y < h; y++)
			for (int x = 0; x < w; x++)
			{
				A(c1,c2) += weight(y,x) * basis[c1](y,x) * basis[c2](y,x);
			}
		}

		for (int c2 = 0; c2 < c1; c2++)
		{
			A(c1,c2) = A(c2,c1);
		}
	}

	const double tol = 1e-20;
	Matrix1D<RFLOAT> x(cc);
	solve(A, b, x, tol);

	std::vector<double> out(cc);

	for (int c = 0; c < cc; c++)
	{
		out[c] = x(c);
	}

	return out;
}

std::vector<double> TiltHelper::optimiseBasis(
		const Image<Complex>& xy,
		const Image<RFLOAT>& weight,
		const std::vector<Image<RFLOAT>>& basis,
		const std::vector<double>& initial)
{
	BasisOptimisation prob(xy, weight, basis, false);

	std::vector<double> opt = NelderMead::optimize(
		initial, prob, 0.01, 0.000001, 100000, 1.0, 2.0, 0.5, 0.5, false);

	return opt;
}

std::vector<double> TiltHelper::optimiseBasis(
		const Image<Complex>& xy,
		const Image<RFLOAT>& weight0,
		const Image<RFLOAT>& Axx,
		const Image<RFLOAT>& Axy,
		const Image<RFLOAT>& Ayy,
		const std::vector<Image<RFLOAT>>& basis,
		const std::vector<double>& initial)
{
	AnisoBasisOptimisation prob(xy, weight0, Axx, Axy, Ayy, basis, false);

	std::vector<double> opt = NelderMead::optimize(
		initial, prob, 0.01, 0.000001, 100000, 1.0, 2.0, 0.5, 0.5, false);

	return opt;
}

void TiltHelper::optimizeAnisoTilt(
    const Image<Complex> &xy,
    const Image<RFLOAT> &weight,
    double Cs, double lambda,
	double angpix, const Matrix2D<RFLOAT>& mag,
    bool L1,
    double shift0_x, double shift0_y,
    double tilt0_x, double tilt0_y,
    double *shift_x, double *shift_y,
    double *tilt_x, double *tilt_y,
    double* tilt_xx, double* tilt_xy, double* tilt_yy,
    Image<RFLOAT> *fit)
{
    TiltOptimization prob(xy, weight, angpix, mag, L1, true);

    double scale = 180.0/(0.360 * Cs * 10000000 * lambda * lambda * 3.141592654);

    std::vector<double> initial{shift0_x, shift0_y, -tilt0_x/scale, -tilt0_y/scale, 0.0, 1.0};

    std::vector<double> opt = NelderMead::optimize(
                initial, prob, 0.01, 0.000001, 100000,
                1.0, 2.0, 0.5, 0.5, false);

    *shift_x = opt[0];
    *shift_y = opt[1];

    *tilt_x = -opt[2]*scale;
    *tilt_y = -opt[3]*scale;

    *tilt_xx = 1.0;
    *tilt_xy = opt[4];
    *tilt_yy = opt[5];

    std::cout << opt[0] << ", " << opt[1] << ", " << opt[2]
        << ", " << opt[3] << ", " << opt[4] << ", " << opt[5] << "\n";

    std::cout << "tilt_x = " << *tilt_x << "\n";
    std::cout << "tilt_y = " << *tilt_y << "\n";
    std::cout << "tilt_xx = " << *tilt_xx << "\n";
    std::cout << "tilt_xy = " << *tilt_xy << "\n";
    std::cout << "tilt_yy = " << *tilt_yy << "\n";

    if (fit != 0)
    {
        const int w = xy.data.xdim;
        const int h = xy.data.ydim;
        const RFLOAT as = (RFLOAT)h * angpix;

        *fit = Image<RFLOAT>(w,h);
        drawPhaseShift(opt[0], opt[1], opt[2], opt[3], 1.0, opt[4], opt[5], w, h, as, mag, fit);
    }
}

void TiltHelper::drawPhaseShift(
        double shift_x, double shift_y,
        double tilt_x, double tilt_y,
        int w, int h, double as,
        const Matrix2D<RFLOAT>& mag,
        Image<RFLOAT> *tgt)
{
	for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++)
    {
        const double x0 = xi;
        const double y0 = yi < w? yi : yi-h;

		const double x = (mag(0,0) * x0 + mag(0,1) * y0) / as;
		const double y = (mag(1,0) * x0 + mag(1,1) * y0) / as;

		const double q = x*x + y*y;

        DIRECT_A2D_ELEM(tgt->data, yi, xi) = x * shift_x + y * shift_y + q * x * tilt_x + q * y * tilt_y;
    }
}

void TiltHelper::drawPhaseShift(
        double shift_x, double shift_y,
        double tilt_x, double tilt_y,
        double tilt_xx, double tilt_xy, double tilt_yy,
        int w, int h, double as,
        const Matrix2D<RFLOAT>& mag,
        Image<RFLOAT> *tgt)
{
    for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++)
    {
		const double x0 = xi;
        const double y0 = yi < w? yi : yi-h;

        const double x = (mag(0,0) * x0 + mag(0,1) * y0) / as;
		const double y = (mag(1,0) * x0 + mag(1,1) * y0) / as;

        const double q = tilt_xx * x * x + 2.0 * tilt_xy * x * y + tilt_yy * y * y;

        DIRECT_A2D_ELEM(tgt->data, yi, xi) = x * shift_x + y * shift_y + q * x * tilt_x + q * y * tilt_y;
    }
}

TiltOptimization::TiltOptimization(
	const Image<Complex> &xy,
    const Image<RFLOAT> &weight,
	double angpix, const Matrix2D<RFLOAT> &mag,
	bool L1, bool anisotropic)
:
	xy(xy),
	weight(weight),
	angpix(angpix),
	mag(mag),
	L1(L1),
	anisotropic(anisotropic)
{
}

double TiltOptimization::f(const std::vector<double> &x, void* tempStorage) const
{
    double out = 0.0;

    const int w = xy.data.xdim;
    const int h = xy.data.ydim;

    const double as = (double)h * angpix;

    for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++)
    {
        const double xd0 = xi/as;
        const double yd0 = yi < w? yi/as : (yi-h)/as;
        double rr;

		const double xd = mag(0,0) * xd0 + mag(0,1) * yd0;
		const double yd = mag(1,0) * xd0 + mag(1,1) * yd0;

        if (anisotropic)
        {
            rr = xd*xd + 2.0*x[4]*xd*yd + x[5]*yd*yd;
        }
        else
        {
            rr = xd*xd + yd*yd;
        }

        const double phi = x[0]*xd + x[1]*yd + rr*x[2]*xd + rr*x[3]*yd;
        const double e = (xy(yi,xi) - Complex(cos(phi), sin(phi))).norm();

        if (L1)
        {
            out += weight(yi,xi) * sqrt(e);
        }
        else
        {
            out += weight(yi,xi) * e;
        }

        /*double d = phi - atan2(xy(yi,xi).imag, xy(yi,xi).real);
        out += weight(yi,xi) * d*d;*/
    }

    return out;
}


BasisOptimisation::BasisOptimisation(
		const Image<Complex> &xy,
		const Image<RFLOAT> &weight,
		const std::vector<Image<RFLOAT> > &basis,
		bool L1)
:	w(xy.data.xdim),
	h(xy.data.ydim),
	cc(basis.size()),
	xy(xy),
	weight(weight),
	basis(basis),
	L1(L1)
{
}

double BasisOptimisation::f(const std::vector<double> &x, void *tempStorage) const
{
	Image<RFLOAT>& recomb = *((Image<RFLOAT>*)tempStorage);
	recomb.data.initZeros();

	for (int c  = 0; c < cc; c++)
	for (int yp = 0; yp < h; yp++)
	for (int xp = 0; xp < w; xp++)
	{
		recomb(yp,xp) += x[c] * basis[c](yp,xp);
	}

	double sum = 0.0;

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		Complex zPred(cos(recomb(y,x)), sin(recomb(y,x)));
		sum += weight(y,x) * (zPred - xy(y,x)).norm();
	}

	return sum;
}

void *BasisOptimisation::allocateTempStorage() const
{
	return new Image<RFLOAT>(w,h);
}

void BasisOptimisation::deallocateTempStorage(void *ts) const
{	
	delete static_cast<Image<RFLOAT>*>(ts);
}

AnisoBasisOptimisation::AnisoBasisOptimisation(
		const Image<Complex> &xy,
		const Image<RFLOAT> &weight0,
		const Image<RFLOAT>& Axx,
		const Image<RFLOAT>& Axy,
		const Image<RFLOAT>& Ayy,
		const std::vector<Image<RFLOAT> > &basis,
		bool L1)
:	w(xy.data.xdim),
	h(xy.data.ydim),
	cc(basis.size()),
	xy(xy),
	weight0(weight0),
	Axx(Axx),
	Axy(Axy),
	Ayy(Ayy),
	basis(basis),
	L1(L1)
{
}

double AnisoBasisOptimisation::f(const std::vector<double> &x, void *tempStorage) const
{
	Image<RFLOAT>& recomb = *((Image<RFLOAT>*)tempStorage);
	recomb.data.initZeros();

	for (int c  = 0; c < cc; c++)
	for (int yp = 0; yp < h; yp++)
	for (int xp = 0; xp < w; xp++)
	{
		recomb(yp,xp) += x[c] * basis[c](yp,xp);
	}

	double sum = 0.0;

	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		d2Vector e(cos(recomb(y,x)) - xy(y,x).real, sin(recomb(y,x)) - xy(y,x).imag);

		sum += weight0(y,x) * (Axx(y,x)*e.x*e.x + 2.0*Axy(y,x)*e.x*e.y + Ayy(y,x)*e.y*e.y);
	}

	return sum;
}

void* AnisoBasisOptimisation::allocateTempStorage() const
{
	return new Image<RFLOAT>(w,h);
}

void AnisoBasisOptimisation::deallocateTempStorage(void *ts) const
{
	delete static_cast<Image<RFLOAT>*>(ts);
}

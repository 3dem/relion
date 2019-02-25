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

#include "magnification_helper.h"

#include <src/jaz/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>

using namespace gravis;

void MagnificationHelper::updateScaleFreq(
		const Image<Complex> &prediction, const Image<Complex> &observation,
		CTF &ctf, double angpix, Volume<Equation2x2> &eqs)
{
	const long w = prediction.data.xdim;
	const long h = prediction.data.ydim;
	
	const double as = (double)h * angpix;
	Volume<gravis::d2Vector> gradReal(w,h,1), gradImg(w,h,1);
	
	FilterHelper::centralGrad2D(prediction, gradReal, gradImg);
	
	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		const double xf = x;
		const double yf = y < w? y : y - h;
		
		Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
		Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);
		
		double c = ctf.getCTF(xf/as, yf/as);
		
		gravis::d2Vector gr = gradReal(x,y,0);
		gravis::d2Vector gi = gradImg(x,y,0);
		
		eqs(x,y,0).Axx += c * c * (gr.x * gr.x + gi.x * gi.x);
		eqs(x,y,0).Axy += c * c * (gr.x * gr.y + gi.x * gi.y);
		eqs(x,y,0).Ayy += c * c * (gr.y * gr.y + gi.y * gi.y);
		
		eqs(x,y,0).bx += c * (gr.x * (c * vx.real - vy.real) + gi.x * (c * vx.imag - vy.imag));
		eqs(x,y,0).by += c * (gr.y * (c * vx.real - vy.real) + gi.y * (c * vx.imag - vy.imag));
	}
}

void MagnificationHelper::updateScaleReal(
		const Image<Complex> &prediction,
		const Image<Complex> &observation,
		const Image<RFLOAT>& snr,
		CTF &ctf, double angpix,
		Volume<Equation2x2> &eqs)
{
	const long ww = 2*(observation.data.xdim - 1);
	const long w = prediction.data.xdim;
	const long h = prediction.data.ydim;
	
	const double as = (double)h * angpix;
	
	Image<Complex> pred2(w,h), obs2(w,h);
	Image<RFLOAT> realPred(ww, h), realObs(ww, h);
	
	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		const double xf = x;
		const double yf = y < w? y : y - h;
		
		Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
		Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);
		Complex sn = DIRECT_A2D_ELEM(snr.data, y, x);
		double c = ctf.getCTF(xf/as, yf/as);
		
		DIRECT_A2D_ELEM(pred2.data, y, x) = sn * c * vx;
		DIRECT_A2D_ELEM(obs2.data, y, x) = sn * vy;
	}
	
	FourierTransformer ft0, ft1;
	ft0.inverseFourierTransform(pred2.data, realPred.data);
	ft1.inverseFourierTransform(obs2.data, realObs.data);
	
	Volume<gravis::d2Vector> grad(ww,h,1);
	FilterHelper::centralGrad2D(realPred, grad);
	
	for (long y = 0; y < h; y++)
	for (long x = 0; x < ww; x++)
	{
		double vx = DIRECT_A2D_ELEM(realPred.data, y, x);
		double vy = DIRECT_A2D_ELEM(realObs.data, y, x);
		
		gravis::d2Vector g = grad(x,y,0);
		
		eqs(x,y,0).Axx += g.x * g.x;
		eqs(x,y,0).Axy += g.x * g.y;
		eqs(x,y,0).Ayy += g.y * g.y;
		
		eqs(x,y,0).bx += g.x * (vx - vy);
		eqs(x,y,0).by += g.y * (vx - vy);
	}
}

void MagnificationHelper::solvePerPixel(
		const Volume<Equation2x2> &eqs,
		Image<RFLOAT> &vx, Image<RFLOAT> &vy)
{
	const long w = eqs.dimx;
	const long h = eqs.dimy;
	
	vx = Image<RFLOAT>(w,h);
	vy = Image<RFLOAT>(w,h);
	
	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		Equation2x2 eq = eqs(x,y,0);
		
		gravis::d2Vector b(eq.bx, eq.by);
		gravis::d2Matrix A;
		A(0,0) = eq.Axx;
		A(0,1) = eq.Axy;
		A(1,0) = eq.Axy;
		A(1,1) = eq.Ayy;
		
		double det = A(0,0)*A(1,1) - A(0,1)*A(1,0);
		
		if (det == 0.0)
		{
			DIRECT_A2D_ELEM(vx.data, y, x) = 0.0;
			DIRECT_A2D_ELEM(vy.data, y, x) = 0.0;
		}
		else
		{
			gravis::d2Matrix Ai = A;
			Ai.invert();
			
			gravis::d2Vector xx = Ai * b;
			
			DIRECT_A2D_ELEM(vx.data, y, x) = xx.x;
			DIRECT_A2D_ELEM(vy.data, y, x) = xx.y;
		}
	}
}

void MagnificationHelper::solveLinearlyFreq(
		const Volume<Equation2x2> &eqs,
		const Image<RFLOAT>& snr,
		Image<RFLOAT> &mat,
		Image<RFLOAT> &vx, Image<RFLOAT> &vy)
{
	const long w = eqs.dimx;
	const long h = eqs.dimy;
	
	vx = Image<RFLOAT>(w,h);
	vy = Image<RFLOAT>(w,h);
	mat = Image<RFLOAT>(2,2);
	
	d4Vector b(0.0, 0.0, 0.0, 0.0);
	
	d4Matrix A(0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0);
	
	for (long yi = 1; yi < h-1; yi++)
	for (long xi = 1; xi < w-1; xi++)
	{
		Equation2x2 eq = eqs(xi,yi,0);
		double wgh = DIRECT_A2D_ELEM(snr.data, yi, xi);
		
		const double x = xi;
		const double y = yi < w? yi : (yi - h);
		
		A(0,0) += wgh * x * x * eq.Axx;
		A(0,1) += wgh * x * y * eq.Axx;
		A(0,2) += wgh * x * x * eq.Axy;
		A(0,3) += wgh * x * y * eq.Axy;
		
		A(1,0) += wgh * x * y * eq.Axx;
		A(1,1) += wgh * y * y * eq.Axx;
		A(1,2) += wgh * x * y * eq.Axy;
		A(1,3) += wgh * y * y * eq.Axy;
		
		A(2,0) += wgh * x * x * eq.Axy;
		A(2,1) += wgh * x * y * eq.Axy;
		A(2,2) += wgh * x * x * eq.Ayy;
		A(2,3) += wgh * x * y * eq.Ayy;
		
		A(3,0) += wgh * x * y * eq.Axy;
		A(3,1) += wgh * y * y * eq.Axy;
		A(3,2) += wgh * x * y * eq.Ayy;
		A(3,3) += wgh * y * y * eq.Ayy;
		
		b[0] += wgh * x * eq.bx;
		b[1] += wgh * y * eq.bx;
		b[2] += wgh * x * eq.by;
		b[3] += wgh * y * eq.by;
	}
	
	d4Matrix Ai = A;
	Ai.invert();
	
	d4Vector opt = Ai * b;
	
	DIRECT_A2D_ELEM(mat.data, 0, 0) = opt[0];
	DIRECT_A2D_ELEM(mat.data, 0, 1) = opt[1];
	DIRECT_A2D_ELEM(mat.data, 1, 0) = opt[2];
	DIRECT_A2D_ELEM(mat.data, 1, 1) = opt[3];
	
	std::cout << opt[0] << ", " << opt[1] << "\n" 
	          << opt[2] << ", " << opt[3] << "\n";
	
	for (long yi = 0; yi < h; yi++)
	for (long xi = 0; xi < w; xi++)
	{
		const double x = xi;
		const double y = yi < w? yi : (yi - h);
		
		DIRECT_A2D_ELEM(vx.data, yi, xi) = opt[0] * x + opt[1] * y;
		DIRECT_A2D_ELEM(vy.data, yi, xi) = opt[2] * x + opt[3] * y;
	}
}

void MagnificationHelper::readEQs(std::string path, Volume<Equation2x2> &eqs)
{
	Image<RFLOAT> Axx, Axy, Ayy, bx, by;
	
	Axx.read(path+"_Axx.mrc");
	Axy.read(path+"_Axy.mrc");
	Ayy.read(path+"_Ayy.mrc");
	bx.read(path+"_bx.mrc");
	by.read(path+"_by.mrc");
	
	const long w = Axx.data.xdim;
	const long h = Axx.data.ydim;
	
	eqs.resize(w,h,1);
	
	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		eqs(x,y,0).Axx = DIRECT_A2D_ELEM(Axx.data, y, x);
		eqs(x,y,0).Axy = DIRECT_A2D_ELEM(Axy.data, y, x);
		eqs(x,y,0).Ayy = DIRECT_A2D_ELEM(Ayy.data, y, x);
		
		eqs(x,y,0).bx = DIRECT_A2D_ELEM(bx.data, y, x);
		eqs(x,y,0).by = DIRECT_A2D_ELEM(by.data, y, x);
	}
}

void MagnificationHelper::writeEQs(const Volume<Equation2x2> &eqs, std::string path)
{
	const long w = eqs.dimx;
	const long h = eqs.dimy;
	
	Image<RFLOAT> Axx(w,h);
	Image<RFLOAT> Axy(w,h);
	Image<RFLOAT> Ayy(w,h);
	Image<RFLOAT> bx(w,h);
	Image<RFLOAT> by(w,h);
	
	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		Equation2x2 eq = eqs(x,y,0);
		
		DIRECT_A2D_ELEM(Axx.data, y, x) = eq.Axx;
		DIRECT_A2D_ELEM(Axy.data, y, x) = eq.Axy;
		DIRECT_A2D_ELEM(Ayy.data, y, x) = eq.Ayy;
		
		DIRECT_A2D_ELEM(bx.data, y, x) = eq.bx;
		DIRECT_A2D_ELEM(by.data, y, x) = eq.by;
	}
	
	Axx.write(path+"_Axx.mrc");
	Axy.write(path+"_Axy.mrc");
	Ayy.write(path+"_Ayy.mrc");
	bx.write(path+"_bx.mrc");
	by.write(path+"_by.mrc");
}

void MagnificationHelper::updatePowSpec(
		const Image<Complex> &prediction,
		const Image<Complex> &observation,
		CTF &ctf, double angpix,
		Image<RFLOAT> &powSpecPred,
		Image<RFLOAT> &powSpecObs)
{
	const long w = prediction.data.xdim;
	const long h = prediction.data.ydim;
	
	const double as = (double)h * angpix;
	
	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		const double xf = x;
		const double yf = y < w? y : y - h;
		
		Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
		Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);
		double c = ctf.getCTF(xf/as, yf/as);
		
		DIRECT_A2D_ELEM(powSpecPred.data, y, x) += (c * vx).abs();
		DIRECT_A2D_ELEM(powSpecObs.data, y, x) += vy.abs();
	}
}

Equation2x2::Equation2x2()
	: Axx(0.0), Axy(0.0), Ayy(0.0),
	  bx(0.0), by(0.0)
{
	
}

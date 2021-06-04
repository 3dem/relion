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

#include <src/jaz/single_particle/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/single_particle/d3x3/dsyev2.h>

using namespace gravis;

Matrix2D<RFLOAT> MagnificationHelper::polarToMatrix(
		double scaleMajor, double scaleMinor, double angleDeg)
{
	// based on definition by T. Nakane

	Matrix2D<RFLOAT> out(2,2);

	const double angle = DEG2RAD(angleDeg);
	const double si = sin(angle), co = cos(angle);
	const double si2 = si * si, co2 = co * co;

	/*
 	    Out = Rot(angle) * Diag(scale_major, scale_minor) * Rot(-angle), where
	    Rot(angle) = [[cos(angle), -sin(angle)], [sin(angle), cos(angle)]]

		[  c s ] [ j 0 ] [ c -s ]
		[ -s c ] [ 0 n ] [ s  c ]
		=
		[  c s ] [ jc -js ]
		[ -s c ] [ ns  nc ]
		=
		[  jcc+nss -jcs+ncs ]
		[ -jcs+ncs  jss+ncc ]
	*/

	out(0, 0) =  scaleMajor * co2 + scaleMinor * si2;
	out(1, 1) =  scaleMajor * si2 + scaleMinor * co2;
	out(0, 1) = (-scaleMajor + scaleMinor) * si * co;
	out(1, 0) = out(0, 1);

	return out;
}

void MagnificationHelper::matrixToPolar(
		const Matrix2D<RFLOAT>& mat,
		RFLOAT& scaleMajor,
		RFLOAT& scaleMinor,
		RFLOAT& angleDeg)
{
	matrixToPolar(
		d2Matrix(mat(0,0), mat(0,1), mat(1,0), mat(1,1)),
		scaleMajor, scaleMinor, angleDeg);
}

void MagnificationHelper::matrixToPolar(
		const d2Matrix& mat,
		RFLOAT& scaleMajor,
		RFLOAT& scaleMinor,
		RFLOAT& angleDeg)
{
	const double m00 = mat(0,0);
	const double m11 = mat(1,1);
	const double m01 = 0.5 * (mat(0,1) + mat(1,0));

	double ev0, ev1, cs, sn;

	dsyev2(m00, m01, m11, &ev0, &ev1, &cs, &sn);

	scaleMajor = ev0;
	scaleMinor = ev1;
	angleDeg = RAD2DEG(atan2(sn,cs));

	return;
}

// deprecated: use the other one!
void MagnificationHelper::updateScaleFreq(
		const Image<Complex> &prediction,
		const Volume<t2Vector<Complex>>& predGradient,
		const Image<Complex> &observation,
		CTF &ctf, double angpix, Volume<Equation2x2> &eqs,
		bool do_ctf_padding)
{
	const long w = prediction.data.xdim;
	const long h = prediction.data.ydim;

	/*Volume<gravis::d2Vector> gradReal(w,h,1), gradImg(w,h,1);

	FilterHelper::centralGrad2D(prediction, gradReal, gradImg);*/

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix, false, false, false, true, do_ctf_padding);

	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
		Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);

		double c = ctfImg(y,x);

		gravis::d2Vector gr(predGradient(x,y,0).x.real, predGradient(x,y,0).y.real);
		gravis::d2Vector gi(predGradient(x,y,0).x.imag, predGradient(x,y,0).y.imag);

		eqs(x,y,0).Axx += c * c * (gr.x * gr.x + gi.x * gi.x);
		eqs(x,y,0).Axy += c * c * (gr.x * gr.y + gi.x * gi.y);
		eqs(x,y,0).Ayy += c * c * (gr.y * gr.y + gi.y * gi.y);

		eqs(x,y,0).bx += c * (gr.x * (vy.real - c * vx.real) + gi.x * (vy.imag - c * vx.imag));
		eqs(x,y,0).by += c * (gr.y * (vy.real - c * vx.real) + gi.y * (vy.imag - c * vx.imag));
	}
}

void MagnificationHelper::updateScale(
		const RawImage<fComplex>& prediction,
		const RawImage<t2Vector<fComplex>>& predGradient,
		const RawImage<fComplex>& observation,
		const RawImage<float>& freqWeights,
		const RawImage<float>& doseWeights,
		CTF& ctf, double angpix,
		RawImage<Equation2x2>& eqs)
{
	const long w = prediction.xdim;
	const long h = prediction.ydim;

	BufferedImage<float> ctfImg(w,h);
	ctf.draw(h, h, angpix, 0, &ctfImg(0,0));

	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		const fComplex vx = prediction(x,y);
		const fComplex vy = observation(x,y);
		const double c = ctfImg(x,y) * doseWeights(x,y);
		const double g = freqWeights(x,y) * freqWeights(x,y);

		d2Vector gr(predGradient(x,y,0).x.real, predGradient(x,y,0).y.real);
		d2Vector gi(predGradient(x,y,0).x.imag, predGradient(x,y,0).y.imag);

		eqs(x,y).Axx += g * c * c * (gr.x * gr.x + gi.x * gi.x);
		eqs(x,y).Axy += g * c * c * (gr.x * gr.y + gi.x * gi.y);
		eqs(x,y).Ayy += g * c * c * (gr.y * gr.y + gi.y * gi.y);

		eqs(x,y).bx += g * c * (gr.x * (vy.real - c * vx.real) + gi.x * (vy.imag - c * vx.imag));
		eqs(x,y).by += g * c * (gr.y * (vy.real - c * vx.real) + gi.y * (vy.imag - c * vx.imag));
	}
}

void MagnificationHelper::updateScaleReal(
		const Image<Complex> &prediction,
		const Image<Complex> &observation,
		const Image<RFLOAT>& snr,
		CTF &ctf, double angpix,
		Volume<Equation2x2> &eqs,
		bool do_ctf_padding)
{
	const long ww = 2*(observation.data.xdim - 1);
	const long w = prediction.data.xdim;
	const long h = prediction.data.ydim;

	Image<Complex> pred2(w,h), obs2(w,h);
	Image<RFLOAT> realPred(ww, h), realObs(ww, h);

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix, false, false, false, true, do_ctf_padding);

	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
		Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);
		Complex sn = DIRECT_A2D_ELEM(snr.data, y, x);
		double c = ctfImg(y,x);

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

		eqs(x,y,0).bx += g.x * (vy - vx);
		eqs(x,y,0).by += g.y * (vy - vx);
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

BufferedImage<double> MagnificationHelper::solvePerPixel(
			const RawImage<Equation2x2>& eqs)
{
	const long w = eqs.xdim;
	const long h = eqs.ydim;

	BufferedImage<double> out(w,h,2);

	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		Equation2x2 eq = eqs(x,y);

		gravis::d2Vector b(eq.bx, eq.by);
		gravis::d2Matrix A;
		A(0,0) = eq.Axx;
		A(0,1) = eq.Axy;
		A(1,0) = eq.Axy;
		A(1,1) = eq.Ayy;

		double det = A(0,0)*A(1,1) - A(0,1)*A(1,0);

		if (det == 0.0)
		{
			out(x,y,0) = 0.0;
			out(x,y,1) = 0.0;
		}
		else
		{
			gravis::d2Matrix Ai = A;
			Ai.invert();

			gravis::d2Vector xx = Ai * b;

			out(x,y,0) = xx.x;
			out(x,y,1) = xx.y;
		}
	}

	return out;
}

Matrix2D<RFLOAT> MagnificationHelper::solveLinearlyFreq(
		const Volume<Equation2x2> &eqs,
		const Image<RFLOAT>& snr,
		Image<RFLOAT> &vx, Image<RFLOAT> &vy)
{
	Matrix2D<RFLOAT> mat(2,2);

	const long w = eqs.dimx;
	const long h = eqs.dimy;

	vx = Image<RFLOAT>(w,h);
	vy = Image<RFLOAT>(w,h);

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

	mat(0,0) = opt[0] + 1.0;
	mat(0,1) = opt[1];
	mat(1,0) = opt[2];
	mat(1,1) = opt[3] + 1.0;

//	std::cout << opt[0] << ", " << opt[1] << "\n"
//	          << opt[2] << ", " << opt[3] << "\n";

	for (long yi = 0; yi < h; yi++)
	for (long xi = 0; xi < w; xi++)
	{
		const double x = xi;
		const double y = yi < w? yi : (yi - h);

		DIRECT_A2D_ELEM(vx.data, yi, xi) = opt[0] * x + opt[1] * y;
		DIRECT_A2D_ELEM(vy.data, yi, xi) = opt[2] * x + opt[3] * y;
	}

	return mat;
}

d2Matrix MagnificationHelper::solveLinearly(
		const RawImage<Equation2x2>& eqs)
{
	d2Matrix mat(0.0, 0.0, 0.0, 0.0);

	const long w = eqs.xdim;
	const long h = eqs.ydim;


	d4Vector b(0.0, 0.0, 0.0, 0.0);

	d4Matrix A(0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0,
			   0.0, 0.0, 0.0, 0.0);

	for (long yi = 1; yi < h-1; yi++)
	for (long xi = 1; xi < w-1; xi++)
	{
		Equation2x2 eq = eqs(xi,yi,0);

		const double x = xi;
		const double y = yi < h/2? yi : (yi - h);

		const double wgh = (x*x + y*y < (w-1)*(w-1))? 1.0 : 0.0;

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

	mat(0,0) = opt[0] + 1.0;
	mat(0,1) = opt[1];
	mat(1,0) = opt[2];
	mat(1,1) = opt[3] + 1.0;

	return mat;
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
		Image<RFLOAT> &powSpecObs,
		bool do_ctf_padding)
{
	const long w = prediction.data.xdim;
	const long h = prediction.data.ydim;

	Image<RFLOAT> ctfImg(w,h);
	ctf.getFftwImage(ctfImg(), h, h, angpix, false, false, false, true, do_ctf_padding);

	for (long y = 0; y < h; y++)
	for (long x = 0; x < w; x++)
	{
		const double xf = x;
		const double yf = y < w? y : y - h;

		Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
		Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);
		double c = ctfImg(y,x);

		DIRECT_A2D_ELEM(powSpecPred.data, y, x) += (c * vx).abs();
		DIRECT_A2D_ELEM(powSpecObs.data, y, x) += vy.abs();
	}
}

void MagnificationHelper::adaptAstigmatism(
		const std::vector<Matrix2D<RFLOAT>>& dMs,
		std::vector<MetaDataTable>& partMdts,
		bool perParticle, ObservationModel* obsModel)
{
	const long int mc = partMdts.size();
	const int ogc = dMs.size();

	std::vector<d2Matrix> M(ogc), Mi(ogc), Mit(ogc);

	for (int og = 0; og < ogc; og++)
	{
		M[og] = d2Matrix(dMs[og](0,0), dMs[og](0,1), dMs[og](1,0), dMs[og](1,1));

		Mi[og] = M[og];
		Mi[og].invert();

		Mit[og] = Mi[og];
		Mit[og].transpose();
	}

	for (long int m = 0; m < mc; m++)
	{
		const int pc = partMdts[m].numberOfObjects();

		std::vector<d2Matrix> A(pc), D(pc), Q(pc);

		for (long int p = 0; p < pc; p++)
		{
			double deltafU, deltafV, phiDeg;

			partMdts[m].getValue(EMDL_CTF_DEFOCUSU, deltafU, p);
			partMdts[m].getValue(EMDL_CTF_DEFOCUSV, deltafV, p);
			partMdts[m].getValue(EMDL_CTF_DEFOCUS_ANGLE, phiDeg, p);

			const double phi = DEG2RAD(phiDeg);

			const double si = sin(phi);
			const double co = cos(phi);

			Q[p] = d2Matrix(co, si, -si, co);
			D[p] = d2Matrix(-deltafU, 0.0, 0.0, -deltafV);
			d2Matrix Qt(co, -si, si, co);

			A[p] = Qt * D[p] * Q[p];
		}

		if (perParticle)
		{
			for (long int p = 0; p < pc; p++)
			{
				int og;
				partMdts[m].getValue(EMDL_IMAGE_OPTICS_GROUP, og, p);
				og--;

				d2Matrix A2 = Mit[og] * A[p] * Mi[og];

				RFLOAT deltafU_neg, deltafV_neg, phiDeg;
				matrixToPolar(A2, deltafU_neg, deltafV_neg, phiDeg);

				partMdts[m].setValue(EMDL_CTF_DEFOCUSU, -deltafU_neg, p);
				partMdts[m].setValue(EMDL_CTF_DEFOCUSV, -deltafV_neg, p);
				partMdts[m].setValue(EMDL_CTF_DEFOCUS_ANGLE, phiDeg, p);
			}
		}
		else // keep difference between deltafU and deltafV, as well as the azimuth angle,
		     // constant for all particles in the same micrograph and optics group
		{
			std::vector<int> optGroups = obsModel->getOptGroupsPresent_oneBased(partMdts[m]);
			const int cc = optGroups.size();

			std::vector<int> groupToIndex(obsModel->numberOfOpticsGroups()+1, -1);

			for (int i = 0; i < cc; i++)
			{
				groupToIndex[optGroups[i]] = i;
			}

			for (int c = 0; c < cc; c++)
			{
				const int og = optGroups[c] - 1;

				d2Matrix A_mean(0.0, 0.0, 0.0, 0.0);

				for (long int p = 0; p < pc; p++)
				{
					int ogp;
					partMdts[m].getValue(EMDL_IMAGE_OPTICS_GROUP, ogp, p);
					ogp--;

					if (ogp == og)
					{
						A_mean += A[p] * (1.0/(double)pc);
					}
				}


				A_mean = Mit[og] * A_mean * Mi[og];

				double deltafU_mean_neg, deltafV_mean_neg, co, si;
				dsyev2(A_mean(0,0), A_mean(1,0), A_mean(1,1),
					   &deltafU_mean_neg, &deltafV_mean_neg, &co, &si);

				d2Matrix Q2(co, si, -si, co);
				d2Matrix Qt2(co, -si, si, co);

				double meanDef_mean = 0.5 * (deltafU_mean_neg + deltafV_mean_neg);

				for (long int p = 0; p < pc; p++)
				{
					int ogp;
					partMdts[m].getValue(EMDL_IMAGE_OPTICS_GROUP, ogp, p);
					ogp--;

					if (ogp == og)
					{
						d2Matrix Ap2 = Mit[og] * A[p] * Mi[og];

						double deltafU_p_neg, deltafV_p_neg, cop, sip;
						dsyev2(Ap2(0,0), Ap2(1,0), Ap2(1,1),
							   &deltafU_p_neg, &deltafV_p_neg, &cop, &sip);

						double meanDef_p = 0.5 * (deltafU_p_neg + deltafV_p_neg);

						d2Matrix Dp2(deltafU_mean_neg - meanDef_mean + meanDef_p, 0.0,
									 0.0, deltafV_mean_neg - meanDef_mean + meanDef_p);

						d2Matrix Apa2 = Qt2 * Dp2 * Q2;

						RFLOAT deltafU_pa_neg, deltafV_pa_neg, phiDeg;
						matrixToPolar(Apa2, deltafU_pa_neg, deltafV_pa_neg, phiDeg);

						partMdts[m].setValue(EMDL_CTF_DEFOCUSU, -deltafU_pa_neg, p);
						partMdts[m].setValue(EMDL_CTF_DEFOCUSV, -deltafV_pa_neg, p);
						partMdts[m].setValue(EMDL_CTF_DEFOCUS_ANGLE, phiDeg, p);
					}
				}
			}
		}
	}
}

Equation2x2::Equation2x2()
	: Axx(0.0), Axy(0.0), Ayy(0.0),
	  bx(0.0), by(0.0)
{

}

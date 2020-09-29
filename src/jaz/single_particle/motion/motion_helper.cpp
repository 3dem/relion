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

#include "motion_helper.h"

#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/interpolation.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/single_particle/Fourier_helper.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/single_particle/t_complex.h>
#include <src/jaz/single_particle/new_ft.h>
#include <omp.h>

using namespace gravis;

std::vector<std::vector<Image<RFLOAT>>> MotionHelper::movieCC(
		const std::vector<std::vector<Image<Complex>>>& movie,
		const std::vector<Image<Complex>>& preds,
		const std::vector<Image<RFLOAT> > &damageWeights,
		double pad, int threads)
{
	const int pc = movie.size();
	const int fc = movie[0].size();

	const int s = movie[0][0]().ydim;
	const int sh = s/2 + 1;

	std::vector<std::vector<Image<RFLOAT>>> out(pc);

	std::vector<Image<dComplex>> ccsFs(threads);
	std::vector<Image<double>> ccsRs(threads);
	
	const int s2 = (int)(pad * s);
	const int sh2 = s2/2 + 1;

	for (int t = 0; t < threads; t++)
	{
		ccsFs[t] = Image<dComplex>(sh2,s2);
		ccsFs[t].data.xinit = 0;
		ccsFs[t].data.yinit = 0;

		ccsRs[t] = Image<double>(s2,s2);
		ccsRs[t].data.xinit = 0;
		ccsRs[t].data.yinit = 0;
	}

	NewFFT::DoublePlan plan(s2,s2,1);

	for (int p = 0; p < pc; p++)
	{
		out[p] = std::vector<Image<RFLOAT>>(fc, Image<RFLOAT>(s2,s2));

#pragma omp parallel for num_threads(threads)
		for (int f = 0; f < fc; f++)
		{
			int t = omp_get_thread_num();

			for (int y = 0; y < s; y++)
				for (int x = 0; x < sh; x++)
				{
					Complex z = movie[p][f](y,x) * damageWeights[f](y,x) * preds[p](y,x).conj();

					const int yy = y < sh? y : s2 - (s - y);

					ccsFs[t](yy,x) = dComplex(z.real, z.imag);
				}

			NewFFT::inverseFourierTransform(ccsFs[t](), ccsRs[t](), plan);

			for (int y = 0; y < s2; y++)
				for (int x = 0; x < s2; x++)
				{
					out[p][f](y,x) = s * s * ccsRs[t](y,x);
				}
		}
	}

	return out;
}

/*std::vector<std::vector<Image<RFLOAT>>> MotionHelper::movieCC(
	Projector& projector0,
	Projector& projector1,
	const ObservationModel &obsModel,
	const MetaDataTable &viewParams,
	const std::vector<std::vector<Image<Complex> > > &movie,
	const std::vector<double> &sigma2,
	const std::vector<Image<RFLOAT> > &damageWeights,
	std::vector<ParFourierTransformer>& fts, int threads)
{
	const int pc = movie.size();
	const int fc = movie[0].size();

	const int s = movie[0][0]().ydim;
	const int sh = s/2 + 1;

	std::vector<std::vector<Image<RFLOAT>>> out(pc);

	std::vector<Image<Complex>> ccsFs(threads);
	std::vector<Image<RFLOAT>> ccsRs(threads);

	for (int t = 0; t < threads; t++)
	{
		ccsFs[t] = Image<Complex>(sh,s);
		ccsFs[t].data.xinit = 0;
		ccsFs[t].data.yinit = 0;

		ccsRs[t] = Image<RFLOAT>(s,s);
		ccsRs[t].data.xinit = 0;
		ccsRs[t].data.yinit = 0;
	}

	Image<Complex> pred;

	for (int p = 0; p < pc; p++)
	{
		out[p] = std::vector<Image<RFLOAT>>(fc);

		for (int f = 0; f < fc; f++)
		{
			out[p][f] = Image<RFLOAT>(s,s);
		}

		int randSubset;
		viewParams.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
		randSubset -= 1;

		if (randSubset == 0)
		{
			pred = obsModel.predictObservation(projector0, viewParams, p, true, true);
		}
		else
		{
			pred = obsModel.predictObservation(projector1, viewParams, p, true, true);
		}

		noiseNormalize(pred, sigma2, pred);

		#pragma omp parallel for num_threads(threads)
		for (int f = 0; f < fc; f++)
		{
			int t = omp_get_thread_num();

			for (int y = 0; y < s; y++)
			for (int x = 0; x < sh; x++)
			{
				ccsFs[t](y,x) = movie[p][f](y,x) * damageWeights[f](y,x) * pred(y,x).conj();
			}

			fts[t].inverseFourierTransform(ccsFs[t](), ccsRs[t]());

			for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)
			{
				out[p][f](y,x) = s * s * ccsRs[t](y,x);
			}
		}
	}

	return out;
}*/

/*std::vector<d2Vector> MotionHelper::getGlobalTrack(
	const std::vector<std::vector<Image<RFLOAT>>>& movieCC)
{
	const int pc = movieCC.size();
	const int fc = movieCC[0].size();

	const int s = movieCC[0][0]().xdim;
	const int sh = s/2 + 1;

	std::vector<d2Vector> out(fc);
	const double eps = 1e-30;

	std::vector<Image<RFLOAT>> e_sum(fc);

	for (int f = 0; f < fc; f++)
	{
		e_sum[f] = Image<RFLOAT>(s, s);
		e_sum[f].data.initZeros();

		for (int p = 0; p < pc; p++)
		{
			for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)
			{
				e_sum[f](y,x) += movieCC[p][f](y,x);
			}
		}

		d2Vector pos = Interpolation::quadraticMaxWrapXY(e_sum[f], eps);

		if (pos.x >= sh) pos.x -= s;
		if (pos.y >= sh) pos.y -= s;

		out[f] = pos;
	}

	return out;
}*/

std::vector<Image<RFLOAT> > MotionHelper::addCCs(
		const std::vector<std::vector<Image<RFLOAT>>> &movieCC)
{
	const int pc = movieCC.size();
	const int fc = movieCC[0].size();

	const int s = movieCC[0][0]().xdim;

	std::vector<Image<RFLOAT>> e_sum(fc);

	for (int f = 0; f < fc; f++)
	{
		e_sum[f] = Image<RFLOAT>(s, s);
		e_sum[f].data.initZeros();

		for (int p = 0; p < pc; p++)
		{
			for (int y = 0; y < s; y++)
				for (int x = 0; x < s; x++)
				{
					e_sum[f](y,x) += movieCC[p][f](y,x);
				}
		}
	}

	return e_sum;
}

std::vector<d2Vector> MotionHelper::getGlobalTrack(
		const std::vector<Image<RFLOAT>> &movieCcSum, double cc_pad)
{
	const int fc = movieCcSum.size();

	const int s = movieCcSum[0]().xdim;
	const int sh = s/2 + 1;

	std::vector<d2Vector> out(fc);
	const double eps = 1e-30;

	for (int f = 0; f < fc; f++)
	{
		d2Vector pos = Interpolation::quadraticMaxWrapXY(movieCcSum[f], eps);

		if (pos.x >= sh) pos.x -= s;
		if (pos.y >= sh) pos.y -= s;

		out[f] = pos/cc_pad;
	}

	return out;
}

std::vector<d2Vector> MotionHelper::getGlobalOffsets(
		const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
		const std::vector<std::vector<gravis::d2Vector>>& initialTracks,
		double cc_pad, double sigma, int wMax, int hMax, int threads)
{
	const int pc = movieCC.size();
	const int fc = movieCC[0].size();
	const int s = movieCC[0][0]().xdim;
	const int sh = s/2 + 1;
	const double eps = 1e-30;

	std::vector<d2Vector> out(pc);
	Image<RFLOAT> weight(s,s);

	for (int y = 0; y < s; y++)
		for (int x = 0; x < s; x++)
		{
			double xx = x >= sh? x - s: x;
			double yy = y >= sh? y - s: y;

			weight(y,x) = exp(-0.5*(xx*xx + yy*yy)/(sigma*sigma));
		}

#pragma omp parallel for num_threads(threads)
	for (int p = 0; p < pc; p++)
	{
		Image<RFLOAT> pSum(s,s);
		pSum.data.initZeros();

		for (int f = 0; f < fc; f++)
		{
			const d2Vector g = initialTracks[p][f];

			for (int y = 0; y < s; y++)
				for (int x = 0; x < s; x++)
				{
					pSum(y,x) += Interpolation::cubicXY(movieCC[p][f], x + g.x, y + g.y, 0, 0, true);
				}
		}

		for (int y = 0; y < s; y++)
			for (int x = 0; x < s; x++)
			{
				pSum(y,x) *= weight(y,x);
			}

		d2Vector out_p = Interpolation::quadraticMaxWrapXY(pSum, eps, wMax, hMax);
		
		if (out_p.x >= sh) out_p.x -= s;
		if (out_p.y >= sh) out_p.y -= s;

		out[p] = out_p/cc_pad;
	}

	return out;
}

void MotionHelper::noiseNormalize(
		const Image<Complex> &img, const std::vector<double> &sigma2, Image<Complex>& dest)
{
	int wf = img().xdim;
	int w = 2*wf - 1;
	int h = img().ydim;

	const double area = 0.25*PI*w*h;

	if (dest.data.xdim != img.data.xdim || dest.data.ydim != img.data.ydim)
	{
		dest.data.reshape(img.data);
	}

	dest.data.xinit = 0;
	dest.data.yinit = 0;

	for (int y = 0; y < h; y++)
		for (int x = 0; x < wf; x++)
		{
			if (x == 0 && y == 0)
			{
				dest(y,x) = Complex(0.0);
				continue;
			}

			const double yy = y < wf? y : y - h;
			const double xx = x;

			const int r = (int) sqrt(xx*xx + yy*yy);

			if (r >= wf || sigma2[r] == 0.0)
			{
				dest(y,x) = Complex(0.0);
			}
			else
			{
				dest(y,x) = DIRECT_A2D_ELEM(img.data, y, x) / sqrt(sigma2[r]*area);
			}
		}
}

void MotionHelper::writeTracks(
		const std::vector<std::vector<d2Vector>>& tracksInPix,
		std::string fn, double angpix)
{
	const int pc = tracksInPix.size();
	const int fc = tracksInPix[0].size();

	std::string path = fn.substr(0, fn.find_last_of('/'));
	mktree(path);

	std::ofstream ofs(fn);
	MetaDataTable mdt;

	mdt.setName("general");
	mdt.setIsList(true);
	mdt.addObject();
	mdt.setValue(EMDL_PARTICLE_NUMBER, pc);

	mdt.write(ofs);
	mdt.clear();

	for (int p = 0; p < pc; p++)
	{
		std::stringstream sts;
		sts << p;
		mdt.setName(sts.str());

		for (int f = 0; f < fc; f++)
		{
			mdt.addObject();
			
			mdt.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, angpix * tracksInPix[p][f].x);
			mdt.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, angpix * tracksInPix[p][f].y);
		}

		mdt.write(ofs);
		mdt.clear();
	}
}

std::vector<std::vector<d2Vector>> MotionHelper::readTracksInPix(std::string fn, double angpix)
{
	std::ifstream ifs(fn);

	if (!ifs)
	{
		REPORT_ERROR("MotionHelper::readTracks: unable to read " + fn + ".");
	}

	MetaDataTable mdt;

	mdt.readStar(ifs, "general");

	int pc;

	if (!mdt.getValue(EMDL_PARTICLE_NUMBER, pc))
	{
		REPORT_ERROR("MotionHelper::readTracks: missing particle number in "+fn+".");
	}

	std::vector<std::vector<d2Vector>> out(pc);
	int fc = 0, lastFc = 0;

	for (int p = 0; p < pc; p++)
	{
		std::stringstream sts;
		sts << p;
		mdt.readStar(ifs, sts.str());

		fc = mdt.numberOfObjects();

		if (p > 0 && fc != lastFc)
		{
			REPORT_ERROR("MotionHelper::readTracks: broken file: "+fn+".");
		}

		lastFc = fc;

		out[p] = std::vector<d2Vector>(fc);

		for (int f = 0; f < fc; f++)
		{
			mdt.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, out[p][f].x, f);
			mdt.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, out[p][f].y, f);
			
			out[p][f] /= angpix;
		}
	}

	return out;
}



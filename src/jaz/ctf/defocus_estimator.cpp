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

#include "defocus_estimator.h"
#include "defocus_helper.h"
#include "ctf_refiner.h"

#include <src/ctf.h>
#include <src/time.h>

#include <src/jaz/obs_model.h>
#include <src/jaz/reference_map.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/gravis/t2Vector.h>

using namespace gravis;

DefocusEstimator::DefocusEstimator()
:	ready(false)
{}

void DefocusEstimator::read(IOParser &parser, int argc, char *argv[])
{
	noGlobAstig = !parser.checkOption("--glob_astig", "Estimate per-micrograph astigmatism");
	fitCs = parser.checkOption("--fit_Cs", "Fit spherical aberration (per micrograph)");
	fitPhase = parser.checkOption("--fit_phase", "Fit phase shift (amplitude contrast) per-micrograph");
	globOnly = parser.checkOption("--glob", "Only perform per-micrograph fit");

	fitAstigmatism = parser.checkOption("--astig", "Estimate independent astigmatism for each particle");

	kmin = textToFloat(parser.getOption("--kmin_defocus",
						"Inner freq. threshold for defocus estimation [Angst]", "30.0"));

	defocusRange = textToFloat(parser.getOption("--range", "Defocus scan range (in A)", "2000."));

	/* #TODO: crash on:
		!glob_astig && (fit_Cs || fit_phase || glob)
	*/
}

void DefocusEstimator::init(
		int verb, int nr_omp_threads,
		bool debug, bool diag, std::string outPath,
		ReferenceMap *reference, ObservationModel *obsModel)
{
	this->verb = verb;
	this->nr_omp_threads = nr_omp_threads;

	this->debug = debug;
	this->diag = diag;
	this->outPath = outPath;

	this->reference = reference;
	this->obsModel = obsModel;

	angpix = obsModel->getPixelSizes();
	obsModel->getBoxSizes(s, sh);

	freqWeights.resize(angpix.size());

	for (int i = 0; i < angpix.size(); i++)
	{
		freqWeights[i] = reference->getHollowWeight(kmin, s[i], angpix[i]);
	}

	ready = true;
}

void DefocusEstimator::processMicrograph(
		long g, MetaDataTable& mdt,
		const std::vector<Image<Complex>>& obs,
		const std::vector<Image<Complex>>& pred)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: DefocusEstimator::processMicrograph: DefocusEstimator not initialized.");
	}

	long pc = obs.size();

	std::stringstream stsg;
	stsg << g;

	if (!noGlobAstig)
	{
		REPORT_ERROR("Per-micrograph CTF-refinement temporarily disabled.");

		/*CTF ctf0;
		ctf0.readByGroup(mdt, obsModel, 0);

		Image<RFLOAT> dataVis;

		if (diag)
		{
			Image<RFLOAT> dotp0(sh,s), cc(sh,s), wgh0(sh,s), wgh1(sh,s), dotp0_full(s,s);
			dotp0.data.initZeros();
			wgh0.data.initZeros();
			wgh1.data.initZeros();

			for (long p = 0; p < pc; p++)
			{
				for (long y = 0; y < s; y++)
				for (long x = 0; x < sh; x++)
				{
					Complex vx = DIRECT_A2D_ELEM(pred[p].data, y, x);
					const Complex vy = DIRECT_A2D_ELEM(obs[p].data, y, x);

					dotp0(y,x) += vy.real*vx.real + vy.imag*vx.imag;
					wgh0(y,x) += vx.norm();
					wgh1(y,x) += vy.norm();
				}
			}

			for (long y = 0; y < s; y++)
			for (long x = 0; x < sh; x++)
			{
				double nrm = sqrt(wgh0(y,x) * wgh1(y,x));
				cc(y,x) = nrm > 0.0? 10.0 * dotp0(y,x) / nrm : 0.0;
			}

			FftwHelper::decenterDouble2D(cc.data, dotp0_full.data);

			ImageLog::write(dotp0_full, outPath+"diag_m"+stsg.str()+"_global_CTF-data");
			dataVis = FilterHelper::polarBlur(dotp0_full, 10.0);
			ImageLog::write(dataVis, outPath+"diag_m"+stsg.str()+"_global_CTF-data-blurred");

			Image<RFLOAT> ctfFit(s,s);
			ctf0.getCenteredImage(ctfFit.data, angpix, false, false, false, false);
			ImageLog::write(ctfFit, outPath+"diag_m"+stsg.str()+"_global_CTF-fit_initial");

			ctfFit.data.xinit = 0;
			ctfFit.data.yinit = 0;

			Image<RFLOAT> vis = FilterHelper::sectorBlend(dataVis, ctfFit, 12);
			ImageLog::write(vis, outPath+"diag_m"+stsg.str()+"_global_CTF-vis_initial");
		}

		double u, v, phi, phase, newCs;

		mdt.getValue(EMDL_CTF_PHASESHIFT, phase, 0);
		mdt.getValue(EMDL_CTF_CS, newCs, 0);

		if (fitCs)
		{
			if (verb > 0)
			{
				std::cout << "initial phi and Cs: " << phase << ", " << newCs << "\n";
			}

			DefocusHelper::findAstigmatismPhaseAndCsNM(
				pred, obs, freqWeight, ctf0, angpix, &u, &v, &phi, &phase, &newCs);

			if (verb > 0)
			{
				std::cout << "final phi and Cs: " << phase << ", " << newCs << "\n";
			}
		}
		else if (fitPhase)
		{
			if (verb > 0)
			{
				std::cout << "initial phase shift: " << phase << "\n";
			}

			DefocusHelper::findAstigmatismAndPhaseNM(
				pred, obs, freqWeight, ctf0, angpix, &u, &v, &phi, &phase);

			if (verb > 0)
			{
				std::cout << "final phase shift: " << phase << "\n";
			}
		}
		else
		{
			DefocusHelper::findAstigmatismNM(
				pred, obs, freqWeight, ctf0, angpix, &u, &v, &phi);
		}

		for (long p = 0; p < pc; p++)
		{
			mdt.setValue(EMDL_CTF_DEFOCUSU, u, p);
			mdt.setValue(EMDL_CTF_DEFOCUSV, v, p);
			mdt.setValue(EMDL_CTF_DEFOCUS_ANGLE, phi, p);

			if (fitPhase) mdt.setValue(EMDL_CTF_PHASESHIFT, phase, p);
			if (fitCs) mdt.setValue(EMDL_CTF_CS, newCs, p);
		}

		if (diag)
		{
			CTF ctf1;
			ctf1.readByGroup(mdt, obsModel, 0);

			Image<RFLOAT> ctfFit(s,s);
			ctf1.getCenteredImage(ctfFit.data, angpix, false, false, false, false);
			ImageLog::write(ctfFit, outPath+"diag_m"+stsg.str()+"_global_CTF-fit_final");
			ctfFit.data.xinit = 0;
			ctfFit.data.yinit = 0;

			Image<RFLOAT> vis = FilterHelper::sectorBlend(dataVis, ctfFit, 12);
			ImageLog::write(vis, outPath+"diag_m"+stsg.str()+"_global_CTF-vis_final");
		}*/
	}

	if (globOnly) return;

	if (diag)
	{
		std::ofstream ofst(outPath+"diag_m"+stsg.str()+"_defocus_cost.dat");
		std::ofstream ofsto(outPath+"diag_m"+stsg.str()+"_defocus_opt.dat");

		for (long p = 0; p < pc; p++)
		{
			const int og = obsModel->getOpticsGroup(mdt, p);

			if (obsModel->getCtfPremultiplied(og))
				REPORT_ERROR("ERROR: you cannot perform defocus estimation on CTF-premultiplied images...");

			CTF ctf0;
			ctf0.readByGroup(mdt, obsModel, p);

			std::vector<d2Vector> cost = DefocusHelper::diagnoseDefocus(
				pred[p], obs[p], freqWeights[og],
				ctf0, angpix[og], defocusRange, 100, nr_omp_threads);

			double cMin = cost[0][1];
			double dOpt = cost[0][0];

			for (int i = 0; i < cost.size(); i++)
			{
				ofst << cost[i][0] << " " << cost[i][1] << "\n";

				if (cost[i][1] < cMin)
				{
					cMin = cost[i][1];
					dOpt = cost[i][0];
				}
			}

			ofsto << dOpt << " " << cMin << "\n";

			ofst << "\n";
		}
	}

	// Parallel loop over all particles in this micrograph
	#pragma omp parallel for num_threads(nr_omp_threads)
	for (long p = 0; p < pc; p++)
	{
		const int og = obsModel->getOpticsGroup(mdt, p);

		std::stringstream stsp;
		stsp << p;

		CTF ctf0;
		ctf0.readByGroup(mdt, obsModel, p);

		if (fitAstigmatism)
		{
			double u, v, phi;
			DefocusHelper::findAstigmatismNM(
				pred[p], obs[p], freqWeights[og], ctf0,
				angpix[og], &u, &v, &phi);

			mdt.setValue(EMDL_CTF_DEFOCUSU, u, p);
			mdt.setValue(EMDL_CTF_DEFOCUSV, v, p);
			mdt.setValue(EMDL_CTF_DEFOCUS_ANGLE, phi, p);
		}
		else
		{
			double u, v;
			DefocusHelper::findDefocus1D(
				pred[p], obs[p], freqWeights[og], ctf0,
				angpix[og], &u, &v, defocusRange);

			mdt.setValue(EMDL_CTF_DEFOCUSU, u, p);
			mdt.setValue(EMDL_CTF_DEFOCUSV, v, p);
		}

	}

	// Output a diagnostic Postscript file
	writeEPS(mdt);

	// Now write out STAR file with optimised values for this micrograph
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	mdt.write(outRoot + "_defocus_fit.star");
}

void DefocusEstimator::writeEPS(const MetaDataTable& mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: DefocusEstimator::writeEPS: DefocusEstimator not initialized.");
	}

	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);

	FileName fn_eps = outRoot + "_ctf-refine_fit.eps";

	CPlot2D plot2D(fn_eps);
	plot2D.SetXAxisSize(600);
	plot2D.SetYAxisSize(600);
	plot2D.SetDrawLegend(false);
	plot2D.SetFlipY(true);

	RFLOAT min_defocus = 99.e10;
	RFLOAT max_defocus = -99.e10;

	const int pc = mdt.numberOfObjects();

	for (int p = 0; p < pc; p++)
	{
		RFLOAT defU, defV;
		mdt.getValue(EMDL_CTF_DEFOCUSU, defU, p);
		mdt.getValue(EMDL_CTF_DEFOCUSV, defV, p);
		defU = (defU + defV) / 2.;

		min_defocus = XMIPP_MIN(min_defocus, defU);
		max_defocus = XMIPP_MAX(max_defocus, defU);
	}

	for (int p = 0; p < pc; p++)
	{
		RFLOAT defU, defV;
		RFLOAT xcoor, ycoor;

		mdt.getValue(EMDL_IMAGE_COORD_X, xcoor, p);
		mdt.getValue(EMDL_IMAGE_COORD_Y, ycoor, p);
		mdt.getValue(EMDL_CTF_DEFOCUSU, defU, p);
		mdt.getValue(EMDL_CTF_DEFOCUSV, defV, p);

		defU = (defU + defV) / 2.;

		RFLOAT val  = (defU - min_defocus) / (max_defocus - min_defocus);
		const RFLOAT eps = 1e-10;
		if (max_defocus - min_defocus < eps) val = 0.5; // to avoid NaN in color

		CDataSet dataSet;

		dataSet.SetDrawMarker(true);
		dataSet.SetDrawLine(false);
		dataSet.SetMarkerSize(10);
		dataSet.SetDatasetColor(val, 0.8 * val, 1.0 - val);

		CDataPoint point(xcoor, ycoor);

		dataSet.AddDataPoint(point);

		plot2D.AddDataSet(dataSet);
	}

	char title[256];
	snprintf(title, 255, "Defocus range from blue to orange: %.0f A", max_defocus - min_defocus);
	plot2D.SetXAxisTitle(title);

	plot2D.OutputPostScriptPlot(fn_eps);
}

bool DefocusEstimator::isFinished(const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: DefocusEstimator::isFinished: DefocusEstimator not initialized.");
	}

	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);

	return exists(outRoot + "_defocus_fit.star");
}

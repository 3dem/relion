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

#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/reference_map.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/single_particle/ctf/modular_ctf_optimisation.h>
#include <src/jaz/optimization/lbfgs.h>

using namespace gravis;

DefocusEstimator::DefocusEstimator()
:	ready(false)
{}

void DefocusEstimator::read(IOParser &parser, int argc, char *argv[])
{
	fittingMode = parser.getOption("--fit_mode",          
                                     "String of 5 characters describing whether to fit the phase shift (1), \n\
                                      defocus (2), astigmatism (3), spherical aberration (4) and B-factors (5) \n\
                                      per particle ('p'), per micrograph ('m') or to keep them fixed ('f')\n\
                                      during the per-micrograph CTF refinement.", "fpmfm");
								   
	max_iters = textToInteger(parser.getOption("--max_defocus_iters", 
						"Maximum number of iterations for CTF refinement.", "100"));
			
	bruteForcePre = parser.checkOption("--bf0", 
                                     "Perform brute-force per-particle defocus search (as in RELION 3.0) prior \n\
                                      to the per-micrograph CTF refinement.");

    bruteForcePost = parser.checkOption("--bf1", 
						"Perform brute-force defocus search after CTF refinement.");
			
	bruteForceOnly = parser.checkOption("--bf_only", 
						"Skip CTF refinement and only perform a brute-force defocus search.");
									  
	defocusRange = textToDouble(parser.getOption("--bf_range", 
						"Defocus scan range (in A) for brute-force search.", "2000."));

	fitAstigmatism = parser.checkOption("--legacy_astig", 
						"Estimate independent per-particle astigmatism (from RELION 3.0)");

	kmin = textToFloat(parser.getOption("--kmin_defocus",
						"Inner freq. threshold for defocus estimation [Angst]", "30.0"));
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
	
	if (verb > 0)
	{
		std::vector<std::string> names{"per particle (p)", "per micrograph (m)", "fixed (f)" };
		
		if (!ModularCtfOptimisation::validateModeString(fittingMode))
		{
			REPORT_ERROR_STR("DefocusEstimator::init: illegal fitting mode string: " << fittingMode);
		}
		
		std::vector<ModularCtfOptimisation::Mode> modes = ModularCtfOptimisation::decodeModes(fittingMode);
		
		std::cout << " + Defocus fitting-mode string: " << fittingMode << "\n";
		std::cout << "     => Estimating:" << "\n";
		std::cout << "            phase-shift:     " << names[(int)modes[ModularCtfOptimisation::Phase]] << "\n";
		std::cout << "            defocus:         " << names[(int)modes[ModularCtfOptimisation::Defocus]] << "\n";
		std::cout << "            astigmatism:     " << names[(int)modes[ModularCtfOptimisation::Astigmatism1]] << "\n";
		std::cout << "            sph. aberration: " << names[(int)modes[ModularCtfOptimisation::SphericalAberration]] << "\n";
		std::cout << "            B/Scale factors: " << names[(int)modes[ModularCtfOptimisation::BFactor]] << std::endl;
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

	if (bruteForcePre || bruteForceOnly) 
	{
		bruteForceFit(g, mdt, obs, pred, "pre");
	}
		
	if (!bruteForceOnly)
	{
		ModularCtfOptimisation mco(mdt, obsModel, obs, pred, freqWeights, fittingMode, nr_omp_threads);
		std::vector<double> x0 = mco.encodeInitial();
			
		std::vector<double> x = LBFGS::optimize(x0, mco, debug, max_iters, 1e-9);
		
		mco.writeToTable(x);
			
		if (bruteForcePost) 
		{
			bruteForceFit(g, mdt, obs, pred, "post");
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

void DefocusEstimator::bruteForceFit(
		long g, MetaDataTable &mdt, 
		const std::vector<Image<Complex> > &obs, 
		const std::vector<Image<Complex> > &pred,
		std::string tag)
{
	long pc = obs.size();
	
	std::stringstream stsg;
	stsg << g;
	
	if (diag)
	{
		std::ofstream ofst(outPath+"diag_m"+stsg.str()+"_bf-defocus-"+tag+"_cost.dat");
		std::ofstream ofsto(outPath+"diag_m"+stsg.str()+"_bf-defocus-"+tag+"_opt.dat");

		for (long p = 0; p < pc; p++)
		{
			const int og = obsModel->getOpticsGroup(mdt, p);
			
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

			/*if (debug)
			{
				double u0, v0;
				mdt.getValue(EMDL_CTF_DEFOCUSU, u0, p);
				mdt.getValue(EMDL_CTF_DEFOCUSV, v0, p);
				
				std::cout << u0 << " -> " << u << ", " << v0 << " -> " << v << "\n";
			}*/
			
			mdt.setValue(EMDL_CTF_DEFOCUSU, u, p);
			mdt.setValue(EMDL_CTF_DEFOCUSV, v, p);
		}

	}	
}

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

#include "tilt_estimator.h"
#include "tilt_helper.h"
#include "ctf_refiner.h"

#include <src/jaz/obs_model.h>
#include <src/jaz/reference_map.h>
#include <src/jaz/image_op.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/gravis/t2Vector.h>

#include <src/args.h>
#include <src/ctf.h>

#include <set>

#include <omp.h>

using namespace gravis;


TiltEstimator::TiltEstimator()
	:	ready(false)
{}

void TiltEstimator::read(IOParser &parser, int argc, char *argv[])
{
	kmin = textToFloat(parser.getOption("--kmin_tilt", 
						"Inner freq. threshold for beamtilt estimation [Angst]", "20.0"));
		
	aniso = parser.checkOption("--anisotropic_tilt", "Use anisotropic coma model");
}

void TiltEstimator::init(
		int verb, int s, int nr_omp_threads, 
		bool debug, bool diag, std::string outPath, 
		MetaDataTable& mdt0,
		ReferenceMap* reference, ObservationModel* obsModel)
{
	this->verb = verb;
	this->s = s;
	sh = s/2 + 1;
	this->nr_omp_threads = nr_omp_threads;
	
	this->debug = debug;
	this->diag = diag;
	this->outPath = outPath;
	
	this->reference = reference;
	this->obsModel = obsModel;
	
	angpix = obsModel->angpix;
	
	tiltClasses.clear();
	
	if (mdt0.containsLabel(EMDL_PARTICLE_BEAM_TILT_CLASS))
	{
		std::set<int> tiltClassSet;
	
		const long pc = mdt0.numberOfObjects();
		
		for (long p = 0; p < pc; p++)
		{
			int s = 0;
			mdt0.getValue(EMDL_PARTICLE_BEAM_TILT_CLASS, s, p);
			
			tiltClassSet.insert(s);			
		}
		
		if (verb > 0)
		{
			std::cout << " + " << tiltClassSet.size() << " beam tilt classes found." << std::endl;
		}
		
		for (std::set<int>::iterator it = tiltClassSet.begin(); it != tiltClassSet.end(); it++)
		{
			tiltClasses.push_back(*it);
			classNameToIndex[*it] = tiltClasses.size()-1;
			
			std::cout << tiltClasses[tiltClasses.size()-1] << "\n";
		}
	}
	else
	{
		tiltClasses.push_back(0);
		classNameToIndex[0] = 0;
	}
	
	ready = true;
}

void TiltEstimator::processMicrograph(
		long g, MetaDataTable& mdt, 
		const std::vector<Image<Complex>>& obs,
		const std::vector<Image<Complex>>& pred)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::processMicrograph: TiltEstimator not initialized.");
	}
	
	const int pc = mdt.numberOfObjects();
	const int cc = tiltClasses.size();
		
	std::vector<Image<Complex>> xyAcc(nr_omp_threads*cc);
	std::vector<Image<RFLOAT>> wAcc(nr_omp_threads*cc);
	
	for (int i = 0; i < nr_omp_threads*cc; i++)
	{
		xyAcc[i] = Image<Complex>(sh,s);
		xyAcc[i].data.initZeros();
		
		wAcc[i] = Image<RFLOAT>(sh,s);
		wAcc[i].data.initZeros();
	}
	
	std::vector<int> partsPerClass(cc, 0);
	
	#pragma omp parallel for num_threads(nr_omp_threads)
	for (long p = 0; p < pc; p++)
	{
		CTF ctf;
		ctf.read(mdt, mdt, p);
		
		int threadnum = omp_get_thread_num();
		
		int cn = 0;
		
		if (cc > 1)
		{
			mdt.getValue(EMDL_PARTICLE_BEAM_TILT_CLASS, cn, p);
		}
		
		const int ci = classNameToIndex[cn];
		
		partsPerClass[ci]++;
		
		TiltHelper::updateTiltShift(
			pred[p], obs[p], ctf, obsModel->angpix, 
			xyAcc[cc*threadnum + ci], 
			wAcc[cc*threadnum + ci]);
	}
	
	// Combine the accumulated weights from all threads for this subset, 
	// store weighted sums in xyAccSum and wAccSum
	
	for (int ci = 0; ci < cc; ci++)
	{
		if (partsPerClass[ci] == 0) continue;
		
		Image<Complex> xyAccSum(sh,s);
		Image<RFLOAT> wAccSum(sh,s);
		
		for (int threadnum = 0; threadnum < nr_omp_threads; threadnum++)
		{
			ImageOp::linearCombination(xyAccSum, xyAcc[cc*threadnum + ci], 1.0, 1.0, xyAccSum);
			ImageOp::linearCombination(wAccSum, wAcc[cc*threadnum + ci], 1.0, 1.0, wAccSum);
		}
		
		// Write out the intermediate results per-micrograph:
		
		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
		
		std::stringstream sts;
		sts << tiltClasses[ci];
		
		ComplexIO::write(xyAccSum(), outRoot + "_xyAcc_class_" + sts.str(), ".mrc");
		wAccSum.write(outRoot+"_wAcc_class_" + sts.str() + ".mrc");
	}
}

void TiltEstimator::parametricFit(
		const std::vector<MetaDataTable>& mdts, 
		double Cs, double lambda, 
		MetaDataTable& mdtOut)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::parametricFit: TiltEstimator not initialized.");
	}
	
	if (verb > 0)
	{
		std::cout << " + Fitting beam tilt ..." << std::endl;
	}
	
	const int gc = mdts.size();	
	const int cc = tiltClasses.size();
	
	std::vector<d2Vector> tiltPerClass(cc);
	
	for (int ci = 0; ci < cc; ci++)
	{	
		std::stringstream sts;
		sts << tiltClasses[ci];
		std::string cns = sts.str();
		
		Image<Complex> xyAccSum(sh,s);
		Image<RFLOAT> wAccSum(sh,s);
		
		xyAccSum.data.initZeros();
		wAccSum.data.initZeros();
			
		for (long g = 0; g < gc; g++)
		{
			std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);
			
			if (   exists(outRoot+"_xyAcc_class_"+cns+"_real.mrc")
				&& exists(outRoot+"_xyAcc_class_"+cns+"_imag.mrc")
				&& exists(outRoot+ "_wAcc_class_"+cns+".mrc"))
			{
				Image<Complex> xyAcc;
				Image<RFLOAT> wAcc;
				
				wAcc.read(outRoot+"_wAcc_class_"+cns+".mrc");
				ComplexIO::read(xyAcc, outRoot+"_xyAcc_class_"+cns, ".mrc");
				
				xyAccSum() += xyAcc();
				wAccSum()  +=  wAcc();
			}	
		}
			
		Image<RFLOAT> wgh, phase, fit, phaseFull, fitFull;
		
		FilterHelper::getPhase(xyAccSum, phase);
		
		Image<Complex> xyNrm(sh,s);
	
		double kmin_px = obsModel->angToPix(kmin, s);
				 
		Image<RFLOAT> wgh0 = reference->getHollowWeight(kmin_px);
	
		FilterHelper::multiply(wAccSum, wgh0, wgh);
		
		for (int y = 0; y < s; y++)
		for (int x = 0; x < sh; x++)
		{
			xyNrm(y,x) = wAccSum(y,x) > 0.0? xyAccSum(y,x)/wAccSum(y,x) : Complex(0.0, 0.0);
		}
		
		Image<RFLOAT> wghFull;
		FftwHelper::decenterDouble2D(wgh(), wghFull());
		
		if (debug)
		{
			ImageLog::write(wghFull, outPath + "beamtilt_weight-full_class_"+cns);
		}
		
		FftwHelper::decenterUnflip2D(phase.data, phaseFull.data);
			
		ImageLog::write(phaseFull, outPath + "beamtilt_delta-phase_per-pixel_class_"+cns);
		
		double shift_x, shift_y, tilt_x, tilt_y;
		
		TiltHelper::fitTiltShift(
			phase, wgh, Cs, lambda, angpix,
			&shift_x, &shift_y, &tilt_x, &tilt_y, &fit);
			
		FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
		
		ImageLog::write(fitFull, outPath + "beamtilt_delta-phase_lin-fit_class_"+cns);
		
		std::ofstream os(outPath+"beamtilt_lin-fit_class_"+cns+".txt");
		os << "beamtilt_x = " << tilt_x << "\n";
		os << "beamtilt_y = " << tilt_y << "\n";
		os.close();
		
		double tilt_xx, tilt_xy, tilt_yy;
		
		if (aniso)
		{
			TiltHelper::optimizeAnisoTilt(
				xyNrm, wgh, Cs, lambda, angpix, false,
				shift_x, shift_y, tilt_x, tilt_y,
				&shift_x, &shift_y, &tilt_x, &tilt_y,
				&tilt_xx, &tilt_xy, &tilt_yy, &fit);
		}
		else
		{
			TiltHelper::optimizeTilt(
				xyNrm, wgh, Cs, lambda, angpix, false,
				shift_x, shift_y, tilt_x, tilt_y,
				&shift_x, &shift_y, &tilt_x, &tilt_y, &fit);
		}
		
		FftwHelper::decenterUnflip2D(fit.data, fitFull.data);
		
		ImageLog::write(fitFull, outPath+"beamtilt_delta-phase_iter-fit_class_"+cns);
		
		std::ofstream os2(outPath+"beamtilt_iter-fit_class_"+cns+".txt");
		os2 << "beamtilt_x = " << tilt_x << "\n";
		os2 << "beamtilt_y = " << tilt_y << "\n";
		os2.close();
		
		tiltPerClass[ci] = d2Vector(tilt_x, tilt_y);
	}
	
	// Now write the beamtilt into mdtOut	
	const long tpc = mdtOut.numberOfObjects();
	
	for (long p = 0; p < tpc; p++)
	{
		int cn = 0;
		
		if (cc > 1)
		{
			mdtOut.getValue(EMDL_PARTICLE_BEAM_TILT_CLASS, cn, p);
		}
		
		const d2Vector t = tiltPerClass[classNameToIndex[cn]];		
		
		mdtOut.setValue(EMDL_IMAGE_BEAMTILT_X, t.x, p);
		mdtOut.setValue(EMDL_IMAGE_BEAMTILT_Y, t.y, p);
	}
}

bool TiltEstimator::isFinished(const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::isFinished: DefocusEstimator not initialized.");
	}
	
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	
	bool allDone = true;
	
	const int cc = tiltClasses.size();
	
	for (int ci = 0; ci < cc; ci++)
	{	
		std::stringstream sts;
		sts << tiltClasses[ci];
		std::string cns = sts.str();
		
		if (   !exists(outRoot+"_xyAcc_class_"+cns+"_real.mrc")
			|| !exists(outRoot+"_xyAcc_class_"+cns+"_imag.mrc")
			|| !exists(outRoot+"_wAcc_class_"+cns+".mrc"))
		{
			allDone = false;
			break;
		}
	}
	
	return allDone;
}

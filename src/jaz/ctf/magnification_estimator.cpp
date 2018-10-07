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

#include "magnification_estimator.h"
#include "magnification_helper.h"
#include "ctf_refiner.h"
#include "equation2x2.h"

#include <src/jaz/reference_map.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/image_log.h>

#include <omp.h>


using namespace gravis;

MagnificationEstimator::MagnificationEstimator()
{
	
}

void MagnificationEstimator::read(
		IOParser &parser, int argc, char *argv[])
{
	kmin = textToFloat(parser.getOption("--kmin_mag", 
				"Inner freq. threshold for anisotropic magnification estimation [Angst]", "20.0"));
	
	adaptAstig = !parser.checkOption("--keep_astig", "Do not translate astigmatism into new coordinates");
	perMgAstig = !parser.checkOption("--part_astig", "Allow astigmatism to vary among the particles of a micrograph");
}

void MagnificationEstimator::init(
		int verb, int s, int nr_omp_threads, 
		bool debug, bool diag, 
		std::string outPath, 
		ReferenceMap* reference, 
		ObservationModel* obsModel)
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
	
	angpix = obsModel->getPixelSize(0);
	
	ready = true;
}

void MagnificationEstimator::processMicrograph(
		long g, MetaDataTable& mdt, 
		const std::vector<Image<Complex>>& obs, 
		const std::vector<Image<Complex>>& pred,
		const std::vector<Volume<t2Vector<Complex>>>& predGradient)
{
	if (!ready)
	{
		REPORT_ERROR_STR("ERROR: MagnificationEstimator::processMicrograph: "
						 << "MagnificationEstimator not initialized.");
	}
	
	const int pc = mdt.numberOfObjects();
	
	std::vector<int> optGroups = obsModel->getOptGroupsPresent_oneBased(mdt);	
	const int cc = optGroups.size();
	
	std::vector<int> groupToIndex(obsModel->numberOfOpticsGroups()+1, -1);
	
	for (int i = 0; i < cc; i++)
	{
		groupToIndex[optGroups[i]] = i;
	}
	
	std::vector<Volume<Equation2x2>> magEqs(nr_omp_threads*cc);
	
	for (int i = 0; i < nr_omp_threads*cc; i++)
	{
		magEqs[i] = Volume<Equation2x2>(sh,s,1);
	}
	
	#pragma omp parallel for num_threads(nr_omp_threads)
	for (long p = 0; p < pc; p++)
	{
		CTF ctf;
		ctf.readByGroup(mdt, obsModel, p);
		
		int threadnum = omp_get_thread_num();
		
		int og;
		mdt.getValue(EMDL_IMAGE_OPTICS_GROUP, og, p);

		const int ci = groupToIndex[og];
		
		MagnificationHelper::updateScaleFreq(
			pred[p], predGradient[p], obs[p], ctf, angpix, magEqs[cc*threadnum + ci]);
	}
	
	for (int ci = 0; ci < cc; ci++)
	{
		Volume<Equation2x2> magEq(sh,s,1);
				
		for (int threadnum = 0; threadnum < nr_omp_threads; threadnum++)
		{
			magEq += magEqs[cc*threadnum + ci];
		}
		
		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
		
		std::stringstream sts;
		sts << optGroups[ci];
		
		MagnificationHelper::writeEQs(magEq, outRoot+"_mag_optics-group_" + sts.str());
	}
}

void MagnificationEstimator::parametricFit(
		std::vector<MetaDataTable> &mdts, MetaDataTable &optOut)
{
	if (!ready)
	{
		REPORT_ERROR_STR("ERROR: MagnificationEstimator::parametricFit: "
					 << "MagnificationEstimator not initialized.");
	}
	
	if (verb > 0)
	{
		std::cout << " + Fitting anisotropic magnification ..." << std::endl;
	}
	
	const int gc = mdts.size();
	const int ogc = obsModel->numberOfOpticsGroups();
	
	std::vector<Matrix2D<RFLOAT>> mat_by_optGroup(ogc);
	
	for (int og = 0; og < ogc; og++)
	{
		Volume<Equation2x2> magEqs(sh,s,1), magEqsG(sh,s,1);
		
		std::stringstream sts;
		sts << (og+1);
		
		for (long g = 0; g < gc; g++)
		{
			std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);
			
			try
			{
				MagnificationHelper::readEQs(outRoot + "_mag_optics-group_" + sts.str(), magEqsG);
			
				magEqs += magEqsG;
			}
			catch (RelionError e)
			{}
		}
		
		Image<RFLOAT> flowx, flowy;
		MagnificationHelper::solvePerPixel(magEqs, flowx, flowy);
			
		Image<RFLOAT> flowxFull, flowyFull;
		FftwHelper::decenterUnflip2D(flowx.data, flowxFull.data);
		FftwHelper::decenterUnflip2D(flowy.data, flowyFull.data);
		
		ImageLog::write(flowxFull, outPath + "mag_disp_x_optics-group_" + sts.str());
		ImageLog::write(flowyFull, outPath + "mag_disp_y_optics-group_" + sts.str());
		
		Image<RFLOAT> freqWght = reference->getHollowWeight(kmin, s, angpix);
		
		Matrix2D<RFLOAT> mat = MagnificationHelper::solveLinearlyFreq(magEqs, freqWght, flowx, flowy);
		
		FftwHelper::decenterUnflip2D(flowx.data, flowxFull.data);
		FftwHelper::decenterUnflip2D(flowy.data, flowyFull.data);
		
		ImageLog::write(flowxFull, outPath + "mag_disp_x_fit_optics-group_" + sts.str());
		ImageLog::write(flowyFull, outPath + "mag_disp_y_fit_optics-group_" + sts.str());
		
		std::ofstream os(outPath + "mag_matrix_optics-group_" + sts.str() + ".txt");
		os << mat(0,0) << " " << mat(0,1) << "\n";
		os << mat(1,0) << " " << mat(1,1) << "\n";
		os.close();
		
		mat_by_optGroup[og] = mat;
		
		// @TODO: Add support for optics groups!
		
		Matrix2D<RFLOAT> mat0(2,2);
		mat0.initIdentity();
		
		optOut.getValue(EMDL_IMAGE_MAG_MATRIX_00, mat0(0,0), 0);
		optOut.getValue(EMDL_IMAGE_MAG_MATRIX_01, mat0(0,1), 0);
		optOut.getValue(EMDL_IMAGE_MAG_MATRIX_10, mat0(1,0), 0);
		optOut.getValue(EMDL_IMAGE_MAG_MATRIX_11, mat0(1,1), 0);
		
		Matrix2D<RFLOAT> mat1 = mat * mat0;
		
		optOut.setValue(EMDL_IMAGE_MAG_MATRIX_00, mat1(0,0), 0);
		optOut.setValue(EMDL_IMAGE_MAG_MATRIX_01, mat1(0,1), 0);
		optOut.setValue(EMDL_IMAGE_MAG_MATRIX_10, mat1(1,0), 0);
		optOut.setValue(EMDL_IMAGE_MAG_MATRIX_11, mat1(1,1), 0);		
	}
	
	if (adaptAstig)
	{
		MagnificationHelper::adaptAstigmatism(mat_by_optGroup, mdts, !perMgAstig, obsModel);
	}
}

bool MagnificationEstimator::isFinished(
		const MetaDataTable &mdt)
{
	if (!ready)
	{
		REPORT_ERROR("ERROR: TiltEstimator::isFinished: DefocusEstimator not initialized.");
	}
	
	std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
	
	return exists(outRoot+"_mag_Axx.mrc")
	    && exists(outRoot+"_mag_Axy.mrc")
		&& exists(outRoot+"_mag_Ayy.mrc")
		&& exists(outRoot+"_mag_bx.mrc")
		&& exists(outRoot+"_mag_by.mrc");
}

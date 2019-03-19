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
#include <src/jaz/img_proc/color_helper.h>

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
		int verb, int nr_omp_threads, 
		bool debug, bool diag, 
		std::string outPath, 
		ReferenceMap* reference, 
		ObservationModel* obsModel)
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
	
	std::vector<std::pair<int, std::vector<int>>> particlesByOpticsGroup 
			= obsModel->splitParticlesByOpticsGroup(mdt);
	
	for (int pog = 0; pog < particlesByOpticsGroup.size(); pog++)
	{			
		const int og = particlesByOpticsGroup[pog].first;
		const std::vector<int>& partIndices = particlesByOpticsGroup[pog].second;
		
		const int pc = partIndices.size();
		
		std::vector<Volume<Equation2x2>> magEqs(nr_omp_threads);
		
		for (int i = 0; i < nr_omp_threads; i++)
		{
			magEqs[i] = Volume<Equation2x2>(sh[og],s[og],1);
		}
		
		#pragma omp parallel for num_threads(nr_omp_threads)
		for (long pp = 0; pp < pc; pp++)
		{
			const int p = partIndices[pp];
			
			CTF ctf;
			ctf.readByGroup(mdt, obsModel, p);
			
			int threadnum = omp_get_thread_num();
			
			MagnificationHelper::updateScaleFreq(
				pred[p], predGradient[p], obs[p], ctf, angpix[og], magEqs[threadnum]);
		}
		
		Volume<Equation2x2> magEq(sh[og], s[og],1);
				
		for (int threadnum = 0; threadnum < nr_omp_threads; threadnum++)
		{
			magEq += magEqs[threadnum];
		}
		
		std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdt, outPath);
		
		std::stringstream sts;
		sts << (og+1);
		
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
	
	bool hasMagMatrices = optOut.labelExists(EMDL_IMAGE_MAG_MATRIX_00)
	                   && optOut.labelExists(EMDL_IMAGE_MAG_MATRIX_01)
	                   && optOut.labelExists(EMDL_IMAGE_MAG_MATRIX_10)
	                   && optOut.labelExists(EMDL_IMAGE_MAG_MATRIX_11);
			
	std::vector<Matrix2D<RFLOAT>> mat_by_optGroup(ogc);

	#pragma omp parallel for num_threads(nr_omp_threads)	
	for (int og = 0; og < ogc; og++)
	{
		Volume<Equation2x2> magEqs(sh[og],s[og],1), magEqsG(sh[og],s[og],1);
		
		std::stringstream sts;
		sts << (og+1);
		
		bool groupPresent = false;
		
		for (long g = 0; g < gc; g++)
		{
			std::string outRoot = CtfRefiner::getOutputFilenameRoot(mdts[g], outPath);
			
			std::string fn = outRoot + "_mag_optics-group_" + sts.str();
			
			if (exists(fn+"_Axx.mrc")
			 || exists(fn+"_Axy.mrc")
			 || exists(fn+"_Ayy.mrc")
			 || exists(fn+"_bx.mrc")
			 || exists(fn+"_by.mrc"))
			{
				try
				{
					MagnificationHelper::readEQs(fn, magEqsG);
					magEqs += magEqsG;
					groupPresent = true;
				}
				catch (RelionError e)
				{}			
			}
		}
		
		if (!groupPresent) 
		{
			mat_by_optGroup[og] = Matrix2D<RFLOAT>(2,2);
			mat_by_optGroup[og].initIdentity();
			continue;
		}
		
		Image<RFLOAT> flowx, flowy;
		MagnificationHelper::solvePerPixel(magEqs, flowx, flowy);
			
		Image<RFLOAT> flowxFull, flowyFull;
		FftwHelper::decenterUnflip2D(flowx.data, flowxFull.data);
		FftwHelper::decenterUnflip2D(flowy.data, flowyFull.data);
		
		ImageLog::write(flowxFull, outPath + "mag_disp_x_optics-group_" + sts.str());
		ImageLog::write(flowyFull, outPath + "mag_disp_y_optics-group_" + sts.str());
		ColorHelper::writeSignedToPNG(flowxFull, 
			outPath + "mag_disp_x_optics-group_" + sts.str() + "_[-1,+1]");
		ColorHelper::writeSignedToPNG(flowyFull, 
			outPath + "mag_disp_y_optics-group_" + sts.str() + "_[-1,+1]");
		
		Image<RFLOAT> freqWght = reference->getHollowWeight(kmin, s[og], angpix[og]);
		
		Matrix2D<RFLOAT> mat = MagnificationHelper::solveLinearlyFreq(magEqs, freqWght, flowx, flowy);
		
		FftwHelper::decenterUnflip2D(flowx.data, flowxFull.data);
		FftwHelper::decenterUnflip2D(flowy.data, flowyFull.data);
		
		ImageLog::write(flowxFull, outPath + "mag_disp_x_fit_optics-group_" + sts.str());
		ImageLog::write(flowyFull, outPath + "mag_disp_y_fit_optics-group_" + sts.str());
		ColorHelper::writeSignedToPNG(flowxFull, 
			outPath + "mag_disp_x_fit_optics-group_" + sts.str() + "_[-1,+1]");
		ColorHelper::writeSignedToPNG(flowyFull, 
			outPath + "mag_disp_y_fit_optics-group_" + sts.str() + "_[-1,+1]");
		
		std::ofstream os(outPath + "mag_matrix_optics-group_" + sts.str() + ".txt");
		os << mat(0,0) << " " << mat(0,1) << "\n";
		os << mat(1,0) << " " << mat(1,1) << "\n";
		os.close();
		
		mat_by_optGroup[og] = mat;
		
		Matrix2D<RFLOAT> mat0 = obsModel->getMagMatrix(og);
		Matrix2D<RFLOAT> mat1 = mat * mat0;
		
		#pragma omp critical
		{
			optOut.setValue(EMDL_IMAGE_MAG_MATRIX_00, mat1(0,0), og);
			optOut.setValue(EMDL_IMAGE_MAG_MATRIX_01, mat1(0,1), og);
			optOut.setValue(EMDL_IMAGE_MAG_MATRIX_10, mat1(1,0), og);
			optOut.setValue(EMDL_IMAGE_MAG_MATRIX_11, mat1(1,1), og);
		}

		obsModel->setMagMatrix(og, mat1);
	}
	
	obsModel->hasMagMatrices = true;
	
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

	bool allThere = true;
	
	std::vector<int> ogp = obsModel->getOptGroupsPresent_zeroBased(mdt);
	
	for (int pog = 0; pog < ogp.size(); pog++)
	{
		const int og = ogp[pog];
		
		std::stringstream sts;
		sts << (og+1);
		std::string ogstr = sts.str();
		
		if(  !exists(outRoot+"_mag_optics-group_"+ogstr+"_Axx.mrc")
		  || !exists(outRoot+"_mag_optics-group_"+ogstr+"_Axy.mrc")
		  || !exists(outRoot+"_mag_optics-group_"+ogstr+"_Ayy.mrc")
		  || !exists(outRoot+"_mag_optics-group_"+ogstr+"_bx.mrc")
		  || !exists(outRoot+"_mag_optics-group_"+ogstr+"_by.mrc"))
		{
			allThere = false;
			break;
		}
	}
	
	return allThere;
}

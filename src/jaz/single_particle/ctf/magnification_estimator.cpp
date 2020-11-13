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
#include "ctf_refiner.h"

#include <src/jaz/single_particle/reference_map.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/image/color_helper.h>
#include <src/jaz/optics/magnification_helper.h>
#include <src/jaz/math/equation2x2.h>

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
		const std::vector<Volume<t2Vector<Complex>>>& predGradient,
		bool do_ctf_padding)
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

		// TODO: SHWS 29mar2018: when data is CTF-premultiplied: do we need to change updateScaleFreq??
		if (obsModel->getCtfPremultiplied(og))
			std::cerr << "TODO: check magnification correction with CTF-premultiplied data!!" << std::endl;

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
				pred[p], predGradient[p], obs[p], ctf, angpix[og], magEqs[threadnum], do_ctf_padding);
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
		std::vector<MetaDataTable> &mdts, MetaDataTable &optOut, std::vector<FileName> &fn_eps)
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

		std::vector<Image<RFLOAT> > imgs_for_eps;
		std::vector<double> scales;
		std::vector<std::string> labels;

		Image<RFLOAT> flowx, flowy;
		MagnificationHelper::solvePerPixel(magEqs, flowx, flowy);

		Image<RFLOAT> flowxFull, flowyFull;
		FftwHelper::decenterUnflip2D(flowx.data, flowxFull.data);
		FftwHelper::decenterUnflip2D(flowy.data, flowyFull.data);

		ImageLog::write(flowxFull, outPath + "mag_disp_x_optics-group_" + sts.str());
		ImageLog::write(flowyFull, outPath + "mag_disp_y_optics-group_" + sts.str());

		imgs_for_eps.push_back(flowxFull);
		scales.push_back(1.);
		labels.push_back("X-disp obs [-1,+1] "+obsModel->getGroupName(og));

		imgs_for_eps.push_back(flowyFull);
		scales.push_back(1.);
		labels.push_back("Y-disp obs [-1,+1] "+obsModel->getGroupName(og));


		Image<RFLOAT> freqWght = reference->getHollowWeight(kmin, s[og], angpix[og]);

		Matrix2D<RFLOAT> mat = MagnificationHelper::solveLinearlyFreq(magEqs, freqWght, flowx, flowy);

		FftwHelper::decenterUnflip2D(flowx.data, flowxFull.data);
		FftwHelper::decenterUnflip2D(flowy.data, flowyFull.data);

		ImageLog::write(flowxFull, outPath + "mag_disp_x_fit_optics-group_" + sts.str());
		ImageLog::write(flowyFull, outPath + "mag_disp_y_fit_optics-group_" + sts.str());


		imgs_for_eps.push_back(flowxFull);
		scales.push_back(1.);
		labels.push_back("X-disp fit [-1,+1] "+obsModel->getGroupName(og));

		imgs_for_eps.push_back(flowyFull);
		scales.push_back(1.);
		labels.push_back("Y-disp fit [-1,+1] "+obsModel->getGroupName(og));

#ifdef DEBUG
		std::ofstream os(outPath + "mag_matrix_optics-group_" + sts.str() + ".txt");
		os << mat(0,0) << " " << mat(0,1) << "\n";
		os << mat(1,0) << " " << mat(1,1) << "\n";
		os.close();
#endif
		mat_by_optGroup[og] = mat;

		Matrix2D<RFLOAT> mat0 = obsModel->getMagMatrix(og);
		Matrix2D<RFLOAT> mat1 = mat * mat0;

		Matrix2D<RFLOAT> u, vh;
		Matrix1D<RFLOAT> eig;
		svdcmp(mat1, u, eig, vh);
		const RFLOAT mean_mag = (eig(0) + eig(1)) / 2;
		const RFLOAT aniso_mag = fabs(eig(0) - eig(1));// / mean_mag;
		
		#pragma omp critical
		{
			optOut.setValue(EMDL_IMAGE_MAG_MATRIX_00, mat1(0,0), og);
			optOut.setValue(EMDL_IMAGE_MAG_MATRIX_01, mat1(0,1), og);
			optOut.setValue(EMDL_IMAGE_MAG_MATRIX_10, mat1(1,0), og);
			optOut.setValue(EMDL_IMAGE_MAG_MATRIX_11, mat1(1,1), og);

			std::cout << "   - Magnification anisotropy of optics group #" << (og + 1) << " named '" << obsModel->getGroupName(og) << "': " << (aniso_mag * 100.0) <<  " %" << std::endl;
			if (fabs(mean_mag - 1) > 0.005)
			{
				std::cerr << "WARNING: Overall magnification of optics group #" << (og + 1) << " (" << obsModel->getGroupName(og) << ") differs from the nominal pixel size by " << ((mean_mag - 1) * 100) << " %.\n";
				std::cerr << "WARNING: This overall difference changes the actual pixel size of the reconstruction!" << std::endl;
			}
		}

		obsModel->setMagMatrix(og, mat1);

		FileName fn_root = outPath + "mag_disp_optics-group_"+ sts.str();
		ColorHelper::writeSignedToEPS(fn_root, 2, imgs_for_eps, scales, labels);
		fn_eps.push_back(fn_root+".eps");
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

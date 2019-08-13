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

#include "frame_recombiner.h"
#include "motion_refiner.h"
#include "motion_helper.h"

#include <src/jaz/micrograph_handler.h>
#include <src/jaz/obs_model.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/vtk_helper.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/filename.h>

using namespace gravis;

FrameRecombiner::FrameRecombiner()
{}

void FrameRecombiner::read(IOParser& parser, int argc, char* argv[])
{
	parser.addSection("Combine frames options");

	doCombineFrames = parser.checkOption("--combine_frames", "Combine movie frames into polished particles.");
	scale_arg = textToInteger(parser.getOption("--scale", "Re-scale the particles to this size (by default read from particles star file)", "-1"));
	box_arg = textToInteger(parser.getOption("--window", "Re-window the particles to this size (in movie-pixels; by default read from particles star file)", "-1"));
	crop_arg = textToInteger(parser.getOption("--crop", "Crop the scaled particles to this size after CTF pre-multiplication", "-1"));
	do_ctf_multiply = parser.checkOption("--ctf_multiply", "Premultiply by CTF.");
	k0a = textToDouble(parser.getOption("--bfac_minfreq", "Min. frequency used in B-factor fit [Angst]", "20"));
	k1a = textToDouble(parser.getOption("--bfac_maxfreq", "Max. frequency used in B-factor fit [Angst]", "-1"));
	bfacFn = parser.getOption("--bfactors", "A .star file with external B/k-factors", "");
	bfac_diag = parser.checkOption("--diag_bfactor", "Write out B/k-factor diagnostic data");
	suffix = parser.getOption("--suffix", "Add this suffix to shiny MRCS and STAR files", "");

	if (box_arg > 0 || scale_arg > 0)
	{
		std::cerr << "WARNING: Changing the box size (--window and/or --scale) might "
		          << "invalidate the current particle offsets.\nPlease remember to "
		          << "run relion_refine again." << std::endl;
	}

	if (box_arg > 0 && box_arg % 2 != 0)
	{
		REPORT_ERROR_STR("The window size (--window) has to be an even number.\n");
	}

	if (scale_arg > 0 && scale_arg % 2 != 0)
	{
		REPORT_ERROR_STR("The rescaled window size (--scale) has to be an even number.\n");
	}

	if (crop_arg > 0 && !do_ctf_multiply)
	{
		REPORT_ERROR("--crop is meaningless without --ctf_multiply");
	}

	if (box_arg > 0 && box_arg % 2 != 0)
	{
		REPORT_ERROR("--window must be an even number");
	}

	if (scale_arg > 0 && scale_arg % 2 != 0)
	{
		REPORT_ERROR("--scale must be an even number");
	}

	if (crop_arg > 0 && crop_arg % 2 != 0)
	{
		REPORT_ERROR("--crop must be an even number");
	}
}

void FrameRecombiner::init(
	const std::vector<MetaDataTable>& allMdts,
	int verb, int s_ref, int fc, double maxFreq, double angpix_ref,
	int nr_omp_threads, std::string outPath, bool debug,
	ReferenceMap* reference,
	ObservationModel* obsModel,
	MicrographHandler* micrographHandler)
{
	this->verb = verb;
	this->s_ref = s_ref;
	this->sh_ref = s_ref/2 + 1;
	this->fc = fc;
	this->nr_omp_threads = nr_omp_threads;
	this->outPath = outPath;
	this->debug = debug;
	this->reference = reference;
	this->obsModel = obsModel;
	this->micrographHandler = micrographHandler;
	this->angpix_ref = angpix_ref;
	this->maxFreq = maxFreq;

	/*
	  OLD:

		neither window nor scale provided:

			angpix_out = angpix_ref
			s_out = s_ref

		only window provided:

			angpix_out = angpix_ref
			s_out = window * angpix_mov / angpix_ref

		only scale provided:

			window = s_ref * angpix_mov / angpix_ref
			angpix_out = angpix_mov * window / scale
			s_out = scale

		both provided:

			angpix_out = angpix_mov * window / scale
			s_out = scale

	  NEW:

		neither window nor scale provided:

			angpix_out = angpix(particle)
			s_out = box_size(particle)

		only window provided:

			angpix_out = angpix(particle)
			s_out = window * angpix_mov / angpix(particle)

		only scale provided:

			window_mov = box_size(particle) * angpix(particle) / angpix_mov
			angpix_out = angpix_mov * window_mov / scale
			s_out = scale

		both provided:

			angpix_out = angpix_mov * window / scale
			s_out = scale
	*/

	const int nog = obsModel->numberOfOpticsGroups();

	s_mov.resize(nog);
	s_out.resize(nog);
	sh_out.resize(nog);
	angpix_out.resize(nog);
	freqWeights.resize(nog);

	const double angpix_mov = micrographHandler->movie_angpix; // TODO: TAKANORI: make sure the movie_angpix is the same for all micrographs

	// The optics group is used only to account for different pixel sizes.
	// Dose weighting uses information from all optics groups.
	for (int og = 0; og < nog; og++)
	{
		if (box_arg > 0) s_mov[og] = box_arg;
		else s_mov[og] = 2 * (int)(0.5 * s_ref * angpix_ref / angpix_mov + 0.5);

		if (scale_arg > 0)
		{
			s_out[og] = scale_arg;
			angpix_out[og] = angpix_mov * s_mov[og] / (double) scale_arg;
		}
		else
		{
			s_out[og] = 2 * (int)(0.5 * s_mov[og] * angpix_mov / angpix_ref + 0.5);
			angpix_out[og] = angpix_ref;
		}

		sh_out[og] = s_out[og]/2 + 1;

		if (debug)
		{
			std::cout << "optics group " << og << "\n";
			std::cout << "s_out: " << s_out[og] << "\n";
			std::cout << "s_mov: " << s_mov[og] << "\n";
			std::cout << "angpix_out: " << angpix_out[og] << "\n";
		}

		if (s_out[og] > s_mov[og])
		{
			REPORT_ERROR_STR("Images can only be scaled down, not up!\n"
			                 << "You are trying to extract squares of size " << s_mov[og] << " px from the movies and "
			                 << "scale them up to " << s_out[og] << " px\n");
		}

		// Either calculate weights from FCC or from user-provided B-factors
		const bool hasBfacs = bfacFn != "";
		std::stringstream sts;
		sts << "optics-group-" << (og + 1);

		if (!hasBfacs)
		{
			freqWeights[og] = weightsFromFCC(allMdts, s_out[og], angpix_out[og], sts.str());
		}
		else
		{
			freqWeights[og] = weightsFromBfacs(allMdts, s_out[og], angpix_out[og]);
		}
	}
}

void FrameRecombiner::process(const std::vector<MetaDataTable>& mdts, long g_start, long g_end)
{
	int barstep;
	int my_nr_micrographs = g_end - g_start + 1;

	if (verb > 0)
	{
		std::cout << " + Combining frames for all micrographs ... " << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
	}

	std::vector<ParFourierTransformer> fts(nr_omp_threads);

	int pctot = 0;
	long nr_done = 0;

	for (long g = g_start; g <= g_end; g++)
	{
		// Abort through the pipeline_control system, TODO: check how this goes with MPI....
		if (pipeline_control_check_abort_job())
			exit(RELION_EXIT_ABORTED);

		const int pc = mdts[g].numberOfObjects();
		if (pc == 0) continue;

		pctot += pc;

		// optics group representative of this micrograph
		// (only the pixel and box sizes have to be identical)
		int ogmg = obsModel->getOpticsGroup(mdts[g], 0);

		if (!obsModel->allPixelAndBoxSizesIdentical(mdts[g]))
		{
			std::cerr << "WARNING: varying pixel or box sizes detected in "
			          << MotionRefiner::getOutputFileNameRoot(outPath, mdts[g])
			          << " - skipping micrograph." << std::endl;

			continue;
		}


		FileName fn_root = MotionRefiner::getOutputFileNameRoot(outPath, mdts[g]);
		std::vector<std::vector<d2Vector>> shift0;
		shift0 = MotionHelper::readTracksInPix(fn_root + "_tracks.star", angpix_out[ogmg]);

		std::vector<std::vector<d2Vector>> shift = shift0;

		std::vector<std::vector<Image<Complex>>> movie;

		// loadMovie() will extract squares around the value of shift0 rounded in movie coords,
		// and return the remainder in shift (in output coordinates)
		movie = micrographHandler->loadMovie(mdts[g], s_out[ogmg], angpix_out[ogmg], fts, &shift0, &shift);

		const int out_size = crop_arg > 0 ? crop_arg : s_out[ogmg];
		Image<RFLOAT> stack(out_size, out_size, 1, pc);

		#pragma omp parallel for num_threads(nr_omp_threads)
		for (int p = 0; p < pc; p++)
		{
			int threadnum = omp_get_thread_num();

			Image<Complex> sum(sh_out[ogmg], s_out[ogmg]);
			sum.data.initZeros();

			Image<Complex> obs(sh_out[ogmg], s_out[ogmg]);

			for (int f = 0; f < fc; f++)
			{
				shiftImageInFourierTransform(movie[p][f](), obs(), s_out[ogmg],
				                             -shift[p][f].x, -shift[p][f].y);

				for (int y = 0; y < s_out[ogmg]; y++)
				for (int x = 0; x < sh_out[ogmg]; x++)
				{
					sum(y,x) += freqWeights[ogmg][f](y,x) * obs(y,x);
				}
			}

			Image<RFLOAT> real(s_out[ogmg], s_out[ogmg]);

			// Premultiply by CTF
			if (do_ctf_multiply)
			{
				CTF ctf;
				ctf.readByGroup(mdts[g], obsModel, p);
				int og = obsModel->getOpticsGroup(mdts[g], p);
	 			if (obsModel->getBoxSize(og) != s_out[og])
				{
					obsModel->setBoxSize(og, s_out[og]);
				}

				MultidimArray<RFLOAT> Fctf;
				Fctf.resize(YSIZE(sum()), XSIZE(sum()));
				ctf.getFftwImage(Fctf, s_out[og], s_out[og], angpix_out[og], false, false, false, true, false);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sum())
				{
					 DIRECT_MULTIDIM_ELEM(sum(), n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}
			}

			fts[threadnum].inverseFourierTransform(sum(), real());
			real().setXmippOrigin();

			const int half_out = out_size / 2;
			for (int y = 0; y < out_size; y++)
			for (int x = 0; x < out_size; x++)
			{
				DIRECT_NZYX_ELEM(stack(), p, 0, y, x) = real(y - half_out, x - half_out); // Image() is logical access
			}
		}

		stack.write(fn_root+"_shiny" + suffix + ".mrcs");

		if (debug)
		{
			VtkHelper::writeTomoVTK(
				stack, fn_root+"_shiny" + suffix + ".vtk", false,
				angpix_out[ogmg],
				-angpix_out[ogmg] * s_out[ogmg] * 0.5 * d3Vector(1,1,0));
		}

		MetaDataTable mdtOut = mdts[g];

		for (int p = 0; p < pc; p++)
		{
			std::stringstream sts;
			sts << (p+1);
			mdtOut.setValue(EMDL_IMAGE_NAME, sts.str() + "@" + fn_root+"_shiny" + suffix + ".mrcs", p);
		}

		mdtOut.write(fn_root+"_shiny" + suffix + ".star");

		nr_done++;

		if (verb > 0 && nr_done % barstep == 0)
		{
			progress_bar(nr_done);
		}
	}

	if (verb > 0)
	{
		progress_bar(my_nr_micrographs);
	}
}

std::vector<Image<RFLOAT>> FrameRecombiner::weightsFromFCC(
		const std::vector<MetaDataTable>& allMdts,
		int s, double angpix, std::string og_name)
{
	if (debug && verb > 0)
	{
		std::cout << " + Summing up FCCs..." << std::endl;
	}

	Image<RFLOAT> fccData, fccWgh0, fccWgh1;
	Image<RFLOAT> fccDataMg, fccWgh0Mg, fccWgh1Mg;

	bool first = true;

	// Compute B/k-factors from all available FCCs (allMdts),
	// even if only a subset of micrographs (chosenMdts) is being recombined.
	for (long g = 0; g < allMdts.size(); g++)
	{
		FileName fn_root = MotionRefiner::getOutputFileNameRoot(outPath, allMdts[g]);

		if (!( exists(fn_root + "_FCC_cc.mrc")
			&& exists(fn_root + "_FCC_w0.mrc")
			&& exists(fn_root + "_FCC_w1.mrc")))
		{
			continue;
		}

		fccDataMg.read(fn_root + "_FCC_cc.mrc");
		fccWgh0Mg.read(fn_root + "_FCC_w0.mrc");
		fccWgh1Mg.read(fn_root + "_FCC_w1.mrc");

		if (first)
		{
			sh_ref = fccDataMg.data.xdim;
			s_ref = 2 * (sh_ref-1);

			fc = fccDataMg.data.ydim;

			fccData = Image<RFLOAT>(sh_ref,fc);
			fccWgh0 = Image<RFLOAT>(sh_ref,fc);
			fccWgh1 = Image<RFLOAT>(sh_ref,fc);

			fccData.data.initZeros();
			fccWgh0.data.initZeros();
			fccWgh1.data.initZeros();

			first = false;
		}

		for (int y = 0; y < fc; y++)
		for (int x = 0; x < sh_ref; x++)
		{
			if (fccDataMg(y,x) == fccDataMg(y,x)) fccData(y,x) += fccDataMg(y,x);
			if (fccWgh0Mg(y,x) == fccWgh0Mg(y,x)) fccWgh0(y,x) += fccWgh0Mg(y,x);
			if (fccWgh1Mg(y,x) == fccWgh1Mg(y,x)) fccWgh1(y,x) += fccWgh1Mg(y,x);
		}
	}

	Image<RFLOAT> fcc(sh_ref,fc);

	for (int y = 0; y < fc; y++)
	for (int x = 0; x < sh_ref; x++)
	{
		const double wgh = sqrt(fccWgh0(y,x) * fccWgh1(y,x));

		if (wgh > 0.0)
		{
			fcc(y,x) = fccData(y,x) / wgh;
		}
		else
		{
			fcc(y,x) = 0.0;
		}
	}

	if (debug) std::cout << "done\n";

	k0 = (int) reference->angToPix(k0a);

	if (!outerFreqKnown())
	{
		k1a = maxFreq;
	}

	k1 = (int) reference->angToPix(k1a);

	if (verb > 0)
	{
		std::cout << " + Fitting B/k-factors for " << og_name << " using FCCs from all particles between " << k0 << " and " << k1 << " pixels, or "
		          << k0a << " and " << k1a << " Angstrom ..." << std::endl;
	}

	std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs = DamageHelper::fitBkFactors(fcc, k0, k1);

	// sigmas (bkFacs[f].x) are given in pixels;
	// rescale if a different box size is to be extracted
	std::vector<d2Vector> bkFacsRescaled(fc);

	if (s == s_ref)
	{
		bkFacsRescaled = bkFacs.first;
	}
	else
	{
		for (int f = 0; f < fc; f++)
		{
			bkFacsRescaled[f].x = (s * angpix) * bkFacs.first[f].x / (double) (s_ref * angpix_ref);
			bkFacsRescaled[f].y = bkFacs.first[f].y;
		}
	}

	const int sh = s/2 + 1;

	std::vector<Image<RFLOAT>> freqWeights;
	freqWeights = DamageHelper::computeWeights(bkFacsRescaled, sh);

	const double cf = 8.0 * angpix_ref * angpix_ref * sh_ref * sh_ref;

	if (bfac_diag)
	{
		mktree(outPath + "/bfacs");

		Image<RFLOAT> bfacFit = DamageHelper::renderBkFit(bkFacs, sh_ref, fc);
		Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs, sh_ref, fc, true);

		ImageLog::write(bfacFit, outPath + "/bfacs/glob_Bk-fit");
		ImageLog::write(bfacFitNoScale, outPath + "/bfacs/glob_Bk-fit_noScale");
		ImageLog::write(fcc, outPath + "/bfacs/glob_Bk-data");
		ImageLog::write(freqWeights, outPath + "/bfacs/freqWeights");

		std::ofstream bfacsDat(outPath + "/bfacs/Bfac.dat");
		std::ofstream kfacsDat(outPath + "/bfacs/kfac.dat");

		for (int i = 0; i < fc; i++)
		{
			double sig = bkFacs.first[i].x;
			double b = -cf/(sig*sig);

			bfacsDat << i << " " << b << std::endl;
			kfacsDat << i << " " << log(bkFacs.first[i].y) << std::endl;
		}

		bfacsDat.close();
		kfacsDat.close();
	}

	MetaDataTable mdt;
	mdt.setName("perframe_bfactors");

	for (int f = 0; f < fc; f++ )
	{
		double sig = bkFacs.first[f].x;
		double b = -cf/(sig*sig);
		double k = log(bkFacs.first[f].y);

		mdt.addObject();
		mdt.setValue(EMDL_IMAGE_FRAME_NR, f);
		mdt.setValue(EMDL_POSTPROCESS_BFACTOR, b);
		mdt.setValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, k);
	}

	mdt.write(outPath + "/bfactors.star");

	// Also write out EPS plots of the B-factors and scale factors
	CPlot2D plot2D("Polishing B-factors");
	plot2D.SetXAxisSize(600);
	plot2D.SetYAxisSize(400);
	plot2D.SetDrawLegend(false);
	plot2D.SetXAxisTitle("movie frame");
	plot2D.SetYAxisTitle("B-factor");
	mdt.addToCPlot2D(&plot2D, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_BFACTOR);
	plot2D.OutputPostScriptPlot(outPath + "bfactors.eps");

	CPlot2D plot2Db("Polishing scale-factors");
	plot2Db.SetXAxisSize(600);
	plot2Db.SetYAxisSize(400);
	plot2Db.SetDrawLegend(false);
	plot2Db.SetXAxisTitle("movie frame");
	plot2Db.SetYAxisTitle("Scale-factor");
	mdt.addToCPlot2D(&plot2Db, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT);
	plot2Db.OutputPostScriptPlot(outPath + "scalefactors.eps");

	return freqWeights;
}

std::vector<Image<RFLOAT>> FrameRecombiner::weightsFromBfacs(
		const std::vector<MetaDataTable>& allMdts,
		int s, double angpix)
{
	const int sh = s/2 + 1;

	// initialization on the first line to avoid copying of return value
	std::vector<Image<RFLOAT>> freqWeights;

	MetaDataTable mdt;
	mdt.read(bfacFn);

	fc = mdt.numberOfObjects();

	std::vector<d2Vector> bkFacs(fc);

	double bfacOff = 0.0;

	for (int f = 0; f < fc; f++)
	{
		double b;
		mdt.getValue(EMDL_POSTPROCESS_BFACTOR, b, f);

		if (b > bfacOff) bfacOff = b;
	}

	const double cf = 8.0 * angpix_ref * angpix_ref * sh * sh;

	for (int f = 0; f < fc; f++)
	{
		double b, k;
		mdt.getValue(EMDL_POSTPROCESS_BFACTOR, b, f);
		mdt.getValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, k, f);

		bkFacs[f] = d2Vector(sqrt(-cf/(b-bfacOff-1)), exp(k));
	}

	freqWeights = DamageHelper::computeWeights(bkFacs, sh);

	if (bfac_diag)
	{
		mktree(outPath + "bfacs");

		std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs2;
		bkFacs2.first = bkFacs;
		bkFacs2.second = std::vector<double>(fc, 1.0);

		Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs2, sh, fc, true);

		ImageLog::write(bfacFitNoScale, outPath + "/bfacs/glob_Bk-fit_noScale");
		ImageLog::write(freqWeights, outPath + "/bfacs/freqWeights");

		std::ofstream bfacsDat(outPath + "/bfacs/Bfac.dat");
		std::ofstream kfacsDat(outPath + "/bfacs/kfac.dat");

		for (int i = 0; i < fc; i++)
		{
			double sig = bkFacs[i].x;
			double b = -cf/(sig*sig);

			bfacsDat << i << " " << b << std::endl;
			kfacsDat << i << " " << log(bkFacs[i].y) << std::endl;
		}

		bfacsDat.close();
		kfacsDat.close();
	}

	return freqWeights;
}

bool FrameRecombiner::doingRecombination()
{
	return doCombineFrames;
}

bool FrameRecombiner::outerFreqKnown()
{
	return k1a > 0.0;
}

std::vector<MetaDataTable> FrameRecombiner::findUnfinishedJobs(
		const std::vector<MetaDataTable> &mdts, std::string path)
{
	std::vector<MetaDataTable> out(0);

	const int gc = mdts.size();

	for (int g = 0; g < gc; g++)
	{
		std::string fn_root = MotionRefiner::getOutputFileNameRoot(path, mdts[g]);

		if (!isJobFinished(fn_root))
		{
			out.push_back(mdts[g]);
		}
	}

	return out;
}

double FrameRecombiner::getOutputPixelSize(int opticsGroup)
{
	return angpix_out[opticsGroup];
}

int FrameRecombiner::getOutputBoxSize(int opticsGroup)
{
	if (crop_arg > 0)
		return crop_arg;
	else
		return s_out[opticsGroup];
}

bool FrameRecombiner::isCtfMultiplied(int opticsGroup)
{
	return do_ctf_multiply;
}

std::string FrameRecombiner::getOutputSuffix()
{
	return suffix;
}

bool FrameRecombiner::isJobFinished(std::string filenameRoot)
{
	return exists(filenameRoot+"_shiny" + suffix + ".mrcs")
	    && exists(filenameRoot+"_shiny" + suffix + ".star");
}

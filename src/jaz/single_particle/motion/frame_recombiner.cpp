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

#include <src/jaz/single_particle/micrograph_handler.h>
#include <src/jaz/single_particle/obs_model.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/single_particle/damage_helper.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/image/translation.h>
#include <src/jaz/image/padding.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/util/zio.h>
#include <src/filename.h>

using namespace gravis;

FrameRecombiner::FrameRecombiner()
{}

void FrameRecombiner::read(IOParser& parser, int argc, char* argv[])
{
	parser.addSection("Combine frames options");

	doCombineFrames = parser.checkOption("--combine_frames", "Combine movie frames into polished particles.");
	write_float16  = parser.checkOption("--float16", "Write in half-precision 16 bit floating point numbers (MRC mode 12), instead of 32 bit (MRC mode 0).");
	scale_arg = textToInteger(parser.getOption("--scale", "Re-scale the particles to this size (by default read from particles star file)", "-1"));
	box_arg = textToInteger(parser.getOption("--window", "Re-window the particles to this size (in movie-pixels; by default read from particles star file)", "-1"));
	crop_arg = textToInteger(parser.getOption("--crop", "Crop the scaled particles to this size after CTF pre-multiplication", "-1"));
	do_ctf_multiply = parser.checkOption("--ctf_multiply", "Premultiply by CTF.");
	k0a = textToDouble(parser.getOption("--bfac_minfreq", "Min. frequency used in B-factor fit [Angst]", "20"));
	k1a = textToDouble(parser.getOption("--bfac_maxfreq", "Max. frequency used in B-factor fit [Angst]", "-1"));
	bfacFn = parser.getOption("--bfactors", "A .star file with external B/k-factors", "");
	bfac_diag = parser.checkOption("--diag_bfactor", "Write out B/k-factor diagnostic data");
	suffix = parser.getOption("--suffix", "Add this suffix to shiny MRCS and STAR files", "");

	do_recenter = parser.checkOption("--recenter", "Re-center particle according to rlnOriginX/Y in --reextract_data_star STAR file");
	recenter_x = textToFloat(parser.getOption("--recenter_x", "X-coordinate (in pixel inside the reference) to recenter re-extracted data on", "0."));
	recenter_y = textToFloat(parser.getOption("--recenter_y", "Y-coordinate (in pixel inside the reference) to recenter re-extracted data on", "0."));
	recenter_z = textToFloat(parser.getOption("--recenter_z", "Z-coordinate (in pixel inside the reference) to recenter re-extracted data on", "0."));

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
	data_angpix.resize(nog);
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
		data_angpix[og] = obsModel->getPixelSize(og);

		if (debug)
		{
			std::cout << "optics group " << og << "\n";
			std::cout << "s_out: " << s_out[og] << "\n";
			std::cout << "s_mov: " << s_mov[og] << "\n";
			std::cout << "angpix_out: " << angpix_out[og] << "\n";
			std::cout << "data_angpix: " << data_angpix[og] << "\n";
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
	const RFLOAT ref_angpix = reference->angpix;
	const RFLOAT coords_angpix = micrographHandler->coords_angpix;

	if (verb > 0)
	{
		std::cout << " + Combining frames for micrographs ... " << std::endl;
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
		{
			exit(RELION_EXIT_ABORTED);
		}

		const int pc = mdts[g].numberOfObjects();
		if (pc == 0) continue;

		pctot += pc;

		MetaDataTable mdtOut = mdts[g];

		// optics group representative of this micrograph
		// (only the pixel and box sizes have to be identical)
		int ogmg = obsModel->getOpticsGroup(mdtOut, 0);

		if (!obsModel->allPixelAndBoxSizesIdentical(mdtOut))
		{
			std::cerr << "WARNING: varying pixel or box sizes detected in "
			          << MotionRefiner::getOutputFileNameRoot(outPath, mdtOut)
			          << " - skipping micrograph." << std::endl;

			continue;
		}


		FileName fn_root = MotionRefiner::getOutputFileNameRoot(outPath, mdtOut);
		std::vector<std::vector<d2Vector>> priorShift;
		priorShift = MotionHelper::readTracksInPix(fn_root + "_tracks.star", angpix_out[ogmg]);

		std::vector<std::vector<d2Vector>> shift(pc);
			
		for (int p = 0; p < pc; p++)
		{
			int og = obsModel->getOpticsGroup(mdtOut, p);
			
			if (obsModel->getBoxSize(og) != s_out[og])
			{
				obsModel->setBoxSize(og, s_out[og]);
			}
			
			if (obsModel->getPixelSize(og) != angpix_out[og])
			{
				obsModel->setPixelSize(og, angpix_out[og]);
			}
			
			shift[p] = {d2Vector(0)};
		}
		
		if (do_recenter)
		{
			recenterParticles(mdtOut, ref_angpix, coords_angpix);
		}
		
		const int out_size = crop_arg > 0 ? crop_arg : s_out[ogmg];
		
		
		BufferedImage<Complex> sumStack(sh_out[ogmg], s_out[ogmg], pc);
		sumStack.fill(Complex(0,0));
		
		
		for (int f = 0; f < fc; f++)
		{
			std::vector<std::vector<gravis::d2Vector>> priorShift_f(pc);
			
			for (int p = 0; p < pc; p++)
			{
				priorShift_f[p] = {priorShift[p][f]};
			}
			
			std::vector<std::vector<Image<Complex>>> fullFrame = micrographHandler->loadMovie(
						mdtOut, s_out[ogmg], angpix_out[ogmg], fts, 
						&priorShift_f, &shift, data_angpix[ogmg], f);
			
			#pragma omp parallel for num_threads(nr_omp_threads)
			for (int p = 0; p < pc; p++)
			{
				RawImage<Complex> obs(fullFrame[p][0]);
				
				Translation::shiftInFourierSpace2D(obs, -shift[p][0].x, -shift[p][0].y);
				
				for (int y = 0; y < s_out[ogmg]; y++)
				for (int x = 0; x < sh_out[ogmg]; x++)
				{
					sumStack(x,y,p) += freqWeights[ogmg][f](y,x) * obs(x,y);
				}
			}
		}
		
		Image<RFLOAT> outStack_xmipp(out_size, out_size, 1, pc);
		RawImage<RFLOAT> outStack(outStack_xmipp);

		#pragma omp parallel for num_threads(nr_omp_threads)
		for (int p = 0; p < pc; p++)
		{
			BufferedImage<Complex> sumCopy = sumStack.getSliceRef(p);
			
			if (do_ctf_multiply)
			{
				CTF ctf;
				ctf.readByGroup(mdtOut, obsModel, p);
				int og = obsModel->getOpticsGroup(mdtOut, p);
				
				MultidimArray<RFLOAT> ctfImg_xmipp;
				ctfImg_xmipp.resize(s_out[og], sh_out[og]);
				
				ctf.getFftwImage(
						ctfImg_xmipp, s_out[og], s_out[og], angpix_out[og], 
						false,  // do_abs
						false,  // do_only_flip_phases
						false,  // do_intact_until_first_peak
						true,   // do_damping
						false );// do_ctf_padding
				
				RawImage<RFLOAT> ctfImg(ctfImg_xmipp);
				
				sumCopy *= ctfImg;
			}
			
			RFLOAT scale2 = StackHelper::computePower(sumCopy, false);
			sumCopy /= out_size * sqrt(scale2);
			
			BufferedImage<RFLOAT> real(s_out[ogmg], s_out[ogmg]);
			FFT::inverseFourierTransform(sumCopy, real);
			
			RawImage<RFLOAT> outSlice = outStack.getSliceRef(p);
			Padding::copyUnpaddedCenter2D_full(real, outSlice);
		}
		
		std::string stackFn = fn_root + "_shiny" + suffix + ".mrcs";
		
		outStack_xmipp.setSamplingRateInHeader(angpix_out[ogmg]);
		outStack_xmipp.write(stackFn, -1, true, WRITE_OVERWRITE, write_float16 ? Float16: Float);

		for (int p = 0; p < pc; p++)
		{
			mdtOut.setValue(EMDL_IMAGE_NAME, ZIO::itoa(p+1) + "@" + stackFn, p);
		}

		mdtOut.write(fn_root + "_shiny" + suffix + ".star");

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

void FrameRecombiner::recenterParticles(
		MetaDataTable& mdtOut,
		RFLOAT ref_angpix,
		RFLOAT coords_angpix)
{
	const int pc = mdtOut.numberOfObjects();
	
	for (int p = 0; p < pc; p++)
	{
		// FIXME: code duplication from preprocess.cpp
		RFLOAT xoff, yoff, xcoord, ycoord;
		Matrix1D<RFLOAT> my_projected_center(3);
		my_projected_center.initZeros();

		mdtOut.getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, xoff, p); // in A
		mdtOut.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, p);

		xoff /= ref_angpix; // Now in reference pixels
		yoff /= ref_angpix;

		if (fabs(recenter_x) > 0. || fabs(recenter_y) > 0. || fabs(recenter_z) > 0.)
		{
			RFLOAT rot, tilt, psi;
			mdtOut.getValue(EMDL_ORIENT_ROT, rot, p);
			mdtOut.getValue(EMDL_ORIENT_TILT, tilt, p);
			mdtOut.getValue(EMDL_ORIENT_PSI, psi, p);

			// Project the center-coordinates
			Matrix1D<RFLOAT> my_center(3);
			Matrix2D<RFLOAT> A3D;
			XX(my_center) = recenter_x; // in reference pixels
			YY(my_center) = recenter_y;
			ZZ(my_center) = recenter_z;
			Euler_angles2matrix(rot, tilt, psi, A3D, false);
			my_projected_center = A3D * my_center;
		}

		xoff -= XX(my_projected_center);
		yoff -= YY(my_projected_center);
		xoff = xoff * ref_angpix / coords_angpix; // Now in (possibly binned) micrograph's pixel
		yoff = yoff * ref_angpix / coords_angpix;

		mdtOut.getValue(EMDL_IMAGE_COORD_X, xcoord, p);
		mdtOut.getValue(EMDL_IMAGE_COORD_Y, ycoord, p);

		xcoord -= ROUND(xoff);
		ycoord -= ROUND(yoff);
		xoff -= ROUND(xoff);
		yoff -= ROUND(yoff);
		
		mdtOut.setValue(EMDL_IMAGE_COORD_X, xcoord, p);
		mdtOut.setValue(EMDL_IMAGE_COORD_Y, ycoord, p);

		mdtOut.setValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, coords_angpix * xoff, p);
		mdtOut.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, coords_angpix * yoff, p);
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
	// TODO: BUG: og_name is not used.
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

std::vector<bool> FrameRecombiner::findUnfinishedJobs(
		const std::vector<MetaDataTable> &mdts, std::string path)
{
	const int gc = mdts.size();

	std::vector<bool> out(gc);

	for (int g = 0; g < gc; g++)
	{
		std::string fn_root = MotionRefiner::getOutputFileNameRoot(path, mdts[g]);

		out[g] = !isJobFinished(fn_root);
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

int FrameRecombiner::getVerbosity()
{
	return verb;
}

void FrameRecombiner::setVerbosity(int v)
{
	verb = v;
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

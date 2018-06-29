/***************************************************************************
 *
 * Authors: "Jasenko Zivanov & Sjors H.W. Scheres"
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

#include "ctf_refiner.h"
#include "tilt_helper.h"

#include <src/jaz/image_log.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/resampling_helper.h>
#include <src/jaz/ctf_helper.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/image_op.h>
#include <src/jaz/parallel_ft.h>

#include <src/ctf.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/time.h>

#include <omp.h>


using namespace gravis;


CtfRefiner::CtfRefiner()
{}

void CtfRefiner::read(int argc, char **argv)
{
	IOParser parser;

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	starFn = parser.getOption("--i", "Input STAR file");

	reference.read(parser, argc, argv);

	outPath = parser.getOption("--o", "Output directory, e.g. CtfRefine/job041/");
	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Skip those steps for which output files already exist.");

	diag = parser.checkOption("--diag", "Write out diagnostic data (slower)");

	int fit_section = parser.addSection("Defocus fit options");
	do_defocus_fit = parser.checkOption("--fit_defocus", "Perform refinement of per-particle defocus values?");
	defocusEstimator.read(parser, argc, argv);

	int tilt_section = parser.addSection("Beam-tilt options");
	do_tilt_fit = parser.checkOption("--fit_beamtilt", "Perform refinement of beamtilt for micrograph groups?");
	tiltEstimator.read(parser, argc, argv);

	int aniso_section = parser.addSection("Anisotropic magnification options");
	do_mag_fit = parser.checkOption("--fit_aniso", "Estimate anisotropic magnification?");
	magnificationEstimator.read(parser, argc, argv);

	int comp_section = parser.addSection("Computational options");
	nr_omp_threads = textToInteger(parser.getOption("--j", "Number of (OMP) threads", "1"));
	minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
	maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));

	debug = parser.checkOption("--debug", "Write debugging data");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));


	int expert_section = parser.addSection("Expert options");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix) - read from STAR file by default", "-1"));
	Cs = textToFloat(parser.getOption("--Cs", "Spherical aberration - read from STAR file by default", "-1"));
	kV = textToFloat(parser.getOption("--kV", "Electron energy (keV) - read from STAR file by default", "-1"));
	beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in X-direction (in mrad)", "0."));
	beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in Y-direction (in mrad)", "0."));
	clTilt = ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.;
	beamtilt_xx = textToFloat(parser.getOption("--beamtilt_xx", "Anisotropic beamtilt, XX-coefficient", "1."));
	beamtilt_xy = textToFloat(parser.getOption("--beamtilt_xy", "Anisotropic beamtilt, XY-coefficient", "0."));
	beamtilt_yy = textToFloat(parser.getOption("--beamtilt_yy", "Anisotropic beamtilt, YY-coefficient", "1."));

	JazConfig::writeMrc = !debug;
	JazConfig::writeVtk = debug;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	// Make sure outPath ends with a slash
	if (outPath[outPath.length()-1] != '/')
	{
		outPath += "/";
	}
}

void CtfRefiner::init()
{
	if (verb > 0)
	{
		std::cout << " + Reading " << starFn << "..." << std::endl;
	}

	mdt0.read(starFn);

	if (!ObservationModel::containsAllNeededColumns(mdt0))
	{
		REPORT_ERROR_STR(starFn << " does not contain all of the required columns ("
			<< "rlnOriginX, rlnOriginY, rlnAngleRot, rlnAngleTilt, rlnAnglePsi and rlnRandomSubset)");
	}

	if (Cs < 0.0)
	{
		mdt0.getValue(EMDL_CTF_CS, Cs, 0);

		if (verb > 0)
		{
			std::cout << "   - Using spherical aberration from the input STAR file: " << Cs << std::endl;
		}
	}
	else
	{
		for (int i = 0; i < mdt0.numberOfObjects(); i++)
		{
			mdt0.setValue(EMDL_CTF_CS, Cs, i);
		}
	}

	if (kV < 0.0)
	{
		mdt0.getValue(EMDL_CTF_VOLTAGE, kV, 0);

		if (verb > 0)
		{
			std::cout << "   - Using voltage from the input STAR file: " << kV << " kV" << std::endl;
		}
	}
	else
	{
		for (int i = 0; i < mdt0.numberOfObjects(); i++)
		{
			mdt0.setValue(EMDL_CTF_VOLTAGE, kV, i);
		}
	}

	if (angpix <= 0.0)
	{
		RFLOAT mag, dstep;
		mdt0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
		mdt0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
		angpix = 10000 * dstep / mag;

		if (verb > 0)
		{
			std::cout << "   - Using pixel size calculated from magnification and detector "
			          << "pixel size in the input STAR file: " << angpix << std::endl;
		}
	}

	allMdts = StackHelper::splitByMicrographName(&mdt0);

	// Only work on a user-specified subset of the micrographs
	if (maxMG < 0 || maxMG >= allMdts.size())
	{
		maxMG = allMdts.size()-1;
	}

	if (minMG < 0 || minMG >= allMdts.size())
	{
		minMG = 0;
	}

	if (minMG > 0 || maxMG < allMdts.size()-1)
	{
		if (verb > 0)
		{
			std::cout << "   - Will only process micrographs in range: ["
					  << minMG << "-" << maxMG << "]"  << std::endl;
		}

		std::vector<MetaDataTable> todo_mdts;

		for (long int g = minMG; g <= maxMG; g++ )
		{
			todo_mdts.push_back(allMdts[g]);
		}

		allMdts = todo_mdts;
	}

	// Make sure output directory ends in a '/'
	if (outPath[outPath.length()-1] != '/')
	{
		outPath+="/";
	}

	RFLOAT V = kV * 1e3;
	lambda = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));
	obsModel = ObservationModel(angpix);

	if (clTilt)
	{
		obsModel = ObservationModel(angpix, Cs, kV * 1e3, beamtilt_x, beamtilt_y);

		if (anisoTilt)
		{
			obsModel.setAnisoTilt(beamtilt_xx, beamtilt_xy, beamtilt_yy);
		}
	}

	if (verb > 0)
	{
		std::cout << " + Reading references ..." << std::endl;
	}

	reference.load(verb, debug);

	// Get dimensions
	s = reference.s;
	sh = s/2 + 1;

	tiltEstimator.init(verb, s, nr_omp_threads, debug, diag, outPath, &reference, &obsModel);
	defocusEstimator.init(verb, s, nr_omp_threads, debug, diag, outPath, &reference, &obsModel);
	magnificationEstimator.init(verb, s, nr_omp_threads, debug, diag, outPath, &reference, &obsModel);

	// check whether output files exist and if they do, then skip this micrograph
	if (only_do_unfinished)
	{
		for (long int g = minMG; g <= maxMG; g++ )
		{
			bool is_done =
				   (!do_tilt_fit || tiltEstimator.isFinished(allMdts[g]))
				&& (!do_defocus_fit || defocusEstimator.isFinished(allMdts[g]))
				&& (!do_mag_fit || magnificationEstimator.isFinished(allMdts[g]));

			if (!is_done)
			{
				unfinishedMdts.push_back(allMdts[g]);
			}
		}
		if (verb > 0)
		{
			std::cout << "   - Will only process " << unfinishedMdts.size()
					  << " unfinished (out of " << allMdts.size()
					  << ") micrographs" << std::endl;
		}
	}
	else
	{
		unfinishedMdts = allMdts;
	}
}

void CtfRefiner::processSubsetMicrographs(long g_start, long g_end)
{
	int barstep;
	int my_nr_micrographs = g_end - g_start + 1;

	if (verb > 0)
	{
		std::cout << " + Performing loop over all micrographs ... " << std::endl;
		init_progress_bar(my_nr_micrographs);
		barstep = XMIPP_MAX(1, my_nr_micrographs/ 60);
	}

	std::vector<ParFourierTransformer> fts(nr_omp_threads);

	long nr_done = 0;
	FileName prevdir = "";

	for (long g = g_start; g <= g_end; g++)
	{
		std::vector<Image<Complex> > obs;

		// both defocus_tit and tilt_fit need the same observations
		obs = StackHelper::loadStackFS(&unfinishedMdts[g], "", nr_omp_threads, &fts, true);

		// Make sure output directory exists
		FileName newdir = getOutputFilenameRoot(unfinishedMdts[g], outPath);
		newdir = newdir.beforeLastOf("/");

		if (newdir != prevdir)
		{
			std::string command = " mkdir -p " + newdir;
			int res = system(command.c_str());
		}

		std::vector<Image<Complex>> predSame, predOpp;

		// use prediction from same half-set for defocus estimation (overfitting danger):
		if (do_defocus_fit)
		{
			predSame = reference.predictAll(
				unfinishedMdts[g], obsModel, ReferenceMap::Own, nr_omp_threads,
				false, true, false);
		}

		// use prediction from opposite half-set otherwise:
		if (do_tilt_fit || do_mag_fit)
		{
			predOpp = reference.predictAll(
				unfinishedMdts[g], obsModel, ReferenceMap::Opposite, nr_omp_threads,
				false, false, false);
		}

		if (do_defocus_fit)
		{
			defocusEstimator.processMicrograph(g, unfinishedMdts[g], obs, predSame);
		}

		if (do_tilt_fit)
		{
			tiltEstimator.processMicrograph(g, unfinishedMdts[g], obs, predOpp);
		}

		if (do_mag_fit)
		{
			magnificationEstimator.processMicrograph(g, unfinishedMdts[g], obs, predOpp);
		}

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

void CtfRefiner::run()
{
	if (do_defocus_fit || do_tilt_fit || do_mag_fit)
	{
		// The subsets will be used in openMPI parallelisation:
		// instead of over g0->gc, they will be over smaller subsets
		processSubsetMicrographs(0, unfinishedMdts.size()-1);
	}

	finalise();
}

void CtfRefiner::finalise()
{
	MetaDataTable mdtOut = mdt0;

	// Read back from disk the metadata-tables and eps-plots for the defocus fit
	// Note: only micrographs for which the defoci were estimated (either now or before)
	// will end up in mdtOut - micrographs excluded through min_MG and max_MG will not.
	if (do_defocus_fit)
	{
		defocusEstimator.merge(allMdts, mdtOut);
	}
	else
	{
		mdtOut = mdt0;
	}

	// Sum up the per-pixel beamtilt fits of all micrographs and fit a parametric model to them.
	// Then, write the beamtilt parameters into mdtOut
	if (do_tilt_fit)
	{
		tiltEstimator.parametricFit(allMdts, Cs, lambda, mdtOut);
	}

	// Do the equivalent for mag. fit
	if (do_mag_fit)
	{
		magnificationEstimator.parametricFit(allMdts, mdtOut);
	}

	mdtOut.write(outPath + "particles_ctf_refine.star");

	if (verb > 0)
	{
		std::cout << " + Done! Written out : " << outPath << "particles_ctf_refine.star" << std::endl;
	}
}

int CtfRefiner::getVerbosityLevel()
{
	return verb;
}

FileName CtfRefiner::getOutputFilenameRoot(const MetaDataTable &mdt, std::string outPath)
{
	FileName fn_mic;
	mdt.getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);
	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);

	return outPath + fn_post.withoutExtension();
}

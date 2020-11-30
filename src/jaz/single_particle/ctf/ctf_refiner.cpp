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

#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/complex_io.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/resampling_helper.h>
#include <src/jaz/single_particle/ctf_helper.h>
#include <src/jaz/single_particle/refinement_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/parallel_ft.h>

#include <src/jaz/util/zio.h>

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
	starFn = parser.getOption("--i", "Input STAR file containing the particles");

	reference.read(parser, argc, argv);

	outPath = parser.getOption("--o", "Output directory, e.g. CtfRefine/job041/");
	only_do_unfinished = parser.checkOption("--only_do_unfinished",
		"Skip those steps for which output files already exist.");

	do_ctf_padding = parser.checkOption("--ctf_pad", "Use larger box to calculate CTF and then downscale to mimic boxing operation in real space");
	diag = parser.checkOption("--diag", "Write out diagnostic data (slower)");

	int fit_section = parser.addSection("Defocus fit options");
	do_defocus_fit = parser.checkOption("--fit_defocus",
		"Perform refinement of per-particle defocus values?");

	defocusEstimator.read(parser, argc, argv);

	int bfac_section = parser.addSection("B-factor options");
	do_bfac_fit = parser.checkOption("--fit_bfacs",
		"Estimate CTF B-factors");

	bfactorEstimator.read(parser, argc, argv);

	int tilt_section = parser.addSection("Beam-tilt options");
	do_tilt_fit = parser.checkOption("--fit_beamtilt",
		"Perform refinement of beamtilt");

	tiltEstimator.read(parser, argc, argv);

	int aberr_section = parser.addSection("Symmetric aberrations options");
	do_aberr_fit = parser.checkOption("--fit_aberr",
		"Estimate symmetric aberrations");

	aberrationEstimator.read(parser, argc, argv);

	int aniso_section = parser.addSection("Anisotropic magnification options");
	do_mag_fit = parser.checkOption("--fit_aniso",
		"Estimate anisotropic magnification");

	magnificationEstimator.read(parser, argc, argv);

	int comp_section = parser.addSection("Computational options");
	nr_omp_threads = textToInteger(parser.getOption("--j", "Number of (OMP) threads", "1"));
	minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
	maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));

	debug = parser.checkOption("--debug", "Write debugging data");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	JazConfig::writeMrc = !debug;
	JazConfig::writeVtk = debug;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
	{
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	}

	// Make sure outPath ends with a slash and exists
	outPath = ZIO::prepareSpaOutputDirectory(outPath);
}

void CtfRefiner::init()
{
	if (verb > 0)
	{
		std::cout << " + Reading " << starFn << "..." << std::endl;
	}

	// Make sure output directory ends in a '/'
	if (outPath[outPath.length()-1] != '/')
	{
		outPath += "/";
	}

	ObservationModel::loadSafely(starFn, obsModel, mdt0);

	if (!ObservationModel::containsAllColumnsNeededForPrediction(mdt0))
	{
		REPORT_ERROR_STR(starFn << " does not contain all columns needed for view prediction: \n"
						 << "rlnOriginXAngst, rlnOriginYAngst, "
						 << "rlnAngleRot, rlnAngleTilt, rlnAnglePsi and rlnRandomSubset");
	}


/* // TAKANORI: TODO Put it somewhere
	if (Cs <= 0.1 && verb > 0)
	{
		std::cerr << "WARNING: Your Cs value is very small. Beam tilt refinement might be unnecessary. Sometimes it gives unrealistically large tilts." << std::endl;
	}
*/

	// after all the necessary changes to mdt0 have been applied
	// in ObservationModel::loadSafely(), split it by micrograph

	allMdts = StackHelper::splitByMicrographName(mdt0);

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


	if (verb > 0)
	{
		std::cout << " + Reading references ..." << std::endl;
	}

	reference.load(verb, debug);

	// Get dimensions
	int s = reference.s;

	tiltEstimator.init(verb, nr_omp_threads, debug, diag, outPath, &reference, &obsModel);
	aberrationEstimator.init(verb, nr_omp_threads, debug, diag, outPath, &reference, &obsModel);
	defocusEstimator.init(verb && do_defocus_fit, nr_omp_threads, debug, diag, outPath, &reference, &obsModel);
	bfactorEstimator.init(verb, nr_omp_threads, debug, diag, outPath, &reference, &obsModel);
	magnificationEstimator.init(verb, nr_omp_threads, debug, diag, outPath, &reference, &obsModel);

	// check whether output files exist and skip the micrographs for which they do
	if (only_do_unfinished)
	{
		for (long int g = minMG; g <= maxMG; g++ )
		{
			bool is_done =
				   (!do_defocus_fit || defocusEstimator.isFinished(allMdts[g]))
				&& (!do_bfac_fit    || bfactorEstimator.isFinished(allMdts[g]))
				&& (!do_tilt_fit    || tiltEstimator.isFinished(allMdts[g]))
				&& (!do_aberr_fit   || aberrationEstimator.isFinished(allMdts[g]))
				&& (!do_mag_fit     || magnificationEstimator.isFinished(allMdts[g]));

			if (!is_done)
			{
				unfinishedMdts.push_back(allMdts[g]);
			}
		}

		if (verb > 0)
		{
			if (unfinishedMdts.size() < allMdts.size())
			{
				std::cout << "   - Will only process " << unfinishedMdts.size()
					  << " unfinished (out of " << allMdts.size()
					  << ") micrographs" << std::endl;
			}
			else
			{
				std::cout << "   - Will process all " << unfinishedMdts.size()
					  << " micrographs" << std::endl;
			}
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
		// Abort through the pipeline_control system, TODO: check how this goes with MPI....
		if (pipeline_control_check_abort_job())
		{
			exit(RELION_EXIT_ABORTED);
		}

		std::vector<Image<Complex> > obs;

		// all CTF-refinement programs need the same observations
		obs = StackHelper::loadStackFS(unfinishedMdts[g], "", nr_omp_threads, true, &obsModel);

		// Make sure output directory exists
		FileName newdir = getOutputFilenameRoot(unfinishedMdts[g], outPath);
		newdir = newdir.beforeLastOf("/");

		if (newdir != prevdir)
		{
			std::string command = " mkdir -p " + newdir;
			int res = system(command.c_str());
		}

		std::vector<Image<Complex>>
				predSameT, // phase-demodulated (defocus)
				predOppNT, // not phase-demodulated (tilt)
				predOppT;  // phase-demodulated (mag and aberr)
		// applyMtf is always true

		// Four booleans in predictAll are applyCtf, applyTilt, applyShift, applyMtf.
		if (do_defocus_fit || do_bfac_fit)
		{
			predSameT = reference.predictAll(
				unfinishedMdts[g], obsModel, ReferenceMap::Own, nr_omp_threads,
				false, true, false, true, do_ctf_padding);
		}

		if (do_tilt_fit)
		{
			predOppNT = reference.predictAll(
				unfinishedMdts[g], obsModel, ReferenceMap::Own, nr_omp_threads,
				false, false, false, true, do_ctf_padding);
		}

		if (do_aberr_fit || do_mag_fit)
		{
			predOppT = reference.predictAll(
				unfinishedMdts[g], obsModel, ReferenceMap::Own, nr_omp_threads,
				false, true, false, true, do_ctf_padding);
		}

		if (do_defocus_fit)
		{
			defocusEstimator.processMicrograph(g, unfinishedMdts[g], obs, predSameT);
		}

		// B-factor fit is always performed after the defocus fit (so it can use the optimal CTFs)
		// The prediction is *not* CTF-weighted, so an up-to-date CTF can be used internally
		if (do_bfac_fit)
		{
			bfactorEstimator.processMicrograph(g, unfinishedMdts[g], obs, predSameT, do_ctf_padding);
		}

		if (do_tilt_fit)
		{
			tiltEstimator.processMicrograph(g, unfinishedMdts[g], obs, predOppNT, do_ctf_padding);
		}

		if (do_aberr_fit)
		{
			aberrationEstimator.processMicrograph(g, unfinishedMdts[g], obs, predOppT);
		}

		if (do_mag_fit)
		{
			std::vector<Volume<t2Vector<Complex>>> predGradient =
				reference.predictAllComplexGradients(
					unfinishedMdts[g], obsModel, ReferenceMap::Opposite, nr_omp_threads,
					false, true, false, true, do_ctf_padding);

			magnificationEstimator.processMicrograph(g, unfinishedMdts[g], obs, predOppT, predGradient, do_ctf_padding);
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
	if (do_defocus_fit || do_bfac_fit || do_tilt_fit || do_aberr_fit || do_mag_fit)
	{
		// The subsets will be used in openMPI parallelisation:
		// instead of over g0->gc, they will be over smaller subsets
		processSubsetMicrographs(0, unfinishedMdts.size()-1);
	}

	finalise();
}

void CtfRefiner::finalise()
{
	std::vector<MetaDataTable> mdtOut;
	std::vector<FileName> fn_eps, fn_eps_earlier, fn_eps_later;

	// Read back from disk the metadata-tables and eps-plots for the B-factor or defocus fit.
	// Note: only micrographs for which the defoci or B-factors were estimated (either now or before)
	// will end up in mdtOut - micrographs excluded through min_MG and max_MG will not.

	if (do_bfac_fit || do_defocus_fit)
	{
		mdtOut = merge(allMdts, fn_eps_later);
	}
	else
	{
		mdtOut = allMdts;
	}

	// ...and for the magnification fit
	if (do_mag_fit)
	{
		magnificationEstimator.parametricFit(mdtOut, obsModel.opticsMdt, fn_eps_earlier);
	}

	// Sum up the per-pixel beamtilt fits of all micrographs and fit a parametric model to them.
	// Then, write the beamtilt parameters into obsModel.opticsMdt
	if (do_tilt_fit)
	{
		tiltEstimator.parametricFit(mdtOut, obsModel.opticsMdt, fn_eps_earlier);
	}

	// Do the equivalent for the symmetrical aberrations...
	if (do_aberr_fit)
	{
		aberrationEstimator.parametricFit(mdtOut, obsModel.opticsMdt, fn_eps_earlier);
	}

	// Make sure to have the EPS from fn_eps_earlier before the ones from fn_eps_later
	// Sort earlier fn_eps, because openMP may have resulted in random order in optics groups
	std::sort(fn_eps_earlier.begin(), fn_eps_earlier.end());
	fn_eps.reserve( fn_eps_earlier.size() + fn_eps_later.size() ); // preallocate memory
	fn_eps.insert( fn_eps.end(), fn_eps_earlier.begin(), fn_eps_earlier.end() );
	fn_eps.insert( fn_eps.end(), fn_eps_later.begin(), fn_eps_later.end() );
	if (fn_eps.size() > 0)
	{
		joinMultipleEPSIntoSinglePDF(outPath + "logfile.pdf", fn_eps);
		if (verb > 0)
		{
			std::cout << " + Written out : " << outPath << "logfile.pdf" << std::endl;
		}
	}

	MetaDataTable mdtOutAll = StackHelper::merge(mdtOut);

	obsModel.save(mdtOutAll, outPath + "particles_ctf_refine.star");

	if (verb > 0)
	{
		std::cout << " + Done! Written out : " << outPath << "particles_ctf_refine.star" << std::endl;
	}
}

std::vector<MetaDataTable> CtfRefiner::merge(const std::vector<MetaDataTable>& mdts, std::vector<FileName> &fn_eps )
{
	int gc = mdts.size();
	int barstep;

	if (verb > 0)
	{
		std::cout << " + Combining data for all micrographs " << std::endl;
		init_progress_bar(gc);
		barstep = 1;
	}

	std::vector<MetaDataTable> mdtOut;

	for (long g = 0; g < gc; g++)
	{
		std::string outRoot = getOutputFilenameRoot(mdts[g], outPath);

		MetaDataTable mdt;

		// If a B-factor fit has been performed, then this has been done after a potential defocus fit,
		// so the B-factor fit files are always more up-to-date.
		if (do_bfac_fit)
		{
			// Read in STAR file with B-factor fit data
			mdt.read(outRoot+"_bfactor_fit.star");
		}
		else if (do_defocus_fit)
		{
			// Read in STAR file with defocus fit data
			mdt.read(outRoot+"_defocus_fit.star");
		}

		mdtOut.push_back(mdt);

		if (exists(outRoot+"_ctf-refine_fit.eps"))
		{
			fn_eps.push_back(outRoot+"_ctf-refine_fit.eps");
		}

		if (exists(outRoot+"_bfactor_fit.eps"))
		{
			fn_eps.push_back(outRoot+"_bfactor_fit.eps");
		}

		if (verb > 0)
		{
			progress_bar(g);
		}
	}

	if (verb > 0)
	{
		progress_bar(gc);
	}

	return mdtOut;
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

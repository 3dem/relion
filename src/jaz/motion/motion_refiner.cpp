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

#include "motion_refiner.h"

#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/image_log.h>
#include <src/jaz/slice_helper.h>
#include <src/jaz/spectral_helper.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/backprojection_helper.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/fsc_helper.h>
#include <src/jaz/image_op.h>
#include <src/jaz/parallel_ft.h>

#include "gp_motion_fit.h"
#include "motion_helper.h"

using namespace gravis;

MotionRefiner::MotionRefiner()
:   motionParamEstimator(),
    motionEstimator(),
    frameRecombiner()
{
}

void MotionRefiner::read(int argc, char **argv)
{
    IOParser parser;
    parser.setCommandLine(argc, argv);

    parser.addSection("General options");
	// TODO: fn_opt = parser.getOption("--opt", "optimiser STAR file from a previous 3D auto-refinement");

    starFn = parser.getOption("--i", "Input STAR file");	
    outPath = parser.getOption("--o", "Output directory, e.g. MotionFit/job041/");

    reference.read(parser, argc, argv);

    micrographHandler.firstFrame = textToInteger(parser.getOption("--first_frame", "", "1")) - 1;
    micrographHandler.lastFrame = textToInteger(parser.getOption("--last_frame", "", "-1")) - 1;
	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Skip those steps for which output files already exist.");
    verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

    motionEstimator.read(parser, argc, argv);
    motionParamEstimator.read(parser, argc, argv);
    frameRecombiner.read(parser, argc, argv);

    parser.addSection("Computational options");

    nr_omp_threads = textToInteger(parser.getOption("--j", "Number of (OMP) threads", "1"));
    minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
    maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));

    micrographHandler.saveMem = parser.checkOption("--sbs", "Load movies slice-by-slice to save memory (slower)");

    parser.addSection("Expert options");

	angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix) - read from STAR file by default", "-1"));
    debug = parser.checkOption("--debug", "Write debugging data");

    micrographHandler.movie_path = parser.getOption("--mov", "Path to movies", "");
    micrographHandler.corrMicFn = parser.getOption("--corr_mic", "List of uncorrected micrographs (e.g. corrected_micrographs.star)", "");
    micrographHandler.preextracted = parser.checkOption("--preex", "Preextracted movie stacks");
    micrographHandler.meta_path = parser.getOption("--meta", "Path to per-movie metadata star files", "");
    micrographHandler.gain_path = parser.getOption("--gain_path", "Path to gain references", "");
    micrographHandler.movie_ending = parser.getOption("--mov_end", "Ending of movie filenames", "");
    micrographHandler.movie_angpix = textToFloat(parser.getOption("--mps", "Pixel size of input movies (Angst/pix)", "-1"));
    micrographHandler.coords_angpix = textToFloat(parser.getOption("--cps", "Pixel size of particle coordinates in star-file (Angst/pix)", "-1"));
    micrographHandler.hotCutoff = textToFloat(parser.getOption("--hot", "Clip hot pixels to this max. value (-1 = off, TIFF only)", "-1"));
    micrographHandler.debug = parser.checkOption("--debug_mov", "Write debugging data for movie loading");

    movie_toReplace = parser.getOption("--mov_toReplace", "Replace this string in micrograph names...", "");
    movie_replaceBy = parser.getOption("--mov_replaceBy", "..by this one", "");

	// Check for errors in the command-line option
	if (parser.checkForErrors())
    {
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
    }

    if (argc == 1) parser.writeUsage(std::cerr);
}

void MotionRefiner::init()
{
    if (micrographHandler.movie_angpix <= 0 && micrographHandler.corrMicFn == "")
    {
        REPORT_ERROR("ERROR: Movie pixel size (--mps) is required unless a corrected_micrographs.star (--corr_mic) is provided.");
    }

    if (micrographHandler.coords_angpix <= 0 && micrographHandler.corrMicFn == "")
    {
        REPORT_ERROR("ERROR: Coordinates pixel size (--cps) is required unless a corrected_micrographs.star (--corr_mic) is provided.");
    }

	if (outPath[outPath.length()-1] != '/')
    {
        outPath += "/";
    }

    if (verb > 0)
    {
		std::cout << " + Reading " << starFn << "...\n";
    }

    mdt0.read(starFn);

    adaptMovieNames();

	if (angpix <= 0.0)
	{
		RFLOAT mag, dstep;
		mdt0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
		mdt0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
		angpix = 10000 * dstep / mag;

		if (verb > 0)
        {
            std::cout << "   - Using pixel size calculated from magnification and "
                      << "detector pixel size in the input STAR file: " << angpix << "\n";
        }
	}

    {
        double fractDose;
        std::string metaFn = "";

        // initialise corrected/uncorrected micrograph dictionary, then load the header
        // of the first movie (or read corrected_micrographs.star) to obtain:
        //  frame count, micrograph size and the fractional dose
        micrographHandler.init(
            mdt0, angpix, verb, nr_omp_threads, // in
            fc, fractDose, metaFn); // out

        // metaFn is only needed for console output
        // verb is needed since motionEstimator has not been initialised yet
        motionEstimator.proposeDosePerFrame(fractDose, metaFn, verb);
    }

    allMdts = StackHelper::splitByMicrographName(&mdt0);

    if (minMG >= allMdts.size())
    {
        std::stringstream sts0, sts1;
        sts0 << minMG;
        sts0 << allMdts.size();

        REPORT_ERROR("ERROR: Cannot start with micrograph "+sts0.str()
            +" (--min_MG); only "+sts1.str()+" micrographs defined in "+starFn+".");
    }

    if (minMG < 0)
    {
        minMG = 0;
    }

	// Only work on a user-specified subset of the micrographs
    if (maxMG < 0 || maxMG >= allMdts.size())
    {
        maxMG = allMdts.size()-1;
    }

    if (minMG > 0 || maxMG < allMdts.size()-1)
	{
		if (verb > 0)
        {
            std::cout << "   - Will only process micrographs in range: ["
                      << minMG << "-" << maxMG << "]"  << std::endl;
        }

        chosenMdts.clear();

        for (long int g = minMG; g <= maxMG; g++)
		{
            chosenMdts.push_back(allMdts[g]);
        }
	}
    else
    {
        chosenMdts = allMdts;
    }

    // make sure we don't try to load too many frames
    micrographHandler.findLowestFrameCount(chosenMdts, verb);

	if (only_do_unfinished)
    {
        motionMdts.clear();
        recombMdts.clear();

        motionMdts = MotionEstimator::findUnfinishedJobs(chosenMdts, outPath);
        recombMdts = FrameRecombiner::findUnfinishedJobs(chosenMdts, outPath);

		if (verb > 0)
        {
            if (motionMdts.size() > 0)
            {
                if (motionMdts.size() < chosenMdts.size())
                {
                    std::cout << "   - Will only estimate motion for " << motionMdts.size()
                          << " unfinished micrographs\n";
                }
                else
                {
                    std::cout << "   - Will estimate motion for all " << motionMdts.size()
                          << " micrographs - none are finished\n";
                }
            }
            else
            {
                std::cout << "   - Motion has already been estimated for all micrographs\n";
            }

            if (frameRecombiner.doingRecombination())
            {
                if (recombMdts.size() > 0)
                {
                    if (recombMdts.size() < chosenMdts.size())
                    {
                        std::cout << "   - Will only recombine frames for " << recombMdts.size()
                              << " unfinished micrographs\n";
                    }
                    else
                    {
                        std::cout << "   - Will recombine frames for all " << recombMdts.size()
                              << " micrographs - none are finished\n";
                    }
                }
                else
                {
                    std::cout << "   - Frames have already been recombined for all micrographs; "
                              << "a new STAR file will be generated\n";
                }
            }
		}
	}
    else
    {
        motionMdts = chosenMdts;
        recombMdts = chosenMdts;
    }

    estimateParams = motionParamEstimator.anythingToDo();
    estimateMotion = motionMdts.size() > 0;
    recombineFrames = frameRecombiner.doingRecombination() && (recombMdts.size() > 0);
    generateStar = frameRecombiner.doingRecombination();

    bool doAnything = estimateParams || estimateMotion || recombineFrames;
    bool needsReference = estimateParams || estimateMotion;

    if (doAnything)
	{
        obsModel = ObservationModel(angpix);

        // Read the first reference
        // (even if there is no motion to estimate - only to learn the image size)
        // @TODO: replace this once the data is tree-structured
        Image<RFLOAT> map0;
        map0.read(reference.reconFn0, false);

		// Get dimensions
        s = map0.data.xdim;
        sh = s/2 + 1;

        if (needsReference)
        {
            if (verb > 0)
            {
                std::cout << " + Reading references ...\n";
            }

            reference.load(verb);
        }

    }

    if (estimateMotion || estimateParams)
    {
        motionEstimator.init(verb, s, fc, nr_omp_threads, debug, outPath,
                             &reference, &obsModel, &micrographHandler);
    }

    if (estimateParams)
    {
        motionParamEstimator.init(
            verb, nr_omp_threads, debug,
            s, fc, chosenMdts, &motionEstimator, &reference, &obsModel);
    }
}

void MotionRefiner::run()
{
    if (estimateParams)
    {
        motionParamEstimator.run();

        return;
        // @TODO: apply the optimized parameters, then continue with motion estimation
    }

    // The subsets will be used in openMPI parallelisation: instead of over g0->gc,
    // they will be over smaller subsets
    if (estimateMotion)
    {
        motionEstimator.process(motionMdts, 0, motionMdts.size()-1);
    }

    if (recombineFrames)
    {
        frameRecombiner.init(
            allMdts, verb, s, fc, nr_omp_threads, outPath, debug,
            &obsModel, &micrographHandler);

        frameRecombiner.process(recombMdts, 0, recombMdts.size()-1);
	}

    if (generateStar)
    {
        combineEPSAndSTARfiles();
    }
}

int MotionRefiner::getVerbosityLevel()
{
    return verb;
}

/*
void MotionRefiner::calculateSingleFrameReconstruction(int iframe)
{

	FileName fn_vol1, fn_vol2;
	fn_vol1.compose(outPath + "frame", iframe, "", 3);
	fn_vol1 += "_half1_unfil.mrc";
	fn_vol2.compose(outPath + "frame", iframe, "", 3);
	fn_vol2 += "_half2_unfil.mrc";
	if (only_do_unfinished && exists(fn_vol1) && exists(fn_vol2))
	{
		if (verb > 0)
			std::cout << std::endl << " + " << fn_vol1 << " and" << fn_vol2 << " already exist: skipping per-frame reconstruction." << std::endl;
		return;
	}

	// Set current_size
	int current_size;
	if (perframe_highres > 0.)
		current_size = 2 * ROUND(ori_size * angpix / perframe_highres);
	else
		current_size = ori_size;

	BackProjector BP1(ori_size, 3, fn_sym, TRILINEAR, 1); // use only pad1 to save memory...
	BP1.initZeros(current_size);
	BackProjector BP2(ori_size, 3, fn_sym, TRILINEAR, 1);
	BP2.initZeros(current_size);

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    // Loop over all individual micrographs
	long nr_done = 0;
	long int pctot;
	for (long g = 0; g <= mdts.size(); g++)
	{

        const int pc = mdts[g].numberOfObjects();
        if (pc == 0) continue;

        pctot += pc;

        std::vector<std::vector<Image<Complex>>> movie;
        movie = loadMovie(g, pc, fts, iframe);

		FileName fn_root = getOutputFileNameRoot(g);
        std::vector<std::vector<d2Vector>> shift;
        shift = MotionHelper::readTracks(fn_root+"_tracks.star");


        // Loop over all particles in this micrograph
        for (int p = 0; p < pc; p++)
        {
    		CTF ctf;
    		std::string dum;
    		Matrix2D<RFLOAT> A3D;
    		MultidimArray<Complex > Faux, F2D, F2Dp, Fsub;
    		MultidimArray<RFLOAT> Fweight, Fctf;

    		RFLOAT xtrans, ytrans;
    		RFLOAT rot, tilt, psi;
    		int i_half;

    		windowFourierTransform(movie[p][0](), F2D, current_size);


			mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, i_half, p);
			mdts[g].getValue(EMDL_ORIENT_ROT_PRIOR, rot, p);
			mdts[g].getValue(EMDL_ORIENT_TILT_PRIOR, tilt, p);
			mdts[g].getValue(EMDL_ORIENT_PSI_PRIOR, psi, p);
			Euler_angles2matrix(rot, tilt, psi, A3D);

			// Translations
			xtrans=ytrans=0.;
			mdts[g].getValue(EMDL_ORIENT_ORIGIN_X, xtrans, p);
			mdts[g].getValue(EMDL_ORIENT_ORIGIN_Y, ytrans, p);

			// And fitted tracks
			xtrans -= shift[p][iframe].x;
			ytrans -= shift[p][iframe].y;

			shiftImageInFourierTransform(F2D, F2D, s, xtrans, ytrans);


			// CTF
			Fctf.resize(F2D);
			Fctf.initConstant(1.);
			if (do_ctf)
			{
				ctf.read(mdts[g], mdts[g], p);
				ctf.getFftwImage(Fctf, s, s, angpix, false, false, false, true);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
				{
					DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
				}
			}

			if (i_half == 1)
				BP1.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
			else if (i_half == 2)
				BP2.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
			else
				REPORT_ERROR("ERROR: invalid rlnRandomSubset value");
        }
   	}

	BP1.symmetrise(helical_nr_asu, helical_twist, helical_rise / angpix);
	BP2.symmetrise(helical_nr_asu, helical_twist, helical_rise / angpix);

	// Now do the reconstructions
	MultidimArray<RFLOAT> dummy;
	Image<RFLOAT> vol;
	BP1.reconstruct(vol(), 10, false, 1., dummy, dummy, dummy, dummy, dummy);
	vol.write(fn_vol1);
	BP2.reconstruct(vol(), 10, false, 1., dummy, dummy, dummy, dummy, dummy);
	vol.write(fn_vol2);

}
*/

// combine all EPS files into one logfile.pdf
void MotionRefiner::combineEPSAndSTARfiles()
{
    std::vector<FileName> fn_eps;

    if (verb > 0)
    {
    	std::cout << " + Combining all EPS and STAR files ... " << std::endl;
    }

    MetaDataTable mdtAll;

    if (frameRecombiner.doingRecombination())
    {
    	if (exists(outPath+"bfactors.eps"))
        {
    		fn_eps.push_back(outPath+"bfactors.eps");
        }

    	if (exists(outPath+"scalefactors.eps"))
        {
    		fn_eps.push_back(outPath+"scalefactors.eps");
        }
    }

    for (long g = 0; g < allMdts.size(); g++)
	{
        FileName fn_root = getOutputFileNameRoot(outPath, allMdts[g]);

		if (exists(fn_root+"_tracks.eps"))
        {
			fn_eps.push_back(fn_root+"_tracks.eps");
        }

        if (frameRecombiner.doingRecombination() && exists(fn_root+"_shiny.star"))
		{
			MetaDataTable mdt;
            mdt.read(fn_root+"_shiny.star");
            mdtAll.append(mdt);
		}
	}

	if (fn_eps.size() > 0)
	{
		joinMultipleEPSIntoSinglePDF(outPath + "logfile.pdf ", fn_eps);
	}

    if (frameRecombiner.doingRecombination())
    {
		mdtAll.write(outPath+"shiny.star");
    }

	if (verb > 0)
	{
		std::cout << " + Done! " << std::endl;
		std::cout << " + Written logfile in " << outPath << "logfile.pdf" << std::endl;

        if (frameRecombiner.doingRecombination())
        {
            std::cout << " + Written new particle STAR file in "
                      << outPath << "shiny.star" << std::endl;
        }
	}
}

// Get the output filename from the micrograph filename
FileName MotionRefiner::getOutputFileNameRoot(std::string outPath, const MetaDataTable& mdt)
{
    FileName fn_mic, fn_pre, fn_jobnr, fn_post;
    mdt.getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);
    decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
	return outPath + fn_post.withoutExtension();
}

void MotionRefiner::adaptMovieNames()
{
    if (movie_toReplace != "")
    {
        std::string name;

        for (int i = 0; i < mdt0.numberOfObjects(); i++)
        {
            mdt0.getValue(EMDL_MICROGRAPH_NAME, name, i);

            if (i == 0 && verb > 0)
            {
                std::cout << "   - Replacing: " << name << "\n";
            }

            std::string::size_type pos0 = name.find(movie_toReplace);

            if (pos0 != std::string::npos)
            {
                std::string::size_type pos1 = pos0 + movie_toReplace.length();

                std::string before = name.substr(0, pos0);
                std::string after = pos1 < name.length()? name.substr(pos1) : "";

                name = before + movie_replaceBy + after;
            }

            if (i == 0 && verb > 0)
            {
                std::cout << "                -> " << name << "\n";
            }

            if (micrographHandler.movie_path != "")
            {
                name = micrographHandler.movie_path + name.substr(name.find_last_of('/')+1);

                if (i == 0 && verb > 0)
                {
                    std::cout << "                -> " << name << "\n";
                }
            }

            if (micrographHandler.movie_ending != "")
            {
                name = name.substr(0, name.find_last_of('.')) + micrographHandler.movie_ending;

                if (i == 0 && verb > 0)
                {
                    std::cout << "                -> " << name << "\n";
                }
            }

            mdt0.setValue(EMDL_MICROGRAPH_NAME, name, i);
        }
    }
}



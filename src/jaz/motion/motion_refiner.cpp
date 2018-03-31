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
:   motionParamEstimator(*this),
    motionEstimator(*this),
    frameRecombiner(*this)
{
}

void MotionRefiner::read(int argc, char **argv)
{
    IOParser parser;
    parser.setCommandLine(argc, argv);

    parser.addSection("General options");
	// TODO: fn_opt = parser.getOption("--opt", "optimiser STAR file from a previous 3D auto-refinement");

    starFn = parser.getOption("--i", "Input STAR file");
	reconFn0 = parser.getOption("--m1", "Reference map, half 1");
    reconFn1 = parser.getOption("--m2", "Reference map, half 2");
    maskFn = parser.getOption("--mask", "Reference mask", "");
    fscFn = parser.getOption("--f", "Input STAR file with the FSC of the reference");
    outPath = parser.getOption("--o", "Output directory, e.g. MotionFit/job041/");
    firstFrame = textToInteger(parser.getOption("--first_frame", "", "1")) - 1;
    lastFrame = textToInteger(parser.getOption("--last_frame", "", "-1")) - 1;
	only_do_unfinished = parser.checkOption("--only_do_unfinished", "Skip those steps for which output files already exist.");
    verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

    motionEstimator.read(parser, argc, argv);
    motionParamEstimator.read(parser, argc, argv);
    frameRecombiner.read(parser, argc, argv);

    parser.addSection("Computational options");

    nr_omp_threads = textToInteger(parser.getOption("--j", "Number of (OMP) threads", "1"));
    paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
    minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
    maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));

    saveMem = parser.checkOption("--sbs", "Load movies slice-by-slice to save memory (slower)");

    parser.addSection("Expert options");

	angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix) - read from STAR file by default", "-1"));
    movie_path = parser.getOption("--mov", "Path to movies", "");
    corrMicFn = parser.getOption("--corr_mic", "List of uncorrected micrographs (e.g. corrected_micrographs.star)", "");
    preextracted = parser.checkOption("--preex", "Preextracted movie stacks");
    meta_path = parser.getOption("--meta", "Path to per-movie metadata star files", "");
    gain_path = parser.getOption("--gain_path", "Path to gain references", "");
    movie_ending = parser.getOption("--mov_end", "Ending of movie filenames", "");
    movie_toReplace = parser.getOption("--mov_toReplace", "Replace this string in micrograph names...", "");
    movie_replaceBy = parser.getOption("--mov_replaceBy", "..by this one", "");
    movie_angpix = textToFloat(parser.getOption("--mps", "Pixel size of input movies (Angst/pix)", "-1"));
    coords_angpix = textToFloat(parser.getOption("--cps", "Pixel size of particle coordinates in star-file (Angst/pix)", "-1"));
    hotCutoff = textToFloat(parser.getOption("--hot", "Clip hot pixels to this max. value (-1 = off, TIFF only)", "-1"));
    debug = parser.checkOption("--debug", "Write debugging data");
    debugMov = parser.checkOption("--debug_mov", "Write debugging data for movie loading");

	// Check for errors in the command-line option
	if (parser.checkForErrors())
    {
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
    }

    if (argc == 1) parser.writeUsage(std::cerr);
}

void MotionRefiner::init()
{
    if (movie_angpix <= 0 && corrMicFn == "")
    {
        REPORT_ERROR("ERROR: Movie pixel size (--mps) is required unless a corrected_micrographs.star (--corr_mic) is provided.");
    }

    if (coords_angpix <= 0 && corrMicFn == "")
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

    // initialise corrected/uncorrected micrograph dictionary
	if (corrMicFn != "")
	{
		MetaDataTable corrMic;
		corrMic.read(corrMicFn);

		mic2meta.clear();

		std::string micName, metaName;

		for (int i = 0; i < corrMic.numberOfObjects(); i++)
		{
			corrMic.getValueToString(EMDL_MICROGRAPH_NAME, micName, i);
			corrMic.getValueToString(EMDL_MICROGRAPH_METADATA_NAME, metaName, i);

			// remove the pipeline job prefix
			FileName fn_pre, fn_jobnr, fn_post;
			decomposePipelineFileName(micName, fn_pre, fn_jobnr, fn_post);

			mic2meta[fn_post] = metaName;
		}

		hasCorrMic = true;
	}
    else
    {
        hasCorrMic = false;
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

    loadInitialMovie();

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

		for (long int g = minMG; g <= maxMG; g++ )
		{
            chosenMdts.push_back(allMdts[g]);
        }
	}
    else
    {
        chosenMdts = allMdts;
    }

    // Check which motion-fit and recombination output files already exist.
	if (only_do_unfinished)
    {
        motionMdts.clear();
        recombMdts.clear();

        for (long int g = 0; g < chosenMdts.size(); g++ )
		{
            std::string fnRoot = getOutputFileNameRoot(chosenMdts[g]);

            bool motionDone = MotionEstimator::isFinished(fnRoot);
            bool recombDone = FrameRecombiner::isFinished(fnRoot);

            if (!motionDone)
            {
                motionMdts.push_back(chosenMdts[g]);
            }

            if (frameRecombiner.doCombineFrames && (!motionDone || !recombDone))
            {
                recombMdts.push_back(chosenMdts[g]);
            }
        }

		if (verb > 0)
        {
            if (motionMdts.size() > 0)
            {
                std::cout << "   - Will only estimate motion for " << motionMdts.size()
                          << " unfinished micrographs\n";
            }
            else
            {
                std::cout << "   - Will not estimate motion for any unfinished micrographs\n";
            }

            if (recombMdts.size() > 0) // implies frameRecombiner.doCombineFrames
            {
                std::cout << "   - Will only recombine frames for " << recombMdts.size()
                          << " unfinished micrographs\n";
            }
            else if (frameRecombiner.doCombineFrames)
            {
                std::cout << "   - Will not combine frames for any unfinished micrographs, "
                          << "just generate a STAR file\n";
            }
		}
	}

    estimateParams = motionParamEstimator.estim2 ||  motionParamEstimator.estim3;
    estimateMotion = motionMdts.size() > 0;
    recombineFrames = frameRecombiner.doCombineFrames && (recombMdts.size() > 0);
    generateStar = frameRecombiner.doCombineFrames;


    bool doAnything = estimateParams || estimateMotion || recombineFrames;
    bool needsReference = estimateParams || estimateMotion;

    if (doAnything)
	{
		obsModel = ObservationModel(angpix);

		if (verb > 0)
        {
			std::cout << " + Reading references ...\n";
        }

        // Read in the first reference
        // (even if there is no motion to estimate - only to learn the image size)
        // TODO: replace once there is global information available
		maps[0].read(reconFn0);

		if (maps[0].data.xdim != maps[0].data.ydim || maps[0].data.ydim != maps[0].data.zdim)
        {
			REPORT_ERROR(reconFn0 + " is not cubical.\n");
        }

		// Get dimensions
		s = maps[0].data.xdim;
        sh = s/2 + 1;

        if (needsReference)
		{
			// Read in the second reference
			maps[1].read(reconFn1);

			if (maps[1].data.xdim != maps[1].data.ydim || maps[1].data.ydim != maps[1].data.zdim)
            {
				REPORT_ERROR(reconFn1 + " is not cubical.\n");
            }

			if (   maps[0].data.xdim != maps[1].data.xdim
				|| maps[0].data.ydim != maps[1].data.ydim
				|| maps[0].data.zdim != maps[1].data.zdim)
            {
				REPORT_ERROR(reconFn0 + " and " + reconFn1 + " are of unequal size.\n");
            }

			if (maskFn != "")
			{
                if (verb > 0) std::cout << " + Masking references ...\n";

				Image<RFLOAT> mask, maskedRef;

				mask.read(maskFn);

				ImageOp::multiply(mask, maps[0], maskedRef);
				maps[0] = maskedRef;

				ImageOp::multiply(mask, maps[1], maskedRef);
				maps[1] = maskedRef;
			}


            if (verb > 0) std::cout << " + Transforming references ...\n";

			projectors[0] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
			projectors[0].computeFourierTransformMap(maps[0].data, powSpec[0].data, maps[0].data.xdim);
			projectors[1] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
			projectors[1].computeFourierTransformMap(maps[1].data, powSpec[1].data, maps[1].data.xdim);

			if (fscFn != "")
			{
				MetaDataTable fscMdt;
				fscMdt.read(fscFn, "fsc");

				if (!fscMdt.containsLabel(EMDL_SPECTRAL_IDX))
                {
                    REPORT_ERROR(fscFn + " does not contain a value for "
                                 + EMDL::label2Str(EMDL_SPECTRAL_IDX));
                }

				if (!fscMdt.containsLabel(EMDL_POSTPROCESS_FSC_TRUE))
                {
                    REPORT_ERROR(fscFn + " does not contain a value for "
                                 + EMDL::label2Str(EMDL_POSTPROCESS_FSC_TRUE));
                }

				RefinementHelper::drawFSC(&fscMdt, freqWeight1D, freqWeight);
			}
			else
			{
				freqWeight1D = std::vector<double>(sh,1.0);
				freqWeight = Image<RFLOAT>(sh,s);
				freqWeight.data.initConstant(1.0);
            }

			int k_out = sh;

			for (int i = 1; i < sh; i++)
			{
				if (freqWeight1D[i] <= 0.0)
				{
					k_out = i;
					break;
				}
			}

			if (verb > 0)
            {
                std::cout << " + maximum frequency to consider: "
                    << (s * angpix)/(RFLOAT)k_out << " A (" << k_out << " px)\n";
            }
		}
    } // if (doAnything)

    if (estimateMotion || estimateParams)
    {
        motionEstimator.init();
    }

    if (estimateParams)
    {
        motionParamEstimator.init();
    }

    if (recombineFrames)
    {
        frameRecombiner.init();
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
        frameRecombiner.process(recombMdts, 0, recombMdts.size()-1);
	}

    if (generateStar)
    {
        combineEPSAndSTARfiles();
    }
}

double MotionRefiner::angToPix(double a)
{
    return s * angpix / a;
}

double MotionRefiner::pixToAng(double p)
{
    return s * angpix / p;
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

    if (frameRecombiner.doCombineFrames)
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
        FileName fn_root = getOutputFileNameRoot(allMdts[g]);

		if (exists(fn_root+"_tracks.eps"))
        {
			fn_eps.push_back(fn_root+"_tracks.eps");
        }

        if (frameRecombiner.doCombineFrames && exists(fn_root+"_shiny.star"))
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

    if (frameRecombiner.doCombineFrames)
    {
		mdtAll.write(outPath+"shiny.star");
    }

	if (verb > 0)
	{
		std::cout << " + Done! " << std::endl;
		std::cout << " + Written logfile in " << outPath << "logfile.pdf" << std::endl;

        if (frameRecombiner.doCombineFrames)
        {
			std::cout << " + Written new particle STAR file in " << outPath << "shiny.star" << std::endl;
        }
	}
}

// Get the coordinate filename from the micrograph filename
FileName MotionRefiner::getOutputFileNameRoot(const MetaDataTable& mdt)
{
    FileName fn_mic, fn_pre, fn_jobnr, fn_post;
    mdt.getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);
    decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
	return outPath + fn_post.withoutExtension();
}

std::vector<std::vector<Image<Complex>>> MotionRefiner::loadMovie(
    const MetaDataTable& mdt, std::vector<ParFourierTransformer>& fts)
{
    std::vector<std::vector<Image<Complex>>> movie;

    if (preextracted)
    {
        movie = StackHelper::loadMovieStackFS(
            &mdt, "", false, nr_omp_threads, &fts,
            firstFrame, lastFrame);
    }
    else
    {
        std::string mgFn;
        mdt.getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);
        FileName fn_pre, fn_jobnr, fn_post;
        decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

        if (hasCorrMic)
        {
            //@TODO: make safe:
            std::string metaFn = mic2meta[fn_post];

            if (meta_path != "")
            {
                metaFn = meta_path + "/" + metaFn.substr(metaFn.find_last_of("/")+1);
            }

            micrograph = Micrograph(metaFn);

            std::string mgFn = micrograph.getMovieFilename();
            std::string gainFn = micrograph.getGainFilename();

            if (movie_ending != "")
            {
                mgFn.substr(0, mgFn.find_last_of(".")+1) + movie_ending;
            }

            if (movie_path != "")
            {
                mgFn = movie_path + "/" + mgFn.substr(mgFn.find_last_of("/")+1);
            }

            bool mgHasGain = false;

            if (gainFn != "")
            {
                if (gain_path != "")
                {
                    gainFn = gain_path + "/" + gainFn.substr(gainFn.find_last_of("/")+1);
                }

                if (gainFn != last_gainFn)
                {
                    lastGainRef.read(gainFn);
                    last_gainFn = gainFn;
                }

                mgHasGain = true;
            }

            movie = StackHelper::extractMovieStackFS(
                &mdt, mgHasGain? &lastGainRef : 0,
                mgFn, angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, firstFrame, lastFrame,
                hotCutoff, debugMov, saveMem);
        }
        else
        {
            // gain no longer supported without a corrected_micrographs.star:
            /*movie = StackHelper::extractMovieStackFS(
                &mdts[g], meta_path, movie_path, movie_ending,
                angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, firstFrame, lastFrame, hotCutoff, debugMov);*/

            movie = StackHelper::extractMovieStackFS(
                &mdt, 0,
                mgFn, angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, firstFrame, lastFrame,
                hotCutoff, debugMov, saveMem);
        }

        const int pc = movie.size();

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            StackHelper::varianceNormalize(movie[p], false);
        }
    }

    // Really warn each time?
    if (angpix < coords_angpix)
    {
        std::cerr << "WARING: pixel size (--angpix) is greater than the AutoPick pixel size (--coords_angpix)\n";

        if (coords_angpix < angpix + 0.01)
        {
            std::cerr << "        This is probably a rounding error. It is recommended to set --angpix ("
                      << angpix << ") to at least " << coords_angpix << "\n";

        }
    }

    if (angpix < movie_angpix)
    {
        std::cerr << "WARING: pixel size (--angpix) is greater than the movie pixel size (--movie_angpix)\n";

        if (movie_angpix < angpix + 0.01)
        {
            std::cerr << "        This is probably a rounding error. It is recommended to set --angpix ("
                      << angpix << ") to at least " << movie_angpix << "\n";

        }
    }

    return movie;
}

std::vector<std::vector<Image<Complex>>> MotionRefiner::loadInitialMovie()
{
    if (preextracted)
    {
        if (lastFrame < 0)
        {
             std::string name;
             allMdts[0].getValue(EMDL_MICROGRAPH_NAME, name, 0);

             Image<RFLOAT> stack0;
             stack0.read(name, false);

             const int pc0 = allMdts[0].numberOfObjects();
             const bool zstack = stack0.data.zdim > 1;
             const int stackSize = zstack? stack0.data.zdim : stack0.data.ndim;
             fc = stackSize / pc0 - firstFrame;
        }
        else
        {
            fc = lastFrame - firstFrame + 1;
        }

        FileName fn_mic;
        mdt0.getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);

        Image<RFLOAT> dum;
        dum.read(fn_mic, false);
        micrograph_xsize = XSIZE(dum());
        micrograph_ysize = YSIZE(dum());
    }
    else
    {
        if (hasCorrMic)
        {
            std::string mgFn;
            allMdts[0].getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

            // remove the pipeline job prefix
            FileName fn_pre, fn_jobnr, fn_post;
            decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

            //@TODO: make safe:
            std::string metaFn = mic2meta[fn_post];

            if (meta_path != "")
            {
                metaFn = meta_path + "/" + metaFn.substr(metaFn.find_last_of("/")+1);
            }

            micrograph = Micrograph(metaFn);

            if (movie_angpix <= 0)
            {
                movie_angpix = micrograph.angpix;

                if (verb > 0)
                {
                    std::cout << " + Using movie pixel size from " << metaFn << ": " << movie_angpix << " A\n";
                }
            }
            else
            {
                if (verb > 0)
                {
                    std::cout << " + Using movie pixel size from command line: " << movie_angpix << " A\n";
                }
            }

            if (coords_angpix <= 0)
            {
                coords_angpix = micrograph.angpix * micrograph.getBinningFactor();

                if (verb > 0)
                {
                    std::cout << " + Using coord. pixel size from " << metaFn << ": " << coords_angpix << " A\n";
                }
            }
            else
            {
                if (verb > 0)
                {
                    std::cout << " + Using coord. pixel size from command line: " << coords_angpix << " A\n";
                }
            }

            // this is safe to do - motionEstimator.read() has been called already
            if (motionEstimator.dosePerFrame < 0)
            {
                motionEstimator.dosePerFrame = micrograph.dose_per_frame;

                if (verb > 0)
                {
                    std::cout << " + Using dose per frame from " << metaFn << ": "
                              << motionEstimator.dosePerFrame << " A\n";
                }
            }

            micrograph_xsize = micrograph.getWidth();
            micrograph_ysize = micrograph.getHeight();

            if (lastFrame < 0)
            {
                fc = micrograph.getNframes() - firstFrame;
            }
            else
            {
                fc = lastFrame - firstFrame + 1;
            }
        }
        else
        {
            if (lastFrame < 0)
            {
                FileName mgFn;
                mdt0.getValue(EMDL_MICROGRAPH_NAME, mgFn, 0);

                Image<RFLOAT> dum;
                dum.read(mgFn, false);
                micrograph_xsize = XSIZE(dum());
                micrograph_ysize = YSIZE(dum());

                fc = (dum().zdim > 1? dum().zdim : dum().ndim) - firstFrame;
            }
            else
            {
                fc = lastFrame - firstFrame + 1;
            }
        }
    }
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

            if (movie_path != "")
            {
                name = movie_path + name.substr(name.find_last_of('/')+1);

                if (i == 0 && verb > 0)
                {
                    std::cout << "                -> " << name << "\n";
                }
            }

            if (movie_ending != "")
            {
                name = name + name.substr(0, name.find_last_of('.'));

                if (i == 0 && verb > 0)
                {
                    std::cout << "                -> " << name << "\n";
                }
            }

            mdt0.setValue(EMDL_MICROGRAPH_NAME, name, i);
        }
    }
}



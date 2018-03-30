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

#include "src/motion_fitter.h"

void MotionFitter::read(int argc, char **argv)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
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

	int fit_section = parser.addSection("Motion fit options");
    dmga = textToFloat(parser.getOption("--dmg_a", "Damage model, parameter a", " 3.40"));
    dmgb = textToFloat(parser.getOption("--dmg_b", "                        b", "-1.06"));
    dmgc = textToFloat(parser.getOption("--dmg_c", "                        c", "-0.54"));
    dosePerFrame = textToFloat(parser.getOption("--fdose", "Electron dose per frame (in e^-/A^2)", "-1"));
    sig_vel = textToFloat(parser.getOption("--s_vel", "Velocity sigma [Angst/dose]", "1.6"));
    sig_div = textToFloat(parser.getOption("--s_div", "Divergence sigma [Angst]", "500.0"));
    sig_acc = textToFloat(parser.getOption("--s_acc", "Acceleration sigma [Angst/dose]", "-1.0"));
    global_init = parser.checkOption("--gi", "Initialize with global trajectories instead of loading them from metadata file");
    expKer = parser.checkOption("--exp_k", "Use exponential kernel instead of sq. exponential");
    // TODO: provide k_cut in Angstroms!
    k_cutoff = textToFloat(parser.getOption("--k_cut", "Freq. cutoff (in pixels) for parameter estimation", "-1.0"));
    maxIters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "10000"));
    maxStep = textToFloat(parser.getOption("--max_step", "Maximum step size", "0.05"));
    unregGlob = parser.checkOption("--unreg_glob", "Do not regularize global component of motion");
    noGlobOff = parser.checkOption("--no_glob_off", "Do not compute initial per-particle offsets");
    debugOpt = parser.checkOption("--debug_opt", "Write optimization debugging info");
    diag = parser.checkOption("--diag", "Write out diagnostic data");

    /*
    parser.addSection("Parameter estimation");
    paramEstim2 = parser.checkOption("--params2", "Estimate 2 parameters instead of motion");
    paramEstim3 = parser.checkOption("--params3", "Estimate 3 parameters instead of motion");
    param_rV = textToFloat(parser.getOption("--r_vel", "Test s_vel +/- r_vel * s_vel", "0.5"));
    param_rD = textToFloat(parser.getOption("--r_div", "Test s_div +/- r_div * s_div", "0.5"));
    param_rA = textToFloat(parser.getOption("--r_acc", "Test s_acc +/- r_acc * s_acc", "0.5"));
    paramEstimIters = textToInteger(parser.getOption("--par_iters", "Parameter estimation is iterated this many times, each time halving the search range", "3"));
    paramEstimSteps = textToInteger(parser.getOption("--par_steps", "Parameter estimation takes max. this many steps before halving the range", "10"));
	*/

	int combine_section = parser.addSection("Combine frames options");
	doCombineFrames = parser.checkOption("--combine_frames", "Combine movie frames into polished particles.");
	k0a = textToFloat(parser.getOption("--bfac_minfreq", "Min. frequency used in B-factor fit [Angst]", "20"));
    k1a = textToFloat(parser.getOption("--bfac_maxfreq", "Max. frequency used in B-factor fit [Angst]", "-1"));
    bfacFn = parser.getOption("--bfactors", "A .star file with external B/k-factors", "");
    bfac_debug = parser.checkOption("--debug_bfactor", "Write out B/k-factor diagnostic data");

	int comp_section = parser.addSection("Computational options");
    nr_omp_threads = textToInteger(parser.getOption("--j", "Number of (OMP) threads", "1"));
    paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
    minMG = textToInteger(parser.getOption("--min_MG", "First micrograph index", "0"));
    maxMG = textToInteger(parser.getOption("--max_MG", "Last micrograph index (default is to process all)", "-1"));
    imgPath = parser.getOption("--img", "Path to images - read from STAR file by default", "");
    debug = parser.checkOption("--debug", "Write debugging data");
    debugMov = parser.checkOption("--debug_mov", "Write debugging data for movie loading");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));
    saveMem = parser.checkOption("--sbs", "Load movies slice-by-slice to save memory (slower)");

	int expert_section = parser.addSection("Expert options");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (angst/pix) - read from STAR file by default", "-1"));
	imgPath = parser.getOption("--mov", "Path to movies", "");
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

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	// Make sure outPath ends with a slash
	if (outPath[outPath.length()-1] != '/')
		outPath += "/";

}

void MotionFitter::usage()
{
	parser.writeUsage(std::cout);
}

void MotionFitter::initialise()
{

    if ((paramEstim2 || paramEstim3) && k_cutoff < 0)
        REPORT_ERROR("ERROR: Parameter estimation requires a freq. cutoff (--k_cut).");

    if (!global_init && corrMicFn == "")
    {
    	if (verb > 0)
    		std::cerr << "\nWarning: in the absence of a corrected_micrographs.star file (--corr_mic), global paths are used for initialization.\n";
        global_init = true;
    }

    if (paramEstim2 && paramEstim3)
    	 REPORT_ERROR("ERROR: Only 2 or 3 parameters can be estimated (--params2 or --params3), not both.");

    if (movie_angpix <= 0 && corrMicFn == "")
        REPORT_ERROR("ERROR: Movie pixel size (--mps) is required unless a corrected_micrographs.star (--corr_mic) is provided.");

    if (coords_angpix <= 0 && corrMicFn == "")
        REPORT_ERROR("ERROR: Coordinates pixel size (--cps) is required unless a corrected_micrographs.star (--corr_mic) is provided.");

	// Make sure output directory ends in a '/'
	if (outPath[outPath.length()-1] != '/')
		outPath+="/";

    if (verb > 0)
		std::cout << " + Reading " << starFn << "...\n";
	mdt0.read(starFn);

	// Get micrograph_xsize and micrograph_ysize for EPS plots
    FileName fn_mic;
    mdt0.getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);
    Image<RFLOAT> dum;
    dum.read(fn_mic, false);
    micrograph_xsize = XSIZE(dum());
    micrograph_ysize = YSIZE(dum());

    if (movie_toReplace != "")
    {
        std::string name;
        for (int i = 0; i < mdt0.numberOfObjects(); i++)
        {
            mdt0.getValue(EMDL_MICROGRAPH_NAME, name, i);

            if (i == 0 && verb > 0)
            	std::cout << "   - Replacing: " << name << " -> ";

            std::string::size_type pos0 = name.find(movie_toReplace);

            if (pos0 != std::string::npos)
            {
                std::string::size_type pos1 = pos0 + movie_toReplace.length();

                std::string before = name.substr(0, pos0);
                std::string after = pos1 < name.length()? name.substr(pos1) : "";

                name = before + movie_replaceBy + after;
            }

            if (i == 0 && verb > 0)
            	std::cout << name << "\n";

            mdt0.setValue(EMDL_MICROGRAPH_NAME, name, i);
        }
    }

	if (angpix <= 0.0)
	{
		RFLOAT mag, dstep;
		mdt0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
		mdt0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
		angpix = 10000 * dstep / mag;

		if (verb > 0)
			std::cout << "   - Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << "\n";
	}

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

	mdts = StackHelper::splitByMicrographName(&mdt0);

	// Jasenko's LoadInitialMovieValues
	if (preextracted)
	{

		if (lastFrame < 0)
		{
			std::string name, fullName, movieName;
			 mdts[0].getValue(EMDL_IMAGE_NAME, fullName, 0);
			 mdts[0].getValue(EMDL_MICROGRAPH_NAME, movieName, 0);
			 name = fullName.substr(fullName.find("@")+1);

			 std::string finName;

			 if (imgPath == "")
			 {
				 finName = name;
			 }
			 else
			 {
				 finName = imgPath + "/" + movieName.substr(movieName.find_last_of("/")+1);
			 }

			 Image<RFLOAT> stack0;
			 stack0.read(finName, false);

			 const int pc0 = mdts[0].numberOfObjects();
			 const bool zstack = stack0.data.zdim > 1;
			 const int stackSize = zstack? stack0.data.zdim : stack0.data.ndim;
			 fc = stackSize / pc0 - firstFrame;
		}
		else fc = lastFrame - firstFrame + 1;
	}
	else
	{
		if (hasCorrMic)
		{
			std::string mgFn;
			mdts[0].getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);

			// remove the pipeline job prefix
			FileName fn_pre, fn_jobnr, fn_post;
			decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);
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
					std::cout << " + Using movie pixel size from " << metaFn << ": " << movie_angpix << " A\n";
			}
			else
			{
				if (verb > 0)
					std::cout << " + Using movie pixel size from command line: " << movie_angpix << " A\n";
			}

			if (coords_angpix <= 0)
			{
				coords_angpix = micrograph.angpix * micrograph.getBinningFactor();
				if (verb > 0)
					std::cout << " + Using coord. pixel size from " << metaFn << ": " << coords_angpix << " A\n";
			}
			else
			{
				if (verb > 0)
					std::cout << " + Using coord. pixel size from command line: " << coords_angpix << " A\n";
			}

			if (dosePerFrame <= 0)
			{
				dosePerFrame = micrograph.dose_per_frame;
				if (verb > 0)
					std::cout << " + Using dose per frame from " << metaFn << ": " << dosePerFrame << " A\n";

			}

			if (lastFrame < 0) fc = micrograph.getNframes() - firstFrame;
			else fc = lastFrame - firstFrame + 1;

		}
		else
		{
			if (lastFrame < 0)
			{
				std::vector<std::vector<Image<Complex>>> movie = StackHelper::extractMovieStackFS(
						&mdts[0], meta_path, imgPath, movie_ending, coords_angpix, angpix, movie_angpix, s,
						nr_omp_threads, false, firstFrame, lastFrame, hotCutoff, debugMov);

				fc = movie[0].size() - firstFrame;
			}
			else fc = lastFrame - firstFrame + 1;
		}
	}

	// Only work on a user-specified subset of the micrographs
	if (maxMG < 0 || maxMG >= mdts.size())
		maxMG = mdts.size()-1;
	if (minMG < 0 || minMG >= mdts.size())
		minMG = 0;
	if (minMG > 0 || maxMG < mdts.size()-1)
	{
		if (verb > 0)
			std::cout << "   - Will only process micrographs in range: [" << minMG << "-" << maxMG << "]"  << std::endl;

		std::vector<MetaDataTable> todo_mdts;
		for (long int g = minMG; g <= maxMG; g++ )
		{
			todo_mdts.push_back(mdts[g]);
		}
		mdts = todo_mdts;
	}

	// check whether motion_fit output files exist and if they do, then skip this micrograph
	if (only_do_unfinished)
	{
		std::vector<MetaDataTable> unfinished_mdts;
		for (long int g = minMG; g <= maxMG; g++ )
		{
			bool is_done = true;
			if (!exists(getOutputFileNameRoot(g)+"_tracks.star"))
				is_done = false;
			if (!exists(getOutputFileNameRoot(g)+"_FCC_w1.mrc"))
				is_done = false;
			if (!is_done)
				unfinished_mdts.push_back(mdts[g]);
		}
		mdts = unfinished_mdts;
		if (verb > 0)
		{
			if (mdts.size() > 0)
				std::cout << "   - Will only estimate motion for " << mdts.size() << " unfinished micrographs" << std::endl;
			else
				std::cout << "   - Will not estimate motion for any unfinished micrographs" << std::endl;
		}
	}

	// Check whether there is something to do at all...
	if (mdts.size() > 0 || doCombineFrames)
	{

		obsModel = ObservationModel(angpix);

		// Only create projectors if there are (still) micrographs to process
		// Read in the first reference
		if (verb > 0)
			std::cout << " + Reading references ...\n";
		maps[0].read(reconFn0);

		if (maps[0].data.xdim != maps[0].data.ydim || maps[0].data.ydim != maps[0].data.zdim)
			REPORT_ERROR(reconFn0 + " is not cubical.\n");

		// Get dimensions
		s = maps[0].data.xdim;
		sh = s/2 + 1;

		// Only read and initialise the rest for unfinisged motion_fit
		if (mdts.size() > 0)
		{
			// Read in the second reference
			maps[1].read(reconFn1);

			if (maps[1].data.xdim != maps[1].data.ydim || maps[1].data.ydim != maps[1].data.zdim)
				REPORT_ERROR(reconFn1 + " is not cubical.\n");

			if (   maps[0].data.xdim != maps[1].data.xdim
				|| maps[0].data.ydim != maps[1].data.ydim
				|| maps[0].data.zdim != maps[1].data.zdim)
				REPORT_ERROR(reconFn0 + " and " + reconFn1 + " are of unequal size.\n");

			if (maskFn != "")
			{
				if (verb > 0)
					std::cout << " + Masking references ...\n";

				Image<RFLOAT> mask, maskedRef;

				mask.read(maskFn);

				ImageOp::multiply(mask, maps[0], maskedRef);
				maps[0] = maskedRef;

				ImageOp::multiply(mask, maps[1], maskedRef);
				maps[1] = maskedRef;
			}


			if (verb > 0)
				std::cout << " + Transforming references ...\n";


			projectors[0] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
			projectors[0].computeFourierTransformMap(maps[0].data, powSpec[0].data, maps[0].data.xdim);
			projectors[1] = Projector(s, TRILINEAR, paddingFactor, 10, 2);
			projectors[1].computeFourierTransformMap(maps[1].data, powSpec[1].data, maps[1].data.xdim);

			if (fscFn != "")
			{
				MetaDataTable fscMdt;
				fscMdt.read(fscFn, "fsc");

				if (!fscMdt.containsLabel(EMDL_SPECTRAL_IDX))
					REPORT_ERROR(fscFn + " does not contain a value for " + EMDL::label2Str(EMDL_SPECTRAL_IDX));
				if (!fscMdt.containsLabel(EMDL_POSTPROCESS_FSC_TRUE))
					REPORT_ERROR(fscFn + " does not contain a value for " + EMDL::label2Str(EMDL_POSTPROCESS_FSC_TRUE));

				RefinementHelper::drawFSC(&fscMdt, freqWeight1D, freqWeight);
			}
			else
			{
				freqWeight1D = std::vector<double>(sh,1.0);
				freqWeight = Image<RFLOAT>(sh,s);
				freqWeight.data.initConstant(1.0);
			}

			dmgWeight = DamageHelper::damageWeights(s, angpix, firstFrame, fc, dosePerFrame, dmga, dmgb, dmgc);
			int k_out = sh;

			for (int i = 1; i < sh; i++)
			{
				if (freqWeight1D[i] <= 0.0)
				{
					k_out = i;
					break;
				}
			}

			// TODO: output k_out in Angstroms
			if (verb > 0)
				std::cout << " + maximum frequency to be consider: = " << (s * angpix)/(RFLOAT)k_out << " Angstrom" << std::endl;

			for (int f = 0; f < fc; f++)
			{
				dmgWeight[f].data.xinit = 0;
				dmgWeight[f].data.yinit = 0;

				if (k_cutoff > 0.0)
				{
					std::stringstream stsf;
					stsf << f;
					dmgWeight[f] = FilterHelper::ButterworthEnvFreq2D(dmgWeight[f], k_cutoff-1, k_cutoff+1);

					ImageOp::multiplyBy(dmgWeight[f], freqWeight);
				}
			}
		}

	}

}

void MotionFitter::initialiseCombineFrames()
{

	// Split again, as a subset might have been done before for only_do_unfinished...
    mdts.clear();
    mdts = StackHelper::splitByMicrographName(&mdt0);

    // check whether combine_frames output stack exist and if they do, then skip this micrograph
	if (only_do_unfinished)
	{
		std::vector<MetaDataTable> unfinished_mdts;
		for (long int g = minMG; g <= maxMG; g++ )
		{
			bool is_done = true;
			if (!exists(getOutputFileNameRoot(g)+"_shiny.mrcs"))
				is_done = false;
			if (!exists(getOutputFileNameRoot(g)+"_shiny.star"))
				is_done = false;
			if (!is_done)
				unfinished_mdts.push_back(mdts[g]);
		}
		mdts = unfinished_mdts;
		if (verb > 0)
		{
			if (mdts.size() > 0)
				std::cout << "   - Will only combine frames for " << mdts.size() << " unfinished micrographs" << std::endl;
			else
				std::cout << "   - Will not combine frames for any unfinished micrographs, just generate a STAR file" << std::endl;
		}
	}


}

void MotionFitter::run()
{

	// The subsets will be used in openMPI parallelisation: instead of over g0->gc, they will be over smaller subsets
	if (mdts.size() > 0)
		processSubsetMicrographs(0, mdts.size()-1);

	if (doCombineFrames)
	{
		initialiseCombineFrames();
		if (mdts.size() > 0)
			combineFramesSubsetMicrographs(0, mdts.size()-1);
	}

	combineEPSAndSTARfiles();
}

/// Helper functions

void MotionFitter::processSubsetMicrographs(long g_start, long g_end)
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
    std::vector<Image<RFLOAT>>
            tables(nr_omp_threads),
            weights0(nr_omp_threads),
            weights1(nr_omp_threads);

    for (int i = 0; i < nr_omp_threads; i++)
    {
        FscHelper::initFscTable(sh, fc, tables[i], weights0[i], weights1[i]);
    }

    const double sig_vel_nrm = dosePerFrame * sig_vel / angpix;
    const double sig_acc_nrm = dosePerFrame * sig_acc / angpix;
    const double sig_div_nrm = dosePerFrame * sig_div / coords_angpix;

    int pctot = 0;

	long nr_done = 0;
	FileName prevdir = "";
	for (long g = g_start; g <= g_end; g++)
	{

        const int pc = mdts[g].numberOfObjects();
        if (pc == 0) continue;

		// Make sure output directory exists
		FileName newdir = getOutputFileNameRoot(g);
		newdir = newdir.beforeLastOf("/");
		if (newdir != prevdir)
		{
			std::string command = " mkdir -p " + newdir;
			int res = system(command.c_str());
		}


        std::vector<std::vector<Image<Complex>>> movie;
        std::vector<std::vector<Image<RFLOAT>>> movieCC;
        std::vector<d2Vector> positions;
        std::vector<std::vector<d2Vector>> initialTracks;
        std::vector<d2Vector> globComp;

        prepMicrograph(g, fts, dmgWeight,
                movie, movieCC, positions, initialTracks, globComp);

        pctot += pc;

        std::vector<std::vector<gravis::d2Vector>> tracks = optimize(
                movieCC, initialTracks,
                sig_vel_nrm, sig_acc_nrm, sig_div_nrm,
                positions, globComp,
                maxStep, 1e-9, 1e-9, maxIters, 0.0);

        updateFCC(movie, tracks, mdts[g], tables, weights0, weights1);

        writeOutput(tracks, tables, weights0, weights1, positions, outPath, g, 30.0);

        for (int i = 0; i < nr_omp_threads; i++)
        {
            tables[i].data.initZeros();
            weights0[i].data.initZeros();
            weights1[i].data.initZeros();
        }

		nr_done++;
		if (verb > 0 && nr_done % barstep == 0)
			progress_bar(nr_done);

	}


    if (verb > 0)
		progress_bar(my_nr_micrographs);

}

void MotionFitter::combineFramesSubsetMicrographs(long g_start, long g_end)
{

    std::vector<Image<RFLOAT>> freqWeights;

    // Either calculate weights from FCC or from user-provided B-factors
    hasBfacs = bfacFn != "";
    if (!hasBfacs)
    {
        freqWeights = weightsFromFCC();
    }
    else
    {
        freqWeights = weightsFromBfacs();
    }

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
	FileName prevdir = "";
	for (long g = g_start; g <= g_end; g++)
	{

        const int pc = mdts[g].numberOfObjects();
        if (pc == 0) continue;

        pctot += pc;

        std::vector<std::vector<Image<Complex>>> movie;
        movie = loadMovie(g, pc, fts);

		FileName fn_root = getOutputFileNameRoot(g);
        std::vector<std::vector<d2Vector>> shift;
        shift = MotionRefinement::readTracks(fn_root+"_tracks.star");

        Image<RFLOAT> stack(s,s,1,pc);

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            int threadnum = omp_get_thread_num();

            Image<Complex> sum(sh,s);
            sum.data.initZeros();

            Image<Complex> obs(sh,s);

            for (int f = 0; f < fc; f++)
            {
            	shiftImageInFourierTransform(movie[p][f](), obs(), s, -shift[p][f].x, -shift[p][f].y);

                for (int y = 0; y < s; y++)
                for (int x = 0; x < sh; x++)
                {
                    sum(y,x) += freqWeights[f](y,x) * obs(y,x);
                }
            }

            Image<RFLOAT> real(s,s);

            fts[threadnum].inverseFourierTransform(sum(), real());

            for (int y = 0; y < s; y++)
            for (int x = 0; x < s; x++)
            {
                DIRECT_NZYX_ELEM(stack(), p, 0, y, x) = real(y,x);
            }
        }

        stack.write(fn_root+"_shiny.mrcs");

        if (debug)
        {
            VtkHelper::writeTomoVTK(stack, fn_root+"_shiny.vtk");
        }

        for (int p = 0; p < pc; p++)
        {
            std::stringstream sts;
            sts << (p+1);
            mdts[g].setValue(EMDL_IMAGE_NAME, sts.str() + "@" + fn_root+"_shiny.mrcs", p);
        }

        mdts[g].write(fn_root+"_shiny.star");

		nr_done++;
		if (verb > 0 && nr_done % barstep == 0)
			progress_bar(nr_done);

	}

    if (verb > 0)
		progress_bar(my_nr_micrographs);

}


/*
void MotionFitter::calculateSingleFrameReconstruction(int iframe)
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
        shift = MotionRefinement::readTracks(fn_root+"_tracks.star");


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
void MotionFitter::combineEPSAndSTARfiles()
{
    std::vector<FileName> fn_eps;
    if (verb > 0)
    	std::cout << " + Combining all EPS and STAR files ... " << std::endl;


    MetaDataTable mdtAll;

    if (doCombineFrames)
    {
    	if (exists(outPath+"bfactors.eps"))
    		fn_eps.push_back(outPath+"bfactors.eps");
    	if (exists(outPath+"scalefactors.eps"))
    		fn_eps.push_back(outPath+"scalefactors.eps");

    }

	// Split again, as a subset might have been done before for only_do_unfinished...
    mdts.clear();
    mdts = StackHelper::splitByMicrographName(&mdt0);

    for (long g = 0; g < mdts.size(); g++)
	{
		FileName fn_root = getOutputFileNameRoot(g);
		if (exists(fn_root+"_tracks.eps"))
			fn_eps.push_back(fn_root+"_tracks.eps");

		if (doCombineFrames)
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

	if (doCombineFrames)
		mdtAll.write(outPath+"shiny.star");

	if (verb > 0)
	{
		std::cout << " + Done! " << std::endl;
		std::cout << " + Written logfile in " << outPath << "logfile.pdf" << std::endl;
		if (doCombineFrames)
			std::cout << " + Written new particle STAR file in " << outPath << "shiny.star" << std::endl;
	}

}


// Get the coordinate filename from the micrograph filename
FileName MotionFitter::getOutputFileNameRoot(long int g)
{

    FileName fn_mic, fn_pre, fn_jobnr, fn_post;
	mdts[g].getValue(EMDL_MICROGRAPH_NAME, fn_mic, 0);
    decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
	return outPath + fn_post.withoutExtension();

}

std::string MotionFitter::getMicrographTag(long g)
{
    std::string tag;
    mdts[g].getValue(EMDL_IMAGE_NAME, tag, 0);
    tag = tag.substr(0,tag.find_last_of('.'));
    tag = tag.substr(tag.find_first_of('@')+1);

    return tag;
}

std::vector<std::vector<Image<Complex>>> MotionFitter::loadMovie(
        long g, int pc, std::vector<ParFourierTransformer>& fts, int only_this_frame)
{
    std::vector<std::vector<Image<Complex>>> movie;

    int myFirstFrame = (only_this_frame < 0) ? firstFrame : only_this_frame;
    int myLastFrame = (only_this_frame < 0) ? lastFrame : only_this_frame;

    if (preextracted)
    {
        movie = StackHelper::loadMovieStackFS(
            &mdts[g], imgPath, false, nr_omp_threads, &fts,
            myFirstFrame, myLastFrame);
    }
    else
    {
        std::string mgFn;
        mdts[g].getValueToString(EMDL_MICROGRAPH_NAME, mgFn, 0);
        FileName fn_pre, fn_jobnr, fn_post;
        decomposePipelineFileName(mgFn, fn_pre, fn_jobnr, fn_post);

        if (hasCorrMic)
        {
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

            if (imgPath != "")
            {
                mgFn = imgPath + "/" + mgFn.substr(mgFn.find_last_of("/")+1);
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
                &mdts[g], mgHasGain? &lastGainRef : 0,
                mgFn, angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, myFirstFrame, myLastFrame, hotCutoff, debugMov, saveMem);
        }
        else
        {
            movie = StackHelper::extractMovieStackFS(
                &mdts[g], meta_path, imgPath, movie_ending,
                angpix, coords_angpix, movie_angpix, s,
                nr_omp_threads, true, myFirstFrame, myLastFrame, hotCutoff, debugMov);
        }

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            StackHelper::varianceNormalize(movie[p], false);
        }
    }

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

void MotionFitter::prepMicrograph(
        long g, std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<RFLOAT>>& dmgWeight,
        std::vector<std::vector<Image<Complex>>>& movie,
        std::vector<std::vector<Image<RFLOAT>>>& movieCC,
        std::vector<d2Vector>& positions,
        std::vector<std::vector<d2Vector>>& initialTracks,
        std::vector<d2Vector>& globComp)
{
    const int pc = mdts[g].numberOfObjects();

    movie = loadMovie(g, pc, fts); // throws exceptions
    std::vector<double> sigma2 = StackHelper::powerSpectrum(movie);

    #pragma omp parallel for num_threads(nr_omp_threads)
    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        MotionRefinement::noiseNormalize(movie[p][f], sigma2, movie[p][f]);
    }

    positions = std::vector<gravis::d2Vector>(pc);

    for (int p = 0; p < pc; p++)
    {
        mdts[g].getValue(EMDL_IMAGE_COORD_X, positions[p].x, p);
        mdts[g].getValue(EMDL_IMAGE_COORD_Y, positions[p].y, p);
    }

    movieCC = MotionRefinement::movieCC(
            projectors[0], projectors[1], obsModel, mdts[g], movie,
            sigma2, dmgWeight, fts, nr_omp_threads);

    initialTracks = std::vector<std::vector<d2Vector>>(pc);

    if (global_init)
    {
        std::vector<Image<RFLOAT>> ccSum = MotionRefinement::addCCs(movieCC);
        std::vector<gravis::d2Vector> globTrack = MotionRefinement::getGlobalTrack(ccSum);
        std::vector<gravis::d2Vector> globOffsets;

        if (noGlobOff)
        {
            globOffsets = std::vector<d2Vector>(pc, d2Vector(0,0));
        }
        else
        {
            globOffsets = MotionRefinement::getGlobalOffsets(
                    movieCC, globTrack, 0.25*s, nr_omp_threads);
        }

        if (diag)
        {
            ImageLog::write(ccSum, getOutputFileNameRoot(g) + "_CCsum.mrc", CenterXY);
        }

        for (int p = 0; p < pc; p++)
        {
            initialTracks[p] = std::vector<d2Vector>(fc);

            for (int f = 0; f < fc; f++)
            {
                if (unregGlob)
                {
                    initialTracks[p][f] = globOffsets[p];
                }
                else
                {
                    initialTracks[p][f] = globTrack[f] + globOffsets[p];
                }
            }
        }

        globComp = unregGlob? globTrack : std::vector<d2Vector>(fc, d2Vector(0,0));
    }
    else
    {
        const d2Vector inputScale(
                coords_angpix / (movie_angpix * micrograph.getWidth()),
                coords_angpix / (movie_angpix * micrograph.getHeight()));

        const double outputScale = movie_angpix / angpix;

        globComp = std::vector<d2Vector>(fc, d2Vector(0,0));

        if (unregGlob)
        {
            for (int f = 0; f < fc; f++)
            {
                RFLOAT sx, sy;
                micrograph.getShiftAt(f+1, 0, 0, sx, sy, false);

                globComp[f] = -outputScale * d2Vector(sx, sy);
            }
        }

        for (int p = 0; p < pc; p++)
        {
            initialTracks[p] = std::vector<d2Vector>(fc);

            for (int f = 0; f < fc; f++)
            {
                d2Vector in(inputScale.x * positions[p].x - 0.5,
                            inputScale.y * positions[p].y - 0.5);

                RFLOAT sx, sy;

                micrograph.getShiftAt(f+1, in.x, in.y, sx, sy, true);

                initialTracks[p][f] = -outputScale * d2Vector(sx,sy) - globComp[f];
            }
        }
    }
}

std::vector<std::vector<d2Vector>> MotionFitter::optimize(
        const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
        const std::vector<std::vector<d2Vector>>& inTracks,
        double sig_vel_px, double sig_acc_px, double sig_div_px,
        const std::vector<d2Vector>& positions,
        const std::vector<d2Vector>& globComp,
        double step, double minStep, double minDiff,
        long maxIters, double inertia)
{
    const double eps = 1e-20;

    if (sig_vel_px < eps)
    {
        std::cerr << "Warning: sig_vel < " << eps << " px. Setting to " << eps << ".\n";
        sig_vel_px = eps;
    }

    if (sig_div_px < eps)
    {
        std::cerr << "Warning: sig_div < " << eps << " px. Setting to " << eps << ".\n";
        sig_div_px = eps;
    }

    const int pc = inTracks.size();

    if (pc == 0) return std::vector<std::vector<d2Vector>>(0);

    const int fc = inTracks[0].size();

    GpMotionFit gpmf(movieCC, sig_vel_px, sig_div_px, sig_acc_px,
                     pc, positions,
                     globComp, nr_omp_threads, expKer);


    std::vector<double> initialCoeffs(2*(pc + pc*(fc-1)));

    gpmf.posToParams(inTracks, initialCoeffs);

    std::vector<double> optPos = GradientDescent::optimize(
            initialCoeffs, gpmf, step, minStep, minDiff, maxIters, inertia, debugOpt);

    std::vector<std::vector<d2Vector>> out(pc, std::vector<d2Vector>(fc));
    gpmf.paramsToPos(optPos, out);

    return out;
}


void MotionFitter::updateFCC(
        const std::vector<std::vector<Image<Complex>>>& movie,
        const std::vector<std::vector<d2Vector>>& tracks,
        const MetaDataTable& mdt,
        std::vector<Image<RFLOAT>>& tables,
        std::vector<Image<RFLOAT>>& weights0,
        std::vector<Image<RFLOAT>>& weights1)
{
    const int pc = mdt.numberOfObjects();

    #pragma omp parallel for num_threads(nr_omp_threads)
    for (int p = 0; p < pc; p++)
    {
        int threadnum = omp_get_thread_num();

        Image<Complex> pred;
        std::vector<Image<Complex>> obs = movie[p];

        for (int f = 0; f < fc; f++)
        {
            shiftImageInFourierTransform(obs[f](), obs[f](), s, -tracks[p][f].x, -tracks[p][f].y);
        }

        int randSubset;
        mdt.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
        randSubset -= 1;

        if (randSubset == 0)
        {
            pred = obsModel.predictObservation(projectors[1], mdt, p, true, true);
        }
        else
        {
            pred = obsModel.predictObservation(projectors[0], mdt, p, true, true);
        }

        FscHelper::updateFscTable(obs, pred, tables[threadnum],
                                  weights0[threadnum], weights1[threadnum]);
    }
}

void MotionFitter::writeOutput(
        const std::vector<std::vector<d2Vector>>& tracks,
        const std::vector<Image<RFLOAT>>& fccData,
        const std::vector<Image<RFLOAT>>& fccWeight0,
        const std::vector<Image<RFLOAT>>& fccWeight1,
        const std::vector<d2Vector>& positions,
        std::string outPath, long g,
        double visScale)
{
    const int pc = tracks.size();

    if (pc == 0) return;

    const int fc = tracks[0].size();

    FileName fn_root = getOutputFileNameRoot(g);
    MotionRefinement::writeTracks(tracks, fn_root + "_tracks.star");

    Image<RFLOAT> fccDataSum(sh,fc), fccWeight0Sum(sh,fc), fccWeight1Sum(sh,fc);
    fccDataSum.data.initZeros();
    fccWeight0Sum.data.initZeros();
    fccWeight1Sum.data.initZeros();

    for (int i = 0; i < fccData.size(); i++)
    {
        for (int y = 0; y < fc; y++)
        for (int x = 0; x < sh; x++)
        {
            fccDataSum(y,x) += fccData[i](y,x);
            fccWeight0Sum(y,x) += fccWeight0[i](y,x);
            fccWeight1Sum(y,x) += fccWeight1[i](y,x);
        }
    }

    fccDataSum.write(fn_root + "_FCC_cc.mrc");
    fccWeight0Sum.write(fn_root + "_FCC_w0.mrc");
    fccWeight1Sum.write(fn_root + "_FCC_w1.mrc");


    // plot EPS graph with all observed and fitted tracks
    std::vector<std::vector<gravis::d2Vector>> visTracks(pc);

    for (int p = 0; p < pc; p++)
    {
    	visTracks[p] = std::vector<gravis::d2Vector>(fc);
    }

    std::vector<gravis::d2Vector> globalTrack(fc);

    for (int f = 0; f < fc; f++)
    {
        globalTrack[f] = d2Vector(0,0);

        for (int p = 0; p < pc; p++)
        {
            globalTrack[f] += tracks[p][f];
        }

        globalTrack[f] /= pc;
        for (int p = 0; p < pc; p++)
        {
            visTracks[p][f] = positions[p] + visScale * tracks[p][f];
        }
    }

	// Make a postscript with the tracks
	FileName fn_eps = fn_root + "_tracks.eps";
	CPlot2D *plot2D=new CPlot2D(fn_eps);
 	plot2D->SetXAxisSize(600);
 	plot2D->SetYAxisSize(600);
	plot2D->SetDrawLegend(false);

 	// Global track in the middle
	CDataSet dataSet;
	dataSet.SetDrawMarker(false);
	dataSet.SetDatasetColor(0.0,0.0,1.0);
	dataSet.SetLineWidth(1.);
	RFLOAT xshift, yshift;
	const RFLOAT xcenter =  micrograph_xsize / 2.0;
	const RFLOAT ycenter =  micrograph_ysize / 2.0;
	for (int f = 0; f < fc; f++)
    {
		CDataPoint point(xcenter + visScale * globalTrack[f].x, ycenter + visScale * globalTrack[f].y);
		dataSet.AddDataPoint(point);
	}
	plot2D->AddDataSet(dataSet);

	// Mark starting point global track
	CDataSet dataSetStart;
	dataSetStart.SetDrawMarker(true);
	dataSetStart.SetMarkerSize(2);
	dataSetStart.SetDatasetColor(1.0,0.0,0.0);
	CDataPoint point2(xcenter + visScale * globalTrack[0].x, ycenter + visScale * globalTrack[0].y);
	dataSetStart.AddDataPoint(point2);
	plot2D->AddDataSet(dataSetStart);

	// Now loop over all particles for local tracks
	for (int p = 0; p < pc; p++)
    {
		CDataSet fit;
		fit.SetDrawMarker(false);
		fit.SetDatasetColor(0.0,0.0,0.0);
		fit.SetLineWidth(1);

		for (int f = 0; f < fc; f++)
	    {
			CDataPoint point(visTracks[p][f].x, visTracks[p][f].y);
			fit.AddDataPoint(point);
		}
		plot2D->AddDataSet(fit);

		// Mark start of each track
		CDataSet patch_start;
		patch_start.SetDrawMarker(true);
		patch_start.SetMarkerSize(2);
		patch_start.SetDatasetColor(1.0,0.0,0.0);
		CDataPoint point3(visTracks[p][0].x, visTracks[p][0].y);
		patch_start.AddDataPoint(point3);
		plot2D->AddDataSet(patch_start);

    }


	char title[256];
	snprintf(title, 255, "X (in pixels; trajectory scaled by %.0f)", visScale);
	plot2D->SetXAxisTitle(title);
	title[0] = 'Y';
	plot2D->SetYAxisTitle(title);

	plot2D->OutputPostScriptPlot(fn_eps);


	// Compatibility with Jasenko's diagnostic .dat files
	// TODO: remove this
	if (!diag) return;

    std::ofstream rawOut(fn_root + "_tracks.dat");
    std::ofstream visOut(fn_root + "_visTracks.dat");
    std::ofstream visOut15(fn_root + "_visTracks_first15.dat");

    for (int p = 0; p < pc; p++)
    {
        rawOut << "#particle " << p << "\n";
        visOut << "#particle " << p << "\n";
        visOut15 << "#particle " << p << "\n";

        for (int f = 0; f < fc; f++)
        {
            rawOut << tracks[p][f].x << " " << tracks[p][f].y << "\n";
            visOut << visTracks[p][f].x << " " << visTracks[p][f].y << "\n";

            if (f < 15) visOut15 << visTracks[p][f].x << " " << visTracks[p][f].y << "\n";
        }

        rawOut << "\n";
        visOut << "\n";
        visOut15 << "\n";
    }

    std::ofstream glbOut(fn_root + "_globTrack.dat");

    for (int f = 0; f < fc; f++)
    {
        glbOut << globalTrack[f].x << " " << globalTrack[f].y << "\n";
    }
}

std::vector<Image<RFLOAT>> MotionFitter::weightsFromFCC()
{
    if (debug && verb > 0)
    	std::cout << " + Summing up FCCs...\n";

    Image<RFLOAT> fccData, fccWgh0, fccWgh1;
    Image<RFLOAT> fccDataMg, fccWgh0Mg, fccWgh1Mg;

    bool first = true;

    for (long g = 0; g < mdts.size(); g++)
    {
        FileName fn_root = getOutputFileNameRoot(g);

		fccDataMg.read(fn_root + "_FCC_cc.mrc");
		fccWgh0Mg.read(fn_root + "_FCC_w0.mrc");
		fccWgh1Mg.read(fn_root + "_FCC_w1.mrc");

		if (first)
		{
			sh = fccDataMg.data.xdim;
			s = 2 * (sh-1);
			fc = fccDataMg.data.ydim;

			fccData = Image<RFLOAT>(sh,fc);
			fccWgh0 = Image<RFLOAT>(sh,fc);
			fccWgh1 = Image<RFLOAT>(sh,fc);

			fccData.data.initZeros();
			fccWgh0.data.initZeros();
			fccWgh1.data.initZeros();

			first = false;
		}

        for (int y = 0; y < fc; y++)
        for (int x = 0; x < sh; x++)
        {
            fccData(y,x) += fccDataMg(y,x);
            fccWgh0(y,x) += fccWgh0Mg(y,x);
            fccWgh1(y,x) += fccWgh1Mg(y,x);
        }
    }

    Image<RFLOAT> fcc(sh,fc);

    for (int y = 0; y < fc; y++)
    for (int x = 0; x < sh; x++)
    {
        const double wgh = sqrt(fccWgh0Mg(y,x) * fccWgh1Mg(y,x));

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

    k0 = (int) (s*angpix / k0a);
    k1 = k1a > 0.0? (int) (s*angpix / k1a) : sh;

    if (verb > 0)
    	std::cout << " + Fitting B/k-factors between " << k0 << " and " << k1 << " pixels, or "
    	<< k0a << " and " << k1a << " Angstrom ...\n";

    std::pair<std::vector<d2Vector>,std::vector<double>> bkFacs
            = DamageHelper::fitBkFactors(fcc, k0, k1);

    std::vector<Image<RFLOAT>> freqWeights;
    freqWeights = DamageHelper::computeWeights(bkFacs.first, sh);

    const double cf = 8.0 * angpix*angpix * sh*sh;

    if (bfac_debug)
    {
        mktree(outPath + "/bfacs");

        Image<RFLOAT> bfacFit = DamageHelper::renderBkFit(bkFacs, sh, fc);
        Image<RFLOAT> bfacFitNoScale = DamageHelper::renderBkFit(bkFacs, sh, fc, true);

        ImageLog::write(bfacFit, outPath + "/bfacs/glob_Bk-fit");
        ImageLog::write(bfacFitNoScale, outPath + "/bfacs/glob_Bk-fit_noScale");
        ImageLog::write(fcc, outPath + "/bfacs/glob_Bk-data");
        ImageLog::write(freqWeights, outPath + "/bfacs/freqWeights");

        std::ofstream bfacsDat(outPath + "/bfacs/Bfac.dat");
        std::ofstream kfacsDat(outPath + "/bfacs/kfac.dat");

        for (int i = 0; i < fc; i++)
        {
            double s = bkFacs.first[i].x;
            double b = -cf/(s*s);

            bfacsDat << i << " " << b << "\n";
            kfacsDat << i << " " << log(bkFacs.first[i].y) << "\n";
        }

        bfacsDat.close();
        kfacsDat.close();
    }

    MetaDataTable mdt;
    mdt.setName("perframe_bfactors");

    for (int f = 0; f < fc; f++ )
    {
        double s = bkFacs.first[f].x;
        double b = -cf/(s*s);
        double k = log(bkFacs.first[f].y);

        mdt.addObject();
        mdt.setValue(EMDL_IMAGE_FRAME_NR, f);
        mdt.setValue(EMDL_POSTPROCESS_BFACTOR, b);
        mdt.setValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, k);
    }

    mdt.write(outPath + "/bfactors.star");

    // Also write out EPS plots of the B-factors and scale factors
    CPlot2D *plot2D=new CPlot2D("Polishing B-factors");
	plot2D->SetXAxisSize(600);
	plot2D->SetYAxisSize(400);
	plot2D->SetDrawLegend(false);
	plot2D->SetXAxisTitle("movie frame");
	plot2D->SetYAxisTitle("B-factor");
	mdt.addToCPlot2D(plot2D, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_BFACTOR);
	plot2D->OutputPostScriptPlot(outPath + "bfactors.eps");

	CPlot2D *plot2Db=new CPlot2D("Polishing scale-factors");
	plot2Db->SetXAxisSize(600);
	plot2Db->SetYAxisSize(400);
	plot2Db->SetDrawLegend(false);
	plot2Db->SetXAxisTitle("movie frame");
	plot2Db->SetYAxisTitle("Scale-factor");
	mdt.addToCPlot2D(plot2Db, EMDL_IMAGE_FRAME_NR, EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT);
	plot2Db->OutputPostScriptPlot(outPath + "scalefactors.eps");

    return freqWeights;
}

std::vector<Image<RFLOAT>> MotionFitter::weightsFromBfacs()
{
    MetaDataTable mdt;
    mdt.read(bfacFn);

    fc = mdt.numberOfObjects();

    std::vector<d2Vector> bkFacs(fc);

    const double cf = 8.0 * angpix*angpix * sh*sh;

    for (int f = 0; f < fc; f++)
    {
        int ff;
        mdt.getValue(EMDL_IMAGE_FRAME_NR, ff);

        double b, k;
        mdt.getValue(EMDL_POSTPROCESS_BFACTOR, b, f);
        mdt.getValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, k, f);

        bkFacs[f] = d2Vector(sqrt(-cf/b), exp(k));
    }

    std::vector<Image<RFLOAT>> freqWeights;
    freqWeights = DamageHelper::computeWeights(bkFacs, sh);

    if (bfac_debug)
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
            double s = bkFacs[i].x;
            double b = -cf/(s*s);

            bfacsDat << i << " " << b << "\n";
            kfacsDat << i << " " << log(bkFacs[i].y) << "\n";
        }

        bfacsDat.close();
        kfacsDat.close();
    }

    return freqWeights;
}



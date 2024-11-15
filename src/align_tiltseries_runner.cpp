/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
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
#include "src/align_tiltseries_runner.h"

void AlignTiltseriesRunner::read(int argc, char **argv, int rank)
{
	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "STAR file with all input tomograms, or a unix wildcard to all tomogram files, e.g. \"mics/*.mrc\"");
	fn_out = parser.getOption("--o", "Directory, where all output files will be stored", "AlignTiltSeries/");
    tomogram_thickness = textToFloat(parser.getOption("--tomogram_thickness", "Estimated tomogram thickness (in nm)", "300"));
	continue_old = parser.checkOption("--only_do_unfinished", "Only estimate CTFs for those tomograms for which there is not yet a logfile with Final values.");
    do_at_most = textToInteger(parser.getOption("--do_at_most", "Only process up to this number of (unprocessed) tomograms.", "-1"));
    fn_batchtomo_exe = parser.getOption("--batchtomo_exe", "IMOD's batchruntomo executable (default is set through $RELION_BATCHTOMO_EXECUTABLE)", "");
    fn_adoc_template = parser.getOption("--adoc_template", "Filename of a template for the IMOD directive file (by default, use hard-coded defaults)", "");

    int fid_section = parser.addSection("IMOD fiducial-based alignment options");
    do_imod_fiducials = parser.checkOption("--imod_fiducials", "Use IMOD's fiducial-based alignment method");
    fiducial_diam = textToFloat(parser.getOption("--fiducial_diameter", "Diameter of the fiducials (in nm)", "10"));

    int pat_section = parser.addSection("IMOD patch-tracking alignment options");
    do_imod_patchtrack = parser.checkOption("--imod_patchtrack", "OR: Use IMOD's patrick-tracking alignment method");
    patch_overlap = textToFloat(parser.getOption("--patch_overlap", "Overlap between the patches (in %)", "10."));
    patch_size = textToFloat(parser.getOption("--patch_size", "Patch size (in nm)", "100."));

    int aretomo2_section = parser.addSection("AreTomo2 alignment options");
    do_aretomo = parser.checkOption("--aretomo2", "OR: Use AreTomo2 alignment method");
    fn_aretomo_exe = parser.getOption("--aretomo_exe", "AreTomo executable (can be set through $RELION_ARETOMO_EXECUTABLE, defaults to AreTomo2)", "");
    do_aretomo_tiltcorrect = parser.checkOption("--aretomo_tiltcorrect", "Specify to correct the tilt angle offset in the tomogram (AreTomo -TiltCor option; default=false)");
    aretomo_tilcorrect_angle = textToFloat(parser.getOption("--aretomo_tiltcorrect_angle", "User-specified tilt angle correction (value > 180, means estimate automatically", "999."));
    do_aretomo_ctf = parser.checkOption("--aretomo_ctf", "Perform CTF estimation in AreTomo? (default=false)");
    do_aretomo_phaseshift = parser.checkOption("--aretomo_phaseshift", "Perform CTF estimation in AreTomo? (default=false)");
    gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread, e.g 0:1:2:3", "");

    int exp_section = parser.addSection("Expert options");
    other_wrapper_args  = parser.getOption("--other_wrapper_args", "Additional command-line arguments that will be passed onto the wrapper.", "");

    // Initialise verb for non-parallel execution
	verb = 1;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void AlignTiltseriesRunner::usage()
{
	parser.writeUsage(std::cout);
}

void AlignTiltseriesRunner::initialise(bool is_leader)
{
	// Get the Imod wrapper executable
	if (fn_batchtomo_exe == "")
	{
		char *penv;
		penv = getenv("RELION_BATCHTOMO_EXECUTABLE");
		if (penv != NULL)
        {
            fn_batchtomo_exe = (std::string)penv;
        }
        else
        {
            fn_batchtomo_exe = "batchruntomo";
        }
	}

    if (fn_aretomo_exe == "")
    {
        char *penv;
        penv = getenv("RELION_ARETOMO_EXECUTABLE");
        if (penv != NULL)
        {
            fn_aretomo_exe = (std::string)penv;
        }
        else
        {
            fn_aretomo_exe = "AreTomo2";
        }
    }


    if (fn_adoc_template != "")
    {
        std::ifstream in(fn_adoc_template.data(), std::ios_base::in);
        if (in.fail()) REPORT_ERROR("ERROR: cannot open directives template file: " + fn_adoc_template);

        in.seekg(0);
        std::string line;
        my_adoc_template = "";
        while (getline(in, line, '\n'))
        {
            my_adoc_template += line + "\n";
        }
        in.close();
    }
    else if (do_imod_fiducials)
    {
        my_adoc_template = fiducial_directive;
    }
    else if (do_imod_patchtrack)
    {
        my_adoc_template = patchtrack_directive;
    }


    int i = 0;
    if (do_imod_fiducials) i++;
    if (do_imod_patchtrack) i++;
    if (do_aretomo) i++;
    if (i != 1) REPORT_ERROR("ERROR: you need to specify one of these options: --imod_fiducials or --imod_patchtrack or --aretomo");

	// Make sure fn_out ends with a slash
	if (fn_out[fn_out.length()-1] != '/')
		fn_out += "/";

    // Check if this is a TomographyExperiment starfile, if not raise an error
    if (!tomogramSet.read(fn_in, 1))
    {
        REPORT_ERROR("ERROR: the input file is not a valid tilt series star file");
    }

	idx_tomograms_all.clear();
	idx_tomograms.clear();
	bool warned = false;
	for (long int itomo = 0; itomo < tomogramSet.size(); itomo++)
	{
        FileName fn_star;
        tomogramSet.globalTable.getValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_star, itomo);
        FileName fn_newstar = getOutputFileWithNewUniqueDate(fn_star, fn_out);
        tomogramSet.globalTable.setValue(EMDL_TOMO_TILT_SERIES_STARFILE, fn_newstar, itomo);

		bool process_this = true;
        bool ignore_this = false;
        if (continue_old)
		{
			if (checkResults(itomo))
			{
				process_this = false; // already done
			}
		}

		if (do_at_most >= 0 && idx_tomograms.size() >= do_at_most)
		{
			if (process_this) {
				ignore_this = true;
				process_this = false;
				if (!warned)
				{
					warned = true;
					std::cout << "NOTE: processing of some tomograms will be skipped as requested by --do_at_most" << std::endl;
				}
			}
			// If this tomogram has already been processed, the result should be included in the output.
			// So ignore_this remains false.
		}

		if (process_this)
		{
			idx_tomograms.push_back(itomo);
		}

		if (!ignore_this)
		{
			idx_tomograms_all.push_back(itomo);
		}
	}

	if (is_leader && do_at_most >= 0 )
	{
		std::cout << tomogramSet.size() << " tomograms were given in the input tomogram set, but we process only ";
		std::cout  << do_at_most << " tomograms as specified in --do_at_most." << std::endl;
	}

    if (do_aretomo)
    {
        if (gpu_ids.length() > 0)
            untangleDeviceIDs(gpu_ids, allThreadIDs);
        else if (verb>0)
            std::cout << "WARNING: --gpu_ids not specified, threads will automatically be mapped to devices."<< std::endl;
    }

    if (verb > 0)
	{
        if (do_aretomo)
            std::cout << " Using AreTomo executable in: " << fn_aretomo_exe << std::endl;
        else
            std::cout << " Using batchruntomo executable in: " << fn_batchtomo_exe << std::endl;
		std::cout << " to align tilt series for the following tomograms: " << std::endl;
		if (continue_old)
			std::cout << " (skipping all tomograms for output files with tilt series alignment parameters already exists)" << std::endl;
		for (unsigned  int  i = 0; i < idx_tomograms.size(); ++i)
			std::cout << "  * " << tomogramSet.getTomogramName(idx_tomograms[i]) << std::endl;
	}
}

void AlignTiltseriesRunner::run()
{

    int barstep;
    if (verb > 0)
    {
        std::cout << " Aligning tilt series ..." << std::endl;
        init_progress_bar(idx_tomograms.size());
        barstep = XMIPP_MAX(1, idx_tomograms.size() / 60);
    }

    std::vector<std::string> alltomonames;
    for (long int itomo = 0; itomo < idx_tomograms.size(); itomo++)
    {

        // Abort through the pipeline_control system
        if (pipeline_control_check_abort_job())
            exit(RELION_EXIT_ABORTED);

        if (do_aretomo)
        {
            executeAreTomo(idx_tomograms[itomo]);
        }
        else if (do_imod_fiducials || do_imod_patchtrack)
        {
            //executeImodWrapper(idx_tomograms[itomo]);
            executeIMOD(idx_tomograms[itomo]);
        }

        if (verb > 0 && itomo % barstep == 0)
            progress_bar(itomo);
    }

    if (verb > 0)
        progress_bar(idx_tomograms.size());

}

bool AlignTiltseriesRunner::checkResults(long idx_tomo)
{

    std::string tomoname = tomogramSet.getTomogramName(idx_tomo);
    FileName fn_dir = fn_out + "external/" + tomoname + '/';

    if (do_aretomo)
    {
        // check that .aln (and _ctf.txt if do_aretomo_ctf) file(s) has been written out
        FileName fn_aln = fn_dir + tomoname + ".aln";
        FileName fn_ctf = fn_dir + tomoname + "_ctf.txt";

        if (do_aretomo_ctf)
        {
            return (exists(fn_ctf) && exists(fn_aln));
        }
        else
        {
            return exists(fn_aln);
        }

    }
    else
    {
        // check that .xf and .tlt files have been written out
        FileName fn_xf = fn_dir + tomoname + ".xf";
        FileName fn_tlt = fn_dir + tomoname + ".tlt";

        return (exists(fn_xf) && exists(fn_tlt));

    }

    return false;
}

void AlignTiltseriesRunner::generateMRCStackAndRawTiltFile(long idx_tomo, bool is_aretomo)
{
    FileName fn_dir = fn_out + "external/" + tomogramSet.getTomogramName(idx_tomo) + '/';
    FileName fn_tilt = fn_dir + tomogramSet.getTomogramName(idx_tomo) + ".rawtlt";
    FileName fn_series = fn_dir + tomogramSet.getTomogramName(idx_tomo) + ".mrc";

    std::ofstream  fh;
    fh.open((fn_tilt).c_str(), std::ios::out);
    Image<RFLOAT> Imic, Iseries;
    std::vector<int> frame_dose_order = tomogramSet.getFrameDoseOrder(idx_tomo);

    int fc = tomogramSet.tomogramTables[idx_tomo].numberOfObjects();
    for (int f = 0; f < fc; f++)
    {
        FileName fn_mic;
        RFLOAT tiltangle;
        tomogramSet.tomogramTables[idx_tomo].getValue(EMDL_MICROGRAPH_NAME, fn_mic, f);
        tomogramSet.tomogramTables[idx_tomo].getValue(EMDL_TOMO_NOMINAL_TILT_STAGE_ANGLE, tiltangle, f);
        RFLOAT dose;
        tomogramSet.tomogramTables[idx_tomo].getValue(EMDL_MICROGRAPH_PRE_EXPOSURE, dose, f);
        // IMOD: 1-column and AreTomo2 2-column raw tilt angle file with tilt angles (and order of acquisition).
        fh << tiltangle;
        if (is_aretomo) fh << " " << frame_dose_order[f];
        fh << std::endl;
        Imic.read(fn_mic);
        if (f == 0) Iseries().resize(fc, YSIZE(Imic()), XSIZE(Imic()));
        Iseries().setSlice(f, Imic());

    }

    RFLOAT angpix = tomogramSet.getTiltSeriesPixelSize(idx_tomo);
    tomogramSet.globalTable.getDouble(EMDL_TOMO_IMPORT_FRACT_DOSE, idx_tomo);
    if (is_aretomo) Iseries.setSamplingRateInHeader(angpix);
    Iseries.write(fn_series);
    fh.close();




}

void AlignTiltseriesRunner::executeIMOD(long idx_tomo, int rank)
{

    // Generate external output directory and write input files for AreTomo
    std::string tomoname = tomogramSet.getTomogramName(idx_tomo);
    RFLOAT mythickness = tomogram_thickness;
    if (tomogramSet.globalTable.containsLabel(EMDL_TOMO_TOMOGRAM_THICKNESS))
        tomogramSet.globalTable.getValueSafely(EMDL_TOMO_TOMOGRAM_THICKNESS, mythickness, idx_tomo);
    FileName fn_dir = fn_out + "external/" + tomoname + '/';
    mktree(fn_dir);

    FileName fn_series = fn_dir + tomoname + ".mrc";
    FileName fn_tilt = fn_dir + tomoname + ".rawtlt";
    FileName fn_adoc = fn_dir + "batchDirective.adoc";
    FileName fn_log = fn_dir + "log.txt"; // for compatbility with Alister's wrapper
    FileName fn_com = fn_dir + tomoname + ".com";

    std::ofstream  fhadoc;
    fhadoc.open((fn_adoc).c_str(), std::ios::out);

    RFLOAT pixel_size = tomogramSet.getTiltSeriesPixelSize(idx_tomo);
    RFLOAT rotangle = tomogramSet.tomogramTables[idx_tomo].getDouble(EMDL_TOMO_NOMINAL_TILT_AXIS_ANGLE, 0);
    if (do_imod_fiducials)
    {
        fhadoc << my_adoc_template;
        fhadoc << "setupset.copyarg.rotation = " << rotangle << std::endl;
        fhadoc << "setupset.copyarg.pixel = " << pixel_size/10. << std::endl;
        fhadoc << "setupset.copyarg.gold = " << fiducial_diam << std::endl;
    }
    else if (do_imod_patchtrack)
    {
        // Find best power-of-2 binning factor that brings pixel size to 10A
        int mybinning = 1, mybestbinning = 1;
        RFLOAT mindiff = fabs(pixel_size * mybinning - 10.); // 10A is the target binned pixel size
        for (int i = 0; i < 8; i++)
        {
            mybinning *= 2;
            //std::cerr << " mybinning= " << mybinning << " diff =" << fabs(pixel_size * mybinning - 10.) <<" mindiff= " << mindiff <<" bestbin= " << mybestbinning << std::endl;
            if (fabs(pixel_size * mybinning - 10.) < mindiff)
            {
                mindiff = fabs(pixel_size * mybinning - 10.);
                mybestbinning = mybinning;
            }
        }
        int binned_patch_size = ROUND( (10. * patch_size) / (pixel_size * mybestbinning) );

        fhadoc << my_adoc_template;
        fhadoc << "setupset.copyarg.rotation = " << rotangle << std::endl;
        fhadoc << "setupset.copyarg.pixel = " << pixel_size/10. << std::endl;
        fhadoc << "comparam.prenewst.newstack.BinByFactor = " << mybestbinning << std::endl;
        fhadoc << "comparam.xcorr_pt.tiltxcorr.SizeOfPatchesXandY = " << binned_patch_size << "," << binned_patch_size << std::endl;
        fhadoc << "comparam.xcorr_pt.tiltxcorr.OverlapOfPatchesXandY = " << patch_overlap/100. << "," << patch_overlap/100. << std::endl;
        fhadoc << "comparam.tilt.tilt.THICKNESS = " << tomogram_thickness * 10. << std::endl;
    }
    else
    {
        REPORT_ERROR("ERROR: either do_imod_fiducials or do_imod_patchtrack should be true.");
    }
    fhadoc.close();


    // Make sure metadata table is sorted on rlnTomoNominalStageTiltAngle (it should be, but anyways...)
    tomogramSet.tomogramTables[idx_tomo].newSort(EMDL_TOMO_NOMINAL_TILT_STAGE_ANGLE);

    generateMRCStackAndRawTiltFile(idx_tomo, false);

    // Now run the actual IMOD command
    std::string command = fn_batchtomo_exe;
    command += " -DirectiveFile " + fn_adoc;
    command += " -CurrentLocation " + fn_dir;
    command += " -RootName " + tomoname;
    command += " -EndingStep 6";

    if (other_wrapper_args.length() > 0)
        command += " " + other_wrapper_args;

    command += " > " + fn_log + " 2>&1 ";

    // Write the command to a .com file
    std::ofstream  fhc;
    fhc.open((fn_com).c_str(), std::ios::out);
    fhc << command << std::endl;
    fhc.close();

    if (system(command.c_str()))
    {
        std::cerr << "WARNING: there was an error in executing: " << command << std::endl;
    }


}

void AlignTiltseriesRunner::executeAreTomo(long idx_tomo, int rank)
{

    // Generate external output directory and write input files for AreTomo
    std::string tomoname = tomogramSet.getTomogramName(idx_tomo);
    FileName fn_dir = fn_out + "external/" + tomoname + '/';
    mktree(fn_dir);

    FileName fn_series = fn_dir + tomoname + ".mrc";
    FileName fn_tilt = fn_dir + tomoname + ".rawtlt";
    FileName fn_ali = fn_dir + tomoname + "_aligned.mrc";
    FileName fn_log = fn_dir + tomoname + ".log";
    FileName fn_com = fn_dir + tomoname + ".com";

    std::ofstream  fh;
    fh.open((fn_tilt).c_str(), std::ios::out);
    Image<RFLOAT> Imic, Iseries;
    int fc = tomogramSet.tomogramTables[idx_tomo].numberOfObjects();
    std::vector<int> frame_dose_order = tomogramSet.getFrameDoseOrder(idx_tomo);

    generateMRCStackAndRawTiltFile(idx_tomo, true);

    RFLOAT frac_dose = tomogramSet.globalTable.getDouble(EMDL_TOMO_IMPORT_FRACT_DOSE, idx_tomo);
    RFLOAT pixel_size = tomogramSet.getTiltSeriesPixelSize(idx_tomo);
    // Tomogram_thickness is in nm, should become unbinned pixels
    RFLOAT mythickness = tomogram_thickness;
    if (tomogramSet.globalTable.containsLabel(EMDL_TOMO_TOMOGRAM_THICKNESS))
        tomogramSet.globalTable.getValueSafely(EMDL_TOMO_TOMOGRAM_THICKNESS, mythickness, idx_tomo);
    RFLOAT thickness_pix = mythickness*10./pixel_size;

    // Now run the actual AreTomo command
    std::string command = fn_aretomo_exe + " ";
    command += " -InMrc " + fn_series;
    command += " -AngFile " + fn_tilt;
    command += " -OutMrc " + fn_ali;
    command += " -ImgDose " + floatToString(frac_dose);
    // Tomogram thickness should be in unbinned pixels
    command += " -AlignZ " + floatToString(thickness_pix);
    // Skip reconstruction of the tomogram in AreTomo...
    command += " -volZ 0";

    if (tomogramSet.tomogramTables[idx_tomo].containsLabel(EMDL_TOMO_NOMINAL_TILT_AXIS_ANGLE))
    {
        RFLOAT tiltaxis_angle = tomogramSet.tomogramTables[idx_tomo].getDouble(EMDL_TOMO_NOMINAL_TILT_AXIS_ANGLE, 0);
        command += " -TiltAxis " + floatToString(tiltaxis_angle);
    }

    if (do_aretomo_tiltcorrect)
    {
        command += " -TiltCor 1 ";
        if (aretomo_tilcorrect_angle < 180.)
            command += floatToString(aretomo_tilcorrect_angle);
    }
    else
    {
        command += " -TiltCor -1 ";
    }


    if (do_aretomo_ctf)
    {
        // Also estimate CTF parameters in AreTomo
        RFLOAT kV, Cs, Q0;
        RFLOAT angpix = tomogramSet.getTiltSeriesPixelSize(idx_tomo);
        tomogramSet.globalTable.getValue(EMDL_CTF_VOLTAGE, kV, idx_tomo);
        tomogramSet.globalTable.getValue(EMDL_CTF_CS, Cs, idx_tomo);
        tomogramSet.globalTable.getValue(EMDL_CTF_Q0, Q0, idx_tomo);
        command += " -Kv " + floatToString(kV);
        command += " -Cs " + floatToString(Cs);
        command += " -AmpContrast " + floatToString(Q0);
        command += " -PixSize " + floatToString(angpix);
        if (do_aretomo_phaseshift)
            command += " -ExtPhase 90 180";
    }

    if (gpu_ids.length() > 0)
    {
        if (rank >= allThreadIDs.size())
            REPORT_ERROR("ERROR: not enough MPI nodes specified for the GPU IDs.");

        command += " -Gpu " ;
        for (int igpu = 0; igpu < allThreadIDs[rank].size(); igpu++)
        {
            command += allThreadIDs[rank][igpu] + " ";
        }
    }

    if (other_wrapper_args.length() > 0)
        command += " " + other_wrapper_args;

    command += " > " + fn_log + " 2>&1 ";

    // Write the command to a .com file
    std::ofstream  fhc;
    fhc.open((fn_com).c_str(), std::ios::out);
    fhc << command << std::endl;
    fhc.close();

    if (system(command.c_str()))
    {
        std::cerr << "WARNING: there was an error in executing: " << command << std::endl;
    }

}

void AlignTiltseriesRunner::makePerTiltSeriesEPSFiles(long idx_tomo, bool do_ctf)
{

    std::string tomoname = tomogramSet.getTomogramName(idx_tomo);
    FileName fn_dir = fn_out + "external/" + tomoname + '/';
    FileName fn_eps_tilt = fn_dir + tomoname + "_tilts.eps";
    FileName fn_eps_shift = fn_dir + tomoname + "_shifts.eps";
    FileName fn_eps_defocus = fn_dir + tomoname + "_defocus.eps";
    FileName fn_eps_ctffom = fn_dir + tomoname + "_ctfscore.eps";

    if (continue_old && exists(fn_eps_tilt) && exists(fn_eps_shift)) return;

    // Nominal and refined tilt values for all images in this tiltseries
    CPlot2D *plot2Da=new CPlot2D("nominal and refined tilt angles for " + tomoname);
    tomogramSet.tomogramTables[idx_tomo].addToCPlot2D(plot2Da, EMDL_UNDEFINED, EMDL_TOMO_NOMINAL_TILT_STAGE_ANGLE, 0., 0., 1.);
    tomogramSet.tomogramTables[idx_tomo].addToCPlot2D(plot2Da, EMDL_UNDEFINED, EMDL_TOMO_YTILT, 1., 0., 0.);
    plot2Da->SetDrawLegend(true);
    plot2Da->OutputPostScriptPlot(fn_eps_tilt);
    delete plot2Da;

    // Refined X and Y-shift values for all images in this tiltseries
    CPlot2D *plot2Db=new CPlot2D("X and Y shifts for " + tomoname);
    tomogramSet.tomogramTables[idx_tomo].addToCPlot2D(plot2Db, EMDL_UNDEFINED, EMDL_TOMO_XSHIFT_ANGST, 0., 0., 1.);
    tomogramSet.tomogramTables[idx_tomo].addToCPlot2D(plot2Db, EMDL_UNDEFINED, EMDL_TOMO_YSHIFT_ANGST, 1., 0., 0.);
    plot2Db->SetDrawLegend(true);
    plot2Db->OutputPostScriptPlot(fn_eps_shift);
    delete plot2Db;

    if (do_ctf)
    {
        CPlot2D *plot2Dc=new CPlot2D("Defocus U and V for " + tomoname);
        tomogramSet.tomogramTables[idx_tomo].addToCPlot2D(plot2Dc, EMDL_UNDEFINED, EMDL_CTF_DEFOCUSU, 0., 0., 1.);
        tomogramSet.tomogramTables[idx_tomo].addToCPlot2D(plot2Dc, EMDL_UNDEFINED, EMDL_CTF_DEFOCUSV, 1., 0., 0.);
        plot2Dc->SetDrawLegend(true);
        plot2Dc->OutputPostScriptPlot(fn_eps_defocus);
        delete plot2Dc;

        CPlot2D *plot2Dd=new CPlot2D("CTF score for " + tomoname);
        tomogramSet.tomogramTables[idx_tomo].addToCPlot2D(plot2Dd, EMDL_UNDEFINED, EMDL_CTF_FOM, 0., 0., 1.);
        plot2Dd->OutputPostScriptPlot(fn_eps_ctffom);
        delete plot2Dd;

    }


}


bool AlignTiltseriesRunner::readIMODResults(long idx_tomo, std::string &error_message)
{
    std::string tomoname = tomogramSet.getTomogramName(idx_tomo);
    FileName fn_dir = fn_out + "external/" + tomoname + '/';
    FileName fn_rawtlt = fn_dir + tomoname + ".rawtlt";
    FileName fn_xf = fn_dir + tomoname + ".xf";
    FileName fn_tlt = fn_dir + tomoname + ".tlt";
    FileName fn_edf = fn_dir + tomoname + ".edf";
    FileName fn_align = fn_dir + "align.log";

    // In case of MPI, the tomogramTable may not yet have been sorted on the tilt angle. Do it again now
    tomogramSet.tomogramTables[idx_tomo].newSort(EMDL_TOMO_NOMINAL_TILT_STAGE_ANGLE);

    int fc = tomogramSet.tomogramTables[idx_tomo].numberOfObjects();
    RFLOAT angpix = tomogramSet.getTiltSeriesPixelSize(idx_tomo);

    std::string line;
    std::vector<std::string> words;

    // 1. Get tiltangle_offset from the align.log file
    RFLOAT tiltangle_offset = 0., mean_error = 0., stddev_error = 0., leaveout_error = 0.;
    std::ifstream in0(fn_align.data(), std::ios_base::in);
    if (in0.fail())
    {
        error_message = " ERROR: cannot open IMOD's align.log file: " + fn_align;
        return false;
    }

    in0.seekg(0);
    while (getline(in0, line, '\n'))
    {
        if (line.find("AngleOffset = ") != std::string::npos)
        {
            tokenize(line, words);
            if (words.size() < 3) REPORT_ERROR("ERROR: fewer than 3 columns on this line from align.log file: " + line);
            tiltangle_offset = textToFloat(words[2]);
        }
        else if (line.find("Residual error mean and sd:") != std::string::npos)
        {
            tokenize(line, words);
            if (words.size() < 7) REPORT_ERROR("ERROR: fewer than 7 columns on this line from align.log file: " + line);
            mean_error = textToFloat(words[5]);
            stddev_error = textToFloat(words[6]);
        }
        else if (line.find("Global leave-out error") != std::string::npos)
        {
            tokenize(line, words);
            if (words.size() < 6) REPORT_ERROR("ERROR: fewer than 6 columns on this line from align.log file: " + line);
            leaveout_error = textToFloat(words[5]);
        }

    }
    in0.close();

    // 2. Get specimen_shifts and zrot from the .xf file
    std::ifstream in(fn_xf.data(), std::ios_base::in);
    if (in.fail())
    {
        error_message = " ERROR: cannot open IMOD's .xf file: " + fn_xf;
        return false;
    }

    in.seekg(0);
    int f = 0;
    while (getline(in, line, '\n'))
    {
        Matrix2D<RFLOAT> A(2,2);
        Matrix1D<RFLOAT> s(2), image_shifts(2);
        tokenize(line, words);
        if (words.size() != 6) REPORT_ERROR("ERROR: did not find 6 columns in xf file: " + fn_xf);
        A(0, 0) = textToFloat(words[0]);
        A(0, 1) = textToFloat(words[1]);
        A(1, 0) = textToFloat(words[2]);
        A(1, 1) = textToFloat(words[3]);
        s(0) =  textToFloat(words[4]);
        s(1) =  textToFloat(words[5]);

        // get specimen shifts
        image_shifts = A.inv() * s;
        tomogramSet.tomogramTables[idx_tomo].setValue(EMDL_TOMO_XSHIFT_ANGST, -image_shifts(0) * angpix, f);
        tomogramSet.tomogramTables[idx_tomo].setValue(EMDL_TOMO_YSHIFT_ANGST, -image_shifts(1) * angpix, f);

        // get in-plane rotation (zrot)
        RFLOAT myzrot = RAD2DEG(acos(A(0, 0)));
        RFLOAT inizrot;
        tomogramSet.tomogramTables[idx_tomo].getValue(EMDL_TOMO_NOMINAL_TILT_AXIS_ANGLE, inizrot, f);
        RFLOAT diff = fabs(inizrot - myzrot);
        RFLOAT flippeddiff = fabs(-inizrot - myzrot);
        if (flippeddiff < diff) myzrot = -myzrot;
        tomogramSet.tomogramTables[idx_tomo].setValue(EMDL_TOMO_ZROT, myzrot, f);

        f++;
    }
    in.close();

    if (f != fc) REPORT_ERROR("ERROR: found " + integerToString(f) + " entries in xf file " + fn_xf + ", but expected " +
                                                 integerToString(fc) + " entries...");

    // 3. Get ytilt angles straight from the .tlt file (and apply tiltangle_offset)
    std::ifstream in2(fn_tlt.data(), std::ios_base::in);
    if (in2.fail())
    {
        error_message = " ERROR: cannot open IMOD's .tlt file: " + fn_tlt;
        return false;
    }

    in2.seekg(0);
    f = 0;
    while (getline(in2, line, '\n'))
    {
        tokenize(line, words);
        if (words.size() != 1) REPORT_ERROR("ERROR: did not find 1 column in tlt file: " + fn_tlt);
        tomogramSet.tomogramTables[idx_tomo].setValue(EMDL_TOMO_XTILT, 0., f);
        // Subtract tilt_angle_offset
        tomogramSet.tomogramTables[idx_tomo].setValue(EMDL_TOMO_YTILT, textToFloat(words[0]) - tiltangle_offset, f);

        f++;
    }
    in2.close();

    if (f != fc) REPORT_ERROR("ERROR: found " + integerToString(f) + " entries in tlt file " + fn_tlt + ", but expected " +
                              integerToString(fc) + " entries...");

    // Also set the edf filename in the overall tilt_series STAR file
    tomogramSet.globalTable.setValue(EMDL_TOMO_ETOMO_DIRECTIVE_FILE, fn_edf, idx_tomo);
    tomogramSet.globalTable.setValue(EMDL_TOMO_IMOD_ERROR_MEAN, mean_error, idx_tomo);
    tomogramSet.globalTable.setValue(EMDL_TOMO_IMOD_ERROR_STDDEV, stddev_error, idx_tomo);
    tomogramSet.globalTable.setValue(EMDL_TOMO_IMOD_LEAVEOUT_ERROR, leaveout_error, idx_tomo);

    // Also make per-tiltseries EPS files
    makePerTiltSeriesEPSFiles(idx_tomo);

    return true;

}

bool AlignTiltseriesRunner::readAreTomoResults(long idx_tomo, std::string &error_message)
{

    std::string tomoname = tomogramSet.getTomogramName(idx_tomo);
    FileName fn_dir = fn_out + "external/" + tomoname + '/';
    FileName fn_aln = fn_dir + tomoname + ".aln";
    FileName fn_ctf_img = fn_dir + tomoname + "_ctf.mrc";
    FileName fn_ctf = fn_dir + tomoname + "_ctf.txt";
    FileName fn_log = fn_dir + tomoname + ".log";

    int fc = tomogramSet.tomogramTables[idx_tomo].numberOfObjects();
    RFLOAT angpix = tomogramSet.getTiltSeriesPixelSize(idx_tomo);

    std::string line;
    std::vector<std::string> words;


    // 0. Get the best tilt axis score the log file
    RFLOAT tiltaxis_score = -0.999;
    std::ifstream in0(fn_log.data(), std::ios_base::in);
    if (in0.fail())
    {
        error_message = " ERROR: cannot open the AreTomo2 log file: " + fn_log;
        return false;
    }

    in0.seekg(0);
    while (getline(in0, line, '\n'))
    {
        if (line.find("Best tilt axis:") != std::string::npos)
        {
            tokenize(line, words);
            if (words.size() < 6) REPORT_ERROR("ERROR: fewer than 6 columns on this line from logfile: " + line);
            if (textToFloat(words[5]) > tiltaxis_score)
                tiltaxis_score = textToFloat(words[5]);
        }
    }
    in0.close();

    // Set the best score in the global table
    tomogramSet.globalTable.setValue(EMDL_TOMO_ARETOMO_TILTAXIS_SCORE, tiltaxis_score, idx_tomo);


    // 1. Get tiltseries alignment parameters from the .aln file
    std::ifstream in(fn_aln.data(), std::ios_base::in);
    if (in.fail())
    {
        error_message = " ERROR: cannot open alignment file: " + fn_aln;
        return false;
    }

    std::vector<RFLOAT> rot, tilt, tx, ty;
    std::vector<RFLOAT> defU, defV, defAngle, phaseShift, corr, maxres;
    std::vector<int> indices, dark_frames;

    in.seekg(0);
    while (getline(in, line, '\n'))
    {
        // See if any dark frames were excluded by AreTomo
        if (line.find("# DarkFrame =") == 0)
        {
            tokenize(line, words);
            int idx = textToInteger(words[4]); // column 4 is the order in the original input starfile, column 3 is ordered by tilt angle
            //std::cerr << " dark idx= " << idx << std::endl;
            dark_frames.push_back(idx);
        }

        // find the header, which is the last line starting with a '#'
        if (line.find("#") != 0)
        {
            // Data lines are all lines without a leading #
            tokenize(line, words);
            int idx =  textToInteger(words[0]);
            if (idx < 0 || idx >= fc) REPORT_ERROR("BUG: idx= " + integerToString(idx) + " fc= " + integerToString(fc) + " from .aln file: " + fn_aln);
            indices.push_back(idx);
            rot.push_back(textToFloat(words[1]));
            tx.push_back(angpix * textToFloat(words[3]));
            ty.push_back(angpix * textToFloat(words[4]));
            tilt.push_back(textToFloat(words[9]));
        }
    }
    in.close();

    if (rot.size() != fc - dark_frames.size())
    {
        error_message = " ERROR_: unexpected number of data rows in parameter file: " + fn_aln + " : " +
                integerToString(rot.size()) + " (expected: " + integerToString(fc) + " - " + integerToString(dark_frames.size()) + " dark frames )";
        return false;
    }

    // 1. Get CTF estimation parameters from the _ctf.txt file
    if (do_aretomo_ctf)
    {

        // Set this for later writing out to power_spectra star file

        // Get CTF parameters from the _ctf.txt file
        std::ifstream in2(fn_ctf.data(), std::ios_base::in);
        if (in2.fail())
        {
            error_message = " ERROR: cannot open CTF parameter file: " + fn_ctf;
            return false;
        }

        in2.seekg(0);
        while (getline(in2, line, '\n'))
        {
            if (line.find("#") != 0)
            {

                tokenize(line, words);

                defU.push_back(textToFloat(words[1]));
                defV.push_back(textToFloat(words[2]));
                defAngle.push_back(textToFloat(words[3]));
                if (do_aretomo_phaseshift)
                {
                    phaseShift.push_back(textToFloat(words[4]));
                }
                corr.push_back(textToFloat(words[5]));
                maxres.push_back(textToFloat(words[6]));
            }
        }

        if (defU.size() != rot.size())
        {
            error_message = " ERROR_: unexpected number of data rows in CTF parameter file " + fn_ctf + ": "
                    + integerToString(defU.size()) + " (expected from .aln file = " + integerToString(rot.size()) + ")";
            return false;
        }
    }


    MetaDataTable MDnew;
    for (int i = 0; i < rot.size(); i++)
    {
        int f = indices[i];

        MDnew.addObject(tomogramSet.tomogramTables[idx_tomo].getObject(f));

        MDnew.setValue(EMDL_TOMO_XTILT, 0.);
        MDnew.setValue(EMDL_TOMO_YTILT, tilt[i]);
        MDnew.setValue(EMDL_TOMO_ZROT, rot[i]);
        // For AreTomo, also position one of the refined Zrot values in the globalTable, so it can be plotted for the logfile
        if (i==0) tomogramSet.globalTable.setValue(EMDL_TOMO_ZROT, rot[i], idx_tomo);
        MDnew.setValue(EMDL_TOMO_XSHIFT_ANGST, tx[i]);
        MDnew.setValue(EMDL_TOMO_YSHIFT_ANGST, ty[i]);

        if (do_aretomo_ctf)
        {

            MDnew.setValue(EMDL_CTF_DEFOCUSU, defU[i]);
            MDnew.setValue(EMDL_CTF_DEFOCUSV, defV[i]);
            MDnew.setValue(EMDL_CTF_DEFOCUS_ANGLE, defAngle[i]);
            if (do_aretomo_phaseshift) MDnew.setValue(EMDL_CTF_PHASESHIFT, phaseShift[i]);
            MDnew.setValue(EMDL_CTF_FOM, corr[i]);
            MDnew.setValue(EMDL_CTF_MAXRES, maxres[i]);
            FileName fn_img;
            fn_img.compose(i+1, fn_ctf_img);
            MDnew.setValue(EMDL_CTF_IMAGE, fn_img+":mrcs");

        }

    }

    MDnew.sort(EMDL_TOMO_NOMINAL_TILT_STAGE_ANGLE);
    MDnew.setName(tomogramSet.tomogramTables[idx_tomo].getName());
    tomogramSet.tomogramTables[idx_tomo] = MDnew;

    // Also make per-tiltseries EPS files
    makePerTiltSeriesEPSFiles(idx_tomo, do_aretomo_ctf);

    return true;

}

void AlignTiltseriesRunner::joinResults()
{
    // Check again the STAR file exists and has the right labels
    // Also check for the presence of any eTomoDirective files

    MetaDataTable MDout;
    MDout.setName("global");

    // Fill all the individual tilt series STAR files with the parameters from the AreTomo .aln files
    std::string error_message;
    MetaDataTable MDpower;
    std::vector<std::string> failed_tomograms;
    for (long itomo = 0; itomo < tomogramSet.size(); itomo++)
    {
        if (do_aretomo)
        {
            if (!readAreTomoResults(itomo, error_message))
            {
                std::string myname = tomogramSet.getTomogramName(itomo);
                std::cerr << " Error for reading AreTomo results from tomogram: " << myname << ":" << std::endl;
                std::cerr << error_message << std::endl;
                failed_tomograms.push_back(myname);
            }
            else
            {
                if (do_aretomo_ctf)
                {
                    FOR_ALL_OBJECTS_IN_METADATA_TABLE(tomogramSet.tomogramTables[itomo])
                    {
                        MDpower.addObject(tomogramSet.tomogramTables[itomo].getObject(current_object));
                    }
                }
            }
        }
        else if (do_imod_fiducials || do_imod_patchtrack)
        {
            if (!readIMODResults(itomo, error_message))
            {
                std::string myname = tomogramSet.getTomogramName(itomo);
                std::cerr << " Error for reading IMOD results from tomogram: " << myname << ":" << std::endl;
                std::cerr << error_message << std::endl;
                failed_tomograms.push_back(myname);
            }
        }
    }

    if (failed_tomograms.size() > 0)
    {
        if (failed_tomograms.size() == tomogramSet.size())
            REPORT_ERROR("ERROR: all tomograms failed alignment, exiting now... ");

        std::cout << " !!! WARNING: there have been " << failed_tomograms.size() << " tomograms for which alignment has failed. " << std::endl;
        std::cout << " !!! WARNING: the failed tomograms are: " << std::endl;
        for (int i =0; i < failed_tomograms.size(); i++)
            std::cout << " !!! WARNING:  - " + failed_tomograms[i] << std::endl;
        std::cout << " !!! WARNING: these failed tomograms will not be part of the output STAR file..." << std::endl;
        std::cout << " !!! WARNING: you may want to see whether you can solve the errors above in order not to loose these data." << std::endl;

        for (int i =0; i < failed_tomograms.size(); i++)
        {
            tomogramSet.removeTomogram(failed_tomograms[i]);
        }
    }

    tomogramSet.write(fn_out+"aligned_tilt_series.star");

    if (do_aretomo && do_aretomo_ctf)
    {
        if (verb > 0) std::cout << " Saving a file called " << fn_out << "power_spectra_fits.star for visualisation of Thon ring fits..." << std::endl;
        MDpower.deactivateLabel(EMDL_MICROGRAPH_NAME);
        MDpower.deactivateLabel(EMDL_MICROGRAPH_MOVIE_NAME);
        MDpower.write(fn_out+"power_spectra_fits.star");
    }

    if (verb > 0)
    {
        std::cout << " Written out: " << fn_out <<  "aligned_tilt_series.star" << std::endl;
    }



    if (verb > 0)
    {
        std::cout << " Now generating logfile.pdf ... " << std::endl;
    }

    std::vector<EMDLabel> plot_labels;
    if (do_imod_fiducials || do_imod_patchtrack)
    {
        plot_labels.push_back(EMDL_TOMO_IMOD_LEAVEOUT_ERROR);
        plot_labels.push_back(EMDL_TOMO_IMOD_ERROR_MEAN);
        plot_labels.push_back(EMDL_TOMO_IMOD_ERROR_STDDEV);
    }
    else if (do_aretomo)
    {
        plot_labels.push_back(EMDL_TOMO_ZROT);
        plot_labels.push_back(EMDL_TOMO_ARETOMO_TILTAXIS_SCORE);
        if (do_aretomo_tiltcorrect && aretomo_tilcorrect_angle > 180.)
            plot_labels.push_back(EMDL_TOMO_ARETOMO_TILTANGLE_OFFSET);
    }
    FileName fn_eps;
    std::vector<FileName> all_fn_eps;
    for (int i = 0; i < plot_labels.size(); i++)
    {
        EMDLabel label = plot_labels[i];
        if (tomogramSet.globalTable.containsLabel(label))
        {
            // Values for all micrographs
            CPlot2D *plot2Db=new CPlot2D(EMDL::label2Str(label) + " for all tiltseries");
            tomogramSet.globalTable.addToCPlot2D(plot2Db, EMDL_UNDEFINED, label, 1.);
            plot2Db->SetDrawLegend(false);
            fn_eps = fn_out + "all_" + EMDL::label2Str(label) + ".eps";
            plot2Db->OutputPostScriptPlot(fn_eps);
            all_fn_eps.push_back(fn_eps);
            delete plot2Db;
            if ( tomogramSet.globalTable.numberOfObjects() > 3)
            {
                // Histogram
                std::vector<RFLOAT> histX, histY;
                CPlot2D *plot2D=new CPlot2D("");
                tomogramSet.globalTable.columnHistogram(label,histX,histY,0, plot2D);
                fn_eps = fn_out + "hist_" + EMDL::label2Str(label) + ".eps";
                plot2D->OutputPostScriptPlot(fn_eps);
                all_fn_eps.push_back(fn_eps);
                delete plot2D;
            }
        }
    }

    // Also add all the EPS plots of the refined tilt series alignments
    for (long itomo = 0; itomo < tomogramSet.size(); itomo++)
    {
        std::string tomoname = tomogramSet.getTomogramName(itomo);
        FileName fn_dir = fn_out + "external/" + tomoname + '/';
        FileName fn_eps_tilt = fn_dir + tomoname + "_tilts.eps";
        FileName fn_eps_shift = fn_dir + tomoname + "_shifts.eps";
        all_fn_eps.push_back(fn_eps_tilt);
        all_fn_eps.push_back(fn_eps_shift);
        if (do_aretomo && do_aretomo_ctf)
        {
            FileName fn_eps_defocus = fn_dir + tomoname + "_defocus.eps";
            FileName fn_eps_ctffom = fn_dir + tomoname + "_ctfscore.eps";
            all_fn_eps.push_back(fn_eps_defocus);
            all_fn_eps.push_back(fn_eps_ctffom);
        }
    }

    joinMultipleEPSIntoSinglePDF(fn_out + "logfile.pdf", all_fn_eps);

    if (verb > 0 )
    {
        std::cout << " Done! Written out: " << fn_out << "logfile.pdf" << std::endl;
    }

}

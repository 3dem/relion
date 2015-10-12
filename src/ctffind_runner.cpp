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
#include "src/ctffind_runner.h"

void CtffindRunner::read(int argc, char **argv, int rank)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	int ctf_section = parser.addSection("CTF estimation");
	fn_in = parser.getOption("--i", "STAR file with all input micrographs, or a unix wildcard to all micrograph files, e.g. \"mics/*.mrc\"");
	fn_out = parser.getOption("--o", "Name for the STAR file with CTF params for each micrograph", "micrographs_ctf.star");
	do_only_join_results = parser.checkOption("--only_make_star", "Don't run CTFFIND, only join the logfile results in a STAR file");
	continue_old = parser.checkOption("--only_do_unfinished", "Only run CTFFIND for those micrographs for which there is not yet a logfile with Final values.");

	// Use a smaller squared part of the micrograph to estimate CTF (e.g. to avoid film labels...)
	ctf_win =  textToInteger(parser.getOption("--ctfWin", "Size (in pixels) of a centered, squared window to use for CTF-estimation", "-1"));

	fn_ctffind_exe = parser.getOption("--ctffind_exe","Location of ctffind executable (or through RLN_CTFFIND_EXECUTABLE environment variable)","");

	// First parameter line in CTFFIND
	Cs = textToFloat(parser.getOption("--CS", "Spherical Aberration (mm) ","2.0"));
	Voltage = textToFloat(parser.getOption("--HT", "Voltage (kV)","300"));
	AmplitudeConstrast = textToFloat(parser.getOption("--AmpCnst", "Amplitude constrast", "0.1"));
	Magnification = textToFloat(parser.getOption("--XMAG", "Magnification", "60000"));
	PixelSize = textToFloat(parser.getOption("--DStep", "Detector pixel size (um)", "14"));
	// Second parameter line in CTFFIND
	box_size = textToFloat(parser.getOption("--Box", "Size of the boxes to calculate FFTs", "512"));
	resol_min = textToFloat(parser.getOption("--ResMin", "Minimum resolution (in A) to include in calculations", "100"));
	resol_max = textToFloat(parser.getOption("--ResMax", "Maximum resolution (in A) to include in calculations", "7"));
	min_defocus = textToFloat(parser.getOption("--dFMin", "Minimum defocus value (in A) to search", "10000"));
	max_defocus = textToFloat(parser.getOption("--dFMax", "Maximum defocus value (in A) to search", "50000"));
	step_defocus = textToFloat(parser.getOption("--FStep", "defocus step size (in A) for search", "250"));
	amount_astigmatism  = textToFloat(parser.getOption("--dAst", "amount of astigmatism (in A)", "0"));

	// Initialise verb for non-parallel execution
	verb = 1;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void CtffindRunner::usage()
{
	parser.writeUsage(std::cerr);
}

void CtffindRunner::initialise()
{

	// Get the CTFFIND executable
	if (fn_ctffind_exe == "")
	{
		char * penv;
		penv = getenv ("RLN_CTFFIND_EXECUTABLE");
		if (penv!=NULL)
			fn_ctffind_exe = (std::string)penv;
	}

	// Set up which micrographs to estimate CTFs from
	if (fn_in.isStarFile())
	{
		MetaDataTable MDin;
		MDin.read(fn_in);
		fn_micrographs.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			FileName fn_mic;
			MDin.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			fn_micrographs.push_back(fn_mic);
		}
	}
	else
	{
		fn_in.globFiles(fn_micrographs);
	}

	// If we're continuing an old run, see which micrographs have not been finished yet...
	if (continue_old)
	{
		std::vector<FileName> fns_todo;
		for (long int imic = 0; imic < fn_micrographs.size(); imic++)
		{
			FileName fn_microot = fn_micrographs[imic].without(".mrc");
			RFLOAT defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep;
			if (!getCtffindResults(fn_microot, defU, defV, defAng, CC,
					HT, CS, AmpCnst, XMAG, DStep, false)) // false: dont die if not found Final values
				fns_todo.push_back(fn_micrographs[imic]);
		}

		fn_micrographs = fns_todo;

	}


	if (verb > 0)
	{
		std::cout << " Using CTFFINDs executable in: " << fn_ctffind_exe << std::endl;
		std::cout << " to estimate CTF parameters for the following micrographs: " << std::endl;
		if (continue_old)
			std::cout << " (skipping all micrographs for which a logfile with Final values already exists " << std::endl;
		for(unsigned  int  i = 0; i < fn_micrographs.size(); ++i)
			std::cout << "  * " << fn_micrographs[i] << std::endl;
	}
}

void CtffindRunner::run()
{

	if (!do_only_join_results)
	{
		int barstep;
		if (verb > 0)
		{
			std::cout << " Estimating CTF parameters using Niko Grigorieff's CTFFIND ..." << std::endl;
			init_progress_bar(fn_micrographs.size());
			barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
		}

		for (long int imic = 0; imic < fn_micrographs.size(); imic++)
		{
			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);

			executeCtffind(fn_micrographs[imic]);
		}

		if (verb > 0)
			progress_bar(fn_micrographs.size());
	}

	joinCtffindResults();

}

void CtffindRunner::joinCtffindResults()
{
	MetaDataTable MDctf;
	for (long int imic = 0; imic < fn_micrographs.size(); imic++)
    {
		FileName fn_microot = fn_micrographs[imic].without(".mrc");
		RFLOAT defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep;
		bool has_this_ctf = getCtffindResults(fn_microot, defU, defV, defAng, CC,
				HT, CS, AmpCnst, XMAG, DStep);

		if (!has_this_ctf)
			REPORT_ERROR("CtffindRunner::joinCtffindResults ERROR; cannot get CTF values for");

		FileName fn_ctf = fn_microot + ".ctf:mrc";
		MDctf.addObject();
		MDctf.setValue(EMDL_MICROGRAPH_NAME, fn_micrographs[imic]);
	    MDctf.setValue(EMDL_CTF_IMAGE, fn_ctf);
		MDctf.setValue(EMDL_CTF_DEFOCUSU, defU);
	    MDctf.setValue(EMDL_CTF_DEFOCUSV, defV);
	    MDctf.setValue(EMDL_CTF_DEFOCUS_ANGLE, defAng);
	    MDctf.setValue(EMDL_CTF_VOLTAGE, HT);
	    MDctf.setValue(EMDL_CTF_CS, CS);
	    MDctf.setValue(EMDL_CTF_Q0, AmpCnst);
	    MDctf.setValue(EMDL_CTF_MAGNIFICATION, XMAG);
	    MDctf.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, DStep);
	    MDctf.setValue(EMDL_CTF_FOM, CC);
    }
	MDctf.write(fn_out);
}

void CtffindRunner::executeCtffind(FileName fn_mic)
{

	FileName fn_root = fn_mic.withoutExtension();
	FileName fn_script = fn_root + "_ctffind3.com";
	FileName fn_log = fn_root + "_ctffind3.log";
	FileName fn_ctf = fn_root + ".ctf";
    FileName fn_mic_win;

	std::ofstream  fh;
	fh.open((fn_script).c_str(), std::ios::out);
	if (!fh)
	 REPORT_ERROR( (std::string)"CtffindRunner::execute_ctffind cannot create file: " + fn_script);

        // If given, then put a square window of ctf_win on the micrograph for CTF estimation
        if (ctf_win > 0)
        {
            // Window micrograph to a smaller, squared sub-micrograph to estimate CTF on
            fn_mic_win = fn_root + "_win.mrc";
            // Read in micrograph, window and write out again
            Image<RFLOAT> I;
            I.read(fn_mic);
            I().setXmippOrigin();
            I().window(FIRST_XMIPP_INDEX(ctf_win), FIRST_XMIPP_INDEX(ctf_win), LAST_XMIPP_INDEX(ctf_win), LAST_XMIPP_INDEX(ctf_win));
            // Calculate mean, stddev, min and max
            RFLOAT avg, stddev, minval, maxval;
            I().computeStats(avg, stddev, minval, maxval);
            I.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
            I.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
            I.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
            I.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);
            I.write(fn_mic_win);
        }
        else
            fn_mic_win = fn_mic;

	// Write script to run ctffind
	fh << "#!/usr/bin/env csh"<<std::endl;
	fh << fn_ctffind_exe << " > " << fn_log << " << EOF"<<std::endl;
	fh << fn_mic_win << std::endl;
	fh << fn_ctf << std::endl;
	// CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]
	fh << Cs << ", " << Voltage << ", " << AmplitudeConstrast << ", " << Magnification << ", " << PixelSize<< std::endl;
	// Box, ResMin[A], ResMax[A], dFMin[A], dFMax[A], FStep[A], dAst[A]
	fh << box_size << ", " << resol_min << ", " << resol_max << ", " << min_defocus << ", " << max_defocus << ", " << step_defocus << ", " << amount_astigmatism << std::endl;
	fh <<"EOF"<<std::endl;
	fh.close();

	// Execute ctffind
	if (!system(NULL))
	 REPORT_ERROR("There is a problem with the system call to run ctffind");
	FileName fn_cmnd = "csh "+ fn_script;
	system ( fn_cmnd.c_str() );

	// Remove windowed file again
	if (ctf_win > 0)
	{
		if( remove( fn_mic_win.c_str() ) != 0 )
			REPORT_ERROR( "Error deleting windowed micrograph file..." );
	}

}

bool getCtffindResults(FileName fn_microot, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
		RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep, bool die_if_not_found)
{

	FileName fn_log = fn_microot + "_ctffind3.log";
	std::ifstream in(fn_log.data(), std::ios_base::in);
    if (in.fail())
    	return false;

    // Start reading the ifstream at the top
    in.seekg(0);

    // Proceed until the next "Final values" statement
    // The loop statement may be necessary for data blocks that have a list AND a table inside them
    bool Final_is_found = false;
    bool Cs_is_found = false;
    std::string line;
    std::vector<std::string> words;
    while (getline(in, line, '\n'))
    {
        // Find data_ lines

    	 if (line.find("CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]") != std::string::npos)
    	 {
    		 Cs_is_found = true;
    		 getline(in, line, '\n');
    		 tokenize(line, words);
    		 if (words.size() < 5)
    			 REPORT_ERROR("ERROR: Unexpected number of words on data line with CS[mm], HT[kV], etc in " + fn_log);
    		 CS = textToFloat(words[0]);
    		 HT = textToFloat(words[1]);
    		 AmpCnst = textToFloat(words[2]);
    		 XMAG = textToFloat(words[3]);
    		 DStep = textToFloat(words[4]);
    	 }

    	if (line.find("Final Values") != std::string::npos)
        {
        	Final_is_found = true;
            tokenize(line, words);
            if (words.size() < 6)
            	REPORT_ERROR("ERROR: Unexpected number of words on Final values line in " + fn_log);
            defU = textToFloat(words[0]);
            defV = textToFloat(words[1]);
            defAng = textToFloat(words[2]);
            CC = textToFloat(words[3]);
        }
    }
    if (!Cs_is_found && die_if_not_found)
    	REPORT_ERROR("ERROR: cannot find line with Cs[mm], HT[kV], etc values in " + fn_log);
    if (!Final_is_found && die_if_not_found)
    	REPORT_ERROR("ERROR: cannot find line with Final values in " + fn_log);

    in.close();

    return Final_is_found;

}



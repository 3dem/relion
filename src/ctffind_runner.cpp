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
	fn_out = parser.getOption("--o", "Directory, where all output files will be stored", "CtfEstimate/");
	do_only_join_results = parser.checkOption("--only_make_star", "Don't estimate any CTFs, only join all logfile results in a STAR file");
	continue_old = parser.checkOption("--only_do_unfinished", "Only estimate CTFs for those micrographs for which there is not yet a logfile with Final values.");
	// Use a smaller squared part of the micrograph to estimate CTF (e.g. to avoid film labels...)
	ctf_win =  textToInteger(parser.getOption("--ctfWin", "Size (in pixels) of a centered, squared window to use for CTF-estimation", "-1"));

	int mic_section = parser.addSection("Microscopy parameters");
	// First parameter line in CTFFIND
	Cs = textToFloat(parser.getOption("--CS", "Spherical Aberration (mm) ","2.0"));
	Voltage = textToFloat(parser.getOption("--HT", "Voltage (kV)","300"));
	AmplitudeConstrast = textToFloat(parser.getOption("--AmpCnst", "Amplitude constrast", "0.1"));
	Magnification = textToFloat(parser.getOption("--XMAG", "Magnification", "60000"));
	PixelSize = textToFloat(parser.getOption("--DStep", "Detector pixel size (um)", "14"));

	int ctffind_section = parser.addSection("CTFFIND parameters");

	// Second parameter line in CTFFIND
	fn_ctffind_exe = parser.getOption("--ctffind_exe","Location of ctffind executable (or through RELION_CTFFIND_EXECUTABLE environment variable)","");
	box_size = textToFloat(parser.getOption("--Box", "Size of the boxes to calculate FFTs", "512"));
	resol_min = textToFloat(parser.getOption("--ResMin", "Minimum resolution (in A) to include in calculations", "100"));
	resol_max = textToFloat(parser.getOption("--ResMax", "Maximum resolution (in A) to include in calculations", "7"));
	min_defocus = textToFloat(parser.getOption("--dFMin", "Minimum defocus value (in A) to search", "10000"));
	max_defocus = textToFloat(parser.getOption("--dFMax", "Maximum defocus value (in A) to search", "50000"));
	step_defocus = textToFloat(parser.getOption("--FStep", "defocus step size (in A) for search", "250"));
	amount_astigmatism  = textToFloat(parser.getOption("--dAst", "amount of astigmatism (in A)", "0"));

	int gctf_section = parser.addSection("Gctf parameters");
	do_use_gctf = parser.checkOption("--use_gctf", "Use Gctf instead of CTFFIND to estimate the CTF parameters");
	fn_gctf_exe = parser.getOption("--gctf_exe","Location of Gctf executable (or through RELION_GCTF_EXECUTABLE environment variable)","");
	angpix = textToFloat(parser.getOption("--angpix", "Magnified pixel size in Angstroms", "1."));
	do_ignore_ctffind_params = parser.checkOption("--ignore_ctffind_params", "Use Gctf default parameters instead of CTFFIND parameters");
	do_EPA = parser.checkOption("--EPA", "Use equi-phase averaging to calculate Thon rinds in Gctf");
	do_validation = parser.checkOption("--do_validation", "Use validation inside Gctf to analyse quality of the fit?");

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
		penv = getenv ("RELION_CTFFIND_EXECUTABLE");
		if (penv!=NULL)
			fn_ctffind_exe = (std::string)penv;
	}
	// Get the GCTF executable
	if (do_use_gctf && fn_gctf_exe == "")
	{
		char * penv;
		penv = getenv ("RELION_GCTF_EXECUTABLE");
		if (penv!=NULL)
			fn_gctf_exe = (std::string)penv;
	}

	if (do_use_gctf && ctf_win>0)
		REPORT_ERROR("CtffindRunner::initialise ERROR: Running Gctf together with --ctfWin is not implemented, please use CTFFIND instead.");

	// Make sure fn_out ends with a slash
	if (fn_out[fn_out.length()-1] != '/')
		fn_out += "/";

	// Set up which micrographs to estimate CTFs from
	if (fn_in.isStarFile())
	{
		MetaDataTable MDin;
		MDin.read(fn_in);
		fn_micrographs_all.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			FileName fn_mic;
			MDin.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			fn_micrographs_all.push_back(fn_mic);
		}
	}
	else
	{
		fn_in.globFiles(fn_micrographs_all);
	}

	// If we're continuing an old run, see which micrographs have not been finished yet...
	if (continue_old)
	{
		fn_micrographs.clear();
		for (long int imic = 0; imic < fn_micrographs_all.size(); imic++)
		{
			FileName fn_microot = fn_micrographs_all[imic].without(".mrc");
			RFLOAT defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep, maxres=-1., bfac = -1., valscore = -1.;
			if (!getCtffindResults(fn_microot, defU, defV, defAng, CC,
					HT, CS, AmpCnst, XMAG, DStep, maxres, bfac, valscore, false)) // false: dont die if not found Final values
				fn_micrographs.push_back(fn_micrographs_all[imic]);
		}
	}
	else
		fn_micrographs = fn_micrographs_all;

	// Make symbolic links of the input micrographs in the output directory because ctffind and gctf write output files alongside the input micropgraph
    char temp [180];
    getcwd(temp, 180);
    currdir = std::string(temp);
    // Make sure fn_out ends with a slash
	if (currdir[currdir.length()-1] != '/')
		currdir += "/";
	FileName prevdir="";
	for (size_t i = 0; i < fn_micrographs.size(); i++)
	{
		// Remove the UNIQDATE part of the filename if present
		FileName output = getOutputFileWithNewUniqueDate(fn_micrographs[i], fn_out);
		// Create output directory if neccesary
		FileName newdir = output.beforeLastOf("/");
		if (newdir != prevdir)
		{
			std::string command = " mkdir -p " + newdir;
			int res = system(command.c_str());
		}
		symlink((currdir+fn_micrographs[i]).c_str(), output.c_str());
	}

	if (do_use_gctf && fn_micrographs.size()>0)
	{
		// Find the dimensions of the first micrograph, to later on ensure all micrographs are the same size
		Image<double> Itmp;
		Itmp.read(fn_micrographs[0], false); // false means only read header!
		xdim = XSIZE(Itmp());
		ydim = YSIZE(Itmp());
	}

	if (verb > 0)
	{
		if (do_use_gctf)
			std::cout << " Using Gctf executable in: " << fn_gctf_exe << std::endl;
		else
			std::cout << " Using CTFFIND executable in: " << fn_ctffind_exe << std::endl;
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
			if (do_use_gctf)
				std::cout << " Estimating CTF parameters using Kai Zhang's Gctf ..." << std::endl;
			else
				std::cout << " Estimating CTF parameters using Niko Grigorieff's CTFFIND ..." << std::endl;
			init_progress_bar(fn_micrographs.size());
			barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
		}

		std::vector<std::string> allmicnames;
		for (long int imic = 0; imic < fn_micrographs.size(); imic++)
		{

			if (do_use_gctf)
			{
				//addToGctfJobList(imic, allmicnames);
				executeGctf(imic, allmicnames, imic+1==fn_micrographs.size());
			}
			else
			{
				executeCtffind(imic);
			}

			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);
		}

		//if (do_use_gctf)
		//	executeGctf(allmicnames);

		if (verb > 0)
			progress_bar(fn_micrographs.size());
	}

	joinCtffindResults();

}

void CtffindRunner::joinCtffindResults()
{

	MetaDataTable MDctf;
	for (long int imic = 0; imic < fn_micrographs_all.size(); imic++)
    {
		FileName fn_microot = fn_micrographs_all[imic].without(".mrc");
		RFLOAT defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep, maxres=-1., bfac = -1., valscore = -1.;
		bool has_this_ctf = getCtffindResults(fn_microot, defU, defV, defAng, CC,
				HT, CS, AmpCnst, XMAG, DStep, maxres, bfac, valscore);

		if (!has_this_ctf)
			REPORT_ERROR("CtffindRunner::joinCtffindResults ERROR; cannot get CTF values for " + fn_micrographs_all[imic] );

		FileName fn_ctf = fn_out + fn_microot + ".ctf:mrc";
		MDctf.addObject();
		MDctf.setValue(EMDL_MICROGRAPH_NAME, fn_micrographs_all[imic]);
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
	    if (maxres > 0.)
	    	MDctf.setValue(EMDL_CTF_MAXRES, maxres);
	    if (bfac > 0.)
	    	MDctf.setValue(EMDL_CTF_BFACTOR, bfac);
	    if (valscore > 0.)
	    	MDctf.setValue(EMDL_CTF_VALIDATIONSCORE, valscore);

    }
	MDctf.write(fn_out+"micrographs_ctf.star");
	std::cout << " Done! Written out: " << fn_out <<  "micrographs_ctf.star" << std::endl;

	if (do_use_gctf)
	{
		FileName fn_gctf_junk = "micrographs_all_gctf";
		if (exists(fn_gctf_junk))
			remove(fn_gctf_junk.c_str());
		fn_gctf_junk = "extra_micrographs_all_gctf";
		if (exists(fn_gctf_junk))
			remove(fn_gctf_junk.c_str());
	}


}

/*
void CtffindRunner::addToGctfJobList(long int i, std::vector<std::string>  &allmicnames)
{

	Image<double> Itmp;
	FileName outputfile = getOutputFileWithNewUniqueDate(fn_micrographs[i], fn_out);
	Itmp.read(outputfile, false); // false means only read header!
	if (XSIZE(Itmp()) != xdim || YSIZE(Itmp()) != ydim)
		REPORT_ERROR("CtffindRunner::executeGctf ERROR: Micrographs do not all have the same size! " + fn_micrographs[i] + " is different from the first micrograph!");
	if (ZSIZE(Itmp()) > 1 || NSIZE(Itmp()) > 1)
		REPORT_ERROR("CtffindRunner::executeGctf ERROR: No movies or volumes allowed for " + fn_micrographs[i]);

	allmicnames.push_back(outputfile);

}


void CtffindRunner::executeGctf(std::vector<std::string> &allmicnames)
{

	std::string command = fn_gctf_exe;
	//command +=  " --ctfstar " + fn_out + "tt_micrographs_ctf.star";
	command +=  " --apix " + floatToString(angpix);
	command +=  " --cs " + floatToString(Cs);
	command +=  " --kV " + floatToString(Voltage);
	command +=  " --ac " + floatToString(AmplitudeConstrast);
	if (!do_ignore_ctffind_params)
	{
		command += " --boxsize " + floatToString(box_size);
		command += " --resL " + floatToString(resol_min);
		command += " --resH " + floatToString(resol_max);
		command += " --defL " + floatToString(min_defocus);
		command += " --defH " + floatToString(max_defocus);
		command += " --defS " + floatToString(step_defocus);
		command += " --astm " + floatToString(amount_astigmatism);
	}

	if (do_EPA)
		command += " --do_EPA ";

	if (do_validation)
		command += " --do_validation ";

	for (size_t i = 0; i<allmicnames.size(); i++)
		command += allmicnames[i];

	// Redirect all gctf output
	command += " >& " + fn_out + "gctf.out ";

	int res = system(command.c_str());

	// Cleanup all the symbolic links again
	for (size_t i = 0; i < fn_micrographs.size(); i++)
	{
		FileName output = getOutputFileWithNewUniqueDate(fn_micrographs[i], fn_out);
		remove(output.c_str());
	}


}
*/

void CtffindRunner::executeGctf(long int imic, std::vector<std::string> &allmicnames, bool is_last, int rank)
{

	// Always add the new micrograph to the TODO list
	Image<double> Itmp;
	FileName outputfile = getOutputFileWithNewUniqueDate(fn_micrographs[imic], fn_out);
	Itmp.read(outputfile, false); // false means only read header!
	if (XSIZE(Itmp()) != xdim || YSIZE(Itmp()) != ydim)
		REPORT_ERROR("CtffindRunner::executeGctf ERROR: Micrographs do not all have the same size! " + fn_micrographs[imic] + " is different from the first micrograph!");
	if (ZSIZE(Itmp()) > 1 || NSIZE(Itmp()) > 1)
		REPORT_ERROR("CtffindRunner::executeGctf ERROR: No movies or volumes allowed for " + fn_micrographs[imic]);

	allmicnames.push_back(outputfile);

	// Execute Gctf every 20 images, and always for the last one
	if (imic+1%20 || is_last)
	{
		std::string command = fn_gctf_exe;
		//command +=  " --ctfstar " + fn_out + "tt_micrographs_ctf.star";
		command +=  " --apix " + floatToString(angpix);
		command +=  " --cs " + floatToString(Cs);
		command +=  " --kV " + floatToString(Voltage);
		command +=  " --ac " + floatToString(AmplitudeConstrast);
		if (!do_ignore_ctffind_params)
		{
			command += " --boxsize " + floatToString(box_size);
			command += " --resL " + floatToString(resol_min);
			command += " --resH " + floatToString(resol_max);
			command += " --defL " + floatToString(min_defocus);
			command += " --defH " + floatToString(max_defocus);
			command += " --defS " + floatToString(step_defocus);
			command += " --astm " + floatToString(amount_astigmatism);
		}

		if (do_EPA)
			command += " --do_EPA ";

		if (do_validation)
			command += " --do_validation ";

		for (size_t i = 0; i<allmicnames.size(); i++)
			command += allmicnames[i];

		// TODO: better control over which GPU to use. For now, gid = rank!
		command += " --gid " + integerToString(rank);

		// Redirect all gctf output
		command += " >> " + fn_out + "gctf" + integerToString(rank)+".out  2>> " + fn_out + "gctf" + integerToString(rank)+".err";
		int res = system(command.c_str());

		// Cleanup all the symbolic links again
		for (size_t i = 0; i < allmicnames.size(); i++)
			remove(allmicnames[i].c_str());

		// Re-set the allmicnames vector
		allmicnames.clear();
	}


}

void CtffindRunner::executeCtffind(long int imic)
{

	FileName fn_mic = getOutputFileWithNewUniqueDate(fn_micrographs[imic], fn_out);
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
	int res = system ( fn_cmnd.c_str() );

	// Remove windowed file again
	if (ctf_win > 0)
	{
		if( remove( fn_mic_win.c_str() ) != 0 )
			REPORT_ERROR( "Error deleting windowed micrograph file..." );
	}

}

bool CtffindRunner::getCtffindResults(FileName fn_microot, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
		RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep,
		RFLOAT &maxres, RFLOAT &bfac, RFLOAT &valscore, bool die_if_not_found)
{

	FileName fn_log = fn_out + fn_microot + "_ctffind3.log";
	if (do_use_gctf && !exists(fn_log)) // also test _gctf.log file
		fn_log = fn_out + fn_microot + "_gctf.log";

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

    	if (do_use_gctf && line.find("Resolution limit estimated by EPA:") != std::string::npos)
    	{
            tokenize(line, words);
            maxres = textToFloat(words[words.size()-1]);
    	}

    	if (do_use_gctf && line.find("Estimated Bfactor:") != std::string::npos)
    	{
            tokenize(line, words);
             if (words.size() < 4)
             	REPORT_ERROR("ERROR: Unexpected number of words on Resolution limit line in " + fn_log);
             bfac = textToFloat(words[3]);
    	}

    	if (do_use_gctf && line.find("OVERALL_VALIDATION_SCORE:") != std::string::npos)
    	{
            tokenize(line, words);
            valscore = textToFloat(words[words.size()-1]);
    	}
    }

    if (!Cs_is_found && die_if_not_found)
    	REPORT_ERROR("ERROR: cannot find line with Cs[mm], HT[kV], etc values in " + fn_log);
    if (!Final_is_found && die_if_not_found)
    	REPORT_ERROR("ERROR: cannot find line with Final values in " + fn_log);

    in.close();

    return Final_is_found;

}



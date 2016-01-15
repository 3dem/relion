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
#include "src/motioncorr_runner.h"

void MotioncorrRunner::read(int argc, char **argv, int rank)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "STAR file with all input micrographs, or a Linux wildcard with all micrographs to operate on");
	fn_out = parser.getOption("--o", "Name for the output directory", "MotionCorr");
	fn_movie = parser.getOption("--movie", "Rootname to identify movies", "movie");
	continue_old = parser.checkOption("--only_do_unfinished", "Only run MOTIONCORR for those micrographs for which there is not yet an output micrograph.");
	do_save_movies  = parser.checkOption("--save_movies", "Also save the motion-corrected movies.");

	// Use a smaller squared part of the micrograph to estimate CTF (e.g. to avoid film labels...)
	bin_factor =  textToInteger(parser.getOption("--bin_factor", "Binning factor (integer) for scaling inside MOTIONCORR", "1"));
	first_frame_ali =  textToInteger(parser.getOption("--first_frame_ali", "First movie frame used in alignment (start at 1)", "1"));
	last_frame_ali =  textToInteger(parser.getOption("--last_frame_ali", "Last movie frame used in alignment (0: use all)", "0"));
	first_frame_sum =  textToInteger(parser.getOption("--first_frame_sum", "First movie frame used in output sum (start at 1)", "1"));
	last_frame_sum =  textToInteger(parser.getOption("--last_frame_sum", "Last movie frame used in output sum (0: use all)", "0"));
	fn_other_args = parser.getOption("--other_motioncorr_args", "Additional arguments to MOTIONCORR", "");

	fn_motioncorr_exe = parser.getOption("--motioncorr_exe","Location of MOTIONCORR executable (or through RELION_MOTIONCORR_EXECUTABLE environment variable)","");

	// Initialise verb for non-parallel execution
	verb = 1;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void MotioncorrRunner::usage()
{
	parser.writeUsage(std::cerr);
}

void MotioncorrRunner::initialise()
{

	// Get the CTFFIND executable
	if (fn_motioncorr_exe == "")
	{
		char * penv;
		penv = getenv ("RELION_MOTIONCORR_EXECUTABLE");
		if (penv!=NULL)
			fn_motioncorr_exe = (std::string)penv;
	}

	MDout.clear();

	FileName fn_avg, fn_mov;

	// Set up which micrograph movies to run MOTIONCORR on
	if (fn_in.isStarFile())
	{
		MetaDataTable MDin;
		MDin.read(fn_in);
		fn_micrographs.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			FileName fn_mic;
			MDin.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_mic);
			fn_micrographs.push_back(fn_mic);

			// For output STAR file
			getOutputFileNames(fn_mic, fn_avg, fn_mov);
			MDout.addObject(MDin.getObject());
			MDout.setValue(EMDL_MICROGRAPH_NAME, fn_avg);
		}
	}
	else
	{
		fn_in.globFiles(fn_micrographs);

		// For output STAR file
		for (size_t imic = 0; imic < fn_micrographs.size(); imic++)
		{
			getOutputFileNames(fn_micrographs[imic], fn_avg, fn_mov);
			MDout.addObject();
			MDout.setValue(EMDL_MICROGRAPH_NAME, fn_avg);
		}
	}

	// If we're continuing an old run, see which micrographs have not been finished yet...
	if (continue_old)
	{
		std::vector<FileName> fns_todo;
		for (long int imic = 0; imic < fn_micrographs.size(); imic++)
		{
			FileName fn_avg, fn_mov;
			getOutputFileNames(fn_micrographs[imic], fn_avg, fn_mov);
			if (!exists(fn_avg) || (do_save_movies && !exists(fn_mov)) )
				fns_todo.push_back(fn_micrographs[imic]);
		}
		fn_micrographs = fns_todo;
	}

	// Make sure fn_out ends with a slash
	if (fn_out[fn_out.length()-1] != '/')
		fn_out += "/";

	// Make all output directories if necessary
	FileName prevdir="";
	for (size_t i = 0; i < fn_micrographs.size(); i++)
	{
		FileName newdir = fn_micrographs[i].beforeLastOf("/");
		if (newdir != prevdir)
		{
			std::string command = " mkdir -p " + fn_out + newdir;
			int res = system(command.c_str());
		}
	}

	// Motioncorr starts counting frames at 0:
	first_frame_ali -= 1;
	first_frame_sum -= 1;
	if (last_frame_ali != 0)
		last_frame_ali -= 1;
	if (last_frame_sum != 0)
		last_frame_sum -= 1;

	if (verb > 0)
	{
		std::cout << " Using MOTIONCORR executable in: " << fn_motioncorr_exe << std::endl;
		std::cout << " to correct beam-induced motion for the following micrographs: " << std::endl;
		if (continue_old)
			std::cout << " (skipping all micrographs for which a corrected movie already exists) " << std::endl;
		for(unsigned  int  i = 0; i < fn_micrographs.size(); ++i)
			std::cout << "  * " << fn_micrographs[i] << std::endl;
	}

}

void MotioncorrRunner::getOutputFileNames(FileName fn_mic, FileName &fn_avg, FileName &fn_mov)
{
	// If there are any dots in the filename, replace them by underscores
	FileName fn_root = fn_mic.withoutExtension();
	// If fn_root already contains "_movie", then remove that from fn_root
	fn_root = fn_root.without("_"+fn_movie);

	size_t pos = 0;
	while (true)
	{
		pos = fn_root.find(".");
		if (pos == std::string::npos)
			break;
		fn_root.replace(pos, 1, "_");
	}

	fn_avg = fn_out + fn_root + ".mrc";
	fn_mov = fn_out + fn_root + "_" + fn_movie + ".mrcs";
}


void MotioncorrRunner::run()
{

	int barstep;
	if (verb > 0)
	{
		std::cout << " Correcting beam-induced motions using UCSF's MOTIONCORR ..." << std::endl;
		init_progress_bar(fn_micrographs.size());
		barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
	}

	for (long int imic = 0; imic < fn_micrographs.size(); imic++)
	{
		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		executeMotioncorr(fn_micrographs[imic]);
	}

	if (verb > 0)
		progress_bar(fn_micrographs.size());

	// Write out STAR file at the end
	MDout.write(fn_out + "/corrected_micrographs.star");

}


void MotioncorrRunner::executeMotioncorr(FileName fn_mic)
{

	FileName fn_avg, fn_mov;
	getOutputFileNames(fn_mic, fn_avg, fn_mov);


	FileName fn_out = fn_avg.withoutExtension() + ".out";
	FileName fn_log = fn_avg.withoutExtension() + ".log";
	FileName fn_err = fn_avg.withoutExtension() + ".err";
	FileName fn_cmd = fn_avg.withoutExtension() + ".com";

	std::string command = fn_motioncorr_exe + " ";

	command += fn_mic + " -fcs " + fn_avg;
	command += " -flg " + fn_log;
	command += " -nst " + integerToString(first_frame_ali) + " -nss " + integerToString(first_frame_sum);
	command += " -ned " + integerToString(last_frame_ali) + " -nes " + integerToString(last_frame_sum);

	if (do_save_movies)
		command += " -dsp 0 -ssc 1 -fct " + fn_mov;

	if (bin_factor > 1)
		command += " -bin " + integerToString(bin_factor);

	if (fn_other_args.length() > 0)
		command += " " + fn_other_args;


	command += " >> " + fn_out + " 2>> " + fn_err;

	// Save the command that was executed
	std::ofstream fh;
	fh.open(fn_cmd.c_str(), std::ios::out);
	fh << command << std::endl;
	fh.close();

	int res = system(command.c_str());

}

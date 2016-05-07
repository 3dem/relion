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
#ifdef CUDA
#include "src/gpu_utils/cuda_mem_utils.h"
#endif

void MotioncorrRunner::read(int argc, char **argv, int rank)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "STAR file with all input micrographs, or a Linux wildcard with all micrographs to operate on");
	fn_out = parser.getOption("--o", "Name for the output directory", "MotionCorr");
	fn_movie = parser.getOption("--movie", "Rootname to identify movies", "movie");
	continue_old = parser.checkOption("--only_do_unfinished", "Only run mottion correctiob for those micrographs for which there is not yet an output micrograph.");
	do_save_movies  = parser.checkOption("--save_movies", "Also save the motion-corrected movies.");
	gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread, e.g 0:1:2:3", "");

	// Use a smaller squared part of the micrograph to estimate CTF (e.g. to avoid film labels...)
	int motioncorr_section = parser.addSection("MOTIONCORR options");
	bin_factor =  textToInteger(parser.getOption("--bin_factor", "Binning factor (integer) for scaling inside MOTIONCORR", "1"));
	bfactor =  textToFloat(parser.getOption("--bfactor", "B-factor (in pix^2) that will be used inside MOTIONCORR", "150"));
	first_frame_ali =  textToInteger(parser.getOption("--first_frame_ali", "First movie frame used in alignment (start at 1)", "1"));
	last_frame_ali =  textToInteger(parser.getOption("--last_frame_ali", "Last movie frame used in alignment (0: use all)", "0"));
	first_frame_sum =  textToInteger(parser.getOption("--first_frame_sum", "First movie frame used in output sum (start at 1)", "1"));
	last_frame_sum =  textToInteger(parser.getOption("--last_frame_sum", "Last movie frame used in output sum (0: use all)", "0"));
	fn_other_motioncorr_args = parser.getOption("--other_motioncorr_args", "Additional arguments to MOTIONCORR", "");
	fn_motioncorr_exe = parser.getOption("--motioncorr_exe","Location of MOTIONCORR executable (or through RELION_MOTIONCORR_EXECUTABLE environment variable)","");

	int unblur_section = parser.addSection("UNBLUR/SUMMOVIE options");
	do_unblur = parser.checkOption("--use_unblur", "Use Niko Grigorieff's UNBLUR instead of MOTIONCORR.");
	fn_unblur_exe = parser.getOption("--unblur_exe","Location of UNBLUR (v1.0.2) executable (or through RELION_UNBLUR_EXECUTABLE environment variable)","");
	nr_threads = textToInteger(parser.getOption("--j","Number of threads in the unblur executable","1"));

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

	if (do_unblur)
	{
		// Get the CTFFIND executable
		if (fn_unblur_exe == "")
		{
			char * penv;
			penv = getenv ("RELION_UNBLUR_EXECUTABLE");
			if (penv!=NULL)
				fn_unblur_exe = (std::string)penv;
		}
	}
	else
	{
		// Get the CTFFIND executable
		if (fn_motioncorr_exe == "")
		{
			char * penv;
			penv = getenv ("RELION_MOTIONCORR_EXECUTABLE");
			if (penv!=NULL)
				fn_motioncorr_exe = (std::string)penv;
		}
	}

	MDavg.clear();
	MDmov.clear();

#ifdef CUDA
	if (gpu_ids.length() > 0)
		untangleDeviceIDs(gpu_ids, allThreadIDs);
	else if (!do_unblur && verb>0)
		std::cout << "gpu-ids not specified, threads will automatically be mapped to devices (incrementally)."<< std::endl;

	HANDLE_ERROR(cudaGetDeviceCount(&devCount));
#endif

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
			MDmov.addObject();
			MDmov.setValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_mov);
			MDavg.addObject();
			MDavg.setValue(EMDL_MICROGRAPH_NAME, fn_avg);
		}
	}
	else
	{
		fn_in.globFiles(fn_micrographs);

		// For output STAR file
		for (size_t imic = 0; imic < fn_micrographs.size(); imic++)
		{
			getOutputFileNames(fn_micrographs[imic], fn_avg, fn_mov);
			MDmov.addObject();
			MDmov.setValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_mov);
			MDavg.addObject();
			MDavg.setValue(EMDL_MICROGRAPH_NAME, fn_avg);
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
		if (do_unblur)
			std::cout << " Using UNBLUR executable in: " << fn_unblur_exe << std::endl;
		else
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
		if (do_unblur)
			std::cout << " Correcting beam-induced motions using Niko Grigorieff's UNBLUR ..." << std::endl;
		else
			std::cout << " Correcting beam-induced motions using UCSF's MOTIONCORR ..." << std::endl;
		init_progress_bar(fn_micrographs.size());
		barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
	}

	for (long int imic = 0; imic < fn_micrographs.size(); imic++)
	{
		std::vector<float> xshifts, yshifts;

		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		if (do_unblur)
			executeUnblur(fn_micrographs[imic], xshifts, yshifts);
		else
			executeMotioncorr(fn_micrographs[imic], xshifts, yshifts);

		plotShifts(fn_micrographs[imic], xshifts, yshifts);
	}

	if (verb > 0)
		progress_bar(fn_micrographs.size());

	// Make a logfile with the shifts in pdf format
	generateLogFilePDF();

	// Write out STAR files at the end
	MDavg.write(fn_out + "corrected_micrographs.star");
	MDmov.write(fn_out + "corrected_micrograph_movies.star");

}


void MotioncorrRunner::executeMotioncorr(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts, int rank)
{


	FileName fn_avg, fn_mov;
	getOutputFileNames(fn_mic, fn_avg, fn_mov);

	FileName fn_out = fn_avg.withoutExtension() + ".out";
	FileName fn_log = fn_avg.withoutExtension() + ".log";
	FileName fn_err = fn_avg.withoutExtension() + ".err";
	FileName fn_cmd = fn_avg.withoutExtension() + ".com";

	for (int ipass = 0; ipass < 3; ipass++)
	{

		std::string command = fn_motioncorr_exe + " ";

		command += fn_mic + " -fcs " + fn_avg;
		command += " -flg " + fn_log;
		command += " -nst " + integerToString(first_frame_ali) + " -nss " + integerToString(first_frame_sum);
		command += " -ned " + integerToString(last_frame_ali) + " -nes " + integerToString(last_frame_sum);
		command += " -bft " + floatToString(bfactor);

		if (do_save_movies)
			command += " -dsp 0 -ssc 1 -fct " + fn_mov;

		if (bin_factor > 1)
			command += " -bin " + integerToString(bin_factor);


		if (fn_other_motioncorr_args.length() > 0)
			command += " " + fn_other_motioncorr_args;

		if ( allThreadIDs.size() == 0)
		{
			// Automated mapping
			command += " -gpu " + integerToString(rank % devCount);
		}
		else
		{
			if (rank >= allThreadIDs.size())
				REPORT_ERROR("ERROR: not enough MPI nodes specified for the GPU IDs.");
			command += " -gpu " + allThreadIDs[rank][0];
		}

		command += " >> " + fn_out + " 2>> " + fn_err;

		// Save the command that was executed
		std::ofstream fh;
		fh.open(fn_cmd.c_str(), std::ios::out);
		fh << command << std::endl;
		fh.close();

		if (system(command.c_str()))
			REPORT_ERROR("ERROR executing: " + command);

		// After motion-correction, check for all-zero average micrographs
		if (exists(fn_avg))
		{
			Image<RFLOAT> Itest;
			Itest.read(fn_avg, false);
			RFLOAT avg, stddev;
			Itest.MDMainHeader.getValue(EMDL_IMAGE_STATS_STDDEV, stddev);
			Itest.MDMainHeader.getValue(EMDL_IMAGE_STATS_AVG, avg);
			if (fabs(stddev) > 0.00001 || fabs(avg) > 0.00001)
			{
				break;
			}
		}
		else if (ipass == 2)
		{
			std::cerr << " WARNING: " << fn_avg << " still did not exist or had zero mean and variance after 3 attempts! " << std::endl;
		}
	}

	// Also analyse the shifts
	getShiftsMotioncorr(fn_log, xshifts, yshifts);

}

void MotioncorrRunner::getShiftsMotioncorr(FileName fn_log, std::vector<float> &xshifts, std::vector<float> &yshifts)
{

	std::ifstream in(fn_log.data(), std::ios_base::in);
	if (in.fail())
		return;

	xshifts.clear();
	yshifts.clear();

    std::string line;

    // Start reading the ifstream at the top
    in.seekg(0);

    // Read throught the shifts file
    int i = 0;
    bool have_found_final = false;
    while (getline(in, line, '\n'))
    {
    	// ignore all commented lines, just read first two lines with data
    	if (line.find("Final shift") != std::string::npos)
    	{
    		have_found_final = true;
    	}
    	else if (have_found_final)
    	{
    		if (line.find("Shift") != std::string::npos)
    		{
        		std::vector<std::string> words;
        		tokenize(line, words);
        		if (words.size() < 7)
        		{
        			std::cerr << " fn_log= " << fn_log << std::endl;
        			REPORT_ERROR("ERROR: unexpected number of words on line from MOTIONCORR logfile: " + line);
        		}
        		xshifts.push_back(textToFloat(words[5]));
    			yshifts.push_back(textToFloat(words[6]));
    		}
    		else
    		{
    			// Stop now
    			break;
    		}
    	}
    }
    in.close();


    if (xshifts.size() != yshifts.size())
    	REPORT_ERROR("ERROR: got an unequal number of x and yshifts from " + fn_log);

}

void MotioncorrRunner::executeUnblur(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts)
{

	FileName fn_avg, fn_mov;
	getOutputFileNames(fn_mic, fn_avg, fn_mov);
	FileName fn_root = fn_avg.withoutExtension();
	FileName fn_log = fn_root + "_unblur.log";
	FileName fn_com = fn_root + "_unblur.com";
	FileName fn_shifts = fn_root + "_shifts.txt";

	FileName fn_tmp_mic;
	// Unblur cannot handle .mrcs extensions
	if (fn_mic.getExtension() == "mrcs")
	{
		fn_tmp_mic = fn_out + fn_mic.withoutExtension() + "_in.mrc";

		// See how many directories deep is the output name
		size_t ndir = std::count(fn_tmp_mic.begin(), fn_tmp_mic.end(), '/');
		FileName fn_link = "";
		for (size_t i =0; i < ndir; i++)
		{
			fn_link += "../";
		}
		// Make the symbolic link relative to the project dir
		fn_link += fn_mic;
		symlink(fn_link.c_str(), fn_tmp_mic.c_str());
	}
	else
	{
		fn_tmp_mic = fn_mic;
	}
	FileName fn_tmp_mov = fn_mov.withoutExtension() + ".mrc";


	std::ofstream  fh;
	fh.open((fn_com).c_str(), std::ios::out);
	if (!fh)
	 REPORT_ERROR( (std::string)"executeUnblur cannot create file: " + fn_com);

	// Write script to run ctffind
	fh << "#!/usr/bin/env csh"<<std::endl;
	fh << "setenv  OMP_NUM_THREADS " << integerToString(nr_threads)<<std::endl;
	fh << fn_unblur_exe << " > " << fn_log << "  << EOF"<<std::endl;
	fh << fn_tmp_mic << std::endl;
	fh << std::endl;
	fh << fn_avg << std::endl;
	fh << fn_shifts << std::endl;
	fh << "1.0" << std::endl;
	fh << "NO" << std::endl; // no dose filtering
	if (do_save_movies)
	{
		fh << "YES" << std::endl; // save movie frames
		fh << fn_tmp_mov << std::endl;
	}
	else
	{
		fh << "NO" << std::endl; // dont set expert options
	}
	fh << "NO" << std::endl; // dont set expert options
	fh <<"EOF"<<std::endl;
	fh.close();

	// Execute unblur
	std::string command = "csh "+ fn_com;
	if (system(command.c_str()))
		REPORT_ERROR("ERROR in executing: " + command);

	// Also analyse the shifts
	getShiftsUnblur(fn_shifts, xshifts, yshifts);

	// Move movie .mrc to new .mrcs filename
	std::rename(fn_tmp_mov.c_str(), fn_mov.c_str());

}

void MotioncorrRunner::getShiftsUnblur(FileName fn_shifts, std::vector<float> &xshifts, std::vector<float> &yshifts)
{

	std::ifstream in(fn_shifts.data(), std::ios_base::in);
	if (in.fail())
		return;

	xshifts.clear();
	yshifts.clear();

	std::string line, token;

    // Start reading the ifstream at the top
    in.seekg(0);

    // Read throught the shifts file
    int i = 0;
    while (getline(in, line, '\n'))
    {
    	// ignore all commented lines, just read first two lines with data
    	if (line[0] != '#')
    	{
    		if (i>1)
    			REPORT_ERROR("ERROR: reading more than 2 data lines from " + fn_shifts);

    		std::vector<std::string> words;
    	    tokenize(line, words);
    		for (int j = 0; j < words.size(); j++)
    		{
    			float sh = textToFloat(words[j]);
    			if (i==0)
    				xshifts.push_back(sh);
    			else if (i==1)
    				yshifts.push_back(sh);
    		}
    		i++;
    	}
    }
    in.close();

    if (xshifts.size() != yshifts.size())
    	REPORT_ERROR("ERROR: got an unequal number of x and yshifts from " + fn_shifts);
}

// Plot the shifts
void MotioncorrRunner::plotShifts(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts)
{

	if (xshifts.size() == 0)
		return;

	FileName fn_eps = fn_out + fn_mic.withoutExtension() + "_shifts.eps";
	CPlot2D *plot2D=new CPlot2D(fn_eps);
 	plot2D->SetXAxisSize(600);
 	plot2D->SetYAxisSize(600);

 	CDataSet dataSet;
	dataSet.SetDrawMarker(false);
	dataSet.SetDatasetColor(0.0,0.0,0.0);
	for (int j = 0; j < xshifts.size(); j++)
	{
		CDataPoint point(xshifts[j], yshifts[j]);
		dataSet.AddDataPoint(point);
	}
	plot2D->AddDataSet(dataSet);

	// Different starting point
	CDataSet dataSetStart;
	dataSetStart.SetDrawMarker(true);
	dataSetStart.SetDatasetColor(1.0,0.0,0.0);
	CDataPoint point2(xshifts[0], yshifts[0]);
	dataSetStart.AddDataPoint(point2);
	plot2D->AddDataSet(dataSetStart);

	plot2D->SetXAxisTitle("X-shift (in pixels)");
	plot2D->SetYAxisTitle("Y-shift (in pixels)");
	plot2D->OutputPostScriptPlot(fn_eps);


}


void MotioncorrRunner::generateLogFilePDF()
{

	if (fn_micrographs.size() > 0)
	{

		std::vector<FileName> fn_eps;


		FileName fn_prev="";
		for (long int i = 0; i < fn_micrographs.size(); i++)
		{
			if (fn_prev != fn_micrographs[i].beforeLastOf("/"))
			{
				fn_prev = fn_micrographs[i].beforeLastOf("/");
				fn_eps.push_back(fn_out + fn_prev+"/*.eps");
			}
		}

		joinMultipleEPSIntoSinglePDF(fn_out + "logfile.pdf ", fn_eps);

	}
}


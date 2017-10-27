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
	continue_old = parser.checkOption("--only_do_unfinished", "Only run motion correction for those micrographs for which there is not yet an output micrograph.");
	do_save_movies  = parser.checkOption("--save_movies", "Also save the motion-corrected movies.");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms","-1"));
	first_frame_ali =  textToInteger(parser.getOption("--first_frame_ali", "First movie frame used in alignment (start at 1)", "1"));
	last_frame_ali =  textToInteger(parser.getOption("--last_frame_ali", "Last movie frame used in alignment (0: use all)", "0"));
	first_frame_sum =  textToInteger(parser.getOption("--first_frame_sum", "First movie frame used in output sum (start at 1)", "1"));
	last_frame_sum =  textToInteger(parser.getOption("--last_frame_sum", "Last movie frame used in output sum (0: use all)", "0"));

	int motioncor2_section = parser.addSection("MOTIONCOR2 options");
	do_motioncor2 = parser.checkOption("--use_motioncor2", "Use Shawn Zheng's MOTIONCOR2 instead of UNBLUR.");
	fn_motioncor2_exe = parser.getOption("--motioncor2_exe","Location of MOTIONCOR2 executable (or through RELION_MOTIONCOR2_EXECUTABLE environment variable)","");
	bin_factor =  textToFloat(parser.getOption("--bin_factor", "Binning factor (can be non-integer)", "1"));
	bfactor =  textToFloat(parser.getOption("--bfactor", "B-factor (in pix^2) that will be used inside MOTIONCOR2", "150"));
	fn_gain_reference = parser.getOption("--gainref","Location of MRC file with the gain reference to be applied","");
	patch_x = textToInteger(parser.getOption("--patch_x", "Patching in X-direction for MOTIONCOR2", "1"));
	patch_y = textToInteger(parser.getOption("--patch_y", "Patching in Y-direction for MOTIONCOR2", "1"));
	group = textToInteger(parser.getOption("--group_frames", "Average together this many frames before calculating the beam-induced shifts", "1"));
	fn_defect = parser.getOption("--defect_file","Location of a MOTIONCOR2-style detector defect file","");
	fn_archive = parser.getOption("--archive","Location of the directory for archiving movies in 4-byte MRC format","");

	fn_other_motioncor2_args = parser.getOption("--other_motioncor2_args", "Additional arguments to MOTIONCOR2", "");
	gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread, e.g 0:1:2:3", "");

	int unblur_section = parser.addSection("UNBLUR/SUMMOVIE options");
	do_unblur = parser.checkOption("--use_unblur", "Use Niko Grigorieff's UNBLUR instead of MOTIONCOR2.");
	fn_unblur_exe = parser.getOption("--unblur_exe","Location of UNBLUR (v1.0.2) executable (or through RELION_UNBLUR_EXECUTABLE environment variable)","");
	fn_summovie_exe = parser.getOption("--summovie_exe","Location of SUMMOVIE(v1.0.2) executable (or through RELION_SUMMOVIE_EXECUTABLE environment variable)","");
	nr_threads = textToInteger(parser.getOption("--j","Number of threads in the unblur executable","1"));

	int doseweight_section = parser.addSection("Dose-weighting options");
	do_dose_weighting = parser.checkOption("--dose_weighting", "Use MOTIONCOR2s or UNBLURs dose-weighting scheme");
	voltage = textToFloat(parser.getOption("--voltage","Voltage (in kV) for dose-weighting inside MOTIONCOR2/UNBLUR","300"));
	dose_per_frame = textToFloat(parser.getOption("--dose_per_frame", "Electron dose (in electrons/A2/frame) for dose-weighting inside MOTIONCOR2/UNBLUR", "0"));
	pre_exposure = textToFloat(parser.getOption("--preexposure", "Pre-exposure (in electrons/A2) for dose-weighting inside UNBLUR", "0"));

	// Initialise verb for non-parallel execution
	verb = 1;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void MotioncorrRunner::usage()
{
	parser.writeUsage(std::cout);
}

void MotioncorrRunner::initialise()
{

	if (do_unblur)
	{
		// Get the UNBLUR executable
		if (fn_unblur_exe == "")
		{
			char * penv;
			penv = getenv ("RELION_UNBLUR_EXECUTABLE");
			if (penv!=NULL)
				fn_unblur_exe = (std::string)penv;
		}
		// Get the SUMMOVIE executable
		if (fn_summovie_exe == "")
		{
			char * penv;
			penv = getenv ("RELION_SUMMOVIE_EXECUTABLE");
			if (penv!=NULL)
				fn_summovie_exe = (std::string)penv;
		}

		if (angpix < 0)
			REPORT_ERROR("ERROR: For Unblur it is mandatory to provide the pixel size in Angstroms through --angpix.");
	}
	else if (do_motioncor2)
	{
		// Get the MOTIONCOR2 executable
		if (fn_motioncor2_exe == "")
		{
			char * penv;
			penv = getenv ("RELION_MOTIONCOR2_EXECUTABLE");
			if (penv!=NULL)
				fn_motioncor2_exe = (std::string)penv;
		}
	}
	else
		REPORT_ERROR(" ERROR: You have to specify which programme to use through either --use_motioncor2 or --use_unblur");

	if (do_dose_weighting)
	{
		if (!(do_unblur || do_motioncor2))
			REPORT_ERROR("ERROR: Dose-weighting can only be done by UNBLUR or MOTIONCOR2.");
		if (angpix < 0)
			REPORT_ERROR("ERROR: For dose-weighting it is mandatory to provide the pixel size in Angstroms through --angpix.");

	}

#ifdef CUDA
	if (!do_unblur)
	{
		if (gpu_ids.length() > 0)
			untangleDeviceIDs(gpu_ids, allThreadIDs);
		else if (verb>0)
			std::cout << "gpu-ids not specified, threads will automatically be mapped to devices (incrementally)."<< std::endl;
		HANDLE_ERROR(cudaGetDeviceCount(&devCount));
	}
#endif

	// Set up which micrograph movies to run MOTIONCOR2 on
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
		}
	}
	else
	{
		fn_in.globFiles(fn_micrographs);
	}

	// If we're continuing an old run, see which micrographs have not been finished yet...
	fn_ori_micrographs = fn_micrographs;
	if (continue_old)
	{
		fn_micrographs.clear();
		for (long int imic = 0; imic < fn_ori_micrographs.size(); imic++)
		{
			FileName fn_avg, fn_mov;
			getOutputFileNames(fn_ori_micrographs[imic], fn_avg, fn_mov);
			if (!exists(fn_avg) || (do_save_movies && !exists(fn_mov)) )
				fn_micrographs.push_back(fn_ori_micrographs[imic]);
		}
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

	if (verb > 0)
	{
		if (do_unblur)
			std::cout << " Using UNBLUR executable in: " << fn_unblur_exe << std::endl;
		else
			std::cout << " Using MOTIONCOR2 executable in: " << fn_motioncor2_exe << std::endl;
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
			std::cout << " Correcting beam-induced motions using Tim Grant's UNBLUR ..." << std::endl;
		else if (do_motioncor2)
			std::cout << " Correcting beam-induced motions using Shawn Zheng's MOTIONCOR2 ..." << std::endl;
		else
			REPORT_ERROR("Bug: by now it should be clear whether to use MotionCor2 or Unblur...");

		init_progress_bar(fn_micrographs.size());
		barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
	}

	for (long int imic = 0; imic < fn_micrographs.size(); imic++)
	{
		std::vector<float> xshifts, yshifts;

		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);


		bool result = false;
		if (do_unblur)
			result = executeUnblur(fn_micrographs[imic], xshifts, yshifts);
		else if (do_motioncor2)
			result = executeMotioncor2(fn_micrographs[imic], xshifts, yshifts);
		else
			REPORT_ERROR("Bug: by now it should be clear whether to use MotionCor2 or Unblur...");

		if (result)
			plotShifts(fn_micrographs[imic], xshifts, yshifts);

	}

	if (verb > 0)
		progress_bar(fn_micrographs.size());

	// Make a logfile with the shifts in pdf format and write output STAR files
	generateLogFilePDFAndWriteStarFiles();


}

bool MotioncorrRunner::executeMotioncor2(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts, int rank)
{


	FileName fn_avg, fn_mov;
	getOutputFileNames(fn_mic, fn_avg, fn_mov);

	FileName fn_out = fn_avg.withoutExtension() + ".out";
	FileName fn_log = fn_avg.withoutExtension() + ".log";
	FileName fn_err = fn_avg.withoutExtension() + ".err";
	FileName fn_cmd = fn_avg.withoutExtension() + ".com";

	std::string command = fn_motioncor2_exe + " ";

	if (fn_mic.getExtension() == "tif" || fn_mic.getExtension() == "tiff")
		command += " -InTiff " + fn_mic;
	else
		command += " -InMrc " + fn_mic;

	command += " -OutMrc " + fn_avg;
        command += " -LogFile " + fn_avg.withoutExtension();
	command += " -Bft " + floatToString(bfactor);
	command += " -PixSize " + floatToString(angpix);

	if (do_save_movies)
		command += " -OutStack 1";

	command += " -Patch " + integerToString(patch_x) + " " + integerToString(patch_y);

	if (group > 1)
		command += " -Group " + integerToString(group);

	if (fn_gain_reference != "")
		command += " -Gain " + fn_gain_reference;

	// Throw away first few frames?
	if (first_frame_sum > 1)
		command += " -Throw " + integerToString(first_frame_sum - 1);

	// Throw away last few frames?
	if (last_frame_sum > 0)
	{
		// Cannot read TIFFs here...
		if (fn_mic.getExtension() == "tif" || fn_mic.getExtension() == "tiff")
			REPORT_ERROR("ERROR: you cannot use --last_frame_sum > 0 when reading TIFFs. Try using -Trunk in the --other_motioncor2_args instead.");

		// Read in header of the movie, to see how many frames it has
		Image<RFLOAT> Ihead;
		Ihead.read(fn_mic, false);
		int n_frames = NSIZE(Ihead());

		int trunc = n_frames - last_frame_sum;
		if (trunc < 0)
			REPORT_ERROR("ERROR: movie " + fn_mic + " does not have enough frames.");
		command += " -Trunc " + integerToString(trunc);
	}

	if (bin_factor > 1)
		command += " -FtBin " + floatToString(bin_factor);

	if (do_dose_weighting)
	{
		command += " -Kv " + floatToString(voltage);
		command += " -FmDose " + floatToString(dose_per_frame);
		command += " -InitDose " + floatToString(pre_exposure);
	}

	if (fn_defect != "")
		command += " -DefectFile " + fn_defect;

	if (fn_archive != "")
		command += " -ArcDir " + fn_archive;

	if (fn_other_motioncor2_args.length() > 0)
		command += " " + fn_other_motioncor2_args;

	if ( allThreadIDs.size() == 0)
	{
		// Automated mapping
		command += " -Gpu " + integerToString(rank % devCount);
	}
	else
	{
		if (rank >= allThreadIDs.size())
			REPORT_ERROR("ERROR: not enough MPI nodes specified for the GPU IDs.");
		command += " -Gpu " + allThreadIDs[rank][0];
	}

	command += " >> " + fn_out + " 2>> " + fn_err;

	// Save the command that was executed
	std::ofstream fh;
	fh.open(fn_cmd.c_str(), std::ios::out);
	fh << command << std::endl;
	fh.close();

	if (system(command.c_str()))
	{
		std::cerr << " WARNING: there was an error executing: " << command << std::endl;
		return false;
	}
	else
	{
		if (do_dose_weighting)
		{
			// Move .mrc to _noDW.mrc filename
			FileName fn_tmp = fn_avg.withoutExtension() + "_noDW.mrc";
			if (std::rename(fn_avg.c_str(), fn_tmp.c_str()))
			{
				std::cerr << "ERROR in renaming: " << fn_avg << " to " << fn_tmp <<std::endl;
				return false;
			}
			// Move _DW.mrc to .mrc filename
			fn_tmp = fn_avg.withoutExtension() + "_DW.mrc";
			if (std::rename(fn_tmp.c_str(), fn_avg.c_str()))
			{
				std::cerr << "ERROR in renaming: " << fn_tmp << " to " << fn_avg <<std::endl;
				return false;
			}
		}

		if (do_save_movies)
		{
			// Move movie .mrc to new .mrcs filename
			FileName fn_tmp = fn_avg.withoutExtension() + "_Stk.mrc";
			if (std::rename(fn_tmp.c_str(), fn_mov.c_str()))
			{
				std::cerr << "ERROR in renaming: " << fn_tmp << " to " << fn_mov <<std::endl;
				return false;
			}
		}

		// Also analyse the shifts
		getShiftsMotioncor2(fn_out, xshifts, yshifts);
	}

	// Success!
	return true;
}

void MotioncorrRunner::getShiftsMotioncor2(FileName fn_log, std::vector<float> &xshifts, std::vector<float> &yshifts)
{

	std::ifstream in(fn_log.data(), std::ios_base::in);
	if (in.fail())
		return;

	xshifts.clear();
	yshifts.clear();

    std::string line;

    // Start reading the ifstream at the top
    in.seekg(0);

    // Read through the shifts file
    int i = 0;
    bool have_found_final = false;
    while (getline(in, line, '\n'))
    {
    	// ignore all commented lines, just read first two lines with data
    	if (line.find("Full-frame alignment shift") != std::string::npos)
    	{
    		have_found_final = true;
    	}
    	else if (have_found_final)
    	{
    		size_t shiftpos = line.find("shift:");
    		if (shiftpos != std::string::npos)
    		{
    			std::vector<std::string> words;
    			tokenize(line.substr(shiftpos+7), words);;
        		if (words.size() < 2)
        		{
        			std::cerr << " fn_log= " << fn_log << std::endl;
        			REPORT_ERROR("ERROR: unexpected number of words on line from MOTIONCORR logfile: " + line);
        		}
        		xshifts.push_back(textToFloat(words[0]));
    			yshifts.push_back(textToFloat(words[1]));
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
bool MotioncorrRunner::executeUnblur(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts)
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
		int res = symlink(fn_link.c_str(), fn_tmp_mic.c_str());
	}
	else
	{
		fn_tmp_mic = fn_mic;
	}
	FileName fn_tmp_mov = fn_mov.withoutExtension() + ".mrc";

	Image<RFLOAT> Itest;
	Itest.read(fn_mic, false);
	int Nframes = NSIZE(Itest());

	std::ofstream  fh;
	fh.open((fn_com).c_str(), std::ios::out);
	if (!fh)
	 REPORT_ERROR( (std::string)"executeUnblur cannot create file: " + fn_com);

	// Write script to run ctffind
	fh << "#!/usr/bin/env csh"<<std::endl;
	fh << "setenv  OMP_NUM_THREADS " << integerToString(nr_threads)<<std::endl;
	fh << fn_unblur_exe << " > " << fn_log << "  << EOF"<<std::endl;
	fh << fn_tmp_mic << std::endl;
	fh << Nframes << std::endl;
	fh << fn_avg << std::endl;
	fh << fn_shifts << std::endl;
	fh << angpix << std::endl; // pixel size
	if (do_dose_weighting)
	{
		fh << "YES" << std::endl; // apply dose weighting
		fh << dose_per_frame << std::endl;
		fh << voltage << std::endl;
		fh << pre_exposure << std::endl;
	}
	else
	{
		fh << "NO" << std::endl; // no dose filtering
	}
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
	{
		std::cerr << "ERROR in executing: " << command << std::endl;
		return false;
    }

	// Also analyse the shifts
	getShiftsUnblur(fn_shifts, xshifts, yshifts);

	// If the requested sum is only a subset, then use summovie to make the average
	int mylastsum = (last_frame_sum == 0) ? Nframes : last_frame_sum;
	if (first_frame_sum != 1 || mylastsum != Nframes)
	{
		FileName fn_com2 = fn_root + "_summovie.com";
		FileName fn_log2 = fn_root + "_summovie.log";
		FileName fn_frc = fn_root + "_frc.txt";

		std::ofstream  fh2;
		fh2.open((fn_com2).c_str(), std::ios::out);
		if (!fh2)
		 REPORT_ERROR( (std::string)"executeUnblur cannot create file: " + fn_com2);

		// Write script to run ctffind
		fh2 << "#!/usr/bin/env csh"<<std::endl;
		fh2 << "setenv  OMP_NUM_THREADS " << integerToString(nr_threads)<<std::endl;
		fh2 << fn_summovie_exe << " > " << fn_log2 << "  << EOF"<<std::endl;
		fh2 << fn_tmp_mic << std::endl;
		fh2 << Nframes << std::endl;
		fh2 << fn_avg << std::endl;
		fh2 << fn_shifts << std::endl;
		fh2 << fn_frc << std::endl;
		fh2 << first_frame_sum << std::endl;
		fh2 << mylastsum << std::endl;
		fh2 << angpix << std::endl; // pixel size
		fh2 << "NO" << std::endl; // dont set expert options
		fh2 <<"EOF"<<std::endl;
		fh2.close();

		// Execute summovie
		std::string command2 = "csh "+ fn_com2;
		if (system(command2.c_str()))
		{
			std::cerr << "ERROR in executing: " << command2 <<std::endl;
			return false;
		}

		// Plot the FRC
		plotFRC(fn_frc);
	}

	// Move movie .mrc to new .mrcs filename
	if (do_save_movies)
	{
		if (std::rename(fn_tmp_mov.c_str(), fn_mov.c_str()))
		{
			std::cerr << "ERROR in renaming: " << fn_tmp_mov << " to " << fn_mov <<std::endl;
			return false;
		}
	}

	// remove symbolic link
	std::remove(fn_tmp_mic.c_str());

	// Success!
	return true;
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

// Plot the UNBLUR FRC curve
void MotioncorrRunner::plotFRC(FileName fn_frc)
{

	std::ifstream in(fn_frc.data(), std::ios_base::in);
	if (in.fail())
		return;

	FileName fn_eps = fn_frc.withoutExtension() + ".eps";
	CPlot2D *plot2D=new CPlot2D(fn_eps);
 	plot2D->SetXAxisSize(600);
 	plot2D->SetYAxisSize(600);
 	CDataSet dataSet;
	dataSet.SetDrawMarker(false);
	dataSet.SetDatasetColor(0.0,0.0,0.0);

 	std::string line, token;
	// Read through the frc file
	while (getline(in, line, '\n'))
	{
		// ignore all commented lines, just read first two lines with data
		if (line[0] != '#')
		{
			std::vector<std::string> words;
			tokenize(line, words);
			if (words.size() < 2)
			{
				break;
			}
			else
			{
				CDataPoint point(textToFloat(words[0]), textToFloat(words[1]));
			}
		}
	}
	in.close();

	plot2D->SetXAxisTitle("Resolution (1/Angstrom)");
	plot2D->SetYAxisTitle("FRC");
	plot2D->OutputPostScriptPlot(fn_eps);

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

	if (do_unblur)
	{
		plot2D->SetXAxisTitle("X-shift (in Angstroms)");
		plot2D->SetYAxisTitle("Y-shift (in Angstroms)");
	}
	else
	{
		plot2D->SetXAxisTitle("X-shift (in pixels)");
		plot2D->SetYAxisTitle("Y-shift (in pixels)");
	}
	plot2D->OutputPostScriptPlot(fn_eps);


}


void MotioncorrRunner::generateLogFilePDFAndWriteStarFiles()
{

	if (fn_ori_micrographs.size() > 0)
	{

		std::vector<FileName> fn_eps;

		FileName fn_prev="";
		for (long int i = 0; i < fn_ori_micrographs.size(); i++)
		{
			if (fn_prev != fn_ori_micrographs[i].beforeLastOf("/"))
			{
				fn_prev = fn_ori_micrographs[i].beforeLastOf("/");
				fn_eps.push_back(fn_out + fn_prev+"/*.eps");
			}
		}

		joinMultipleEPSIntoSinglePDF(fn_out + "logfile.pdf ", fn_eps);

	}

	// Also write out the output star files
	MDavg.clear();
	MDmov.clear();

	for (long int imic = 0; imic < fn_ori_micrographs.size(); imic++)
	{
		// For output STAR file
		FileName fn_avg, fn_mov;
		getOutputFileNames(fn_ori_micrographs[imic], fn_avg, fn_mov);
		if (exists(fn_avg))
		{
			MDavg.addObject();
			if (do_dose_weighting && do_motioncor2)
			{
				FileName fn_avg_wodose = fn_avg.withoutExtension() + "_noDW.mrc";
				MDavg.setValue(EMDL_MICROGRAPH_NAME_WODOSE, fn_avg_wodose);
			}
			MDavg.setValue(EMDL_MICROGRAPH_NAME, fn_avg);
			if (do_save_movies && exists(fn_mov))
			{
				MDmov.addObject();
				MDmov.setValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_mov);
				MDmov.setValue(EMDL_MICROGRAPH_NAME, fn_avg);
			}
		}
	}

	// Write out STAR files at the end
	MDavg.write(fn_out + "corrected_micrographs.star");
	if (do_save_movies)
		MDmov.write(fn_out + "corrected_micrograph_movies.star");

}


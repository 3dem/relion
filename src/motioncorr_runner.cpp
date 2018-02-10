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
#include "src/micrograph_model.h"
#include "src/fftw.h"
#include "src/matrix2d.h"
#include "src/matrix1d.h"
#include <omp.h>

#define TIMING
#ifdef TIMING
	#define RCTIC(label) (timer.tic(label))
	#define RCTOC(label) (timer.toc(label))

	Timer timer;
	int TIMING_READ_GAIN = timer.setNew("read gain");
	int TIMING_READ_MOVIE = timer.setNew("read movie");
	int TIMING_APPLY_GAIN = timer.setNew("apply gain");
	int TIMING_INITIAL_SUM = timer.setNew("initial sum");
	int TIMING_DETECT_HOT = timer.setNew("detect hot pixels");
	int TIMING_GLOBAL_FFT = timer.setNew("global FFT");
	int TIMING_GLOBAL_ALIGNMENT = timer.setNew("global alignment");
	int TIMING_GLOBAL_IFFT = timer.setNew("global iFFT");
	int TIMING_PREP_PATCH = timer.setNew("prepare patch");
	int TIMING_CLIP_PATCH = timer.setNew("prepare patch - real space clip");
	int TIMING_PATCH_FFT = timer.setNew("prepare patch - FFT");
	int TIMING_PATCH_ALIGN = timer.setNew("patch alignment");
	int TIMING_PREP_WEIGHT = timer.setNew("align - prep weight");
	int TIMING_MAKE_REF = timer.setNew("align - make reference");
	int TIMING_CCF_CALC = timer.setNew("align - calc CCF");
	int TIMING_CCF_IFFT = timer.setNew("align - iFFT CCF");
	int TIMING_CCF_RECENTRE = timer.setNew("align - recentre CCF");
	int TIMING_CCF_FIND_MAX = timer.setNew("align - argmax CCF");
	int TIMING_FOURIER_SHIFT = timer.setNew("align - shift in Fourier space");
	int TIMING_FIT_POLYNOMIAL = timer.setNew("fit polynomial");
	int TIMING_DOSE_WEIGHTING = timer.setNew("dose weighting");
	int TIMING_DW_WEIGHT = timer.setNew("dose weighting - calc weight");
	int TIMING_DW_IFFT = timer.setNew("dose weighting - iFFT");
	int TIMING_REAL_SPACE_INTERPOLATION = timer.setNew("real space interpolation");
	int TIMING_BINNING = timer.setNew("binning");
//	int TIMING_ = timer.setNew("");

#else
	#define RCTIC(timer,label)
	#define RCTOC(timer,label)
#endif

void MotioncorrRunner::read(int argc, char **argv, int rank)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "STAR file with all input micrographs, or a Linux wildcard with all micrographs to operate on");
	fn_out = parser.getOption("--o", "Name for the output directory", "MotionCorr");
	n_threads = textToInteger(parser.getOption("--j", "Number of threads per movie (= process)", "1"));
	fn_movie = parser.getOption("--movie", "Rootname to identify movies", "movie");
	continue_old = parser.checkOption("--only_do_unfinished", "Only run motion correction for those micrographs for which there is not yet an output micrograph.");
	do_save_movies  = parser.checkOption("--save_movies", "Also save the motion-corrected movies.");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms","-1"));
	first_frame_sum =  textToInteger(parser.getOption("--first_frame_sum", "First movie frame used in output sum (start at 1)", "1"));
	if (first_frame_sum < 1) first_frame_sum = 1;
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

	do_own = parser.checkOption("--use_own","Use our own implementation of motion correction");
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
	else if (do_own) {
		std::cout << "     !!! WARNING !!!" << std::endl << " Our own implementation of motion correction is under development!" << std::endl;
	} else {
		REPORT_ERROR(" ERROR: You have to specify which programme to use through either --use_motioncor2 or --use_unblur");
	}

	if (do_dose_weighting)
	{
		if (!(do_unblur || do_motioncor2 || do_own))
			REPORT_ERROR("ERROR: Dose-weighting can only be done by UNBLUR or MOTIONCOR2.");
		if (angpix < 0)
			REPORT_ERROR("ERROR: For dose-weighting it is mandatory to provide the pixel size in Angstroms through --angpix.");

	}

#ifdef CUDA
	if (do_motioncor2)
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
		if (do_own)
			std::cout << " Correcting beam-induced motions using our own implementation ..." << std::endl;
		else if (do_unblur)
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
		if (do_own)
			result = executeOwnMotionCorrection(fn_micrographs[imic], xshifts, yshifts);
		else if (do_unblur)
			result = executeUnblur(fn_micrographs[imic], xshifts, yshifts);
		else if (do_motioncor2)
			result = executeMotioncor2(fn_micrographs[imic], xshifts, yshifts);
		else
			REPORT_ERROR("Bug: by now it should be clear whether to use MotionCor2 or Unblur...");

		if (result) {
			saveModel(fn_micrographs[imic], xshifts, yshifts);
			plotShifts(fn_micrographs[imic], xshifts, yshifts);
		}
	}

	if (verb > 0)
		progress_bar(fn_micrographs.size());

	// Make a logfile with the shifts in pdf format and write output STAR files
	generateLogFilePDFAndWriteStarFiles();

#ifdef TIMING
        timer.printTimes(false);
#endif
#ifdef TIMING_FFTW
	timer_fftw.printTimes(false);
#endif
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

	// Always take the first frame to be aligned as the origin of the motion
	command += " -FmRef 0";

	// Throw away last few frames?
	if (last_frame_sum > 0)
	{
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

	// Write script to run Unblur
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

void MotioncorrRunner::saveModel(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts) {
	Micrograph m(fn_mic, fn_gain_reference, bin_factor);

	FileName fn_avg, fn_mov;
        getOutputFileNames(fn_mic, fn_avg, fn_mov);

	for (int i = 0, ilim = xshifts.size(); i < ilim; i++) {
		int frame = i + 1;

		// UNBLUR processes all frames, but MotionCor2 not.
		// So we have to adjust...
		if (do_motioncor2) frame += (first_frame_sum - 1);

		m.setGlobalShift(frame, xshifts[i], yshifts[i]);
	}

	m.write(fn_avg.withoutExtension() + ".star");
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

// TODO:
// - defect
// - grouping

bool MotioncorrRunner::executeOwnMotionCorrection(FileName fn_mic, std::vector<float> &xshifts, std::vector<float> &yshifts) {
	std::cout << std::endl;
	std::cout << "Now working on " << fn_mic << " nthreads = " << n_threads << std::endl;
	omp_set_num_threads(n_threads);

	Image<RFLOAT> Ihead, Igain, Iref;
	std::vector<MultidimArray<Complex> > Fframes;
	std::vector<Image<RFLOAT> > Iframes;
	std::vector<int> frames;
	std::vector<FourierTransformer> transformers(n_threads);
	FourierTransformer transformer;

	const int hotpixel_sigma = 6;

	// Check image size
	Ihead.read(fn_mic, false);
	const int nx = XSIZE(Ihead()), ny = YSIZE(Ihead()), nn = NSIZE(Ihead());

	// Which frame to use?
	std::cout << "Movie X = " << nx << " Y = " << ny << " N = " << nn << std::endl;
	std::cout << "Frames to be used:";
	for (int i = 0; i < nn; i++) {
		// For users, all numbers are 1-indexed. Internally they are 0-indexed.
		int frame = i + 1;
		if (frame < first_frame_sum) continue;
		if (last_frame_sum > 0 && frame > last_frame_sum) continue;
		frames.push_back(i);
		std::cout << " " << frame;
	}
	std::cout << std::endl;
	
	const int n_frames = frames.size();
	Iframes.resize(n_frames);
	Fframes.resize(n_frames);
	xshifts.resize(n_frames);
	yshifts.resize(n_frames);

	RCTIC(TIMING_READ_GAIN);
	if (fn_gain_reference != "") {
		Igain.read(fn_gain_reference);
		if (XSIZE(Igain()) != nx || YSIZE(Igain()) != ny) {
			REPORT_ERROR("The size of the image and the size of the gain reference do not match.");
		}
	}
	RCTOC(TIMING_READ_GAIN);

	// Read images
	RCTIC(TIMING_READ_MOVIE);
	#pragma omp parallel for
	for (int iframe = 0; iframe < n_frames; iframe++) {
		Iframes[iframe].read(fn_mic, true, frames[iframe]);
	}
	RCTOC(TIMING_READ_MOVIE);

	// Apply gain
	RCTIC(TIMING_APPLY_GAIN);
	if (fn_gain_reference != "") {
		for (int iframe = 0; iframe < n_frames; iframe++) {
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Igain()) {
				DIRECT_MULTIDIM_ELEM(Iframes[iframe](), n) *= DIRECT_MULTIDIM_ELEM(Igain(), n);
			}
		}
	}
	RCTOC(TIMING_APPLY_GAIN);

	MultidimArray<RFLOAT> Iframe(ny, nx);
	Iframe.initZeros();

	// Simple unaligned sum
	RCTIC(TIMING_INITIAL_SUM);
	for (int iframe = 0; iframe < n_frames; iframe++) {
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe) {
			DIRECT_MULTIDIM_ELEM(Iframe, n) += DIRECT_MULTIDIM_ELEM(Iframes[iframe](), n);
		}
	}
	RCTOC(TIMING_INITIAL_SUM);

	// Hot pixel
	RCTIC(TIMING_DETECT_HOT);
	RFLOAT mean = 0, std = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe) {
		mean += DIRECT_MULTIDIM_ELEM(Iframe, n);
	}
	mean /=  YXSIZE(Iframe);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe) {
		RFLOAT d = (DIRECT_MULTIDIM_ELEM(Iframe, n) - mean);
		std += d * d;
	}
	std = std::sqrt(std / YXSIZE(Iframe));
	const RFLOAT threshold = mean + hotpixel_sigma * std;
	std::cout << "Mean = " << mean << " Std = " << std << ", Hotpixel threshold = " << threshold << std::endl;

	MultidimArray<bool> bBad(ny, nx);
	bBad.initZeros();
	int n_bad = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe) {
		if (DIRECT_MULTIDIM_ELEM(Iframe, n) > threshold) {
			DIRECT_MULTIDIM_ELEM(bBad, n) = true;
			n_bad++;
		}
	}
	std::cout << "Detected " << n_bad << " hot pixels to be corrected. (BUT correction not implemented yet)" << std::endl;
	RCTOC(TIMING_DETECT_HOT);
	// TODO: fix defects

	// FFT
	RCTIC(TIMING_GLOBAL_FFT);
	#pragma omp parallel for
	for (int iframe = 0; iframe < n_frames; iframe++) {
//		std::cout << iframe << " thread " << omp_get_thread_num() << std::endl;
		transformers[omp_get_thread_num()].FourierTransform(Iframes[iframe](), Fframes[iframe]);
	}
	RCTOC(TIMING_GLOBAL_FFT);

	// Global alignment
	RCTIC(TIMING_GLOBAL_ALIGNMENT);
	alignPatch(Fframes, nx, ny, xshifts, yshifts);
	RCTOC(TIMING_GLOBAL_ALIGNMENT);

	std::cout << "Global alignment done." << std::endl;	
	Iref().reshape(Iframes[0]());
	Iref().initZeros();
	RCTIC(TIMING_GLOBAL_IFFT);
	#pragma omp parallel for
	for (int iframe = 0; iframe < n_frames; iframe++) {
		transformers[omp_get_thread_num()].inverseFourierTransform(Fframes[iframe], Iframes[iframe]());
	}
	RCTOC(TIMING_GLOBAL_IFFT);

	// Patch based alignment
	std::cout << "Full Size: X = " << nx << " Y = " << ny << std::endl;
	std::cout << "Patches: X = " << patch_x << " Y = " << patch_y << std::endl;
	const int patch_nx = nx / patch_x, patch_ny = ny / patch_y, n_patches = patch_x * patch_y;
	std::vector<RFLOAT> patch_xshifts, patch_yshifts, patch_frames, patch_xs, patch_ys;
	std::vector<MultidimArray<Complex> > Fpatches(n_frames);

	int ipatch = 1;
	for (int iy = 0; iy < patch_y; iy++) {
		for (int ix = 0; ix < patch_x; ix++) {
			int x_start = ix * patch_nx, y_start = iy * patch_ny;
			int x_end = x_start + patch_nx, y_end = y_start + patch_ny;
			if (x_end > nx) x_end = nx;
			if (y_end > ny) y_end = ny;

			int x_center = (x_start + x_end - 1) / 2, y_center = (y_start + y_end - 1) / 2;
			std::cout << "Patch (" << iy + 1 << ", " << ix + 1 << ") " << ipatch << " / " << patch_x * patch_y;
			std::cout << ", X range = [" << x_start << ", " << x_end << "), Y range = [" << y_start << ", " << y_end << ")";
			std::cout << ", Center = (" << x_center << ", " << y_center << ")" << std::endl;
			ipatch++;

			std::vector<float> local_xshifts(n_frames), local_yshifts(n_frames);
			RCTIC(TIMING_PREP_PATCH);
			MultidimArray<RFLOAT> Iframe(y_end - y_start + 1, x_end - x_start + 1);
			for (int iframe = 0; iframe < n_frames; iframe++) {
				RCTIC(TIMING_CLIP_PATCH);
				for (int ipy = y_start; ipy < y_end; ipy++) {
					for (int ipx = x_start; ipx < x_end; ipx++) {
						DIRECT_A2D_ELEM(Iframe, ipy - y_start, ipx - x_start) = DIRECT_A2D_ELEM(Iframes[iframe](), ipy, ipx);
					}
				}
				RCTOC(TIMING_CLIP_PATCH);

				RCTIC(TIMING_PATCH_FFT);
				transformer.FourierTransform(Iframe, Fpatches[iframe]);
				RCTOC(TIMING_PATCH_FFT);
			}
			RCTOC(TIMING_PREP_PATCH);
			
			RCTIC(TIMING_PATCH_ALIGN);
			bool converged = alignPatch(Fpatches, x_end - x_start + 1, y_end - y_start + 1, local_xshifts, local_yshifts);
			RCTOC(TIMING_PATCH_ALIGN);
			if (!converged) continue;

			for (int iframe = 0; iframe < n_frames; iframe++) {
				patch_xshifts.push_back(local_xshifts[iframe]);
				patch_yshifts.push_back(local_yshifts[iframe]);
				patch_frames.push_back(iframe);
				patch_xs.push_back(x_center);
				patch_ys.push_back(y_center);
			}
		}
	}
	Fpatches.clear();

	// Fit polynomial model

	// TODO: outlier rejection
	RCTIC(TIMING_FIT_POLYNOMIAL);
	const int n_obs = patch_frames.size();
	const int n_params = 18;
	Matrix2D <RFLOAT> matA(n_obs, n_params);
	Matrix1D <RFLOAT> vecX(n_obs), vecY(n_obs), coeffX(n_params), coeffY(n_params);
	for (int i = 0; i < n_obs; i++) {
		VEC_ELEM(vecX, i) = patch_xshifts[i]; VEC_ELEM(vecY, i) = patch_yshifts[i];

		const RFLOAT x = patch_xs[i] / nx - 0.5;
		const RFLOAT y = patch_ys[i] / ny - 0.5;
		const RFLOAT z = patch_frames[i];
		const RFLOAT x2 = x * x, y2 = y * y, xy = x * y, z2 = z * z;
		const RFLOAT z3 = z2 * z;		

		MAT_ELEM(matA, i, 0)  =      z;
		MAT_ELEM(matA, i, 1)  =      z2; 
		MAT_ELEM(matA, i, 2)  =      z3;
 
		MAT_ELEM(matA, i, 3)  = x  * z;
		MAT_ELEM(matA, i, 4)  = x  * z2;
		MAT_ELEM(matA, i, 5)  = x  * z3;

		MAT_ELEM(matA, i, 6)  = x2 * z;
		MAT_ELEM(matA, i, 7)  = x2 * z2;
		MAT_ELEM(matA, i, 8)  = x2 * z3;

		MAT_ELEM(matA, i, 9)  = y  * z;
		MAT_ELEM(matA, i, 10) = y  * z2;
		MAT_ELEM(matA, i, 11) = y  * z3;

		MAT_ELEM(matA, i, 12) = y2 * z;
		MAT_ELEM(matA, i, 13) = y2 * z2;
		MAT_ELEM(matA, i, 14) = y2 * z3;

		MAT_ELEM(matA, i, 15) = xy * z;
		MAT_ELEM(matA, i, 16) = xy * z2;
		MAT_ELEM(matA, i, 17) = xy * z3;
	}

	const RFLOAT EPS = 1e-10;
	solve(matA, vecX, coeffX, EPS);
	solve(matA, vecY, coeffY, EPS);

#ifdef DEBUG
	std::cout << "Polynomial fitting coefficients for X and Y:" << std::endl;
	for (int i = 0; i < n_params; i++) {
		std::cout << i << " " << coeffX(i) << " " << coeffY(i) << std::endl;
	}
#endif

	RFLOAT rms_x = 0, rms_y = 0;
	RFLOAT x_fitted, y_fitted;
        for (int i = 0; i < n_obs; i++) {
		const RFLOAT x = patch_xs[i] / nx - 0.5;
		const RFLOAT y = patch_ys[i] / ny - 0.5;
		const RFLOAT z = patch_frames[i];

		getFittedXY(x, y, z, coeffX, coeffY, x_fitted, y_fitted);
		rms_x += (patch_xshifts[i] - x_fitted) * (patch_xshifts[i] - x_fitted);
		rms_y += (patch_yshifts[i] - y_fitted) * (patch_yshifts[i] - y_fitted);
#ifdef DEBUG
		std::cout << "x = " << x << " y = " << y << " z = " << z;
		std::cout << ", Xobs = " << patch_xshifts[i] << " Xfit = " << x_fitted;
		std::cout << ", Yobs = " << patch_yshifts[i] << " Yfit = " << y_fitted << std::endl;
#endif
	}
	rms_x = std::sqrt(rms_x / n_obs); rms_y = std::sqrt(rms_y / n_obs);
	std::cout << "Polynomial fit RMSD X = " << rms_x << " px, Y = " << rms_y << " px" << std::endl;
	RCTOC(TIMING_FIT_POLYNOMIAL);

	// Dose weighting
	RCTIC(TIMING_DOSE_WEIGHTING);
	if (do_dose_weighting) {
		if (std::abs(voltage - 300) > 2 && std::abs(voltage - 200) > 2) {
			REPORT_ERROR("Sorry, dose weighting is supported only for 300 kV or 200 kV");
		}

		std::vector <RFLOAT> doses(n_frames);
		for (int iframe = 0; iframe < n_frames; iframe++) {
			// dose AFTER each frame.
			doses[iframe] = pre_exposure + dose_per_frame * (frames[iframe] + 1);
			if (std::abs(voltage = 200) <= 5) {
				doses[iframe] /= 0.8; // 200 kV electron is more damaging.
			}
		}

		RCTIC(TIMING_DW_WEIGHT);
		doseWeighting(Fframes, doses);
		RCTOC(TIMING_DW_WEIGHT);

		// Update real space images
		RCTIC(TIMING_DW_IFFT);
		#pragma omp parallel for
		for (int iframe = 0; iframe < n_frames; iframe++) {
			transformers[omp_get_thread_num()].inverseFourierTransform(Fframes[iframe], Iframes[iframe]());
		}
		RCTOC(TIMING_DW_IFFT);
	}
	RCTOC(TIMING_DOSE_WEIGHTING);

	RCTIC(TIMING_REAL_SPACE_INTERPOLATION);
	std::cout << "Real space interpolation: ";
	Iref().initZeros();
	for (int iframe = 0; iframe < n_frames; iframe++) {
		std::cout << "." << std::flush;
		#pragma omp parallel for schedule(static)
		for (int ix = 0; ix < nx; ix++) {
			for (int iy = 0; iy < ny; iy++) {
				bool valid = true;
				const RFLOAT x = (float)ix / nx - 0.5;
				const RFLOAT y = (float)iy / ny - 0.5;
				getFittedXY(x, y, iframe, coeffX, coeffY, x_fitted, y_fitted);
				x_fitted = ix - x_fitted; y_fitted = iy - y_fitted;

				int x0 = FLOOR(x_fitted);
				int y0 = FLOOR(y_fitted);
				const int x1 = x0 + 1;
				const int y1 = y0 + 1;

				if (x0 < 0) {x0 = 0; valid = false;}
				if (y0 < 0) {y0 = 0; valid = false;}
				if (x1 >= nx) {x0 = nx - 1; valid = false;}
				if (y1 >= ny) {y0 = ny - 1; valid = false;}
				if (!valid) {
					DIRECT_A2D_ELEM(Iref(), iy, ix) += DIRECT_A2D_ELEM(Iframes[iframe](), y0, x0);
					if (std::isnan(DIRECT_A2D_ELEM(Iref(), iy, ix))) {
						std::cout << "ix = " << ix << " xfit = " << x_fitted << " iy = " << iy << " ifit = " << y_fitted << std::endl;
					}
					continue;
				}

				const RFLOAT fx = x_fitted - x0;
				const RFLOAT fy = y_fitted - y0;

//				std::cout << "ix = " << ix << " xfit = " << x_fitted << " iy = " << iy << " ifit = " << y_fitted << std::endl;
				const RFLOAT d00 = DIRECT_A2D_ELEM(Iframes[iframe](), y0, x0);
				const RFLOAT d01 = DIRECT_A2D_ELEM(Iframes[iframe](), y0, x1);
				const RFLOAT d10 = DIRECT_A2D_ELEM(Iframes[iframe](), y1, x0);
				const RFLOAT d11 = DIRECT_A2D_ELEM(Iframes[iframe](), y1, x1);

				const RFLOAT dx0 = LIN_INTERP(fx, d00, d01);
				const RFLOAT dx1 = LIN_INTERP(fx, d10, d11);
				const RFLOAT val = LIN_INTERP(fy, dx0, dx1);
#ifdef DEBUG
				if (std::isnan(val)) {
					std::cout << "ix = " << ix << " xfit = " << x_fitted << " iy = " << iy << " ifit = " << y_fitted << " d00 " << d00 << " d01 " << d01 << " d10 " << d10 << " d11 " << d11 << " dx0 " << dx0 << " dx1 " << dx1 << std::endl;
				}
#endif
				DIRECT_A2D_ELEM(Iref(), iy, ix) += val;
			}
		}
	}
	std::cout << " done" << std::endl;
	RCTOC(TIMING_REAL_SPACE_INTERPOLATION);
	
	// Apply binning
	RCTIC(TIMING_BINNING);
	if (bin_factor != 1) {
		binNonSquareImage(Iref, bin_factor);
	}
	RCTOC(TIMING_BINNING);
	
	// Final output
	FileName fn_avg, fn_mov;
	getOutputFileNames(fn_mic, fn_avg, fn_mov);
	Iref.write(fn_avg);
	std::cout << "Written aligned sum to " << fn_avg << std::endl;
	return true;
}

bool MotioncorrRunner::alignPatch(std::vector<MultidimArray<Complex> > &Fframes, const int pnx, const int pny, std::vector<float> &xshifts, std::vector<float> &yshifts) {
	std::vector<Image<RFLOAT> > Iccs(n_threads);
	MultidimArray<Complex> Fref;
	std::vector<MultidimArray<Complex> > Fccs(n_threads);
	MultidimArray<RFLOAT> weight;
	std::vector<RFLOAT> cur_xshifts, cur_yshifts;
	std::vector<FourierTransformer> transformers(n_threads);
	bool converged = false;

	// Parameters TODO: make an option
	const int max_iter = 5;
	const int search_range = 50; // px TODO: bound check
	const RFLOAT tolerance = 0.5; // px
	const RFLOAT EPS = 1e-15;

	// Shifts within an iteration
	const int n_frames = xshifts.size();
	cur_xshifts.resize(n_frames);
	cur_yshifts.resize(n_frames);

	const int nfx = XSIZE(Fframes[0]), nfy = YSIZE(Fframes[0]);
	const int nfy_half = nfy / 2;
	Fref.reshape(nfy, nfx);
	for (int i = 0; i < n_threads; i++) {
		Iccs[i]().reshape(pny, pnx);
		Fccs[i].reshape(Fref);
	}

#ifdef DEBUG
	std::cout << "Patch Size X = " << pnx << " Y  = " << pny << std::endl;
	std::cout << "Fframes X = " << nfx << " Y = " << nfy << std::endl;
	std::cout << "Trajectory size: " << xshifts.size() << std::endl;
#endif

	// Initialize B factor weight
	weight.reshape(Fref);
	RCTIC(TIMING_PREP_WEIGHT);
	#pragma omp parallel for schedule(static)
	for (int y = 0; y < nfy; y++) {
		int ly = y;
		if (y > nfy_half) ly = y - nfy;
		RFLOAT ly2 = ly * (RFLOAT)ly / (nfy * (RFLOAT)nfy);

		for (int x = 0; x < nfx; x++) {
			RFLOAT dist2 = ly2 + x * (RFLOAT)x / (nfx * (RFLOAT)nfx);
			DIRECT_A2D_ELEM(weight, y, x) = exp(- 2 * dist2 * bfactor); // 2 for Fref and Fframe
		}
	}
	RCTOC(TIMING_PREP_WEIGHT);

	for (int iter = 1; iter	<= max_iter; iter++) {
		RCTIC(TIMING_MAKE_REF);
		Fref.initZeros();
		
		for (int iframe = 0; iframe < n_frames; iframe++) {
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref) {
				DIRECT_MULTIDIM_ELEM(Fref, n) += DIRECT_MULTIDIM_ELEM(Fframes[iframe], n);
			}
		}
#ifdef DEBUG
		transformers[tid].inverseFourierTransform(Fref, Icc());
		Icc.write("ref.spi");
		std::cout << "Done Fref." << std::endl;
#endif
		RCTOC(TIMING_MAKE_REF);

		#pragma omp parallel for
		for (int iframe = 0; iframe < n_frames; iframe++) {
			const int tid = omp_get_thread_num();
			
			RCTIC(TIMING_CCF_CALC);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref) {
				DIRECT_MULTIDIM_ELEM(Fccs[tid], n) = (DIRECT_MULTIDIM_ELEM(Fref, n) - DIRECT_MULTIDIM_ELEM(Fframes[iframe], n)) * 
				                                     DIRECT_MULTIDIM_ELEM(Fframes[iframe], n).conj() * DIRECT_MULTIDIM_ELEM(weight, n);
			}
			RCTOC(TIMING_CCF_CALC);

			RCTIC(TIMING_CCF_IFFT);
			transformers[tid].inverseFourierTransform(Fccs[tid], Iccs[tid]());
			RCTOC(TIMING_CCF_IFFT);
			
			RCTIC(TIMING_CCF_FIND_MAX);
			RFLOAT maxval = -9999;
			int posx, posy;
			for (int y = -search_range; y < search_range; y++) {
				int iy = y;
				if (y < 0) iy = pny + y;

				for (int x = -search_range; x < search_range; x++) {
					int ix = x;
					if (x < 0) ix = pnx + x;
					RFLOAT val = DIRECT_A2D_ELEM(Iccs[tid](), iy, ix);
//					std::cout << "(x, y) = " << x << ", " << y << ", (ix, iy) = " << ix << " , " << iy << " val = " << val << std::endl;
					if (val > maxval) {
						posx = x; posy = y;
						maxval = val;
					}
				}
			}

			int ipx_n = posx - 1, ipx = posx, ipx_p = posx + 1, ipy_n = posy - 1, ipy = posy, ipy_p = posy + 1;
			if (ipx_n < 0) ipx_n = pnx + ipx_n;
			if (ipx < 0) ipx = pnx + ipx;
			if (ipx_p < 0) ipx_p = pnx + ipx_p;
			if (ipy_n < 0) ipy_n = pny + ipy_n;
			if (ipy < 0) ipy = pny + ipy;
			if (ipy_p < 0) ipy_p = pny + ipy_p;

			// Quadratic interpolation by Jasenko
			RFLOAT vp, vn;
			vp = DIRECT_A2D_ELEM(Iccs[tid](), ipy, ipx_p);
			vn = DIRECT_A2D_ELEM(Iccs[tid](), ipy, ipx_n);
 			if (std::abs(vp + vn - 2.0 * maxval) > EPS) {
				cur_xshifts[iframe] = posx - 0.5 * (vp - vn) / (vp + vn - 2.0 * maxval);
			} else {
				cur_xshifts[iframe] = posx;
			}

			vp = DIRECT_A2D_ELEM(Iccs[tid](), ipy_p, ipx);
			vn = DIRECT_A2D_ELEM(Iccs[tid](), ipy_n, ipx);
 			if (std::abs(vp + vn - 2.0 * maxval) > EPS) {
				cur_yshifts[iframe] = posy - 0.5 * (vp - vn) / (vp + vn - 2.0 * maxval);
			} else {
				cur_yshifts[iframe] = posy;
			}
#ifdef DEBUG
			std::cout << "tid " << tid << " Frame " << 1 + iframe << ": raw shift x = " << posx << " y = " << posy << " cc = " << maxval << " interpolated x = " << cur_xshifts[iframe] << " y = " << cur_yshifts[iframe] << std::endl;
#endif
			RCTOC(TIMING_CCF_FIND_MAX);
		}

		// Set origin
		RFLOAT x_sumsq = 0, y_sumsq = 0;
		for (int iframe = n_frames - 1; iframe >= 0; iframe--) { // do frame 0 last!
			cur_xshifts[iframe] -= cur_xshifts[0];
			cur_yshifts[iframe] -= cur_yshifts[0];
			x_sumsq += cur_xshifts[iframe] * cur_xshifts[iframe];
			y_sumsq += cur_yshifts[iframe] * cur_yshifts[iframe];
		}
		cur_xshifts[0] = 0; cur_yshifts[0] = 0;

		for (int iframe = 0; iframe < n_frames; iframe++) {
			xshifts[iframe] += cur_xshifts[iframe];
			yshifts[iframe] += cur_yshifts[iframe];
//			std::cout << "Shift for Frame " << iframe << ": delta_x = " << cur_xshifts[iframe] << " delta_y = " << cur_yshifts[iframe] << std::endl;
		}

		// Apply shifts
		// Since the image is not necessarily square, we cannot use the method in fftw.cpp
		RCTIC(TIMING_FOURIER_SHIFT);
		#pragma omp parallel for
		for (int iframe = 1; iframe < n_frames; iframe++) {
			shiftNonSquareImageInFourierTransform(Fframes[iframe], -cur_xshifts[iframe] / pnx, -cur_yshifts[iframe] / pny);
		}
		RCTOC(TIMING_FOURIER_SHIFT);

		// Test convergence
		RFLOAT rmsd = std::sqrt((x_sumsq + y_sumsq) / n_frames);
		std::cout << "Iteration " << iter << ": RMSD = " << rmsd << " px" << std::endl;
		if (rmsd < tolerance) {
			converged = true;
			break;
		}
	}

#ifdef DEBUG
	for (int iframe = 0; iframe < n_frames; iframe++) {
		std::cout << iframe << " " << xshifts[iframe] << " " << yshifts[iframe] << std::endl;
	}
#endif

	return converged;
}

// dose is equivalent dose at 300 kV at the END of the frame.
// This implements the model by Timothy Grant & Nikolaus Grigorieff on eLife, 2015
// doi: 10.7554/eLife.06980
void MotioncorrRunner::doseWeighting(std::vector<MultidimArray<Complex> > &Fframes, std::vector<RFLOAT> doses) {
	const int nfx = XSIZE(Fframes[0]), nfy = YSIZE(Fframes[0]);
	const int nfy_half = nfy / 2;
	const RFLOAT nfy2 = (RFLOAT)nfy * nfy;
	const RFLOAT nfx2 = (RFLOAT)(nfx - 1) * (nfx - 1) * 4; // assuming nx is even
	const int n_frames= Fframes.size();
	const RFLOAT A = 0.245, B = -1.665, C = 2.81;

	#pragma omp parallel for schedule(static)
	for (int y = 0; y < nfy; y++) {
		int ly = y;
		if (y > nfy_half) ly = y - nfy;

		const RFLOAT ly2 = (RFLOAT)ly * ly / nfy2;
		for (int x = 0; x < nfx; x++) {
			const RFLOAT dinv2 = ly2 + (RFLOAT)x * x / nfx2;
			const RFLOAT dinv = std::sqrt(dinv2) / angpix; // d = N * angpix / dist, thus dinv = dist / N / angpix
			const RFLOAT Ne = (A * std::pow(dinv, B) + C) * 2; // Eq. 3. 2 comes from Eq. 5
			RFLOAT sum_weight_sq = 0;

			for (int iframe = 0; iframe < n_frames; iframe++) {
				const RFLOAT weight = std::exp(- doses[iframe] / Ne); // Eq. 5. 0.5 is factored out to Ne.
				if (isnan(weight)) {
					std::cout << "dose = " <<  doses[iframe] << " Ne = " << Ne << " frm = " << iframe << " lx = " << x << " ly = " << ly << " reso = " << 1 / dinv << " weight = " << weight << std::endl;
				}
				sum_weight_sq += weight * weight;
				DIRECT_A2D_ELEM(Fframes[iframe], y, x) *= weight;
			}

			sum_weight_sq = std::sqrt(sum_weight_sq);
			if (isnan(sum_weight_sq)) {
				std::cout << " Ne = " << Ne << " lx = " << x << " ly = " << ly << " reso = " << 1 / dinv << " sum_weight_sq NaN" << std::endl;
				REPORT_ERROR("Shouldn't happen.");
			}
			for (int iframe = 0; iframe < n_frames; iframe++) {
				DIRECT_A2D_ELEM(Fframes[iframe], y, x) /= sum_weight_sq; // Eq. 9
			}
		}
	}
}

// shiftx, shifty is relative to the (real space) image size
void MotioncorrRunner::shiftNonSquareImageInFourierTransform(MultidimArray<Complex> &frame, RFLOAT shiftx, RFLOAT shifty) {
	const int nfx = XSIZE(frame), nfy = YSIZE(frame);
	const int nfy_half = nfy / 2;
	RFLOAT twoPI = 2 * PI;

	for (int y = 0; y < nfy; y++) {
		int ly = y;
		if (y > nfy_half) ly = y - nfy;

		for (int x = 0; x < nfx; x++) {
			RFLOAT phase_shift = twoPI * (x * shiftx + ly * shifty);
			RFLOAT a, b, c, d, ac, bd, ab_cd;
#ifdef RELION_SINGLE_PRECISION
			SINCOSF(phase_shift, &b, &a);
#else
			SINCOS(phase_shift, &b, &a);
#endif
			c = DIRECT_A2D_ELEM(frame, y, x).real;
			d = DIRECT_A2D_ELEM(frame, y, x).imag;
			ac = a * c;
			bd = b * d;
			ab_cd = (a + b) * (c + d); // (ab_cd-ac-bd = ad+bc : but needs 4 multiplications)
			DIRECT_A2D_ELEM(frame, y, x) = Complex(ac - bd, ab_cd - ac - bd);	
		}
	}
}

// Iwork is overwritten
void MotioncorrRunner::binNonSquareImage(Image<RFLOAT> &Iwork, RFLOAT bin_factor) {
	FourierTransformer transformer;

	const int nx = XSIZE(Iwork()), ny = YSIZE(Iwork());
	int new_nx = nx / bin_factor, new_ny = ny / bin_factor;
	new_nx -= new_nx % 2; new_ny -= new_ny % 2; // force it to be even
//	std::cout << "Binning from X = " << nx << " Y = " << ny << " to X = " << new_nx << " Y = " << new_ny << std::endl;

	const int half_new_ny = new_ny / 2, new_nfx = new_nx / 2 + 1;
	MultidimArray<Complex> Fref, Fbinned(new_ny, new_nfx);
	transformer.FourierTransform(Iwork(), Fref);

	for (int y = 0; y <= half_new_ny; y++) {
		for (int x = 0; x < new_nfx; x++) {
			DIRECT_A2D_ELEM(Fbinned, y, x) =  DIRECT_A2D_ELEM(Fref, y, x);
		}
	}
	for (int y = half_new_ny + 1; y < new_ny; y++) {
		for (int x = 0; x < new_nfx; x++) {
			DIRECT_A2D_ELEM(Fbinned, y, x) =  DIRECT_A2D_ELEM(Fref, ny - 1 - new_ny + y, x);
		}
	}

	Iwork().reshape(new_ny, new_nx);
	transformer.inverseFourierTransform(Fbinned, Iwork());
}

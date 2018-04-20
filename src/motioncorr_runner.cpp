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
#include "src/acc/cuda/cuda_mem_utils.h"
#endif
#include "src/micrograph_model.h"
#include "src/matrix2d.h"
#include "src/matrix1d.h"
#include <omp.h>

//#define TIMING
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
	int TIMING_CLIP_PATCH = timer.setNew("prep patch - clip (in thread)");
	int TIMING_PATCH_FFT = timer.setNew("prep patch - FFT (in thread)");
	int TIMING_PATCH_ALIGN = timer.setNew("patch alignment");
	int TIMING_PREP_WEIGHT = timer.setNew("align - prep weight");
	int TIMING_MAKE_REF = timer.setNew("align - make reference");
	int TIMING_CCF_CALC = timer.setNew("align - calc CCF (in thread)");
	int TIMING_CCF_IFFT = timer.setNew("align - iFFT CCF (in thread)");
	int TIMING_CCF_FIND_MAX = timer.setNew("align - argmax CCF (in thread)");
	int TIMING_FOURIER_SHIFT = timer.setNew("align - shift in Fourier space");
	int TIMING_FIT_POLYNOMIAL = timer.setNew("fit polynomial");
	int TIMING_DOSE_WEIGHTING = timer.setNew("dose weighting");
	int TIMING_DW_WEIGHT = timer.setNew("dw - calc weight");
	int TIMING_DW_IFFT = timer.setNew("dw - iFFT");
	int TIMING_REAL_SPACE_INTERPOLATION = timer.setNew("real space interpolation");
	int TIMING_BINNING = timer.setNew("binning");
//	int TIMING_ = timer.setNew("");

#else
	#define RCTIC(label)
	#define RCTOC(label)
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

	int doseweight_section = parser.addSection("Dose-weighting options");
	do_dose_weighting = parser.checkOption("--dose_weighting", "Use MOTIONCOR2s or UNBLURs dose-weighting scheme");
	voltage = textToFloat(parser.getOption("--voltage","Voltage (in kV) for dose-weighting inside MOTIONCOR2/UNBLUR","300"));
	dose_per_frame = textToFloat(parser.getOption("--dose_per_frame", "Electron dose (in electrons/A2/frame) for dose-weighting inside MOTIONCOR2/UNBLUR", "0"));
	pre_exposure = textToFloat(parser.getOption("--preexposure", "Pre-exposure (in electrons/A2) for dose-weighting inside UNBLUR", "0"));

	do_own = parser.checkOption("--use_own", "Use our own implementation of motion correction");
	save_noDW = parser.checkOption("--save_noDW", "Save aligned but non dose weighted micrograph. (This is always ON for MOTIONCOR2)");
	interpolate_shifts = parser.checkOption("--interpolate_shifts", "(EXPERIMENTAL) Interpolate shifts");
	ccf_downsample = textToFloat(parser.getOption("--ccf_downsample", "(EXPERT) Downsampling rate of CC map. default = 0 = automatic based on B factor", "0"));
	early_binning = parser.checkOption("--early_binning", "(EXPERT) Do binning before alignment to reduce memory usage. This might dampen signal near Nyquist.");
	if (ccf_downsample > 1) REPORT_ERROR("--ccf_downsample cannot exceed 1.");
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
		if (fn_defect != "") {
			std::cerr << "WARNING: The defect file is not yet supported in our own implementation of motion correction. The defect file is ignored." << std::endl;
		}
		if (patch_x <= 2 || patch_y <= 2) {
			std::cerr << "The number of patches is too small (<= 2). Patch based alignment will be skipped." << std::endl;
		}
		if (group > 1) {
			std::cerr << "WARNING: Frame grouping is under development in our own implementation of motion correction." << std::endl;
		}
		if (do_save_movies) {
			std::cerr << "WARNING: In our own implementation of motion correction, we do not save aligned movies. --save_movies was ignored." << std::endl;
			do_save_movies = false;
		}
	} else {
		REPORT_ERROR(" ERROR: You have to specify which programme to use through either --use_motioncor2 or --use_unblur");
	}

	if (do_dose_weighting)
	{
		if (do_motioncor2) save_noDW = true; // MOTIONCOR2 always writes non-dose weighted micrographs
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
		else if (do_motioncor2)
			std::cout << " Using MOTIONCOR2 executable in: " << fn_motioncor2_exe << std::endl;
		else if (do_own)
			std::cout << " Using our own implementation based on MOTIONCOR2 algorithm" << std::endl;
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
		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);


		Micrograph mic(fn_micrographs[imic], fn_gain_reference, bin_factor);

		bool result = false;
		if (do_own)
			result = executeOwnMotionCorrection(mic);
		else if (do_unblur)
			result = executeUnblur(mic);
		else if (do_motioncor2)
			result = executeMotioncor2(mic);
		else
			REPORT_ERROR("Bug: by now it should be clear whether to use MotionCor2 or Unblur...");

		if (result) {
			saveModel(mic);
			plotShifts(fn_micrographs[imic], mic);
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

bool MotioncorrRunner::executeMotioncor2(Micrograph &mic, int rank)
{
	FileName fn_mic = mic.getMovieFilename();
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
		Ihead.read(fn_mic, false, -1, false, true); // select_img -1, mmap false, is_2D true
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
		getShiftsMotioncor2(fn_out, mic);
	}

	// Success!
	return true;
}

void MotioncorrRunner::getShiftsMotioncor2(FileName fn_log, Micrograph &mic)
{
	std::ifstream in(fn_log.data(), std::ios_base::in);
	if (in.fail())
		return;

	std::string line;

	// Start reading the ifstream at the top
	in.seekg(0);

	// Read through the shifts file
	int frame = first_frame_sum; // 1-indexed
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
 		   		mic.setGlobalShift(frame, textToFloat(words[0]), textToFloat(words[1]));
				frame++;
			}
			else
			{
				// Stop now
				break;
			}
		}
	}
	in.close();

	mic.first_frame = first_frame_sum;

	// Read local shifts
	FileName fn_patch = fn_log.withoutExtension() + "0-Patch-Patch.log";
	if (!exists(fn_patch)) {
		std::cerr << "warning: failed to load local shifts from MotionCor2 logfile, but this does not affect further processing: " << fn_patch << std::endl;
		return;
	}

	in.open(fn_patch.c_str(), std::ios_base::in);
	in.seekg(0);

	while (getline(in, line, '\n'))
	{
		if (line.find("#") != std::string::npos) continue;

		std::vector<std::string> words;
		tokenize(line, words);
    		if (words.size() != 7) continue;

		mic.patchZ.push_back(textToFloat(words[0]) + first_frame_sum - 1); // 1-indexed
		mic.patchX.push_back(textToFloat(words[1]));
		mic.patchY.push_back(textToFloat(words[2]));
		mic.localShiftX.push_back(textToFloat(words[3]));
		mic.localShiftY.push_back(textToFloat(words[4]));
		mic.localFitX.push_back(textToFloat(words[5]));
		mic.localFitY.push_back(textToFloat(words[6]));
	}
}

bool MotioncorrRunner::executeUnblur(Micrograph &mic)
{
	FileName fn_mic = mic.getMovieFilename();
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
	Itest.read(fn_mic, false, -1, false, true); // select_img -1, mmap false, is_2D true
	int Nframes = NSIZE(Itest());

	std::ofstream  fh;
	fh.open((fn_com).c_str(), std::ios::out);
	if (!fh)
	 REPORT_ERROR( (std::string)"executeUnblur cannot create file: " + fn_com);

	// Write script to run Unblur
	fh << "#!/usr/bin/env csh"<<std::endl;
	fh << "setenv  OMP_NUM_THREADS " << integerToString(n_threads)<<std::endl;
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
	getShiftsUnblur(fn_shifts, mic);

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
		fh2 << "setenv  OMP_NUM_THREADS " << integerToString(n_threads)<<std::endl;
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

void MotioncorrRunner::getShiftsUnblur(FileName fn_shifts, Micrograph &mic)
{

	std::ifstream in(fn_shifts.data(), std::ios_base::in);
	if (in.fail())
		return;

	std::vector<RFLOAT> xshifts, yshifts;
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

	for (int i = 0, ilim = xshifts.size(); i < ilim; i++) {
		int frame = i + first_frame_sum; // make 1-indexed
		// shifts from Unblur are in angstrom not pixel
		mic.setGlobalShift(frame, xshifts[i] / angpix, yshifts[i] / angpix);
	}
	mic.first_frame = first_frame_sum;
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
void MotioncorrRunner::plotShifts(FileName fn_mic, Micrograph &mic)
{
	const RFLOAT SCALE = 40;
	RFLOAT shift_scale = SCALE;
	if (do_unblur) {
		shift_scale = 1.0 / angpix; // convert from A to pix, no scaling
	}

	// Global shift
	FileName fn_eps = fn_out + fn_mic.withoutExtension() + "_shifts.eps";
	CPlot2D *plot2D=new CPlot2D(fn_eps);
 	plot2D->SetXAxisSize(600);
 	plot2D->SetYAxisSize(600);
	plot2D->SetDrawLegend(false);

 	CDataSet dataSet;
	dataSet.SetDrawMarker(false);
	dataSet.SetDatasetColor(0.0,0.0,1.0);
	RFLOAT xshift, yshift;
	// UNBLUR does not provide local trajectories, so start the global trajectory from the origin.
	const RFLOAT xcenter = (!do_unblur) ? mic.getWidth() / 2.0 : 0;
	const RFLOAT ycenter = (!do_unblur) ? mic.getHeight() / 2.0 : 0;
	for (int j = mic.first_frame, jlim = mic.getNframes(); j <= jlim; j++) // 1-indexed
	{
		mic.getShiftAt(j, 0, 0, xshift, yshift, false); // no local model
		if (xshift != Micrograph::NOT_OBSERVED) {
			CDataPoint point(xcenter + shift_scale * xshift, ycenter + shift_scale * yshift);
			dataSet.AddDataPoint(point);
		}
	}
	plot2D->AddDataSet(dataSet);

	// Mark starting point
	mic.getShiftAt(mic.first_frame, 0, 0, xshift, yshift, false);
	if (xshift != Micrograph::NOT_OBSERVED) {
		CDataSet dataSetStart;
		dataSetStart.SetDrawMarker(true);
		dataSetStart.SetMarkerSize(5);
		dataSetStart.SetDatasetColor(1.0,0.0,0.0);
		CDataPoint point2(xcenter + shift_scale * xshift, ycenter + shift_scale * yshift);
		dataSetStart.AddDataPoint(point2);
		plot2D->AddDataSet(dataSetStart);
	}

	// Local shifts
	const int n_local = mic.localShiftX.size();
	int i = 0;
	while (i < n_local) {
		CDataSet patch_start;
		patch_start.SetDrawMarker(true);
		patch_start.SetMarkerSize(5);
		patch_start.SetDatasetColor(1.0,0.0,0.0);

		CDataSet fit, obs;
		const int start = i;
		fit.SetDrawMarker(false);
		fit.SetDatasetColor(0.0,0.0,0.0);
		obs.SetDrawMarker(false);
		obs.SetDatasetColor(0.5,0.5,0.5);
		while (i < n_local) {
		//	std::cout << mic.patchX[i] << " " << mic.patchY[i] << " " << mic.patchZ[i] << " " << mic.localShiftX[i] << " " << mic.localShiftY[i] << std::endl;
			if ((mic.patchX[start] != mic.patchX[i]) || (mic.patchY[start] != mic.patchY[i])) {
		//		std::cout << "End of a trace" << std::endl;
				break; // start again from this frame
			}
			CDataPoint p_obs(mic.patchX[start] + shift_scale * mic.localShiftX[i], mic.patchY[start] + shift_scale * mic.localShiftY[i]);
			obs.AddDataPoint(p_obs);
			CDataPoint p_fit(mic.patchX[start] + shift_scale * mic.localFitX[i], mic.patchY[start] + shift_scale * mic.localFitY[i]);
			fit.AddDataPoint(p_fit);
			if (i == start) patch_start.AddDataPoint(p_fit);
			i++;
		}
		plot2D->AddDataSet(obs);
		plot2D->AddDataSet(fit);

		plot2D->AddDataSet(patch_start);
	}

	char title[256];
	snprintf(title, 255, "X (in pixels; trajectory scaled by %.0f)", shift_scale);
	plot2D->SetXAxisTitle(title);
	title[0] = 'Y';
	plot2D->SetYAxisTitle(title);

	plot2D->OutputPostScriptPlot(fn_eps);
}

void MotioncorrRunner::saveModel(Micrograph &mic) {
	if (do_dose_weighting) {
		mic.angpix = angpix;
		mic.voltage = voltage;
		mic.dose_per_frame = dose_per_frame;
		mic.pre_exposure = pre_exposure;
	}

	FileName fn_avg, fn_mov;
	getOutputFileNames(mic.getMovieFilename(), fn_avg, fn_mov);

	mic.write(fn_avg.withoutExtension() + ".star");
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
			if (do_dose_weighting && save_noDW)
			{
				FileName fn_avg_wodose = fn_avg.withoutExtension() + "_noDW.mrc";
                MDavg.setValue(EMDL_MICROGRAPH_NAME_WODOSE, fn_avg_wodose);
			}
			MDavg.setValue(EMDL_MICROGRAPH_NAME, fn_avg);
			MDavg.setValue(EMDL_MICROGRAPH_METADATA_NAME, fn_avg.withoutExtension() + ".star");
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

bool MotioncorrRunner::executeOwnMotionCorrection(Micrograph &mic) {
	FileName fn_mic = mic.getMovieFilename();
	FileName fn_avg, fn_mov;
	getOutputFileNames(fn_mic, fn_avg, fn_mov);
	FileName fn_avg_noDW = fn_avg.withoutExtension() + "_noDW.mrc";
	FileName fn_log = fn_avg.withoutExtension() + ".log";
	std::ofstream logfile;
	logfile.open(fn_log);

	logfile << "Working on " << fn_mic << " with " << n_threads << " thread(s)." << std::endl << std::endl;

	Image<float> Ihead, Igain, Iref;
	std::vector<MultidimArray<fComplex> > Fframes;
	std::vector<Image<float> > Iframes;
	std::vector<int> frames;

	RFLOAT prescaling = 1;

	const int hotpixel_sigma = 6;
	const int fit_rmsd_threshold = 10; // px

	// Check image size
	Ihead.read(fn_mic, false, -1, false, true); // select_img -1, mmap false, is_2D true
	int nx = XSIZE(Ihead()), ny = YSIZE(Ihead()), nn = NSIZE(Ihead());

	// Which frame to use?
	logfile << "Movie size: X = " << nx << " Y = " << ny << " N = " << nn << std::endl;
	logfile << "Frames to be used:";
	for (int i = 0; i < nn; i++) {
		// For users, all numbers are 1-indexed. Internally they are 0-indexed.
		int frame = i + 1;
		if (frame < first_frame_sum) continue;
		if (last_frame_sum > 0 && frame > last_frame_sum) continue;
		frames.push_back(i);
		logfile << " " << frame;
	}
	logfile << std::endl;

	const int n_frames = frames.size();
	Iframes.resize(n_frames);
	Fframes.resize(n_frames);

	std::vector<RFLOAT> xshifts(n_frames), yshifts(n_frames);

	// Setup grouping
	logfile << "Frame grouping: n_frames = " << n_frames << ", requested group size = " << group << std::endl;
	const int n_groups = n_frames / group;
	if (n_groups < 3) REPORT_ERROR("Too few frames (< 3) after grouping.");
	int n_remaining = n_frames % group;
	std::vector<int> group_start(n_groups, 0), group_size(n_groups, group);
	while (n_remaining > 0) {
		for (int i = n_groups - 1; i >= 1 && n_remaining > 0; i--) {
			// Do not expand the first group, where the motion is largest.
			group_size[i]++;
			n_remaining--;
		}
	}
	for (int i = 1; i < n_groups; i++) {
		group_start[i] = group_start[i - 1] + group_size[i - 1];
	}
	logfile << " | ";
	for (int i = 0, igroup = 0; i < n_frames; i++) {
		logfile << frames[i] + 1 << " "; // make 1-indexed
		if (i == group_start[igroup] + group_size[igroup] - 1) {
			logfile << "| ";
			igroup++;
		}
	}
	logfile << std::endl;
	if (n_frames % group != 0) {
		logfile << "Some groups contain more than the requested number of frames (" << group << ") because the number of frames (" << n_frames << ") was not divisible." << std::endl;
		logfile << "If you want to ignore remaining frame(s) instead, use --last_frame_sum to discard last frame(s)." << std::endl;
	}
	logfile << "interpolate_shifts = " << interpolate_shifts << std::endl;
	logfile << std::endl;

	// Read gain reference
	RCTIC(TIMING_READ_GAIN);
	if (fn_gain_reference != "") {
		Igain.read(fn_gain_reference);
		if (XSIZE(Igain()) != nx || YSIZE(Igain()) != ny) {
			std::cerr << "fn_mic: " << fn_mic << std::endl;
			REPORT_ERROR("The size of the image and the size of the gain reference do not match. Make sure the gain reference has been rotated if necessary.");
		}
	}
	RCTOC(TIMING_READ_GAIN);

	// Read images
	RCTIC(TIMING_READ_MOVIE);
	#pragma omp parallel for num_threads(n_threads)
	for (int iframe = 0; iframe < n_frames; iframe++) {
		Iframes[iframe].read(fn_mic, true, frames[iframe], false, true); // mmap false, is_2D true
	}
	RCTOC(TIMING_READ_MOVIE);

	// Apply gain
	RCTIC(TIMING_APPLY_GAIN);
	if (fn_gain_reference != "") {
		for (int iframe = 0; iframe < n_frames; iframe++) {
			#pragma omp parallel for num_threads(n_threads)
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Igain()) {
				DIRECT_MULTIDIM_ELEM(Iframes[iframe](), n) *= DIRECT_MULTIDIM_ELEM(Igain(), n);
			}
		}
	}
	RCTOC(TIMING_APPLY_GAIN);

	/*
	MultidimArray<float> Iframe(ny, nx);
	Iframe.initZeros();

	// Simple unaligned sum
	RCTIC(TIMING_INITIAL_SUM);
	for (int iframe = 0; iframe < n_frames; iframe++) {
		#pragma omp parallel for num_threads(n_threads)
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
	logfile << "Mean = " << mean << " Std = " << std << ", Hotpixel threshold = " << threshold << std::endl;

	MultidimArray<bool> bBad(ny, nx);
	bBad.initZeros();
	int n_bad = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iframe) {
		if (DIRECT_MULTIDIM_ELEM(Iframe, n) > threshold) {
			DIRECT_MULTIDIM_ELEM(bBad, n) = true;
			n_bad++;
		}
	}
	logfile << "Detected " << n_bad << " hot pixels to be corrected. (BUT correction not implemented yet)" << std::endl << std::endl;
	RCTOC(TIMING_DETECT_HOT);
	// TODO: fix defects
	*/

	if (early_binning) {
		nx /= bin_factor; ny /= bin_factor;
		if (nx % 2 != 0 || ny % 2 != 0) {
			// Do not adjust the size because it might lead to non-square pixels
			REPORT_ERROR("The dimensions of the image after binning must be even");
		}
		prescaling = bin_factor;
		logfile << "Image size after binning: X = " << nx << " Y = " << ny << std::endl;
	}

	// FFT
	RCTIC(TIMING_GLOBAL_FFT);
	#pragma omp parallel for num_threads(n_threads)
	for (int iframe = 0; iframe < n_frames; iframe++) {
		if (!early_binning) {
			NewFFT::FourierTransform(Iframes[iframe](), Fframes[iframe]);
		} else {
			MultidimArray<fComplex> Fframe;
			Fframes[iframe].reshape(ny, nx / 2 + 1);
			cropInFourierSpace(Fframe, Fframes[iframe]);
		}
		Iframes[iframe].clear(); // save some memory (global alignment use the most memory)
	}
	RCTOC(TIMING_GLOBAL_FFT);

	// Global alignment
	// TODO: Consider frame grouping in global alignment.
	logfile << std::endl << "Global alignment:" << std::endl;
	RCTIC(TIMING_GLOBAL_ALIGNMENT);
	alignPatch(Fframes, nx, ny, bfactor / (prescaling * prescaling), xshifts, yshifts, logfile);
	RCTOC(TIMING_GLOBAL_ALIGNMENT);
	for (int i = 0, ilim = xshifts.size(); i < ilim; i++) {
		// Should be in the original pixel size
		mic.setGlobalShift(frames[i] + 1, xshifts[i] * prescaling, yshifts[i] * prescaling); // 1-indexed
        }

	Iref().reshape(ny, nx);
	Iref().initZeros();
	RCTIC(TIMING_GLOBAL_IFFT);
	#pragma omp parallel for num_threads(n_threads)
	for (int iframe = 0; iframe < n_frames; iframe++) {
		Iframes[iframe]().reshape(ny, nx);
		NewFFT::inverseFourierTransform(Fframes[iframe], Iframes[iframe]());
		// Unfortunately, we cannot deallocate Fframes here because of dose-weighting
	}
	RCTOC(TIMING_GLOBAL_IFFT);

	// Patch based alignment
	logfile << std::endl << "Local alignments:" << std::endl;
	logfile << "Patches: X = " << patch_x << " Y = " << patch_y << std::endl;
	bool do_local = (patch_x > 2) && (patch_y > 2);
	if (!do_local) {
		logfile << "Too few patches to do local alignments." << std::endl;
	}

	if (do_local) {
		const int patch_nx = nx / patch_x, patch_ny = ny / patch_y, n_patches = patch_x * patch_y;
		std::vector<RFLOAT> patch_xshifts, patch_yshifts, patch_frames, patch_xs, patch_ys;
		std::vector<MultidimArray<fComplex> > Fpatches(n_groups);

		int ipatch = 1;
		for (int iy = 0; iy < patch_y; iy++) {
			for (int ix = 0; ix < patch_x; ix++) {
				int x_start = ix * patch_nx, y_start = iy * patch_ny; // Inclusive
				int x_end = x_start + patch_nx, y_end = y_start + patch_ny; // Exclusive
				if (x_end > nx) x_end = nx;
				if (y_end > ny) y_end = ny;
				// make patch size even
				if ((x_end - x_start) % 2 == 1) {
					if (x_end == nx) x_start++;
					else x_end--;
				}
				if ((y_end - y_start) % 2 == 1) {
					if (y_end == ny) y_start++;
					else y_end--;
				}

				int x_center = (x_start + x_end - 1) / 2, y_center = (y_start + y_end - 1) / 2;
				logfile << "Patch (" << iy + 1 << ", " << ix + 1 << "): " << ipatch << " / " << patch_x * patch_y;
				logfile << ", X range = [" << x_start << ", " << x_end << "), Y range = [" << y_start << ", " << y_end << ")";
				logfile << ", Center = (" << x_center << ", " << y_center << ")" << std::endl;
				ipatch++;

				std::vector<RFLOAT> local_xshifts(n_groups), local_yshifts(n_groups);
				RCTIC(TIMING_PREP_PATCH);
				std::vector<MultidimArray<float> >Ipatches(n_threads);
				#pragma omp parallel for num_threads(n_threads)
				for (int igroup = 0; igroup < n_groups; igroup++) {
					const int tid = omp_get_thread_num();
					Ipatches[tid].reshape(y_end - y_start, x_end - x_start); // end is not included
					RCTIC(TIMING_CLIP_PATCH);
					for (int iframe = group_start[igroup]; iframe < group_start[igroup] + group_size[igroup]; iframe++) {
						for (int ipy = y_start; ipy < y_end; ipy++) {
							for (int ipx = x_start; ipx < x_end; ipx++) {
								DIRECT_A2D_ELEM(Ipatches[tid], ipy - y_start, ipx - x_start) = DIRECT_A2D_ELEM(Iframes[iframe](), ipy, ipx);
							}
						}
					}
					RCTOC(TIMING_CLIP_PATCH);

					RCTIC(TIMING_PATCH_FFT);
					NewFFT::FourierTransform(Ipatches[tid], Fpatches[igroup]);
					RCTOC(TIMING_PATCH_FFT);
				}
				RCTOC(TIMING_PREP_PATCH);

				RCTIC(TIMING_PATCH_ALIGN);
				bool converged = alignPatch(Fpatches, x_end - x_start, y_end - y_start, bfactor / (prescaling * prescaling), local_xshifts, local_yshifts, logfile);
				RCTOC(TIMING_PATCH_ALIGN);
				if (!converged) continue;

				std::vector<RFLOAT> interpolated_xshifts(n_frames), interpolated_yshifts(n_frames);
				interpolateShifts(group_start, group_size, local_xshifts, local_yshifts, n_frames, interpolated_xshifts, interpolated_yshifts);
				if (interpolate_shifts) {
					// Recenter to the first frame
					for (int iframe = 0; iframe < n_frames; iframe++) {
						interpolated_xshifts[iframe] -= interpolated_xshifts[0];
						interpolated_yshifts[iframe] -= interpolated_yshifts[0];
					}
					// Store shifts
					for (int iframe = 0; iframe < n_frames; iframe++) {
						patch_xshifts.push_back(interpolated_xshifts[iframe]);
						patch_yshifts.push_back(interpolated_yshifts[iframe]);
						patch_frames.push_back(iframe);
						patch_xs.push_back(x_center);
						patch_ys.push_back(y_center);
					}
				} else { // only recenter to the center
					for (int igroup = 0; igroup < n_groups; igroup++) {
						patch_xshifts.push_back(local_xshifts[igroup] - interpolated_xshifts[0]);
						patch_yshifts.push_back(local_yshifts[igroup] - interpolated_yshifts[0]);
						RFLOAT middle_frame = group_start[igroup] + group_size[igroup] / 2.0;
						patch_frames.push_back(middle_frame);
						patch_xs.push_back(x_center);
						patch_ys.push_back(y_center);
					}
				}
			}
		}
		Fpatches.clear();

		// Fit polynomial model

		RCTIC(TIMING_FIT_POLYNOMIAL);
		const int n_obs = patch_frames.size();
		const int n_params = 18;

		if (n_obs <= n_params) {
			std::cerr << fn_mic << ": too few valid local trajectories to fit local motion model." << std::endl;
			mic.model = NULL;
			goto skip_fitting; // TODO: Refactor!
		}

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

#ifdef DEBUG_OWN
		std::cout << "Polynomial fitting coefficients for X and Y:" << std::endl;
		for (int i = 0; i < n_params; i++) {
			std::cout << i << " " << coeffX(i) << " " << coeffY(i) << std::endl;
		}
#endif

		ThirdOrderPolynomialModel *model = new ThirdOrderPolynomialModel();
		model->coeffX = coeffX; model->coeffY = coeffY;
		mic.model = model;

#ifdef DEBUG_OWN
		std::cout << "Polynomial Fitting:" << std::endl;
#endif
		RFLOAT rms_x = 0, rms_y = 0;
	        for (int i = 0; i < n_obs; i++) {
			RFLOAT x_fitted, y_fitted;
			const RFLOAT x = patch_xs[i] / nx - 0.5;
			const RFLOAT y = patch_ys[i] / ny - 0.5;
			const RFLOAT z = patch_frames[i];

			model->getShiftAt(z, x, y, x_fitted, y_fitted);
			rms_x += (patch_xshifts[i] - x_fitted) * (patch_xshifts[i] - x_fitted);
			rms_y += (patch_yshifts[i] - y_fitted) * (patch_yshifts[i] - y_fitted);

			// These shifts and RMSDs are for reporting, so should be in the original pixel size
			mic.patchX.push_back(patch_xs[i] * prescaling);
			mic.patchY.push_back(patch_ys[i] * prescaling);
			mic.patchZ.push_back(z + first_frame_sum); // 1-indexed
			mic.localShiftX.push_back(patch_xshifts[i] * prescaling);
			mic.localShiftY.push_back(patch_yshifts[i] * prescaling);
			mic.localFitX.push_back(x_fitted * prescaling);
			mic.localFitY.push_back(y_fitted * prescaling);

#ifdef DEBUG_OWN
			std::cout << " x = " << x << " y = " << y << " z = " << z;
			std::cout << ", Xobs = " << patch_xshifts[i] * prescaling << " Xfit = " << x_fitted * prescaling;
			std::cout << ", Yobs = " << patch_yshifts[i] * prescaling << " Yfit = " << y_fitted * prescaling<< std::endl;
#endif
		}
		rms_x = std::sqrt(rms_x / n_obs) * prescaling; rms_y = std::sqrt(rms_y / n_obs) * prescaling;
		logfile << std::endl << "Polynomial fit RMSD: X = " << rms_x << " px Y = " << rms_y << " px" << std::endl;
		if (rms_x >= fit_rmsd_threshold || rms_y >= fit_rmsd_threshold) {
			logfile << "The polynomial motion model did not explain the observation very well." << std::endl;
			logfile << "Local correction is disabled for this micrograph." << std::endl;
			delete mic.model;
			mic.model = NULL;

			// remove fitted trajectories
			for (int i = 0, ilim = mic.localFitX.size(); i < ilim; i++) {
				mic.localFitX[i] = 0;
				mic.localFitY[i] = 0;
			}
		}
		RCTOC(TIMING_FIT_POLYNOMIAL);
	} else { // !do_local
		mic.model = NULL;
	}

skip_fitting:
	if (!do_dose_weighting || save_noDW) {
		Iref().initZeros(Iframes[0]());

		RCTIC(TIMING_REAL_SPACE_INTERPOLATION);
		logfile << "Summing frames before dose weighting: ";
		realSpaceInterpolation(Iref, Iframes, mic.model, logfile);
		logfile << " done" << std::endl;
		RCTOC(TIMING_REAL_SPACE_INTERPOLATION);

		// Apply binning
		RCTIC(TIMING_BINNING);
		if (!early_binning && bin_factor != 1) {
			binNonSquareImage(Iref, bin_factor);
		}
		RCTOC(TIMING_BINNING);

		// Final output.
		Iref.write(!do_dose_weighting ? fn_avg : fn_avg_noDW);
		logfile << "Written aligned but non-dose weighted sum to " << fn_avg_noDW << std::endl;
	}

	// Dose weighting
	if (do_dose_weighting) {
		RCTIC(TIMING_DOSE_WEIGHTING);
		if (std::abs(voltage - 300) > 2 && std::abs(voltage - 200) > 2) {
			REPORT_ERROR("Sorry, dose weighting is supported only for 300 kV or 200 kV");
		}

		std::vector <RFLOAT> doses(n_frames);
		for (int iframe = 0; iframe < n_frames; iframe++) {
			// dose AFTER each frame.
			doses[iframe] = pre_exposure + dose_per_frame * (frames[iframe] + 1);
			if (std::abs(voltage - 200) <= 5) {
				doses[iframe] /= 0.8; // 200 kV electron is more damaging.
			}
		}

		RCTIC(TIMING_DW_WEIGHT);
		doseWeighting(Fframes, doses, angpix * prescaling);
		RCTOC(TIMING_DW_WEIGHT);

		// Update real space images
		RCTIC(TIMING_DW_IFFT);
		#pragma omp parallel for num_threads(n_threads)
		for (int iframe = 0; iframe < n_frames; iframe++) {
			NewFFT::inverseFourierTransform(Fframes[iframe], Iframes[iframe]());
		}
		RCTOC(TIMING_DW_IFFT);
		RCTOC(TIMING_DOSE_WEIGHTING);

		Iref().initZeros(Iframes[0]());
		RCTIC(TIMING_REAL_SPACE_INTERPOLATION);
		logfile << "Summing frames after dose weighting: ";
		realSpaceInterpolation(Iref, Iframes, mic.model, logfile);
		logfile << " done" << std::endl;
		RCTOC(TIMING_REAL_SPACE_INTERPOLATION);

		// Apply binning
		RCTIC(TIMING_BINNING);
		if (!early_binning & bin_factor != 1) {
			binNonSquareImage(Iref, bin_factor);
		}
		RCTOC(TIMING_BINNING);

		// Final output
		Iref.write(fn_avg);
		logfile << "Written aligned and dose-weighted sum to " << fn_avg << std::endl;
	}

	// Set the start frame for the local motion model.
	mic.first_frame = frames[0] + 1; // NOTE that this is 1-indexed.

	return true;
}

void MotioncorrRunner::interpolateShifts(std::vector<int> &group_start, std::vector<int> &group_size,
                                         std::vector<RFLOAT> &xshifts, std::vector<RFLOAT> &yshifts,
                                         int n_frames,
                                         std::vector<RFLOAT> &interpolated_xshifts, std::vector<RFLOAT> &interpolated_yshifts) {
	const int n_groups = group_start.size();
	if (n_groups < 2) REPORT_ERROR("Assert failed for n_groups >= 2");
	std::vector<RFLOAT> centers(n_groups);

	// Calculate anchor points
	for (int igroup = 0; igroup < n_groups; igroup++) {
		centers[igroup] = group_start[igroup] + group_size[igroup] / 2.0;
#ifdef DEBUG_OWN
		std::cout << "igroup = " << igroup << " center " << centers[igroup] << " x " << xshifts[igroup] << " y " << yshifts[igroup] << std::endl;
#endif
	}

	// Inter-/Extra-polate
	int cur_group = 0;
	for (int iframe = 0; iframe < n_frames; iframe++) {
		if (cur_group < n_groups - 2 && iframe >= centers[cur_group + 1]) cur_group++; // don't go to the last group

		// This formula can be used for both interpolation and extrapolation
		interpolated_xshifts[iframe] = (xshifts[cur_group] * (centers[cur_group + 1] - iframe) +
		                                xshifts[cur_group + 1] * (iframe - centers[cur_group])) /
		                               (centers[cur_group + 1] - centers[cur_group]);
		interpolated_yshifts[iframe] = (yshifts[cur_group] * (centers[cur_group + 1] - iframe) +
		                                yshifts[cur_group + 1] * (iframe - centers[cur_group])) /
		                               (centers[cur_group + 1] - centers[cur_group]);
#ifdef DEBUG_OWN
		std::cout << "iframe = " << iframe << " igroup " << cur_group << " x " << interpolated_xshifts[iframe] << " y " << interpolated_yshifts[iframe] << std::endl;
#endif
	}
}

void MotioncorrRunner::realSpaceInterpolation(Image <float> &Isum, std::vector<Image<float> > &Iframes, MotionModel *model, std::ostream &logfile) {
	int model_version = MOTION_MODEL_NULL;
	if (model != NULL) {
		model_version = model->getModelVersion();
	}

	const int n_frames = Iframes.size();
	if (model_version == MOTION_MODEL_NULL) {
		// Simple sum
		for (int iframe = 0; iframe < n_frames; iframe++) {
			logfile << "." << std::flush;

			#pragma omp parallel for num_threads(n_threads)
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Isum()) {
				DIRECT_MULTIDIM_ELEM(Isum(), n) += DIRECT_MULTIDIM_ELEM(Iframes[iframe](), n);
			}
		}

	} else if (model_version == MOTION_MODEL_THIRD_ORDER_POLYNOMIAL) { // Optimised code
		ThirdOrderPolynomialModel *polynomial_model = (ThirdOrderPolynomialModel*)model;
		realSpaceInterpolation_ThirdOrderPolynomial(Isum, Iframes, *polynomial_model, logfile);
	} else { // general code
		const int nx = XSIZE(Iframes[0]()), ny = YSIZE(Iframes[0]());
		for (int iframe = 0; iframe < n_frames; iframe++) {
			logfile << "." << std::flush;
			const RFLOAT z = iframe;

			#pragma omp parallel for num_threads(n_threads)
			for (int ix = 0; ix < nx; ix++) {
				const RFLOAT x = (RFLOAT)ix / nx - 0.5;
				for (int iy = 0; iy < ny; iy++) {
					bool valid = true;
					const RFLOAT y = (RFLOAT)iy / ny - 0.5;

					RFLOAT x_fitted, y_fitted;
					model->getShiftAt(z, x, y, x_fitted, y_fitted);
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
						DIRECT_A2D_ELEM(Isum(), iy, ix) += DIRECT_A2D_ELEM(Iframes[iframe](), y0, x0);
						if (std::isnan(DIRECT_A2D_ELEM(Isum(), iy, ix))) {
							std::cerr << "ix = " << ix << " xfit = " << x_fitted << " iy = " << iy << " ifit = " << y_fitted << std::endl;
						}
						continue;
					}

					const RFLOAT fx = x_fitted - x0;
					const RFLOAT fy = y_fitted - y0;

					const RFLOAT d00 = DIRECT_A2D_ELEM(Iframes[iframe](), y0, x0);
					const RFLOAT d01 = DIRECT_A2D_ELEM(Iframes[iframe](), y0, x1);
					const RFLOAT d10 = DIRECT_A2D_ELEM(Iframes[iframe](), y1, x0);
					const RFLOAT d11 = DIRECT_A2D_ELEM(Iframes[iframe](), y1, x1);

					const RFLOAT dx0 = LIN_INTERP(fx, d00, d01);
					const RFLOAT dx1 = LIN_INTERP(fx, d10, d11);
					const RFLOAT val = LIN_INTERP(fy, dx0, dx1);
#ifdef DEBUG_OWN
					if (std::isnan(val)) {
						std::cerr << "ix = " << ix << " xfit = " << x_fitted << " iy = " << iy << " ifit = " << y_fitted << " d00 " << d00 << " d01 " << d01 << " d10 " << d10 << " d11 " << d11 << " dx0 " << dx0 << " dx1 " << dx1 << std::endl;
					}
#endif
					DIRECT_A2D_ELEM(Isum(), iy, ix) += val;
				} // y
			} // x
		} // frame
	} // general model
}

void MotioncorrRunner::realSpaceInterpolation_ThirdOrderPolynomial(Image <float> &Isum, std::vector<Image<float> > &Iframes, ThirdOrderPolynomialModel &model, std::ostream &logfile) {
	const int n_frames = Iframes.size();
	const int nx = XSIZE(Iframes[0]()), ny = YSIZE(Iframes[0]());
	const Matrix1D<RFLOAT> coeffX = model.coeffX, coeffY = model.coeffY;

	for (int iframe = 0; iframe < n_frames; iframe++) {
		logfile << "." << std::flush;

		const RFLOAT z = iframe, z2 = iframe * iframe;
		const RFLOAT z3 = z * z2;
		// Common terms
		const RFLOAT x_C0 = coeffX(0)  * z + coeffX(1)  * z2 + coeffX(2)  * z3;
		const RFLOAT x_C1 = coeffX(3)  * z + coeffX(4)  * z2 + coeffX(5)  * z3;
		const RFLOAT x_C2 = coeffX(6)  * z + coeffX(7)  * z2 + coeffX(8)  * z3;
		const RFLOAT x_C3 = coeffX(9)  * z + coeffX(10) * z2 + coeffX(11) * z3;
		const RFLOAT x_C4 = coeffX(12) * z + coeffX(13) * z2 + coeffX(14) * z3;
		const RFLOAT x_C5 = coeffX(15) * z + coeffX(16) * z2 + coeffX(17) * z3;
		const RFLOAT y_C0 = coeffY(0)  * z + coeffY(1)  * z2 + coeffY(2)  * z3;
		const RFLOAT y_C1 = coeffY(3)  * z + coeffY(4)  * z2 + coeffY(5)  * z3;
		const RFLOAT y_C2 = coeffY(6)  * z + coeffY(7)  * z2 + coeffY(8)  * z3;
		const RFLOAT y_C3 = coeffY(9)  * z + coeffY(10) * z2 + coeffY(11) * z3;
		const RFLOAT y_C4 = coeffY(12) * z + coeffY(13) * z2 + coeffY(14) * z3;
		const RFLOAT y_C5 = coeffY(15) * z + coeffY(16) * z2 + coeffY(17) * z3;

		#pragma omp parallel for num_threads(n_threads)
		for (int ix = 0; ix < nx; ix++) {
			const RFLOAT x = (RFLOAT)ix / nx - 0.5;
			for (int iy = 0; iy < ny; iy++) {
				bool valid = true;
				const RFLOAT y = (RFLOAT)iy / ny - 0.5;

				RFLOAT x_fitted = x_C0 + (x_C1 + x_C2 * x) * x + (x_C3 + x_C4 * y + x_C5 * x) * y;
				RFLOAT y_fitted = y_C0 + (y_C1 + y_C2 * x) * x + (y_C3 + y_C4 * y + y_C5 * x) * y;

#ifdef VALIDATE_OPTIMISED_CODE
				RFLOAT x_fitted2, y_fitted2;
				model.getShiftAt(z, x, y, x_fitted2, y_fitted2);
				if (abs(x_fitted - x_fitted2) > 1E-4 || abs(y_fitted - y_fitted2) > 1E-4) {
					std::cout << "error at " << x << ", " << y << " : " << x_fitted << " " << x_fitted2 << " " << y_fitted << " " << y_fitted2 << std::endl;
				}
#endif
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
					DIRECT_A2D_ELEM(Isum(), iy, ix) += DIRECT_A2D_ELEM(Iframes[iframe](), y0, x0);
					if (std::isnan(DIRECT_A2D_ELEM(Isum(), iy, ix))) {
						std::cerr << "ix = " << ix << " xfit = " << x_fitted << " iy = " << iy << " ifit = " << y_fitted << std::endl;
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
#ifdef DEBUG_OWN
				if (std::isnan(val)) {
					std::cerr << "ix = " << ix << " xfit = " << x_fitted << " iy = " << iy << " ifit = " << y_fitted << " d00 " << d00 << " d01 " << d01 << " d10 " << d10 << " d11 " << d11 << " dx0 " << dx0 << " dx1 " << dx1 << std::endl;
				}
#endif
				DIRECT_A2D_ELEM(Isum(), iy, ix) += val;
			}
		}
	}
}

bool MotioncorrRunner::alignPatch(std::vector<MultidimArray<fComplex> > &Fframes, const int pnx, const int pny, const RFLOAT scaled_B, std::vector<RFLOAT> &xshifts, std::vector<RFLOAT> &yshifts, std::ostream &logfile) {
	std::vector<Image<float> > Iccs(n_threads);
	MultidimArray<fComplex> Fref;
	std::vector<MultidimArray<fComplex> > Fccs(n_threads);
	MultidimArray<float> weight;
	std::vector<RFLOAT> cur_xshifts, cur_yshifts;
	bool converged = false;

	// Parameters TODO: make an option
	const int max_iter = 5;
	int search_range = 50; // px 
	const RFLOAT tolerance = 0.5; // px
	const RFLOAT EPS = 1e-15;

	// Shifts within an iteration
	const int n_frames = xshifts.size();
	cur_xshifts.resize(n_frames);
	cur_yshifts.resize(n_frames);

	if (pny % 2 == 1 || pnx % 2 == 1) {
		REPORT_ERROR("Patch size must be even");
	}

	// Calculate the size of down-sampled CCF
	float ccf_requested_scale = ccf_downsample;
	if (ccf_downsample <= 0) {
		ccf_requested_scale = sqrt(-log(1E-8) / (2 * scaled_B)); // exp(-2 B max_dist^2) = 1E-8 
	}
	int ccf_nx = findGoodSize(int(pnx * ccf_requested_scale)), ccf_ny = findGoodSize(int(pny * ccf_requested_scale));
	if (ccf_nx > pnx) ccf_nx = pnx;
	if (ccf_ny > pny) ccf_ny = pny;
	if (ccf_nx % 2 == 1) ccf_nx++;
	if (ccf_ny % 2 == 1) ccf_ny++;
	const int ccf_nfx = ccf_nx / 2 + 1, ccf_nfy = ccf_ny;
	const int ccf_nfy_half = ccf_ny / 2;
	const RFLOAT ccf_scale_x = (RFLOAT)pnx / ccf_nx;
	const RFLOAT ccf_scale_y = (RFLOAT)pny / ccf_ny;
	search_range /= (ccf_scale_x > ccf_scale_y) ? ccf_scale_x : ccf_scale_y; // account for the increase of pixel size in CCF
	if (search_range * 2 + 1 > ccf_nx) search_range = ccf_nx / 2 - 1;
	if (search_range * 2 + 1 > ccf_ny) search_range = ccf_ny / 2 - 1;

	const int nfx = XSIZE(Fframes[0]), nfy = YSIZE(Fframes[0]);
	const int nfy_half = nfy / 2;

	Fref.reshape(ccf_nfy, ccf_nfx);
	for (int i = 0; i < n_threads; i++) {
		Iccs[i]().reshape(ccf_ny, ccf_nx);
		Fccs[i].reshape(Fref);
	}

#ifdef DEBUG
	std::cout << "Patch Size X = " << pnx << " Y  = " << pny << std::endl;
	std::cout << "Fframes X = " << nfx << " Y = " << nfy << std::endl;
	std::cout << "Fccf X = " << ccf_nfx << " Y = " << ccf_nfy << std::endl;
	std::cout << "CCF crop request = " << ccf_requested_scale << ", actual X = " << 1 / ccf_scale_x << " Y = " << 1 / ccf_scale_y << std::endl;
	std::cout << "CCF search range = " << search_range << std::endl;
	std::cout << "Trajectory size: " << xshifts.size() << std::endl;
#endif

	// Initialize B factor weight
	weight.reshape(Fref);
	RCTIC(TIMING_PREP_WEIGHT);
	#pragma omp parallel for num_threads(n_threads) 
	for (int y = 0; y < ccf_nfy; y++) {
		const int ly = (y > ccf_nfy_half) ? (y - ccf_nfy) : y;
		RFLOAT ly2 = ly * (RFLOAT)ly / (nfy * (RFLOAT)nfy);

		for (int x = 0; x < ccf_nfx; x++) {
			RFLOAT dist2 = ly2 + x * (RFLOAT)x / (nfx * (RFLOAT)nfx);
			DIRECT_A2D_ELEM(weight, y, x) = exp(- 2 * dist2 * scaled_B); // 2 for Fref and Fframe
		}
	}
	RCTOC(TIMING_PREP_WEIGHT);

	for (int iter = 1; iter	<= max_iter; iter++) {
		RCTIC(TIMING_MAKE_REF);
		Fref.initZeros();

		#pragma omp parallel for num_threads(n_threads)
		for (int y = 0; y < ccf_nfy; y++) {
			const int ly = (y > ccf_nfy_half) ? (y - ccf_nfy + nfy) : y;
			for (int x = 0; x < ccf_nfx; x++) {
				for (int iframe = 0; iframe < n_frames; iframe++) {
					DIRECT_A2D_ELEM(Fref, y, x) += DIRECT_A2D_ELEM(Fframes[iframe], ly, x);
				}
			}
		}
		RCTOC(TIMING_MAKE_REF);

		#pragma omp parallel for num_threads(n_threads)
		for (int iframe = 0; iframe < n_frames; iframe++) {
			const int tid = omp_get_thread_num();

			RCTIC(TIMING_CCF_CALC);
			for (int y = 0; y < ccf_nfy; y++) {
				const int ly = (y > ccf_nfy_half) ? (y - ccf_nfy + nfy) : y;
				for (int x = 0; x < ccf_nfx; x++) {
					DIRECT_A2D_ELEM(Fccs[tid], y, x) = (DIRECT_A2D_ELEM(Fref, y, x) - DIRECT_A2D_ELEM(Fframes[iframe], ly, x)) *
					                                    DIRECT_A2D_ELEM(Fframes[iframe], ly, x).conj() * DIRECT_A2D_ELEM(weight, y, x);
				}
			}
			RCTOC(TIMING_CCF_CALC);

			RCTIC(TIMING_CCF_IFFT);
			NewFFT::inverseFourierTransform(Fccs[tid], Iccs[tid]());
			RCTOC(TIMING_CCF_IFFT);

			RCTIC(TIMING_CCF_FIND_MAX);
			RFLOAT maxval = -1E30;
			int posx = 0, posy = 0;
			for (int y = -search_range; y <= search_range; y++) {
				const int iy = (y < 0) ? ccf_ny + y : y;

				for (int x = -search_range; x <= search_range; x++) {
					const int ix = (x < 0) ? ccf_nx + x : x;
					RFLOAT val = DIRECT_A2D_ELEM(Iccs[tid](), iy, ix);
//					std::cout << "(x, y) = " << x << ", " << y << ", (ix, iy) = " << ix << " , " << iy << " val = " << val << std::endl;
					if (val > maxval) {
						posx = x; posy = y;
						maxval = val;
					}
				}
			}

			int ipx_n = posx - 1, ipx = posx, ipx_p = posx + 1, ipy_n = posy - 1, ipy = posy, ipy_p = posy + 1;
			if (ipx_n < 0) ipx_n = ccf_nx + ipx_n;
			if (ipx < 0) ipx = ccf_nx + ipx;
			if (ipx_p < 0) ipx_p = ccf_nx + ipx_p;
			if (ipy_n < 0) ipy_n = ccf_ny + ipy_n;
			if (ipy < 0) ipy = ccf_ny + ipy;
			if (ipy_p < 0) ipy_p = ccf_ny + ipy_p;

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
			cur_xshifts[iframe] *= ccf_scale_x;
			cur_yshifts[iframe] *= ccf_scale_y;
#ifdef DEBUG_OWN
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
		#pragma omp parallel for num_threads(n_threads)
		for (int iframe = 1; iframe < n_frames; iframe++) {
			shiftNonSquareImageInFourierTransform(Fframes[iframe], -cur_xshifts[iframe] / pnx, -cur_yshifts[iframe] / pny);
		}
		RCTOC(TIMING_FOURIER_SHIFT);

		// Test convergence
		RFLOAT rmsd = std::sqrt((x_sumsq + y_sumsq) / n_frames);
		logfile << " Iteration " << iter << ": RMSD = " << rmsd << " px" << std::endl;

		if (rmsd < tolerance) {
			converged = true;
			break;
		}
	}

#ifdef DEBUG_OWN
	for (int iframe = 0; iframe < n_frames; iframe++) {
		std::cout << iframe << " " << xshifts[iframe] << " " << yshifts[iframe] << std::endl;
	}
#endif

	return converged;
}

int MotioncorrRunner::findGoodSize(int request) {
	// numbers that do not contain large prime numbers
	const int good_numbers[] = {192, 216, 256, 288, 324, 
	                            384, 432, 486, 512, 576, 648,
	                            768, 800, 864, 972, 1024,
	                            1296, 1536, 1728, 1944,
	                            2048, 2304, 2592, 3072, 3200,
	                            3456, 3888, 4096, 4608, 5000, 5184,
	                            6144, 6250, 6400, 6912, 7776, 8192,
	                            9216, 10240, 12288, 12500, -1};
	//return request; // DEBUG

	for (int i = 0; good_numbers[i] != -1; i++) {
		if (good_numbers[i] < request) continue;
		else return good_numbers[i];
	}
	return request; // return as it is
}
// dose is equivalent dose at 300 kV at the END of the frame.
// This implements the model by Timothy Grant & Nikolaus Grigorieff on eLife, 2015
// doi: 10.7554/eLife.06980
void MotioncorrRunner::doseWeighting(std::vector<MultidimArray<fComplex> > &Fframes, std::vector<RFLOAT> doses, RFLOAT apix) {
	const int nfx = XSIZE(Fframes[0]), nfy = YSIZE(Fframes[0]);
	const int nfy_half = nfy / 2;
	const RFLOAT nfy2 = (RFLOAT)nfy * nfy;
	const RFLOAT nfx2 = (RFLOAT)(nfx - 1) * (nfx - 1) * 4; // assuming nx is even
	const int n_frames= Fframes.size();
	const RFLOAT A = 0.245, B = -1.665, C = 2.81;

	#pragma omp parallel for num_threads(n_threads)
	for (int y = 0; y < nfy; y++) {
		int ly = y;
		if (y > nfy_half) ly = y - nfy;

		const RFLOAT ly2 = (RFLOAT)ly * ly / nfy2;
		for (int x = 0; x < nfx; x++) {
			const RFLOAT dinv2 = ly2 + (RFLOAT)x * x / nfx2;
			const RFLOAT dinv = std::sqrt(dinv2) / apix; // d = N * apix / dist, thus dinv = dist / N / angpix
			const RFLOAT Ne = (A * std::pow(dinv, B) + C) * 2; // Eq. 3. 2 comes from Eq. 5
			RFLOAT sum_weight_sq = 0;

			for (int iframe = 0; iframe < n_frames; iframe++) {
				const RFLOAT weight = std::exp(- doses[iframe] / Ne); // Eq. 5. 0.5 is factored out to Ne.
				if (std::isnan(weight)) {
					std::cerr << "dose = " <<  doses[iframe] << " Ne = " << Ne << " frm = " << iframe << " lx = " << x << " ly = " << ly << " reso = " << 1 / dinv << " weight = " << weight << std::endl;
				}
				sum_weight_sq += weight * weight;
				DIRECT_A2D_ELEM(Fframes[iframe], y, x) *= weight;
			}

			sum_weight_sq = std::sqrt(sum_weight_sq);
			if (std::isnan(sum_weight_sq)) {
				std::cerr << " Ne = " << Ne << " lx = " << x << " ly = " << ly << " reso = " << 1 / dinv << " sum_weight_sq NaN" << std::endl;
				REPORT_ERROR("Shouldn't happen.");
			}
			for (int iframe = 0; iframe < n_frames; iframe++) {
				DIRECT_A2D_ELEM(Fframes[iframe], y, x) /= sum_weight_sq; // Eq. 9
			}
		}
	}
}

// shiftx, shifty is relative to the (real space) image size
void MotioncorrRunner::shiftNonSquareImageInFourierTransform(MultidimArray<fComplex> &frame, RFLOAT shiftx, RFLOAT shifty) {
	const int nfx = XSIZE(frame), nfy = YSIZE(frame);
	const int nfy_half = nfy / 2;
	const RFLOAT twoPI = 2 * PI;

// Reduce calls to SINCOS from nfx * nfy to nfx + nfy
#define USE_TABLE
#ifdef USE_TABLE
	std::vector<RFLOAT> sinx(nfx), cosx(nfx), siny(nfy), cosy(nfy);
	for (int y = 0; y < nfy; y++) {
		int ly = y;
		if (y > nfy_half) ly = y - nfy;

		const RFLOAT phase_y = twoPI * ly * shifty;
		#ifdef RELION_SINGLE_PRECISION
		SINCOSF(phase_y, &siny[y], &cosy[y]);
		#else
		SINCOS(phase_y, &siny[y], &cosy[y]);
		#endif
	}
	for (int x = 0; x < nfx; x++) {
		const RFLOAT phase_x = twoPI * x * shiftx;
		#ifdef RELION_SINGLE_PRECISION
		SINCOSF(phase_x, &sinx[x], &cosx[x]);
		#else
		SINCOS(phase_x, &sinx[x], &cosx[x]);
		#endif
	}
#endif

	for (int y = 0; y < nfy; y++) {
#ifndef USE_TABLE
		int ly = y;
		if (y > nfy_half) ly = y - nfy;
#endif
		for (int x = 0; x < nfx; x++) {
			RFLOAT a, b, c, d, ac, bd, ab_cd;
#ifdef USE_TABLE
			b = sinx[x] * cosy[y] + cosx[x] * siny[y];
			a = cosx[x] * cosy[y] - sinx[x] * siny[y];
#else
			RFLOAT phase_shift = twoPI * (x * shiftx + ly * shifty);
			#ifdef RELION_SINGLE_PRECISION
			SINCOSF(phase_shift, &b, &a);
			#else
			SINCOS(phase_shift, &b, &a);
			#endif
#endif
			c = DIRECT_A2D_ELEM(frame, y, x).real;
			d = DIRECT_A2D_ELEM(frame, y, x).imag;
			ac = a * c;
			bd = b * d;
			ab_cd = (a + b) * (c + d); // (ab_cd-ac-bd = ad+bc : but needs 4 multiplications)
			DIRECT_A2D_ELEM(frame, y, x) = fComplex(ac - bd, ab_cd - ac - bd);
		}
	}
}

// Iwork is overwritten
void MotioncorrRunner::binNonSquareImage(Image<float> &Iwork, RFLOAT bin_factor) {
	const int nx = XSIZE(Iwork()), ny = YSIZE(Iwork());
	int new_nx = nx / bin_factor, new_ny = ny / bin_factor;
	if (new_nx % 2 != 0 || new_ny % 2 != 0) {
		// Do not adjust the size because it might lead to non-square pixels
		REPORT_ERROR("The dimensions of the image after binning must be even");
	}
#ifdef DEBUG_OWN
	std::cout << "Binning from X = " << nx << " Y = " << ny << " to X = " << new_nx << " Y = " << new_ny << std::endl;
#endif
	MultidimArray<fComplex> Fref(ny, nx / 2 + 1), Fbinned(new_ny, new_nx / 2 + 1);
	NewFFT::FourierTransform(Iwork(), Fref);

	cropInFourierSpace(Fref, Fbinned);

	Iwork().reshape(new_ny, new_nx);
	NewFFT::inverseFourierTransform(Fbinned, Iwork());
}

void MotioncorrRunner::cropInFourierSpace(MultidimArray<fComplex> &Fref, MultidimArray<fComplex> &Fbinned) {
	const int nfx = XSIZE(Fref), nfy = YSIZE(Fref);
	const int new_nfx = XSIZE(Fbinned), new_nfy = YSIZE(Fbinned);
	const int half_new_nfy = new_nfy / 2;

	if (new_nfx > nfx || new_nfy > nfy) REPORT_ERROR("Invalid size given to cropInFourierSpace");

	for (int y = 0; y < half_new_nfy; y++) {
		for (int x = 0; x < new_nfx; x++) {
			DIRECT_A2D_ELEM(Fbinned, y, x) =  DIRECT_A2D_ELEM(Fref, y, x);
		}
	}
	for (int y = half_new_nfy; y < new_nfy; y++) {
		for (int x = 0; x < new_nfx; x++) {
			DIRECT_A2D_ELEM(Fbinned, y, x) =  DIRECT_A2D_ELEM(Fref, nfy - new_nfy + y, x);
		}
	}
}

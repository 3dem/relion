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
#include <cmath>

#ifdef _CUDA_ENABLED
#include "src/acc/cuda/cuda_mem_utils.h"
#endif

void CtffindRunner::read(int argc, char **argv, int rank)
{
	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	int ctf_section = parser.addSection("CTF estimation");
	fn_in = parser.getOption("--i", "STAR file with all input micrographs, or a unix wildcard to all micrograph files, e.g. \"mics/*.mrc\"");
	do_use_without_doseweighting = parser.checkOption("--use_noDW", "Estimate CTFs from rlnMicrographNameNoDW instead of rlnMicrographName (only after MotionCor2)");
	fn_out = parser.getOption("--o", "Directory, where all output files will be stored", "CtfEstimate/");
	do_only_join_results = parser.checkOption("--only_make_star", "Don't estimate any CTFs, only join all logfile results in a STAR file");
	continue_old = parser.checkOption("--only_do_unfinished", "Only estimate CTFs for those micrographs for which there is not yet a logfile with Final values.");
	do_at_most = textToInteger(parser.getOption("--do_at_most", "Only process up to this number of (unprocessed) micrographs.", "-1"));
	// Use a smaller squared part of the micrograph to estimate CTF (e.g. to avoid film labels...)
	ctf_win =  textToInteger(parser.getOption("--ctfWin", "Size (in pixels) of a centered, squared window to use for CTF-estimation", "-1"));

	int mic_section = parser.addSection("Microscopy parameters");
	// First parameter line in CTFFIND
	Cs = textToFloat(parser.getOption("--CS", "Spherical Aberration (mm) ","-1"));
	Voltage = textToFloat(parser.getOption("--HT", "Voltage (kV)","-1"));
	AmplitudeConstrast = textToFloat(parser.getOption("--AmpCnst", "Amplitude constrast", "-1"));
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in the input micrographs (A)", "-1"));

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

	int ctffind4_section = parser.addSection("CTFFIND4 parameters");
	is_ctffind4 = parser.checkOption("--is_ctffind4", "The provided CTFFIND executable is CTFFIND4 (version 4.1+)");
	use_given_ps = parser.checkOption("--use_given_ps", "Use pre-calculated power spectra?");
	do_movie_thon_rings = parser.checkOption("--do_movie_thon_rings", "Calculate Thon rings from movie frames?");
	avg_movie_frames = textToInteger(parser.getOption("--avg_movie_frames", "Average over how many movie frames (try to get 4 e-/A2)", "1"));
	movie_rootname = parser.getOption("--movie_rootname", "Rootname plus extension for movies", "_movie.mrcs");
	do_phaseshift = parser.checkOption("--do_phaseshift", "Estimate the phase shift in the images (e.g. from a phase-plate)");
	phase_min  = textToFloat(parser.getOption("--phase_min", "Minimum phase shift (in degrees)", "0."));
	phase_max  = textToFloat(parser.getOption("--phase_max", "Maximum phase shift (in degrees)", "180."));
	phase_step = textToFloat(parser.getOption("--phase_step", "Step in phase shift (in degrees)", "10."));
	nr_threads = textToInteger(parser.getOption("--j", "Number of threads (for CTFIND4 only)", "1"));
	do_fast_search = parser.checkOption("--fast_search", "Disable \"Slower, more exhaustive search\" in CTFFIND4.1 (faster but less accurate)");

	int gctf_section = parser.addSection("Gctf parameters");
	do_use_gctf = parser.checkOption("--use_gctf", "Use Gctf instead of CTFFIND to estimate the CTF parameters");
	fn_gctf_exe = parser.getOption("--gctf_exe","Location of Gctf executable (or through RELION_GCTF_EXECUTABLE environment variable)","");
	angpix = textToFloat(parser.getOption("--angpix", "Magnified pixel size in Angstroms", "1."));
	do_ignore_ctffind_params = parser.checkOption("--ignore_ctffind_params", "Use Gctf default parameters instead of CTFFIND parameters");
	do_EPA = parser.checkOption("--EPA", "Use equi-phase averaging to calculate Thon rinds in Gctf");
	do_validation = parser.checkOption("--do_validation", "Use validation inside Gctf to analyse quality of the fit?");
	additional_gctf_options = parser.getOption("--extra_gctf_options", "Additional options for Gctf", "");
	gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread, e.g 0:1:2:3","");

	// Initialise verb for non-parallel execution
	verb = 1;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void CtffindRunner::usage()
{
	parser.writeUsage(std::cout);
}

void CtffindRunner::initialise(bool is_leader)
{
	// Get the CTFFIND executable
	if (fn_ctffind_exe == "")
	{
		char *penv;
		penv = getenv("RELION_CTFFIND_EXECUTABLE");
		if (penv != NULL)
			fn_ctffind_exe = (std::string)penv;
	}
	// Get the GCTF executable
	if (do_use_gctf && fn_gctf_exe == "")
	{
		char *penv;
		penv = getenv("RELION_GCTF_EXECUTABLE");
		if (penv != NULL)
			fn_gctf_exe = (std::string)penv;
	}

	fn_shell = "/bin/sh";
	char *shell_name;
	shell_name = getenv("RELION_SHELL");
	if (shell_name != NULL)
		fn_shell = (std::string)shell_name;

	if (do_use_gctf && ctf_win>0)
		REPORT_ERROR("ERROR: Running Gctf together with --ctfWin is not implemented, please use CTFFIND instead.");

	if (!do_phaseshift && (additional_gctf_options.find("--phase_shift_L") != std::string::npos ||
	                       additional_gctf_options.find("--phase_shift_H") != std::string::npos ||
	                       additional_gctf_options.find("--phase_shift_S") != std::string::npos))
		REPORT_ERROR("ERROR: Please don't specify --phase_shift_L, H, S in 'Other Gctf options' (--extra_gctf_options). Use 'Estimate phase shifts' (--do_phaseshift) and 'Phase shift - Min, Max, Step' (--phase_min, --phase_max, --phase_step) instead.");

	if (do_use_gctf && use_given_ps)
		REPORT_ERROR("ERROR: --use_given_ps is available only with CTFFIND 4.1");

	if (use_given_ps && do_movie_thon_rings)
		REPORT_ERROR("ERROR: You cannot enable --use_given_ps and --do_movie_thon_rings simultaneously");

	if (use_given_ps)
		do_use_without_doseweighting = false;

	// Make sure fn_out ends with a slash
	if (fn_out[fn_out.length()-1] != '/')
		fn_out += "/";

	// Set up which micrographs to estimate CTFs from
	if (fn_in.isStarFile())
	{
		MetaDataTable MDin;
		ObservationModel::loadSafely(fn_in, obsModel, MDin, "micrographs", verb);

		if (MDin.numberOfObjects() > 0 && !MDin.containsLabel(EMDL_MICROGRAPH_NAME))
			REPORT_ERROR("ERROR: There is no rlnMicrographName label in the input micrograph STAR file.");

		if (do_use_without_doseweighting && MDin.numberOfObjects() > 0 && !MDin.containsLabel(EMDL_MICROGRAPH_NAME_WODOSE))
			REPORT_ERROR("ERROR: You are using --use_noDW, but there is no rlnMicrographNameNoDW label in the input micrograph STAR file.");

		if (use_given_ps && MDin.numberOfObjects() > 0 && !MDin.containsLabel(EMDL_CTF_POWER_SPECTRUM))
			REPORT_ERROR("ERROR: You are using --use_given_ps, but there is no rlnCtfPowerSpectrum label in the input micrograph STAR file.");

		fn_micrographs_all.clear();
		optics_group_micrographs_all.clear();
		fn_micrographs_ctf_all.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			FileName fn_mic;
			MDin.getValue(EMDL_MICROGRAPH_NAME, fn_mic);

			fn_micrographs_all.push_back(fn_mic); // Dose weighted image

			if (do_use_without_doseweighting)
				MDin.getValue(EMDL_MICROGRAPH_NAME_WODOSE, fn_mic);
			else if (use_given_ps)
				MDin.getValue(EMDL_CTF_POWER_SPECTRUM, fn_mic);
			fn_micrographs_ctf_all.push_back(fn_mic); // Image for CTF estsimation

			int optics_group;
			MDin.getValue(EMDL_IMAGE_OPTICS_GROUP, optics_group);
			optics_group_micrographs_all.push_back(optics_group);
		}
	}
	else
	{
		fn_in.globFiles(fn_micrographs_all);
		optics_group_micrographs_all.resize(fn_in.size(), 1);
		obsModel.opticsMdt.clear();
		obsModel.opticsMdt.addObject();
	}

	// Make sure obsModel.opticsMdt has all the necessary information
	// If voltage or pixel size were not in the input STAR file, set them from the command line options
	if (!obsModel.opticsMdt.containsLabel(EMDL_CTF_CS))
	{
		if (Cs < 0.)
		{
			REPORT_ERROR("ERROR: the input STAR file does not contain the spherical aberration, and it is not given through --CS.");
		}
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
		{
			obsModel.opticsMdt.setValue(EMDL_CTF_CS, Cs);
		}
	}
	if (!obsModel.opticsMdt.containsLabel(EMDL_CTF_VOLTAGE))
	{
		if (Voltage < 0.)
		{
			REPORT_ERROR("ERROR: the input STAR file does not contain the acceleration voltage, and it is not given through --HT.");
		}
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
		{
			obsModel.opticsMdt.setValue(EMDL_CTF_VOLTAGE, Voltage);
		}
	}
	if (!obsModel.opticsMdt.containsLabel(EMDL_CTF_Q0))
	{
		if (AmplitudeConstrast < 0.)
		{
			REPORT_ERROR("ERROR: the input STAR file does not contain the amplitude contrast, and it is not given through --AmpCnst.");
		}
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
		{
			obsModel.opticsMdt.setValue(EMDL_CTF_Q0, AmplitudeConstrast);
		}
	}
	if (!obsModel.opticsMdt.containsLabel(EMDL_MICROGRAPH_PIXEL_SIZE))
	{
		if (angpix < 0.)
		{
			REPORT_ERROR("ERROR: the input STAR file does not contain the micrograph pixel size, and it is not given through --angpix.");
		}
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
		{
			obsModel.opticsMdt.setValue(EMDL_MICROGRAPH_PIXEL_SIZE, angpix);
		}
	}

	// First backup the given list of all micrographs
	std::vector<int> optics_group_given_all = optics_group_micrographs_all;
	std::vector<FileName> fn_mic_given_all = fn_micrographs_all;
	std::vector<FileName> fn_mic_ctf_given_all = fn_micrographs_ctf_all;
	// These lists contain those for the output STAR & PDF files
	optics_group_micrographs_all.clear();
	fn_micrographs_all.clear();
	fn_micrographs_ctf_all.clear();
	// These are micrographs to be processed
	optics_group_micrographs.clear();
	fn_micrographs.clear();
	fn_micrographs_ctf.clear();

	bool warned = false;
	for (long int imic = 0; imic < fn_mic_given_all.size(); imic++)
	{
		bool ignore_this = false;
		bool process_this = true;

		if (continue_old)
		{
			FileName fn_microot = fn_mic_ctf_given_all[imic].withoutExtension();
			RFLOAT defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep, maxres=-1., valscore = -1., phaseshift = 0.;
			if (getCtffindResults(fn_microot, defU, defV, defAng, CC,
			     HT, CS, AmpCnst, XMAG, DStep, maxres, valscore, phaseshift, false)) // false: dont warn if not found Final values
			{
				process_this = false; // already done
			}
		}

		if (do_at_most >= 0 && fn_micrographs.size() >= do_at_most)
		{
			if (process_this) {
				ignore_this = true;
				process_this = false;
				if (!warned)
				{
					warned = true;
					std::cout << "NOTE: processing of some micrographs will be skipped as requested by --do_at_most" << std::endl;
				}
			}
			// If this micrograph has already been processed, the result should be included in the output.
			// So ignore_this remains false.
		}

		if (process_this)
		{
			optics_group_micrographs.push_back(optics_group_given_all[imic]);
			fn_micrographs.push_back(fn_mic_given_all[imic]);
			fn_micrographs_ctf.push_back(fn_mic_ctf_given_all[imic]);
		}

		if (!ignore_this)
		{
			optics_group_micrographs_all.push_back(optics_group_given_all[imic]);
			fn_micrographs_all.push_back(fn_mic_given_all[imic]);
			fn_micrographs_ctf_all.push_back(fn_mic_ctf_given_all[imic]);
		}
	}

	if (is_leader)
	{
		std::cout << fn_mic_given_all.size() << " micrographs were given but we process only ";
		std::cout  << do_at_most << " micrographs as specified in --do_at_most." << std::endl;
	}

	// Make symbolic links of the input micrographs in the output directory because ctffind and gctf write output files alongside the input micropgraph
	if (is_leader)
	{
		char temp [180];
		char *cwd = getcwd(temp, 180);
		currdir = std::string(temp);
		// Make sure fn_out ends with a slash
		if (currdir[currdir.length()-1] != '/')
			currdir += "/";
		FileName prevdir="";
		for (size_t i = 0; i < fn_micrographs.size(); i++)
		{
			FileName myname = fn_micrographs_ctf[i];
			if (do_movie_thon_rings)
				myname = myname.withoutExtension() + movie_rootname;
			// Remove the UNIQDATE part of the filename if present
			FileName output = getOutputFileWithNewUniqueDate(myname, fn_out);
			// Create output directory if neccesary
			FileName newdir = output.beforeLastOf("/");
			if (newdir != prevdir)
			{
				mktree(newdir);
			}
			symlink(currdir + myname, output);
		}
	}

	if (do_use_gctf && fn_micrographs.size()>0)
	{
		untangleDeviceIDs(gpu_ids, allThreadIDs);
		if (allThreadIDs[0].size()==0 || (!std::isdigit(*gpu_ids.begin())) )
		{
#ifdef _CUDA_ENABLED
			if (verb>0)
				std::cout << "gpu-ids were not specified, so threads will automatically be mapped to devices (incrementally)."<< std::endl;
			HANDLE_ERROR(cudaGetDeviceCount(&devCount));
#else
			if (verb>0)
				REPORT_ERROR("gpu-ids were not specified, but we could not figure out which GPU to use because RELION was not compiled with CUDA support.");
#endif
		}

		// Find the dimensions of the first micrograph, to later on ensure all micrographs are the same size
		Image<double> Itmp;
		Itmp.read(fn_micrographs[0], false); // false means only read header!
		xdim = XSIZE(Itmp());
		ydim = YSIZE(Itmp());
	}

	if (is_ctffind4 && ctf_win > 0 && do_movie_thon_rings)
		REPORT_ERROR("CtffindRunner::initialise ERROR: You cannot use a --ctfWin operation on movies.");

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
			{
				if (is_ctffind4)
					std::cout << " Estimating CTF parameters using Alexis Rohou's and Niko Grigorieff's CTFFIND4.1 ..." << std::endl;
				else
					std::cout << " Estimating CTF parameters using Niko Grigorieff's CTFFIND ..." << std::endl;
			}
			init_progress_bar(fn_micrographs.size());
			barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
		}

		std::vector<std::string> allmicnames;
		for (long int imic = 0; imic < fn_micrographs.size(); imic++)
		{

			// Abort through the pipeline_control system
			if (pipeline_control_check_abort_job())
				exit(RELION_EXIT_ABORTED);

			// Get angpix and voltage from the optics groups:
			obsModel.opticsMdt.getValue(EMDL_CTF_CS, Cs, optics_group_micrographs[imic]-1);
			obsModel.opticsMdt.getValue(EMDL_CTF_VOLTAGE, Voltage, optics_group_micrographs[imic]-1);
			obsModel.opticsMdt.getValue(EMDL_CTF_Q0, AmplitudeConstrast, optics_group_micrographs[imic]-1);
			obsModel.opticsMdt.getValue(EMDL_MICROGRAPH_PIXEL_SIZE, angpix, optics_group_micrographs[imic]-1);

			if (do_use_gctf)
			{
				executeGctf(imic, allmicnames, imic+1==fn_micrographs.size());
			}
			else if (is_ctffind4)
			{
				executeCtffind4(imic);
			}
			else
			{
				executeCtffind3(imic);
			}

			if (verb > 0 && imic % barstep == 0)
				progress_bar(imic);
		}

		if (verb > 0)
			progress_bar(fn_micrographs.size());
	}

	joinCtffindResults();
}

void CtffindRunner::joinCtffindResults()
{
	long int barstep = XMIPP_MAX(1, fn_micrographs_all.size() / 60);
	if (verb > 0)
	{
		std::cout << " Generating joint STAR file ... " << std::endl;
		init_progress_bar(fn_micrographs_all.size());
	}

	MetaDataTable MDctf;
	for (long int imic = 0; imic < fn_micrographs_all.size(); imic++)
	{
		FileName fn_microot = fn_micrographs_ctf_all[imic].withoutExtension();
		RFLOAT defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep;
		RFLOAT maxres = -999., valscore = -999., phaseshift = -999.;
		bool has_this_ctf = getCtffindResults(fn_microot, defU, defV, defAng, CC,
		                                      HT, CS, AmpCnst, XMAG, DStep, maxres, valscore, phaseshift);

		if (!has_this_ctf)
		{
			std::cerr << " WARNING: skipping, since cannot get CTF values for " << fn_micrographs_all[imic] <<std::endl;
		}
		else
		{
			FileName fn_root = getOutputFileWithNewUniqueDate(fn_microot, fn_out);
			FileName fn_ctf = fn_root + ".ctf:mrc";
			MDctf.addObject();

			if (do_use_without_doseweighting)
				MDctf.setValue(EMDL_MICROGRAPH_NAME_WODOSE, fn_micrographs_ctf_all[imic]);
			MDctf.setValue(EMDL_MICROGRAPH_NAME, fn_micrographs_all[imic]);
			MDctf.setValue(EMDL_IMAGE_OPTICS_GROUP, optics_group_micrographs_all[imic]);
			MDctf.setValue(EMDL_CTF_IMAGE, fn_ctf);
			MDctf.setValue(EMDL_CTF_DEFOCUSU, defU);
			MDctf.setValue(EMDL_CTF_DEFOCUSV, defV);
			MDctf.setValue(EMDL_CTF_ASTIGMATISM, fabs(defU-defV));
			MDctf.setValue(EMDL_CTF_DEFOCUS_ANGLE, defAng);
			if (!std::isfinite(CC)) CC = 0.0; // GCTF might return NaN
			MDctf.setValue(EMDL_CTF_FOM, CC);
			if (fabs(maxres + 999.) > 0.)
			{
				// Put an upper limit on maxres, as gCtf may put 999. now max is 25.
				MDctf.setValue(EMDL_CTF_MAXRES, XMIPP_MIN(25., maxres));
			}
			if (fabs(phaseshift + 999.) > 0.)
				MDctf.setValue(EMDL_CTF_PHASESHIFT, phaseshift);
			if (fabs(valscore + 999.) > 0.)
				MDctf.setValue(EMDL_CTF_VALIDATIONSCORE, valscore);
		}

		if (verb > 0 && imic % 60 == 0) progress_bar(imic);
	}

	obsModel.save(MDctf, fn_out+"micrographs_ctf.star", "micrographs");

	if (verb > 0)
	{
		progress_bar(fn_micrographs_all.size());
		std::cout << " Done! Written out: " << fn_out <<  "micrographs_ctf.star" << std::endl;
	}

	if (verb > 0)
	{
		std::cout << " Now generating logfile.pdf ... " << std::endl;
	}

	std::vector<EMDLabel> plot_labels;
	plot_labels.push_back(EMDL_CTF_DEFOCUSU);
	plot_labels.push_back(EMDL_CTF_DEFOCUS_ANGLE);
	plot_labels.push_back(EMDL_CTF_ASTIGMATISM);
	plot_labels.push_back(EMDL_CTF_MAXRES);
	plot_labels.push_back(EMDL_CTF_PHASESHIFT);
	plot_labels.push_back(EMDL_CTF_FOM);
	plot_labels.push_back(EMDL_CTF_VALIDATIONSCORE);
	FileName fn_eps, fn_eps_root = fn_out+"micrographs_ctf";
	std::vector<FileName> all_fn_eps;
	for (int i = 0; i < plot_labels.size(); i++)
	{
		EMDLabel label = plot_labels[i];
		if (MDctf.containsLabel(label))
		{
			// Values for all micrographs
			CPlot2D *plot2Db=new CPlot2D(EMDL::label2Str(label) + " for all micrographs");
			MDctf.addToCPlot2D(plot2Db, EMDL_UNDEFINED, label, 1.);
			plot2Db->SetDrawLegend(false);
			fn_eps = fn_eps_root + "_all_" + EMDL::label2Str(label) + ".eps";
			plot2Db->OutputPostScriptPlot(fn_eps);
			all_fn_eps.push_back(fn_eps);
			delete plot2Db;
			if (MDctf.numberOfObjects() > 3)
			{
				// Histogram
				std::vector<RFLOAT> histX, histY;
				CPlot2D *plot2D=new CPlot2D("");
				MDctf.columnHistogram(label,histX,histY,0, plot2D);
				fn_eps = fn_eps_root + "_hist_" + EMDL::label2Str(label) + ".eps";
				plot2D->OutputPostScriptPlot(fn_eps);
				all_fn_eps.push_back(fn_eps);
				delete plot2D;
			}
		}
	}
	joinMultipleEPSIntoSinglePDF(fn_out + "logfile.pdf", all_fn_eps);

	if (verb > 0 )
	{
		std::cout << " Done! Written out: " << fn_out << "logfile.pdf" << std::endl;
	}

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

void CtffindRunner::executeGctf(long int imic, std::vector<std::string> &allmicnames, bool is_last, int rank)
{
	// Always add the new micrograph to the TODO list
	Image<double> Itmp;
	FileName outputfile = getOutputFileWithNewUniqueDate(fn_micrographs_ctf[imic], fn_out);
	Itmp.read(outputfile, false); // false means only read header!
	if (XSIZE(Itmp()) != xdim || YSIZE(Itmp()) != ydim)
		REPORT_ERROR("CtffindRunner::executeGctf ERROR: Micrographs do not all have the same size! " + fn_micrographs_ctf[imic] + " is different from the first micrograph!");
	if (ZSIZE(Itmp()) > 1 || NSIZE(Itmp()) > 1)
		REPORT_ERROR("CtffindRunner::executeGctf ERROR: No movies or volumes allowed for " + fn_micrographs_ctf[imic]);

	allmicnames.push_back(outputfile);

	// Execute Gctf every 20 images, and always for the last one
	if ( ((imic+1)%20) == 0 || is_last)
	{
		std::string command = fn_gctf_exe;
		//command +=  " --ctfstar " + fn_out + "tt_micrographs_ctf.star";
		command +=  " --apix " + floatToString(angpix);
		command +=  " --cs " + floatToString(Cs);
		command +=  " --kV " + floatToString(Voltage);
		command +=  " --ac " + floatToString(AmplitudeConstrast);
		command += " --astm " + floatToString(amount_astigmatism);
		command +=  " --logsuffix _gctf.log";

		if (!do_ignore_ctffind_params)
		{
			command += " --boxsize " + floatToString(box_size);
			command += " --resL " + floatToString(resol_min);
			command += " --resH " + floatToString(resol_max);
			command += " --defL " + floatToString(min_defocus);
			command += " --defH " + floatToString(max_defocus);
			command += " --defS " + floatToString(step_defocus);
		}

		if (do_phaseshift)
		{
			command += " --phase_shift_L " + floatToString(phase_min);
			command += " --phase_shift_H " + floatToString(phase_max);
			command += " --phase_shift_S " + floatToString(phase_step);
		}

		if (do_EPA)
			command += " --do_EPA ";

		if (do_validation)
			command += " --do_validation ";

		for (size_t i = 0; i<allmicnames.size(); i++)
			command += " " + allmicnames[i];

		if (allThreadIDs[0].size()==0 || (!std::isdigit(*gpu_ids.begin())) )
		{
			// Automated mapping
			command += " -gid " + integerToString(rank % devCount);
		}
		else
		{
			// User-specified mapping
			command += " -gid " + allThreadIDs[rank][0];
		}

		// extra Gctf options
		command += " " + additional_gctf_options;

		// Redirect all gctf output
		command += " >> " + fn_out + "gctf" + integerToString(rank)+".out  2>> " + fn_out + "gctf" + integerToString(rank)+".err";

		//std::cerr << " command= " << command << std::endl;
		int res = system(command.c_str());

		// Re-set the allmicnames vector
		allmicnames.clear();
	}
}

void CtffindRunner::executeCtffind3(long int imic)
{
	FileName fn_mic = getOutputFileWithNewUniqueDate(fn_micrographs_ctf[imic], fn_out);
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

	std::string ctffind4_options = (is_ctffind4) ? " --omp-num-threads " + integerToString(nr_threads) + " --old-school-input-ctffind4 " : "";

	// Write script to run ctffind
	fh << "#!/usr/bin/env " << fn_shell << std::endl;
	fh << fn_ctffind_exe << ctffind4_options << " > " << fn_log << " << EOF"<<std::endl;
	// line 1: input image
	if (do_movie_thon_rings)
		fh << fn_mic_win.withoutExtension() + movie_rootname << std::endl;
	else
		fh << fn_mic_win << std::endl;
	// line 2: diagnostic .ctf image
	fh << fn_ctf << std::endl;
	// line 3: CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]
	fh << Cs << ", " << Voltage << ", " << AmplitudeConstrast << ", 10000, " << angpix<< std::endl;
	// line 4: Box, ResMin[A], ResMax[A], dFMin[A], dFMax[A], FStep[A], dAst[A]
	fh << box_size << ", " << resol_min << ", " << resol_max << ", " << min_defocus << ", " << max_defocus << ", " << step_defocus << ", " << amount_astigmatism << std::endl;
	if (is_ctffind4)
	{
		// line 4: Movie Thon rings: $input_is_stack_of_frames,$number_of_frames_to_average
		if (do_movie_thon_rings)
			fh << " 1  " <<  integerToString(avg_movie_frames) << std::endl;
		else
			fh << " 0  1" << std::endl;
		// line 5: Phase-shifts: $find_phase_shift,$min_ps,$max_ps,$step_ps (in rads)
		if (do_phaseshift)
			fh << " 1, " << DEG2RAD(phase_min) << ", " << DEG2RAD(phase_max) << ", " << DEG2RAD(phase_step) << std::endl;
		else
			fh << " 0, 0, 3.15, 0.2" << std::endl;
	}
	fh <<"EOF"<<std::endl;
	fh.close();

	// Execute ctffind
	std::string command = fn_shell + " "+ fn_script;
	if (system(command.c_str()))
		std::cerr << "WARNING: there was an error in executing: " << command << std::endl;

	// Remove windowed file again
	if (ctf_win > 0)
	{
		if( remove( fn_mic_win.c_str() ) != 0 )
			std::cerr << "WARNING: there was an error deleting windowed micrograph file " << fn_mic_win << std::endl;
	}
}

void CtffindRunner::executeCtffind4(long int imic)
{
	FileName fn_mic = getOutputFileWithNewUniqueDate(fn_micrographs_ctf[imic], fn_out);
	FileName fn_root = fn_mic.withoutExtension();
	FileName fn_script = fn_root + "_ctffind4.com";
	FileName fn_log = fn_root + "_ctffind4.log";
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

	int ctf_boxsize = box_size;
	RFLOAT ctf_angpix = angpix;
	if (use_given_ps)
	{
		Image<RFLOAT> Ihead;
		Ihead.read(fn_mic_win, false);
		ctf_boxsize = XSIZE(Ihead());
		ctf_angpix = Ihead.samplingRateX();
	}
	//std::string ctffind4_options = " --omp-num-threads " + integerToString(nr_threads);
	std::string ctffind4_options = "";
	if (use_given_ps)
		ctffind4_options += " --amplitude-spectrum-input";

	// Write script to run ctffind
	fh << "#!/usr/bin/env " << fn_shell << std::endl;
	fh << fn_ctffind_exe << ctffind4_options << " > " << fn_log << " << EOF"<<std::endl;
	// line 1: input image
	if (do_movie_thon_rings)
	{
		fh << fn_mic_win.withoutExtension() + movie_rootname << std::endl;
		fh << "yes" << std::endl;
		fh << avg_movie_frames << std::endl;
	}
	else
	{
		fh << fn_mic_win << std::endl;
	}

	// line 2: diagnostic .ctf image
	fh << fn_ctf << std::endl;
	fh << ctf_angpix << std::endl;
	fh << Voltage << std::endl;
	fh << Cs << std::endl;
	fh << AmplitudeConstrast << std::endl;
	fh << ctf_boxsize << std::endl;
	fh << resol_min << std::endl;
	fh << resol_max << std::endl;
	fh << min_defocus << std::endl;
	fh << max_defocus << std::endl;
	fh << step_defocus << std::endl;
	// Do you know what astigmatism is present?
	fh << "no" << std::endl;
	// Slower, more exhaustive search?
	// The default was "no" in CTFFIND 4.1.5, but turned out to be less accurate.
	// The default was changed to "yes" in CTFFIND 4.1.8.
	// Ref: http://grigoriefflab.janelia.org/ctffind4
	// So, we say "yes" regardless of the version unless "--fast_search" is specified.
	if (!do_fast_search)
		fh << "yes" << std::endl;
	else
		fh << "no" << std::endl;
	// Use a restraint on astigmatism?
	fh << "yes" << std::endl;
	// Expected (tolerated) astigmatism
	fh << amount_astigmatism << std::endl;
	if (do_phaseshift)
	{
		fh << "yes" << std::endl;
		fh << DEG2RAD(phase_min) << std::endl;
		fh << DEG2RAD(phase_max) << std::endl;
		fh << DEG2RAD(phase_step) << std::endl;
	}
	else
		fh << "no" << std::endl;
	// Set expert options?
	fh << "no" << std::endl;

	fh <<"EOF"<<std::endl;
	fh << "exit 0" << std::endl;
	fh.close();

	// Execute ctffind
	FileName command = fn_shell + " "+ fn_script;
	if (system(command.c_str()))
		std::cerr << "WARNING: there was an error in executing: " << command << std::endl;

	// Remove windowed file again
	if (ctf_win > 0)
	{
		if( remove( fn_mic_win.c_str() ) != 0 )
			std::cerr << "WARNING: there was an error deleting windowed micrograph file " << fn_mic_win << std::endl;
	}
}

bool CtffindRunner::getCtffindResults(FileName fn_microot, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
		RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep,
		RFLOAT &maxres, RFLOAT &valscore, RFLOAT &phaseshift, bool do_warn)
{
	if (is_ctffind4)
	{
		return getCtffind4Results(fn_microot, defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep,
		                          maxres, phaseshift, do_warn);
	}
	else
	{
		return getCtffind3Results(fn_microot, defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep,
		                          maxres, phaseshift, valscore, do_warn);
	}
}

bool CtffindRunner::getCtffind3Results(FileName fn_microot, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
		RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep,
		RFLOAT &maxres, RFLOAT &phaseshift, RFLOAT &valscore, bool do_warn)
{
	FileName fn_root = getOutputFileWithNewUniqueDate(fn_microot, fn_out);
	FileName fn_log = fn_root + "_ctffind3.log";
	if (do_use_gctf)
		fn_log = fn_root + "_gctf.log";

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

		if (line.find("CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]") != std::string::npos ||
                    line.find("CS[mm], HT[kV], ac, XMAG, DStep[um]") != std::string::npos) // GCTF 1.18 B1 changed the line...
		{
			getline(in, line, '\n');
			tokenize(line, words);
			if (words.size() == 5)
			{
    				Cs_is_found = true;
				CS = textToFloat(words[0]);
				HT = textToFloat(words[1]);
				AmpCnst = textToFloat(words[2]);
				XMAG = textToFloat(words[3]);
				DStep = textToFloat(words[4]);
			}
		}

		int nr_exp_cols = (do_phaseshift) ? 7 : 6;
		if (line.find("Final Values") != std::string::npos)
		{
			tokenize(line, words);
			if (words.size() == nr_exp_cols)
			{
				Final_is_found = true;
				defU = textToFloat(words[0]);
				defV = textToFloat(words[1]);
				defAng = textToFloat(words[2]);
				if (do_use_gctf && do_phaseshift)
				{
					phaseshift = textToFloat(words[3]);
					CC = textToFloat(words[4]);
				}
				else
					CC = textToFloat(words[3]);
			}
		}

		if (do_use_gctf)
		{
			if (line.find("Resolution limit estimated by EPA:") != std::string::npos)
			{
				tokenize(line, words);
				maxres = textToFloat(words[words.size()-1]);
			}

			if (line.find("OVERALL_VALIDATION_SCORE:") != std::string::npos)
			{
				tokenize(line, words);
				valscore = textToFloat(words[words.size()-1]);
			}
		}
	}

	if (!Cs_is_found)
	{
		if (do_warn)
			std::cerr << "WARNING: cannot find line with Cs[mm], HT[kV], etc values in " << fn_log << std::endl;
		return false;
	}
	if (!Final_is_found)
	{
		if (do_warn)
			std::cerr << "WARNING: cannot find line with Final values in " << fn_log << std::endl;
		return false;
	}

	in.close();

	return Final_is_found;
}


bool CtffindRunner::getCtffind4Results(FileName fn_microot, RFLOAT &defU, RFLOAT &defV, RFLOAT &defAng, RFLOAT &CC,
		RFLOAT &HT, RFLOAT &CS, RFLOAT &AmpCnst, RFLOAT &XMAG, RFLOAT &DStep,
		RFLOAT &maxres, RFLOAT &phaseshift, bool do_warn)
{
	FileName fn_root = getOutputFileWithNewUniqueDate(fn_microot, fn_out);
	FileName fn_log = fn_root + "_ctffind4.log";
	std::ifstream in(fn_log.data(), std::ios_base::in);
	if (in.fail())
    		return false;

	// Start reading the ifstream at the top
	in.seekg(0);
	std::string line;
	std::vector<std::string> words;
	bool found_log = false;
	while (getline(in, line, '\n'))
	{
		// Find the file with the summary of the results
		if (line.find("Summary of results") != std::string::npos)
		{
			tokenize(line, words);
			fn_log = words[words.size() - 1];
			found_log = true;
			break;
		}
	}
	in.close();

	if (!found_log)
		return false;

	// Now open the file with the summry of the results
	std::ifstream in2(fn_log.data(), std::ios_base::in);
	if (in2.fail())
		return false;
	bool Final_is_found = false;
	bool Cs_is_found = false;
	while (getline(in2, line, '\n'))
	{
		// Find data_ lines
		if (line.find("acceleration voltage:") != std::string::npos)
		{
			Cs_is_found = true;
			tokenize(line, words);
			if (words.size() < 19)
				REPORT_ERROR("ERROR: Unexpected number of words on data line with acceleration voltage in " + fn_log);
			CS = textToFloat(words[13]);
			HT = textToFloat(words[8]);
			AmpCnst = textToFloat(words[18]);
			DStep = textToFloat(words[3]);
			XMAG = 10000.;
		}
		else if (line.find("Columns: ") != std::string::npos)
		{
			getline(in2, line, '\n');
			tokenize(line, words);
			if (words.size() < 7)
				REPORT_ERROR("ERROR: Unexpected number of words on data line below Columns line in " + fn_log);
			Final_is_found = true;
			defU = textToFloat(words[1]);
			defV = textToFloat(words[2]);
			defAng = textToFloat(words[3]);
			if (do_phaseshift)
				phaseshift = RAD2DEG(textToFloat(words[4]));
			CC = textToFloat(words[5]);
			if (words[6] == "inf")
				maxres= 999.;
			else
				maxres = textToFloat(words[6]);
		}
	}

	if (!Cs_is_found)
	{
		if (do_warn)
			std::cerr << " WARNING: cannot find line with acceleration voltage etc in " << fn_log << std::endl;
		return false;
	}
	if (!Final_is_found)
	{
		if (do_warn)
			std::cerr << "WARNING: cannot find line with Final values in " << fn_log << std::endl;
		return false;
	}

	in2.close();

	return Final_is_found;
}

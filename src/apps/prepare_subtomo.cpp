/***************************************************************************
 *
 * Author: "Shaoda He"
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

#include <src/args.h>
#include <src/helix.h>
#include <src/funcs.h>
#include <src/macros.h>
#include <src/image.h>

//#define DEBUG

class prepare_subtomo
{
public:
	IOParser parser;

	// Directory of IMOD executables 'extractilts' and 'newstack'
	FileName dir_imod;

	// Tomogram list STAR file
	FileName fn_tomo_list;

	// CTFFIND and Gctf executables
	FileName fn_ctffind_exe, fn_gctf_exe;

	// Extracted particle STAR file
	FileName fn_part;

	// Alias of the particle extraction job
	FileName fn_extract_job_alias;

	bool continue_old, do_skip_ctf_correction;

	bool do_use_trials_for_ctffind, do_use_only_lower_tilt_defoci;

	RFLOAT lower_tilt_defoci_limit;

	RFLOAT bfactor;

	bool is_coords_star_file;

	bool show_usage;

	bool dont_check_input_files;

	////// CTFFIND parameters
	// Size of the box to calculate FFTw
	RFLOAT box_size;

	// Minimum and maximum resolution (in A) to be taken into account
	RFLOAT resol_min, resol_max;

	// Defocus search parameters (in A, positive is underfocus)
	RFLOAT min_defocus, max_defocus, step_defocus;

	// Amount of astigmatism (in A)
	RFLOAT amount_astigmatism;

	// Voltage (kV)
	RFLOAT Voltage;

	// Spherical aberration
	RFLOAT Cs;

	// Amplitude contrast (e.g. 0.07)
	RFLOAT AmplitudeConstrast;

	// Magnification
	RFLOAT Magnification;

	// Detector pixel size (um)
	RFLOAT PixelSize;

	// For Gctf: directly provide angpix!
	RFLOAT angpix;

	////// Additional Gctf Parameters
	bool do_use_gctf, do_ignore_ctffind_params, do_EPA, do_validation;

	std::string additional_gctf_options, gpu_ids;

	prepare_subtomo()
	{
		clear();
	};

	~prepare_subtomo()
	{
		clear();
	};

	void usage()
	{
		parser.writeUsage(std::cerr);
	};

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General");
		show_usage = parser.checkOption("--help", "Show usage");
		dont_check_input_files = parser.checkOption("--dont_check", "(Not recommended) Don't check input files in the initialisation step");
		fn_tomo_list = parser.getOption("--i", "Tomogram STAR file", "all_tomograms.star");
		fn_extract_job_alias = parser.getOption("--o_extract", "Extract job alias. Alias of the job given during particle extraction for RELION 2.0", "extract_tomo");
		dir_imod = parser.getOption("--imod_dir", "Directory of IMOD executables", "/public/EM/imod/imod-4.5.8/IMOD/bin");
		continue_old = parser.checkOption("--only_do_unfinished", "Only extract individual frames, estimate CTFs for those micrographs for which there is not yet a logfile with Final values. Only write out .sh commands where CTF .mrc files have not yet been reconstructed.");
		do_skip_ctf_correction = parser.checkOption("--skip_ctf", "Skip CTF correction? The 3D CTF model will have no CTF modulations, but will still use the Tilt and Bfactor weighting.");
		do_use_trials_for_ctffind = parser.checkOption("--use_trials", "Use trials for CTFFIND. Please keep a tomogram.trial stack in the Tomograms directory containing two trials from either side of the record region. Please note that the tilt order of the files should be same as the aligned stack.");
		do_use_only_lower_tilt_defoci = parser.checkOption("--use_low_tilt", "If you don't have extra trials, then maybe you can set an upper limit of abs(tilt), over which the average defocus value from lower tilts is used."); // default true
		lower_tilt_defoci_limit = textToFloat(parser.getOption("--low_tilt_limit", "(See above)", "30."));
		bfactor = textToFloat(parser.getOption("--bfac", "3D CTF model weighting B-factor per e-/A2", "4."));

		int ctffind_section = parser.addSection("CTFFIND parameters (CTFFIND is used as default)");
		fn_ctffind_exe = parser.getOption("--ctffind_exe", "Location of ctffind executable (or through RELION_CTFFIND_EXECUTABLE environment variable)", "/public/EM/ctffind/ctffind.exe");
		Cs = textToFloat(parser.getOption("--CS", "Spherical Aberration (mm) ","2.7"));
		Voltage = textToFloat(parser.getOption("--HT", "Voltage (kV)","300"));
		AmplitudeConstrast = textToFloat(parser.getOption("--AmpCnst", "Amplitude constrast", "0.07"));
		Magnification = textToFloat(parser.getOption("--XMAG", "Magnification", "53000"));
		PixelSize = textToFloat(parser.getOption("--DStep", "Detector pixel size (um)", "11.57"));

		box_size = textToFloat(parser.getOption("--Box", "Size of the boxes to calculate FFTs", "256"));
		resol_min = textToFloat(parser.getOption("--ResMin", "Minimum resolution (in A) to include in calculations", "50"));
		resol_max = textToFloat(parser.getOption("--ResMax", "Maximum resolution (in A) to include in calculations", "8"));
		min_defocus = textToFloat(parser.getOption("--dFMin", "Minimum defocus value (in A) to search", "20000"));
		max_defocus = textToFloat(parser.getOption("--dFMax", "Maximum defocus value (in A) to search", "50000"));
		step_defocus = textToFloat(parser.getOption("--FStep", "defocus step size (in A) for search", "1000"));
		amount_astigmatism = textToFloat(parser.getOption("--dAst", "amount of astigmatism (in A)", "2000"));

		int gctf_section = parser.addSection("Gctf parameters (ignored if CTFFIND is used)");
		do_use_gctf = parser.checkOption("--use_gctf", "Use Gctf instead of CTFFIND to estimate the CTF parameters");
		fn_gctf_exe = parser.getOption("--gctf_exe", "Location of Gctf executable (or through RELION_GCTF_EXECUTABLE environment variable)", "/lmb/home/kzhang/Public/Gctf/bin/Gctf-v0.50_sm_30_cu7.5_x86_64");
		angpix = textToFloat(parser.getOption("--angpix", "Magnified pixel size in Angstroms", "2.18302"));
		do_ignore_ctffind_params = parser.checkOption("--ignore_ctffind_params", "Use Gctf default parameters instead of CTFFIND parameters");
		do_EPA = parser.checkOption("--EPA", "Use equi-phase averaging to calculate Thon rinds in Gctf");
		do_validation = parser.checkOption("--do_validation", "Use validation inside Gctf to analyse quality of the fit?");
		gpu_ids = parser.getOption("--gpu", "(DOUBLE QUOTES NEEDED) Device ids for each MPI-thread, e.g \"0:1:2:3\"", "");
		additional_gctf_options = parser.getOption("--extra_gctf_options", "(DOUBLE QUOTES NEEDED) Additional options for Gctf (e.g. \"--refine_local_astm\")", "");
	};

	void clear()
	{
		parser.clear();
	};

	void showMessages()
	{
		// Messages
		std::cout << std::endl;
		std::cout << " ### RELION 2.0 sub-tomogram averaging - 23:59, FEB 19, 2014 ###" << std::endl;
		std::cout << " # The original python script was written by Tanmay Bharat to support sub-tomogram averaging in RELION." << std::endl;
		std::cout << " # This 'relion_prepare_subtomo' executable was written by Shaoda He in Sjors Scheres' lab." << std::endl;
		std::cout << " # Please ensure that you have provided the directory containing IMOD executables 'extracttilts' and 'newstack'" << std::endl;
		std::cout << " # Please provide either CTFFIND or Gctf executable." << std::endl;
		std::cout << " # Please report bugs and comments to tbharat@mrc-lmb.cam.ac.uk or scheres@mrc-lmb.cam.ac.uk" << std::endl;
		std::cout << " # Please read the documentation on the RELION wiki, several questions are answered there." << std::endl;
		std::cout << " # This version can set defocus values above a certain tilt to the defocus value of the zero degree tilt." << std::endl;
		std::cout << " # This version will write out all the CTF reconstruction commands in the master file." << std::endl;
		std::cout << " # This version supports RELION 2.0 only. For compatibility with older RELION, please use the original python script." << std::endl;
		std::cout << " # This version depends on IMOD executables (extracttilts, newstack) and CTFFIND or Gctf." << std::endl;
		std::cout << std::endl;
		std::cout << " ### RELION 2.0 sub-tomogram averaging - Usage (also refer to RELION wiki) ###" << std::endl;
		std::cout << " # Before running the program: " << std::endl;
		std::cout << " # 1. Create a directory 'Tomogram/tomo\?\?\?' for each reconstructed 3D tomogram." << std::endl;
		std::cout << " # 2. In each of the individual tomogram directories you need:" << std::endl;
		std::cout << " #    a. tomo.mrc	   : the actual reconstructed tomogram." << std::endl;
		std::cout << " #    b. tomo.mrcs   : the aligned tilt series in MRC-stack format (Please rename if they are in .st format!)" << std::endl;
		std::cout << " #    c. tomo.star   : a STAR file with at least 3 columns: _rlnCoordinateX, Y and Z. (e.g. STAR file generated by 'relion_helix_toolbox --interpo')" << std::endl;
		std::cout << " #        OR  (if STAR file exists then .coords file will be ignored)" << std::endl;
		std::cout << " #       tomo.coords : a text file with 3 columns: the X, Y and Z coordinates of each subtomogram (e.g. save this from IMOD)." << std::endl;
		std::cout << " #    d. tomo.order  : a text file with 2 columns: the tilt angle of each image in tomo.mrcs and the accumulated dose in e-/A2 for that image." << std::endl;
		std::cout << " #    e. tomo.tlt    : (OPTIONAL) a text file with the final tilt angles from IMOD. If this is not provided then the extended header of the .mrcs will be read." << std::endl;
		std::cout << " # 3. Run the program. (Input files will be checked in the initialisation step. Please pay attention if error messages pop up.)" << std::endl;
		std::cout << " # 4. Check the contents of 'do_all_reconstruct_ctfs.sh', (split it into multiple files for parallelisation) and run the .sh script (please provide the reconstruction box size)." << std::endl;
		std::cout << " # 5. Process the data with RELION 2.0 GUI." << std::endl;
		std::cout << std::endl;
		if (show_usage)
			REPORT_ERROR("All the available parameters and their default values are listed above.");
	};

	void initialChecks()
	{
		FileName fn1, fn2, fn3;
		std::ifstream fin1;
		MetaDataTable MD_tomo;
		std::string line, command;
		std::vector<std::string> words;
		std::vector<FileName> fns_tomo;
		RFLOAT calc_angpix = 0.;
		int res = 0;

		std::cout << " ###################################################################" << std::endl;
		std::cout << " Checking input data ..." << std::endl;
		if (Magnification < 1.)
			REPORT_ERROR("Invalid magnification!");
		calc_angpix = 10000. * PixelSize / Magnification;
		std::cout << " Calculated pixel size (10000 * DPix / Mag) = " << calc_angpix << " Angstrom(s)" << std::endl;
		if (calc_angpix < 0.001)
			REPORT_ERROR("Calculated pixel size is smaller than 0.001!");

		// Check CTFFIND or Gctf executables
		if (do_skip_ctf_correction)
		{
			Cs = 0.;
			AmplitudeConstrast = 1.;
			do_use_trials_for_ctffind = false;
			do_use_only_lower_tilt_defoci = false;
		}
		else
		{
			if (do_use_gctf)
			{
				//REPORT_ERROR("Gctf is not currently supported!");
				// Get the GCTF executable
				if (fn_gctf_exe == "")
				{
					char* penv = getenv ("RELION_GCTF_EXECUTABLE");
					if (penv != NULL)
						fn_gctf_exe = (std::string)penv;
				}
				if ( (fn_gctf_exe.length() < 2) && (!exists(fn_gctf_exe)) )
					REPORT_ERROR("Cannot find Gctf executable " + fn_gctf_exe);
#ifdef DEBUG
				if (gpu_ids != "")
				{
					std::cout << " str_gpu_ids = " << gpu_ids << std::endl;
					//if (gpu_ids.length() < 2)
					//	REPORT_ERROR("Invalid GPU ids!");
					//if ( (gpu_ids[0] != '\"') || (gpu_ids[gpu_ids.length() - 1] != '\"') )
					//	REPORT_ERROR("GPU ids should come with double quotes outside! (e.g. \"0:1:2:3\")");
				}
				if (additional_gctf_options != "")
				{
					std::cout << " str_additional_gctf_options = " << additional_gctf_options << std::endl;
					//if (additional_gctf_options.length() < 2)
					//	REPORT_ERROR("Invalid additional gctf options!");
					//if ( (additional_gctf_options[0] != '\"') || (additional_gctf_options[additional_gctf_options.length() - 1] != '\"') )
					//	REPORT_ERROR("Additional gctf options should come with double quotes outside! (e.g. \"--refine_local_astm\")");
				}
				return RELION_EXIT_FAILURE;
#endif
				std::cout << " Pixel size used in Gctf = " << angpix << " Angstrom(s)" << std::endl;
				RFLOAT ratio = angpix / calc_angpix; // calc_angpix >= 0.001, no need to check again
				if ( (ratio < 0.99) || (ratio > 1.01) )
					REPORT_ERROR("Calculated and user-defined pixel sizes are different (> 1% of error)!");
			}
			else
			{
				// Get the CTFFIND executable
				if (fn_ctffind_exe == "")
				{
					char* penv = getenv ("RELION_CTFFIND_EXECUTABLE");
					if (penv != NULL)
						fn_ctffind_exe = (std::string)penv;
				}
				if ( (fn_ctffind_exe.length() < 2) && (!exists(fn_ctffind_exe)) )
					REPORT_ERROR("Cannot find CTFFIND executable " + fn_ctffind_exe);
			}
			// TODO: output CTF parameters!
		}

		if (do_use_only_lower_tilt_defoci)
		{
			if (!(lower_tilt_defoci_limit > 0.))
				REPORT_ERROR("Lower tilt defoci limit should be larger than 0.0!");
			// TODO: Report lower_tilt_defoci_limit
		}

		// Check IMOD executables: 'extracttilts' and 'newstack
		fn1 = dir_imod + "/extracttilts";
		if (!exists(fn1))
			REPORT_ERROR("Cannot find IMOD 'extractilts' executable " + fn1);
		fn1 = dir_imod + "/newstack";
		if (!exists(fn1))
			REPORT_ERROR("Cannot find IMOD 'newstack' executable " + fn1);

		// Check tomogram list
		if (!exists(fn_tomo_list))
			REPORT_ERROR("Cannot find the STAR file with all the tomograms " + fn_tomo_list);
		MD_tomo.clear();
		MD_tomo.read(fn_tomo_list);
		if (MD_tomo.numberOfObjects() < 1)
			REPORT_ERROR("Tomogram STAR file " + fn_tomo_list + " is empty!");
		if (!MD_tomo.containsLabel(EMDL_MICROGRAPH_NAME))
			REPORT_ERROR("Tomogram STAR file " + fn_tomo_list + " should contain _rlnMicrographName!");
		// Check whether each tomogram sits in a separate folder
		fns_tomo.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_tomo)
		{
			MD_tomo.getValue(EMDL_MICROGRAPH_NAME, fn1);
			fns_tomo.push_back(fn1);
		}
		if (fns_tomo.size() < 1)
			REPORT_ERROR("Tomogram STAR file " + fn_tomo_list + " is empty!");
		std::stable_sort(fns_tomo.begin(), fns_tomo.end());
		for (size_t idx = 0; idx < (fns_tomo.size() - 1); idx++)
		{
			fn1 = fns_tomo[idx];
			fn2 = fns_tomo[idx + 1];
			if (fn1 == fn2)
				REPORT_ERROR("Tomogram " + fn1 + " appears in the tomogram STAR file more than once!");
			if (fn1.beforeLastOf("/") == fn2.beforeLastOf("/"))
				REPORT_ERROR("Tomograms " + fn1 + " and " + fn2 + " are located in the same folder!");
		}

		// Check dependent files
		MetaDataTable MD_tmp;
		Image<RFLOAT> img;
		int xdim = 0, ydim = 0, zdim = 0;
		long int ndim = 0, nr_frames = 0, nr_lines = 0;
		RFLOAT xx = 0., yy = 0., zz = 0.;
		bool is_star_coords = false, is_txt_coords = false;
		MD_tmp.clear();
		img.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_tomo)
		{
			MD_tomo.getValue(EMDL_MICROGRAPH_NAME, fn1);
			std::cout << " 3D reconstructed tomogram in STAR file: " << fn1 << std::flush;

			// Check 3D reconstructed tomogram
			if (!exists(fn1))
				REPORT_ERROR("Cannot find 3D reconstructed tomogram " + fn1);
			img.read(fn1, false);
			img.getDimensions(xdim, ydim, zdim, ndim);
			std::cout << " , Dimensions XYZN = " << xdim << " * " << ydim << " * " << zdim << " * " << ndim << std::endl;
			if ( (zdim > 1) && (ndim > 1) )
				REPORT_ERROR("Reconstructed 3D tomogram " + fn1 + " is 4D!");
			if ( (xdim < box_size) || (ydim < box_size) )
				REPORT_ERROR("X and/or Y dimensions of reconstructed 3D tomogram " + fn1 + " is smaller than CTF box_size " + integerToString(box_size));
			// TODO: consider Gctf box_size?

			// Check original tilt series
			fn2 = fn1.withoutExtension() + ".mrcs";
			if (!exists(fn2))
				REPORT_ERROR("Cannot find original tilt series (.mrcs) " + fn2);
			std::cout << " Tilt series  : " << fn2 << std::flush;
			img.read(fn2, false);
			img.getDimensions(xdim, ydim, zdim, ndim);
			std::cout << " , Dimensions XYZN = " << xdim << " * " << ydim << " * " << zdim << " * " << ndim << std::endl;
			if ( (zdim > 1) && (ndim > 1) )
				REPORT_ERROR("Tilt series " + fn2 + " is 4D!");
			if ( (xdim < box_size) || (ydim < box_size) )
				REPORT_ERROR("X and/or Y dimensions of tilt series " + fn2 + " is smaller than CTF box_size " + integerToString(box_size));
			nr_frames = zdim * ndim;
			if (nr_frames < 2)
				REPORT_ERROR("Tilt series " + fn2 + " contains less than 2 frames!");

			// Check .tlt (optional)
			fn2 = fn1.withoutExtension() + ".tlt";
			if (exists(fn2))
			{
				std::cout << "  File .tlt    : " << fn2 << "  (optional)" << std::endl;
				nr_lines = 0;
				fin1.open(fn2.c_str(), std::ios_base::in);
				if (fin1.fail())
					REPORT_ERROR("Cannot open .tlt file: " + fn2);
				while (getline(fin1, line, '\n'))
				{
					words.clear();
					tokenize(line, words);
					if (words.size() == 0) // Empty line
						continue;
					if (words.size() != 1) // 1 blocks: tilt angle
						REPORT_ERROR("Invalid .tlt file: " + fn2);
					xx = textToFloat(words[0]);
					nr_lines++; // A valid line
				}
				fin1.close();
				if (nr_lines != nr_frames)
					REPORT_ERROR("Tilt series has " + integerToString(nr_frames) + " frames but .tlt file " + fn2 + " has " + integerToString(nr_lines) + " lines!");
			}

			// Check .order
			fn2 = fn1.withoutExtension() + ".order";
			if (!exists(fn2))
				REPORT_ERROR("Cannot find .order file " + fn2);
			std::cout << "  File .order  : " << fn2 << std::endl;
			nr_lines = 0;
			fin1.open(fn2.c_str(), std::ios_base::in);
			if (fin1.fail())
				REPORT_ERROR("Cannot open .order file: " + fn2);
			while (getline(fin1, line, '\n'))
			{
				words.clear();
				tokenize(line, words);
				if (words.size() == 0) // Empty line
					continue;
				if (words.size() != 2) // 2 blocks: tilt angle, accumulated dose
					REPORT_ERROR("Invalid .order file: " + fn2);
				xx = textToFloat(words[0]);
				yy = textToFloat(words[1]);
				nr_lines++; // A valid line
			}
			fin1.close();
			if (nr_lines != nr_frames)
				REPORT_ERROR("Tilt series has " + integerToString(nr_frames) + " frames but .order file " + fn2 + " has " + integerToString(nr_lines) + " lines!");

			// Check txt or STAR coords
			fn2 = fn1.withoutExtension() + ".coords";
			fn3 = fn1.withoutExtension() + ".star";
			if ( (!exists(fn2)) && (!exists(fn3)) )
				REPORT_ERROR("Cannot find .coord OR .star file " + fn2 + " OR " + fn3);
			if (exists(fn3))
			{
				std::cout << "  Coords STAR  : " << fn3 << std::endl;
				MetaDataTable MD_this_tomo;
				MD_this_tomo.clear();
				MD_this_tomo.read(fn3);
				if (MD_this_tomo.numberOfObjects() < 1)
					REPORT_ERROR("Coordinates STAR file " + fn3 + " is empty!");
				if (MD_tmp.numberOfObjects() < 1) // MD_tmp is empty. 'fn3' is the first STAR file ever read.
					MD_tmp = MD_this_tomo;
				else
				{
					if (!MetaDataTable::compareLabels(MD_tmp, MD_this_tomo))
						REPORT_ERROR("Coordinates STAR file " + fn3 + " has a different set of activeLabels!");
				}
				if ( (!MD_this_tomo.containsLabel(EMDL_IMAGE_COORD_X))
						|| (!MD_this_tomo.containsLabel(EMDL_IMAGE_COORD_Y))
						|| (!MD_this_tomo.containsLabel(EMDL_IMAGE_COORD_Z)) )
					REPORT_ERROR("Coordinates STAR file " + fn3 + " should contain _rlnCoordinateX Y and Z!");

				is_star_coords = true;
				if (is_txt_coords)
					REPORT_ERROR("All coordinates files should have the same extensions (either all .star or all .coords)!");
			}
			else
			{
				std::cout << "  File .coords : " << fn2 << std::endl;
				nr_lines = 0;
				fin1.open(fn2.c_str(), std::ios_base::in);
				if (fin1.fail())
					REPORT_ERROR("Cannot open .coords file: " + fn2);
				while (getline(fin1, line, '\n'))
				{
					words.clear();
					tokenize(line, words);
					if (words.size() == 0) // Empty line
						continue;
					if (words.size() != 3) // 3 blocks: x, y, z
						REPORT_ERROR("Invalid .coords file: " + fn2);
					xx = textToFloat(words[0]);
					yy = textToFloat(words[1]);
					zz = textToFloat(words[2]);
					nr_lines++; // A valid line
				}
				fin1.close();

				is_txt_coords = true;
				if (is_star_coords)
					REPORT_ERROR("All coordinates files should have the same extensions (either all .star or all .coords)!");
			}

			// Check .trial
			if (do_use_trials_for_ctffind)
			{
				fn2 = fn1.withoutExtension() + ".trial";
				if (!exists(fn2))
					REPORT_ERROR("Cannot find .trial file " + fn2);
				std::cout << "  File .trial  : " << fn2 << std::flush;
				fn2 += ":mrcs"; // Open this file as .mrcs stack
				img.read(fn2, false);
				img.getDimensions(xdim, ydim, zdim, ndim);
				std::cout << " , Dimensions XYZN = " << xdim << " * " << ydim << " * " << zdim << " * " << ndim << std::endl;
				if ( (zdim > 1) && (ndim > 1) )
					REPORT_ERROR("Trial series " + fn2 + " is 4D!");
				if ( (xdim < box_size) || (ydim < box_size) )
					REPORT_ERROR("X and/or Y dimensions of trial series " + fn2 + " is smaller than CTF box_size " + integerToString(box_size));
				if ( (zdim * ndim) != (nr_frames * 2))
					REPORT_ERROR("Trial series has " + integerToString(zdim * ndim) + " frames, not 2X the total frames " + integerToString(nr_frames) + " in the tilt series!");
			}

			// Create folders for CTF correction
			fn2 = fn1.beforeLastOf("/") + "/Ctffind/Results";
			std::cout << " Folder containing CTF correction results: " << fn2 << std::endl;
			if (!exists(fn2))
			{
				std::cout << " This folder does not exist. Create it." << std::endl;
				command = "mkdir -p " + fn2;
#ifdef DEBUG
				std::cout << " " << command << std::endl;
#endif
				res = system(command.c_str());
			}

			// Create folders for particle extraction
			fn2 = "Particles/" + fn1.beforeLastOf("/");
			std::cout << " Folder containing extracted particles: " << fn2 << std::endl;
			if (!exists(fn2))
			{
				std::cout << " This folder does not exist. Create it." << std::endl;
				command = "mkdir -p " + fn2;
#ifdef DEBUG
				std::cout << " " << command << std::endl;
#endif
				res = system(command.c_str());
			}
		}

		std::cout << " Input data checked." << std::endl;
		std::cout << " ###################################################################" << std::endl;
	};

	void run()
	{
		std::ifstream fin1;
		std::ofstream fout1;
		std::string line, command;
		std::vector<std::string> words;
		int res = 0;
		FileName fn_ctf_recon = "do_all_reconstruct_ctfs.sh";

		MetaDataTable MD_tomo, MD_part;
		MD_tomo.clear();
		MD_tomo.read(fn_tomo_list);

		// Initialise MetaDataTable for particles
		MD_part.clear();
		MD_part.addLabel(EMDL_MICROGRAPH_NAME);
		MD_part.addLabel(EMDL_IMAGE_COORD_X);
		MD_part.addLabel(EMDL_IMAGE_COORD_Y);
		MD_part.addLabel(EMDL_IMAGE_COORD_Z);
		MD_part.addLabel(EMDL_IMAGE_NAME);
		MD_part.addLabel(EMDL_CTF_IMAGE);
		MD_part.addLabel(EMDL_CTF_MAGNIFICATION);
		MD_part.addLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE);

		// Open output file fn_ctf_recon
		fout1.open(fn_ctf_recon.c_str(), std::ios_base::out);
		if (fout1.fail())
			REPORT_ERROR("Cannot open output file: " + (std::string)(fn_ctf_recon));
		if (continue_old)
		{
			fout1 << "# Option '--only_do_unfinished' is enabled." << std::endl;
			fout1 << "# Commands are commented if reconstructed CTF files exist." << std::endl;
			fout1 << std::endl;
		}

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_tomo)
		{
			std::cout << std::endl;

			FileName fn_tomo;
			MD_tomo.getValue(EMDL_MICROGRAPH_NAME, fn_tomo);
			std::cout << "#### Processing tomogram " << fn_tomo << " ... ####" << std::endl;

			FileName dir_ctf = fn_tomo.beforeLastOf("/") + "/Ctffind";
			FileName dir_ctf_results = fn_tomo.beforeLastOf("/") + "/Ctffind/Results";
			FileName fn_order = fn_tomo.withoutExtension() + ".order";
			FileName fn_stack = fn_tomo.withoutExtension() + ".mrcs";
			FileName fn_trial = fn_tomo.withoutExtension() + ".trial"; // optional
			FileName fn_tilt = fn_tomo.withoutExtension() + ".tlt"; // optional
			FileName fn_coords = fn_tomo.withoutExtension() + ".star";
			if (!exists(fn_coords))
				fn_coords = fn_tomo.withoutExtension() + ".coords";
			bool have_tilt = exists(fn_tilt);

			// ExtractTilt
			if (!have_tilt)
			{
				std::cout << " File .tlt does not exist. Call IMOD extracttilt." << std::endl;
				command = dir_imod + "/extracttilts -InputFile " + fn_stack + " -tilts -OutputFile " + dir_ctf + "/tiltangles.txt > " + dir_ctf + "/tiltangles_scratch.txt";
//#ifdef DEBUG
				std::cout << " " << command << std::endl;
//#endif
				res = system(command.c_str());
			}
			else
			{
				std::cout << " File .tlt file exists. Copy it to CTFFIND folder." << std::endl;
				command = "cp " + fn_tilt + " " + dir_ctf + "/tiltangles.txt";
//#ifdef DEBUG
				std::cout << " " << command << std::endl;
//#endif
				res = system(command.c_str());
			}
			if (do_use_trials_for_ctffind)
			{
				std::cout << " Use trails in CTFFIND. Call IMOD extracttilt." << std::endl;
				command = dir_imod + "/extracttilts -InputFile " + fn_trial + " -tilts -OutputFile " + dir_ctf + "/trial_tiltangles.txt > " + dir_ctf + "/trial_tiltangles_scratch.txt";
//#ifdef DEBUG
				std::cout << " " << command << std::endl;
//#endif
				res = system(command.c_str());
			}
			std::cout << " Tilt values have been extracted." << std::endl;

			// CTF correction
			//FileName fn_ctf = fn_tomo.withoutExtension() + "_images.star"; // OLD
			FileName fn_ctf = fn_tomo.beforeLastOf("/") + "/Ctffind/" + (fn_tomo.afterLastOf("/")).withoutExtension() + "_images.star"; // NEW
#ifdef DEBUG
			std::cout << " fn_ctf = " << fn_ctf << std::endl;
#endif
			FileName fn_tilt_txt = dir_ctf + "/tiltangles.txt";
			FileName fn_trial_tilt_txt = dir_ctf + "/trial_tiltangles.txt";
			std::vector<RFLOAT> tilts, trial_tilts;
			MetaDataTable MD_ctf, MD_ctf_results;

			tilts.clear(); trial_tilts.clear();
			MD_ctf.clear();
			MD_ctf_results.clear();

			MD_ctf.addLabel(EMDL_MICROGRAPH_NAME);
			fin1.open(fn_tilt_txt.c_str(), std::ios_base::in);
			if (fin1.fail())
				REPORT_ERROR("Cannot open input file: " + (std::string)(fn_tilt_txt));

			// Get tilt angles (without trials)
			while (getline(fin1, line, '\n'))
			{
				words.clear();
				tokenize(line, words);
				if (words.size() == 0) // Empty line
					continue;
				if (words.size() != 1)
					REPORT_ERROR("Invalid .tlt file: " + fn_tilt_txt);
				tilts.push_back(textToFloat(words[0]));
			}
			fin1.close();

			// Get tilt angles (with trials)
			if (do_use_trials_for_ctffind)
			{
				fin1.open(fn_trial_tilt_txt.c_str(), std::ios_base::in);
				if (fin1.fail())
					REPORT_ERROR("Cannot open input file: " + (std::string)(fn_trial_tilt_txt));
				while (getline(fin1, line, '\n'))
				{
					words.clear();
					tokenize(line, words);
					if (words.size() == 0) // Empty line
						continue;
					if (words.size() != 1)
						REPORT_ERROR("Invalid .tlt trial file: " + fn_trial_tilt_txt);
					trial_tilts.push_back(textToFloat(words[0]));
				}
				fin1.close();
				if ( (tilts.size() * 2) != trial_tilts.size())
					REPORT_ERROR("Invalid .tlt and/or .tlt trial files: " + fn_tilt_txt + " and " + fn_trial_tilt_txt + " . The trial file should contain 2X lines.");
			}

			// Extract individual tilt frames using IMOD executable 'newstack' or 'relion_image_handler'
			// Note that frame id starts from 0 in' newstack' but 1 in 'relion_image_handler'
			command = "touch " + dir_ctf + "/temp_newstack_out.txt";
#ifdef DEBUG
			std::cout << " " << command << std::endl;
#endif
			for (int ida = 0; ida < tilts.size(); ida++)
			{
				FileName fn_sec;
				if (do_use_trials_for_ctffind)
				{
					// Command 1
					fn_sec = dir_ctf + "/" + (fn_trial.afterLastOf("/")).withoutExtension() + "_image" + floatToString(tilts[ida]) + "_" + integerToString(2 * ida) + ".mrc";
					if ( (continue_old) && (exists(fn_sec)) )
					{}
					else
					{
						// TODO: I want to use relion_image_handler but failed. It does not produce the same set of individual frames.
						command = "relion_image_handler --i " + integerToString(2 * ida + 1) + "@" + fn_trial + ":mrcs --o " + fn_sec;
						//command = dir_imod + "/newstack -secs " + integerToString(2 * ida) + " " + fn_trial + " " + fn_sec + " >> " + dir_ctf + "/temp_newstack_out.txt";
//#ifdef DEBUG
						std::cout << " " << command << std::endl;
//#endif
						res = system(command.c_str());
					}
					MD_ctf.addObject();
					MD_ctf.setValue(EMDL_MICROGRAPH_NAME, fn_sec);

					// Command 2
					fn_sec = dir_ctf + "/" + (fn_trial.afterLastOf("/")).withoutExtension() + "_image" + floatToString(tilts[ida]) + "_" + integerToString(2 * ida + 1) + ".mrc";
					if ( (continue_old) && (exists(fn_sec)) )
					{}
					else
					{
						// TODO: I want to use relion_image_handler but failed. It does not produce the same set of individual frames.
						command = "relion_image_handler --i " + integerToString(2 * ida + 2) + "@" + fn_trial + ":mrcs --o " + fn_sec;
						//command = dir_imod + "/newstack -secs " + integerToString(2 * ida + 1) + " " + fn_trial + " " + fn_sec + " >> " + dir_ctf + "/temp_newstack_out.txt";
//#ifdef DEBUG
						std::cout << " " << command << std::endl;
//#endif
						res = system(command.c_str());
					}
					MD_ctf.addObject();
					MD_ctf.setValue(EMDL_MICROGRAPH_NAME, fn_sec);
				}
				else
				{
					fn_sec = dir_ctf + "/" + (fn_stack.afterLastOf("/")).withoutExtension() + "_image" + floatToString(tilts[ida]) + "_" + integerToString(ida) + ".mrc";
					if ( (continue_old) && (exists(fn_sec)) )
					{}
					else
					{
						// TODO: I want to use relion_image_handler but failed. It does not produce the same set of individual frames.
						command = "relion_image_handler --i " + integerToString(ida + 1) + "@" + fn_stack + ":mrcs --o " + fn_sec;
						//command = dir_imod + "/newstack -secs " + integerToString(ida) + " " + fn_stack + " " + fn_sec + " >> " + dir_ctf + "/temp_newstack_out.txt";
//#ifdef DEBUG
						std::cout << " " << command << std::endl;
//#endif
						res = system(command.c_str());
					}
					MD_ctf.addObject();
					MD_ctf.setValue(EMDL_MICROGRAPH_NAME, fn_sec);
				}
			}
			MD_ctf.write(fn_ctf);

			// Run CTFFIND and store the estimated parameters
			if (do_skip_ctf_correction)
			{
				MD_ctf_results.clear();
				MD_ctf_results = MD_ctf;
				MD_ctf_results.addLabel(EMDL_CTF_DEFOCUSU);
				MD_ctf_results.addLabel(EMDL_CTF_DEFOCUSV);
				//MD_ctf_results.addLabel(EMDL_CTF_DEFOCUS_ANGLE);
				FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_ctf_results)
				{
					MD_ctf_results.setValue(EMDL_CTF_DEFOCUSU, 0.);
					MD_ctf_results.setValue(EMDL_CTF_DEFOCUSV, 0.);
					//MD_ctf_results.setValue(EMDL_CTF_DEFOCUS_ANGLE, 0.);
				}
				MD_ctf_results.write(dir_ctf_results + "/micrographs_ctf.star");
			}
			else
			{
				command = "relion_run_ctffind --i " + fn_ctf + " --o " + dir_ctf_results + "/ --CS " + floatToString(Cs)
						+ " --HT " + floatToString(Voltage) + " --ctfWin -1 --AmpCnst " + floatToString(AmplitudeConstrast)
						+ " --DStep " + floatToString(PixelSize) + " --XMAG " + floatToString(Magnification)
						+ " --Box " + floatToString(box_size) + " --dFMin " + floatToString(min_defocus)
						+ " --dFMax " + floatToString(max_defocus) + " --FStep " + floatToString(step_defocus)
						+ " --dAst " + floatToString(amount_astigmatism) + " --ResMin " + floatToString(resol_min)
						+ " --ResMax " + floatToString(resol_max);
				if (continue_old)
					command += " --only_do_unfinished";

				if (do_use_gctf)
				{
					command += " --use_gctf --gctf_exe " + fn_gctf_exe;
					command += " --angpix " + floatToString(angpix);
					if (do_ignore_ctffind_params)
						command += " --ignore_ctffind_params";
					if (do_EPA)
						command += " --EPA";
					if (do_validation)
						command += " --do_validation";
					//if (gpu_ids.length() < 2)
					//	gpu_ids = "\"\""; // If gpu_ids is empty, put double quotes outside
					command += " --gpu \"" + gpu_ids + "\""; // TODO: User needs to provide double quotes outside
					if (additional_gctf_options.length() > 0)
						command += " --extra_gctf_options \"" + additional_gctf_options + "\""; // TODO: User needs to provide double quotes outside
				}
				else
				{
					command += " --ctffind_exe \"" + fn_ctffind_exe + "  --omp-num-threads 1 --old-school-input\"";
				}
//#ifdef DEBUG
				std::cout << " " << command << std::endl;
//#endif
				res = system(command.c_str());

				// TODO: Support GCTF ?
			}

			// Making .star files for each 3D CTF Volume
			//FileName dir_rec = "Particles/" + fn_tomo.beforeLastOf("/"); // I don't need this...
			//FileName fn_rec = "Particles/" + fn_tomo.withoutExtension() + "_rec_CTF_volumes.sh"; // I don't need this...
			std::vector<RFLOAT> order_tilts, accu_dose, avg_defoci;
			RFLOAT du = 0., dv = 0.;

			order_tilts.clear(); accu_dose.clear(); avg_defoci.clear();
			std::cout << " Reading tilt series order file " + fn_order + " for dose dependent B-Factor weighting..." << std::endl;

			fin1.open(fn_order.c_str(), std::ios_base::in);
			if (fin1.fail())
				REPORT_ERROR("Cannot open input file: " + (std::string)(fn_order));
			while (getline(fin1, line, '\n'))
			{
				words.clear();
				tokenize(line, words);
				if (words.size() == 0) // Empty line
					continue;
				if (words.size() != 2)
					REPORT_ERROR("Invalid input file!");
				order_tilts.push_back(textToFloat(words[0]));
				accu_dose.push_back(textToFloat(words[1]));
			}
			fin1.close();

			MD_ctf_results.clear();
			MD_ctf_results.read(dir_ctf_results + "/micrographs_ctf.star");
			if ( (!MD_ctf_results.containsLabel(EMDL_CTF_DEFOCUSU))
					|| (!MD_ctf_results.containsLabel(EMDL_CTF_DEFOCUSV))
					|| (!MD_ctf_results.containsLabel(EMDL_MICROGRAPH_NAME)) )
				REPORT_ERROR("micrographs_ctf.star should contain _rlnDefocusU, _rlnDefocusV and _rlnMicrographName! Please check whether CTF estimation was done successfully.");
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_ctf_results)
			{
				MD_ctf_results.getValue(EMDL_CTF_DEFOCUSU, du);
				MD_ctf_results.getValue(EMDL_CTF_DEFOCUSV, dv);
				avg_defoci.push_back(du); // TODO: Why just read in defocusU but not with defocusV and defocusAngle ???
			}
			if (do_use_trials_for_ctffind)
			{
#ifdef DEBUG
				std::cout << " trial_tilts.size() = " << trial_tilts.size() << " , order_tilts.size() = " << order_tilts.size() << std::endl;
#endif
				if (trial_tilts.size() != (order_tilts.size() * 2))
					REPORT_ERROR("Tilt series are different in .tlt (or stack header) and .order files: " + fn_trial_tilt_txt + " and " + fn_order);
#ifdef DEBUG
				std::cout << " avg_defoci.size() = " << avg_defoci.size() << " , trial_tilts.size() = " << trial_tilts.size() << std::endl;
#endif
				if (avg_defoci.size() != trial_tilts.size())
					REPORT_ERROR("Tilt series are different in .tlt (or stack header) and micrographs_ctf.star files: " + fn_trial_tilt_txt + " and " + dir_ctf_results + "/micrographs_ctf.star");
				std::vector<RFLOAT> tmp_vec;
				tmp_vec.clear();
				for (int id = 0; id < avg_defoci.size(); id += 2)
					tmp_vec.push_back((avg_defoci[id] + avg_defoci[id + 1]) / 2.);
				avg_defoci.clear();
				avg_defoci = tmp_vec;
			}
			else
			{
				if (tilts.size() != order_tilts.size())
					REPORT_ERROR("Tilt series are different in .tlt (or stack header) and .order files: " + fn_tilt_txt + " and " + fn_order);
				if (avg_defoci.size() != tilts.size())
					REPORT_ERROR("Tilt series are different in .tlt (or stack header) and micrographs_ctf.star files!");
			}
			if (avg_defoci.size() != tilts.size())
				REPORT_ERROR("Tilt series are different in .tlt (or stack header) and micrographs_ctf.star files: " + fn_tilt_txt + " and " + dir_ctf_results + "/micrographs_ctf.star");

			// Deal with lower tilts
			if (do_use_only_lower_tilt_defoci)
			{
				RFLOAT sum_defoci = 0., nr_defoci = 0.;
				for (int id = 0; id < tilts.size(); id++)
				{
					if (fabs(tilts[id]) < lower_tilt_defoci_limit)
					{
						sum_defoci += avg_defoci[id];
						nr_defoci += 1.;
					}
				}
				if (nr_defoci > 0.5)
					sum_defoci /= nr_defoci;
				std::cout << " " << fn_tomo << " : Average defocus from the lower tilt images below " << lower_tilt_defoci_limit << " is " << sum_defoci << std::endl;
				for (int id = 0; id < tilts.size(); id++)
				{
					// TODO:
					if (fabs(tilts[id]) > lower_tilt_defoci_limit) // python script
					{
					//if (!(fabs(tilts[id]) < lower_tilt_defoci_limit)) // Shaoda's idea
						avg_defoci[id] = sum_defoci;
					}
				}
			}

			// TODO: consider Y/Z flipped tomograms ? Maybe not. The coordinates have already been flipped when processing the .mod files.
			Image<RFLOAT> img;
			int xdim = 0, ydim = 0, zdim = 0;
			long int ndim = 0;
			RFLOAT calc_angpix = 10000. * PixelSize / Magnification;
			std::cout << " Calculated pixel size = " << calc_angpix << " Angstrom(s)" << std::endl;
			std::cout << " Extract XYZN dimensions of the tomogram " << fn_tomo << std::endl;
			img.read(fn_tomo, false);
			img.getDimensions(xdim, ydim, zdim, ndim);
			std::cout << " Tomogram XYZN dimensions = " << xdim << " * " << ydim << " * " << zdim << " * " << ndim << std::endl;
			std::cout << " Writing out .star files to make 3D CTF volumes..." << std::endl;

			RFLOAT xx = 0., yy = 0., zz = 0.;
			int nr_subtomo = 0;
			bool write_star_file = false;
			FileName fn_subtomo_star, fn_subtomo_mrc;
			MetaDataTable MD_coords, MD_this_subtomo;

			// Load coordinates into MD_coords
			MD_coords.clear();
#ifdef DEBUG
			std::cout << " fn_coords = " << fn_coords << std::endl;
#endif
			if (fn_coords.getExtension() != "star")
			{
				MD_coords.addLabel(EMDL_IMAGE_COORD_X);
				MD_coords.addLabel(EMDL_IMAGE_COORD_Y);
				MD_coords.addLabel(EMDL_IMAGE_COORD_Z);

				fin1.open(fn_coords.c_str(), std::ios_base::in);
				while (getline(fin1, line, '\n'))
				{
					words.clear();
					tokenize(line, words);
					if (words.size() == 0) // Empty line
						continue;
					if (words.size() != 3)
						REPORT_ERROR("Invalid input file: " + fn_coords);

					MD_coords.addObject();
					MD_coords.setValue(EMDL_IMAGE_COORD_X, textToFloat(words[0]));
					MD_coords.setValue(EMDL_IMAGE_COORD_Y, textToFloat(words[1]));
					MD_coords.setValue(EMDL_IMAGE_COORD_Z, textToFloat(words[2]));
				}
				fin1.close();
			}
			else
			{
				// All MD_coords have the same set of EMDLabels
				MD_coords.read(fn_coords);

				// Append extra columns from MD_coords
				std::vector<EMDLabel> labels = MD_coords.getActiveLabels();
				for (size_t idx = 0; idx < labels.size(); idx++)
				{
					if (!MD_part.containsLabel(labels[idx]))
						MD_part.addLabel(labels[idx]);
				}
			}

			// Loop over every picked 3D point
			if (MD_coords.numberOfObjects() < 1)
				REPORT_ERROR("MD_coords is empty! It reads from .coord or .star file: " + fn_coords);
			nr_subtomo = 0;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD_coords)
			{
				nr_subtomo++;
				MD_coords.getValue(EMDL_IMAGE_COORD_X, xx);
				MD_coords.getValue(EMDL_IMAGE_COORD_Y, yy);
				MD_coords.getValue(EMDL_IMAGE_COORD_Z, zz);

				write_star_file = ((do_skip_ctf_correction) && (nr_subtomo == 1)) || (!do_skip_ctf_correction);
				if (!write_star_file)
					continue; // TODO: check this! OK. I think it is fine.

				// Only do once if do_skip_ctf_correction
				if ( (do_skip_ctf_correction) && (nr_subtomo == 1) )
				{
					fn_subtomo_star = "Particles/" + fn_tomo.withoutExtension() + "_ctf.star";
					fn_subtomo_mrc  = "Particles/" + fn_tomo.withoutExtension() + "_ctf.mrc";
#ifdef DEBUG
					std::cout << " fn_subtomo_star = " << fn_subtomo_star << " , fn_subtomo_mrc = " << fn_subtomo_mrc << std::endl;
#endif
				}
				// For every sub-tomogram
				if (!do_skip_ctf_correction)
				{
					fn_subtomo_star = "Particles/" + fn_tomo.withoutExtension() + "_ctf" + integerToString(nr_subtomo, 6, '0') + ".star";
					fn_subtomo_mrc	= "Particles/" + fn_tomo.withoutExtension() + "_ctf" + integerToString(nr_subtomo, 6, '0') + ".mrc";
#ifdef DEBUG
					std::cout << " fn_subtomo_star = " << fn_subtomo_star << " , fn_subtomo_mrc = " << fn_subtomo_mrc << std::endl;
#endif
				}
				MD_this_subtomo.clear();
				MD_this_subtomo.addLabel(EMDL_CTF_DEFOCUSU);
				MD_this_subtomo.addLabel(EMDL_CTF_VOLTAGE);
				MD_this_subtomo.addLabel(EMDL_CTF_CS);
				MD_this_subtomo.addLabel(EMDL_CTF_Q0);
				MD_this_subtomo.addLabel(EMDL_ORIENT_ROT);
				MD_this_subtomo.addLabel(EMDL_ORIENT_TILT);
				MD_this_subtomo.addLabel(EMDL_ORIENT_PSI);
				MD_this_subtomo.addLabel(EMDL_CTF_BFACTOR);
				MD_this_subtomo.addLabel(EMDL_CTF_SCALEFACTOR);

				// Find minimum and maximum tilts
				RFLOAT min_tilt = 999999., max_tilt = -999999.;
				for (int id = 0; id < tilts.size(); id++)
				{
					if (tilts[id] < min_tilt)
						min_tilt = tilts[id];
					if (tilts[id] > max_tilt)
						max_tilt = tilts[id];
				}

				RFLOAT defoci = 0., tilt_deg = 0., tilt_rad = 0., xxtomo = 0., zztomo = 0.;
				RFLOAT xximg = 0., deltaD = 0., ptcldefocus = 0., tilt_scale = 0.;
				RFLOAT tilt_step = 0., tilt_diff = 0., best_tilt_diff = 0., dose_w = 0.;
				RFLOAT cur_accu_dose = 0.;
				for (int ida = 0; ida < tilts.size(); ida++)
				{
					defoci = avg_defoci[ida];
					tilt_deg = tilts[ida];
					tilt_rad = DEG2RAD(tilt_deg);
					xxtomo = float(xx - (xdim / 2) ) * calc_angpix;
					zztomo = float(zz - (zdim / 2) ) * calc_angpix;

					// Calculating the height difference of the particle from the tilt axis
					xximg = xxtomo * cos(tilt_rad) + zztomo * sin(tilt_rad);
					deltaD = xximg * sin(tilt_rad);
					ptcldefocus = defoci + deltaD;
					if (do_skip_ctf_correction)
						ptcldefocus = 0.; // TODO: Should be 0.000. I think it is fine.

					// Now weighting the 3D CTF model using the tilt dependent scale factor and the dose dependent B-Factor
					tilt_scale = cos(fabs(tilt_rad));
					if (tilts.size() <= 1)
						REPORT_ERROR("Less than 2 tilt angles are found in std::vector<RFLOAT> tilts. Tilt angles are read from " + fn_tilt_txt); // This is checked in the initialisation step.
					tilt_step = (max_tilt - min_tilt) / (RFLOAT(tilts.size()) - 1.); // Ensure that the denominator is always >= 1.0
					best_tilt_diff = tilt_step + 0.5;

					for (int idb = 0; idb < order_tilts.size(); idb++)
					{
						tilt_diff = fabs(tilt_deg - order_tilts[idb]);
						if (tilt_diff < (tilt_step + 0.25))
						{
							if (tilt_diff < best_tilt_diff)
							{
								best_tilt_diff = tilt_diff;
								cur_accu_dose = accu_dose[idb]; // TODO: cur_accu_dose always reset? Copied from the python script. Not a good way of C/C++ coding. But I think it is fine.
							}
						}
					}
					dose_w = cur_accu_dose * bfactor; // TODO: check this bfactor. T think it is fine.

					MD_this_subtomo.addObject();
					MD_this_subtomo.setValue(EMDL_CTF_DEFOCUSU, ptcldefocus);
					MD_this_subtomo.setValue(EMDL_CTF_VOLTAGE, Voltage);
					MD_this_subtomo.setValue(EMDL_CTF_CS, Cs);
					MD_this_subtomo.setValue(EMDL_CTF_Q0, AmplitudeConstrast);
					MD_this_subtomo.setValue(EMDL_ORIENT_ROT, 0.);
					MD_this_subtomo.setValue(EMDL_ORIENT_TILT, tilt_deg);
					MD_this_subtomo.setValue(EMDL_ORIENT_PSI, 0.);
					MD_this_subtomo.setValue(EMDL_CTF_BFACTOR, dose_w);
					MD_this_subtomo.setValue(EMDL_CTF_SCALEFACTOR, tilt_scale);
				}
				MD_this_subtomo.write(fn_subtomo_star);

				// Write a new line to fn_ctf_recon
				command = "";
				if ( (continue_old) && (exists(fn_subtomo_mrc)) ) // If the reconstructed CTF .mrc file exists, comment that line
					command += "# ";
				command += "relion_reconstruct --i " + fn_subtomo_star + " --o " + fn_subtomo_mrc + " --reconstruct_ctf $1 --angpix " + floatToString(calc_angpix);
				fout1 << command << std::endl;

				// Write a new object to MD_part
				FileName fn_extract_part = "Extract/" + fn_extract_job_alias + "/" + fn_tomo.withoutExtension() + integerToString(nr_subtomo, 6, '0') + ".mrc";
				if (fn_coords.getExtension() == "star")
					MD_part.addObject(MD_coords.getObject()); // Append extra information from MD_coords
				else
					MD_part.addObject(); // Otherwise, add an empty object
				MD_part.setValue(EMDL_MICROGRAPH_NAME, fn_tomo);
				MD_part.setValue(EMDL_IMAGE_COORD_X, xx);
				MD_part.setValue(EMDL_IMAGE_COORD_Y, yy);
				MD_part.setValue(EMDL_IMAGE_COORD_Z, zz);
				MD_part.setValue(EMDL_IMAGE_NAME, fn_extract_part);
				MD_part.setValue(EMDL_CTF_IMAGE, fn_subtomo_mrc);
				MD_part.setValue(EMDL_CTF_MAGNIFICATION, Magnification);
				MD_part.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, PixelSize);
			}
			// Close fn_rec here? And chmod u+x? I don't need this.
		}

		MD_part.write("particles_subtomo.star");
		fout1.close();

		command = "chmod u+x " + fn_ctf_recon;
#ifdef DEBUG
		std::cout << " " << command << std::endl;
#endif
		res = system(command.c_str());

		// Check whether all files are closed in time. I have checked. Fine.

		if (do_use_gctf) // Delete Gctf temporary files
		{
			if (exists("micrographs_all_gctf.star"))
			{
				command = "rm -rf micrographs_all_gctf.star";
				res = system(command.c_str());
			}
		}
		std::cout << std::endl;
		std::cout << " All done!" << std::endl;
		std::cout << std::endl;
		std::cout << " Please extract sub-tomograms using the RELION GUI. Remember to use the same subtomoname '--o_extract' as you gave in this script: " << fn_extract_job_alias << std::endl;
		//std::cout << " Please run the 3D CTF model volume reconstructions using the .sh scripts written in the working directory." << std::endl;
		std::cout << " Check, (split into multiple files for parallelisation) and run the script(s) from the command line:" << std::endl;
		std::cout << "  --> do_all_reconstruct_ctfs.sh SubtomogramSize" << std::endl;
		std::cout << "  (If you want to rescale the particle boxes in extraction, for example, from 200 to 100, then use 200 when running the script and rewindow the reconstructed CTF .mrc files with 'relion_image_handler --new_box 100')" << std::endl;
		std::cout << " STAR file to use for refinement (after sub-tomogram extraction and 3D CTF volume reconstruction):  particles_subtomo.star" << std::endl;
		std::cout << std::endl;

		return;
	};
};

int main(int argc, char *argv[])
{

//	time_config();

	prepare_subtomo prm;

	try
	{
		prm.read(argc, argv);
		prm.showMessages();
		if (prm.show_usage)
			return RELION_EXIT_SUCCESS;
		if (!prm.dont_check_input_files)
			prm.initialChecks();
		prm.run();
	}
	catch (RelionError XE)
	{
		//prm.usage();
		std::cout << XE;
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

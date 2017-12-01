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
#include <src/assembly.h>
#include <src/mask.h>
#include <src/fftw.h>
#include <src/image.h>
#include <src/euler.h>
#include <src/healpix_sampling.h>
#include <src/transformations.h>

#define RELION_HELIX_TOOLBOX

#define CART_TO_HELICAL_COORDS true
#define HELICAL_TO_CART_COORDS false

#ifdef RELION_HELIX_TOOLBOX
class helix_bilder_parameters
{
public:
	IOParser parser;

	// Available options
	// PLEASE MAKE SURE THAT ALL THESE OPTIONS ARE INITIALISED IN THE PARSING STEP!
	// ----------------------------------------
	bool show_usage_for_an_option;

	bool do_extract_coords_relion;
	bool do_extract_coords_ximdisp;
	bool do_convert_coords_xim2rln;
	bool do_extract_coords_eman;
	bool do_convert_coords_emn2rln;
	bool do_combine_GCTF_results;
	bool do_apply_spherical_mask_3D;
	bool do_crop_central_Z;
	bool do_create_cylinder_3D;
	bool do_set_default_tilt;
	bool do_remove_segments_with_bad_tilt;
	bool do_remove_segments_with_bad_psi;
	bool do_remove_mics_with_bad_ctf;
	bool do_simulate_helix_3D;
	bool do_impose_helical_symmetry;
	bool do_local_search_helical_symmetry;
	bool do_PDB_helix;
	bool do_divide_star_file;
	bool do_merge_star_files;
	bool do_sort_datastar_tubeID;
	bool do_simulate_helical_segments_2D;
	bool do_cut_out;
	bool do_set_xmipp_origin;
	bool do_impose_helical_symmetry_fourier_space;
	bool do_check_parameters;
	bool do_normalise_segments;
	bool do_interpo_3D_curve;
	bool do_select_3D_subtomo_from_2D_proj;
	bool do_debug;
	// ----------------------------------------

	// Input files
	FileName fn_in, fn_in1, fn_in2;

	// Rootnames of input files
	FileName fn_in_root, fn_in1_root, fn_in2_root;

	// Output files
	FileName fn_out;

	// Rootnames of output files
	FileName fn_out_root;

	// Dimensions
	int Xdim, Ydim, boxdim;

	// Number of helical subunits
	int nr_subunits;

	// Number of helical asymmetrical units
	int nr_asu;

	// Rotational symmetry - Cn
	int sym_Cn;

	// Number of filaments in a helix with seam (>= 2)
	int nr_filaments_helix_with_seam;

	// Helical rise and its local searches
	RFLOAT rise_A, rise_min_A, rise_max_A, rise_inistep_A;

	// Helical twist and its local searches
	RFLOAT twist_deg, twist_min_deg, twist_max_deg, twist_inistep_deg;

	// Pixel size in Angstroms
	RFLOAT pixel_size_A;

	// Width of soft edge
	RFLOAT width_edge_pix;

	// % of box size as the 2D / 3D spherical mask
	RFLOAT sphere_percentage;

	// % of Zdim as the central Z mask for helices
	RFLOAT z_percentage;

	// Inner and outer diameters of Z cylindrical mask
	RFLOAT cyl_inner_diameter_A, cyl_outer_diameter_A;

	// Remove segments of bad tilt angles - Maximum deviation of tilt angles allowed (away from 90 degrees)
	RFLOAT tilt_max_dev_deg;

	// Remove segments of bad psi angles - Maximum deviation of psi angles allowed (away from psi prior)
	RFLOAT psi_max_dev_deg;

	// Translate all atoms in the original PDB file to the center of mass of the molecule?
	bool do_center_of_mass_each_PDB_molecule;

	// Divide one into multiple STAR files - Number of output files
	int nr_outfiles;

	// Simulate helical segments with a STAR file - Number of helical tubes
	int nr_tubes;

	// Simulate helical subtomograms with a STAR file ?
	bool is_3d_tomo;

	// Simulate helical segments / subtomograms with a STAR file - sigma tilt, psi and offset
	RFLOAT sigma_tilt, sigma_psi, sigma_offset;

	// Simulate helical segments / subtomograms - Standard deviation of added white Gaussian noise
	RFLOAT white_noise;

	// Diameter of helical subunits (in Angstroms)
	RFLOAT subunit_diameter_A;

	// Minimum threshold of CTF FOM value, lowest resolution of EPA, minimum and maximum of defocus values
	RFLOAT ctf_fom_min, EPA_lowest_res, df_min, df_max;

	// Do bimoidal searches of tilt and psi angles in 3D helical reconstruction?
	bool do_bimodal_searches;

	// Cut helical tubes into segments?
	bool do_cut_into_segments;

	// Ignore helical symmetry in 3D reconstruction?
	bool ignore_helical_symmetry;

	// Perform local searches of helical symmetry in 3D reconstruction?
	bool do_helical_symmetry_local_refinement;

	// Construct a 3D reference for helical reconstruction with polarity along Z axis?
	bool do_polar_reference;

	// Top-bottom width ratio for construction of polarised helical reference
	RFLOAT topbottom_ratio;

	// Cut out a small part of the helix within this angle (in degrees)
	RFLOAT ang;

	// Binning factor used in manual segment picking
	int binning_factor;

	// Random seed
	int random_seed;

	// Verbosity?
	bool verb;

	void initBoolOptions()
	{
		show_usage_for_an_option = false;

		do_extract_coords_relion = false;
		do_extract_coords_ximdisp = false;
		do_convert_coords_xim2rln = false;
		do_extract_coords_eman = false;
		do_convert_coords_emn2rln = false;
		do_combine_GCTF_results = false;
		do_apply_spherical_mask_3D = false;
		do_crop_central_Z = false;
		do_create_cylinder_3D = false;
		do_set_default_tilt = false;
		do_remove_segments_with_bad_tilt = false;
		do_remove_segments_with_bad_psi = false;
		do_remove_mics_with_bad_ctf = false;
		do_simulate_helix_3D = false;
		do_impose_helical_symmetry = false;
		do_local_search_helical_symmetry = false;
		do_PDB_helix = false;
		do_divide_star_file = false;
		do_merge_star_files = false;
		do_sort_datastar_tubeID = false;
		do_simulate_helical_segments_2D = false;
		do_cut_out = false;
		do_set_xmipp_origin = false;
		do_impose_helical_symmetry_fourier_space = false;
		do_check_parameters = false;
		do_normalise_segments = false;
		do_interpo_3D_curve = false;
		do_select_3D_subtomo_from_2D_proj = false;
		do_debug = false;
	};

	helix_bilder_parameters()
	{
		clear();
	};

	~helix_bilder_parameters()
	{
		clear();
	};

	void usage()
	{
		parser.writeUsage(std::cerr);
	};

	void writeCommand(FileName fn_cmd)
	{
		std::ofstream ofs;
		ofs.open(fn_cmd.c_str(), std::ofstream::out | std::ofstream::app);

		time_t now = time(0);
	    char nodename[64] = "undefined";
	    gethostname(nodename,sizeof(nodename));
	    std::string hostname(nodename);
		ofs << std::endl << " ++++ Executed the following command at host " << hostname << " on " << ctime(&now);
		ofs << "  `which relion_helix_toolbox` " << std::flush;
		parser.writeCommandLine(ofs);
		ofs.close();
	};

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int init_section = parser.addSection("Show usage");
		show_usage_for_an_option = parser.checkOption("--function_help", "Show usage for the selected function (FEB 19, 2017)");

		int options_section = parser.addSection("List of functions (alphabetically ordered)");
		do_check_parameters = parser.checkOption("--check", "Check parameters for 3D helical reconstruction in RELION");
		do_cut_out = parser.checkOption("--cut_out", "Cut out a small part of the helix");
		do_create_cylinder_3D = parser.checkOption("--cylinder", "Create a cylinder as 3D initial reference");
		do_impose_helical_symmetry = parser.checkOption("--impose", "Impose helical symmetry (in real space)");
		do_interpo_3D_curve = parser.checkOption("--interpo", "Interpolate 3D curve for 3D helical sub-tomogram extraction");
		do_normalise_segments = parser.checkOption("--norm", "Normalise 2D/3D helical segments in a STAR file");
		do_PDB_helix = parser.checkOption("--pdb_helix", "Simulate a helix from a single PDB file of protein molecule");
		do_remove_mics_with_bad_ctf = parser.checkOption("--remove_bad_ctf", "Remove micrographs with poor-quality CTF");
		do_remove_segments_with_bad_tilt = parser.checkOption("--remove_bad_tilt", "Remove helical segments with large tilt angle deviation (away from 90 degrees)");
		do_remove_segments_with_bad_psi = parser.checkOption("--remove_bad_psi", "Remove helical segments with large psi angle deviation (away from psi prior)");
		do_local_search_helical_symmetry = parser.checkOption("--search", "Local search of helical symmetry");
		do_select_3D_subtomo_from_2D_proj = parser.checkOption("--select_3dtomo", "Select 3D subtomograms given 2D projections");
		do_simulate_helix_3D = parser.checkOption("--simulate_helix", "Create a helical 3D reference of spheres");
		do_simulate_helical_segments_2D = parser.checkOption("--simulate_segments", "Simulate helical segments using a STAR file");
		do_sort_datastar_tubeID = parser.checkOption("--sort_tube_id", "Sort segments in _data.star file according to helical tube IDs");
		do_apply_spherical_mask_3D = parser.checkOption("--spherical_mask", "Apply soft spherical mask to 3D helical reference");

		int options_old_section = parser.addSection("List of functions which can be called in Relion GUI");
		do_combine_GCTF_results = parser.checkOption("--combine_gctf", "Combine Autopicker priors (tilt and psi) with Gctf local search results");
		do_crop_central_Z = parser.checkOption("--central_mask", "Crop the central part of a helix");
		do_convert_coords_emn2rln = parser.checkOption("--coords_emn2rln", "Convert EMAN2 coordinates of helical segments into RELION STAR format");
		do_convert_coords_xim2rln = parser.checkOption("--coords_xim2rln", "Convert XIMDISP coordinates of helical segments into RELION STAR format");
		do_divide_star_file = parser.checkOption("--divide", "Divide one huge STAR file into many small ones");
		do_extract_coords_eman = parser.checkOption("--extract_emn", "Extract EMAN2 coordinates of helical segments from specified straight tubes");
		do_extract_coords_relion = parser.checkOption("--extract_rln", "Extract RELION coordinates of helical segments from specified straight tubes");
		do_extract_coords_ximdisp = parser.checkOption("--extract_xim", "Extract XIMDISP coordinates of helical segments from specified straight tubes");
		do_impose_helical_symmetry_fourier_space = parser.checkOption("--impose_fourier", "Impose helical symmetry (simulate what is done in 3D reconstruction in Fourier space)");
		do_set_default_tilt = parser.checkOption("--init_tilt", "Set tilt angles to 90 degrees for all helical segments");
		do_merge_star_files = parser.checkOption("--merge", "Merge small STAR files into a huge one");
		do_set_xmipp_origin = parser.checkOption("--set_xmipp_origin", "Set Xmipp origin");
		do_debug = parser.checkOption("--debug", "(Debug only)");

		int params_section = parser.addSection("Parameters (alphabetically ordered)");
		is_3d_tomo = parser.checkOption("--3d_tomo", "Simulate 3D subtomograms using a STAR file?");
		ang = textToFloat(parser.getOption("--ang", "Cut out a small part of the helix within this angle (in degrees)", "91."));
		pixel_size_A = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "1."));
		do_bimodal_searches = parser.checkOption("--bimodal", "Do bimodal searches of tilt and psi angles in 3D helical reconstruction?");
		binning_factor = textToInteger(parser.getOption("--bin", "Binning factor used in manual segment picking", "1"));
		boxdim = textToInteger(parser.getOption("--boxdim", "Box size (in pixels)", "-1"));
		do_center_of_mass_each_PDB_molecule = parser.checkOption("--center_pdb", "Translate all atoms in the original PDB to the center of mass of this molecule?");
		ctf_fom_min = textToFloat(parser.getOption("--ctf_fom_min", "Minimum figure-of-merit - threshold used in removing micrographs with bad CTF", "-999"));
		cyl_inner_diameter_A = textToFloat(parser.getOption("--cyl_inner_diameter", "Inner diameter of the cylindrical mask (in Angstroms)", "-1"));
		cyl_outer_diameter_A = textToFloat(parser.getOption("--cyl_outer_diameter", "Outer diameter of the cylindrical mask (in Angstroms)", "-1"));
		df_min = textToFloat(parser.getOption("--df_min", "Minimum defocus (in Angstroms)", "-999999."));
		df_max = textToFloat(parser.getOption("--df_max", "Maximum defocus (in Angstroms)", "999999."));
		EPA_lowest_res = textToFloat(parser.getOption("--EPA_lowest_res", "Lowest EPA resolution (in Angstroms) - threshold used in removing micrographs with bad CTF", "999"));
		fn_in = parser.getOption("--i", "Input file", "file.in");
		fn_in1 = parser.getOption("--i1", "Input file #1", "file01.in");
		fn_in2 = parser.getOption("--i2", "Input file #2", "file02.in");
		fn_in_root = parser.getOption("--i_root", "Rootname of input files", "_rootnameIn.star");
		fn_in1_root = parser.getOption("--i1_root", "Rootname #1 of input files", "_rootnameIn01.star");
		fn_in2_root = parser.getOption("--i2_root", "Rootname #2 of input files", "_rootnameIn02.star");
		ignore_helical_symmetry = parser.checkOption("--ignore_helical_symmetry", "Ignore helical symmetry in 3D reconstruction?");
		nr_asu = textToInteger(parser.getOption("--nr_asu", "Number of helical asymmetrical units", "1"));
		nr_outfiles = textToInteger(parser.getOption("--nr_outfiles", "Number of output files", "10"));
		nr_subunits = textToInteger(parser.getOption("--nr_subunits", "Number of helical subunits", "-1"));
		nr_tubes = textToInteger(parser.getOption("--nr_tubes", "Number of helical tubes", "-1"));
		fn_out = parser.getOption("--o", "Output file", "file.out");
		fn_out_root = parser.getOption("--o_root", "Rootname of output files", "_rootnameOut.star");
		do_polar_reference = parser.checkOption("--polar", "Construct a 3D reference for helical reconstruction with polarity along Z axis?");
		psi_max_dev_deg = textToFloat(parser.getOption("--psi_max_dev", "Maximum deviation of psi angles allowed (away from psi prior)", "15."));
		random_seed = textToFloat(parser.getOption("--random_seed", "Random seed (set to system time if negative)", "-1"));
		rise_A = textToFloat(parser.getOption("--rise", "Helical rise (in Angstroms)", "-1"));
		rise_inistep_A = textToFloat(parser.getOption("--rise_inistep", "Initial step of helical rise search (in Angstroms)", "-1"));
		rise_min_A = textToFloat(parser.getOption("--rise_min", "Minimum helical rise (in Angstroms)", "-1"));
		rise_max_A = textToFloat(parser.getOption("--rise_max", "Maximum helical rise (in Angstroms)", "-1"));
		nr_filaments_helix_with_seam = textToInteger(parser.getOption("--seam_nr_filaments", "Number of filaments in a helix with seam (>= 2)", "-1"));
		do_helical_symmetry_local_refinement = parser.checkOption("--search_sym", "Perform local searches of helical symmetry in 3D reconstruction?");
		do_cut_into_segments = parser.checkOption("--segments", "Cut helical tubes into segments?");
		sigma_offset = textToFloat(parser.getOption("--sigma_offset", "Sigma of translational offsets (in pixels)", "5."));
		sigma_psi = textToFloat(parser.getOption("--sigma_psi", "Sigma of psi angles (in degrees)", "5."));
		sigma_tilt = textToFloat(parser.getOption("--sigma_tilt", "Sigma of tilt angles (in degrees)", "5."));
		sphere_percentage = textToFloat(parser.getOption("--sphere_percentage", "Diameter of spherical mask divided by the box size (0.10~0.90 or 0.01~0.99)", "0.9"));
		subunit_diameter_A = textToFloat(parser.getOption("--subunit_diameter", "Diameter of helical subunits (in Angstroms)", "-1"));
		sym_Cn = textToInteger(parser.getOption("--sym_Cn", "Rotational symmetry Cn", "1"));
		tilt_max_dev_deg = textToFloat(parser.getOption("--tilt_max_dev", "Maximum deviation of tilt angles allowed (away from +90 degrees)", "15."));
		topbottom_ratio = textToFloat(parser.getOption("--topbottom_ratio", "Top-bottom width ratio for construction of polarised helical reference", "0.5"));
		twist_deg = textToFloat(parser.getOption("--twist", "Helical twist (in degrees, + for right-handedness)", "-1"));
		twist_inistep_deg = textToFloat(parser.getOption("--twist_inistep", "Initial step of helical twist search (in degrees)", "-1"));
		twist_min_deg = textToFloat(parser.getOption("--twist_min", "Minimum helical twist (in degrees, + for right-handedness)", "-1"));
		twist_max_deg = textToFloat(parser.getOption("--twist_max", "Maximum helical twist (in degrees, + for right-handedness)", "-1"));
		verb = parser.checkOption("--verb", "Detailed screen output?");
		white_noise = textToFloat(parser.getOption("--white_noise", "Standard deviation of added white Gaussian noise", "1."));
		width_edge_pix = textToFloat(parser.getOption("--width", "Width of cosine soft edge (in pixels)", "5."));
		Xdim = textToInteger(parser.getOption("--xdim", "Dimension X (in pixels) of the micrographs", "4096"));
		Ydim = textToInteger(parser.getOption("--ydim", "Dimension Y (in pixels) of the micrographs", "4096"));
		z_percentage = textToFloat(parser.getOption("--z_percentage", "Percentage of cropped length (along Z axis, 0.1~0.9)", "0.3"));

		RFLOAT tmp_RFLOAT = 0.;
		if (rise_min_A > rise_max_A)
			SWAP(rise_min_A, rise_max_A, tmp_RFLOAT);
		if (twist_min_deg > twist_max_deg)
			SWAP(twist_min_deg, twist_max_deg, tmp_RFLOAT);

		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
	};

	void clear()
	{
		parser.clear();
		initBoolOptions();
	};

	void displayEmptyLine()
	{
		std::cout << "=========================================================================" << std::endl;
	}

	void run()
	{
		// Check options
		int valid_options = 0;

		valid_options += (do_extract_coords_relion) ? (1) : (0);
		valid_options += (do_extract_coords_ximdisp) ? (1) : (0);
		valid_options += (do_extract_coords_eman) ? (1) : (0);
		valid_options += (do_convert_coords_emn2rln) ? (1) : (0);
		valid_options += (do_convert_coords_xim2rln) ? (1) : (0);
		valid_options += (do_combine_GCTF_results) ? (1) : (0);
		valid_options += (do_apply_spherical_mask_3D) ? (1) : (0);
		valid_options += (do_crop_central_Z) ? (1) : (0);
		valid_options += (do_create_cylinder_3D) ? (1) : (0);
		valid_options += (do_set_default_tilt) ? (1) : (0);
		valid_options += (do_remove_segments_with_bad_tilt) ? (1) : (0);
		valid_options += (do_remove_segments_with_bad_psi) ? (1) : (0);
		valid_options += (do_remove_mics_with_bad_ctf) ? (1) : (0);
		valid_options += (do_simulate_helix_3D) ? (1) : (0);
		valid_options += (do_impose_helical_symmetry) ? (1) : (0);
		valid_options += (do_local_search_helical_symmetry) ? (1) : (0);
		valid_options += (do_PDB_helix) ? (1) : (0);
		valid_options += (do_divide_star_file) ? (1) : (0);
		valid_options += (do_merge_star_files) ? (1) : (0);
		valid_options += (do_sort_datastar_tubeID) ? (1) : (0);
		valid_options += (do_simulate_helical_segments_2D) ? (1) : (0);
		valid_options += (do_cut_out) ? (1) : (0);
		valid_options += (do_set_xmipp_origin) ? (1) : (0);
		valid_options += (do_impose_helical_symmetry_fourier_space) ? (1) : (0);
		valid_options += (do_check_parameters) ? (1) : (0);
		valid_options += (do_normalise_segments) ? (1) : (0);
		valid_options += (do_interpo_3D_curve) ? (1) : (0);
		valid_options += (do_select_3D_subtomo_from_2D_proj) ? (1) : (0);
		valid_options += (do_debug) ? (1) : (0);

		if (valid_options <= 0)
			REPORT_ERROR("Please specify one option!");
		if (valid_options > 1)
			REPORT_ERROR("Only one option can be specified at one time! valid_options = " + integerToString(valid_options));

		if (do_extract_coords_relion || do_extract_coords_ximdisp || do_extract_coords_eman)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Extract coordinates of helical segments from specified straight tubes" << std::endl;
				std::cout << "  USAGE (EMAN2 format)  : --extract_emn --i_root _boxes.txt  --o_root _segments.star --nr_asu 30 --rise 1.408 --angpix 1.126 --xdim 4096 --ydim 4096 --boxdim 320 --bimodal --segments" << std::endl;
				std::cout << "  USAGE (RELION format) : --extract_rln --i_root _tubes.star --o_root _segments.star --nr_asu 30 --rise 1.408 --angpix 1.126 --xdim 4096 --ydim 4096 --boxdim 320 --bimodal --segments" << std::endl;
				std::cout << "  USAGE (XIMDISP format): --extract_xim --i_root .mrc.coords --o_root _segments.star --nr_asu 30 --rise 1.408 --angpix 1.126 --xdim 4096 --ydim 4096 --boxdim 320 --bimodal --segments" << std::endl;
				displayEmptyLine();
				return;
			}

			int format_tag;
			if (do_extract_coords_relion)
				format_tag = RELION_STAR_FORMAT;
			else if (do_extract_coords_ximdisp)
				format_tag = XIMDISP_COORDS_FORMAT;
			else if (do_extract_coords_eman)
				format_tag = EMAN2_FORMAT;
			extractHelicalSegmentsFromTubes_Multiple(
					fn_in_root,
					fn_out_root,
					format_tag,
					nr_asu,
					rise_A,
					pixel_size_A,
					Xdim,
					Ydim,
					boxdim,
					do_bimodal_searches,
					do_cut_into_segments);
		}
		else if (do_convert_coords_emn2rln || do_convert_coords_xim2rln)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Convert EMAN2 / XIMDISP coordinates of helical segments into RELION STAR format" << std::endl;
				std::cout << "  USAGE (EMAN2 format)  : --coords_emn2rln --i_root _helix_ptcl_coords.txt --o_root _segments.star --xdim 4096 --ydim 4096 --boxdim 320 --bimodal" << std::endl;
				std::cout << "  USAGE (XIMDISP format): --coords_xim2rln --i_root .mrc.coords            --o_root _segments.star --xdim 4096 --ydim 4096 --boxdim 320 --bimodal" << std::endl;
				displayEmptyLine();
				return;
			}

			int format_tag;
			if (do_convert_coords_xim2rln)
				format_tag = XIMDISP_COORDS_FORMAT;
			else if (do_convert_coords_emn2rln)
				format_tag = EMAN2_FORMAT;
			convertHelicalSegmentCoordsToStarFile_Multiple(
					fn_in_root,
					fn_out_root,
					format_tag,
					Xdim,
					Ydim,
					boxdim,
					do_bimodal_searches);
		}
		else if (do_combine_GCTF_results)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Combine Autopicker priors (tilt and psi) with Gctf local search results" << std::endl;
				std::cout << "  USAGE: --combine_gctf --i1_root _autopick.star --i2_root _gctf_local.star --o_root _combined.star" << std::endl;
				displayEmptyLine();
				return;
			}

			combineParticlePriorsWithKaiLocalCTF_Multiple(
					fn_in1_root,
					fn_in2_root,
					fn_out_root);
		}
		else if (do_apply_spherical_mask_3D)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Apply soft spherical mask to 3D helical reference" << std::endl;
				std::cout << "  USAGE: --spherical_mask --i in.mrc --o out.mrc (--sphere_percentage 0.9 --width 5)" << std::endl;
				displayEmptyLine();
				return;
			}

			int box_size;
			if ( (sphere_percentage < 0.009) || (sphere_percentage > 0.991) )
				REPORT_ERROR("Diameter of spherical mask divided by the box size should be within range 0.01~0.99!");
			Image<RFLOAT> img;
			img.read(fn_in);
			img().setXmippOrigin();
			box_size = ((XSIZE(img())) < (YSIZE(img()))) ? (XSIZE(img())) : (YSIZE(img()));
			box_size = (box_size < (ZSIZE(img()))) ? box_size : (ZSIZE(img()));
			applySoftSphericalMask(
					img(),
					(RFLOAT(box_size) * sphere_percentage),
					width_edge_pix);
			img.write(fn_out);
		}
		else if (do_crop_central_Z)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Crop the central part of a helix" << std::endl;
				std::cout << "  USAGE: --central_mask --i in.mrc --o out.mrc (--z_percentage 0.3 --width 5)" << std::endl;
				displayEmptyLine();
				return;
			}

			Image<RFLOAT> img;
			img.read(fn_in);
			cutZCentralPartOfSoftMask(
					img(),
					z_percentage,
					width_edge_pix);
			img.write(fn_out);
		}
		else if (do_create_cylinder_3D)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Create a cylinder for 3D initial reference" << std::endl;
				std::cout << "  USAGE: --cylinder --o out.mrc --boxdim 300 (--cyl_inner_diameter -1) --cyl_outer_diameter 200 --angpix 1.34 (--polar --topbottom_ratio 0.5) (--sphere_percentage 0.9 --width 5)" << std::endl;
				displayEmptyLine();
				return;
			}

			if (pixel_size_A < 0.01)
				REPORT_ERROR("Pixel size should be larger than 0!");
			if (boxdim < 20)
				REPORT_ERROR("Box size should be larger than 20 pixels!");
			if ( (sphere_percentage < 0.009) || (sphere_percentage > 0.991) )
				REPORT_ERROR("Diameter of spherical mask divided by the box size should be within range 0.01~0.99!");
			Image<RFLOAT> img;
			if (!do_polar_reference)
				topbottom_ratio = 1.;
			createCylindricalReferenceWithPolarity(
					img(),
					boxdim,
					(cyl_inner_diameter_A / pixel_size_A),
					(cyl_outer_diameter_A / pixel_size_A),
					topbottom_ratio,
					width_edge_pix);
			applySoftSphericalMask(
					img(),
					(RFLOAT(boxdim) * sphere_percentage),
					width_edge_pix);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, pixel_size_A);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, pixel_size_A);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, pixel_size_A);
			img.write(fn_out);
		}
		else if (do_set_default_tilt)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Set tilt angles to 90 degrees for all helical segments" << std::endl;
				std::cout << "  USAGE: --init_tilt --i in.star --o out.star" << std::endl;
				displayEmptyLine();
				return;
			}

			setNullTiltPriorsInDataStar(
					fn_in,
					fn_out);
		}
		else if (do_remove_segments_with_bad_tilt)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Remove helical segments with large tilt angle deviation (away from 90 degrees)" << std::endl;
				std::cout << "  USAGE: --remove_bad_tilt --i in.star --o out.star --tilt_max_dev 15" << std::endl;
				displayEmptyLine();
				return;
			}

			removeBadTiltHelicalSegmentsFromDataStar(
					fn_in,
					fn_out,
					tilt_max_dev_deg);
		}
		else if (do_remove_segments_with_bad_psi)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Remove helical segments with large psi angle deviation (away from psi prior)" << std::endl;
				std::cout << "  USAGE: --remove_bad_psi --i in.star --o out.star --psi_max_dev 15" << std::endl;
				displayEmptyLine();
				return;
			}

			removeBadPsiHelicalSegmentsFromDataStar(
					fn_in,
					fn_out,
					psi_max_dev_deg);
		}
		else if (do_remove_mics_with_bad_ctf)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Remove micrographs with poor-quality CTF" << std::endl;
				std::cout << "  USAGE: --remove_bad_ctf --i in.star --o out.star (--ctf_fom_min 0.1 --EPA_lowest_res 5 --df_min 10000 --df_max 30000)" << std::endl;
				displayEmptyLine();
				return;
			}

			excludeLowCTFCCMicrographs(
					fn_in,
					fn_out,
					ctf_fom_min,
					EPA_lowest_res,
					df_min,
					df_max);
		}
		else if (do_simulate_helix_3D)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Create a helical 3D reference of spheres" << std::endl;
				std::cout << "  USAGE: --simulate_helix --o ref.mrc --subunit_diameter 30 --cyl_outer_diameter 200 --angpix 1.126 --rise 1.408 --twist 22.03 --boxdim 300 (--sym_Cn 1) (--polar --topbottom_ratio 0.5 --cyl_inner_diameter 20) (--sphere_percentage 0.9 --width 5) (--seam_nr_filaments 13)" << std::endl;
				displayEmptyLine();
				return;
			}

			if (pixel_size_A < 0.01)
				REPORT_ERROR("Pixel size should be larger than 0!");
			if (boxdim < 20)
				REPORT_ERROR("Box size should be larger than 20 pixels!");
			if ( (sphere_percentage < 0.009) || (sphere_percentage > 0.991) )
				REPORT_ERROR("Diameter of spherical mask divided by the box size should be within range 0.01~0.99!");

			Image<RFLOAT> img;
			makeHelicalReference3DWithPolarity(
					img(),
					boxdim,
					pixel_size_A,
					twist_deg,
					rise_A,
					cyl_outer_diameter_A,
					subunit_diameter_A,
					(do_polar_reference) ? (cyl_inner_diameter_A) : (subunit_diameter_A),
					(do_polar_reference) ? (topbottom_ratio) : (1.),
					sym_Cn,
					nr_filaments_helix_with_seam);
			applySoftSphericalMask(
					img(),
					(RFLOAT(boxdim) * sphere_percentage),
					width_edge_pix);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, pixel_size_A);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, pixel_size_A);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, pixel_size_A);
			img.write(fn_out);
		}
		else if (do_impose_helical_symmetry)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Impose helical symmetry (in real space)" << std::endl;
				std::cout << "  USAGE: --impose --i in.mrc --o out.mrc (--cyl_inner_diameter -1) --cyl_outer_diameter 200 --angpix 1.126 --rise 1.408 --twist 22.03 (--z_percentage 0.3 --sphere_percentage 0.9 --width 5)" << std::endl;
				displayEmptyLine();
				return;
			}

			int box_size;
			RFLOAT sphere_diameter_A;

			Image<RFLOAT> img;
			img.read(fn_in);

			box_size = ((XSIZE(img())) < (YSIZE(img()))) ? (XSIZE(img())) : (YSIZE(img()));
			box_size = (box_size < (ZSIZE(img()))) ? (box_size) : (ZSIZE(img()));
			sphere_diameter_A = pixel_size_A * sphere_percentage * RFLOAT(box_size);

			img().setXmippOrigin();
			imposeHelicalSymmetryInRealSpace(
					img(),
					pixel_size_A,
					sphere_diameter_A / 2.,
					cyl_inner_diameter_A / 2.,
					cyl_outer_diameter_A / 2.,
					z_percentage,
					rise_A,
					twist_deg,
					width_edge_pix);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, pixel_size_A);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, pixel_size_A);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, pixel_size_A);
			img.write(fn_out);
		}
		else if (do_local_search_helical_symmetry)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Local search of helical symmetry" << std::endl;
				std::cout << "  USAGE: --search --i in.mrc (--cyl_inner_diameter -1) --cyl_outer_diameter 200 --angpix 1.126 --rise_min 1.3 --rise_max 1.5 (--rise_inistep -1) --twist_min 20 --twist_max 24 (--twist_inistep -1) (--z_percentage 0.3) (--verb)" << std::endl;
				displayEmptyLine();
				return;
			}

			int box_size;
			RFLOAT sphere_diameter_A, rise_refined_A, twist_refined_deg;

			Image<RFLOAT> img;
			img.read(fn_in);

			box_size = ((XSIZE(img())) < (YSIZE(img()))) ? (XSIZE(img())) : (YSIZE(img()));
			box_size = (box_size < (ZSIZE(img()))) ? (box_size) : (ZSIZE(img()));
			sphere_diameter_A = pixel_size_A * sphere_percentage * RFLOAT(box_size);

			img().setXmippOrigin();
			localSearchHelicalSymmetry(
					img(),
					pixel_size_A,
					sphere_diameter_A / 2.,
					cyl_inner_diameter_A / 2.,
					cyl_outer_diameter_A / 2.,
					z_percentage,
					rise_min_A,
					rise_max_A,
					rise_inistep_A,
					rise_refined_A,
					twist_min_deg,
					twist_max_deg,
					twist_inistep_deg,
					twist_refined_deg,
					((verb == true) ? (&std::cout) : (NULL)) );
			std::cout << " Done! Refined helical rise = " << rise_refined_A << " Angstroms, twist = " << twist_refined_deg << " degrees." << std::endl;
		}
		else if (do_PDB_helix)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Simulate a helix from a single PDB file of protein molecule" << std::endl;
				std::cout << "  USAGE: --pdb_helix --i in.pdb --o out.pdb --cyl_outer_diameter 50 --rise 1.408 --twist 22.03 --nr_subunits 300 (--center_pdb)" << std::endl;
				displayEmptyLine();
				return;
			}

			if ( (fn_in.getExtension() != "pdb") || (fn_out.getExtension() != "pdb") )
				REPORT_ERROR("Input and output files should be in .pdb format!");
			if (cyl_outer_diameter_A < 0.) // TODO: PLEASE CHECK THIS FOR OTHER OPTIONS !
				cyl_outer_diameter_A = 0.;

			Assembly pdb_ori, pdb_helix;
			pdb_ori.readPDB(fn_in);
			makeSimpleHelixFromPDBParticle(
					pdb_ori,
					pdb_helix,
					cyl_outer_diameter_A / 2.,
					twist_deg,
					rise_A,
					nr_subunits,
					do_center_of_mass_each_PDB_molecule);
			pdb_helix.writePDB(fn_out);
		}
		else if (do_divide_star_file)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Divide one huge STAR file into many small ones" << std::endl;
				std::cout << "  USAGE: --divide --i in.star (--nr_outfiles 10)" << std::endl;
				displayEmptyLine();
				return;
			}

			divideStarFile(fn_in, nr_outfiles);
		}
		else if (do_merge_star_files)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Merge small STAR files into a huge one" << std::endl;
				std::cout << "  USAGE: --merge --i_root _subset.star" << std::endl;
				displayEmptyLine();
				return;
			}

			mergeStarFiles(fn_in_root);
		}
		else if (do_sort_datastar_tubeID)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Sort segments in _data.star file according to helical tube IDs" << std::endl;
				std::cout << "  USAGE: --sort_tube_id --i in.star --o out.star" << std::endl;
				displayEmptyLine();
				return;
			}

			MetaDataTable MD;
			MD.read(fn_in);
			sortHelicalTubeID(MD);
			MD.write(fn_out);
		}
		else if (do_simulate_helical_segments_2D)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Simulate helical segments / subtomograms using a STAR file" << std::endl;
				std::cout << "  USAGE: --simulate_segments --i 3dvol-for-projection.mrc --o segments.star --boxdim 200 --nr_subunits 5000 --nr_asu 5 --nr_tubes 20 --twist 22.03 --rise 1.408 --cyl_outer_diameter 200 --angpix 1.126 (--bimodal --3d_tomo --sigma_tilt 5 --sigma_psi 5 --sigma_offset 5 --white_noise 1 --random_seed 1400014000)" << std::endl;
				std::cout << "  BEWARE: '--boxdim' is the shrunk box size of the simulated output (2D or 3D). It should not be bigger than the box size of the input 3D volume." << std::endl;
				displayEmptyLine();
				return;
			}

			if ( (pixel_size_A < 0.001) || ((rise_A / pixel_size_A) < 0.001) )
				REPORT_ERROR("Helical rise should be larger than 0.001 pixels!");

			simulateHelicalSegments(
					is_3d_tomo,
					fn_in,
					fn_out,
					white_noise,
					boxdim,
					nr_subunits,
					nr_asu,
					nr_tubes,
					do_bimodal_searches,
					cyl_outer_diameter_A,
					pixel_size_A,
					rise_A,
					twist_deg,
					sigma_psi,
					sigma_tilt,
					sigma_offset,
					random_seed);

			std::cout << " WARNING: Please check the output STAR files before you execute the .sh script! Use '*_helical_priors.star' or '*_no_priors.star' as the input particle STAR file!" << std::endl;
		}
		else if (do_cut_out)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Cut out a small part of the helix" << std::endl;
				std::cout << "  USAGE: --cut_out --i in.mrc --o out.mrc (--boxdim 100 --z_percentage 0.3 --ang 30)" << std::endl;
				displayEmptyLine();
				return;
			}

			Image<RFLOAT> img1, img2;
			img1.clear();
			img2.clear();
			img1.read(fn_in);
			cutOutPartOfHelix(img1(), img2(), boxdim, ang, z_percentage);
			img2.write(fn_out);
		}
		else if (do_set_xmipp_origin)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Set Xmipp origin" << std::endl;
				std::cout << "  USAGE: --set_xmipp_origin --i in.mrc --o out.mrc" << std::endl;
				displayEmptyLine();
				return;
			}
			Image<RFLOAT> img;
			img.clear();
			img.read(fn_in);
			img().setXmippOrigin();
			img.write(fn_out);
		}
		else if (do_impose_helical_symmetry_fourier_space)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Impose helical symmetry (Fourier space simulation)" << std::endl;
				std::cout << "  USAGE: --impose_fourier --i in.mrc --o out.mrc --angpix 1.126 --nr_asu 5 --rise 1.408 --twist 22.03" << std::endl;
				displayEmptyLine();
				return;
			}

			if (nr_asu <= 1)
			{
				std::cout << " Number of asymmetrical units is smaller than 1. Nothing is needed to be done..." << std::endl;
				return;
			}

			MultidimArray<RFLOAT> Msum, Maux1;
			Matrix1D<RFLOAT> transZ(3);
			Image<RFLOAT> img;
			long int Xdim, Ydim, Zdim, Ndim;

			img.read(fn_in);
			img().getDimensions(Xdim, Ydim, Zdim, Ndim);
			img().setXmippOrigin();

			if ( (Xdim != Ydim) || (Ydim != Zdim) )
				REPORT_ERROR("Error in the input 3D map: DimX != DimY or DimY != DimZ");

			Msum.clear();
			Msum.initZeros(img());
			Msum.setXmippOrigin();
			int h_min = -nr_asu / 2;
			int h_max = -h_min + nr_asu % 2;
			XX(transZ) = YY(transZ) = 0.;
			for (int hh = h_min; hh < h_max; hh++)
			{
				if (hh == 0)
					Msum += img();
				else
				{
					rotate(img(), Maux1, RFLOAT(hh) * twist_deg);
					ZZ(transZ) = RFLOAT(hh) * rise_A / pixel_size_A;
					selfTranslate(Maux1, transZ, WRAP);
					Msum += Maux1;
				}
			}
			img() = Msum / RFLOAT(nr_asu);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_X, pixel_size_A);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Y, pixel_size_A);
			img.MDMainHeader.setValue(EMDL_IMAGE_SAMPLINGRATE_Z, pixel_size_A);
			img.write(fn_out);
		}
		else if (do_check_parameters)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Check parameters for 3D helical reconstruction in RELION" << std::endl;
				std::cout << "  USAGE: --check --boxdim 300 --angpix 1.126 --sphere_percentage 0.9 (--cyl_inner_diameter 20) --cyl_outer_diameter 240 (--ignore_helical_symmetry --search_sym --z_percentage 0.3 --nr_asu 20 --rise 1.408 --rise_min 1.3 --rise_max 1.5 --twist 22.03 --twist_min 21 --twist_max 23)" << std::endl;
				displayEmptyLine();
				return;
			}
			bool result = checkParametersFor3DHelicalReconstruction(
					ignore_helical_symmetry,
					do_helical_symmetry_local_refinement,
					nr_asu,
					rise_A,
					rise_min_A,
					rise_max_A,
					twist_deg,
					twist_min_deg,
					twist_max_deg,
					boxdim,
					pixel_size_A,
					z_percentage,
					sphere_percentage * boxdim * pixel_size_A,
					cyl_inner_diameter_A,
					cyl_outer_diameter_A,
					true);
			if (result)
				std::cout << " Done! All the parameters seem OK for 3D helical reconstruction in RELION." << std::endl;
		}
		else if (do_normalise_segments)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Normalise 2D/3D helical segments in a STAR file" << std::endl;
				std::cout << "  USAGE: --norm --i imgs_input.star --o_root _norm --angpix 1.126 --cyl_outer_diameter 200" << std::endl;
				displayEmptyLine();
				return;
			}
			normaliseHelicalSegments(
					fn_in,
					fn_out_root,
					cyl_outer_diameter_A,
					pixel_size_A);
		}
		else if (do_interpo_3D_curve)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Interpolate 3D curve for 3D helical sub-tomogram extraction" << std::endl;
				std::cout << "  USAGE: --interpo --i_root Btub_tomo1 --o_root _interpo --nr_asu 1 --rise 52.77 --angpix 2.18 --boxdim 200 --bin 1 (--bimodal)" << std::endl;
				displayEmptyLine();
				return;
			}
			Interpolate3DCurves(
					fn_in_root,
					fn_out_root,
					nr_asu,
					rise_A,
					pixel_size_A,
					boxdim,
					binning_factor,
					do_bimodal_searches);
		}
		else if (do_select_3D_subtomo_from_2D_proj)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " Select 3D subtomograms given 2D projections" << std::endl;
				std::cout << "  USAGE: --select_3dtomo --i1 selected_2d_proj.star --i2 particle_3d_subtomo.star --o selected_3d_subtomo.star" << std::endl;
				displayEmptyLine();
				return;
			}
			if ( (fn_in1.getExtension() != "star")
					|| (fn_in2.getExtension() != "star")
					|| (fn_out.getExtension() != "star") )
				REPORT_ERROR("Input and output files (--i1, --i2, --o) should be in .star format!");

			MetaDataTable MD_2d, MD_3d, MD_out;
			MD_2d.read(fn_in1);
			MD_3d.read(fn_in2);
			select3DsubtomoFrom2Dproj(MD_2d, MD_3d, MD_out);
			MD_out.write(fn_out);
			std::cout << " Done! " << MD_out.numberOfObjects() << " out of " << MD_3d.numberOfObjects() << " subtomograms have been selected." << std::endl;
		}
		else if (do_debug)
		{
			if (show_usage_for_an_option)
			{
				displayEmptyLine();
				std::cout << " (Debug only)" << std::endl;
				displayEmptyLine();
				return;
			}
			//MetaDataTable MD;
			//MD.read(fn_in);
			//setPsiFlipRatioInStarFile(MD);
			//MD.write(fn_out);
			grabParticleCoordinates_Multiple(fn_in, fn_out); // RECOVER THIS !
			//readFileHeader(fn_in, fn_out, 9493);

			//Image<RFLOAT> img;
			//img.read(fn_in);
			//calculateRadialAvg(img(), pixel_size_A);
		}
		else
		{
			REPORT_ERROR("Please specify an option!");
		}
		if ( (!show_usage_for_an_option) && (!do_debug) )
		{
			writeCommand("relion_helix_toolbox.log");
		}
	};
};
#endif




// Just for debugging...
#ifndef RELION_HELIX_TOOLBOX
class helix_bilder_parameters
{
public:
	IOParser parser;
	FileName fn_in, fn_in1, fn_in2, fn_in_ctf, fn_in_ori, fn_out, fn_out1, fn_out2;
	bool do_FT, do_inv_psi, do_inv_tilt;
	RFLOAT rot_in, tilt_in, psi_in, xoff_in, yoff_in, zoff_in;
	int old_sampling, new_sampling;

	void usage()
	{
		//parser.writeUsage(std::cerr);
		return;
	};

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		fn_in = parser.getOption("--i", "Input file", "");
		fn_in1 = parser.getOption("--i1", "Input file 1", "");
		fn_in2 = parser.getOption("--i2", "Input file 2", "");
		fn_out = parser.getOption("--o", "Output file", "");
		fn_out1 = parser.getOption("--o1", "Output file 1", "");
		fn_out2 = parser.getOption("--o2", "Output file 2", "");

		do_FT = parser.checkOption("--doFT", "Do translations in Fourier space instead of real space?");
		do_inv_psi = parser.checkOption("--inv_psi", "Do inv_psi?");
		do_inv_tilt = parser.checkOption("--inv_tilt", "Do inv_tilt?");

		rot_in = textToFloat(parser.getOption("--rot", "Rot angle in degrees", "0."));
		tilt_in = textToFloat(parser.getOption("--tilt", "Tilt angle in degrees", "0."));
		psi_in = textToFloat(parser.getOption("--psi", "Psi angle in degrees", "0."));

		xoff_in = textToFloat(parser.getOption("--x", "Offest X in pixels", "0."));
		yoff_in = textToFloat(parser.getOption("--y", "Offest Y in pixels", "0."));
		zoff_in = textToFloat(parser.getOption("--z", "Offest Z in pixels", "0."));

		fn_in_ctf = parser.getOption("--ictf", "Input file 1", "");
		fn_in_ori = parser.getOption("--iori", "Input file 2", "");

		return;
	};

	void clear()
	{
		//parser.clear();
		fn_in.clear();
		fn_in1.clear();
		fn_in2.clear();
		fn_out.clear();
		fn_out1.clear();
		fn_out2.clear();
		do_FT = false;
		do_inv_psi = false;
		do_inv_tilt = false;
		rot_in = tilt_in = psi_in = 0.;
		return;
	};

	void run()
	{
		MetaDataTable MD1, MD2;
		MD1.read(fn_in1);
		MD2.read(fn_in2);
		bool result = compareLabels(MD1, MD2);
		std::cout << " MetaDataTables are the same? = " << result << std::endl;
		return;

		//std::cout << RAD2DEG(atan2(23.60, 2.73)) << std::endl;
		//std::cout << RAD2DEG(atan2(-21.19, -10.73)) << std::endl;
		//std::cout << RAD2DEG(atan2(0, 0)) << std::endl;
		//return;


		Matrix2D<RFLOAT> Ah_c, Ac_h, Amult, B;
		RFLOAT rot_out, tilt_out, psi_out;

		Euler_angles2matrix(rot_in, tilt_in, psi_in, Ah_c);
		std::cout << " == Ah_c ==" << std::endl;
		std::cout << Ah_c(0, 0) << ", " << Ah_c(0, 1) << ", " << Ah_c(0, 2) << std::endl;
		std::cout << Ah_c(1, 0) << ", " << Ah_c(1, 1) << ", " << Ah_c(1, 2) << std::endl;
		std::cout << Ah_c(2, 0) << ", " << Ah_c(2, 1) << ", " << Ah_c(2, 2) << std::endl;

		//Euler_angles2matrix(rot_in, -tilt_in, -psi_in, Ac_h);
		//std::cout << " == Ac_h ==" << std::endl;
		//std::cout << Ac_h(0, 0) << ", " << Ac_h(0, 1) << ", " << Ac_h(0, 2) << std::endl;
		//std::cout << Ac_h(1, 0) << ", " << Ac_h(1, 1) << ", " << Ac_h(1, 2) << std::endl;
		//std::cout << Ac_h(2, 0) << ", " << Ac_h(2, 1) << ", " << Ac_h(2, 2) << std::endl;

		//Amult = Ah_c * Ac_h;
		//std::cout << " == Amult = h_c * c_h ==" << std::endl;
		//std::cout << Amult(0, 0) << ", " << Amult(0, 1) << ", " << Amult(0, 2) << std::endl;
		//std::cout << Amult(1, 0) << ", " << Amult(1, 1) << ", " << Amult(1, 2) << std::endl;
		//std::cout << Amult(2, 0) << ", " << Amult(2, 1) << ", " << Amult(2, 2) << std::endl;

		//Amult = Ac_h * Ah_c;
		//std::cout << " == Amult = c_h * h_c ==" << std::endl;
		//std::cout << Amult(0, 0) << ", " << Amult(0, 1) << ", " << Amult(0, 2) << std::endl;
		//std::cout << Amult(1, 0) << ", " << Amult(1, 1) << ", " << Amult(1, 2) << std::endl;
		//std::cout << Amult(2, 0) << ", " << Amult(2, 1) << ", " << Amult(2, 2) << std::endl;

		Ac_h = Ah_c.transpose();
		Euler_matrix2angles(Ac_h, rot_out, tilt_out, psi_out);
		std::cout << " out rot, tilt, psi = " << rot_out << ", " << tilt_out << ", " << psi_out << std::endl;

		std::cout << " ############################################# " << std::endl;
		std::cout << " Rot, tilt, psi in = " << rot_in << ", " << tilt_in << ", " << psi_in << std::endl;
		std::cout << " XYZ in = " << xoff_in << ", " << yoff_in << ", " << zoff_in << std::endl;
		transformCartesianAndHelicalCoords(
				xoff_in, yoff_in, zoff_in,
				xoff_in, yoff_in, zoff_in,
				rot_in, tilt_in, psi_in,
				3,
				CART_TO_HELICAL_COORDS);
		std::cout << " XYZ h = " << xoff_in << ", " << yoff_in << ", " << zoff_in << std::endl;
		std::cout << " Now clear Z" << std::endl;
		zoff_in = 0.;
		transformCartesianAndHelicalCoords(
				xoff_in, yoff_in, zoff_in,
				xoff_in, yoff_in, zoff_in,
				rot_in, tilt_in, psi_in,
				3,
				HELICAL_TO_CART_COORDS);
		std::cout << " XYZ c = " << xoff_in << ", " << yoff_in << ", " << zoff_in << std::endl;

		MetaDataTable MD_in, MD_out;
		MD_in.read(fn_in);
		transformCartesianToHelicalCoordsForStarFiles(MD_in, MD_out);
		MD_out.write(fn_out);

		/*
		Image<RFLOAT> img;
		img.read(fn_in);
		softMaskOutsideMapForHelix(
				img(),
				psi_in,
				tilt_in,
				25.,
				20.);
		img.write(fn_out);
		*/




		/*
		Image<RFLOAT> img;
		img.read(fn_in);
		img().setXmippOrigin();
		softMaskOutsideMapForHelix(img(),
				45.,
				90.,
				400. * 0.8 * 0.5,
				180. * 0.5,
				5.);
		img.write(fn_out);
		*/
		/*
		MetaDataTable MD;
		plotLatticePoints(MD,
				1, 5, 2, 4);
		MD.write(fn_out);
		*/



		/*
		FileName fns_in, fa, fb;
		std::vector<FileName> fn_in_list;
		MetaDataTable MD_out;
		std::string cl;

		fns_in = "*" + fn_in;
		fns_in.globFiles(fn_in_list);
		if (fn_in_list.size() < 1)
			REPORT_ERROR("helix_toolbox.cpp: No input files are found!");
		for (int ii = 0; ii < fn_in_list.size(); ii++)
		{
			fb = fa = fn_in_list[ii];
			fb[0] = 'f';
			cl = "cp " + fa + " " + fb;
			int res = system(cl.c_str());
		}
		*/



		//Image<RFLOAT> img;
		//img.read(fn_in);
		//amplitudeOrPhaseMap(img(), img(), PHASE_MAP);
		//img.write(fn_out);
		/*
		MetaDataTable MD;
		int nr_opposite_polarity;
		MD.read(fn_in);
		sortHelicalTubeID(MD);
		updatePriorsForHelicalReconstruction(
									MD,
									nr_opposite_polarity,
									rot_in,
									false,
									true,
									true,
									25,
									25,
									25,
									2);
		//testDataFileTransformXY(MD);
		MD.write(fn_out);
		*/



		//RFLOAT xin = 10., yin = 7., zin = 10., psi = 0., tilt = 87., rot = 0.;
		//std::cout << " x, y, z, psi = " << xin << ", " << yin << ", " << zin << ", " << psi << std::endl;
		//transformCartesianAndHelicalCoords(xin, yin, zin, xin, yin, zin, rot, tilt, psi, 2, HELICAL_TO_CART_COORDS);
		//std::cout << " x, y, z, psi = " << xin << ", " << yin << ", " << zin << ", " << psi << std::endl;
		//transformCartesianAndHelicalCoords(xin, yin, zin, xin, yin, zin, rot, tilt, psi, 2, CART_TO_HELICAL_COORDS);
		//std::cout << " x, y, z, psi = " << xin << ", " << yin << ", " << zin << ", " << psi << std::endl;



		/*
		FourierTransformer transformer;
		MultidimArray<Complex> Faux;;
		Image<RFLOAT> img;
		std::cout << " Loading..." << std::endl;
		img.read(fn_in);
		img().setXmippOrigin();
		CenterFFT(img(), true);
		std::cout << " Doing FT..." << std::endl;
		transformer.FourierTransform(img(), Faux);

        Complex fzero;
        fzero.real = fzero.imag = 0.;
        std::cout << " Modifying..." << std::endl;
	    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux)
	    {
	    	int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
	        if (idx >= XSIZE(Faux))
	        	continue;

	        // Sphere, centered 0, 0, 76
	        //RFLOAT rr1 = sqrt((ABS(kp) - 76) * (ABS(kp) - 76) + ip*ip + jp*jp);
	        //if (rr1 < 5.01)
	        //	FFTW_ELEM(Faux, kp, ip, jp) = fzero;

	        // Sphere, centered 0, 0, 46
	        //RFLOAT rr2 = sqrt((ABS(kp) - 46) * (ABS(kp) - 46) + ip*ip + jp*jp);
	        //if (rr2 < 5.01)
	        //	FFTW_ELEM(Faux, kp, ip, jp) = fzero;

	        // Disk, centered 0, 0, 76
	        //if (ABS(ABS(kp) - 76) < 5)
	        //{
	        //	RFLOAT rr1 = sqrt(ip*ip + jp*jp);
	        //	if (rr1 < 5.)
	        //		FFTW_ELEM(Faux, kp, ip, jp) = fzero;
	        //}

	        // Disk, centered 0, 0, 46
	        //if (ABS(ABS(kp) - 46) < 5)
	        //{
	        //	RFLOAT rr1 = sqrt(ip*ip + jp*jp);
	        //	if (rr1 < 5.)
	        //		FFTW_ELEM(Faux, kp, ip, jp) = fzero;
	        //}

	        if (ABS(kp) > 44)
	        {
	        	RFLOAT rr1 = sqrt(ip*ip + jp*jp);
	        	if (rr1 < 25.)
	        		FFTW_ELEM(Faux, kp, ip, jp) = fzero;
	        }

        	//RFLOAT rr1 = sqrt(ip*ip + jp*jp);
        	//if (rr1 < 5.)
        	//	FFTW_ELEM(Faux, kp, ip, jp) = fzero;

	    }
	    std::cout << " Doing IFT..." << std::endl;
		transformer.inverseFourierTransform(Faux, img());
		CenterFFT(img(), false);
		img.write(fn_out);
		*/

		/*
		char str[1000];
		int dim;
		Matrix1D<RFLOAT> in, out;
		RFLOAT rot_deg = 0., psi_deg = 0., tilt_deg = 0., x = 0., y = 0., z = 0., angpix = 1.;
		bool direction;

		std::cout << " ========== Helical2Cartesian / Cartesian2Helical coordinates ==========" << std::endl;
		std::cout << " To Helical, press h; to Cartesian, press c..." << std::flush;
		std::cin >> str;
		if( (str[0] != 'H') && (str[0] != 'h') && (str[0] != 'C') && (str[0] != 'c') )
		{
			std::cout << "Error!" << std::endl;
			return;
		}
		direction = HELICAL_TO_CART_COORDS;
		if( (str[0] == 'H') && (str[0] == 'h') )
		{
			direction = CART_TO_HELICAL_COORDS;
		}

		std::cout << " For 2D coordinates, press 2; for 3D coordinates, press 3..." << std::flush;
		std::cin >> str;
		if( (str[0] != '2') && (str[0] != '3') )
		{
			std::cout << " Error!" << std::endl;
			return;
		}
		dim = 3;
		in.clear();
		in.resize(3);
		out.clear();
		out.resize(3);
		if(str[0] == '2')
		{
			dim = 2;
			in.clear();
			in.resize(2);
			out.clear();
			out.resize(2);
		}
		std::cout << " Please input pixel size (A / pix): " << std::flush;
		std::cin >> angpix;
		if(angpix < 0.0001)
		{
			std::cout << " Error!" << std::endl;
			return;
		}

		while(1)
		{
			if(dim == 2)
			{
				std::cout << "  Please input psi angle (in degrees): " << std::flush;
				std::cin >> psi_deg;
				std::cout << "  Please input x, y coordinates (in pixels): " << std::flush;
				std::cin >> x >> y;
				XX(in) = x;
				YY(in) = y;
			}
			else
			{
				std::cout << "  Please input psi and tilt angles (in degrees): " << std::flush;
				std::cin >> psi_deg >> tilt_deg;
				std::cout << "  Please input x, y, z coordinates (in pixels): " << std::flush;
				std::cin >> x >> y >> z;
				XX(in) = x;
				YY(in) = y;
				ZZ(in) = z;
			}
			transformCartesianAndHelicalCoords(in, out, rot_deg, tilt_deg, psi_deg, direction);
			out /= angpix;
			if(dim == 2)
			{
				std::cout << "  Output coordinates (x, y) in Angstroms = "
						<< XX(out) << ", " << YY(out) << std::endl;
			}
			else
			{
				std::cout << "  Output coordinates (x, y, z) in Angstroms = "
						<< XX(out) << ", " << YY(out) << ", " << ZZ(out) << std::endl;
			}
		}
		*/


		//RFLOAT blot_val = -50.;
		//RFLOAT blot_r = 3.;
		//RFLOAT psi_deg = 33.5561;
		//img.clear();
		//img.read(fn_in);
		//psi_deg *= -1.;
		//img().setXmippOrigin();
		//A2D_ELEM(img(), 0, 0) = -100; // y, x
		//A2D_ELEM(img(), 10, 20) = -100; // y, x
		//A2D_ELEM(img(), 20, 40) = -100; // y, x
		//A2D_ELEM(img(), 30, 60) = -100; // y, x
		//A2D_ELEM(img(), 40, 80) = -100; // y, x
		//A2D_ELEM(img(), 50, 100) = -100; // y, x
		//A2D_ELEM(img(), 60, 120) = -100; // y, x
		//makeBlot(img(), 0, 0, blot_r);
		//makeBlot(img(), 5, 10, blot_r);
		//makeBlot(img(), 10, 20, blot_r);
		//makeBlot(img(), 15, 30, blot_r);
		//makeBlot(img(), 20, 40, blot_r);
		//makeBlot(img(), 25, 50, blot_r);
		//makeBlot(img(), 20 * sin(DEG2RAD(psi_deg)), 20 * cos((DEG2RAD(psi_deg))), blot_r);
		//makeBlot(img(), 40 * sin(DEG2RAD(psi_deg)), 40 * cos((DEG2RAD(psi_deg))), blot_r);
		//makeBlot(img(), 60 * sin(DEG2RAD(psi_deg)), 60 * cos((DEG2RAD(psi_deg))), blot_r);
		//makeBlot(img(), 80 * sin(DEG2RAD(psi_deg)), 80 * cos((DEG2RAD(psi_deg))), blot_r);
		//makeBlot(img(), 100 * sin(DEG2RAD(psi_deg)), 100 * cos((DEG2RAD(psi_deg))), blot_r);
		//makeBlot(img(), 120 * sin(DEG2RAD(psi_deg)), 120 * cos((DEG2RAD(psi_deg))), blot_r);
		//img.write(fn_out);

		return;
	};

	helix_bilder_parameters()
	{
		clear();
		return;
	};

	~helix_bilder_parameters()
	{
		clear();
		return;
	};

};
#endif






int main(int argc, char *argv[])
{

//	time_config();

	helix_bilder_parameters prm;

	try
    {
		prm.read(argc, argv);
		prm.run();
    }
    catch (RelionError XE)
    {
    	prm.usage();
        std::cout << XE;
        exit(1);
    }

    return 0;
}

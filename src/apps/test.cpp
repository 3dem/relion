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

#define USEFUL_SCRIPTS_FOR_HELIX

#define CART_TO_HELICAL_COORDS true
#define HELICAL_TO_CART_COORDS false

#ifdef USEFUL_SCRIPTS_FOR_HELIX
class helix_bilder_parameters
{
public:
	IOParser parser;
	int option_id;

	// ID = 1
	FileName fn_in_id1, fn_out_id1;
	int nr_asu_id1, Xdim_id1, Ydim_id1, box_size_id1;
	RFLOAT helical_rise_A_id1, pixel_size_A_id1;

	// ID = 2
	FileName fn_in1_id2, fn_in2_id2, fn_out_id2;

	// ID = 3
	FileName fn_in_id3, fn_out_id3;
	int random_seed_id3;

	// ID = 4
	FileName fn_in1_id4, fn_in2_id4, fn_out_id4;

	// ID = 5
	FileName fn_in_id5, fn_out_id5;
	RFLOAT sphere_diameter_percentage_id5, width_edge_pix_id5;

	// ID = 6
	FileName fn_in_id6, fn_out_id6;
	RFLOAT width_edge_pix_id6, z_percentage_id6;

	// ID = 7
	FileName fn_out_id7;
	int box_size_id7;
	RFLOAT pixel_size_A_id7, sphere_diameter_percentage_id7, cyl_inner_diameter_A_id7, cyl_outer_diameter_A_id7, width_edge_pix_id7;

	// ID = 8
	FileName fn_out_id8;
	int box_size_id8;
	RFLOAT particle_diameter_A_id8, cyl_diameter_A_id8, pixel_size_A_id8;

	// ID = 9
	FileName fn_in_id9, fn_out_id9;

	// ID = 10
	FileName fn_in_id10, fn_out_id10;
	RFLOAT tilt_max_dev_deg_id10;

	// ID = 11
	FileName fn_out_id11;
	RFLOAT cyl_diameter_A_id11, pixel_size_A_id11, sphere_diameter_A_id11, helical_rise_A_id11, helical_twist_deg_id11, box_size_id11;
	int sym_Cn_id11;

	// ID = 12
	FileName fn_out_id12, fn_in_id12;
	RFLOAT cyl_inner_diameter_A_id12, cyl_outer_diameter_A_id12, pixel_size_A_id12, sphere_diameter_percentage_id12, helical_rise_A_id12, helical_twist_deg_id12, z_percentage_id12, width_edge_pix_id12;

	// ID = 13
	FileName fn_in_id13;
	RFLOAT cyl_inner_diameter_A_id13, cyl_outer_diameter_A_id13, pixel_size_A_id13, sphere_diameter_percentage_id13, helical_rise_A_id13, helical_twist_deg_id13, helical_rise_dev_percentage_id13, helical_twist_dev_percentage_id13, z_percentage_id13;

	// ID = 14
	FileName fn_in_id14, fn_out_id14;
	RFLOAT cyl_diameter_A_id14, helical_rise_A_id14, helical_twist_deg_id14;
	int nr_particles_id14;
	bool do_center_id14;

	// ID = 15
	FileName fn_in_id15, fn_out_id15;
	RFLOAT helical_rise_A_id15, helical_twist_deg_id15, pixel_size_A_id15, tilt_deg_id15, psi_deg_id15;
	int nr_projections_id15;

	// ID = 16
	FileName fn_in_id16;
	int nr_outfiles;

	// ID = 17
	FileName fn_in_id17;

	// ID = 18
	FileName fn_in_id18, fn_out_id18;

	// ID = 19
	FileName fn_out_id19;
	int nr_subunits_id19, nr_asu_id19, nr_tubes_id19;
	RFLOAT twist_deg_id19, sigma_psi_id19, sigma_tilt_id19;

	void usage()
	{
		parser.writeUsage(std::cerr);
		return;
	};

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
		option_id = textToInteger(parser.getOption("--id", "Option ID", "-1"));

		int section_1 = parser.addSection("Option ID = 1: Extract coordinates of helical segments for specified straight tubes");
		fn_in_id1 = parser.getOption("--i1", "Suffix of input files containing start and end points of helical tubes in each micrograph", "_manualpick.star");
		fn_out_id1 = parser.getOption("--o1", "Suffix of output files containing coordinates of helical segments in each micrograph", "_segments.star");
		nr_asu_id1 = textToInteger(parser.getOption("--h_nr_asu1", "Number of asymmetrical units in each box", "1"));
		helical_rise_A_id1 = textToFloat(parser.getOption("--h_rise1", "Helical rise (in Angstroms)", "-1"));
		pixel_size_A_id1 = textToFloat(parser.getOption("--angpix1", "Pixel size (in Angstroms)", "1.34"));
		Xdim_id1 = textToInteger(parser.getOption("--Xdim1", "Dimension X (in pixels) of the micrographs (assume all micrographs are of the same size)", "4096"));
		Ydim_id1 = textToInteger(parser.getOption("--Ydim1", "Dimension Y (in pixels) of the micrographs (assume all micrographs are of the same size)", "4096"));
		box_size_id1 = textToInteger(parser.getOption("--box_size1", "Box size (in pixels) for extraction of helical segments", "256"));

		int section_2 = parser.addSection("Option ID = 2: Add Autopicker orientational priors to selected segments");
		fn_in1_id2 = parser.getOption("--i2_prior", "Input file from Autopicker", "segments_all_autopick.star");
		fn_in2_id2 = parser.getOption("--i2_data", "Input file with selected segments", "segments_selected.star");
		fn_out_id2 = parser.getOption("--o2", "Output file", "segments_for_Refine3D.star");

		int section_3 = parser.addSection("Option ID = 3: Divide helical segments into two halves for Refine3D");
		fn_in_id3 = parser.getOption("--i3", "Input file", "segments_all_with_priors.star");
		fn_out_id3 = parser.getOption("--o3", "Output file", "segments_all_with_priors_subsets.star");
		random_seed_id3 = textToInteger(parser.getOption("--seed3", "Random seed (no randomisation if set to 0, system time as seed if set to negative values)", "-1"));

		int section_4 = parser.addSection("Option ID = 4: Combine Autopicker tilt and psi priors with Gctf results");
		fn_in1_id4 = parser.getOption("--i4_prior", "Suffix of input files containing tilt and psi priors", "_autopick.star");
		fn_in2_id4 = parser.getOption("--i4_ctf", "Suffix of input files containing local CTF parameters from Gctf", "_local.star");
		fn_out_id4 = parser.getOption("--o4", "Suffix of output files", "_combined.star");

		int section_5 = parser.addSection("Option ID = 5: Apply soft spherical mask to 3D helical reference");
		fn_in_id5 = parser.getOption("--i5", "Input file", "ref_in.mrc");
		fn_out_id5 = parser.getOption("--o5", "Output file", "ref_out.mrc");
		sphere_diameter_percentage_id5 = textToFloat(parser.getOption("--sphere_diameter5", "Diameter of spherical mask divided by the box size", "0.8"));
		width_edge_pix_id5 = textToFloat(parser.getOption("--width5", "Width of cosine soft edge (in pixels)", "5."));

		int section_6 = parser.addSection("Option ID = 6: Crop the central part of a helix");
		fn_in_id6 = parser.getOption("--i6", "Input file", "mask_in.mrc");
		fn_out_id6 = parser.getOption("--o6", "Output file", "mask_out.mrc");
		width_edge_pix_id6 = textToFloat(parser.getOption("--width6", "Width of cosine soft edge (in pixels)", "10."));
		z_percentage_id6 = textToFloat(parser.getOption("--z_percentage6", "Cropped length (along Z axis, 0.1~0.9)", "0.5"));

		int section_7 = parser.addSection("Option ID = 7: Create a cylinder for 3D initial reference");
		fn_out_id7 = parser.getOption("--o7", "Output file", "ref.mrc");
		box_size_id7 = textToInteger(parser.getOption("--box_size7", "Box size (in pixels)", "256"));
		cyl_inner_diameter_A_id7 = textToFloat(parser.getOption("--cyl_inner_diameter7", "Inner diameter of the cylinder (in Angstroms)", "-1"));
		cyl_outer_diameter_A_id7 = textToFloat(parser.getOption("--cyl_outer_diameter7", "Outer diameter of the cylinder (in Angstroms)", "-1"));
		pixel_size_A_id7 = textToFloat(parser.getOption("--angpix7", "Pixel size (in Angstroms)", "1.34"));
		sphere_diameter_percentage_id7 = textToFloat(parser.getOption("--sphere_diameter_percentage7", "Diameter of spherical mask divided by the box size", "0.8"));
		width_edge_pix_id7 = textToFloat(parser.getOption("--width7", "Width of cosine soft edge (in pixels)", "5."));

		int section_8 = parser.addSection("Option ID = 8: Create a tube for 2D autopicking reference");
		fn_out_id8 = parser.getOption("--o8", "Output file", "ref.mrc");
		box_size_id8 = textToInteger(parser.getOption("--box_size8", "Box size (in pixels)", "256"));
		particle_diameter_A_id8 = textToFloat(parser.getOption("--particle_diameter", "Particle diameter (in Angstroms)", "-1"));
		cyl_diameter_A_id8 = textToFloat(parser.getOption("--cyl_diameter", "Cylindrical diameter (in Angstroms)", "-1"));
		pixel_size_A_id8 = textToFloat(parser.getOption("--angpix8", "Pixel size (in Angstroms)", "1.34"));

		int section_9 = parser.addSection("Option ID = 9: Set tilt= 90 degrees for all segments");
		fn_in_id9 = parser.getOption("--i9", "Input file", "in.star");
		fn_out_id9 = parser.getOption("--o9", "Output file", "out.star");

		int section_10 = parser.addSection("Option ID = 10: Remove helical segments with large tilt angle deviation from _data.star file");
		fn_in_id10 = parser.getOption("--i10", "Input file", "in.star");
		fn_out_id10 = parser.getOption("--o10", "Output file", "out.star");
		tilt_max_dev_deg_id10 = textToFloat(parser.getOption("--tilt_max_dev10", "Maximum deviation of tilt angles allowed", "15."));

		int section_11 = parser.addSection("Option ID = 11: Create a 3D reference of spheres");
		fn_out_id11 = parser.getOption("--o11", "Output file", "ref.mrc");
		cyl_diameter_A_id11 = textToFloat(parser.getOption("--cyl_diameter11", "Cylindrical diameter (in Angstroms)", "-1"));
		pixel_size_A_id11 = textToFloat(parser.getOption("--angpix11", "Pixel size (in Angstroms)", "1.34"));
		sphere_diameter_A_id11 = textToFloat(parser.getOption("--sphere_diameter11", "Sphere diameter (in Angstroms)", "-1"));
		helical_rise_A_id11 = textToFloat(parser.getOption("--h_rise11", "Helical rise (in Angstroms)", "-1"));
		helical_twist_deg_id11 = textToFloat(parser.getOption("--h_twist11", "Helical twist (in degrees, + for right-handedness)", "-1"));
		box_size_id11 = textToInteger(parser.getOption("--box_size11", "Box size (in pixels)", "256"));
		sym_Cn_id11 = textToInteger(parser.getOption("--Cn11", "Cn symmetry?", "1"));

		int section_12 = parser.addSection("Option ID = 12: Impose helical symmetry");
		fn_in_id12 = parser.getOption("--i12", "Input file", "in.mrc");
		fn_out_id12 = parser.getOption("--o12", "Output file", "out.mrc");
		cyl_inner_diameter_A_id12 = textToFloat(parser.getOption("--cyl_inner_diameter12", "Inner diameter of the cylinder (in Angstroms)", "-1"));
		cyl_outer_diameter_A_id12 = textToFloat(parser.getOption("--cyl_outer_diameter12", "Outer diameter of the cylinder (in Angstroms)", "-1"));
		pixel_size_A_id12 = textToFloat(parser.getOption("--angpix12", "Pixel size (in Angstroms)", "1.34"));
		sphere_diameter_percentage_id12 = textToFloat(parser.getOption("--sphere_diameter_percentage12", "Diameter of spherical mask divided by the box size", "0.8"));
		helical_rise_A_id12 = textToFloat(parser.getOption("--h_rise12", "Helical rise (in Angstroms)", "-1"));
		helical_twist_deg_id12 = textToFloat(parser.getOption("--h_twist12", "Helical twist (in degrees, + for right-handedness)", "-1"));
		z_percentage_id12 = textToFloat(parser.getOption("--z_percentage12", "Cropped length (along Z axis, 0.1~0.9)", "0.3"));
		width_edge_pix_id12 = textToFloat(parser.getOption("--width12", "Width of cosine soft edge (in pixels)", "5."));

		int section_13 = parser.addSection("Option ID = 13: Local search of helical symmetry");
		fn_in_id13 = parser.getOption("--i13", "Input file", "in.mrc");
		cyl_inner_diameter_A_id13 = textToFloat(parser.getOption("--cyl_inner_diameter13", "Inner diameter of the cylinder (in Angstroms)", "-1"));
		cyl_outer_diameter_A_id13 = textToFloat(parser.getOption("--cyl_outer_diameter13", "Outer diameter of the cylinder (in Angstroms)", "-1"));
		pixel_size_A_id13 = textToFloat(parser.getOption("--angpix13", "Pixel size (in Angstroms)", "1.34"));
		sphere_diameter_percentage_id13 = textToFloat(parser.getOption("--sphere_diameter_percentage13", "Diameter of spherical mask divided by the box size", "0.8"));
		helical_rise_A_id13 = textToFloat(parser.getOption("--h_rise13", "Helical rise (in Angstroms)", "-1"));
		helical_twist_deg_id13 = textToFloat(parser.getOption("--h_twist13", "Helical twist (in degrees, + for right-handedness)", "-1"));
		helical_rise_dev_percentage_id13 = textToFloat(parser.getOption("--h_rise_dev13", "Maximum deviation of helical rise (0.0~0.3)", "0.1"));
		helical_twist_dev_percentage_id13 = textToFloat(parser.getOption("--h_twist_dev13", "Maximum deviation of helical twist (0.0~0.3)", "0.1"));
		z_percentage_id13 = textToFloat(parser.getOption("--z_percentage13", "Cropped length (along Z axis, 0.1~0.9)", "0.3"));

		int section_14 = parser.addSection("Option ID = 14: Make helix from a PDB file");
		fn_in_id14 = parser.getOption("--i14", "Input file", "in.pdb");
		fn_out_id14 = parser.getOption("--o14", "Output file", "out.pdb");
		cyl_diameter_A_id14 = textToFloat(parser.getOption("--cyl_diameter14", "Cylindrical diameter (in Angstroms)", "0."));
		helical_rise_A_id14 = textToFloat(parser.getOption("--h_rise14", "Helical rise (in Angstroms)", "-1"));
		helical_twist_deg_id14 = textToFloat(parser.getOption("--h_twist14", "Helical twist (in degrees, + for right-handedness)", "-1"));
		nr_particles_id14 = textToInteger(parser.getOption("--n14", "Number of particles", "100"));
		do_center_id14 = parser.checkOption("--center14", "Translate all atoms in the original PDB file to the center of mass of the molecule");

		int section_15 = parser.addSection("Option ID = 15: Make helical reconstruction STAR file from a single 2D class average");
		fn_in_id15 = parser.getOption("--i15", "Input file of single 2D class average", "in.mrc");
		fn_out_id15 = parser.getOption("--o15", "Output star file", "out.star");
		helical_rise_A_id15 = textToFloat(parser.getOption("--h_rise15", "Helical rise (in Angstroms)", "-1"));
		helical_twist_deg_id15 = textToFloat(parser.getOption("--h_twist15", "Helical twist (in degrees, + for right-handedness)", "-1"));
		pixel_size_A_id15 = textToFloat(parser.getOption("--angpix15", "Pixel size (in Angstroms)", "1.34"));
		nr_projections_id15 = textToInteger(parser.getOption("--n15", "Number of projections", "100"));
		tilt_deg_id15 = textToFloat(parser.getOption("--tilt15", "Tilt angle", "90."));
		psi_deg_id15 = textToFloat(parser.getOption("--psi15", "Psi angle", "0."));

		int section_16 = parser.addSection("Option ID = 16: Divide one into multiple STAR files");
		fn_in_id16 = parser.getOption("--i16", "Input star file", "in.star");
		nr_outfiles = textToInteger(parser.getOption("--n16", "Number of output files", "10"));

		int section_17 = parser.addSection("Option ID = 17: Combine multiple STAR files");
		fn_in_id17 = parser.getOption("--i17", "Rootname of input star files", "root");

		int section_18 = parser.addSection("Option ID = 18: Sort _data.star file according to rlnHelicalTubeID");
		fn_in_id18 = parser.getOption("--i18", "Input star file", "in.star");
		fn_out_id18 = parser.getOption("--o18", "Output star file", "out.star");

		int section_19 = parser.addSection("Option ID = 19: Simulate helical segments with a STAR file");
		fn_out_id19 = parser.getOption("--o19", "Output star file", "out.star");
		nr_subunits_id19 = textToInteger(parser.getOption("--h_n19", "Number of subunits", "1000"));
		nr_asu_id19 = textToInteger(parser.getOption("--h_nr_asu19", "Number of asymmetrical units", "1"));
		nr_tubes_id19 = textToInteger(parser.getOption("--h_nr_tube19", "Number of helical tubes", "10"));
		twist_deg_id19 = textToFloat(parser.getOption("--h_twist19", "Helical rise (in Angstroms)", "-1."));
		sigma_psi_id19 = textToFloat(parser.getOption("--sigma_psi19", "Helical rise (in Angstroms)", "1."));
		sigma_tilt_id19 = textToFloat(parser.getOption("--sigma_tilt19", "Helical rise (in Angstroms)", "1."));

		return;
	};

	void clear()
	{
		parser.clear();
		return;
	};

	void run()
	{
		if (option_id == 1)
		{
			extractCoordsForAllHelicalSegments_Multiple(
					fn_in_id1,
					fn_out_id1,
					nr_asu_id1,
					helical_rise_A_id1,
					pixel_size_A_id1,
					Xdim_id1,
					Ydim_id1,
					box_size_id1);
		}
		else if (option_id == 2)
		{
			addPriorsToParticleDataFile(fn_in1_id2, fn_in2_id2, fn_out_id2);
		}
		else if (option_id == 3)
		{
			divideHelicalSegmentsFromMultipleMicrographsIntoRandomHalves(fn_in_id3, fn_out_id3, random_seed_id3);
		}
		else if (option_id == 4)
		{
			combineParticlePriorsWithKaiLocalCTF_Multiple(fn_in1_id4, fn_in2_id4, fn_out_id4);
		}
		else if (option_id == 5)
		{
			int box_size;
			if ( (sphere_diameter_percentage_id5 < 0.1) || (sphere_diameter_percentage_id5 > 0.9) )
				REPORT_ERROR("Diameter of spherical mask divided by the box size should be within range 0.1~0.9!");
			Image<RFLOAT> img;
			img.read(fn_in_id5);
			img().setXmippOrigin();
			box_size = ((XSIZE(img())) < (YSIZE(img()))) ? (XSIZE(img())) : (YSIZE(img()));
			box_size = (box_size < (ZSIZE(img()))) ? box_size : (ZSIZE(img()));
			applySoftSphericalMask(img(), (RFLOAT(box_size) * sphere_diameter_percentage_id5), width_edge_pix_id5);
			img.write(fn_out_id5);
		}
		else if (option_id == 6)
		{
			Image<RFLOAT> img;
			img.read(fn_in_id6);
			cutZCentralPartOfSoftMask(img(), z_percentage_id6, width_edge_pix_id6);
			img.write(fn_out_id6);
		}
		else if (option_id == 7)
		{
			if (pixel_size_A_id7 < 0.01)
				REPORT_ERROR("Pixel size should be larger than 0!");
			if (box_size_id7 < 20)
				REPORT_ERROR("Box size should be larger than 20 pixels!");
			if ( (sphere_diameter_percentage_id7 < 0.1) || (sphere_diameter_percentage_id7 > 0.9) )
				REPORT_ERROR("Diameter of spherical mask divided by the box size should be within range 0.1~0.9!");
			Image<RFLOAT> img;
			createCylindricalReference(img(), box_size_id7, (cyl_inner_diameter_A_id7 / pixel_size_A_id7), (cyl_outer_diameter_A_id7 / pixel_size_A_id7), width_edge_pix_id7);
			applySoftSphericalMask(img(), (RFLOAT(box_size_id7) * sphere_diameter_percentage_id7), width_edge_pix_id7);
			img.write(fn_out_id7);
		}
		else if (option_id == 8)
		{
			Image<RFLOAT> img;
			makeHelicalReference2D(img(), box_size_id8, particle_diameter_A_id8, cyl_diameter_A_id8, pixel_size_A_id8);
			img.write(fn_out_id8);
		}
		else if (option_id == 9)
		{
			bool rewrite_tilt = true;
			bool rewrite_psi = false;
			setNullAlignmentPriorsInDataStar(fn_in_id9, fn_out_id9, rewrite_tilt, rewrite_psi);
		}
		else if (option_id == 10)
		{
			removeBadTiltHelicalSegmentsFromDataStar(fn_in_id10, fn_out_id10, tilt_max_dev_deg_id10);
		}
		else if (option_id == 11)
		{
			Image<RFLOAT> img;
			makeHelicalReference3D(
					img(),
					box_size_id11,
					pixel_size_A_id11,
					helical_twist_deg_id11,
					helical_rise_A_id11,
					cyl_diameter_A_id11,
					sphere_diameter_A_id11,
					sym_Cn_id11);
			img.write(fn_out_id11);
		}
		else if (option_id == 12)
		{
			int box_size_id12;
			RFLOAT sphere_diameter_A_id12;

			Image<RFLOAT> img;
			img.read(fn_in_id12);

			box_size_id12 = ((XSIZE(img())) < (YSIZE(img()))) ? (XSIZE(img())) : (YSIZE(img()));
			box_size_id12 = (box_size_id12 < (ZSIZE(img()))) ? (box_size_id12) : (ZSIZE(img()));
			sphere_diameter_A_id12 = pixel_size_A_id12 * sphere_diameter_percentage_id12 * RFLOAT(box_size_id12);

			img().setXmippOrigin();
			makeHelicalReferenceInRealSpace(
					img(),
					pixel_size_A_id12,
					helical_twist_deg_id12,
					helical_rise_A_id12,
					z_percentage_id12,
					sphere_diameter_A_id12 / 2.,
					cyl_inner_diameter_A_id12 / 2.,
					cyl_outer_diameter_A_id12 / 2.,
					width_edge_pix_id12);
			img.write(fn_out_id12);
		}
		else if (option_id == 13)
		{
			int box_size_id13;
			RFLOAT sphere_diameter_A_id13, rise_refined_A, twist_refined_deg;

			Image<RFLOAT> img;
			img.read(fn_in_id13);

			box_size_id13 = ((XSIZE(img())) < (YSIZE(img()))) ? (XSIZE(img())) : (YSIZE(img()));
			box_size_id13 = (box_size_id13 < (ZSIZE(img()))) ? (box_size_id13) : (ZSIZE(img()));
			sphere_diameter_A_id13 = pixel_size_A_id13 * sphere_diameter_percentage_id13 * RFLOAT(box_size_id13);

			img().setXmippOrigin();
			localSearchHelicalSymmetry(
					img(),
					pixel_size_A_id13,
					sphere_diameter_A_id13 / 2.,
					cyl_inner_diameter_A_id13 / 2.,
					cyl_outer_diameter_A_id13 / 2.,
					z_percentage_id13,
					helical_rise_A_id13,
					helical_twist_deg_id13,
					helical_rise_dev_percentage_id13,
					helical_twist_dev_percentage_id13,
					rise_refined_A,
					twist_refined_deg);

			std::cout << " Refined helical rise = " << rise_refined_A
					<< " Angstroms, refined helical twist = " << twist_refined_deg << " degrees." << std::endl;
		}
		else if (option_id == 14)
		{
			if ( (fn_in_id14.getExtension() != "pdb") || (fn_out_id14.getExtension() != "pdb") )
				REPORT_ERROR("relion_test::option_id == 14: Input and output files should be in pdb format!");

			Assembly pdb_ori, pdb_helix;
			pdb_ori.readPDB(fn_in_id14);
			makeSimpleHelixFromPDBParticle(pdb_ori, pdb_helix, cyl_diameter_A_id14 / 2., helical_twist_deg_id14, helical_rise_A_id14, nr_particles_id14, do_center_id14);
			pdb_helix.writePDB(fn_out_id14);
		}
		else if (option_id == 15)
		{
			makeHelicalReconstructionStarFileFromSingle2DClassAverage(
					fn_in_id15,
					fn_out_id15,
					pixel_size_A_id15,
					helical_twist_deg_id15,
					helical_rise_A_id15,
					tilt_deg_id15,
					psi_deg_id15,
					nr_projections_id15);
		}
		else if (option_id == 16)
		{
			divideStarFile(fn_in_id16, nr_outfiles);
		}
		else if (option_id == 17)
		{
			combineStarFiles(fn_in_id17);
		}
		else if (option_id == 18)
		{
			MetaDataTable MD;
			MD.read(fn_in_id18);
			sortHelicalTubeID(MD);
			MD.write(fn_out_id18);
		}
		else if (option_id == 19)
		{
			simulateHelicalSegments(
					fn_out_id19,
					nr_subunits_id19,
					nr_asu_id19,
					nr_tubes_id19,
					twist_deg_id19,
					sigma_psi_id19,
					sigma_tilt_id19);
		}
		else
		{
			REPORT_ERROR("test.cpp::run(): Invalid option ID!");
		}

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





#ifndef USEFUL_SCRIPTS_FOR_HELIX
class helix_bilder_parameters
{
public:
	IOParser parser;
	FileName fn_in, fn_in1, fn_in2, fn_in_ctf, fn_in_ori, fn_out, fn_out1, fn_out2;
	bool do_FT, do_inv_psi, do_inv_tilt;
	RFLOAT rot_in, tilt_in, psi_in;
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
		// STAR files of 3D helical map projections
		int nr_tubes, nr_segments, nr_subunits, nr_asu, tube_id;
		RFLOAT twist_deg, rot, psi, tilt;
		MetaDataTable MD;

		nr_tubes = 30;
		nr_asu = 2;
		twist_deg = 164.97;
		nr_subunits = 5000;

		// TODO: raise error if nr_asu<0 or too big, n too small!
		if ( (nr_tubes < 2) || (nr_subunits < 100) || (nr_asu < 1) || (((nr_subunits / nr_asu) / nr_tubes) < 5) )
			return;

		nr_segments = nr_subunits / nr_asu;

		MD.clear();
	    MD.addLabel(EMDL_ORIENT_ROT);
	    MD.addLabel(EMDL_ORIENT_TILT);
	    MD.addLabel(EMDL_ORIENT_PSI);
	    MD.addLabel(EMDL_ORIENT_ORIGIN_X);
	    MD.addLabel(EMDL_ORIENT_ORIGIN_Y);
	    MD.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);

	    tube_id = 0;
		for (int id = 0; id < nr_segments; id++)
		{
			if ( ( (id % (nr_segments / nr_tubes)) == 0 )
					&& ( (nr_segments - id) >= (nr_segments / nr_tubes) ) )
			{
				tube_id++;
				tilt = rnd_unif(85., 95.);
				rot = rnd_unif(0.01, 359.99);
				psi = rnd_unif(-179.99, 179.99);
			}

			rot += twist_deg * ((RFLOAT)(nr_asu));
			rot = realWRAP(rot, -180., 180.);

			MD.addObject();
	    	MD.setValue(EMDL_ORIENT_ROT, rot);
	    	MD.setValue(EMDL_ORIENT_TILT, realWRAP(tilt + rnd_gaus(0., 1.), 0., 180.) );
	    	MD.setValue(EMDL_ORIENT_PSI, realWRAP(psi + rnd_gaus(0., 1.), -180., 180.) );
	    	MD.setValue(EMDL_ORIENT_ORIGIN_X, 0.);
	    	MD.setValue(EMDL_ORIENT_ORIGIN_Y, 0.);
	    	MD.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, tube_id);
		}
		MD.write(fn_out);


		//Image<RFLOAT> img;
		//RFLOAT rise_A, twist_deg, rise_refined_A, twist_refined_deg;
		//img.clear();

		/*
		extractCoordsForAllHelicalSegments_Multiple(
				fn_in,
				fn_out,
				90,
				1.406,
				0.562,
				7420,
				7420,
				320);
		*/


		//combineParticlePriorsWithKaiLocalCTF_Multiple(fn_in1, fn_in2, fn_out);

		//combineParticlePriorsWithClass2DDataStar(fn_in1, fn_in2, fn_out);

		//divideHelicalSegmentsFromMultipleMicrographsIntoTwoHalves(fn_in, fn_out);


		/*
		createCylindricalReference(
				img(),
				224,
				-1,
				30);
		*/
		//img.read(fn_in);
		//img().setXmippOrigin();
		//applySoftSphericalMask(img(), (1350. / 5.36));
		//img.write(fn_out);
		//std::cout << "Best psi = " << searchPsiFor2DHelicalSegment(img(), 2.68, 1350. / 2., 880. / 2.) << " degrees." << std::endl;
		//std::cout << "Best psi = " << searchPsiFor2DHelicalSegment(img(), 1.126, 300. / 2., 250. / 2.) << " degrees." << std::endl;
		//selfRotate(img(), rot_in, 'Z', DONT_WRAP, DIRECT_A2D_ELEM(img(), 0, 0));
		//img.write(fn_out);

		//setNullAlignmentPriorsInDataStar(fn_in, fn_out, true, true);


		/*
		Image<RFLOAT> img;
		img.read(fn_in);
		img().setXmippOrigin();
		RFLOAT twist_deg = 30;  // 23.4
		RFLOAT rise_A = 1.9;  // 53.1
		RFLOAT max_dev = 0.3;
		RFLOAT pixel_size_A = 1.126;
		RFLOAT z_percentage = 0.2;
		RFLOAT rise_refined_A;
		RFLOAT twist_refined_deg;

		localSearchHelicalSymmetry_NEW(
				img(),
				pixel_size_A,
				340. / 2.,
				-1.,
				250. / 2.,
				z_percentage,
				rise_A,
				twist_deg,
				max_dev,
				max_dev,
				rise_refined_A,
				twist_refined_deg);
		std::cout << "NEW: " << twist_refined_deg << ", " << rise_refined_A << std::endl;
		std::cout << "##########################################################" << std::endl;

		localSearchHelicalSymmetry(
				img(),
				pixel_size_A,
				340. / 2.,
				-1.,
				250. / 2.,
				z_percentage,
				rise_A,
				twist_deg,
				max_dev,
				rise_refined_A,
				twist_refined_deg);
		std::cout << "OLD: " << twist_refined_deg << ", " << rise_refined_A << std::endl;
		std::cout << "##########################################################" << std::endl;
		*/

		//applySoftSphericalMask(img(), 220);
		//img.write(fn_out);
		/*
		fn_in = "_manualpick.star";
		fn_out = "_asu10.star";
		extractCoordsForAllHelicalSegments_Multiple(
				fn_in,
				fn_out,
				10,
				1.408,
				1.126,
				3600,
				3600,
				320);
		*/

		//fn_in = ".coords";
		//fn_out = "_p.star";
		//transformXimdispHelicalCoordsToStarFile_Multiple(fn_in, fn_out, 4096., 4096., 200.);

		//std::string rootname_priors = "_p.star";
		//std::string rootname_local_ctf = "_local.star";
		//combineParticlePriorsWithKaiLocalCTF_Multiple(rootname_priors, rootname_local_ctf);

		//transformXimdispHelicalCoordsToStarFile(fn_in, fn_out, 4096., 4096., 200.);
		//combineParticlePriorsWithKaiLocalCTF(fn_in1, fn_in2, fn_out);
		//setNullAlignmentPriorsInDataStar(fn_in, fn_out, true, false);
		/*
		HealpixSampling sampling;
		std::vector<int> pointer_dir_nonzeroprior, pointer_psi_nonzeroprior;
		std::vector<RFLOAT> directions_prior, psi_prior;
		pointer_dir_nonzeroprior.clear();
		pointer_psi_nonzeroprior.clear();
		directions_prior.clear();
		psi_prior.clear();

		sampling.clear();
		sampling.is_3D = true;
		sampling.is_3d_trans = false;
		sampling.healpix_order = 5; // 1.8 degrees
		sampling.fn_sym = "c1";
		sampling.limit_tilt = -91;
		sampling.psi_step = 1.875;
		sampling.offset_range = 5;
		sampling.offset_step = 2;
		sampling.perturbation_factor = 0.5;
		sampling.initialise(1, 3, false, true, 1.2, 22.03, -1);


		sampling.rot_angles.clear();
		sampling.tilt_angles.clear();
		sampling.psi_angles.clear();

		sampling.rot_angles.push_back(-179.);
		sampling.rot_angles.push_back(156.);
		sampling.rot_angles.push_back(27.);
		sampling.rot_angles.push_back(-133.);
		sampling.rot_angles.push_back(-144.);
		sampling.rot_angles.push_back(-59.);

		sampling.tilt_angles.push_back(77.);
		sampling.tilt_angles.push_back(80.);
		sampling.tilt_angles.push_back(83.);
		sampling.tilt_angles.push_back(97.);
		sampling.tilt_angles.push_back(100.);
		sampling.tilt_angles.push_back(103.);

		sampling.psi_angles.push_back(29.);
		sampling.psi_angles.push_back(31.);
		sampling.psi_angles.push_back(-151.);
		sampling.psi_angles.push_back(-149.);
		sampling.psi_angles.push_back(209.);
		sampling.psi_angles.push_back(211.);
		sampling.psi_angles.push_back(389.);
		sampling.psi_angles.push_back(391.);
		sampling.psi_angles.push_back(-331.);
		sampling.psi_angles.push_back(-329.);

		sampling.selectOrientationsWithNonZeroPriorProbability(77., 80., 30.,
				0., 1.01, 1.01,
				pointer_dir_nonzeroprior, directions_prior, pointer_psi_nonzeroprior, psi_prior,
				true, true, 3.);
		*/


		//setNullAlignmentPriorsInDataStar(fn_in, fn_out, true, false);

		//extractCoordsForAllHelicalSegments(fn_in, fn_out, 30, 5.80818, 1.33);
		//removeBadTiltHelicalSegmentsFromDataStar(fn_in, fn_out, 15.);
		//combineParticlePriorsWithKaiLocalCTF(fn_in1, fn_in2, fn_out);
		//addNullAlignmentPriorsToDataStar(fn_in, fn_out);

		//Image<RFLOAT> img;
		//img.read(fn_in);
		//enlarge3DReference(img(), 2);
		//applySoftSphericalMask(img());
		//img.write(fn_out);

		////Image<RFLOAT> img;
		//img.read(fn_in);
		//createCylindricalReference(img(), 320, -1., 160.);
		//applySoftSphericalMask(img());
		//cutZCentralPartOfSoftMask(img(), 0.4, 15);
		//img.write(fn_out);

		/*
		MetaDataTable MDlocalctf, MD;
		MDlocalctf.read(fn_in_ctf);
		MD.read(fn_in_ori);
		RFLOAT x, y, v, du, dv, da, cs, ps, fom, mag;
		int ii;
		std::vector<RFLOAT> lv, ldu, ldv, lda, lcs, lps, lfom, lmag;
		lv.clear(); ldu.clear(); ldv.clear(); lda.clear(); lcs.clear(); lps.clear(); lfom.clear(); lmag.clear();

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDlocalctf)
		{
			MDlocalctf.getValue(EMDL_CTF_VOLTAGE, v);
			MDlocalctf.getValue(EMDL_CTF_DEFOCUSU, du);
			MDlocalctf.getValue(EMDL_CTF_DEFOCUSV, dv);
			MDlocalctf.getValue(EMDL_CTF_DEFOCUS_ANGLE, da);
			MDlocalctf.getValue(EMDL_CTF_CS, cs);
			MDlocalctf.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, ps);
			MDlocalctf.getValue(EMDL_CTF_FOM, fom);
			MDlocalctf.getValue(EMDL_CTF_MAGNIFICATION, mag);

			lv.push_back(v);
			ldu.push_back(du);
			ldv.push_back(dv);
			lda.push_back(da);
			lcs.push_back(cs);
			lps.push_back(ps);
			lfom.push_back(fom);
			lmag.push_back(mag);
		}
		ii = -1;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			ii++;
			MD.setValue(EMDL_CTF_VOLTAGE, lv[ii]);
			MD.setValue(EMDL_CTF_DEFOCUSU, ldu[ii]);
			MD.setValue(EMDL_CTF_DEFOCUSV, ldv[ii]);
			MD.setValue(EMDL_CTF_DEFOCUS_ANGLE, lda[ii]);
			MD.setValue(EMDL_CTF_CS, lcs[ii]);
			MD.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, lps[ii]);
			MD.setValue(EMDL_CTF_FOM, lfom[ii]);
			MD.setValue(EMDL_CTF_MAGNIFICATION, lmag[ii]);
		}
		MD.write("out.star");
		*/



		/*
		Image<RFLOAT> img, Isolvent;
		// psi_deg # 46 60 65 83 = -90.7041, 33.5561, 17.2951, -138.72
		RFLOAT psi_deg = -27.986498;
		RFLOAT tilt_deg = 4325324.23598324;
		RFLOAT helical_mask_radius_pix = 55.;
		//RFLOAT helical_mask_len_pix = 30;
		RFLOAT cosine_width = 5.;
		RFLOAT shift_x = 30.;
		RFLOAT shift_y = 60.;
		RFLOAT shift_z = 0.;
		int dim;


		TabSine sin0, sin1;
		TabCosine cos0, cos1;
		std::cout << "==============================" << std::endl;
		std::cout << "Fill in old tabfuncs ..." << std::endl;
		sin0.initialise(5000);
		cos0.initialise(5000);
		std::cout << "Old done !!!" << std::endl;
		std::cout << "Fill in new tabfuncs ..." << std::endl;
		sin1.initialise(1000000);
		cos1.initialise(1000000);
		std::cout << "New done !!!" << std::endl;

		std::ofstream fout1, fout2;
		fout1.open(fn_in.c_str(), std::ios::out);
		fout2.open(fn_out.c_str(), std::ios::out);
		int samplings, multiples;
		samplings = 5000;
		multiples = 99999;
		for (int ii = 0; ii < samplings; ii++)
		{
			RFLOAT val;
			val = ii * (2. * PI / samplings);
			val += multiples * 2. * PI;
			fout1 << ii << "	" << sin0(val) << std::endl;
			fout2 << ii << "	" << sin1(val) << std::endl;
		}
		fout1.close();
		fout2.close();
		*/



		/*
		Matrix1D<RFLOAT> dir, dir1, dir2;
		RFLOAT rot0, tilt0, psi0, rot1, tilt1, psi1, rot2, tilt2, psi2, diff_ang;
		rot0 = tilt0 = psi0 = rot1 = tilt1 = psi1 = rot2 = tilt2 = psi2 = 0.;

		bool out = false;
		for (rot0 = -148.82; rot0 < +165.343; rot0 += 0.135)
		{
			out = false;
			for (tilt0 = 7.36; tilt0 < +171.42; tilt0 += 0.183)
			{
				Euler_angles2direction(rot0, tilt0, dir);
				Euler_direction2angles(dir, rot1, tilt1);
				Euler_direction2angles_New(dir, rot2, tilt2);
				diff_ang = fabs(fabs(rot1 - rot2) + fabs(tilt1 - tilt2));
				if (diff_ang > 0.01)
				{
					out = true;
					std::cout << " #    In: " << rot0 << ", " << tilt0 << std::flush;
					std::cout << " # OutOld: " << rot1 << ", " << tilt1 << std::flush;
					std::cout << " # OutNew: " << rot2 << ", " << tilt2 << std::endl;
				}
			}
			if (out == true)
			{
				std::cout << "----------------------------------------------------------------------" << std::endl;
			}
		}
		*/

		/*
		MetaDataTable MDt;
		int ret;
		MDt.clear();
		MDt.read(fn_in);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDt)
		{
			RFLOAT psi, tilt;

			if (do_inv_psi)
			{
				MDt.getValue(EMDL_ORIENT_PSI, psi);
				if (fabs(psi) > 0.01)
				{
					if (psi > 0.)
						psi -= 180.;
					else
						psi += 180.;
				}
				MDt.setValue(EMDL_ORIENT_PSI, psi);
			}

			MDt.getValue(EMDL_ORIENT_TILT, tilt);
			if (do_inv_tilt)
			{
				if ( (tilt < 45. ) || (tilt > 135.) )
				{
					std::cout << "Tilt angle is not within suitable range!" << std::endl;
				}
				tilt = 180. - tilt;
			}
			else
			{
				tilt = 90.;
			}
			MDt.setValue(EMDL_ORIENT_TILT, tilt);
		}
		MDt.write(fn_out);




		MultidimArray<Complex> Faux1, Faux2;
		FourierTransformer transformer;
		Matrix1D<RFLOAT> trans;

		if (do_FT)
		{
			std::cout << "doFT" << std::endl;
			img.clear();
			Faux1.clear();
			Faux2.clear();
			transformer.clear();
			img.read(fn_in);

			dim = img().getDim();

			CenterFFT(img(), true);
			transformer.FourierTransform(img(), Faux1);
			if (dim == 2)
				shiftImageInFourierTransform(Faux1, Faux2, (RFLOAT)(XSIZE(img())), shift_x, shift_y);
			else
				shiftImageInFourierTransform(Faux1, Faux2, (RFLOAT)(XSIZE(img())), shift_x, shift_y, shift_z);
			transformer.inverseFourierTransform(Faux2, img());
			CenterFFT(img(), false);

			img.write(fn_out);
			img.clear();
			Faux1.clear();
			Faux2.clear();
			transformer.clear();
		}
		else
		{
			std::cout << "doRealSpace" << std::endl;
			trans.clear();
			img.clear();
			img.read(fn_in);

			dim = img().getDim();
			trans.resize(2);
			XX(trans) = shift_x;
			YY(trans) = shift_y;
			if (dim == 3)
			{
				trans.resize(3);
				ZZ(trans) = shift_z;
			}

			trans.selfROUND();
			selfTranslate(img(), trans, DONT_WRAP);

			img.write(fn_out);
			trans.clear();
			img.clear();
		}




		img.read(fn_in);
		softMaskOutsideMapForHelix(img(), psi_deg, tilt_deg, helical_mask_radius_pix, cosine_width);
		img.write(fn_out);




		img.read(fn_in);
		Isolvent.clear();
		Isolvent().resize(img());
		img().setXmippOrigin();
		Isolvent().setXmippOrigin();
		RFLOAT sphere, sphere_p, hradius, hradius_p, r, d, val1, val2;
		int boxsize;
		boxsize = (XSIZE(Isolvent()) < YSIZE(Isolvent())) ? (XSIZE(Isolvent())) : (YSIZE(Isolvent()));
		boxsize = (boxsize < ZSIZE(Isolvent())) ? (boxsize) : (ZSIZE(Isolvent()));
		boxsize = boxsize / 2 - ((boxsize + 1) % 2);
		sphere_p = (RFLOAT)(boxsize);
		sphere = sphere_p - cosine_width;
		hradius = helical_mask_radius_pix;
		hradius_p = hradius + cosine_width;
		FOR_ALL_ELEMENTS_IN_ARRAY3D(Isolvent())
		{
			d = (RFLOAT)(i * i + j * j);
			r = d + (RFLOAT)(k * k);
			d = sqrt(d);
			r = sqrt(r);
			if ( (r < sphere) && (d < hradius) )
				A3D_ELEM(Isolvent(), k, i, j) = 1.;
			else if ( (r > sphere_p) || (d > hradius_p) )
				A3D_ELEM(Isolvent(), k, i, j) = 0.;
			else
			{
				val1 = 0.5 - 0.5 * cos(PI * (sphere_p - r) / cosine_width );
				val2 = 0.5 - 0.5 * cos(PI * (hradius_p - d) / cosine_width );
				A3D_ELEM(Isolvent(), k, i, j) = (val1 < val2) ? val1 : val2;
			}
		}
		img() = img() * Isolvent();
		img.write(fn_out);




		Matrix1D<RFLOAT> c(2), h(2);
		RFLOAT psi_deg = 30;
		RFLOAT tilt_deg = 90.;





		XX(c) = 10; // x
		YY(c) = 0; // y
		std::cout << "psi_deg , tilt_deg = " << psi_deg << ", " << tilt_deg << std::endl;
		transformCartesianAndHelicalCoords(c, h, psi_deg, tilt_deg, CART_TO_HELICAL_COORDS);
		std::cout << "CART to HELICAL" << std::endl;
		std::cout << " CART = " << XX(c) << ", " << YY(c) << std::endl;
		std::cout << " HELICAL = " << XX(h) << ", " << YY(h) << std::endl;
		std::cout << "ROTATE BACK TO CART" << std::endl;
		transformCartesianAndHelicalCoords(h, h, psi_deg, tilt_deg, HELICAL_TO_CART_COORDS);
		std::cout << " CART = " << XX(h) << ", " << YY(h) << std::endl;



		XX(h) = 10; // r
		YY(h) = 5.9; // p
		std::cout << "psi_deg , tilt_deg = " << psi_deg << ", " << tilt_deg << std::endl;
		transformCartesianAndHelicalCoords(h, c, psi_deg, tilt_deg, HELICAL_TO_CART_COORDS);
		std::cout << "HELICAL to CART" << std::endl;
		std::cout << " HELICAL = " << XX(h) << ", " << YY(h) << std::endl;
		std::cout << " CART = " << XX(c) << ", " << YY(c) << std::endl;
		std::cout << "ROTATE BACK TO HELICAL" << std::endl;
		transformCartesianAndHelicalCoords(c, c, psi_deg, tilt_deg, CART_TO_HELICAL_COORDS);
		std::cout << " HELICAL = " << XX(c) << ", " << YY(c) << std::endl;



		char str[1000];
		int dim;
		Matrix1D<RFLOAT> in, out;
		RFLOAT psi_deg, tilt_deg, x, y, z, angpix;
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
			transformCartesianAndHelicalCoords(in, out, psi_deg, tilt_deg, direction);
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



		//int ii;
		//std::ofstream fout;
		//Image<RFLOAT> img, img_out;
		//std::vector<RFLOAT> cn_list, cc_list, rise_pix_list, twist_deg_list;
		//std::vector<int> nr_asym_voxels_list;

		//normalise2DImageSlices(fn_in, fn_out, 140);

		//Assembly pdb;
		//pdb.readPDB(fn_in);
		//pdb.writePDB(fn_out);

		// Make helix from PDB


		/*
		Assembly pdb_ori, pdb_helix;
		pdb_ori.readPDB(fn_in);
		//makeSimpleHelixFromPDBParticle(pdb_ori, pdb_helix, -1., 22.03, 1.408, 500);  // TMV
		makeSimpleHelixFromPDBParticle(pdb_ori, pdb_helix, 310., 146.975, 5.716, 100);  // AChR
		//makeSimpleHelixFromPDBParticle(pdb_ori, pdb_helix, 300, 30.0, 0.01, 4);  // (test)
		//makeSimpleHelixFromPDBParticle(pdb_ori, pdb_helix, -1., 24.35, 52.26, 20);  // MamK
		//pdb_helix = pdb_ori;
		pdb_helix.writePDB(fn_out);
		*/



		// Rotate helix


		/*
		Assembly pdb_ori;
		Matrix1D<RFLOAT> shift;
		Matrix2D<RFLOAT> rotational_matrix;
		rotational_matrix.clear();
		rotation3DMatrix(-90., 'Y', rotational_matrix, false);
		shift.resize(3);
		shift.initZeros();
		pdb_ori.readPDB(fn_in);
		pdb_ori.applyTransformation(rotational_matrix, shift);
		pdb_ori.writePDB(fn_out);
		*/


		//img.clear();
		//img.read(fn_in);
		//rotate2DZSliceInFT(img(), 23., 2);
		//shift3DVolumeAlongZAxisInFT(img(), 50.);
		//expandZaxis(img, 90., 20., -16, 6);
		//img.write(fn_out);

		//fout.open(fn_out.c_str(), std::ios::out);
		//searchHelicalSymmetry(img(), 30., 80., 2.0610, 0.0001, 30, 64.740, 0.001, 20, rise_pix_list, twist_deg_list, cc_list, nr_asym_voxels_list, &fout);
		//searchHelicalSymmetry(img(), 0., 1000., 9.4730, 0.0001, 20, 210.510, 0.001, 20, rise_pix_list, twist_deg_list, cc_list, nr_asym_voxels_list, &fout);
		//searchHelicalSymmetry(img(), 0., 1000., 10.880, 0.001, 20, 343.40, 0.01, 20, rise_pix_list, twist_deg_list, cc_list, nr_asym_voxels_list, &fout);
		//searchHelicalSymmetry(img(), 0., 1000., 5.20, 0.01, 20, 165.0, 0.1, 100, rise_pix_list, twist_deg_list, cc_list, nr_asym_voxels_list, &fout);
		//searchHelicalSymmetry(img(), 0., 1000., 15.00, 0.01, 2, 0.0, 1.0, 360, rise_pix_list, twist_deg_list, cc_list, nr_asym_voxels_list, &fout);
		//fout.close();

		//expandZaxis(img, -132.92, 1.3152, -25, 25);

		//imposeCnSymmetry(img, 7);
		//imposeHelicalSymmetry(img, -132.923, 1.3152, -20, 20);
		//img.clear();
		//img.read(fn_in);
		//expandZaxis(img, -132.923, 1.3152, -30, 30);
		//img.write(fn_out);

		//img.clear();
		//img.read(fn_in);
		//softMaskOutsideMapHelicalParticle2D(img(), 30, 37295.3684, 10., NULL);
		//img.write(fn_out);

		//img.clear();
		//img.read(fn_in);

		//makeGoodHelix(img(), 50, -132.923, 1.3152);
		//img.write(fn_out);

		//img.clear();
		//img.read(fn_in);
		//calcRadialAverage(img(), cc_list);
		//for(ii = 0; ii < cc_list.size(); ii++)
		//{
		//	std::cout << "Radius = " << ii << ", avg_pix_val = " << cc_list[ii] << std::endl;
		//}
		//fout.open(fn_out.c_str(), std::ios::out);
		//searchHelicalSymmetry(img(), 30., 80., 1.0700, 0.0001, 300, 146.940, 0.001, 100, rise_pix_list, twist_deg_list, cc_list, nr_asym_voxels_list, &fout);
		//fout.close();

		//img.clear();
		//img.read(fn_in);
		//imposeHelicalSymmetry(img(), -146.98604, 1.09176, -100, 100);
		//imposeHelicalSymmetry(img(), 146.975, 1.0745, -100, 100);
		//imposeHelicalSymmetry(img(), -14.876, 5.373, -100, 100);
		//imposeHelicalSymmetry(img(), 14.876, 5.373, -100, 100);
		//img.write(fn_out);

		//img.clear();
		//img.read(fn_in);
		//makeGoodHelix(img(), 0.2, 80, -146.975, 1.0745, 2);
		//makeGoodHelix(img(), 0.4, 8000, 22.03, (1.408 / 1.124), 2);
		//img.write(fn_out);

		//img.clear();
		//createHelicalCentralSoftMask(img(), 200, 200, 200, 0.2, -0.0001, 60., 3.);
		//createCylinderRef(img(), 200, 200, 200, -0.0001, 60., 3.);
		//img.write(fn_out);

		//img.clear();
		//img.read(fn_in);
		//imposeHelicalSymmetryInRealSpace(img(), 160., -146.98604, (5.80818 / 2.66), -100, 100);
		//makeGoodHelixInRealSpace(img(), 0.3, 104., -22.03, (1.408 / 1.124));
		//img.write(fn_out);

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




/*
 * helix_bilder.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: she
 */



/*
#include <src/args.h>
#include <src/fftw.h>
#include <src/metadata_table.h>
#include <src/matrix1d.h>
#include <src/matrix2d.h>
#include <src/transformations.h>
#include <src/time.h>
#include "src/image.h"
#include "src/macros.h"
#include <iostream>

#define INPUT_FILE_FIRST_FEW_LINES 6
#define CHAR_BUFFER_SIZE 1000
#define VERY_SMALL_RFLOAT (1e-15)

// START AND END COORDINATES MUST BE OF RFLOAT PAIRS!
// EACH OF THE PAIRS MUST BE OF 2 VALUES!

class all_coordinates_picker_straight_line
{
public:
	IOParser parser;
	FileName fn_in, fn_out, fn_out_cylinder, fn_out_mask;
	RFLOAT stepSizeInPixels;
	int cyl_radius, mask_radius, cyl_3Dbox, add_prior_pos;

	RFLOAT k, b;

	RFLOAT x1, x2, y1, y2, delta_x, delta_y;

	void clear()
	{
		parser.clear();
		fn_in.clear();
		fn_out.clear();
		fn_out_cylinder.clear();
		x1 = x2 = y1 = y2 = delta_x = delta_y = 0.;
		cyl_radius = cyl_3Dbox = add_prior_pos = 0;
		return;
	}

	all_coordinates_picker_straight_line()
	{
		clear();
		return;
	}

	~all_coordinates_picker_straight_line()
	{
		clear();
		return;
	}

	void usage()
	{
		parser.writeUsage(std::cerr);
		return;
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		fn_in = parser.getOption("--i", "Input .star file (start and end coordinates)",
				"/fhgfs2/she/GTMV/Micrographs_bin/CS003_TMV_01_bin_manualpick.star");
	    fn_out = parser.getOption("--o", "Output .star file (all picked coordinates)",
	    		"/fhgfs2/she/GTMV/Micrographs_bin/CS003_TMV_01_bin_manualpickall.star");

	    fn_out_cylinder = parser.getOption("--cyl", "Output cylinder .mrc file",
	    		"/fhgfs2/she/AChR175_my/ref_cylinder.mrc");
	    fn_out_mask = parser.getOption("--cyl", "Output cylinder .mrc file",
	    		"/fhgfs2/she/AChR175_my/mask_cylinder.mrc");

	    // TMV: asu = 30, rise = 1.408A, angpix = 1.126A, step = 30*1.408/1.126
	    stepSizeInPixels = textToFloat(parser.getOption("--step", "Step size of segmentation in pixels", "37.513321492"));

	    cyl_radius = textToInteger(parser.getOption("--rcyl", "Reference cylinder - radius - in pixels", "80"));
	    mask_radius = textToInteger(parser.getOption("--maskcyl", "Reference cylindrical mask - radius - in pixels", "100"));
	    cyl_3Dbox = textToInteger(parser.getOption("--boxcyl", "Reference cylinder - 3Dbox diameter - in pixels", "320"));

	    add_prior_pos = textToInteger(parser.getOption("--psipriorpos", "Number of first additional column of _rlnXXX in new .star files", "3"));

	    return;
	}

	void output_cylinder_reference()
	{
		Image<RFLOAT> cyl_img;
		int _r;

		cyl_img.clear();
		cyl_img().resize(cyl_3Dbox, cyl_3Dbox, cyl_3Dbox);
		cyl_img().setXmippOrigin();

		FOR_ALL_ELEMENTS_IN_ARRAY3D(cyl_img())
		{
			A3D_ELEM(cyl_img(), k, i, j) = 0.;

			_r = ROUND(sqrt(i * i + j * j));
			if(_r <= cyl_radius)
			{
				A3D_ELEM(cyl_img(), k, i, j) = 1.;
			}
		}

		cyl_img.write(fn_out_cylinder);

		return;
	}

	void output_cylinder_mask()
	{
		Image<RFLOAT> cyl_img;
		int _r;

		cyl_img.clear();
		cyl_img().resize(cyl_3Dbox, cyl_3Dbox, cyl_3Dbox);
		cyl_img().setXmippOrigin();

		FOR_ALL_ELEMENTS_IN_ARRAY3D(cyl_img())
		{
			A3D_ELEM(cyl_img(), k, i, j) = 0.;

			_r = ROUND(sqrt(i * i + j * j));
			if(_r <= mask_radius)
			{
				A3D_ELEM(cyl_img(), k, i, j) = 1.;
			}
		}

		cyl_img.write(fn_out_mask);

		return;
	}

	void find_all_coordinates()
	{
		int ii;
		const RFLOAT pi = 3.141592653589793238462643383279502884197;
		char charBuf[CHAR_BUFFER_SIZE];
		RFLOAT now_x, now_y;

		ii = 0;
		now_x = now_y = 0.;
		charBuf[0] = '\0';

		// Open a .star input file
		std::ifstream infile(fn_in.data(), std::ios_base::in);
	    if(infile.fail())
	    {
	    	REPORT_ERROR( (std::string) "all_coordinates_picker_straight_line::find_all_coordinates: File " + fn_in.c_str() + " does not exists" );
	    }
	    FileName ext = fn_in.getFileFormat();
	    if(ext != "star")
	    {
	    	REPORT_ERROR("all_coordinates_picker_straight_line::find_all_coordinates ERROR: metadatatable should have .star extension");
	    }

		// Open a .star output file
	    std::ofstream outfile(fn_out.data(), std::ios_base::out);
	    if(infile.fail())
	    {
	    	REPORT_ERROR( (std::string) "all_coordinates_picker_straight_line::find_all_coordinates: Output to file:  " + fn_out.c_str() + " error" );
	    }

	    for(ii = 0; ii < INPUT_FILE_FIRST_FEW_LINES ;ii++)
	   {
	    	infile.getline(charBuf, CHAR_BUFFER_SIZE - 1);
	    //	outfile << charBuf << std::endl;
	    }
	    //outfile << "_rlnAngleTiltPrior #" << add_prior_pos << std::endl;
	    outfile << std::endl;
	    outfile << "data_" << std::endl;
	    outfile << std::endl;
	    outfile << "loop_" << std::endl;
	    outfile << "_rlnCoordinateX #1" << std::endl;
		outfile << "_rlnCoordinateY #2" << std::endl;
		outfile << "_rlnAngleRot #3" << std::endl;
		outfile << "_rlnAngleTilt #4" << std::endl;
		outfile << "_rlnAnglePsi #5" << std::endl;
		outfile << "_rlnOriginX #6" << std::endl;
		outfile << "_rlnOriginY #7" << std::endl;
		RFLOAT rot, tilt, psi, xoff, yoff;
		rot = 0.;
		tilt = 90.;
		psi = 0.;
		xoff = 0.;
		yoff = 0.;

	    ii = 0;
	    while(infile >> x1)
	    {
	    	ii++;

	    	infile >> y1 >> x2 >> y2;
	    	//std::cout << "x1: " << x1 << " y1: " << y1 << " x2: " << x2 << " y2: " << y2 << std::endl;
	    	if(fabs(x2 - x1) > VERY_SMALL_RFLOAT)
	    	{
	    		k = (y2 - y1) / (x2 - x1);
	    	}
	    	else
	    	{
		    	REPORT_ERROR( (std::string) "all_coordinates_picker_straight_line::find_all_coordinates: X coordinates are the same in at least one point pair" );
	    	}

	    	b = y1 - k * x1;

	    	std::cout << "Detected Line " << ii << ": y = (" << k << ") x + (" << b << ")" << std::endl;

	    	delta_x = stepSizeInPixels / (sqrt(k * k + 1.));
	    	if(x1 > x2)
	    	{
	    		delta_x *= -1.;
	    	}
	    	delta_y = k * delta_x;
	    	std::cout << "  Delta x: " << delta_x << "    Delta y: " << delta_y << std::endl;

	    	now_x = x1;
	    	now_y = y1;
	    	while(1)
	    	{
	    		now_x += delta_x;
	    		now_y += delta_y;
	    		if( ((now_x > x1) && (now_x > x2)) || ((now_x < x1) && (now_x < x2))
	    				|| ((now_y > y1) && (now_y > y2)) || ((now_y < y1) && (now_y < y2)) )
	    		{
	    			break;
	    		}
	    		else
	    		{
	    			psi = -180. * atan(k) / pi;

	    			//outfile << "  " << now_x << "  " << now_y << "  90.0000" << std::endl;
	    			outfile << now_x << "	"
	    					<< now_y << "	"
							<< rot << "	"
							<< tilt << "	"
							<< psi << "	"
							<< xoff << "	"
							<< yoff << "	"
							<< std::endl;
	    		}
	    	}

	    }

	    infile.close();
	    outfile.close();

	    return;
	}


};

int main(int argc, char *argv[])
{

//	time_config();

	all_coordinates_picker_straight_line prm;

	try
    {
		prm.read(argc, argv);
		prm.find_all_coordinates();
		//prm.output_cylinder_reference();
		//prm.output_cylinder_mask();

    }
    catch (RelionError XE)
    {
        prm.usage();
        std::cout << XE;
        exit(1);
    }

    return 0;
}
*/



/*
class test_proj
{
public:
	int nr_part, nr_asu, add_comments;
	RFLOAT rise_A, angpix, twist_deg, sigma_offset_pix, sigma_tilt_deg, dU, dUs, dUe, kV, Cs, Q0, detector_pixel_size;
	IOParser parser;
	FileName fn_root, fn_out;

	void usage()
	{
		parser.writeUsage(std::cerr);
		return;
	};

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");

		nr_part = textToInteger(parser.getOption("--n", "Number of total particles", "1500"));
		nr_asu = textToInteger(parser.getOption("--asu", "Helical asymmetrical units in each new box", "3"));
		add_comments = textToInteger(parser.getOption("--comments", "Output a STAR file with comments?", "0"));
		rise_A = textToFloat(parser.getOption("--rise", "Helical rise in Angstroms", "1.408"));
		twist_deg = textToFloat(parser.getOption("--twist", "Helical twist in degrees", "22.03"));
		angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "1.126"));
		sigma_offset_pix = textToFloat(parser.getOption("--sigma_offset", "Sigma for Gaussian offset in pixels", "0."));
		sigma_tilt_deg = textToFloat(parser.getOption("--sigma_tilt", "Sigma for Gaussian tilt in degrees", "2."));
		dUs = textToFloat(parser.getOption("--dUs", "CTF: Minimum defocus U ", "10000."));
		dUe = textToFloat(parser.getOption("--dUe", "CTF: Maximum defocus U ", "30000."));
		kV = textToFloat(parser.getOption("--kV", "CTF: Voltage (kV)", "300"));
		Cs = textToFloat(parser.getOption("--Cs", "CTF: Spherical Aberration (mm)", "2.7"));
		Q0 = textToFloat(parser.getOption("--Q0", "CTF: Amplitude contrast", "0.1"));
		detector_pixel_size = textToFloat(parser.getOption("--detector_pix", "CTF: Detector pixel size (um)", "14."));

		fn_root = parser.getOption("--root", "Output root name",
				"/lmb/home/she/Desktop/");
		fn_out = parser.getOption("--o", "Output STAR file name (without extension)",
				"r0p0tilt2asu3");

		return;
	};

	void clear()
	{
		nr_part = -1;
		nr_asu = 1;
		add_comments = 0;
		rise_A = 0.;
		twist_deg = 0.;
		angpix = 0.;
		sigma_offset_pix = 0.;
		sigma_tilt_deg = 0.;
		dU = 0.;
		dUs = 0.;
		dUe = 0.;
		kV = 0.;
		Cs = 0.;
		Q0 = 0.;
		detector_pixel_size = 0.;
		parser.clear();
		fn_root.clear();
		fn_out.clear();
		return;
	};

	void checkParameters()
	{
		if( (nr_asu < 1) || (nr_part < 10) || ((nr_part / nr_asu) < 10)
				|| (dUs > dUe) || (dUs < 1000.) || (dUe > 100000.)
				|| (fabs(twist_deg) < 0.01) || (angpix < 0.01) || (sigma_offset_pix < 0.) || (sigma_tilt_deg < 0.) || (detector_pixel_size < 0.) )
		{
			REPORT_ERROR("test_proj::checkParameters(): Wrong parameter(s)!");
			return;
		}
		return;
	};

	void run()
	{
		int ii, jj, nr_particles_in_a_group;
		RFLOAT rot_deg, tilt_deg, psi_deg, x_pix, y_pix, r_pix, p_pix, uni_rand, gaus_rand;
		bool do_inv = true;
		Matrix1D<RFLOAT> cartesian, helical;
		FileName fn_temp;

		fn_temp.clear();
		fn_temp = fn_root + fn_out;
		// Open a .star output file
		if(add_comments <= 0)
		{
			fn_temp += ".star";
		}
		else
		{
			fn_temp += "_comments.star";
		}
		fn_out.clear();
		fn_out = fn_temp;
	    std::ofstream outfile(fn_out.data(), std::ios_base::out);

	    //randomize_random_generator();

	    if(outfile.fail())
	    {
	    	REPORT_ERROR( (std::string) "all_coordinates_picker_straight_line::find_all_coordinates: Output to file:  " + fn_out.c_str() + " error" );
	    }

	    outfile << std::endl;
	    outfile << "data_" << std::endl;
	    outfile << std::endl;
	    outfile << "loop_" << std::endl;
	    outfile << "_rlnOriginX #1" << std::endl;
	    outfile << "_rlnOriginY #2" << std::endl;
	    outfile << "_rlnAngleRot #3" << std::endl;
	    outfile << "_rlnAngleTilt #4" << std::endl;
	    outfile << "_rlnAnglePsi #5" << std::endl;
	    outfile << "_rlnDefocusU #6" << std::endl;
	    outfile << "_rlnVoltage #7" << std::endl;
	    outfile << "_rlnSphericalAberration #8" << std::endl;
	    outfile << "_rlnAmplitudeContrast #9" << std::endl;
	    outfile << "_rlnMagnification #10" << std::endl;
	    outfile << "_rlnDetectorPixelSize #11" << std::endl;

	    rise_A = fabs(rise_A);
	    dU = dUs;
	    for(ii = 0; ii < (nr_part / nr_asu); ii++)
	    {
	    	dU = rnd_unif(dUs, dUe);

	    	// Generate r_pix, p_pix
	    	r_pix = rnd_unif((-0.5 * rise_A / angpix), (0.5 * rise_A / angpix));
	    	p_pix = rnd_gaus(0., sigma_offset_pix);

	    	// Generate 3 Euler angles
	    	rot_deg = ii * nr_asu * twist_deg;
	    	rot_deg = rot_deg - 360. * ROUND(rot_deg / 360.);
	    	tilt_deg = 90. + rnd_gaus(0., sigma_tilt_deg);
	    	psi_deg = rnd_unif(-180.0, 180.0);

	    	// Store helical coordinates
	    	helical.clear();
	    	helical.resize(2);
	    	XX(helical) = r_pix;
	    	YY(helical) = p_pix;

	    	// Helical -> Cartesian coordinates
	    	cartesian.clear();
	    	cartesian.resize(2);
	    	transformCartesianAndHelicalCoords(helical, cartesian, psi_deg, tilt_deg, HELICAL_TO_CART_COORDS);

	    	// Output
	    	outfile << XX(cartesian) << "    "
	    			<< YY(cartesian) << "    "
					<< rot_deg << "    "
					<< tilt_deg << "    "
					<< psi_deg << "    "
					<< dU << "    "
					<< kV << "    "
					<< Cs << "    "
					<< Q0 << "    "
					<< (detector_pixel_size * 10000. / angpix) << "    "
					<< detector_pixel_size
					<< std::flush;
	    	if(add_comments > 0)
	    	{
	    		outfile << "            # r, p = " << XX(helical) << ", " << YY(helical)
	    				<< " pix; x, y = " << XX(cartesian) << ", " << YY(cartesian)
						<< " pix; rot, tilt, psi = " << rot_deg << ", " << tilt_deg << ", " << psi_deg << " degrees." << std::flush;
	    	}
	    	outfile << std::endl;
	    }

	    outfile.close();
		return;
	};

	test_proj()
	{
		clear();
		return;
	};

	~test_proj()
	{
		clear();
		return;
	};

};

int main(int argc, char *argv[])
{

//	time_config();

	test_proj prm;

	try
    {
		prm.read(argc, argv);
		prm.checkParameters();
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
*/


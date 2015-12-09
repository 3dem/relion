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
#include "src/preprocessing.h"

void Preprocessing::read(int argc, char **argv, int rank)
{

	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_star_in = parser.getOption("--i", "The STAR file with all (selected) micrographs to extract particles from","");
	fn_coord_suffix = parser.getOption("--coord_suffix", "The suffix for the coordinate files, e.g. \"_picked.star\" or \".box\"","");
	fn_coord_dir = parser.getOption("--coord_dir", "The directory where the coordinate files are (default is same as micrographs)", "ASINPUT");
	fn_part_dir = parser.getOption("--part_dir", "Output directory for particle stacks", "Particles/");
	fn_part_star = parser.getOption("--part_star", "Output STAR file with all particles metadata", "particles.star");

	int extract_section = parser.addSection("Particle extraction");
	do_extract = parser.checkOption("--extract", "Extract all particles from the micrographs");
	do_helical_segments = parser.checkOption("--helical_segments", "Perform operations on helical segments");
	do_premultiply_ctf = parser.checkOption("--premultiply_ctf", "Premultiply the micrograph/frame with its CTF prior to particle extraction");
	do_ctf_intact_first_peak = parser.checkOption("--ctf_intact_first_peak", "When premultiplying with the CTF, leave frequencies intact until the first peak");
	do_phase_flip = parser.checkOption("--phase_flip", "Flip CTF-phases in the micrograph/frame prior to particle extraction");
	angpix  =  textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms (only necessary for phase-flipping if magnification and detector pixel size are not in STAR file)", "1."));
	extract_size = textToInteger(parser.getOption("--extract_size", "Size of the box to extract the particles in (in pixels)", "-1"));
	extract_bias_x  = textToInteger(parser.getOption("--extract_bias_x", "Bias in X-direction of picked particles (this value in pixels will be added to the coords)", "0"));
	extract_bias_y  = textToInteger(parser.getOption("--extract_bias_y", "Bias in Y-direction of picked particles (this value in pixels will be added to the coords)", "0"));
	do_movie_extract = parser.checkOption("--extract_movies", "Extract particles from movie stacks?");
	movie_name = parser.getOption("--movie_name", "Movie-identifier to extract particles from movie stacks (e.g. _movie.mrcs)", "movie");
	avg_n_frames = textToInteger(parser.getOption("--avg_movie_frames", "Average over this number of individual movie frames", "1"));
	movie_first_frame = textToInteger(parser.getOption("--first_movie_frame", "Extract from this movie frame onwards", "1"));
	movie_first_frame--; // (start counting at 0, not 1)
	movie_last_frame = textToInteger(parser.getOption("--last_movie_frame", "Extract until this movie frame (default=all movie frames)", "0"));
	movie_last_frame--; // (start counting at 0, not 1)
	fn_movie = parser.getOption("--movie_rootname", "Common name to relate movies to the single micrographs (e.g. mic001_movie.mrcs related to mic001.mrc)", "movie");

	int perpart_section = parser.addSection("Particle operations");
	do_project_3d = parser.checkOption("--project3d", "Project sub-tomograms along Z to generate 2D particles");
	scale  = textToInteger(parser.getOption("--scale", "Re-scale the particles to this size (in pixels)", "-1"));
	window  = textToInteger(parser.getOption("--window", "Re-window the particles to this size (in pixels)", "-1"));
	do_normalise = parser.checkOption("--norm", "Normalise the background to average zero and stddev one");
	do_ramp = !parser.checkOption("--no_ramp", "Just subtract the background mean in the normalisation, instead of subtracting a fitted ramping background. ");
	bg_radius = textToInteger(parser.getOption("--bg_radius", "Radius of the circular mask that will be used to define the background area (in pixels)", "-1"));
	bg_helical_radius = textToFloat(parser.getOption("--bg_helical_radius", "Radius of the cylindrical mask that will be used to define the background area (in pixels)", "-1."));
	white_dust_stddev = textToFloat(parser.getOption("--white_dust", "Sigma-values above which white dust will be removed (negative value means no dust removal)","-1"));
	black_dust_stddev = textToFloat(parser.getOption("--black_dust", "Sigma-values above which black dust will be removed (negative value means no dust removal)","-1"));
	do_invert_contrast = parser.checkOption("--invert_contrast", "Invert the contrast in the input images");
	fn_operate_in = parser.getOption("--operate_on", "Use this option to operate on an input stack/STAR file", "");
	fn_operate_out = parser.getOption("--operate_out", "Output rootname when operating on an input stack/STAR file", "preprocessed");

	// Initialise verb for non-parallel execution
	verb = 1;

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void Preprocessing::usage()
{
	parser.writeUsage(std::cerr);
}

void Preprocessing::initialise()
{

	if (!do_extract && fn_operate_in == "")
		REPORT_ERROR("Provide either --extract or --operate_on");

	// Set up which coordinate files to extract particles from (or to join STAR file for)
	if (do_extract)
	{
		if (verb > 0)
		{
			if (fn_star_in=="" || fn_coord_suffix=="")
				REPORT_ERROR("Preprocessing::initialise ERROR: please provide --i and --coord_suffix to extract particles");

			if (extract_size < 0)
				REPORT_ERROR("Preprocessing::initialise ERROR: please provide the size of the box to extract particle using --extract_size ");
		}



		// Read in the micrographs STAR file
		MDmics.read(fn_star_in);
		if (!MDmics.containsLabel(EMDL_MICROGRAPH_NAME))
			REPORT_ERROR("Preprocessing::initialise ERROR: Input micrograph STAR file has no rlnMicrographName column!");
		if ((do_phase_flip||do_premultiply_ctf) && !MDmics.containsLabel(EMDL_CTF_DEFOCUSU))
			REPORT_ERROR("Preprocessing::initialise ERROR: No CTF information found in the input micrograph STAR-file");

		star_has_ctf = MDmics.containsLabel(EMDL_CTF_DEFOCUSU);

		// Make sure the directory names end with a '/'
		if (fn_coord_dir != "ASINPUT" && fn_coord_dir[fn_coord_dir.length()-1] != '/')
			fn_coord_dir+="/";
		if (fn_part_dir[fn_part_dir.length()-1] != '/')
			fn_part_dir+="/";

		// Loop over all micrographs in the input STAR file and warn of coordinate file or micrograph file do not exist
		if (verb > 0)
		{
			FileName fn_mic;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmics)
			{
				MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
				FileName fn_coord = getCoordinateFileName(fn_mic);
				if (do_movie_extract)
					fn_mic = fn_mic.withoutExtension() + "_" + movie_name + ".mrcs";
				if (!exists(fn_coord))
					std::cout << "Warning: coordinate file " << fn_coord << " does not exist..." << std::endl;
				if (!exists(fn_mic))
					std::cout << "Warning: micrograph file " << fn_mic << " does not exist..." << std::endl;
			}
		}
	}

	if (do_extract || fn_operate_in != "")
	{
		// Check whether to do re-scaling
		do_rescale = (scale > 0);
		if (do_rescale && scale%2 != 0)
			REPORT_ERROR("ERROR: only re-scaling to even-sized images is allowed in RELION...");

		// Check whether to do re-windowing
		do_rewindow = (window > 0);
		if (do_rewindow && window%2 != 0)
			REPORT_ERROR("ERROR: only re-windowing to even-sized images is allowed in RELION...");

		// Check for bg_radius in case of normalisation
		if (do_normalise && bg_radius < 0)
			REPORT_ERROR("ERROR: please provide a radius for a circle that defines the background area when normalising...");
	}
}

void Preprocessing::run()
{

	if (do_extract)
		runExtractParticles();
	else if (fn_operate_in != "")
		runOperateOnInputFile();

	if (verb > 0)
		std::cout << " Done!" <<std::endl;
}


void Preprocessing::joinAllStarFiles()
{

	MetaDataTable MDout, MDonestack;

	std::cout << " Joining all metadata in one STAR file..." << std::endl;
	FileName fn_mic;
	// Re-scaled detector pixel size
	RFLOAT DStep;

	for (long int current_object1 = MDmics.firstObject();
	              current_object1 != MetaDataTable::NO_MORE_OBJECTS && current_object1 != MetaDataTable::NO_OBJECTS_STORED;
	              current_object1 = MDmics.nextObject())
	{
		// Micrograph filename
		MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);

		// Get the filename of the STAR file for just this micrograph
		FileName fn_star = getOutputFileNameRoot(fn_mic) + ".star";
		MDonestack.read(fn_star);

		for (long int current_object2 = MDonestack.firstObject();
		              current_object2 != MetaDataTable::NO_MORE_OBJECTS && current_object1 != MetaDataTable::NO_OBJECTS_STORED;
		              current_object2 = MDonestack.nextObject())
		{
			MDout.addObject( MDonestack.getObject());
		}
	}

	// Write out the joined star file
	MDout.write(fn_part_star);
	std::cout << " Written out STAR file with all particles in " << fn_part_star<< std::endl;

}

void Preprocessing::runExtractParticles()
{

	long int nr_mics = MDmics.numberOfObjects();

	int barstep;
	if (verb > 0)
	{
		std::cout << " Extracting particles from the micrographs ..." << std::endl;
		init_progress_bar(nr_mics);
		barstep = XMIPP_MAX(1, nr_mics / 60);
	}

	FileName fn_mic, fn_olddir = "";
	long int imic = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmics)
	{
		MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);

		// Check new-style outputdirectory exists and make it if not!
		FileName fn_dir = getOutputFileNameRoot(fn_mic);
		fn_dir = fn_dir.beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			int res = system(("mkdir -p " + fn_dir).c_str());
			fn_olddir = fn_dir;
		}

		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		extractParticlesFromFieldOfView(fn_mic, imic);

		imic++;
	}

	if (verb > 0)
		progress_bar(fn_coords.size());

	// Now combine all metadata in a single STAR file
	joinAllStarFiles();

}


void Preprocessing::readCoordinates(FileName fn_coord, MetaDataTable &MD)
{
    MD.clear();

    bool is_star = (fn_coord.getExtension() == "star");
    bool is_box = (fn_coord.getExtension() == "box");

    if (is_star)
	{
		MD.read(fn_coord);
		if (do_helical_segments)
		{
			if ( (!MD.containsLabel(EMDL_IMAGE_COORD_X))
					|| (!MD.containsLabel(EMDL_IMAGE_COORD_Y))
					|| (!MD.containsLabel(EMDL_ORIENT_TILT))
					|| (!MD.containsLabel(EMDL_ORIENT_PSI)) )
			{
				REPORT_ERROR("Preprocessing::readCoordinates ERROR: Extraction of helical segments: x, y, tilt or psi are missing in STAR file: " + fn_coord);
			}
		}
	}
	else
	{

		if (do_helical_segments)
			REPORT_ERROR("Preprocessing::readCoordinates ERROR: only STAR files are allowed for helical segments.");

		std::ifstream in(fn_coord.data(), std::ios_base::in);
		if (in.fail())
			REPORT_ERROR( (std::string) "Preprocessing::readCoordinates ERROR: File " + fn_coord + " does not exists" );

		// Start reading the ifstream at the top
		in.seekg(0);
		std::string line;
		int n = 0;
		while (getline(in, line, '\n'))
		{
			std::vector<std::string> words;
			tokenize(line, words);

			if (is_box)
			{
				if (words.size() < 4)
					REPORT_ERROR("Preprocessing::readCoordinates Unexpected number of words on data line of " + fn_coord);
				int num1 =  textToInteger(words[0]);
				int num2 =  textToInteger(words[1]);
				int num3 =  textToInteger(words[2]);
				int num4 =  textToInteger(words[3]);
				int xpos = num1 + num3 / 2;
				int ypos = num2 + num4 / 2;

				MD.addObject();
				MD.setValue(EMDL_IMAGE_COORD_X, (RFLOAT)xpos);
				MD.setValue(EMDL_IMAGE_COORD_Y, (RFLOAT)ypos);
			}
			else // Try reading as plain ASCII....
			{
				// The first line might be a special header (as in ximdisp or xmipp2 format)
				// But could also be a data line (as in plain text format)
				if (n==0)
				{

					int num1, num2, num3;
					bool is_data = false;

					// Ignore lines that do not have at least two integer numbers on it (at this point I do not know dimensionality yet....)
					if (words.size() > 1 && sscanf(words[0].c_str(), "%d", &num1) && sscanf(words[1].c_str(), "%d", &num2))
					{
						MD.addObject();
						MD.setValue(EMDL_IMAGE_COORD_X, (RFLOAT)num1);
						MD.setValue(EMDL_IMAGE_COORD_Y, (RFLOAT)num2);

						// It could also be a X,Y,Z coordinate...
						if (words.size() > 2 && sscanf(words[2].c_str(), "%d", &num3))
							MD.setValue(EMDL_IMAGE_COORD_Z, (RFLOAT)num3);
					}


				}
				else
				{
					// All other lines contain x, y as first two entries, and possibly a Z-coordinate as well.
					// Note the Z-coordinate will only be used for 3D micrographs (i.e. tomograms), so the third column for a 2D Ximdisp format will be ignored later on
					if (words.size() < 2)
						REPORT_ERROR("Preprocessing::readCoordinates Unexpected number of words on data line of " + fn_coord);

					int num1 = textToInteger(words[0]);
					int num2 = textToInteger(words[1]);
					int num3;

					MD.addObject();
					MD.setValue(EMDL_IMAGE_COORD_X, (RFLOAT)num1);
					MD.setValue(EMDL_IMAGE_COORD_Y, (RFLOAT)num2);
					// It could also be a X,Y,Z coordinate...
					if (words.size() > 2 && sscanf(words[2].c_str(), "%d", &num3))
						MD.setValue(EMDL_IMAGE_COORD_Z, (RFLOAT)num3);

				}
			}
			n++;

		}
		in.close();
	}
}

void Preprocessing::extractParticlesFromFieldOfView(FileName fn_mic, long int imic)
{
	// Read in the coordinates file
	MetaDataTable MDin, MDout;
	FileName fn_coord = getCoordinateFileName(fn_mic);
	readCoordinates(fn_coord, MDin);

	// Correct for bias in the picked coordinates
	if (ABS(extract_bias_x) > 0 || ABS(extract_bias_y) > 0)
	{
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			RFLOAT xcoor, ycoor;
			MDin.getValue(EMDL_IMAGE_COORD_X, xcoor);
			MDin.getValue(EMDL_IMAGE_COORD_Y, ycoor);
			xcoor += extract_bias_x;
			ycoor += extract_bias_y;
			MDin.setValue(EMDL_IMAGE_COORD_X, xcoor);
			MDin.setValue(EMDL_IMAGE_COORD_Y, ycoor);
		}
	}

	// Warn for small groups
	int npos = MDin.numberOfObjects();
	if (npos < 10)
	{
		std:: cout << "WARNING: there are only " << npos << " particles in " << fn_coord <<". Consider joining multiple micrographs into one group. "<< std::endl;
	}

	// Get movie or normal micrograph name and check it exists
	// This implicitly assumes that the movie is next to the micrograph....
	// TODO: That may not necessarily be true?
	FileName fn_img = (do_movie_extract) ? fn_mic.withoutExtension() + "_" + movie_name + ".mrcs" : fn_mic;
	if (!exists(fn_img))
	{
		std::cout << "WARNING: cannot find micrograph file " << fn_img << " which has " << npos << " particles" << std::endl;
		return;
	}

	// Read the header of the micrograph to see how many frames there are.
	Image<RFLOAT> Imic;
	Imic.read(fn_img, false, -1, false); // readData = false, select_image = -1, mapData= false, is_2D = true);

	int xdim, ydim, zdim;
	long int ndim;
	Imic.getDimensions(xdim, ydim, zdim, ndim);
	dimensionality = (zdim > 1) ? 3 : 2;
	if (dimensionality == 3 || do_movie_extract)
		do_ramp = false;

	// Name of the output particle stack
	// Add the same root as the output STAR file (that way one could extract two "families" of different particle stacks)
	FileName fn_output_img_root = getOutputFileNameRoot(fn_mic);
    // Name of this micrographs STAR file
    FileName fn_star = fn_output_img_root + ".star";

	// Just to be sure...
	if (do_movie_extract && ndim < 2)
		std::cout << "WARNING: movie " << fn_mic << " does not have multiple frames..." << std::endl;

	long int my_current_nr_images = 0;
	RFLOAT all_avg = 0;
	RFLOAT all_stddev = 0;
	RFLOAT all_minval = 99.e99;
	RFLOAT all_maxval = -99.e99;

	// To deal with default movie_last_frame value
	if (movie_last_frame < 0)
		movie_last_frame = ndim - 1;

	int n_frames = movie_last_frame - movie_first_frame + 1;
	// The total number of images to be extracted
	long int my_total_nr_images = npos * n_frames;

	for (long int iframe = movie_first_frame; iframe <= movie_last_frame; iframe += avg_n_frames)
	{
		extractParticlesFromOneFrame(MDin, fn_img, imic, iframe, n_frames, fn_output_img_root, my_current_nr_images, my_total_nr_images,
				all_avg, all_stddev, all_minval, all_maxval);

        MDout.append(MDin);
		// Keep track of total number of images extracted thus far
		my_current_nr_images += npos;

	}
    MDout.setName("images");
    MDout.write(fn_star);

}

// Actually extract particles. This can be from one (average) micrograph or from a single movie frame
void Preprocessing::extractParticlesFromOneFrame(MetaDataTable &MD,
		FileName fn_mic, int imic, int iframe, int n_frames,
		FileName fn_output_img_root, long int &my_current_nr_images, long int my_total_nr_images,
		RFLOAT &all_avg, RFLOAT &all_stddev, RFLOAT &all_minval, RFLOAT &all_maxval)
{

	Image<RFLOAT> Ipart, Imic, Itmp;

	FileName fn_frame;
	// If movies, then average over avg_n_frames
	if (n_frames > 1)
	{
		for (int ii =0; ii < avg_n_frames; ii++)
		{
			int iiframe = iframe + ii;
			// If we run over the size of the movie, then discard these frames
			if (iiframe >= movie_first_frame + n_frames)
				break;
			fn_frame.compose(iiframe + 1, fn_mic);
			if (ii==0)
			{
				Imic.read(fn_frame, true, -1, false, true); // readData = true, select_image = -1, mapData= false, is_2D = true
			}
			else
			{
				Itmp.read(fn_frame, true, -1, false, true); // readData = true, select_image = -1, mapData= false, is_2D = true
				Imic() += Itmp();
				Itmp.clear();
			}
		}
	}
	else
	{
		fn_frame = fn_mic;
		Imic.read(fn_frame);
	}

	CTF ctf;
	if (star_has_ctf || do_phase_flip || do_premultiply_ctf)
		ctf.read(MDmics, MDmics, imic);

	if (do_phase_flip || do_premultiply_ctf)
	{
		if (dimensionality != 2)
			REPORT_ERROR("extractParticlesFromFieldOfView ERROR: cannot do phase flipping as dimensionality is not 2!");

		RFLOAT mag, dstep, my_angpix;
		if (MDmics.containsLabel(EMDL_CTF_MAGNIFICATION) && MDmics.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
		{
			MDmics.getValue(EMDL_CTF_MAGNIFICATION, mag, imic);
			MDmics.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, imic);
			my_angpix = 10000. * dstep / mag;
		}
		else
			my_angpix = angpix;

		// FFTW only does square images
		long int ori_xsize = XSIZE(Imic());
		long int ori_ysize = YSIZE(Imic());

		if (ori_xsize > ori_ysize)
			rewindow(Imic, ori_xsize);
		else if (ori_ysize > ori_xsize)
			rewindow(Imic, ori_ysize);

		// Now do the actual phase flipping or CTF-multiplication
		MultidimArray<Complex> FTmic;
		FourierTransformer transformer;

		transformer.FourierTransform(Imic(), FTmic, false);

		RFLOAT xs = (RFLOAT)XSIZE(Imic()) * my_angpix;
		RFLOAT ys = (RFLOAT)YSIZE(Imic()) * my_angpix;
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(FTmic)
		{
			RFLOAT x = (RFLOAT)jp / xs;
			RFLOAT y = (RFLOAT)ip / ys;
			DIRECT_A2D_ELEM(FTmic, i, j) *= ctf.getCTF(x, y, false, do_phase_flip, do_ctf_intact_first_peak, false); // true=do_only_phase_flip
		}

		transformer.inverseFourierTransform(FTmic, Imic());

		if (ori_xsize != ori_ysize)
			Imic().window(FIRST_XMIPP_INDEX(ori_ysize), FIRST_XMIPP_INDEX(ori_xsize),
					LAST_XMIPP_INDEX(ori_ysize),  LAST_XMIPP_INDEX(ori_xsize));
	}

	// Jun24,2015 - Shaoda, helical segments
	if (do_helical_segments)
	{
		if ( (!MD.containsLabel(EMDL_IMAGE_COORD_X))
				|| (!MD.containsLabel(EMDL_IMAGE_COORD_Y))
				|| (!MD.containsLabel(EMDL_ORIENT_TILT))
				|| (!MD.containsLabel(EMDL_ORIENT_PSI)) )
		{
			REPORT_ERROR("ERROR: Extraction of helical segments: x, y, tilt and psi are missing in STAR files!");
		}
	}

	// Now window all particles from the micrograph
	int ipos = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		RFLOAT dxpos, dypos, dzpos;
		long int xpos, ypos, zpos;
		long int x0, xF, y0, yF, z0, zF;
		MD.getValue(EMDL_IMAGE_COORD_X, dxpos);
		MD.getValue(EMDL_IMAGE_COORD_Y, dypos);
		xpos = (long int)dxpos;
		ypos = (long int)dypos;
		x0 = xpos + FIRST_XMIPP_INDEX(extract_size);
		xF = xpos + LAST_XMIPP_INDEX(extract_size);
		y0 = ypos + FIRST_XMIPP_INDEX(extract_size);
		yF = ypos + LAST_XMIPP_INDEX(extract_size);
		if (dimensionality == 3)
		{
			MD.getValue(EMDL_IMAGE_COORD_Z, dzpos);
			zpos = (long int)dzpos;
			z0 = zpos + FIRST_XMIPP_INDEX(extract_size);
			zF = zpos + LAST_XMIPP_INDEX(extract_size);
		}

		// extract one particle in Ipart
		if (dimensionality == 3)
			Imic().window(Ipart(), z0, y0, x0, zF, yF, xF);
		else
			Imic().window(Ipart(), y0, x0, yF, xF);

		// Discard particles that are completely outside the micrograph and print a warning
		if (yF < 0 || y0 >= YSIZE(Imic()) || xF < 0 || x0 >= XSIZE(Imic()) ||
				(dimensionality==3 &&
						(zF < 0 || z0 >= ZSIZE(Imic()))
				)
			)
		{
			REPORT_ERROR("Preprocessing::extractParticlesFromOneFrame ERROR: particle" + integerToString(ipos+1) + " lies completely outside micrograph " + fn_mic);
		}
		else
		{
			// Check boundaries: fill pixels outside the boundary with the nearest ones inside
			// This will create lines at the edges, rather than zeros
			Ipart().setXmippOrigin();

			// X-boundaries
			if (x0 < 0 || xF >= XSIZE(Imic()) )
			{
				FOR_ALL_ELEMENTS_IN_ARRAY3D(Ipart())
				{
					if (j + xpos < 0)
						A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), k, i, -xpos);
					else if (j + xpos >= XSIZE(Imic()))
						A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), k, i, XSIZE(Imic()) - xpos - 1);
				}
			}

			// Y-boundaries
			if (y0 < 0 || yF >= YSIZE(Imic()))
			{
				FOR_ALL_ELEMENTS_IN_ARRAY3D(Ipart())
				{
					if (i + ypos < 0)
						A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), k, -ypos, j);
					else if (i + ypos >= YSIZE(Imic()))
						A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), k, YSIZE(Imic()) - ypos - 1, j);
				}
			}

			if (dimensionality == 3)
			{
				// Z-boundaries
				if (z0 < 0 || zF >= ZSIZE(Imic()))
				{
					FOR_ALL_ELEMENTS_IN_ARRAY3D(Ipart())
					{
						if (k + zpos < 0)
							A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), -zpos, i, j);
						else if (k + zpos >= ZSIZE(Imic()))
							A3D_ELEM(Ipart(), k, i, j) = A3D_ELEM(Ipart(), ZSIZE(Imic()) - zpos - 1, i, j);
					}
				}
			}

			//
			if (dimensionality == 3 && do_project_3d)
			{
				// Project the 3D sub-tomogram into a 2D particle again
				Image<RFLOAT> Iproj(YSIZE(Ipart()), XSIZE(Ipart()));
				Iproj().setXmippOrigin();
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Ipart())
				{
					DIRECT_A2D_ELEM(Iproj(), i, j) += DIRECT_A3D_ELEM(Ipart(), k, i, j);
				}
				Ipart = Iproj;
			}

			// performPerImageOperations will also append the particle to the output stack in fn_stack
			// Jun24,2015 - Shaoda, extract helical segments
			RFLOAT tilt_deg, psi_deg;
			tilt_deg = psi_deg = 0.;
			if (do_helical_segments)
			{
				MD.getValue(EMDL_ORIENT_TILT, tilt_deg);
				MD.getValue(EMDL_ORIENT_PSI, psi_deg);
			}
			performPerImageOperations(Ipart, fn_output_img_root, n_frames, my_current_nr_images + ipos, my_total_nr_images,
					tilt_deg, psi_deg,
					all_avg, all_stddev, all_minval, all_maxval);

			// Also store all the particles information in the STAR file
			FileName fn_img;
			if (Ipart().getDim() == 3)
				fn_img.compose(fn_output_img_root, my_current_nr_images + ipos + 1, "mrc");
			else
				fn_img.compose(my_current_nr_images + ipos + 1, fn_output_img_root + ".mrcs"); // start image counting in stacks at 1!
			if (do_movie_extract)
			{
				FileName fn_part, fn_mic;
				long int dum;
				fn_frame.decompose(dum, fn_mic);
				FileName fn_mic_ori = fn_output_img_root.without("_"+movie_name);
				fn_part.compose(ipos + 1,  fn_output_img_root + ".mrcs"); // start image counting in stacks at 1!
				// for automated re-alignment of particles in relion_refine: have rlnParticleName equal to rlnImageName in non-movie star file
				MD.setValue(EMDL_PARTICLE_ORI_NAME, fn_part);
			}
			MD.setValue(EMDL_IMAGE_NAME, fn_img);
			MD.setValue(EMDL_MICROGRAPH_NAME, fn_frame);

			// Also fill in the CTF parameters
			if (star_has_ctf)
			{
				ctf.write(MD);
				RFLOAT mag, dstep;
				if (MDmics.containsLabel(EMDL_CTF_MAGNIFICATION))
				{
					MDmics.getValue(EMDL_CTF_MAGNIFICATION, mag, imic);
					MD.setValue(EMDL_CTF_MAGNIFICATION, mag);
				}
				if (MDmics.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
				{
					MDmics.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, imic);
					if (do_rescale)
						dstep *= (RFLOAT)extract_size/(RFLOAT)scale;
					MD.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
				}
			}

		}
		ipos++;
	}


}


void Preprocessing::runOperateOnInputFile()
{
	Image<RFLOAT> Ipart, Iout;
	MetaDataTable MD;
	long int Nimg;

	FileName fn_stack = fn_operate_out.withoutExtension()+".mrcs";
	FileName fn_star = fn_operate_out.withoutExtension()+".star";

	if (fn_operate_in.isStarFile())
	{
		// Readt STAR file and get total number of images
		MD.read(fn_operate_in);
		Nimg = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			Nimg++;
		}
		MD.firstObject(); // reset pointer to the first object in the table

		// Jun25,2015 - Shaoda, helical segments
		if (do_helical_segments)
		{
			if ( (!MD.containsLabel(EMDL_ORIENT_TILT))
					|| (!MD.containsLabel(EMDL_ORIENT_PSI)) )
			{
				REPORT_ERROR("ERROR: Operations on helical segments: only STAR files with tilt and psi angles are supported!");
			}
		}
	}
	else
	{
		// Jun25,2015 - Shaoda, helical segments
		if (do_helical_segments)
			REPORT_ERROR("ERROR: Operations on helical segments: only STAR files are supported!");

		// Read the header of the stack to see how many images there
		Iout.read(fn_operate_in, false);
		Nimg = NSIZE(Iout());
	}

	RFLOAT all_avg = 0;
	RFLOAT all_stddev = 0;
	RFLOAT all_minval = 99.e99;
	RFLOAT all_maxval = -99.e99;
	init_progress_bar(Nimg);
	int barstep = XMIPP_MAX(1, Nimg / 120);
	for (long int i = 0; i < Nimg; i++)
	{
		FileName fn_tmp;

		// Read in individual miages from the stack
		Ipart.clear();
		if (fn_operate_in.isStarFile())
		{
			MD.getValue(EMDL_IMAGE_NAME, fn_tmp);
			Ipart.read(fn_tmp);

			// Set the new name at this point in the MDtable, e.g. as 000001@out.mrcs
			fn_tmp.compose(i+1,fn_stack);
			MD.setValue(EMDL_IMAGE_NAME, fn_tmp);
			if (i < (Nimg - 1))
				MD.nextObject();
		}
		else
		{
			Ipart.read(fn_operate_in, true, i);
			// Set the new name at this point in the MDtable, e.g. as 000001@out.mrcs
			fn_tmp.compose(i+1,fn_stack);
			MD.addObject();
			MD.setValue(EMDL_IMAGE_NAME, fn_tmp);
		}

		// Jun24,2015 - Shaoda, extract helical segments
		RFLOAT tilt_deg, psi_deg;
		tilt_deg = psi_deg = 0.;
		if (do_helical_segments)
		{
			MD.getValue(EMDL_ORIENT_TILT, tilt_deg);
			MD.getValue(EMDL_ORIENT_PSI, psi_deg);
		}
		performPerImageOperations(Ipart, fn_stack, 1, i, Nimg,
				tilt_deg, psi_deg,
				all_avg, all_stddev, all_minval, all_maxval);

		// progress bar
		if (i % barstep == 0) progress_bar(i);

	}
	progress_bar(Nimg);

	std::cout << " Done writing to " << fn_stack << std::endl;
	MD.setName("images");
	MD.write(fn_star);
	std::cout << " Also written a STAR file with the image names as " << fn_star << std::endl;

}


void Preprocessing::performPerImageOperations(
		Image<RFLOAT> &Ipart,
		FileName fn_output_img_root,
		int nframes,
		long int image_nr,
		long int nr_of_images,
		RFLOAT tilt_deg,
		RFLOAT psi_deg,
		RFLOAT &all_avg,
		RFLOAT &all_stddev,
		RFLOAT &all_minval,
		RFLOAT &all_maxval)
{

	Ipart().setXmippOrigin();

	if (do_rescale) rescale(Ipart, scale);

	if (do_rewindow) rewindow(Ipart, window);

	// Jun24,2015 - Shaoda, helical segments
	if (do_normalise)
	{
		normalise(Ipart, bg_radius, white_dust_stddev, black_dust_stddev, do_ramp,
				do_helical_segments, bg_helical_radius, tilt_deg, psi_deg);
	}

	if (do_invert_contrast) invert_contrast(Ipart);

	// For movies: multiple the image intensities by sqrt(nframes) so the stddev in the average of the normalised frames is again 1
	if (nframes > 1)
		Ipart() *= sqrt((RFLOAT)nframes/(RFLOAT)avg_n_frames);

	// Calculate mean, stddev, min and max
	RFLOAT avg, stddev, minval, maxval;
	Ipart().computeStats(avg, stddev, minval, maxval);

	if (Ipart().getDim() == 3)
	{
		Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, minval);
		Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, maxval);
		Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, avg);
		Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, stddev);

		// Write one mrc file for every subtomogram
		FileName fn_img;
		fn_img.compose(fn_output_img_root, image_nr + 1, "mrc");
		Ipart.write(fn_img);

	}
	else
	{
		// Keep track of overall statistics
		all_minval = XMIPP_MIN(minval, all_minval);
		all_maxval = XMIPP_MAX(maxval, all_maxval);
		all_avg	+= avg;
		all_stddev += stddev*stddev;

		// Last particle: reset the min, max, avg and stddev values in the main header
		if (image_nr == (nr_of_images - 1))
		{
			all_avg /= nr_of_images;
			all_stddev = sqrt(all_stddev/nr_of_images);
			Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MIN, all_minval);
			Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, all_maxval);
			Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_AVG, all_avg);
			Ipart.MDMainHeader.setValue(EMDL_IMAGE_STATS_STDDEV, all_stddev);
		}

		// Write this particle to the stack on disc
		// First particle: write stack in overwrite mode, from then on just append to it
		if (image_nr == 0)
			Ipart.write(fn_output_img_root+".mrcs", -1, (nr_of_images > 1), WRITE_OVERWRITE);
		else
			Ipart.write(fn_output_img_root+".mrcs", -1, false, WRITE_APPEND);
	}

}

// Get the coordinate filename from the micrograph filename
FileName Preprocessing::getCoordinateFileName(FileName fn_mic)
{

	FileName uniqdate;
	size_t slashpos = findUniqueDateSubstring(fn_mic, uniqdate);
	FileName fn_mic_nouniqdate = (slashpos!= std::string::npos) ? fn_mic.substr(slashpos+15) : fn_mic;
	fn_mic_nouniqdate = fn_mic_nouniqdate.withoutExtension();
	FileName fn_coord = fn_coord_dir + fn_mic_nouniqdate + fn_coord_suffix;
	return fn_coord;
}

// Get the coordinate filename from the micrograph filename
FileName Preprocessing::getOutputFileNameRoot(FileName fn_mic)
{

	FileName uniqdate;
	size_t slashpos = findUniqueDateSubstring(fn_mic, uniqdate);
	FileName fn_mic_nouniqdate = (slashpos!= std::string::npos) ? fn_mic.substr(slashpos+15) : fn_mic;
	fn_mic_nouniqdate = fn_mic_nouniqdate.withoutExtension();
	FileName fn_part = fn_part_dir + fn_mic_nouniqdate;
	return fn_part;
}


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
	fn_out = parser.getOption("--o", "Output rootname", "particles");
	// Dont allow for directories here!
	fn_out = fn_out.removeDirectories();

	int particle_section = parser.addSection("Particle selection");
	fns_coords_in = parser.getOption("--coord_files", "The coordinate files for all particles to be output (may contain wildcards e.g. \"mics/*.box\" or \"mics/???.star\" )","");
	fn_star_in = parser.getOption("--mic_star", "The STAR file with all (selected) micrographs to extract particles from","");
	fn_pick_suffix = parser.getOption("--coord_suffix", "The suffix for the coordinate files, e.g. \"_picked.star\" or \".box\"","");

	int extract_section = parser.addSection("Particle extraction");
	do_extract = parser.checkOption("--extract", "Extract all particles from the micrographs");
	extract_size = textToInteger(parser.getOption("--extract_size", "Size of the box to extract the particles in (in pixels)", "-1"));
	extract_bias_x  = textToInteger(parser.getOption("--extract_bias_x", "Bias in X-direction of picked particles (this value in pixels will be added to the coords)", "0"));
	extract_bias_y  = textToInteger(parser.getOption("--extract_bias_y", "Bias in Y-direction of picked particles (this value in pixels will be added to the coords)", "0"));
	do_movie_extract = parser.checkOption("--extract_movies", "Extract particles from movie stacks (e.g. from DDDs)");
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
	bg_radius = textToInteger(parser.getOption("--bg_radius", "Radius of the circular mask that will be used to define the background area (in pixels)", "-1"));
	white_dust_stddev = textToFloat(parser.getOption("--white_dust", "Sigma-values above which white dust will be removed (negative value means no dust removal)","-1"));
	black_dust_stddev = textToFloat(parser.getOption("--black_dust", "Sigma-values above which black dust will be removed (negative value means no dust removal)","-1"));
	do_invert_contrast = parser.checkOption("--invert_contrast", "Invert the contrast in the input images");
	fn_operate_in = parser.getOption("--operate_on", "The operations above are applied to all extracted particles. Use this option to operate on an input stack/STAR file", "");

	// Initialise verb for non-parallel execution
	verb = 1;

	if (!( checkParameter(argc, argv, "--o") || checkParameter(argc, argv, "--operate_on") ))
		REPORT_ERROR("Provide either --o or --operate_on");

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

	if (fns_coords_in == "" && fn_star_in == "" && fn_operate_in == "" && !do_join_starfile)
		REPORT_ERROR("Provide either --extract, --join_starfile or --operate_on");

	// Set up which coordinate files to extract particles from (or to join STAR file for)
	do_join_starfile = false;
	if (fns_coords_in != "" || fn_star_in != "")
	{
		do_join_starfile = true;
		if (do_extract && verb > 0)
		{
			if (extract_size < 0)
				REPORT_ERROR("Preprocessing::initialise ERROR: please provide the size of the box to extract particle using --extract_size ");

			std::cout << " Extract particles based on the following coordinate files: " << std::endl;
		}
		else if (!do_extract && verb > 0)
		{
			std::cout << " Creating output STAR file for particles based on the following coordinate files: " << std::endl;
		}

		// Get the filenames of all micrographs to be processed by CTFFIND
		if (fns_coords_in != "")
		{
			fns_coords_in.globFiles(fn_coords);
		}
		else if (fn_star_in != "")
		{
			MetaDataTable MDmics;
			MDmics.read(fn_star_in);
			if (!MDmics.containsLabel(EMDL_MICROGRAPH_NAME))
				REPORT_ERROR("Preprocessing::initialise ERROR: Input micrograph STAR file has no rlnMicrographName column!");
			fn_coords.clear();
			FileName fn_mic;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmics)
			{
				MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
				fn_mic = fn_mic.withoutExtension() + fn_pick_suffix;
				if (exists(fn_mic))
					fn_coords.push_back(fn_mic);
				else if (verb > 0)
					std::cout << "Warning: coordinate file " << fn_mic << " does not exist..." << std::endl;
			}
		}

		if (verb > 0)
			for(unsigned  int  i = 0; i < fn_coords.size(); ++i)
				std::cout << "  * " << fn_coords[i] << std::endl;
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

	if (do_join_starfile)
		joinAllStarFiles();

	if (fn_operate_in != "")
		runOperateOnInputFile(fn_operate_in);

	if (verb > 0)
		std::cout << " Done!" <<std::endl;
}


void Preprocessing::joinAllStarFiles()
{

	MetaDataTable MDout, MDonestack;

	// Fix order of the labels in the output file
	MDout.addLabel(EMDL_MICROGRAPH_NAME);
	MDout.addLabel(EMDL_IMAGE_COORD_X);
	MDout.addLabel(EMDL_IMAGE_COORD_Y);
	MDout.addLabel(EMDL_IMAGE_NAME);

	std::cout << " Joining all metadata in one STAR file..." << std::endl;
	bool has_other_ctfs, has_this_ctf;
	double defU, defV, defAng, CC, HT, CS, AmpCnst, XMAG, DStep;
	has_other_ctfs = false;
	FileName prev_fn_mic="";
	for (long int ipos = 0; ipos < fn_coords.size(); ipos++)
    {
		FileName fn_mic = getMicrographNameFromRootName((fn_coords[ipos]).withoutExtension());
		// If the micrograph did not exist, particles were not extracted: just continue with the next one
		if (fn_mic == "")
			continue;

		if (fn_mic != prev_fn_mic)
		{
			FileName fn_microot = getRootNameFromMicrographName(fn_mic);
			// Gather the results from ctffind
			has_this_ctf = getCtffindResults(fn_microot, defU, defV, defAng, CC,
					HT, CS, AmpCnst, XMAG, DStep);

		}
		prev_fn_mic = fn_mic;

			// Re-scaled detector pixel size
		if (has_this_ctf && do_rescale)
                    DStep *= (double)extract_size/(double)scale;

		if (ipos == 0 && has_this_ctf)
		{
			// Set has_other_ctfs to true if the first micrograph has a logfile
			// In that case add all CTF labels to the output MDtable
			has_other_ctfs = true;
			MDout.addLabel(EMDL_CTF_DEFOCUSU);
			MDout.addLabel(EMDL_CTF_DEFOCUSV);
			MDout.addLabel(EMDL_CTF_DEFOCUS_ANGLE);
			MDout.addLabel(EMDL_CTF_VOLTAGE);
			MDout.addLabel(EMDL_CTF_CS);
			MDout.addLabel(EMDL_CTF_Q0);
			MDout.addLabel(EMDL_CTF_MAGNIFICATION);
			MDout.addLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE);
			MDout.addLabel(EMDL_CTF_FOM);
		}

		if (!has_this_ctf && has_other_ctfs)
			REPORT_ERROR("joinAllStarFiles%ERROR: Exiting because of missing CTFFIND logfiles for micrograph " + fn_mic);

		if (has_this_ctf && !has_other_ctfs)
			REPORT_ERROR("joinAllStarFiles%ERROR: Exiting because of missing CTFFIND logfiles ...");


		FileName fn_star = "Particles/" + fn_mic.withoutExtension() + "_" + fn_out + ".star";
		MDonestack.read(fn_star);

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDonestack)
		{
			// This was double, right?
			//MDonestack.setValue(EMDL_MICROGRAPH_NAME, fn_mic);

			if (has_this_ctf)
			{
				MDonestack.setValue(EMDL_CTF_DEFOCUSU, defU);
				MDonestack.setValue(EMDL_CTF_DEFOCUSV, defV);
				MDonestack.setValue(EMDL_CTF_DEFOCUS_ANGLE, defAng);
				MDonestack.setValue(EMDL_CTF_VOLTAGE, HT);
				MDonestack.setValue(EMDL_CTF_Q0, AmpCnst);
				MDonestack.setValue(EMDL_CTF_CS, CS);
				MDonestack.setValue(EMDL_CTF_MAGNIFICATION, XMAG);
				MDonestack.setValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, DStep);
				MDonestack.setValue(EMDL_CTF_FOM, CC);
			}

			// Add the entire line to the joined STAR-file
			MDout.addObject(MDonestack.getObject());

			// Remove the individual star file to clean up
			//remove(fn_star.c_str());
		}

    }

	// Write out the joined star file
	FileName fn_tmp;
	if (do_movie_extract)
		fn_tmp = fn_out + "_" + fn_movie + ".star";
	else
		fn_tmp = fn_out + ".star";
	MDout.write(fn_tmp);
	std::cout << " Written out STAR file with all particles in " << fn_tmp << std::endl;

}

void Preprocessing::runExtractParticles()
{

	int barstep;
	if (verb > 0)
	{
		std::cout << " Extracting particles from the micrographs ..." << std::endl;
		init_progress_bar(fn_coords.size());
		barstep = XMIPP_MAX(1, fn_coords.size() / 60);
	}

	FileName fn_olddir = "";
	for (long int ipos = 0; ipos < fn_coords.size(); ipos++)
    {

		FileName fn_dir = "Particles/" + fn_coords[ipos].beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			system(("mkdir -p " + fn_dir).c_str());
			fn_olddir = fn_dir;
		}

		if (verb > 0 && ipos % barstep == 0)
			progress_bar(ipos);

    	extractParticlesFromFieldOfView(fn_coords[ipos]);
	}

	if (verb > 0)
		progress_bar(fn_coords.size());

}


void Preprocessing::readCoordinates(FileName fn_coord, MetaDataTable &MD)
{
    MD.clear();

    bool is_star = (fn_coord.getExtension() == "star");
    bool is_box = (fn_coord.getExtension() == "box");

    if (is_star)
	{
		MD.read(fn_coord);
		/*
		bool has_z = MD.containsLabel(EMDL_IMAGE_COORD_Z);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			Matrix1D<int> onepos(3);
			MD.getValue(EMDL_IMAGE_COORD_X, XX(onepos));
			MD.getValue(EMDL_IMAGE_COORD_Y, YY(onepos));
			if (has_z)
				MD.getValue(EMDL_IMAGE_COORD_Z, ZZ(onepos));
			all_pos.push_back(onepos);
		}
		*/
	}
	else
	{
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
				MD.setValue(EMDL_IMAGE_COORD_X, (double)xpos);
				MD.setValue(EMDL_IMAGE_COORD_Y, (double)ypos);
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
						MD.setValue(EMDL_IMAGE_COORD_X, (double)num1);
						MD.setValue(EMDL_IMAGE_COORD_Y, (double)num2);

						// It could also be a X,Y,Z coordinate...
						if (words.size() > 2 && sscanf(words[2].c_str(), "%d", &num3))
							MD.setValue(EMDL_IMAGE_COORD_Z, (double)num3);
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
					MD.setValue(EMDL_IMAGE_COORD_X, (double)num1);
					MD.setValue(EMDL_IMAGE_COORD_Y, (double)num2);
					// It could also be a X,Y,Z coordinate...
					if (words.size() > 2 && sscanf(words[2].c_str(), "%d", &num3))
						MD.setValue(EMDL_IMAGE_COORD_Z, (double)num3);

				}
			}
			n++;

		}
		in.close();
	}
}

void Preprocessing::extractParticlesFromFieldOfView(FileName fn_coord)
{
	MetaDataTable MDin, MDout;

	// Read in the coordinates file
	readCoordinates(fn_coord, MDin);

	// Correct for bias in the picked coordinates
	if (ABS(extract_bias_x) > 0 || ABS(extract_bias_y) > 0)
	{
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
		{
			double xcoor, ycoor;
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

	// Check the micrograph exists
	FileName fn_mic;
	fn_mic = getMicrographNameFromRootName(fn_coord.withoutExtension());
	// Return if the micrograph does not exist
	if (fn_mic == "")
	{
		std::cout << "WARNING: cannot find micrograph for coordinate file " << fn_coord << " with " << npos << " particles" << std::endl;
		return;
	}

	// Read the header of the micrograph to see how many frames there are.
	Image<double> Imic;
	Imic.read(fn_mic, false, -1, false); // readData = false, select_image = -1, mapData= false, is_2D = true);

	int xdim, ydim, zdim;
	long int ndim;
	Imic.getDimensions(xdim, ydim, zdim, ndim);
	dimensionality = (zdim > 1) ? 3 : 2;

	// Name of the output stack
	// Add the same root as the output STAR file (that way one could extract two "families" of different particle stacks)
	FileName fn_output_img_root = "Particles/" + fn_mic.withoutExtension() + "_" + fn_out;
	// Name of this micrographs STAR file
	FileName fn_star = fn_output_img_root + ".star";

	// Just to be sure...
	if (do_movie_extract && ndim < 2)
		std::cout << "WARNING: movie " << fn_mic << " does not have multiple frames..." << std::endl;

	long int my_current_nr_images = 0;
	double all_avg = 0;
	double all_stddev = 0;
	double all_minval = 99.e99;
	double all_maxval = -99.e99;

	// To deal with default movie_last_frame value
	if (movie_last_frame < 0)
		movie_last_frame = ndim - 1;

	int n_frames = movie_last_frame - movie_first_frame + 1;
	// The total number of images to be extracted
	long int my_total_nr_images = npos * n_frames;

	for (long int iframe = movie_first_frame; iframe <= movie_last_frame; iframe += avg_n_frames)
	{
		extractParticlesFromOneFrame(MDin, fn_mic, iframe, n_frames, fn_output_img_root, my_current_nr_images, my_total_nr_images,
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
		FileName fn_mic, int iframe, int n_frames,
		FileName fn_output_img_root, long int &my_current_nr_images, long int my_total_nr_images,
		double &all_avg, double &all_stddev, double &all_minval, double &all_maxval)
{

	Image<double> Ipart, Imic, Itmp;

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

	// Now window all particles from the micrograph
	int ipos = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{
		double dxpos, dypos, dzpos;
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
				Image<double> Iproj(YSIZE(Ipart()), XSIZE(Ipart()));
				Iproj().setXmippOrigin();
				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Ipart())
				{
					DIRECT_A2D_ELEM(Iproj(), i, j) += DIRECT_A3D_ELEM(Ipart(), k, i, j);
				}
				Ipart = Iproj;
			}

			// performPerImageOperations will also append the particle to the output stack in fn_stack
			performPerImageOperations(Ipart, fn_output_img_root, n_frames, my_current_nr_images + ipos, my_total_nr_images, all_avg, all_stddev, all_minval, all_maxval);

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
				fn_part.compose(ipos + 1,  "Particles/" + getRootNameFromMicrographName(fn_mic)); // start image counting in stacks at 1!
				// for automated re-alignment of particles in relion_refine: have rlnParticleName equal to rlnImageName in non-movie star file
				fn_part += "_" + fn_out + ".mrcs";
				MD.setValue(EMDL_PARTICLE_ORI_NAME, fn_part);
			}
			MD.setValue(EMDL_IMAGE_NAME, fn_img);
			MD.setValue(EMDL_MICROGRAPH_NAME, fn_frame);
		}
		ipos++;
	}


}


void Preprocessing::runOperateOnInputFile(FileName fn_operate_on)
{
	Image<double> Ipart, Iout;
	MetaDataTable MD;
	long int Nimg;

	FileName fn_stack = fn_out+".mrcs";
	FileName fn_star = fn_out+".star";

	if (fn_operate_on.isStarFile())
	{
		// Readt STAR file and get total number of images
		MD.read(fn_operate_on);
		Nimg = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			Nimg++;
		}
		MD.firstObject(); // reset pointer to the first object in the table
	}
	else
	{
		// Read the header of the stack to see how many images there
		Iout.read(fn_operate_on, false);
		Nimg = NSIZE(Iout());
	}

	double all_avg = 0;
	double all_stddev = 0;
	double all_minval = 99.e99;
	double all_maxval = -99.e99;
	init_progress_bar(Nimg);
	int barstep = XMIPP_MAX(1, Nimg / 120);
	for (long int i = 0; i < Nimg; i++)
	{
		FileName fn_tmp;

		// Read in individual miages from the stack
		Ipart.clear();
		if (fn_operate_on.isStarFile())
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
			Ipart.read(fn_operate_on, true, i);
			// Set the new name at this point in the MDtable, e.g. as 000001@out.mrcs
			fn_tmp.compose(i+1,fn_stack);
			MD.addObject();
			MD.setValue(EMDL_IMAGE_NAME, fn_tmp);
		}

		performPerImageOperations(Ipart, fn_stack, 1, i, Nimg, all_avg, all_stddev, all_minval, all_maxval);

		// progress bar
		if (i % barstep == 0) progress_bar(i);

	}
	progress_bar(Nimg);

	std::cout << " Done writing to " << fn_stack << std::endl;
	MD.setName("images");
	MD.write(fn_star);
	std::cout << " Also written a STAR file with the image names as " << fn_star << std::endl;

}


void Preprocessing::performPerImageOperations(Image<double> &Ipart, FileName fn_output_img_root, int nframes, long int image_nr, long int nr_of_images,
		double &all_avg, double &all_stddev, double &all_minval, double &all_maxval)
{

	Ipart().setXmippOrigin();

	if (do_rescale) rescale(Ipart, scale);

	if (do_rewindow) rewindow(Ipart, window);

	if (do_normalise) normalise(Ipart, bg_radius, white_dust_stddev, black_dust_stddev);

	if (do_invert_contrast) invert_contrast(Ipart);

	// For movies: multiple the image intensities by sqrt(nframes) so the stddev in the average of the normalised frames is again 1
	if (nframes > 1)
		Ipart() *= sqrt((double)nframes/(double)avg_n_frames);

	// Calculate mean, stddev, min and max
	double avg, stddev, minval, maxval;
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

FileName Preprocessing::getMicrographNameFromRootName(FileName fn_root)
{
	FileName fn_mic, fn_mic_nos;
	//Search the unique part of the micrograph name (i.e. coord name may have additional text after the micrograph name without extension...
	// e.g. fn_root="mic001_all.pos" may correspond to mic001.mrc
	for (int i = 0; i < fn_root.length(); i++)
	{
		if (do_movie_extract)
		{
			fn_mic = fn_root.substr(0, fn_root.length() - i) + "_" + fn_movie + ".mrcs";
			// For movies also allow name without the s of .mrcs
			fn_mic_nos = fn_root.substr(0, fn_root.length() - i) + "_" + fn_movie + ".mrc";
		}
		else
			fn_mic = fn_root.substr(0, fn_root.length() - i) + ".mrc";
		std::vector<FileName> fn_mics;
		if (fn_mic.globFiles(fn_mics) == 1)
		{
			fn_mic = fn_mics[0];
			break;
		}
		else if (do_movie_extract && fn_mic_nos.globFiles(fn_mics) == 1)
		{
			// For movies also allow name without the s of .mrcs
			fn_mic = fn_mics[0];
			break;
		}
		if (i == fn_root.length() - 1)
		{
			return "";
		}
	}

	return fn_mic;

}

FileName Preprocessing::getRootNameFromMicrographName(FileName fn_mic)
{
	if (do_movie_extract)
	{
		if (fn_mic.contains(".mrcs"))
			return fn_mic.without("_" + fn_movie + ".mrcs");
		else if (fn_mic.contains(".mrc"))
			return fn_mic.without("_" + fn_movie + ".mrc");
		else
			REPORT_ERROR("Preprocessing::getRootNameFromMicrographName ERROR: movie name does not end in \"_\" + movie-identifier + \".mrc\" or \".mrcs\"");
	}
	else
		return fn_mic.without(".mrc");


}

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


#include "helix_inimodel2d.h"
//#define DEBUG


void HelixAlignerModel::initialise(int nr_classes, int ydim, int xdim)
{

	MultidimArray<RFLOAT> tmp;
	tmp.initZeros(ydim, xdim);
	tmp.setXmippOrigin();
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		Aref.push_back(tmp);
		Asum.push_back(tmp);
		Asumw.push_back(tmp);
		pdf.push_back(0);
	}
	tmp.initZeros(ydim,ydim);
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		Arec.push_back(tmp);
	}
}

void HelixAlignerModel::initZeroSums()
{
	for (int iclass = 0; iclass < Asum.size(); iclass++)
	{
		Asum[iclass].initZeros();
		Asumw[iclass].initZeros();
		pdf[iclass] = 0.;
	}
}

// Cleaning up
void HelixAlignerModel::clear()
{
	Aref.clear();
	Arec.clear();
	Asum.clear();
	Asumw.clear();
	pdf.clear();
}

void HelixAligner::clear()
{
	// TODO, clean up..

}

void HelixAligner::usage()
{
	parser.writeUsage(std::cout);
}

void HelixAligner::parseInitial(int argc, char **argv)
{

	parser.setCommandLine(argc, argv);

	// General optimiser I/O stuff
    int general_section = parser.addSection("General options");

    fn_out = parser.getOption("--o", "Output rootname","");
	fn_imgs = parser.getOption("--i", " STAR file with the input images and orientation parameters","");
    // deactivate fn_mics approach: never really worked...
	fn_mics = "";
	/*
	fn_mics = parser.getOption("--mic", "OR: STAR file with the input micrographs","");
	fn_coord_suffix = parser.getOption("--coord_suffix", "The suffix for the start-end coordinate files, e.g. \"_picked.star\" or \".box\"","");
	fn_coord_dir = parser.getOption("--coord_dir", "The directory where the coordinate files are (default is same as micrographs)", "ASINPUT");
	extract_width = textToInteger(parser.getOption("--extract_width", "Width (in pixels) of the images for the helices to be extracted ", "100"));
	*/
    int param_section = parser.addSection("Parameters");
    crossover_distance = textToFloat(parser.getOption("--crossover_distance", "Distance in Angstroms between 2 cross-overs",""));
    nr_iter = textToInteger(parser.getOption("--iter", "Maximum number of iterations to perform", "10"));
    nr_classes = textToInteger(parser.getOption("--K", "Number of classes", "1"));
    angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms (default take from STAR file)", "-1"));
    maxres = textToFloat(parser.getOption("--maxres", "Limit calculations to approximately this resolution in Angstroms", "-1"));
    max_shift_A = textToFloat(parser.getOption("--search_shift", "How many Angstroms to search translations perpendicular to helical axis?", "0"));
    max_rotate = textToFloat(parser.getOption("--search_angle", "How many degrees to search in-plane rotations?", "0"));
    step_rotate = textToFloat(parser.getOption("--step_angle", "The step size (in degrees) of the rotational searches", "1"));
    fn_inimodel = parser.getOption("--iniref", "An initial model to starting optimisation path", "");
    symmetry = textToInteger(parser.getOption("--sym", "Order of symmetry in the 2D xy-slice?", "1"));
    max_smear = textToInteger(parser.getOption("--smear", "Smear out each image along X to ensure continuity", "0"));
    random_seed = textToInteger(parser.getOption("--random_seed", "Random seed (default is with clock)", "-1"));
    search_size = textToInteger(parser.getOption("--search_size", "Search this many pixels up/down of the target downscaled size to fit best crossover distance", "5"));
    mask_diameter = textToFloat(parser.getOption("--mask_diameter", "The diameter (A) of a mask to be aplpied to the 2D reconstruction", "-1"));
    nr_threads = textToInteger(parser.getOption("--j", "Number of (openMP) threads", "1"));
    do_only_make_3d = parser.checkOption("--only_make_3d", "Take the iniref image, and create a 3D model from that without any alignment of the input images");

    verb = 1;

    if (parser.checkForErrors(verb))
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}


// Run multiple iterations


void HelixAligner::initialise()
{

	// Randomise the order of the particles
	if (random_seed == -1) random_seed = time(NULL);
    // Also randomize random-number-generator for perturbations on the angles
    init_random_generator(random_seed);

	if (fn_imgs!= "")
	{
		// Get the image size
		MetaDataTable MD;
		MD.read(fn_imgs);
		MD.firstObject();
		FileName fn_img;
		Image<RFLOAT> img;
		if (MD.containsLabel(EMDL_IMAGE_NAME))
			MD.getValue(EMDL_IMAGE_NAME, fn_img);
		else if (MD.containsLabel(EMDL_MLMODEL_REF_IMAGE))
			MD.getValue(EMDL_MLMODEL_REF_IMAGE, fn_img);
		else
			REPORT_ERROR("ERROR: input STAR file does not contain rlnImageName or rlnReferenceImage!");
		img.read(fn_img, false); // only read the header
		int xdim=XSIZE(img());
		int ydim = YSIZE(img());
		ori_size = xdim;
		if (XSIZE(img()) != YSIZE(img()) || ZSIZE(img()) != 1)
			REPORT_ERROR("ERROR: only squared 2D images are allowed.");

		// Get the pixel size
		if (MD.containsLabel(EMDL_CTF_MAGNIFICATION) && MD.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
		{
			RFLOAT mag, dstep, my_angpix;
			MD.getValue(EMDL_CTF_MAGNIFICATION, mag);
			MD.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
			my_angpix = 10000. * dstep / mag;
			std::cout << " Using pixel size from the input STAR file: " << my_angpix << std::endl;
			angpix = my_angpix;
		}

	}
	else if (fn_mics != "")
	{
		// Read in the micrographs STAR file
		MDmics.read(fn_mics,"micrographs");

		// Get the pixel size
		if (MDmics.containsLabel(EMDL_CTF_MAGNIFICATION) && MDmics.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
		{
			RFLOAT mag, dstep, my_angpix;
			MDmics.getValue(EMDL_CTF_MAGNIFICATION, mag);
			MDmics.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
			my_angpix = 10000. * dstep / mag;
			std::cout << " Using pixel size from the input STAR file: " << my_angpix << std::endl;
			angpix = my_angpix;
		}

		// Make sure the coordinate file directory names end with a '/'
		if (fn_coord_dir != "ASINPUT" && fn_coord_dir[fn_coord_dir.length()-1] != '/')
			fn_coord_dir+="/";

		// Loop over all micrographs in the input STAR file and warn of coordinate file or micrograph file do not exist
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmics)
		{
			FileName fn_mic;
			MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			FileName fn_pre, fn_jobnr, fn_post;
			decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
			FileName fn_coord = fn_coord_dir + fn_post.withoutExtension() + fn_coord_suffix;
			if (!exists(fn_coord))
				std::cerr << "Warning: coordinate file " << fn_coord << " does not exist..." << std::endl;
			if (!exists(fn_mic))
				std::cerr << "Warning: micrograph file " << fn_mic << " does not exist..." << std::endl;
		}

		ori_size = extract_width;
	}
	else if (do_only_make_3d && fn_inimodel != "")
	{
		Image<RFLOAT> img;
		img.read(fn_inimodel);
		img().setXmippOrigin();
		if (angpix < 0.)
		{
			img.MDMainHeader.getValue(EMDL_IMAGE_SAMPLINGRATE_X, angpix);
			std::cout << " Using pixel size from the input file header: " << angpix << std::endl;
		}
		ori_size = XSIZE(img());

		// The 3D reconstruction
		float deg_per_pixel = 180. * angpix / (crossover_distance);
		Image<RFLOAT> vol;
		vol().resize(ori_size, ori_size, ori_size);
		vol().setXmippOrigin();
		for (long int k=STARTINGZ(vol()); k<=FINISHINGZ(vol()); k++)
		{
			float ang = deg_per_pixel * k;
			Matrix2D<RFLOAT> Arot;
			rotation2DMatrix(ang, Arot);

			MultidimArray<RFLOAT> Mrot;
			Mrot.initZeros(img());
			Mrot.setXmippOrigin();
			applyGeometry(img(), Mrot, Arot, true, false);
			FOR_ALL_ELEMENTS_IN_ARRAY2D(Mrot)
				A3D_ELEM(vol(), k, i, j) = A2D_ELEM(Mrot, i, j);
		}
		vol.setSamplingRateInHeader(angpix);
		vol.write(fn_out + ".mrc");
		std::cout << " * Written " << fn_out << ".mrc" << std::endl;
		exit(RELION_EXIT_SUCCESS);
	}
	else
	{
		REPORT_ERROR("ERROR: provide --i, -mic, or --only_make_3d and --iniref");
	}

	if (angpix < 0.)
	{
		REPORT_ERROR("ERROR: provide pixel size through --angpix or through the magnification and detectorpixel size in the input STAR file.");
	}

	if (maxres < 0. || maxres < 2. * angpix)
	{
		maxres = 2. * angpix;
		std::cout << " Setting maximum resolution to " << maxres << std::endl;
	}

	down_size = ori_size * angpix * 2. / maxres;

	// Make sure that the crossover distance is close to an integer times the (downsized) pixel size of the model!
	float best_fit = 1.;
	down_angpix = 0.;
	int best_size = 0;
	for (int delta_size = -search_size; delta_size <= search_size; delta_size+= 2)
	{
		int mysize = down_size + delta_size;
		mysize -= mysize%2; //make even in case it is not already
		if (mysize <= ori_size)
		{
			float myangpix = angpix * (float)ori_size/(float)mysize;
			float mydiv = (2. * crossover_distance) / myangpix;
			// Also want even number of pixels in rectangle!
			float myfit = fmod(mydiv, 2);
			if (myfit > 1.)
				myfit -= 2.;
			myfit = fabs(myfit);
			if (myfit < best_fit)
			{
				best_fit = myfit;
				down_angpix = myangpix;
				best_size = mysize;
			}
			std::cout << " *   mydiv= " << mydiv << " myangpix= " << myangpix << " myfit= " << myfit << std::endl;
		}
	}
	std::cout <<     " *** best_angpix= " << down_angpix << " rectangles xsize= " << (2. * crossover_distance)/down_angpix << std::endl;

	down_size = best_size;
	yrect  = ROUND(ori_size * angpix/down_angpix);
	// Make even
	yrect -= yrect%2;
	xrect  = ROUND((2. * crossover_distance)/down_angpix);
	model.initialise(nr_classes, yrect, xrect);
	max_shift = CEIL(max_shift_A / down_angpix);
	mask_radius_pix = (mask_diameter > 0) ? CEIL(mask_diameter / (2. * down_angpix)) : yrect/2 - 2;
	std::cout << " maxres= " << maxres << " angpix= " << angpix << " down_size= " << down_size << std::endl;
	std::cout << " xrect= " << xrect << " yrect= " << yrect << " down_angpix= " << down_angpix << std::endl;
	std::cout << " max_shift= " << max_shift << " mask_radius_pix= "<< mask_radius_pix<< std::endl;

	// Now read in all images
	if (fn_mics == "")
		readImages();
	else
		getHelicesFromMics();

	initialiseClasses();

}

// Read in all the images
void HelixAligner::readImages()
{

	MD.read(fn_imgs);

	if (verb > 0)
	{
		std::cout << " Reading in all images ..." << std::endl;
		init_progress_bar(MD.numberOfObjects());
	}

	std::vector<MultidimArray<RFLOAT> > dummy;


	long int ipart=0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
	{

		FileName fn_img;
		Image<RFLOAT> img;
		if (MD.containsLabel(EMDL_IMAGE_NAME))
			MD.getValue(EMDL_IMAGE_NAME, fn_img);
		else if (MD.containsLabel(EMDL_MLMODEL_REF_IMAGE))
			MD.getValue(EMDL_MLMODEL_REF_IMAGE, fn_img);
		else
			REPORT_ERROR("ERROR: input STAR file does not contain rlnImageName or rlnReferenceImage!");
		img.read(fn_img);
		img().setXmippOrigin();
		// Rethink this when expanding program to 3D!
		RFLOAT yoff = 0.;
		RFLOAT psi = 0.;
		if (MD.containsLabel(EMDL_ORIENT_ORIGIN_Y_ANGSTROM))
			MD.getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff);
		if (MD.containsLabel(EMDL_ORIENT_PSI))
			MD.getValue(EMDL_ORIENT_PSI, psi);
		ori_psis.push_back(psi);
		ori_yoffs.push_back(yoff);
		// Apply the actual transformation
		Matrix2D<RFLOAT> A;
		rotation2DMatrix(psi, A);
		MAT_ELEM(A,1, 2) = -yoff / angpix;
		selfApplyGeometry(img(), A, IS_INV, DONT_WRAP);


		Xrects.push_back(dummy);

		// Calculate all rotated versions
		if (ipart==0) psis.clear();
		for (int iflip =0; iflip < 2; iflip++)
		{
			for (RFLOAT ang = 0; ang <= max_rotate; ang += step_rotate)
			{
				Matrix2D<RFLOAT> Arot;
				MultidimArray<RFLOAT> Irot;
				RFLOAT myang = (iflip == 1) ? ang + 180. : ang;
				Irot.initZeros(img());
				rotation2DMatrix(myang, Arot);
				applyGeometry(img(), Irot, Arot, true, false);
				resizeMap(Irot, down_size);
				Irot.setXmippOrigin();
				Xrects[Xrects.size()-1].push_back(Irot);
				if (ipart==0) psis.push_back(myang);
				if (ang > 0.)
				{
					// Also rotate in the opposite direction
					Irot.initZeros(img());
					applyGeometry(img(), Irot, Arot, false, false);
					resizeMap(Irot, down_size);
					Irot.setXmippOrigin();
					Xrects[Xrects.size()-1].push_back(Irot);
					if (ipart==0) psis.push_back(-myang);
				}
			}
		}

		ipart++;
		if (verb > 0 && ipart%50 == 0)
			progress_bar(ipart);

//#define DEBUG_READIMAGES
#ifdef DEBUG_READIMAGES
				FileName fnt;
				fnt.compose("helixnew", Xrects.size(),"spi",3);
				Image<RFLOAT> It;
				It()=Xrects[Xrects.size()-1][3];
				It.write(fnt);
#endif

	}

	if (verb > 0)
		progress_bar(MD.numberOfObjects());

#ifdef DEBUG
	std::cerr << "done readImages" << std::endl;
#endif

}

void HelixAligner::getHelicesFromMics()
{
	if (verb > 0)
	{
		std::cout << " Reading in all micrographs ..." << std::endl;
		init_progress_bar(MDmics.numberOfObjects());
	}

	// Loop over all micrographs in the input STAR file and warn of coordinate file or micrograph file do not exist
	long int imic = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmics)
	{
		imic++;
		FileName fn_mic;
		MDmics.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
		FileName fn_pre, fn_jobnr, fn_post;
		decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
		FileName fn_coord = fn_coord_dir + fn_post.withoutExtension() + fn_coord_suffix;
		if (!exists(fn_mic) || !exists(fn_coord))
		{
			if (!exists(fn_mic))
				std::cerr << "Warning: micrograph file " << fn_mic << " does not exist..." << std::endl;
			if (!exists(fn_coord))
				std::cerr << "Warning: coordinate file " << fn_coord << " does not exist..." << std::endl;
		}
		else
		{
			Image<RFLOAT> Imic;
			Imic.read(fn_mic);
			RFLOAT avg = Imic().computeAvg();

			// Read in the coordinate files
			MetaDataTable MDcoords;
			MDcoords.read(fn_coord);
			if (MDcoords.numberOfObjects()%2 == 1)
			{
				std::cerr << " ERROR: not an even number of entries in " << fn_coord << "! Skipping this micrograph... " << std::endl;
				continue;
			}

			// Get all start-end coordinate pairs
			std::vector<RFLOAT> x1_coord_list, y1_coord_list, x2_coord_list, y2_coord_list, pitch_list;
			RFLOAT xp, yp;
			int MDobj_id = 0;
			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDcoords)
		    {
		    	MDobj_id++;
		    	MDcoords.getValue(EMDL_IMAGE_COORD_X, xp);
		    	MDcoords.getValue(EMDL_IMAGE_COORD_Y, yp);
		    	if (MDobj_id % 2)
		    	{
		    		x1_coord_list.push_back(xp);
		    		y1_coord_list.push_back(yp);
		    	}
		    	else
		    	{
		    		x2_coord_list.push_back(xp);
		    		y2_coord_list.push_back(yp);
		    	}
		    }

			// Now extract the images: make all helices stand upright... Y becomes helical axis, X becomes helix width
			// For that we need to do interpolations...
			for (int ipair = 0; ipair < x1_coord_list.size(); ipair++)
			{
				std::vector<MultidimArray<RFLOAT> > dummy;
				Xrects.push_back(dummy);

				// Calculate all rotated versions
				int oldxsize, oldysize, oldsize;
				bool do_set_oldsize = true;
				for (RFLOAT ang = 0.; ang <= max_rotate; ang += step_rotate)
				{
					RFLOAT x1,x2,y1,y2,xcen,ycen;
					x1=x1_coord_list[ipair];
					x2=x2_coord_list[ipair];
					y1=y1_coord_list[ipair];
					y2=y2_coord_list[ipair];
					xcen = x1 + (x2-x1)/2;
					ycen = y1 + (y2-y1)/2;

					int xsize = FLOOR(sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1)));
					RFLOAT phi = RAD2DEG(atan(RFLOAT(y2-y1)/RFLOAT(x2-x1)));
					MultidimArray<RFLOAT> Ihelix;
					Ihelix.resize(extract_width, xsize);
					Ihelix.setXmippOrigin();

					int nrots = (ang > 0.) ? 2 : 1;
					for (int irot = 0; irot < nrots; irot++)
					{
						Matrix2D<RFLOAT> Arot(3,3);
						if (irot == 0)
							rotation2DMatrix(phi+ang, Arot);
						else
							rotation2DMatrix(phi-ang, Arot);
						Arot(0,2) = xcen;
						Arot(1,2) = ycen;

						int m1, n1, m2, n2;
						RFLOAT x, y, xp, yp;
						RFLOAT minxp, minyp, maxxp, maxyp;
						int cen_x, cen_y, cen_xp, cen_yp;
						RFLOAT wx, wy;
						int Xdim, Ydim;

						// Find center and limits of image
						cen_y  = (int)(YSIZE(Ihelix) / 2);
						cen_x  = (int)(XSIZE(Ihelix) / 2);
						cen_yp = (int)(YSIZE(Imic()) / 2);
						cen_xp = (int)(XSIZE(Imic()) / 2);
						minxp  = 0;
						minyp  = 0;
						maxxp  = XSIZE(Imic()) - 1;
						maxyp  = YSIZE(Imic()) - 1;
						Xdim   = XSIZE(Imic());
						Ydim   = YSIZE(Imic());

						for (int i = 0; i < YSIZE(Ihelix); i++)
						{
							// Calculate position of the beginning of the row in the output image
							x = -cen_x;
							y = i - cen_y;

							// Calculate this position in the input image according to the
							// geometrical transformation
							// they are related by
							// coords_output(=x,y) = A * coords_input (=xp,yp)
							xp = x * Arot(0, 0) + y * Arot(0, 1) + Arot(0, 2);
							yp = x * Arot(1, 0) + y * Arot(1, 1) + Arot(1, 2);

							for (int j = 0; j < XSIZE(Ihelix); j++)
							{
								bool interp;
								RFLOAT tmp;

								// If the point is outside the image, apply a periodic extension
								// of the image, what exits by one side enters by the other
								interp = true;
								if (xp < minxp  || xp > maxxp)
									interp = false;

								if (yp < minyp  || yp > maxyp)
									interp = false;

								if (interp)
								{
										// Linear interpolation

										// Calculate the integer position in input image, be careful
										// that it is not the nearest but the one at the top left corner
										// of the interpolation square. Ie, (0.7,0.7) would give (0,0)
										// Calculate also weights for point m1+1,n1+1
										wx = xp;// + cen_xp;
										m1 = (int) wx;
										wx = wx - m1;
										m2 = m1 + 1;
										wy = yp;// + cen_yp;
										n1 = (int) wy;
										wy = wy - n1;
										n2 = n1 + 1;

										// Perform interpolation
										// if wx == 0 means that the rightest point is useless for this
										// interpolation, and even it might not be defined if m1=xdim-1
										// The same can be said for wy.
										tmp  = (RFLOAT)((1 - wy) * (1 - wx) * DIRECT_A2D_ELEM(Imic(), n1, m1));

										if (m2 < Xdim)
											tmp += (RFLOAT)((1 - wy) * wx * DIRECT_A2D_ELEM(Imic(), n1, m2));

										if (n2 < Ydim)
										{
											tmp += (RFLOAT)(wy * (1 - wx) * DIRECT_A2D_ELEM(Imic(), n2, m1));

											if (m2 < Xdim)
												tmp += (RFLOAT)(wy * wx * DIRECT_A2D_ELEM(Imic(), n2, m2));
										}

										dAij(Ihelix, i, j) = tmp;

								} // if interp
								else
									dAij(Ihelix, i, j) = avg;

								// Compute new point inside input image
								xp += Arot(0, 0);
								yp += Arot(1, 0);
							}
						}
//#define DEBUG_GETHELICESFROMMICS
#ifdef DEBUG_GETHELICESFROMMICS
				FileName fntt;
				fntt.compose("helixnew1_beforedown", Xrects.size(),"spi",3);
				Image<RFLOAT> Itt;
				Itt()=Ihelix;
				Itt.write(fntt);
#endif

						// Downscale if needed
						MultidimArray<RFLOAT> Idown = Ihelix;
						if (down_angpix > angpix)
						{
							RFLOAT avg = Idown.computeAvg();

							int oldxsize = XSIZE(Idown);
							int oldysize = YSIZE(Idown);
							int oldsize = oldxsize;

							if ( oldxsize != oldysize )
							{
								oldsize = XMIPP_MAX( oldxsize, oldysize );
								Idown.setXmippOrigin();
								Idown.window(FIRST_XMIPP_INDEX(oldsize), FIRST_XMIPP_INDEX(oldsize),
											  LAST_XMIPP_INDEX(oldsize),  LAST_XMIPP_INDEX(oldsize), avg);
							}
							int newsize = ROUND(oldsize * angpix / down_angpix);
							newsize -= newsize%2; //make even in case it is not already
							resizeMap(Idown, newsize);

							if ( oldxsize != oldysize )
							{
								int newxsize = ROUND(oldxsize * angpix / down_angpix);
								int newysize = ROUND(oldysize * angpix / down_angpix);
								newxsize -= newxsize%2; //make even in case it is not already
								newysize -= newysize%2; //make even in case it is not already
								Idown.setXmippOrigin();
								Idown.window(FIRST_XMIPP_INDEX(newysize), FIRST_XMIPP_INDEX(newxsize),
											  LAST_XMIPP_INDEX(newysize),  LAST_XMIPP_INDEX(newxsize));
							}

						}

						// adhoc normalisation of images
						RFLOAT avg,stddev,min,max;
						Idown.computeStats(avg,stddev,min,max);
						Idown -= avg;
						Idown /= -stddev; // divide by minus stddev to flip contrast to white...

						Xrects[Xrects.size()-1].push_back(Idown);
					} // end for irot (for positive and negative rotation
				} // end for over rotations

				if (verb > 0)
					progress_bar(imic);


//#define DEBUG_GETHELICESFROMMICS2
#ifdef DEBUG_GETHELICESFROMMICS2
				FileName fnt;
				fnt.compose("helixnew1", Xrects.size(),"spi",3);
				Image<RFLOAT> It;
				It()=Xrects[Xrects.size()-1][1];
				It.write(fnt);
#endif
			}
		}
	}
	if (verb > 0)
		progress_bar(MDmics.numberOfObjects());

}

void HelixAligner::initialiseClasses()
{
#ifdef DEBUG
	std::cerr << "Entering initialiseClasses" << std::endl;
#endif
	if (model.Aref.size() == 0)
		REPORT_ERROR("BUG: non-initialised model!");

	if (verb > 0)
		std::cout << " Initialising reference(s) ..." << std::endl;

	if (fn_inimodel != "")
	{

		if (nr_classes > 1)
			REPORT_ERROR("ERROR: can only use initial reference for single-class!");

		Image<RFLOAT> Iref;
		Iref.read(fn_inimodel);
		resizeMap(Iref(), YSIZE(model.Aref[0]));
		Iref().setXmippOrigin();
		std::cerr << " model.Arec.size()= " << model.Arec.size() << std::endl;
		model.Arec[0] = Iref();
	    // Now project the reconstruction back out into the model.Aref[iclass]
	    Projector PP(YSIZE(model.Aref[0]), TRILINEAR, 2, 1, 1);
	    // Set the FT of img inside the Projector
	    MultidimArray<RFLOAT> dummy;
	    PP.computeFourierTransformMap(Iref(), dummy, YSIZE(model.Aref[0]), 1);

	    // Calculate all projected lines
	    for (int j = 0; j < XSIZE(model.Aref[0]); j++)
	    {
	    	Matrix2D<RFLOAT> A2D;
			MultidimArray<RFLOAT> myline(YSIZE(model.Aref[0]));
			MultidimArray<Complex> myFline(YSIZE(model.Aref[0])/2 + 1);
			FourierTransformer transformer;

			RFLOAT rot = (RFLOAT)j*360./(XSIZE(model.Aref[0]));
	    	rotation2DMatrix(rot, A2D);
	    	PP.get2DFourierTransform(myFline, A2D);
	    	transformer.inverseFourierTransform(myFline,myline);
	    	// Shift the image back to the center...
	    	myline.setXmippOrigin();
	    	CenterFFT(myline, false);
			for (int i = 0; i < YSIZE(model.Aref[0]); i++)
				DIRECT_A2D_ELEM(model.Aref[0], i, j) = DIRECT_A1D_ELEM(myline, i);
	    }

#define DEBUGREC2D
#ifdef DEBUGREC2D
	    Image<RFLOAT> It;
	    It()=model.Aref[0];
	    It.write("after_reproject.spi");
#endif

	}
	else
	{
		// Randomly position all particles along the X-direction
		model.initZeroSums();
		// Loop over all particles
		if (verb > 0)
			init_progress_bar(Xrects.size());

		for (int ipart = 0; ipart < Xrects.size(); ipart++)
		{
			// Set into a random class
			int myclass = (int)(rnd_unif() * nr_classes);
			int random_xoffset = (int)(rnd_unif() * xrect);
			for (int j_smear = -max_smear; j_smear <= max_smear; j_smear++)
			{

				double smearw = (max_smear ==0 ) ? 1 : gaussian1D((double)j_smear, (double)max_smear/3);
				FOR_ALL_ELEMENTS_IN_ARRAY2D(Xrects[ipart][0])
				{
					int jp = j + random_xoffset + j_smear;
					while (jp < STARTINGX(model.Aref[myclass]))
						jp += xrect;
					while (jp > FINISHINGX(model.Aref[myclass]))
						jp -= xrect;

					// this places the original image in the offset-translated center of the rectangle
					A2D_ELEM(model.Asum[myclass], i, jp) += smearw * A2D_ELEM(Xrects[ipart][0], i, j);
					A2D_ELEM(model.Asumw[myclass], i, jp) += smearw;

					// This places the Y-flipped image at half a cross-over distance from the first one
					int ip = -i;
					if (ip >= STARTINGY(Xrects[ipart][0]) && ip <= FINISHINGY(Xrects[ipart][0]))
					{
						int jjp = jp + xrect/2;
						while (jjp < STARTINGX(model.Aref[myclass]))
							jjp += xrect;
						while (jjp > FINISHINGX(model.Aref[myclass]))
							jjp -= xrect;
						A2D_ELEM(model.Asum[myclass], ip, jjp) += smearw * A2D_ELEM(Xrects[ipart][0], i, j);
						A2D_ELEM(model.Asumw[myclass], ip, jjp) += smearw;

					}
				}
			}
			model.pdf[myclass] += 1.;
			if (verb > 0)
				progress_bar(ipart);
		}

		if (verb > 0)
			progress_bar(Xrects.size());

		// After all images have been set, maximise the references in the model
		maximisation();
	}
#ifdef DEBUG
	std::cerr << "Leaving initialiseClasses" << std::endl;
#endif
}

void HelixAligner::expectationOneParticleNoFFT(long int ipart)
{

	int twostarty = 2 * FIRST_XMIPP_INDEX(yrect);
	double maxccf = -100.;
	int best_class = -1;
	int best_k_rot = -1;
	int best_i_offset = -1;
	int best_j_offset = -1;
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{

		int k_rot = 0;
		for (int k_rot = 0; k_rot < Xrects[ipart].size(); k_rot++)
		{
			for (int i_offset = -max_shift; i_offset <= max_shift; i_offset++)
			{
				for (int j_offset = 0; j_offset < xrect; j_offset++)
				{
					double ccf_xa = 0;
					double ccf_x2 = 0;
					double ccf_a2 = 0;
				    for (long int i=STARTINGY(Xrects[ipart][k_rot]); i<=FINISHINGY(Xrects[ipart][k_rot]); i++) \
				    {

						int ip = i + i_offset;
						if (ip < -mask_radius_pix || ip > mask_radius_pix)
							continue;

						/*
						while (ip < STARTINGY(model.Aref[iclass]))
							ip += yrect;
						while (ip > FINISHINGY(model.Aref[iclass]))
							ip -= yrect;
						*/

						for (long int j=STARTINGX(Xrects[ipart][k_rot]); j<=FINISHINGX(Xrects[ipart][k_rot]); j++)
				    	{
							int jp = j + j_offset;
							while (jp < STARTINGX(model.Aref[iclass]))
								jp += xrect;
							while (jp > FINISHINGX(model.Aref[iclass]))
								jp -= xrect;

							// This places the Y-flipped image at half a cross-over distance from the first one
							int ipp = -ip;
							// Don't let the image run out of the height of the box
							if (ipp >= STARTINGY(Xrects[ipart][k_rot]) && ipp <= FINISHINGY(Xrects[ipart][k_rot]))
							{
								int jpp = jp + xrect/2;
								while (jpp < STARTINGX(model.Aref[iclass]))
									jpp += xrect;
								while (jpp > FINISHINGX(model.Aref[iclass]))
									jpp -= xrect;

								// this places the original image in the offset-translated center of the rectangle
								ccf_xa += A2D_ELEM(model.Aref[iclass], ip, jp) * A2D_ELEM(Xrects[ipart][k_rot], i, j);
								ccf_a2 += A2D_ELEM(model.Aref[iclass], ip, jp) * A2D_ELEM(model.Aref[iclass], ip, jp);

								ccf_xa += A2D_ELEM(model.Aref[iclass], ipp, jpp) * A2D_ELEM(Xrects[ipart][k_rot], i, j);
								ccf_a2 += A2D_ELEM(model.Aref[iclass], ipp, jpp) * A2D_ELEM(model.Aref[iclass], ipp, jpp);

								ccf_x2 += 2. * A2D_ELEM(Xrects[ipart][k_rot], i, j) * A2D_ELEM(Xrects[ipart][k_rot], i, j);

							}
						} // end loop j
				    } // end loop i

					double ccf = (ccf_x2 > 0. && ccf_a2 > 0.) ? ccf_xa/(sqrt(ccf_x2) * sqrt(ccf_a2)) : 0.;

					// Find the best fit
					if (ccf > maxccf)
					{
						maxccf = ccf;
						best_class = iclass;
						best_k_rot = k_rot;
						best_i_offset = i_offset;
						best_j_offset = j_offset;
					}
				} // end for j_offset
			} // end for i_offset
		} // end for k_rot
	}


	if (maxccf < -1.)
		REPORT_ERROR("BUG: not found maxccf!");

	// Now set the optimal Y-translations and rotations in the output STAR file
	RFLOAT yoff, psi;
	psi = ori_psis[ipart] + psis[best_k_rot];
	yoff = ori_yoffs[ipart] + best_i_offset * down_angpix;

	MD.setValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, yoff, ipart);
	MD.setValue(EMDL_ORIENT_PSI, psi, ipart);

	// Now add the image to that class reference
	// To ensure continuity in the reference: smear out every image along X
	#pragma omp critical
	{
		for (int j_smear = -max_smear; j_smear <= max_smear; j_smear++)
		{
			double smearw = (max_smear< XMIPP_EQUAL_ACCURACY) ? 1 : gaussian1D((double)j_smear, (double)max_smear/3);
			FOR_ALL_ELEMENTS_IN_ARRAY2D(Xrects[ipart][best_k_rot])
			{

				int jp = j + best_j_offset + j_smear;
				while (jp < STARTINGX(model.Aref[best_class]))
					jp += xrect;
				while (jp > FINISHINGX(model.Aref[best_class]))
					jp -= xrect;

				int ip = i + best_i_offset;
				while (ip < STARTINGY(model.Aref[best_class]))
					ip += yrect;
				while (ip > FINISHINGY(model.Aref[best_class]))
					ip -= yrect;

				// this places the original image in the offset-translated center of the rectangle
				A2D_ELEM(model.Asum[best_class], ip, jp) += smearw * A2D_ELEM(Xrects[ipart][best_k_rot], i, j);
				A2D_ELEM(model.Asumw[best_class], ip, jp) += smearw;

				// This places the Y-flipped image at half a cross-over distance from the first one
				int ipp = -ip;
				if (ipp >= STARTINGY(Xrects[ipart][best_k_rot]) && ipp <= FINISHINGY(Xrects[ipart][best_k_rot]))
				{
					int jpp = jp + xrect/2;
					while (jpp > FINISHINGX(model.Aref[best_class]))
						jpp -= xrect;
					A2D_ELEM(model.Asum[best_class], ipp, jpp) += smearw * A2D_ELEM(Xrects[ipart][best_k_rot], i, j);
					A2D_ELEM(model.Asumw[best_class], ipp, jpp) += smearw;
				}
			}
		}
		model.pdf[best_class] += 1.;
	}

}

void HelixAligner::expectation()
{

	// Initialise the wsum_model to zeros
	model.initZeroSums();


	if (verb > 0)
	{
		init_progress_bar(Xrects.size());
	}

	#pragma omp parallel for num_threads(nr_threads)
	for (long int ipart = 0; ipart < Xrects.size(); ipart++)
	{

		expectationOneParticleNoFFT(ipart);

		if (ipart%nr_threads==0)
			progress_bar(ipart);
	}

	progress_bar(Xrects.size());

}

void HelixAligner::maximisation()
{
#ifdef DEBUGREC2D
	Image<RFLOAT> It;
	It()=model.Asumw[0];
	It.write("Asumw.spi");
	It()=model.Asum[0];
	It.write("Asum.spi");
#endif

	// Update the references
	double allsum = 0.;
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		for (int i = 0; i < yrect; i++)
		{
			for (int j = 0; j < xrect; j++)
			{
				if (DIRECT_A2D_ELEM(model.Asumw[iclass], i, j) > 0.)
					DIRECT_A2D_ELEM(model.Aref[iclass], i, j) =  DIRECT_A2D_ELEM(model.Asum[iclass], i, j) / DIRECT_A2D_ELEM(model.Asumw[iclass], i, j);
				else
					DIRECT_A2D_ELEM(model.Aref[iclass], i, j) = 0.;

				// Also store  sum of classes in Asum for writeOut
				DIRECT_A2D_ELEM(model.Asum[iclass], i, j)  = DIRECT_A2D_ELEM(model.Aref[iclass], i, j);
			}
		}
		allsum += model.pdf[iclass];

		reconstruct2D(iclass);
	}
	for (int iclass = 0; iclass < nr_classes; iclass++)
		model.pdf[iclass]/=allsum;

}

void HelixAligner::reconstruct2D(int iclass)
{
#ifdef DEBUG
	std::cerr << "Entering reconstruct2D" << std::endl;
#endif
#ifdef DEBUGREC2D
	Image<RFLOAT> It;
	It()=model.Aref[iclass];
	It.write("before_reproject.spi");
#endif

	// Loop over the length of the helix to get the transforms of all 1D images
	std::vector<MultidimArray<Complex> > myFlines;
	for (int j = 0; j < XSIZE(model.Aref[iclass]); j++)
	{
		MultidimArray<RFLOAT> myline(YSIZE(model.Aref[iclass]));
		MultidimArray<Complex> myFline;
		FourierTransformer transformer;

		for (int i = 0; i < YSIZE(model.Aref[iclass]); i++)
			DIRECT_A1D_ELEM(myline, i) = DIRECT_A2D_ELEM(model.Aref[iclass], i, j);

		CenterFFT(myline, true);
		transformer.FourierTransform(myline, myFline, false);
		myFlines.push_back(myFline);
	}

	// Then reconstruct
    BackProjector BP(YSIZE(model.Aref[iclass]), 2, "C1", TRILINEAR, 2, 1, 0, 1.9, 15, 1, false);
	BP.initialiseDataAndWeight(YSIZE(model.Aref[iclass]));

    for (int j = 0; j < myFlines.size(); j++)
    {
    	Matrix2D<RFLOAT> A2D;
    	RFLOAT rot = (RFLOAT)j*360./(XSIZE(model.Aref[iclass]));
    	rotation2DMatrix(rot, A2D);
    	BP.set2DFourierTransform(myFlines[j], A2D);
    }
    MultidimArray<RFLOAT> dummy;
    model.Arec[iclass].initZeros();
    BP.reconstruct(model.Arec[iclass], 10, false, dummy);

    if (symmetry > 1)
    {

#ifdef DEBUGREC2D
		It()=model.Arec[iclass];
		resizeMap(It(), ori_size);
		It.write("rec_beforesym.spi");
#endif
		MultidimArray<RFLOAT> Asum = model.Arec[iclass];
		for (int i = 1; i < symmetry; i++)
		{
			RFLOAT ang = i*360./(RFLOAT)symmetry;
		   	Matrix2D<RFLOAT> A2D;
	    	rotation2DMatrix(ang, A2D);
	    	MultidimArray<RFLOAT> Arot;
	    	Arot.initZeros(model.Arec[iclass]);
	    	applyGeometry(model.Arec[iclass], Arot, A2D, false, false);
	    	Asum += Arot;
		}
		model.Arec[iclass] = Asum / (RFLOAT)symmetry;
    }

    if (mask_diameter > 0.)
    {
    	RFLOAT pixel_radius = mask_diameter/(2.*down_angpix);
    	softMaskOutsideMap(model.Arec[iclass], pixel_radius, 0.);
    }


#ifdef DEBUGREC2D
	It()=model.Arec[iclass];
	resizeMap(It(), ori_size);
	It.write("rec.spi");
#endif

    // Now project the reconstruction back out into the model.Aref[iclass]
    Projector PP(YSIZE(model.Aref[iclass]), TRILINEAR, 2, 1, 1);
    // Set the FT of img inside the Projector
    PP.computeFourierTransformMap(model.Arec[iclass], dummy, YSIZE(model.Aref[iclass]), 1);

    // Calculate all projected lines
    for (int j = 0; j < myFlines.size(); j++)
    {
    	Matrix2D<RFLOAT> A2D;
		MultidimArray<RFLOAT> myline(YSIZE(model.Aref[iclass]));
		FourierTransformer transformer;

		RFLOAT rot = (RFLOAT)j*360./(XSIZE(model.Aref[iclass]));
    	rotation2DMatrix(rot, A2D);
    	myFlines[j].initZeros();
    	PP.get2DFourierTransform(myFlines[j], A2D);
    	transformer.inverseFourierTransform(myFlines[j],myline);
    	// Shift the image back to the center...
    	CenterFFT(myline, false);

		for (int i = 0; i < YSIZE(model.Aref[iclass]); i++)
    		 DIRECT_A2D_ELEM(model.Aref[iclass], i, j) = DIRECT_A1D_ELEM(myline, i);

    }
#ifdef DEBUGREC2D
	It()=model.Aref[iclass];
	It.write("after_reproject.spi");
#endif
#ifdef DEBUG
	std::cerr << "Leaving reconstruct2D" << std::endl;
#endif

}

void HelixAligner::writeOut(int iter)
{

	//std::cout << " **** Model for iteration " << iter << std::endl;

#ifdef DEBUG
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		FileName fn_class = fn_out + "_it" + integerToString(iter, 3) + "_class" + integerToString(iclass+1, 3)+".spi";
		Image<RFLOAT> Ic;
		Ic()=model.Aref[iclass];
		Ic.write(fn_class);
		std::cout << " * Written " << fn_class << std::endl;
		fn_class = fn_out + "_it" + integerToString(iter, 3) + "_class" + integerToString(iclass+1, 3)+"_reconstructed.spi";
		Ic()=model.Arec[iclass];
		Ic.write(fn_class);
		std::cout << " * Written " << fn_class << std::endl;
	}
#else
	FileName fn_iter = fn_out + "_it" + integerToString(iter, 3);
	MD.write(fn_iter+".star");
	Image<RFLOAT> Aimg(xrect, yrect, 1, nr_classes);

	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(model.Aref[iclass])
		{
			DIRECT_NZYX_ELEM(Aimg(), iclass, 0, i, j) = DIRECT_A2D_ELEM(model.Aref[iclass], i, j);
		}
    }
	Aimg.write(fn_iter + "_reprojections.mrcs");

	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(model.Asum[iclass])
		{
			DIRECT_NZYX_ELEM(Aimg(), iclass, 0, i, j) = DIRECT_A2D_ELEM(model.Asum[iclass], i, j);
		}
    }
	Aimg.write(fn_iter + "_summed_classes.mrcs");

	Image<RFLOAT> Aimg2(yrect, yrect, 1, nr_classes);
	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(model.Arec[iclass])
		{
			DIRECT_NZYX_ELEM(Aimg2(), iclass, 0, i, j) = DIRECT_A2D_ELEM(model.Arec[iclass], i, j);
		}
	}
	Aimg2.write(fn_iter + "_reconstructed.mrcs");
#endif

	if (nr_classes > 1)
	{
		for (int iclass = 0; iclass < nr_classes; iclass++)
		{
			std:: cout << " * Fraction class " << iclass+1 << " = " << model.pdf[iclass] << std::endl;
		}
	}

}

void HelixAligner::reconstruct3D()
{

	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		FileName fn_class = fn_out + "_class" + integerToString(iclass+1, 3)+"_projections.spi";
		Image<RFLOAT> Ic;
		Ic()=model.Aref[iclass];
		Ic.setSamplingRateInHeader(angpix);
		Ic.write(fn_class);
		std::cout << " * Written " << fn_class << std::endl;
		fn_class = fn_out + "_class" + integerToString(iclass+1, 3)+"_rec2d.spi";
		Ic()=model.Arec[iclass];
		resizeMap(Ic(), ori_size);
		Ic.setSamplingRateInHeader(angpix);
		Ic.write(fn_class);
		MultidimArray<RFLOAT> Mori = Ic();
		std::cout << " * Written " << fn_class << std::endl;

		// The 3D reconstruction
		float deg_per_pixel = 180. * angpix / (crossover_distance);
		Ic().resize(ori_size, ori_size, ori_size);
		for (int k = 0; k < ZSIZE(Ic()); k++)
		{
			float ang = deg_per_pixel * k;
			Matrix2D<RFLOAT> Arot;
			rotation2DMatrix(ang, Arot);

			MultidimArray<RFLOAT> Mrot;
			Mrot.initZeros(Mori);
			applyGeometry(Mori, Mrot, Arot, true, false);
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Mrot)
				DIRECT_A3D_ELEM(Ic(), k, i, j) = DIRECT_A2D_ELEM(Mrot, i, j);
		}
		fn_class = fn_out + "_class" + integerToString(iclass+1, 3)+"_rec3d.mrc";
		Ic.setSamplingRateInHeader(angpix);
		Ic.write(fn_class);
		std::cout << " * Written " << fn_class << std::endl;
	}

}

// Run multiple iterations
void HelixAligner::run()
{

	// Write out the starting model as well
	writeOut(0);

	int decrease_smear = ROUND((float)max_smear/(float)(nr_iter+5));
	for (int iter = 1; iter <= nr_iter; iter++)
	{

		if (verb > 0)
		{
			std::cout << " Iteration " << iter <<" of " << nr_iter << std::endl;
			if (max_smear > 0)
				std::cout << "  = smearing references by " << max_smear << " downsampled pixels along helical axis " << std::endl;
		}

		expectation();

		maximisation();

		writeOut(iter);

		if (max_smear > 0)
			max_smear -= decrease_smear;

	}

	// Reconstruct the final solution in 3D
	reconstruct3D();


}


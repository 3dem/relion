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
#include "src/autopicker.h"

//#define DEBUG
//#define DEBUG_HELIX

void ccfPeak::clear()
{
	id = ref = nr_peak_pixel = -1;
	x = y = r = area_percentage = fom_max = psi = dist = fom_thres = (-1.);
	ccf_pixel_list.clear();
}

bool ccfPeak::isValid() const
{
	// Invalid parameters
	if ( (r < 0.) || (area_percentage < 0.) || (ccf_pixel_list.size() < 1) )
		return false;
	// TODO: check ccf values in ccf_pixel_list?
	for (int id = 0; id < ccf_pixel_list.size(); id++)
	{
		if (ccf_pixel_list[id].fom > fom_thres)
			return true;
	}
	return false;
}

bool ccfPeak::operator<(const ccfPeak& b) const
{
	if (fabs(r - b.r) < 0.01)
		return (fom_max < b.fom_max);
	return (r < b.r);
}

bool ccfPeak::refresh()
{
	RFLOAT x_avg, y_avg;
	int nr_valid_pixel;

	area_percentage = (-1.);

	if (ccf_pixel_list.size() < 1)
		return false;

	fom_max = (-99.e99);
	nr_valid_pixel = 0;
	x_avg = y_avg = 0.;
	for (int id = 0; id < ccf_pixel_list.size(); id++)
	{
		if (ccf_pixel_list[id].fom > fom_thres)
		{
			nr_valid_pixel++;

			if (ccf_pixel_list[id].fom > fom_max)
				fom_max = ccf_pixel_list[id].fom;

			x_avg += ccf_pixel_list[id].x;
			y_avg += ccf_pixel_list[id].y;
		}
	}
	nr_peak_pixel = nr_valid_pixel;

	if (nr_valid_pixel < 1)
		return false;

	x = x_avg / (RFLOAT)(nr_valid_pixel);
	y = y_avg / (RFLOAT)(nr_valid_pixel);
	area_percentage = (RFLOAT)(nr_valid_pixel) / ccf_pixel_list.size();

	return true;
};

void AutoPicker::read(int argc, char **argv)
{

	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "Micrograph STAR file OR filenames from which to autopick particles, e.g. \"Micrographs/*.mrc\"");
	fn_out = parser.getOption("--pickname", "Rootname for coordinate STAR files", "autopick");
	fn_odir = parser.getOption("--odir", "Output directory for coordinate files (default is to store next to micrographs)", "AutoPick/");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size of the micrographs in Angstroms", "1"));
	particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms, default=automatic)", "-1"));
	decrease_radius = textToInteger(parser.getOption("--shrink_particle_mask", "Shrink the particle mask by this many pixels (to detect Einstein-from-noise classes)", "2"));
	outlier_removal_zscore= textToFloat(parser.getOption("--outlier_removal_zscore", "Remove pixels that are this many sigma away from the mean", "8."));
	do_write_fom_maps = parser.checkOption("--write_fom_maps", "Write calculated probability-ratio maps to disc (for re-reading in subsequent runs)");
	no_fom_limit = parser.checkOption("--no_fom_limit", "Ignore default maximum limit of 30 fom maps being written","false");
	do_read_fom_maps = parser.checkOption("--read_fom_maps", "Skip probability calculations, re-read precalculated maps from disc");
	do_optimise_scale = !parser.checkOption("--skip_optimise_scale", "Skip the optimisation of the micrograph scale for better prime factors in the FFTs. This runs slower, but at exactly the requested resolution.");
	do_only_unfinished = parser.checkOption("--only_do_unfinished", "Only autopick those micrographs for which the coordinate file does not yet exist");
	do_gpu = parser.checkOption("--gpu", "Use GPU acceleration when availiable");
	gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread","default");
#ifndef CUDA
	if(do_gpu)
	{
		std::cerr << "+ WARNING : Relion was compiled without CUDA of at least version 7.0 - you do NOT have support for GPUs" << std::endl;
		do_gpu = false;
	}
#endif
	int ref_section = parser.addSection("References options");
	fn_ref = parser.getOption("--ref", "STAR file with the reference names, or an MRC stack with all references, or \"gauss\" for blob-picking");
	angpix_ref = textToFloat(parser.getOption("--angpix_ref", "Pixel size of the references in Angstroms (default is same as micrographs)", "-1"));
	do_invert = parser.checkOption("--invert", "Density in micrograph is inverted w.r.t. density in template");
	psi_sampling = textToFloat(parser.getOption("--ang", "Angular sampling (in degrees); use 360 for no rotations", "10"));
	lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter in Angstroms for the references (prevent Einstein-from-noise!)","-1"));
	highpass = textToFloat(parser.getOption("--highpass", "Highpass filter in Angstroms for the micrographs","-1"));
	do_ctf = parser.checkOption("--ctf", "Perform CTF correction on the references?");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
	gauss_max_value = textToFloat(parser.getOption("--gauss_max", "Value of the peak in the Gaussian blob reference","0.1"));

	int helix_section = parser.addSection("Helix options");
	autopick_helical_segments = parser.checkOption("--helix", "Are the references 2D helical segments? If so, in-plane rotation angles (psi) are estimated for the references.");
	helical_tube_curvature_factor_max = textToFloat(parser.getOption("--helical_tube_kappa_max", "Factor of maximum curvature relative to that of a circle", "0.25"));
	helical_tube_diameter = textToFloat(parser.getOption("--helical_tube_outer_diameter", "Tube diameter in Angstroms", "-1"));
	helical_tube_length_min = textToFloat(parser.getOption("--helical_tube_length_min", "Minimum tube length in Angstroms", "-1"));

	int peak_section = parser.addSection("Peak-search options");
	min_fraction_expected_Pratio = textToFloat(parser.getOption("--threshold", "Fraction of expected probability ratio in order to consider peaks?", "0.25"));
	min_particle_distance = textToFloat(parser.getOption("--min_distance", "Minimum distance (in A) between any two particles (default is half the box size)","-1"));
	max_stddev_noise = textToFloat(parser.getOption("--max_stddev_noise", "Maximum standard deviation in the noise area to use for picking peaks (default is no maximum)","-1"));
	autopick_skip_side = textToInteger(parser.getOption("--skip_side", "Keep this many extra pixels (apart from particle_size/2) away from the edge of the micrograph ","0"));

	int expert_section = parser.addSection("Expert options");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));
	random_seed = textToInteger(parser.getOption("--random_seed", "Number for the random seed generator", "1"));
	workFrac = textToFloat(parser.getOption("--shrink", "Reduce micrograph to this fraction size, during correlation calc (saves memory and time)", "1.0"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	if (autopick_helical_segments)
	{
		if ( (helical_tube_curvature_factor_max < 0.0001) || (helical_tube_curvature_factor_max > 1.0001) )
			REPORT_ERROR("Error: Maximum curvature factor should be 0~1!");
		if (!(min_particle_distance > 0.))
			REPORT_ERROR("Error: Helical rise and the number of asymmetrical units between neighbouring helical segments should be positive!");
	}
}

void AutoPicker::usage()
{
	parser.writeUsage(std::cout);
}

void AutoPicker::initialise()
{

#ifdef TIMING
	TIMING_A0  =           timer.setNew("Initialise()");
	TIMING_A1  =           timer.setNew("--Init");
	TIMING_A2  =           timer.setNew("--Read Reference(s)");
	TIMING_A3  =           timer.setNew("--Read Micrograph(s)");
	TIMING_A4  =           timer.setNew("--Prep projectors");
	TIMING_A5  =           timer.setNew("autoPickOneMicrograph()");
	TIMING_A6  =           timer.setNew("--Read Micrographs(s)");
	TIMING_A7  =           timer.setNew("--Micrograph computestats");
	TIMING_A8  =           timer.setNew("--CTF-correct micrograph");
	TIMING_A9  =           timer.setNew("--Resize CCF and PSI-maps");
	TIMING_B1  =           timer.setNew("--FOM prep");
	TIMING_B2  =           timer.setNew("--Read reference(s) via FOM");
	TIMING_B3  =           timer.setNew("--Psi-dep correlation calc");
	TIMING_B4  =           timer.setNew("----ctf-correction");
	TIMING_B5  =           timer.setNew("----first psi");
	TIMING_B6  =           timer.setNew("----rest of psis");
	TIMING_B7  =           timer.setNew("----write fom maps");
	TIMING_B8  =           timer.setNew("----peak-prune/-search");
	TIMING_B9  =           timer.setNew("--final peak-prune");
#endif

#ifdef TIMING
		timer.tic(TIMING_A0);
		timer.tic(TIMING_A1);
#endif

	if (random_seed == -1) random_seed = time(NULL);

	if (fn_in.isStarFile())
	{
		MDmic.read(fn_in);
		fn_micrographs.clear();
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmic)
		{
			FileName fn_mic;
			MDmic.getValue(EMDL_MICROGRAPH_NAME, fn_mic);
			fn_micrographs.push_back(fn_mic);
		}

        RFLOAT mag, dstep;
        if (MDmic.containsLabel(EMDL_CTF_MAGNIFICATION) && MDmic.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
        {
			MDmic.goToObject(0);
        	MDmic.getValue(EMDL_CTF_MAGNIFICATION, mag);
			MDmic.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
			angpix = 10000. * dstep / mag;
			if (verb > 0)
				std::cout << " + Using (micrograph) pixel size from input STAR file of " << angpix << " Angstroms" << std::endl;
        }
        else if (verb > 0 )
        {
        	std::cout << " + Warning: input (micrograph) STAR file does not contain information about pixel size!" << std::endl;
        	std::cout << " + Warning: use --angpix_ref to provide the correct value. Now using " << angpix << " Angstroms" << std::endl;
        }
	}
	else
	{

		if (do_ctf)
			REPORT_ERROR("AutoPicker::initialise ERROR: use an input STAR file with the CTF information when using --ctf");

		fn_in.globFiles(fn_micrographs);
		if (fn_micrographs.size() == 0)
			REPORT_ERROR("Cannot find any micrograph called: "+fns_autopick);
	}

	// If we're continuing an old run, see which micrographs have not been finished yet...
	if (do_only_unfinished)
	{
		if (verb > 0)
			std::cout << " + Skipping those micrographs for which coordinate file already exists" << std::endl;
		std::vector<FileName> fns_todo;
		for (long int imic = 0; imic < fn_micrographs.size(); imic++)
		{

			FileName fn_tmp = getOutputRootName(fn_micrographs[imic]) + "_" + fn_out + ".star";
			if (!exists(fn_tmp))
				fns_todo.push_back(fn_micrographs[imic]);
		}

		fn_micrographs = fns_todo;
	}

	// If there is nothing to do, then go out of initialise
	todo_anything = true;
	if (fn_micrographs.size() == 0)
	{
		if (verb > 0)
			std::cout << " + No new micrographs to do, so exiting autopicking ..." << std::endl;
		todo_anything = false;
		return;
	}

	if (verb > 0)
	{
		if((fn_micrographs.size()>30 && do_write_fom_maps) && !no_fom_limit)
		{
			REPORT_ERROR("\n If you really want to write this many (" + integerToString(fn_micrographs.size()) + ") FOM-maps, add --no_fom_limit");
		}
		std::cout << " + Run autopicking on the following micrographs: " << std::endl;
		for(unsigned  int  i = 0; i < fn_micrographs.size(); ++i)
			std::cout << "    * " << fn_micrographs[i] << std::endl;
	}
#ifdef TIMING
	timer.toc(TIMING_A1);
#endif
#ifdef TIMING
	timer.tic(TIMING_A2);
#endif
	// Make sure that psi-sampling is even around the circle
	RFLOAT old_sampling = psi_sampling;
	int n_sampling = ROUND(360. / psi_sampling);
	psi_sampling = 360. / (RFLOAT) n_sampling;
	if (verb > 0 && fabs(old_sampling - psi_sampling) > 1e-3)
		std::cout << " + Changed psi-sampling rate to: " << psi_sampling << std::endl;

	// Read in the references
	Mrefs.clear();
	if (fn_ref == "gauss")
	{
		if (verb > 0)
			std::cout << " + Will use Gaussian blob as reference, with peak value of " << gauss_max_value << std::endl;

		if(particle_diameter<=0)
			CRITICAL(ERR_GAUSSBLOBSIZE);

		// Set particle boxsize to be 1.5x bigger than circle with particle_diameter
		particle_size =  1.5 * ROUND(particle_diameter/angpix);
		particle_size += particle_size%2;
		psi_sampling = 360.;
		do_ctf = false;

		Image<RFLOAT> Iref;
		Iref().initZeros(particle_size, particle_size);
		Iref().setXmippOrigin();
		// Make a Gaussian reference. sigma is 1/6th of the particle size, such that 3 sigma is at the image edge
		RFLOAT normgauss = gaussian1D(0., particle_size/6., 0.);
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Iref())
		{
			double r = sqrt((RFLOAT)(i*i + j*j));
			A2D_ELEM(Iref(), i, j) = gauss_max_value * gaussian1D(r, particle_size/6., 0.) / normgauss;
		}
		Mrefs.push_back(Iref());

	}
	else if (fn_ref.isStarFile())
	{
		MetaDataTable MDref;
		MDref.read(fn_ref);
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDref)
		{
			// Get all reference images and their names
			Image<RFLOAT> Iref;

			FileName fn_img;
			if (!MDref.getValue(EMDL_MLMODEL_REF_IMAGE, fn_img))
			{
				if (!MDref.getValue(EMDL_IMAGE_NAME, fn_img))
					REPORT_ERROR("AutoPicker::initialise ERROR: either provide rlnReferenceImage or rlnImageName in the reference STAR file!");
			}
#ifdef DEBUG
			std::cerr << " Reference fn= " << fn_img << std::endl;
#endif
			Iref.read(fn_img);
			Iref().setXmippOrigin();
			Mrefs.push_back(Iref());
		}
	}
	else
	{
		Image<RFLOAT> Istk, Iref;
		Istk.read(fn_ref);
		for (int n = 0; n < NSIZE(Istk()); n++)
		{
			Istk().getImage(n, Iref());
			Iref().setXmippOrigin();
			Mrefs.push_back(Iref());
		}
	}
#ifdef TIMING
	timer.toc(TIMING_A2);
#endif
#ifdef TIMING
	timer.tic(TIMING_A3);
#endif


	// Re-scale references if necessary
	if (angpix_ref < 0)
		angpix_ref = angpix;

	// Automated determination of bg_radius (same code as in particle_sorter.cpp!)
	if (particle_diameter < 0.)
	{
		RFLOAT sumr=0.;
		for (int iref = 0; iref < Mrefs.size(); iref++)
		{
			RFLOAT cornerval = DIRECT_MULTIDIM_ELEM(Mrefs[iref], 0);
			// Look on the central X-axis, which first and last values are NOT equal to the corner value
			bool has_set_first=false;
			bool has_set_last=false;
			int last_corner=FINISHINGX(Mrefs[iref]), first_corner=STARTINGX(Mrefs[iref]);
			for (long int j=STARTINGX(Mrefs[iref]); j<=FINISHINGX(Mrefs[iref]); j++)
			{
				if (!has_set_first)
				{
					if (fabs(A3D_ELEM(Mrefs[iref], 0,0,j) - cornerval) > 1e-6)
					{
						first_corner = j;
						has_set_first = true;
					}
				}
				else if (!has_set_last)
				{
					if (fabs(A3D_ELEM(Mrefs[iref], 0,0,j) - cornerval) < 1e-6)
					{
						last_corner = j - 1;
						has_set_last = true;
					}
				}
			}
			sumr += (last_corner - first_corner);
		}
		particle_diameter = sumr / Mrefs.size();
		// diameter is in Angstroms
		particle_diameter *= angpix_ref;
		if (verb>0)
		{
			std::cout << " + Automatically set the background diameter to " << particle_diameter << " Angstrom" << std::endl;
			std::cout << " + You can override this by providing --particle_diameter (in Angstroms)" << std::endl;
		}
	}

	if (fabs(angpix_ref - angpix) > 1e-3)
	{
		int halfoldsize = XSIZE(Mrefs[0]) / 2;
		int newsize = ROUND(halfoldsize * (angpix_ref/angpix));
		newsize *= 2;
		RFLOAT rescale_greyvalue = 1.;
		// If the references came from downscaled particles, then those were normalised differently
		// (the stddev is N times smaller after downscaling N times)
		// This needs to be corrected again
		RFLOAT rescale_factor = 1.;
		if (newsize > XSIZE(Mrefs[0]))
			rescale_factor *= (RFLOAT)(XSIZE(Mrefs[0]))/(RFLOAT)newsize;
		for (int iref = 0; iref < Mrefs.size(); iref++)
		{
			resizeMap(Mrefs[iref], newsize);
			Mrefs[iref] *= rescale_factor;
			Mrefs[iref].setXmippOrigin();
		}
	}

	// Get particle boxsize from the input reference images
	particle_size = XSIZE(Mrefs[0]);

	if (particle_diameter > particle_size * angpix_ref)
	{
		std::cerr << " particle_diameter (A): " << particle_diameter << " box_size (pix): " << particle_size << " pixel size (A): " << angpix_ref << std::endl;
		REPORT_ERROR("ERROR: the particle diameter is larger than the size of the box.");
	}


	if ( (verb > 0) && (autopick_helical_segments))
	{
		std::cout << " + Helical tube diameter = " << helical_tube_diameter << " Angstroms " << std::endl;
		std::cout << " + Helical tube diameter should be smaller than the particle (background) diameter" << std::endl;
	}
	if ( (autopick_helical_segments) && (helical_tube_diameter > particle_diameter) )
		REPORT_ERROR("Error: Helical tube diameter should be smaller than the particle (background) diameter!");


	// Get the squared particle radius (in integer pixels)
	particle_radius2 = ROUND(particle_diameter/(2. * angpix));
	particle_radius2 -= decrease_radius;
	particle_radius2*= particle_radius2;
#ifdef DEBUG
	std::cerr << " particle_size= " << particle_size << " sqrt(particle_radius2)= " << sqrt(particle_radius2) << std::endl;
#endif
	// Invert references if necessary (do this AFTER automasking them!)
	if (do_invert)
	{
		for (int iref = 0; iref < Mrefs.size(); iref++)
		{
			Mrefs[iref] *= -1.;
		}
	}

	// Get micrograph_size
	Image<RFLOAT> Imic;
	Imic.read(fn_micrographs[0], false);
	micrograph_xsize = XSIZE(Imic());
	micrograph_ysize = YSIZE(Imic());
	micrograph_size = (micrograph_xsize != micrograph_ysize) ? XMIPP_MAX(micrograph_xsize, micrograph_ysize) : micrograph_xsize;

	if (lowpass < 0.)
	{
		downsize_mic = micrograph_size;
	}
	else
	{
		downsize_mic = 2 * ROUND(micrograph_size * angpix / lowpass);
	}

	/*
	 * Here we set the size of the micrographs during cross-correlation calculation. The final size is still the same size as
	 * the input micrographs, we simply adjust the frequencies used in fourier space by cropping the frequency-space images in
	 * intermediate calculations.
	 */


	if(workFrac>1) // set size directly
	{
		int tempFrac = (int)ROUND(workFrac);
		tempFrac -= tempFrac%2;
		if(tempFrac<micrograph_size)
		{
			workSize = getGoodFourierDims(tempFrac,micrograph_size);
		}
		else
			REPORT_ERROR("workFrac larger than micrograph_size (--shrink) cannot be used. Choose a fraction 0<frac<1  OR  size<micrograph_size");
	}
	else if(workFrac<=1) // set size as fraction of original
	{

		if(workFrac>0)
		{
			int tempFrac = (int)ROUND(workFrac*(RFLOAT)micrograph_size);
			tempFrac -= tempFrac%2;
			workSize = getGoodFourierDims(tempFrac,micrograph_size);
		}
		else if(workFrac==0)
		{
			workSize = getGoodFourierDims((int)downsize_mic,micrograph_size);
		}
		else
			REPORT_ERROR("negative workFrac (--shrink) cannot be used. Choose a fraction 0<frac<1  OR size<micrograph_size");
	}
	workSize -= workSize%2; //make even in case it is not already

	if ( verb > 0 && workSize < downsize_mic)
	{
		std::cout << " + WARNING: The calculations will be done at a lower resolution than requested." << std::endl;
	}

	if ( verb > 0 && (autopick_helical_segments) && ((float(workSize) / float(micrograph_size)) < 0.4999) )
	{
		std::cerr << " + WARNING: Please consider using a shrink value 0.5~1 for picking helical segments. Smaller values may lead to poor results." << std::endl;
	}

	//printf("workSize = %d, corresponding to a resolution of %g for these settings. \n", workSize, 2*(((RFLOAT)micrograph_size*angpix)/(RFLOAT)workSize));

	if (min_particle_distance < 0)
	{
		min_particle_distance = particle_size * angpix / 2.;
	}
#ifdef TIMING
	timer.toc(TIMING_A3);
#endif
#ifdef TIMING
	timer.tic(TIMING_A4);
#endif

	// Pre-calculate and store Projectors for all references at the right size
	if (!do_read_fom_maps)
	{
		if (verb > 0)
		{
			std::cout << " Initialising FFTs for the references and masks ... " << std::endl;
		}

		// Calculate a circular mask based on the particle_diameter and then store its FT
		FourierTransformer transformer;
		MultidimArray<RFLOAT> Mcirc_mask(particle_size, particle_size);
		MultidimArray<RFLOAT> Maux(micrograph_size, micrograph_size);
		Mcirc_mask.setXmippOrigin();
		Maux.setXmippOrigin();

		// For squared difference, need the mask of the background to locally normalise the micrograph
		nr_pixels_circular_invmask = 0;
		Mcirc_mask.initZeros();
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Mcirc_mask)
		{
			if (i*i + j*j >= particle_radius2)
			{
				A2D_ELEM(Mcirc_mask, i, j) = 1.;
				nr_pixels_circular_invmask++;
			}
		}

		// Now set the mask in the large square and store its FFT
		Maux.initZeros();
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Mcirc_mask)
		{
			A2D_ELEM(Maux, i, j ) = A2D_ELEM(Mcirc_mask, i, j);
		}
		CenterFFT(Maux, true);
		transformer.FourierTransform(Maux, Finvmsk);

		// Also get the particle-area mask
		nr_pixels_circular_mask = 0;
		Mcirc_mask.initZeros();
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Mcirc_mask)
		{
			if (i*i + j*j < particle_radius2)
			{
				A2D_ELEM(Mcirc_mask, i, j) = 1.;
				nr_pixels_circular_mask++;
			}
		}

#ifdef DEBUG
		std::cerr << " min_particle_distance= " << min_particle_distance << " micrograph_size= " << micrograph_size << " downsize_mic= " << downsize_mic << std::endl;
		std::cerr << " nr_pixels_circular_mask= " << nr_pixels_circular_mask << " nr_pixels_circular_invmask= " << nr_pixels_circular_invmask << std::endl;
#endif


		PPref.clear();
		if (verb > 0)
			init_progress_bar(Mrefs.size());

		Projector PP(micrograph_size);
		MultidimArray<RFLOAT> dummy;

		for (int iref = 0; iref < Mrefs.size(); iref++)
		{

			// (Re-)apply the mask to the references
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mrefs[iref])
			{
				DIRECT_MULTIDIM_ELEM(Mrefs[iref], n) *= DIRECT_MULTIDIM_ELEM(Mcirc_mask, n);
			}

			// Set reference in the large box of the micrograph
			Maux.initZeros();
			Maux.setXmippOrigin();
			FOR_ALL_ELEMENTS_IN_ARRAY2D(Mrefs[iref])
			{
				A2D_ELEM(Maux, i, j) = A2D_ELEM(Mrefs[iref], i, j);
			}

			// And compute its Fourier Transform inside the Projector
			PP.computeFourierTransformMap(Maux, dummy, downsize_mic, 1, false);
			PPref.push_back(PP);

			if (verb > 0)
				progress_bar(iref+1);

		}

		if (verb > 0)
			progress_bar(Mrefs.size());

	}
#ifdef TIMING
	timer.toc(TIMING_A4);
	timer.toc(TIMING_A0);
#endif
#ifdef DEBUG
	std::cerr << "Finishing initialise" << std::endl;
#endif

}

#ifdef CUDA
int AutoPicker::deviceInitialise()
{
	int devCount;
	cudaGetDeviceCount(&devCount);

	std::vector < std::vector < std::string > > allThreadIDs;
	untangleDeviceIDs(gpu_ids, allThreadIDs);

	// Sequential initialisation of GPUs on all ranks
	int dev_id;
	if (!std::isdigit(*gpu_ids.begin()))
		dev_id = 0;
	else
		dev_id = textToInteger((allThreadIDs[0][0]).c_str());

	if (verb>0)
	{
		std::cout << " + Using GPU device " << dev_id << std::endl;
	}

	return(dev_id);
}
#endif

void AutoPicker::run()
{
	int barstep;
	if (verb > 0)
	{
		std::cout << " Autopicking ..." << std::endl;
		init_progress_bar(fn_micrographs.size());
		barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
	}


	FileName fn_olddir="";
	for (long int imic = 0; imic < fn_micrographs.size(); imic++)
	{
		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		// Check new-style outputdirectory exists and make it if not!
		FileName fn_dir = getOutputRootName(fn_micrographs[imic]);
		fn_dir = fn_dir.beforeLastOf("/");
		if (fn_dir != fn_olddir)
		{
			// Make a Particles directory
			int res = system(("mkdir -p " + fn_dir).c_str());
			fn_olddir = fn_dir;
		}
#ifdef TIMING
		timer.tic(TIMING_A5);
#endif
		autoPickOneMicrograph(fn_micrographs[imic], imic);
#ifdef TIMING
		timer.toc(TIMING_A5);
#endif
	}

	if (verb > 0)
		progress_bar(fn_micrographs.size());

}

void AutoPicker::pickCCFPeaks(
		const MultidimArray<RFLOAT>& Mccf,
		const MultidimArray<int>& Mclass,
		RFLOAT threshold_value,
		int peak_r_min,
		RFLOAT particle_diameter_pix,
		std::vector<ccfPeak>& ccf_peak_list,
		MultidimArray<RFLOAT>& Mccfplot,
		int skip_side, float scale)
{
	MultidimArray<int> Mrec;
	std::vector<ccfPixel> ccf_pixel_list;
	ccfPeak ccf_peak_small, ccf_peak_big;
	std::vector<ccfPeak> ccf_peak_list_aux;
	int new_micrograph_xsize = (int)((float)micrograph_xsize*scale);
	int new_micrograph_ysize = (int)((float)micrograph_ysize*scale);
	int nr_pixels;
	RFLOAT ratio;

	// Rescale skip_side and particle_diameter_pix
	skip_side = (int)((float)skip_side*scale);
	particle_diameter_pix *= scale;
	//int micrograph_core_size = XMIPP_MIN(micrograph_xsize, micrograph_ysize) - skip_side * 2 - 2;

	if ( (NSIZE(Mccf) != 1) || (ZSIZE(Mccf) != 1) || (YSIZE(Mccf) != XSIZE(Mccf)) )
		REPORT_ERROR("autopicker.cpp::pickCCFPeaks: The micrograph should be a 2D square!");
	if ( (XSIZE(Mccf) < new_micrograph_xsize) || (YSIZE(Mccf) < new_micrograph_ysize) )
		REPORT_ERROR("autopicker.cpp::pickCCFPeaks: Invalid dimensions for Mccf!");
	//if (micrograph_core_size < 100*scale)
	//	REPORT_ERROR("autopicker.cpp::pickCCFPeaks: Size of the micrograph is too small compared to that of the particle box!");
	if ( (STARTINGY(Mccf) != FIRST_XMIPP_INDEX(YSIZE(Mccf))) || (STARTINGX(Mccf) != FIRST_XMIPP_INDEX(XSIZE(Mccf))) )
		REPORT_ERROR("autopicker.cpp::pickCCFPeaks: The origin of input 3D MultidimArray is not at the center (use v.setXmippOrigin() before calling this function)!");
	if (Mccf.sameShape(Mclass) == false)
		REPORT_ERROR("autopicker.cpp::pickCCFPeaks: Mccf and Mclass should have the same shape!");
	if (peak_r_min < 1)
		REPORT_ERROR("autopicker.cpp::pickCCFPeaks: Radii of peak should be positive!");
	if (particle_diameter_pix < 5. * scale)
		REPORT_ERROR("autopicker.cpp::pickCCFPeaks: Particle diameter should be larger than 5 pixels!");

	// Init output
	ccf_peak_list.clear();
	Mccfplot.clear();
	Mccfplot.resize(Mccf);
	Mccfplot.initZeros();
	Mccfplot.setXmippOrigin();

	RFLOAT avg0, stddev0, maxccf0, minccf0;
	Mccf.computeStats(avg0, stddev0, minccf0, maxccf0);

	// Collect all high ccf pixels
	Mrec.clear();
	Mrec.resize(Mccf);
	Mrec.initConstant(0);
	Mrec.setXmippOrigin();
	nr_pixels = 0;
	for (int ii = FIRST_XMIPP_INDEX(new_micrograph_ysize) + skip_side; ii <= LAST_XMIPP_INDEX(new_micrograph_ysize) - skip_side; ii++)
	{
		for (int jj = FIRST_XMIPP_INDEX(new_micrograph_xsize) + skip_side; jj <= LAST_XMIPP_INDEX(new_micrograph_xsize) - skip_side; jj++)
		{
			RFLOAT fom = A2D_ELEM(Mccf, ii, jj);
			nr_pixels++;
			if (fom > threshold_value)
			{
				A2D_ELEM(Mrec, ii, jj) = 1;
				ccf_pixel_list.push_back(ccfPixel(jj, ii, fom));
			}
		}
	}
	std::sort(ccf_pixel_list.begin(), ccf_pixel_list.end());
#ifdef DEBUG_HELIX
	std::cout << " nr_high_ccf_pixels= " << ccf_pixel_list.size() << std::endl;
#endif

	// Do not do anything if nr_high_ccf_pixels is too small or too large! Thres is restricted to 0.01%-10% beforehand.
	if ( (nr_pixels < 100) || (ccf_pixel_list.size() < 10) )
	{
		ccf_peak_list.clear();
		return;
	}
	ratio = ((RFLOAT)(ccf_pixel_list.size())) / ((RFLOAT)(nr_pixels));
	if (ratio > 0.1)
	{
		ccf_peak_list.clear();
		return;
	}

	// Find all peaks! (From the highest fom values)
	ccf_peak_list.clear();
	for (int id = ccf_pixel_list.size() - 1; id >= 0; id--)
	{
		int x_new, y_new, x_old, y_old, rmax, rmax2, iref;
		int rmax_min = peak_r_min;
		int rmax_max;
		int iter_max = 3;
		RFLOAT fom_max;
		RFLOAT area_percentage_min = 0.8;

		// Deal with very small shrink values. But it still not performs well if workFrac < 0.2
		if ( (scale < 0.5) && (scale > 0.2) )
			area_percentage_min = 0.2 + (2. * (scale - 0.2));
		else if (scale < 0.2)
			area_percentage_min = 0.2;

		// Check if this ccf pixel is covered by another peak
		x_old = x_new = ROUND(ccf_pixel_list[id].x);
		y_old = y_new = ROUND(ccf_pixel_list[id].y);
		if (A2D_ELEM(Mrec, y_new, x_new) == 0)
			continue;

		iref = A2D_ELEM(Mclass, y_new, x_new);
		fom_max = A2D_ELEM(Mccf, y_new, x_new);

		// Pick a peak starting from this ccf pixel
		ccf_peak_small.clear();
		ccf_peak_big.clear();
		rmax_max = ROUND(particle_diameter_pix / 2.); // Sep29,2015 ????????????
		if (rmax_max < 100)
			rmax_max = 100;
		for (rmax = rmax_min; rmax < rmax_max; rmax++)
		{
			// Record the smaller peak
			ccf_peak_small = ccf_peak_big;

			//std::cout << " id= " << id << ", rmax= " << rmax << ", p= " << ccf_peak_small.area_percentage << std::endl;

			// 5 iterations to guarantee convergence??????????????
			// Require 5 iterations for stablising the center of this peak under this rmax
			for (int iter = 0; iter < iter_max; iter++)
			{
				// Empty this peak
				ccf_peak_big.clear();

				// New rmax
				rmax2 = rmax * rmax;

				// Get all ccf pixels within this rmax
				for (int dx = -rmax; dx <= rmax; dx++)
				{
					for (int dy = -rmax; dy <= rmax; dy++)
					{
						// Boundary checks
						if ( (dx * dx + dy * dy) > rmax2)
							continue;

						x_new = x_old + dx;
						y_new = y_old + dy;
						if ( (x_new < (FIRST_XMIPP_INDEX(new_micrograph_xsize) + skip_side + 1))
								|| (x_new > (LAST_XMIPP_INDEX(new_micrograph_xsize) - skip_side - 1))
								|| (y_new < (FIRST_XMIPP_INDEX(new_micrograph_ysize) + skip_side + 1))
								|| (y_new > (LAST_XMIPP_INDEX(new_micrograph_ysize) - skip_side - 1)) )
							continue;

						// Push back all ccf pixels within this rmax
						RFLOAT ccf = A2D_ELEM(Mccf, y_new, x_new);
						if (A2D_ELEM(Mrec, y_new, x_new) == 0)
							ccf = minccf0;
						ccf_peak_big.ccf_pixel_list.push_back(ccfPixel(x_new, y_new, ccf));
					}
				}
				// Check ccf_peak.ccf_pixel_list.size() below!

				// Refresh
				ccf_peak_big.r = rmax;
				ccf_peak_big.fom_thres = threshold_value;
				if (ccf_peak_big.refresh() == false)
				{
					//std::cout << " x_old, y_old = " << x_old << ", " << y_old << std::endl;
					//REPORT_ERROR("autopicker.cpp::pickCCFPeaks(): BUG No ccf pixels found within the small circle!");
					break;
				}
				x_new = ROUND(ccf_peak_big.x);
				y_new = ROUND(ccf_peak_big.y);

				// Out of range
				if ( (x_new < (FIRST_XMIPP_INDEX(new_micrograph_xsize) + skip_side + 1))
						|| (x_new > (LAST_XMIPP_INDEX(new_micrograph_xsize) - skip_side - 1))
						|| (y_new < (FIRST_XMIPP_INDEX(new_micrograph_ysize) + skip_side + 1))
						|| (y_new > (LAST_XMIPP_INDEX(new_micrograph_ysize) - skip_side - 1)) )
					break;

				// Converge
				if ( (x_old == x_new) && (y_old == y_new) )
					break;

				x_old = x_new;
				y_old = y_new;

			} // iter++ ends

			// Peak finding is over if peak area does not expand
			if (ccf_peak_big.area_percentage < area_percentage_min)
				break;

		} // rmax++ ends

		// A peak is found
		if (ccf_peak_small.isValid())
		{
			for (int ii = 0; ii < ccf_peak_small.ccf_pixel_list.size(); ii++)
			{
				x_new = ROUND(ccf_peak_small.ccf_pixel_list[ii].x);
				y_new = ROUND(ccf_peak_small.ccf_pixel_list[ii].y);
				A2D_ELEM(Mrec, y_new, x_new) = 0;
			}
			// TODO: if r > ...? do not include this peak?
			ccf_peak_small.ref = iref;
			ccf_peak_small.fom_max = fom_max;
			ccf_peak_list.push_back(ccf_peak_small);
			//std::cout << ccf_peak_list.size() << ", "<< std::flush;
		}
	}// id-- ends
	// Sort the peaks from the weakest to the strongest
	std::sort(ccf_peak_list.begin(), ccf_peak_list.end());
#ifdef DEBUG_HELIX
	std::cout << " nr_peaks= " << ccf_peak_list.size() << std::endl;
#endif

	// Remove too close peaks (retain the stronger ones while remove the weaker)
	Mrec.clear();
	Mrec.resize(Mccf);
	Mrec.initConstant(-1);
	Mrec.setXmippOrigin();
	// Sort the peaks from the weakest to the strongest
	for (int new_id = 0; new_id < ccf_peak_list.size(); new_id++)
	{
		int x, y, peak_r, old_id;
		RFLOAT dist2;
		RFLOAT peak_r2 = ccf_peak_list[new_id].r * ccf_peak_list[new_id].r;
		if (ccf_peak_list[new_id].r > 0.)
			peak_r = CEIL(ccf_peak_list[new_id].r);
		else
			peak_r = (-1.);

		// Remove peaks with too small/big radii!
		if ( (peak_r <= 1) || (peak_r > (particle_diameter_pix / 2.)) )
		{
			ccf_peak_list[new_id].r = (-1.);
			continue;
		}
		for (int dx = -peak_r; dx <= peak_r; dx++)
		{
			for (int dy = -peak_r; dy <= peak_r; dy++)
			{
				dist2 = (RFLOAT)(dx * dx + dy * dy);
				if (dist2 > peak_r2)
					continue;

				x = dx + ROUND(ccf_peak_list[new_id].x);
				y = dy + ROUND(ccf_peak_list[new_id].y);

				// Out of range
				if ( (x < (FIRST_XMIPP_INDEX(new_micrograph_xsize) + skip_side + 1))
						|| (x > (LAST_XMIPP_INDEX(new_micrograph_xsize) - skip_side - 1))
						|| (y < (FIRST_XMIPP_INDEX(new_micrograph_ysize) + skip_side + 1))
						|| (y > (LAST_XMIPP_INDEX(new_micrograph_ysize) - skip_side - 1)) )
					continue;

				old_id = A2D_ELEM(Mrec, y, x);
				if (old_id >= 0)
					ccf_peak_list[old_id].r = (-1.);
				A2D_ELEM(Mrec, y, x) = new_id;
			}
		}
	}

	// Collect all valid peaks
	ccf_peak_list_aux.clear();
	for (int id = 0; id < ccf_peak_list.size(); id++)
	{
		if (ccf_peak_list[id].isValid())
			ccf_peak_list_aux.push_back(ccf_peak_list[id]);
	}
	ccf_peak_list.clear();
	ccf_peak_list = ccf_peak_list_aux;
	ccf_peak_list_aux.clear();
#ifdef DEBUG_HELIX
	std::cout << " nr_peaks_pruned= " << ccf_peak_list.size() << std::endl;
#endif

	// TODO: Remove all discrete peaks (one peak should have at least two neighbouring peaks within r=particle_radius)

	// Plot
	for (int ii = 0; ii < ccf_peak_list.size(); ii++)
	{
		for (int jj = 0; jj < ccf_peak_list[ii].ccf_pixel_list.size(); jj++)
		{
			int x, y;

			if (ccf_peak_list[ii].ccf_pixel_list[jj].fom < ccf_peak_list[ii].fom_thres)
				continue;

			x = ROUND(ccf_peak_list[ii].ccf_pixel_list[jj].x);
			y = ROUND(ccf_peak_list[ii].ccf_pixel_list[jj].y);
			A2D_ELEM(Mccfplot, y, x) = 1.;
		}
	}

	return;
}

void AutoPicker::extractHelicalTubes(
		std::vector<ccfPeak>& peak_list,
		std::vector<std::vector<ccfPeak> >& tube_coord_list,
		std::vector<RFLOAT>& tube_len_list,
		std::vector<std::vector<ccfPeak> >& tube_track_list,
		RFLOAT particle_diameter_pix,
		RFLOAT curvature_factor_max,
		RFLOAT interbox_distance_pix,
		RFLOAT tube_diameter_pix,
		float scale)
{
	std::vector<int> is_peak_on_other_tubes;
	std::vector<int> is_peak_on_this_tube;
	int tube_id;
	RFLOAT curvature_max;

	tube_coord_list.clear();
	tube_len_list.clear();
	tube_track_list.clear();

	//Rescaling
	particle_diameter_pix *= scale;
	interbox_distance_pix *= scale;
	tube_diameter_pix *= scale;

	if (particle_diameter_pix < 5. * scale)
		REPORT_ERROR("autopicker.cpp::extractHelicalTubes: Particle diameter should be larger than 5 pixels!");
	if ( (curvature_factor_max < 0.0001) || (curvature_factor_max > 1.0001) )
		REPORT_ERROR("autopicker.cpp::extractHelicalTubes: Factor of curvature should be 0~1!");
	if ( (interbox_distance_pix < 0.9999) || (interbox_distance_pix > particle_diameter_pix) )
		REPORT_ERROR("autopicker.cpp::extractHelicalTubes: Interbox distance should be > 1 pixel and < particle diameter!");
	if ( (tube_diameter_pix < 1.) || (tube_diameter_pix > particle_diameter_pix) )
		REPORT_ERROR("autopicker.cpp::extractHelicalTubes: Tube diameter should be > 1 pixel and < particle diameter!");
	if (peak_list.size() < 5)
		return;

	// Calculate the maximum curvature
	curvature_max = curvature_factor_max / (particle_diameter_pix / 2.);
	//curvature_max = (sqrt(1. / scale)) * curvature_factor_max / (particle_diameter_pix / 2.);

	// Sort the peaks from the weakest to the strongest
	std::sort(peak_list.begin(), peak_list.end());

	is_peak_on_other_tubes.resize(peak_list.size());
	is_peak_on_this_tube.resize(peak_list.size());
	for (int peak_id0 = 0; peak_id0 < is_peak_on_other_tubes.size(); peak_id0++)
		is_peak_on_other_tubes[peak_id0] = is_peak_on_this_tube[peak_id0] = -1;

	// Traverse peaks from the strongest to the weakest
	tube_id = 0;
	for (int peak_id0 = peak_list.size() - 1; peak_id0 >= 0; peak_id0--)
	{
		RFLOAT rmax2;
		std::vector<ccfPeak> selected_peaks;

		// Check whether this peak is included on other tubes
		if (is_peak_on_other_tubes[peak_id0] > 0)
			continue;

		// Probably a new tube
		tube_id++;
		is_peak_on_other_tubes[peak_id0] = tube_id;
		for (int peak_id1 = 0; peak_id1 < peak_list.size(); peak_id1++)
			is_peak_on_this_tube[peak_id1] = -1;
		is_peak_on_this_tube[peak_id0] = tube_id;

		// Gather all neighboring peaks around
		selected_peaks.clear(); // don't push itself in? No do not push itself!!!
		rmax2 = particle_diameter_pix * particle_diameter_pix / 4.;
		for (int peak_id1 = 0; peak_id1 < peak_list.size(); peak_id1++)
		{
			if (peak_id0 == peak_id1)
				continue;
			if (is_peak_on_other_tubes[peak_id1] > 0)
				continue;

			RFLOAT dx, dy, dist2;
			dx = peak_list[peak_id1].x - peak_list[peak_id0].x;
			dy = peak_list[peak_id1].y - peak_list[peak_id0].y;
			dist2 = dx * dx + dy * dy;
			if (dist2 < rmax2)
			{
				ccfPeak myPeak = peak_list[peak_id1];
				myPeak.dist = sqrt(dist2);
				if ( (fabs(dx) < 0.01) && (fabs(dy) < 0.01) )
					myPeak.psi = 0.;
				else
					myPeak.psi = RAD2DEG(atan2(dy, dx));
				selected_peaks.push_back(myPeak);
			}
		}

		// Sep29,2015 ????????????
		// If less than 3 neighboring peaks are found, this is not a peak along a helical tube!
		if (selected_peaks.size() <= 2)
			continue;

		// This peak has >=2 neighboring peaks! Try to find an orientation!
		RFLOAT local_psi, local_dev, best_local_psi, best_local_dev, dev0, dev1, dev_weights;
		RFLOAT local_psi_sampling = 0.1;
		std::vector<ccfPeak> selected_peaks_dir1, selected_peaks_dir2, helical_track_dir1, helical_track_dir2, helical_track, helical_segments;
		RFLOAT psi_dir1, psi_dir2, len_dir1, len_dir2;

		selected_peaks_dir1.clear();
		selected_peaks_dir2.clear();

		// Find the averaged psi
		best_local_psi = -1.;
		best_local_dev = (1e30);
		// Traverse every possible value of local_psi and calculate the dev
		for (local_psi = 0.; local_psi < 180.; local_psi += local_psi_sampling)
		{
			local_dev = 0.;
			dev_weights = 0.;
			for (int peak_id1 = 0; peak_id1 < selected_peaks.size(); peak_id1++)
			{
				dev0 = ABS(selected_peaks[peak_id1].psi - local_psi);
				if (dev0 > 180.)
					dev0 = ABS(dev0 - 360.);
				if (dev0 > 90.)
					dev0 = ABS(dev0 - 180.);

				RFLOAT pixel_count = (RFLOAT)(selected_peaks[peak_id1].nr_peak_pixel);
				if (pixel_count < 1.)
					pixel_count = 1.;
				local_dev += dev0 * pixel_count;
				dev_weights += pixel_count;
			}
			local_dev /= dev_weights;

			// Refresh if a better local psi is found
			if (local_dev < best_local_dev)
			{
				best_local_psi = local_psi;
				best_local_dev = local_dev;
			}
		}
		// Sort all peaks into dir1, dir2 and others
		psi_dir1 = psi_dir2 = 0.;
		for (int peak_id1 = 0; peak_id1 < selected_peaks.size(); peak_id1++)
		{
			dev0 = ABS(selected_peaks[peak_id1].psi - best_local_psi);
			dev1 = dev0;
			if (dev1 > 180.)
				dev1 = ABS(dev1 - 360.);
			if (dev1 > 90.)
				dev1 = ABS(dev1 - 180.);
			RFLOAT curvature1 = DEG2RAD(dev1) / selected_peaks[peak_id1].dist;

			// Cannot fall into the estimated direction
			if (curvature1 > curvature_max)
				continue;

			// Psi direction or the opposite direction
			if (fabs(dev1 - dev0) < 0.1)
			{
				selected_peaks_dir2.push_back(selected_peaks[peak_id1]);
				psi_dir2 += selected_peaks[peak_id1].psi;
			}
			else
			{
				selected_peaks_dir1.push_back(selected_peaks[peak_id1]);
				psi_dir1 += selected_peaks[peak_id1].psi;
			}
		}

		RFLOAT xc, yc, xc_new, yc_new, xc_old, yc_old, dist_max, nr_psi_within_range;

		//std::cout << " nr Dir1 peaks= " << selected_peaks_dir1.size() << std::endl;
		// Dir1
		if (selected_peaks_dir1.size() >= 1)
		{
			// Init
			psi_dir1 /= selected_peaks_dir1.size();
			dist_max = -1.;
			for (int peak_id1 = 0; peak_id1 < selected_peaks_dir1.size(); peak_id1++)
			{
				if (selected_peaks_dir1[peak_id1].dist > dist_max)
					dist_max = selected_peaks_dir1[peak_id1].dist;
			}
			len_dir1 = 0.;
			xc_old = peak_list[peak_id0].x;
			yc_old = peak_list[peak_id0].y;
			helical_track_dir1.clear();

			while(1)
			{
				// A new center along helical track dir1 is found, record it
				xc_new = xc_old + cos(DEG2RAD(psi_dir1)) * dist_max;
				yc_new = yc_old + sin(DEG2RAD(psi_dir1)) * dist_max;
				len_dir1 += dist_max;

				ccfPeak myPeak;
				myPeak.x = xc_new;
				myPeak.y = yc_new;
				myPeak.psi = psi_dir1;
				helical_track_dir1.push_back(myPeak);
				//std::cout << " Dir1 new center: x, y, psi= " << xc << ", " << yc << ", " << psi_dir1 << std::endl;
				// TODO: other parameters to add?

				// TODO: mark peaks along helical tracks
				xc = (xc_old + xc_new) / 2.;
				yc = (yc_old + yc_new) / 2.;
				rmax2 = ((dist_max + tube_diameter_pix) / 2.) * ((dist_max + tube_diameter_pix) / 2.);
				bool is_new_peak_found = false;
				bool is_combined_with_another_tube = true;
				for (int peak_id1 = 0; peak_id1 < peak_list.size(); peak_id1++)
				{
					RFLOAT dx, dy, dist, dist2, dpsi, h, r;
					dx = peak_list[peak_id1].x - xc;
					dy = peak_list[peak_id1].y - yc;
					dist2 = dx * dx + dy * dy;

					if (dist2 > rmax2)
						continue;

					if ( (fabs(dx) < 0.01) && (fabs(dy) < 0.01) ) // atan2(0,0)
						dpsi = 0.;
					else
						dpsi = RAD2DEG(atan2(dy, dx)) - psi_dir1;
					dist = sqrt(dist2);
					h = dist * fabs(cos(DEG2RAD(dpsi)));
					r = dist * fabs(sin(DEG2RAD(dpsi)));

					if ( (h < ((dist_max + tube_diameter_pix) / 2.)) && (r < (tube_diameter_pix / 2.)) )
					{
						if (is_peak_on_this_tube[peak_id1] < 0)
						{
							is_new_peak_found = true;
							is_peak_on_this_tube[peak_id1] = tube_id;
							if (is_peak_on_other_tubes[peak_id1] < 0)
							{
								is_combined_with_another_tube = false;
								is_peak_on_other_tubes[peak_id1] = tube_id;
							}
						}
					}
				}
				if ( (is_new_peak_found == false) || (is_combined_with_another_tube == true) )
				{
					// TODO: delete the end of this track list or not?
					//helical_track_dir1.pop_back();
					break;
				}

				// TODO: try to find another new center if possible
				xc_old = xc_new;
				yc_old = yc_new;
				rmax2 = particle_diameter_pix * particle_diameter_pix / 4.;
				selected_peaks_dir1.clear();
				for (int peak_id1 = 0; peak_id1 < peak_list.size(); peak_id1++)
				{
					if (is_peak_on_this_tube[peak_id1] > 0)
						continue;

					RFLOAT dx, dy, dist, dist2, dpsi, h, r;
					dx = peak_list[peak_id1].x - xc_old;
					dy = peak_list[peak_id1].y - yc_old;
					dist2 = dx * dx + dy * dy;
					if (dist2 < rmax2)
					{
						myPeak = peak_list[peak_id1];
						myPeak.dist = sqrt(dist2);
						if ( (fabs(dx) < 0.01) && (fabs(dy) < 0.01) ) // atan2(0,0)
							myPeak.psi = 0.;
						else
							myPeak.psi = RAD2DEG(atan2(dy, dx));
						selected_peaks_dir1.push_back(myPeak);
					}
				}

				dist_max = -1.;
				RFLOAT psi_sum = 0.;
				RFLOAT psi_weights = 0.;
				nr_psi_within_range = 0.;
				int id_peak_dist_max;
				for (int peak_id1 = 0; peak_id1 < selected_peaks_dir1.size(); peak_id1++)
				{
					//std::cout << "  Peak id " << selected_peaks_dir1[ii].id << " x, y, r, psi, psidir1= " << selected_peaks_dir1[ii].x << ", " << selected_peaks_dir1[ii].y
					//		<< ", " << selected_peaks_dir1[ii].r << ", " << selected_peaks_dir1[ii].psi << ", " << psi_dir1 << std::endl;

					RFLOAT curvature = DEG2RAD(ABS(selected_peaks_dir1[peak_id1].psi - psi_dir1)) / selected_peaks_dir1[peak_id1].dist;
					if (curvature < curvature_max)
					{
						nr_psi_within_range += 1.;

						RFLOAT pixel_count = (RFLOAT)(selected_peaks_dir1[peak_id1].nr_peak_pixel);
						if (pixel_count < 1.)
							pixel_count = 1.;
						psi_sum += selected_peaks_dir1[peak_id1].psi * pixel_count;
						psi_weights += pixel_count;

						if (selected_peaks_dir1[peak_id1].dist > dist_max)
						{
							dist_max = selected_peaks_dir1[peak_id1].dist;
							id_peak_dist_max = peak_id1;
						}
					}
				}

				// If no peaks are found in this round, the helical track stops, exit
				if (nr_psi_within_range < 0.5)
					break;
				psi_dir1 = psi_sum / psi_weights;
			}
		}

		//std::cout << " nr Dir2 peaks= " << selected_peaks_dir2.size() << std::endl;
		// Dir2
		// ================================================================================================
		if (selected_peaks_dir2.size() >= 1)
		{
			// Init
			psi_dir2 /= selected_peaks_dir2.size();
			dist_max = -1.;
			for (int peak_id1 = 0; peak_id1 < selected_peaks_dir2.size(); peak_id1++)
			{
				if (selected_peaks_dir2[peak_id1].dist > dist_max)
					dist_max = selected_peaks_dir2[peak_id1].dist;
			}
			len_dir2 = 0.;
			xc_old = peak_list[peak_id0].x;
			yc_old = peak_list[peak_id0].y;
			helical_track_dir2.clear();

			while(1)
			{
				// A new center along helical track dir1 is found, record it
				xc_new = xc_old + cos(DEG2RAD(psi_dir2)) * dist_max;
				yc_new = yc_old + sin(DEG2RAD(psi_dir2)) * dist_max;
				len_dir2 += dist_max;

				ccfPeak myPeak;
				myPeak.x = xc_new;
				myPeak.y = yc_new;
				myPeak.psi = psi_dir2;
				helical_track_dir2.push_back(myPeak);
				//std::cout << " Dir1 new center: x, y, psi= " << xc << ", " << yc << ", " << psi_dir1 << std::endl;
				// TODO: other parameters to add?

				// TODO: mark peaks along helical tracks
				xc = (xc_old + xc_new) / 2.;
				yc = (yc_old + yc_new) / 2.;
				rmax2 = ((dist_max + tube_diameter_pix) / 2.) * ((dist_max + tube_diameter_pix) / 2.);
				bool is_new_peak_found = false;
				bool is_combined_with_another_tube = true;
				for (int peak_id1 = 0; peak_id1 < peak_list.size(); peak_id1++)
				{
					RFLOAT dx, dy, dist, dist2, dpsi, h, r;
					dx = peak_list[peak_id1].x - xc;
					dy = peak_list[peak_id1].y - yc;
					dist2 = dx * dx + dy * dy;

					if (dist2 > rmax2)
						continue;

					if ( (fabs(dx) < 0.01) && (fabs(dy) < 0.01) ) // atan2(0,0)
						dpsi = 0.;
					else
						dpsi = RAD2DEG(atan2(dy, dx)) - psi_dir2;
					dist = sqrt(dist2);
					h = dist * fabs(cos(DEG2RAD(dpsi)));
					r = dist * fabs(sin(DEG2RAD(dpsi)));

					if ( (h < ((dist_max + tube_diameter_pix) / 2.)) && (r < (tube_diameter_pix / 2.)) )
					{
						if (is_peak_on_this_tube[peak_id1] < 0)
						{
							is_new_peak_found = true;
							is_peak_on_this_tube[peak_id1] = tube_id;
							if (is_peak_on_other_tubes[peak_id1] < 0)
							{
								is_combined_with_another_tube = false;
								is_peak_on_other_tubes[peak_id1] = tube_id;
							}
						}
					}
				}
				if ( (is_new_peak_found == false) || (is_combined_with_another_tube == true) )
				{
					// TODO: delete the end of this track list or not?
					break;
				}

				// TODO: try to find another new center if possible
				xc_old = xc_new;
				yc_old = yc_new;
				rmax2 = particle_diameter_pix * particle_diameter_pix / 4.;
				selected_peaks_dir2.clear();
				for (int peak_id1 = 0; peak_id1 < peak_list.size(); peak_id1++)
				{
					if (is_peak_on_this_tube[peak_id1] > 0)
						continue;

					RFLOAT dx, dy, dist, dist2, dpsi, h, r;
					dx = peak_list[peak_id1].x - xc_old;
					dy = peak_list[peak_id1].y - yc_old;
					dist2 = dx * dx + dy * dy;
					if (dist2 < rmax2)
					{
						myPeak = peak_list[peak_id1];
						myPeak.dist = sqrt(dist2);
						if ( (fabs(dx) < 0.01) && (fabs(dy) < 0.01) ) // atan2(0,0)
							myPeak.psi = 0.;
						else
							myPeak.psi = RAD2DEG(atan2(dy, dx));
						selected_peaks_dir2.push_back(myPeak);
					}
				}

				dist_max = -1.;
				RFLOAT psi_sum = 0.;
				RFLOAT psi_weights = 0.;
				nr_psi_within_range = 0.;
				int id_peak_dist_max;
				for (int peak_id1 = 0; peak_id1 < selected_peaks_dir2.size(); peak_id1++)
				{
					//std::cout << "  Peak id " << selected_peaks_dir2[ii].id << " x, y, r, psi, psidir2= " << selected_peaks_dir2[ii].x << ", " << selected_peaks_dir2[ii].y
					//		<< ", " << selected_peaks_dir2[ii].r << ", " << selected_peaks_dir2[ii].psi << ", " << psi_dir2 << std::endl;

					RFLOAT curvature = DEG2RAD(ABS(selected_peaks_dir2[peak_id1].psi - psi_dir2)) / selected_peaks_dir2[peak_id1].dist;
					if (curvature < curvature_max)
					{
						nr_psi_within_range += 1.;

						RFLOAT pixel_count = (RFLOAT)(selected_peaks_dir2[peak_id1].nr_peak_pixel);
						if (pixel_count < 1.)
							pixel_count = 1.;
						psi_sum += selected_peaks_dir2[peak_id1].psi * pixel_count;
						psi_weights += pixel_count;

						if (selected_peaks_dir2[peak_id1].dist > dist_max)
						{
							dist_max = selected_peaks_dir2[peak_id1].dist;
							id_peak_dist_max = peak_id1;
						}
					}
				}

				// If no peaks are found in this round, the helical track stops, exit
				if (nr_psi_within_range < 0.5)
					break;
				psi_dir2 = psi_sum / psi_weights;
			}
		}


		// Get a full track
		helical_track.clear();
		for (int ii = helical_track_dir2.size() - 1; ii >= 0; ii--)
			helical_track.push_back(helical_track_dir2[ii]);
		helical_track.push_back(peak_list[peak_id0]);
		for (int ii = 0; ii < helical_track_dir1.size(); ii++)
			helical_track.push_back(helical_track_dir1[ii]);

		// TODO: check below !!!
		if ( (helical_track.size() < 3)
				|| ((len_dir1 + len_dir2) < particle_diameter_pix)
				|| ((len_dir1 + len_dir2) < interbox_distance_pix) )
		{
			helical_track.clear();
		}
		else
		{
			ccfPeak newSegment;
			RFLOAT dist_left, len_total;

			helical_segments.clear();

			// Get the first segment
			newSegment.x = helical_track[0].x;
			newSegment.y = helical_track[0].y;
			newSegment.psi = RAD2DEG(atan2(helical_track[1].y - helical_track[0].y, helical_track[1].x - helical_track[0].x));
			newSegment.ref = helical_track[0].ref;
			helical_segments.push_back(newSegment);

			// Get segments along the track
			dist_left = 0.;
			for (int inext = 1; inext < helical_track.size(); inext++)
			{
				RFLOAT x0, y0, dx, dy, dist, dist_total, psi, nr_segments_float;
				int nr_segments_int;

				x0 = helical_track[inext - 1].x;
				y0 = helical_track[inext - 1].y;
				dx = helical_track[inext].x - helical_track[inext - 1].x;
				dy = helical_track[inext].y - helical_track[inext - 1].y;
				psi = RAD2DEG(atan2(dy, dx));
				dist_total = sqrt(dx * dx + dy * dy);

				nr_segments_float = (dist_left + dist_total) / interbox_distance_pix;
				nr_segments_int = FLOOR(nr_segments_float);
				if (nr_segments_int >= 1)
				{
					for (int iseg = 1; iseg <= nr_segments_int; iseg++)
					{
						dist = (RFLOAT)(iseg) * interbox_distance_pix - dist_left;
						dx = cos(DEG2RAD(psi)) * dist;
						dy = sin(DEG2RAD(psi)) * dist;

						newSegment.x = x0 + dx;
						newSegment.y = y0 + dy;
						newSegment.psi = psi;
						if ( (iseg * 2) < nr_segments_int)
							newSegment.ref = helical_track[inext - 1].ref;
						else
							newSegment.ref = helical_track[inext].ref;
						helical_segments.push_back(newSegment);
					}
				}

				dist_left = (dist_left + dist_total) - ((RFLOAT)(nr_segments_int) * interbox_distance_pix);
			}

			// Get the last segment and mark it as invalid (different from what I did for the first segment)
			int last_id = helical_track.size();
			last_id -= 1;
			newSegment.x = helical_track[last_id].x;
			newSegment.y = helical_track[last_id].y;
			newSegment.psi = (1e30);
			newSegment.ref = helical_track[last_id].ref;
			helical_segments.push_back(newSegment);

			len_total = len_dir1 + len_dir2;
			tube_coord_list.push_back(helical_segments);
			tube_len_list.push_back(len_total);
			tube_track_list.push_back(helical_track);

			// DEBUG
#ifdef DEBUG_HELIX
			for (int ii = 0; ii < helical_track.size(); ii++)
				std::cout << "Track point x, y, psi= " << helical_track[ii].x << ", " << helical_track[ii].y << ", " << helical_track[ii].psi << std::endl;
			std::cout << " Track length= " << (len_dir1 + len_dir2) << std::endl;
#endif
		}
	}

	return;
}

void AutoPicker::exportHelicalTubes(
		const MultidimArray<RFLOAT>& Mccf,
		MultidimArray<RFLOAT>& Mccfplot,
		const MultidimArray<int>& Mclass,
		std::vector<std::vector<ccfPeak> >& tube_coord_list,
		std::vector<std::vector<ccfPeak> >& tube_track_list,
		std::vector<RFLOAT>& tube_len_list,
		FileName& fn_mic_in,
		FileName& fn_star_out,
		RFLOAT particle_diameter_pix,
		RFLOAT tube_length_min_pix,
		int skip_side, float scale)
{
	// Rescale particle_diameter_pix, tube_length_min_pix, skip_side
	tube_length_min_pix *= scale;
	particle_diameter_pix *= scale;
	skip_side = (int)((float)skip_side*scale);

	if ( (tube_coord_list.size() != tube_track_list.size()) || (tube_track_list.size() != tube_len_list.size()) )
		REPORT_ERROR("autopicker.cpp::exportHelicalTubes: BUG tube_coord_list.size() != tube_track_list.size() != tube_len_list.size()");
	if ( (STARTINGY(Mccf) != FIRST_XMIPP_INDEX(YSIZE(Mccf))) || (STARTINGX(Mccf) != FIRST_XMIPP_INDEX(XSIZE(Mccf))) )
		REPORT_ERROR("autopicker.cpp::exportHelicalTubes: The origin of input 3D MultidimArray is not at the center (use v.setXmippOrigin() before calling this function)!");
	if (particle_diameter_pix < 5.) // TODO: 5?
		REPORT_ERROR("autopicker.cpp::exportHelicalTubes: Particle diameter should be larger than 5 pixels!");

	// Mark tracks on Mccfplot
	Mccfplot.setXmippOrigin();
	for (int itrack = 0; itrack < tube_track_list.size(); itrack++)
	{
		for (int icoord = 1; icoord < tube_track_list[itrack].size(); icoord++)
		{
			RFLOAT x0, y0, x1, y1, dx, dy, psi_rad, dist_total;
			int x_int, y_int;

			x0 = tube_track_list[itrack][icoord - 1].x;
			y0 = tube_track_list[itrack][icoord - 1].y;
			x1 = tube_track_list[itrack][icoord].x;
			y1 = tube_track_list[itrack][icoord].y;
			dx = x1 - x0;
			dy = y1 - y0;
			if ( (fabs(dx) < 0.1) && (fabs(dy) < 0.1) )
				psi_rad = 0.;
			psi_rad = atan2(dy, dx);

			dist_total = sqrt(dx * dx + dy * dy);
			if (dist_total < 2.)
				continue;

			for (RFLOAT fdist = 1.; fdist < dist_total; fdist += 1.)
			{
				dx = cos(psi_rad) * fdist;
				dy = sin(psi_rad) * fdist;
				x1 = x0 + dx;
				y1 = y0 + dy;
				x_int = ROUND(x1);
				y_int = ROUND(y1);

				if ( (x_int < (FIRST_XMIPP_INDEX(micrograph_xsize) + 1))
						|| (x_int > (LAST_XMIPP_INDEX(micrograph_xsize) - 1))
						|| (y_int < (FIRST_XMIPP_INDEX(micrograph_ysize) + 1))
						|| (y_int > (LAST_XMIPP_INDEX(micrograph_ysize) - 1)) )
					continue;

				A2D_ELEM(Mccfplot, y_int, x_int) = 1.;
			}
		}
	}

	// Detect crossovers
	RFLOAT dist2_min = particle_diameter_pix * particle_diameter_pix / 4.;
	for (int itube1 = 0; (itube1 + 1) < tube_coord_list.size(); itube1++)
	{
		for (int icoord1 = 0; icoord1 < tube_coord_list[itube1].size(); icoord1++)
		{
			// Coord1 selected
			for (int itube2 = (itube1 + 1); itube2 < tube_coord_list.size(); itube2++)
			{
				for (int icoord2 = 0; icoord2 < tube_coord_list[itube2].size(); icoord2++)
				{
					// Coord2 selected
					RFLOAT x1, y1, x2, y2, dist2;
					x1 = tube_coord_list[itube1][icoord1].x;
					y1 = tube_coord_list[itube1][icoord1].y;
					x2 = tube_coord_list[itube2][icoord2].x;
					y2 = tube_coord_list[itube2][icoord2].y;
					dist2 = (x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1);

					// If this point is around the crossover
					if (dist2 < dist2_min)
						tube_coord_list[itube1][icoord1].psi = tube_coord_list[itube2][icoord2].psi = (1e30);
				}
			}
		}
	}

	// Cancel segments close to the ends of tubes
	/*
	for (int itube = 0; itube < tube_coord_list.size(); itube++)
	{
		if (tube_track_list[itube].size() < 2)
			continue;

		RFLOAT x_start, y_start, x_end, y_end, particle_radius_pix2;
		int last_id;

		last_id = tube_track_list[itube].size();
		last_id -= 1;

		x_start = tube_track_list[itube][0].x;
		y_start = tube_track_list[itube][0].y;
		x_end = tube_track_list[itube][last_id].x;
		y_end = tube_track_list[itube][last_id].y;
		particle_radius_pix2 = particle_diameter_pix * particle_diameter_pix / 4.;

		for (int icoord = 0; icoord < tube_coord_list[itube].size(); icoord++)
		{
			if (fabs(tube_coord_list[itube][icoord].psi) > 360.)
				continue;

			RFLOAT x, y, dx1, dy1, dx2, dy2, dist21, dist22;

			x = tube_coord_list[itube][icoord].x;
			y = tube_coord_list[itube][icoord].y;
			dx1 = x - x_start;
			dy1 = y - y_start;
			dx2 = x - x_end;
			dy2 = y - y_end;
			dist21 = dx1 * dx1 + dy1 * dy1;
			dist22 = dx2 * dx2 + dy2 * dy2;

			if ( (dist21 < particle_radius_pix2) || (dist22 < particle_radius_pix2) )
				tube_coord_list[itube][icoord].psi = (1e30);
		}
	}
	*/

	// Write out a STAR file with the coordinates
	FileName fn_tmp;
	MetaDataTable MDout;
	int helical_tube_id;
	RFLOAT helical_tube_len;

	// Only output STAR header if there are no tubes...
	MDout.clear();
	MDout.addLabel(EMDL_IMAGE_COORD_X);
	MDout.addLabel(EMDL_IMAGE_COORD_Y);
	MDout.addLabel(EMDL_PARTICLE_CLASS);
	MDout.addLabel(EMDL_PARTICLE_AUTOPICK_FOM);
	MDout.addLabel(EMDL_PARTICLE_HELICAL_TUBE_ID);
	MDout.addLabel(EMDL_ORIENT_TILT_PRIOR);
	MDout.addLabel(EMDL_ORIENT_PSI_PRIOR);
	MDout.addLabel(EMDL_PARTICLE_HELICAL_TRACK_LENGTH);
	MDout.addLabel(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO);

	helical_tube_id = 0;
	for (int itube = 0; itube < tube_coord_list.size(); itube++)
	{
		if (tube_length_min_pix > particle_diameter_pix)
		{
			if (tube_len_list[itube] < tube_length_min_pix)
				continue;
		}
		helical_tube_id++;
		helical_tube_len = 0.;
		for (int icoord = 0; icoord < tube_coord_list[itube].size(); icoord++)
		{
			int x_int, y_int, iref;
			RFLOAT fom;

			if (icoord > 0)
			{
				RFLOAT dx = ((RFLOAT)(tube_coord_list[itube][icoord].x)) - ((RFLOAT)(tube_coord_list[itube][icoord - 1].x));
				RFLOAT dy = ((RFLOAT)(tube_coord_list[itube][icoord].y)) - ((RFLOAT)(tube_coord_list[itube][icoord - 1].y));
				helical_tube_len += sqrt(dx * dx + dy * dy);
			}

			// Invalid psi (crossover)
			if (fabs(tube_coord_list[itube][icoord].psi) > 360.)
				continue;

			x_int = ROUND(tube_coord_list[itube][icoord].x);
			y_int = ROUND(tube_coord_list[itube][icoord].y);

			// Out of range
			if ( (x_int < (FIRST_XMIPP_INDEX(micrograph_xsize) + skip_side + 1))
					|| (x_int > (LAST_XMIPP_INDEX(micrograph_xsize) - skip_side - 1))
					|| (y_int < (FIRST_XMIPP_INDEX(micrograph_ysize) + skip_side + 1))
					|| (y_int > (LAST_XMIPP_INDEX(micrograph_ysize) - skip_side - 1)) )
				continue;

			iref = A2D_ELEM(Mclass, y_int, x_int);
			fom = A2D_ELEM(Mccf, y_int, x_int);

			MDout.addObject();
			RFLOAT xval = (tube_coord_list[itube][icoord].x / scale) - (RFLOAT)(FIRST_XMIPP_INDEX(micrograph_xsize));
			RFLOAT yval = (tube_coord_list[itube][icoord].y / scale) - (RFLOAT)(FIRST_XMIPP_INDEX(micrograph_ysize));
			MDout.setValue(EMDL_IMAGE_COORD_X, xval);
			MDout.setValue(EMDL_IMAGE_COORD_Y, yval);
			MDout.setValue(EMDL_PARTICLE_CLASS, iref + 1); // start counting at 1
			MDout.setValue(EMDL_PARTICLE_AUTOPICK_FOM, fom);
			MDout.setValue(EMDL_PARTICLE_HELICAL_TUBE_ID, helical_tube_id);
			MDout.setValue(EMDL_ORIENT_TILT_PRIOR, 90.);
			MDout.setValue(EMDL_ORIENT_PSI_PRIOR, (-1.) * (tube_coord_list[itube][icoord].psi)); // Beware! Multiplied by -1!
			MDout.setValue(EMDL_PARTICLE_HELICAL_TRACK_LENGTH, helical_tube_len);
			MDout.setValue(EMDL_ORIENT_PSI_PRIOR_FLIP_RATIO, BIMODAL_PSI_PRIOR_FLIP_RATIO);
		}
	}

	fn_tmp = getOutputRootName(fn_mic_in) + "_" + fn_star_out + ".star";
	MDout.write(fn_tmp);

	return;
}

void AutoPicker::autoPickOneMicrograph(FileName &fn_mic, long int imic)
{

	Image<RFLOAT> Imic;
	MultidimArray<Complex > Faux, Faux2, Fmic;
	MultidimArray<RFLOAT> Maux, Mstddev, Mmean, Mdiff2, MsumX2, Mccf_best, Mpsi_best, Fctf, Mccf_best_combined;
	MultidimArray<int> Mclass_best_combined;
	FourierTransformer transformer;
	RFLOAT sum_ref_under_circ_mask, sum_ref2_under_circ_mask;
	int my_skip_side = autopick_skip_side + particle_size/2;
	CTF ctf;

	int min_distance_pix = ROUND(min_particle_distance / angpix);
	float scale = (float)workSize / (float)micrograph_size;

	// Always use the same random seed
	init_random_generator(random_seed + imic);

#ifdef DEBUG
	Image<RFLOAT> tt;
	tt().resize(micrograph_size, micrograph_size);
	std::cerr << " fn_mic= " << fn_mic << std::endl;
#endif
#ifdef TIMING
	timer.tic(TIMING_A6);
#endif
	// Read in the micrograph
	Imic.read(fn_mic);
	Imic().setXmippOrigin();
#ifdef TIMING
	timer.toc(TIMING_A6);
#endif
	// Let's just check the square size again....
	RFLOAT my_size, my_xsize, my_ysize;
	my_xsize = XSIZE(Imic());
	my_ysize = YSIZE(Imic());
	my_size = (my_xsize != my_ysize) ? XMIPP_MAX(my_xsize, my_ysize) : my_xsize;

	if (my_size != micrograph_size || my_xsize != micrograph_xsize || my_ysize != micrograph_ysize)
	{
		Imic().printShape();
		std::cerr << " micrograph_size= " << micrograph_size << " micrograph_xsize= " << micrograph_xsize << " micrograph_ysize= " << micrograph_ysize << std::endl;
		REPORT_ERROR("AutoPicker::autoPickOneMicrograph ERROR: No differently sized micrographs are allowed in one run, sorry you will have to run separately for each size...");
	}
#ifdef TIMING
	timer.tic(TIMING_A7);
#endif
	// Set mean to zero and stddev to 1 to prevent numerical problems with one-sweep stddev calculations....
    RFLOAT avg0, stddev0, minval0, maxval0;
	Imic().computeStats(avg0, stddev0, minval0, maxval0);
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imic())
	{
		// Remove pixel values that are too far away from the mean
		if ( ABS(DIRECT_MULTIDIM_ELEM(Imic(), n) - avg0) / stddev0 > outlier_removal_zscore)
			DIRECT_MULTIDIM_ELEM(Imic(), n) = avg0;

		DIRECT_MULTIDIM_ELEM(Imic(), n) = (DIRECT_MULTIDIM_ELEM(Imic(), n) - avg0) / stddev0;
	}

	if (micrograph_xsize != micrograph_ysize)
	{
		// Window non-square micrographs to be a square with the largest side
		rewindow(Imic, micrograph_size);

		// Fill region outside the original window with white Gaussian noise to prevent all-zeros in Mstddev
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Imic())
		{
			if (i < FIRST_XMIPP_INDEX(micrograph_ysize)
					|| i > LAST_XMIPP_INDEX(micrograph_ysize)
					|| j < FIRST_XMIPP_INDEX(micrograph_xsize)
					|| j > LAST_XMIPP_INDEX(micrograph_xsize) )
				A2D_ELEM(Imic(), i, j) = rnd_gaus(0.,1.);
		}
	}
#ifdef TIMING
	timer.toc(TIMING_A7);
#endif
#ifdef TIMING
	timer.tic(TIMING_A8);
#endif
	// Read in the CTF information if needed
	if (do_ctf)
	{
		// Search for this micrograph in the metadata table
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmic)
		{
			FileName fn_tmp;
			MDmic.getValue(EMDL_MICROGRAPH_NAME, fn_tmp);
			if (fn_tmp==fn_mic)
			{
				ctf.read(MDmic, MDmic);
				Fctf.resize(downsize_mic, downsize_mic/2 + 1);
				ctf.getFftwImage(Fctf, micrograph_size, micrograph_size, angpix, false, false, intact_ctf_first_peak, true);
				break;
			}
		}
#ifdef DEBUG
		std::cerr << " Read CTF info from" << fn_mic.withoutExtension()<<"_ctf.star" << std::endl;
		Image<RFLOAT> Ictf;
		Ictf()=Fctf;
		Ictf.write("Mmic_ctf.spi");
#endif
	}
#ifdef TIMING
	timer.toc(TIMING_A8);
#endif
#ifdef TIMING
	timer.tic(TIMING_A9);
#endif
	Mccf_best.resize(workSize, workSize);
	Mpsi_best.resize(workSize, workSize);
#ifdef TIMING
	timer.toc(TIMING_A9);
#endif
#ifdef TIMING
	timer.tic(TIMING_B1);
#endif
	//Sjors 18apr2016
	RFLOAT normfft = (RFLOAT)(micrograph_size * micrograph_size) / (RFLOAT)nr_pixels_circular_mask;
	if (do_read_fom_maps)
	{
		FileName fn_tmp=getOutputRootName(fn_mic)+"_"+fn_out+"_stddevNoise.spi";
		Image<RFLOAT> It;
		It.read(fn_tmp);
		Mstddev = It();
	}
	else
	{
		/*
		 * Squared difference FOM:
		 * Sum ( (X-mu)/sig  - A )^2 =
		 *  = Sum((X-mu)/sig)^2 - 2 Sum (A*(X-mu)/sig) + Sum(A)^2
		 *  = (1/sig^2)*Sum(X^2) - (2*mu/sig^2)*Sum(X) + (mu^2/sig^2)*Sum(1) - (2/sig)*Sum(AX) + (2*mu/sig)*Sum(A) + Sum(A^2)
		 *
		 * However, the squared difference with an "empty" ie all-zero reference is:
		 * Sum ( (X-mu)/sig)^2
		 *
		 * The ratio of the probabilities thereby becomes:
		 * P(ref) = 1/sqrt(2pi) * exp (( (X-mu)/sig  - A )^2 / -2 )   // assuming sigma = 1!
		 * P(zero) = 1/sqrt(2pi) * exp (( (X-mu)/sig )^2 / -2 )
		 *
		 * P(ref)/P(zero) = exp(( (X-mu)/sig  - A )^2 / -2) / exp ( ( (X-mu)/sig )^2 / -2)
		 *                = exp( (- (2/sig)*Sum(AX) + (2*mu/sig)*Sum(A) + Sum(A^2)) / - 2 )
		 *
		 *                Therefore, I do not need to calculate (X-mu)/sig beforehand!!!
		 *
		 */

        // Fourier Transform (and downscale) Imic()
		CenterFFT(Imic(), true);
		transformer.FourierTransform(Imic(), Fmic);

		if (highpass > 0.)
        {
			lowPassFilterMap(Fmic, micrograph_size, highpass, angpix, 2, true); // true means highpass instead of lowpass!
        	transformer.inverseFourierTransform(Fmic, Imic()); // also calculate inverse transform again for squared calculation below
        }

		// Also calculate the FFT of the squared micrograph
		Maux.resize(micrograph_size,micrograph_size);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Maux)
		{
			DIRECT_MULTIDIM_ELEM(Maux, n) = DIRECT_MULTIDIM_ELEM(Imic(), n) * DIRECT_MULTIDIM_ELEM(Imic(), n);
		}
		MultidimArray<Complex> Fmic2;
		transformer.FourierTransform(Maux, Fmic2);

		Maux.resize(workSize,workSize);

#ifdef DEBUG
		std::cerr << " nr_pixels_circular_invmask= " << nr_pixels_circular_invmask << std::endl;
		std::cerr << " nr_pixels_circular_mask= " << nr_pixels_circular_mask << std::endl;
		windowFourierTransform(Finvmsk, Faux2, micrograph_size);
		transformer.inverseFourierTransform(Faux2, tt());
		CenterFFT(tt(), false);
		tt.write("Minvmask.spi");
#endif

		// The following calculate mu and sig under the solvent area at every position in the micrograph
		calculateStddevAndMeanUnderMask(Fmic, Fmic2, Finvmsk, nr_pixels_circular_invmask, Mstddev, Mmean);

		if (do_write_fom_maps)
		{
			// TMP output
			FileName fn_tmp=getOutputRootName(fn_mic)+"_"+fn_out+"_stddevNoise.spi";
			Image<RFLOAT> It;
			It() = Mstddev;
			It.write(fn_tmp);
		}

		// From now on use downsized Fmic, as the cross-correlation with the references can be done at lower resolution
		windowFourierTransform(Fmic, Faux, downsize_mic);
		Fmic = Faux;

	}// end if do_read_fom_maps
#ifdef TIMING
	timer.toc(TIMING_B1);
#endif

	// Now start looking for the peaks of all references
	// Clear the output vector with all peaks
	std::vector<Peak> peaks;
	peaks.clear();

	if (autopick_helical_segments)
	{
		if (do_read_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It_float;
			Image<int> It_int;

			fn_tmp = getOutputRootName(fn_mic)+"_"+fn_out+"_combinedCCF.spi";
			It_float.read(fn_tmp);
			Mccf_best_combined = It_float();

			fn_tmp = getOutputRootName(fn_mic)+"_"+fn_out+"_combinedCLASS.spi";
			It_int.read(fn_tmp);
			Mclass_best_combined = It_int();
		}
		else
		{
			Mccf_best_combined.clear();
			Mccf_best_combined.resize(workSize, workSize);
			Mccf_best_combined.initConstant(-99.e99);
			Mclass_best_combined.clear();
			Mclass_best_combined.resize(workSize, workSize);
			Mclass_best_combined.initConstant(-1);
		}
	}

	for (int iref = 0; iref < Mrefs.size(); iref++)
	{
		RFLOAT expected_Pratio; // the expectedFOM for this (ctf-corrected) reference
		if (do_read_fom_maps)
		{
#ifdef TIMING
			timer.tic(TIMING_B2);
#endif
			if (!autopick_helical_segments)
			{
				FileName fn_tmp;
				Image<RFLOAT> It;

				fn_tmp.compose(getOutputRootName(fn_mic)+"_"+fn_out+"_ref", iref,"_bestCCF.spi");
				It.read(fn_tmp);
				Mccf_best = It();
				It.MDMainHeader.getValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);  // Retrieve expected_Pratio from the header of the image

				fn_tmp.compose(getOutputRootName(fn_mic)+"_"+fn_out+"_ref", iref,"_bestPSI.spi");
				It.read(fn_tmp);
				Mpsi_best = It();
			}
#ifdef TIMING
			timer.toc(TIMING_B2);
#endif
		} //end else if do_read_fom_maps
		else
		{
#ifdef TIMING
			timer.tic(TIMING_B3);
#endif
			Mccf_best.initConstant(-LARGE_NUMBER);
			bool is_first_psi = true;
			for (RFLOAT psi = 0. ; psi < 360.; psi+=psi_sampling)
			{
				// Get the Euler matrix
				Matrix2D<RFLOAT> A(3,3);
				Euler_angles2matrix(0., 0., psi, A);

				// Now get the FT of the rotated (non-ctf-corrected) template
				Faux.initZeros(downsize_mic, downsize_mic/2 + 1);
				PPref[iref].get2DFourierTransform(Faux, A, IS_NOT_INV);

#ifdef DEBUG
				std::cerr << " psi= " << psi << std::endl;
				windowFourierTransform(Faux, Faux2, micrograph_size);
				tt().resize(micrograph_size, micrograph_size);
				transformer.inverseFourierTransform(Faux2, tt());
				CenterFFT(tt(), false);
				tt.write("Mref_rot.spi");

				windowFourierTransform(Fmic, Faux2, micrograph_size);
				transformer.inverseFourierTransform(Faux2, tt());
				CenterFFT(tt(), false);
				tt.write("Mmic.spi");

#endif
#ifdef TIMING
	timer.tic(TIMING_B4);
#endif
				// Apply the CTF on-the-fly (so same PPref can be used for many different micrographs)
				if (do_ctf)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
					{
						DIRECT_MULTIDIM_ELEM(Faux, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
#ifdef TIMING
	timer.toc(TIMING_B4);
#endif
#ifdef DEBUG
				MultidimArray<RFLOAT> ttt(micrograph_size, micrograph_size);
				windowFourierTransform(Faux, Faux2, micrograph_size);
				transformer.inverseFourierTransform(Faux2, ttt);
				CenterFFT(ttt, false);
				ttt.setXmippOrigin();
				tt().resize(particle_size, particle_size);
				tt().setXmippOrigin();
				FOR_ALL_ELEMENTS_IN_ARRAY2D(tt())
				{
					A2D_ELEM(tt(), i, j) = A2D_ELEM(ttt, i, j);
				}
				tt.write("Mref_rot_ctf.spi");
#endif
				}

				if (is_first_psi)
				{
#ifdef TIMING
	timer.tic(TIMING_B5);
#endif
					// Calculate the expected ratio of probabilities for this CTF-corrected reference
					// and the sum_ref_under_circ_mask and sum_ref_under_circ_mask2
					// Do this also if we're not recalculating the fom maps...
					// This calculation needs to be done on an "non-shrinked" micrograph, in order to get the correct I^2 statistics
					windowFourierTransform(Faux, Faux2, micrograph_size);
					Maux.resize(micrograph_size, micrograph_size);
					transformer.inverseFourierTransform(Faux2, Maux);
					CenterFFT(Maux, false);
					Maux.setXmippOrigin();
//#ifdef DEBUG
					Image<RFLOAT> ttt;
					ttt()=Maux;
					ttt.write("Maux.spi");
//#endif
					sum_ref_under_circ_mask = 0.;
					sum_ref2_under_circ_mask = 0.;
					RFLOAT suma2 = 0.;
					RFLOAT sumn = 1.;
					MultidimArray<RFLOAT> Mctfref(particle_size, particle_size);
					Mctfref.setXmippOrigin();
					FOR_ALL_ELEMENTS_IN_ARRAY2D(Mctfref) // only loop over smaller Mctfref, but take values from large Maux!
					{
						if (i*i + j*j < particle_radius2)
						{
							suma2 += A2D_ELEM(Maux, i, j) * A2D_ELEM(Maux, i, j);
							suma2 += 2. * A2D_ELEM(Maux, i, j) * rnd_gaus(0., 1.);
							sum_ref_under_circ_mask += A2D_ELEM(Maux, i, j);
							sum_ref2_under_circ_mask += A2D_ELEM(Maux, i, j) * A2D_ELEM(Maux, i, j);
							sumn += 1.;
						}
#ifdef DEBUG
						A2D_ELEM(Mctfref, i, j) = A2D_ELEM(Maux, i, j);
#endif
					}
					sum_ref_under_circ_mask /= sumn;
					sum_ref2_under_circ_mask /= sumn;
					expected_Pratio = exp(suma2 / (2. * sumn));
#ifdef DEBUG
					std::cerr << " expected_Pratio["<<iref<<"]= " << expected_Pratio << std::endl;
					tt()=Mctfref;
					tt.write("Mctfref.spi");
					std::cerr << "suma2 " << suma2<< " sumn " << sumn << " suma2/2sumn="<< suma2 / (2. * sumn) << std::endl;
					std::cerr << " nr_pixels_under_mask= " << nr_pixels_circular_mask << " nr_pixels_under_invmask= " << nr_pixels_circular_invmask << std::endl;
					std::cerr << "sum_ref_under_circ_mask " << sum_ref_under_circ_mask << std::endl;
					std::cerr << "sum_ref2_under_circ_mask " << sum_ref2_under_circ_mask << std::endl;
					std::cerr << "expected_Pratio " << expected_Pratio << std::endl;
#endif

					// Maux goes back to the workSize
					Maux.resize(workSize, workSize);
#ifdef TIMING
	timer.toc(TIMING_B5);
#endif
				}

#ifdef TIMING
	timer.tic(TIMING_B6);
#endif
				// Now multiply template and micrograph to calculate the cross-correlation
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
				{
					DIRECT_MULTIDIM_ELEM(Faux, n) = conj(DIRECT_MULTIDIM_ELEM(Faux, n)) * DIRECT_MULTIDIM_ELEM(Fmic, n);
				}

				// If we're not doing shrink, then Faux is bigger than Faux2!
				windowFourierTransform(Faux, Faux2, workSize);
				transformer.inverseFourierTransform(Faux2, Maux);
				CenterFFT(Maux, false);
#ifdef DEBUG
				tt()=Maux*normfft;
				tt.write("Mcc.spi");
#endif

				// Calculate ratio of prabilities P(ref)/P(zero)
				// Keep track of the best values and their corresponding iref and psi

				// So now we already had precalculated: Mdiff2 = 1/sig*Sum(X^2) - 2/sig*Sum(X) + mu^2/sig*Sum(1)
				// Still to do (per reference): - 2/sig*Sum(AX) + 2*mu/sig*Sum(A) + Sum(A^2)
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Maux)
				{
					RFLOAT diff2 = - 2. * normfft * DIRECT_MULTIDIM_ELEM(Maux, n);
					diff2 += 2. * DIRECT_MULTIDIM_ELEM(Mmean, n) * sum_ref_under_circ_mask;
					if (DIRECT_MULTIDIM_ELEM(Mstddev, n) > 1E-10)
						diff2 /= DIRECT_MULTIDIM_ELEM(Mstddev, n);
					diff2 += sum_ref2_under_circ_mask;
					diff2 = exp(- diff2 / 2.); // exponentiate to reflect the Gaussian error model. sigma=1 after normalization, 0.4=1/sqrt(2pi)

					// Store fraction of (1 - probability-ratio) wrt  (1 - expected Pratio)
					diff2 = (diff2 - 1.) / (expected_Pratio - 1.);
#ifdef DEBUG
					DIRECT_MULTIDIM_ELEM(Maux, n) = diff2;
#endif
					if (diff2 > DIRECT_MULTIDIM_ELEM(Mccf_best, n))
					{
						DIRECT_MULTIDIM_ELEM(Mccf_best, n) = diff2;
						DIRECT_MULTIDIM_ELEM(Mpsi_best, n) = psi;
					}
				}
#ifdef DEBUG
				std::cerr << " Maux.computeMax()= " << Maux.computeMax() << std::endl;
				tt()=Maux;
				tt.write("Mccf.spi");
			    std::cerr << " Press any key to continue... "  << std::endl;
			    char c;
			    std::cin >> c;
#endif
			    is_first_psi = false;
#ifdef TIMING
	timer.toc(TIMING_B6);
#endif
			} // end for psi
#ifdef TIMING
	timer.toc(TIMING_B3);
#endif
#ifdef TIMING
	timer.tic(TIMING_B7);
#endif
			if (do_write_fom_maps && !autopick_helical_segments)
			{
				FileName fn_tmp;
				Image<RFLOAT> It;

				It() = Mccf_best;
				It.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);  // Store expected_Pratio in the header of the image
				fn_tmp.compose(getOutputRootName(fn_mic)+"_"+fn_out+"_ref", iref,"_bestCCF.spi");
				It.write(fn_tmp);

				It() = Mpsi_best;
				fn_tmp.compose(getOutputRootName(fn_mic)+"_"+fn_out+"_ref", iref,"_bestPSI.spi");
				It.write(fn_tmp);

//				for (long int n=0; n<((Mccf_best).nzyxdim/10); n+=1)
//				{
//					std::cerr << DIRECT_MULTIDIM_ELEM(Mccf_best, n) << std::endl;
//				}
//				exit(0);
			} // end if do_write_fom_maps
#ifdef TIMING
	timer.toc(TIMING_B7);
#endif
		} // end if do_read_fom_maps
#ifdef TIMING
	timer.tic(TIMING_B8);
#endif
		if (autopick_helical_segments)
		{
			if (!do_read_fom_maps)
			{
				// Combine Mccf_best and Mpsi_best from all refs
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mccf_best)
				{
					RFLOAT new_ccf = DIRECT_MULTIDIM_ELEM(Mccf_best, n);
					RFLOAT old_ccf = DIRECT_MULTIDIM_ELEM(Mccf_best_combined, n);
					if (new_ccf > old_ccf)
					{
						DIRECT_MULTIDIM_ELEM(Mccf_best_combined, n) = new_ccf;
						DIRECT_MULTIDIM_ELEM(Mclass_best_combined, n) = iref;
					}
				}
			}
		}
		else
		{
			// Now that we have Mccf_best and Mpsi_best, get the peaks
			std::vector<Peak> my_ref_peaks;

			Mstddev.setXmippOrigin();
			Mccf_best.setXmippOrigin();
			Mpsi_best.setXmippOrigin();

			peakSearch(Mccf_best, Mpsi_best, Mstddev, iref, my_skip_side, my_ref_peaks, scale);
			prunePeakClusters(my_ref_peaks, min_distance_pix, scale);
			peaks.insert(peaks.end(), my_ref_peaks.begin(), my_ref_peaks.end());  // append the peaks of this reference to all the other peaks
		}
#ifdef TIMING
	timer.toc(TIMING_B8);
#endif
	} // end for iref


	if (autopick_helical_segments)
	{
		RFLOAT thres = min_fraction_expected_Pratio;
		int peak_r_min = 1;
		std::vector<ccfPeak> ccf_peak_list;
		std::vector<std::vector<ccfPeak> > tube_coord_list, tube_track_list;
		std::vector<RFLOAT> tube_len_list;
		MultidimArray<RFLOAT> Mccfplot;

		Mccf_best_combined.setXmippOrigin();
		Mclass_best_combined.setXmippOrigin();
		pickCCFPeaks(Mccf_best_combined, Mclass_best_combined, thres, peak_r_min, (particle_diameter / angpix),
				ccf_peak_list, Mccfplot, my_skip_side, scale);
		extractHelicalTubes(ccf_peak_list, tube_coord_list, tube_len_list, tube_track_list,
				(particle_diameter / angpix), helical_tube_curvature_factor_max,
				(min_particle_distance / angpix), (helical_tube_diameter / angpix), scale);
		exportHelicalTubes(Mccf_best_combined, Mccfplot, Mclass_best_combined,
					tube_coord_list, tube_track_list, tube_len_list,
					fn_mic, fn_out,
					(particle_diameter / angpix),
					(helical_tube_length_min / angpix),
					my_skip_side, scale);

		if (do_write_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It_float;
			Image<int> It_int;

			It_float() = Mccf_best_combined;
			fn_tmp = getOutputRootName(fn_mic) + "_" + fn_out + "_combinedCCF.spi";
			It_float.write(fn_tmp);

			It_int() = Mclass_best_combined;
			fn_tmp = getOutputRootName(fn_mic) + + "_" + fn_out + "_combinedCLASS.spi";
			It_int.write(fn_tmp);
		} // end if do_write_fom_maps

		if (do_write_fom_maps || do_read_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It;

			It() = Mccfplot;
			fn_tmp =  getOutputRootName(fn_mic) + "_" + fn_out + "_combinedPLOT.spi";
			It.write(fn_tmp);
		}
	}
	else
	{
#ifdef TIMING
	timer.tic(TIMING_B9);
#endif
		//Now that we have done all references, prune the list again...
		prunePeakClusters(peaks, min_distance_pix, scale);
		// And remove all too close neighbours
		removeTooCloselyNeighbouringPeaks(peaks, min_distance_pix, scale);
		// Write out a STAR file with the coordinates
		MetaDataTable MDout;
		for (int ipeak =0; ipeak < peaks.size(); ipeak++)
		{
			MDout.addObject();
			MDout.setValue(EMDL_IMAGE_COORD_X, (RFLOAT)(peaks[ipeak].x) / scale);
			MDout.setValue(EMDL_IMAGE_COORD_Y, (RFLOAT)(peaks[ipeak].y) / scale);
			MDout.setValue(EMDL_PARTICLE_CLASS, peaks[ipeak].ref + 1); // start counting at 1
			MDout.setValue(EMDL_PARTICLE_AUTOPICK_FOM, peaks[ipeak].fom);
			MDout.setValue(EMDL_ORIENT_PSI, peaks[ipeak].psi);
		}
		FileName fn_tmp = getOutputRootName(fn_mic) + "_" + fn_out + ".star";
		MDout.write(fn_tmp);
#ifdef TIMING
	timer.toc(TIMING_B9);
#endif
	}
}

FileName AutoPicker::getOutputRootName(FileName fn_mic)
{

	FileName fn_pre, fn_jobnr, fn_post;
	decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
	return fn_odir + fn_post.withoutExtension();

}

void AutoPicker::calculateStddevAndMeanUnderMask(const MultidimArray<Complex > &_Fmic, const MultidimArray<Complex > &_Fmic2,
		MultidimArray<Complex > &_Fmsk, int nr_nonzero_pixels_mask, MultidimArray<RFLOAT> &_Mstddev, MultidimArray<RFLOAT> &_Mmean)
{

	MultidimArray<Complex > Faux, Faux2;
	MultidimArray<RFLOAT> Maux(workSize, workSize);
	FourierTransformer transformer;

	_Mstddev.initZeros(workSize, workSize);
	RFLOAT normfft = (RFLOAT)(micrograph_size * micrograph_size) / (RFLOAT)nr_nonzero_pixels_mask;

	// Calculate convolution of micrograph and mask, to get average under mask at all points
	Faux.resize(_Fmic);
#ifdef DEBUG
	Image<RFLOAT> tt;
#endif

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
	{
		DIRECT_MULTIDIM_ELEM(Faux, n) = DIRECT_MULTIDIM_ELEM(_Fmic, n) * conj(DIRECT_MULTIDIM_ELEM(_Fmsk, n));
	}
	windowFourierTransform(Faux, Faux2, workSize);
	transformer.inverseFourierTransform(Faux2, Maux);
	Maux *= normfft;
	_Mmean = Maux;
	CenterFFT(_Mmean, false);

#ifdef DEBUG
	tt()=Maux;
	CenterFFT(tt(), false);
	tt.write("Mavg_mic.spi");
#endif

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(_Mstddev)
	{
		// store minus average-squared already in _Mstddev
		DIRECT_MULTIDIM_ELEM(_Mstddev, n) = -DIRECT_MULTIDIM_ELEM(Maux, n) * DIRECT_MULTIDIM_ELEM(Maux, n);
	}

	// Calculate convolution of micrograph-squared and mask
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
	{
		DIRECT_MULTIDIM_ELEM(Faux, n) = DIRECT_MULTIDIM_ELEM(_Fmic2, n) * conj(DIRECT_MULTIDIM_ELEM(_Fmsk, n));
	}
	windowFourierTransform(Faux, Faux2, workSize);
	transformer.inverseFourierTransform(Faux2, Maux);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(_Mstddev)
	{
		// we already stored minus average-squared in _Mstddev
		DIRECT_MULTIDIM_ELEM(_Mstddev, n) += normfft * DIRECT_MULTIDIM_ELEM(Maux, n);
		if (DIRECT_MULTIDIM_ELEM(_Mstddev, n) > (RFLOAT)1E-10)
			DIRECT_MULTIDIM_ELEM(_Mstddev, n) = sqrt(DIRECT_MULTIDIM_ELEM(_Mstddev, n) );
		else
			DIRECT_MULTIDIM_ELEM(_Mstddev, n) = 1.;
	}

	CenterFFT(_Mstddev, false);

#ifdef DEBUG
	tt()=_Mstddev;
	tt.write("Msig_mic.spi");
#endif
}

void AutoPicker::peakSearch(const MultidimArray<RFLOAT> &Mfom, const MultidimArray<RFLOAT> &Mpsi, const MultidimArray<RFLOAT> &Mstddev, int iref,
		int skip_side, std::vector<Peak> &peaks, float scale)
{

	peaks.clear();
	Peak peak;
	peak.ref = iref;

	skip_side = (int)((float)skip_side*scale);

	// Skip the pixels along the side of the micrograph!
	// At least 1, so dont have to check for the borders!
	skip_side = XMIPP_MAX(1, skip_side);
	for (int i = FIRST_XMIPP_INDEX((int)((float)micrograph_ysize*scale)) + skip_side; i <= LAST_XMIPP_INDEX((int)((float)micrograph_ysize*scale)) - skip_side; i++)
	{
		for (int j = FIRST_XMIPP_INDEX((int)((float)micrograph_xsize*scale)) + skip_side; j <= LAST_XMIPP_INDEX((int)((float)micrograph_xsize*scale)) - skip_side; j++)
		{

			RFLOAT myval = A2D_ELEM(Mfom, i, j);
			// check if this element is above the threshold
			if (myval  >= min_fraction_expected_Pratio)
			{

				// Only check stddev in the noise areas if max_stddev_noise is positive!
				if (max_stddev_noise > 0. && A2D_ELEM(Mstddev, i, j) > max_stddev_noise)
					continue;

				if (scale < 1.)
				{
					// When we use shrink, then often peaks aren't 5 pixels big anymore....
					if (A2D_ELEM(Mfom, i-1, j) > myval )
						continue;
					if (A2D_ELEM(Mfom, i+1, j) > myval )
						continue;
					if (A2D_ELEM(Mfom, i, j-1) > myval )
						continue;
					if (A2D_ELEM(Mfom, i, j+1) > myval )
						continue;
				}
				else
				{
					// This is a peak if all four neighbours are also above the threshold, AND have lower values than myval
					if (A2D_ELEM(Mfom, i-1, j) < min_fraction_expected_Pratio || A2D_ELEM(Mfom, i-1, j) > myval )
						continue;
					if (A2D_ELEM(Mfom, i+1, j) < min_fraction_expected_Pratio || A2D_ELEM(Mfom, i+1, j) > myval )
						continue;
					if (A2D_ELEM(Mfom, i, j-1) < min_fraction_expected_Pratio || A2D_ELEM(Mfom, i, j-1) > myval )
						continue;
					if (A2D_ELEM(Mfom, i, j+1) < min_fraction_expected_Pratio || A2D_ELEM(Mfom, i, j+1) > myval )
						continue;
				}
				peak.x = j - FIRST_XMIPP_INDEX((int)((float)micrograph_xsize*scale));
				peak.y = i - FIRST_XMIPP_INDEX((int)((float)micrograph_ysize*scale));
				peak.psi = A2D_ELEM(Mpsi, i, j);
				peak.fom = A2D_ELEM(Mfom, i, j);
				peak.relative_fom = myval;
				peaks.push_back(peak);
			}
		}
	}

}

void AutoPicker::prunePeakClusters(std::vector<Peak> &peaks, int min_distance, float scale)
{
	int mind2 = (float)(min_distance*min_distance)*scale*scale;
	int nclus = 0;

	std::vector<Peak> pruned_peaks;
	while (peaks.size() > 0)
	{
		nclus++;
		std::vector<Peak> cluster;
		cluster.push_back(peaks[0]);
		peaks.erase(peaks.begin());
		for (int iclus = 0; iclus < cluster.size(); iclus++)
		{
			int my_x = cluster[iclus].x;
			int my_y = cluster[iclus].y;
			for (int ipeakp = 0; ipeakp < peaks.size(); ipeakp++)
			{
				int dx = my_x - peaks[ipeakp].x;
				int dy = my_y - peaks[ipeakp].y;
				if (dx*dx + dy*dy < ( (float)(particle_radius2)*scale*scale ))
				{
					// Put ipeakp in the cluster, and remove from the peaks list
					cluster.push_back(peaks[ipeakp]);
					peaks.erase(peaks.begin()+ipeakp);
					ipeakp--;
				}
			}
		}

		// Now search for the peak from the cluster with the best ccf.
		// Then search again if there are any other peaks in the cluster that are further than particle_diameter apart from the selected peak
		// If so, again search for the maximum
		int ipass = 0;
		while (cluster.size() > 0)
		{
			RFLOAT best_relative_fom=-1.;
			Peak bestpeak;
			for (int iclus = 0; iclus < cluster.size(); iclus++)
			{
				if ( cluster[iclus].relative_fom > best_relative_fom)
				{
					best_relative_fom = cluster[iclus].relative_fom;
					bestpeak = cluster[iclus];
				}
			}

			// Store this peak as pruned
			pruned_peaks.push_back(bestpeak);

			// Remove all peaks within mind2 from the clusters
			for (int iclus = 0; iclus < cluster.size(); iclus++)
			{
				int dx = cluster[iclus].x - bestpeak.x;
				int dy = cluster[iclus].y - bestpeak.y;
				if (dx*dx + dy*dy < mind2)
				{
					cluster.erase(cluster.begin()+iclus);
					iclus--;
				}
			}
			ipass++;
		}
	} // end while peaks.size > 0

	// Set the pruned peaks back into the input vector
	peaks = pruned_peaks;

}

void AutoPicker::removeTooCloselyNeighbouringPeaks(std::vector<Peak> &peaks, int min_distance, float scale)
{
	// Now only keep those peaks that are at least min_particle_distance number of pixels from any other peak
	std::vector<Peak> pruned_peaks;
	int mind2 = (float)(min_distance*min_distance)*scale*scale;
	for (int ipeak = 0; ipeak < peaks.size(); ipeak++)
	{
		int my_x = peaks[ipeak].x;
		int my_y = peaks[ipeak].y;
		int my_mind2 = 99999;
		for (int ipeakp = 0; ipeakp < peaks.size(); ipeakp++)
		{
			if (ipeakp != ipeak)
			{
				int dx = peaks[ipeakp].x - my_x;
				int dy = peaks[ipeakp].y - my_y;
				int d2 = dx*dx + dy*dy;
				if ( d2 < (int)(((float)my_mind2)*scale*scale) )
					my_mind2 = d2;
			}
		}
		if (my_mind2 > mind2)
			pruned_peaks.push_back(peaks[ipeak]);
	}

	// Set the pruned peaks back into the input vector
	peaks = pruned_peaks;

}

int AutoPicker::largestPrime(int query)
{
	int i(2), primeF(query);
	while (i*i<=primeF)
	{
		if(primeF%i!=0)
			i+=1;
		else
			primeF /= i;
	}
	return primeF;
}

int AutoPicker::getGoodFourierDims(int requestedSizeRealX, int lim)
{

	if (!do_optimise_scale)
		return requestedSizeRealX;

	int inputPrimeF =XMIPP_MAX(largestPrime(requestedSizeRealX),largestPrime(requestedSizeRealX/2+1));
	if(inputPrimeF<=LARGEST_ACCEPTABLE_PRIME)
	{
		if (verb > 0)
			std::cout << " + Will use micrographs scaled to " << requestedSizeRealX << " pixels as requested. The largest prime factor in FFTs is " << inputPrimeF << std::endl;
		return requestedSizeRealX;
	}

	int S_up = LARGEST_ACCEPTABLE_PRIME;
	int S_down = LARGEST_ACCEPTABLE_PRIME;

	// Search upwards - can take a long time if unlucky and/or small LARGEST_ACCEPTABLE_PRIME
	int currentU = requestedSizeRealX;
	S_up = 			 largestPrime(currentU);
	S_up = XMIPP_MAX(largestPrime(currentU/2+1),S_up);
	while(S_up>=LARGEST_ACCEPTABLE_PRIME  && currentU<=(lim+2))
	{
		currentU += 2;
		S_up = 			 largestPrime(currentU);
		S_up = XMIPP_MAX(largestPrime(currentU/2+1),S_up);
	}


	// Search downwards - guaranteed to find in reasonable time
	int currentD = requestedSizeRealX;
	S_down = 		   largestPrime(currentD);
	S_down = XMIPP_MAX(largestPrime(currentD/2+1),S_down);
	while(S_down>=LARGEST_ACCEPTABLE_PRIME)
	{
		currentD -= 2;
		S_down = 		   largestPrime(currentD);
		S_down = XMIPP_MAX(largestPrime(currentD/2+1),S_down);
	}


	if (verb > 0)
	{
		std::cout << " + WARNING: Requested rescale of micrographs is " << requestedSizeRealX << " pixels. The largest prime factor in FFTs is " << inputPrimeF << std::endl;
	}
	if((currentU-requestedSizeRealX)>(requestedSizeRealX-currentD) || (currentU>lim))
	{
		if (verb > 0)
		{
			std::cout << " + WARNING: Will change rescaling of micrographs to " << currentD << " pixels, because the prime factor then becomes " <<  S_down << std::endl;
			std::cout << " + WARNING: add --skip_optimise_scale to your autopick command to prevent rescaling " << std::endl;
		}
		return currentD;
	}
	else
	{
		if (verb > 0)
		{
			std::cout << " + WARNING: Will change rescaling of micrographs to " << currentU << " pixels, because the prime factor then becomes " <<  S_up << std::endl;
			std::cout << " + WARNING: add --skip_optimise_scale to your autopick command to prevent rescaling " << std::endl;
		}
		return currentU;
	}

}

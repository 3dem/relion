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

#include "src/postprocessing.h"

void Postprocessing::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);
	int gen_section = parser.addSection("General options");
	fn_I1 = parser.getOption("--i", "Input name of half1, e.g. run_half1_class001_unfil.mrc", "");
	fn_I2 = parser.getOption("--i2", "Input name of half2, (default replaces half1 from --i with half2)", "");
	fn_OS = parser.getOption("--ios", "Input tomo optimiser set file. It is used to set --i if not provided. Updated output optimiser set is created.", "");
	fn_out = parser.getOption("--o", "Output rootname", "postprocess");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms", "-1"));
	write_halfmaps = parser.checkOption("--half_maps", "Write post-processed half maps for validation");
	mtf_angpix  = textToFloat(parser.getOption("--mtf_angpix", "Pixel size in the original micrographs/movies (in Angstroms)", "-1."));
	molweight = textToFloat(parser.getOption("--molweight", "Molecular weight (in kDa) of ordered protein mass", "-1"));

	int mask_section = parser.addSection("Masking options");
	do_auto_mask = parser.checkOption("--auto_mask", "Perform automated masking, based on a density threshold");
	ini_mask_density_threshold = textToFloat(parser.getOption("--inimask_threshold", "Density at which to threshold the map for the initial seed mask", "0.02"));
	extend_ini_mask = textToFloat(parser.getOption("--extend_inimask", "Number of pixels to extend the initial seed mask", "3."));
	width_soft_mask_edge  = textToFloat(parser.getOption("--width_mask_edge", "Width for the raised cosine soft mask edge (in pixels)", "6."));
	fn_mask = parser.getOption("--mask", "Filename of a user-provided mask (1=protein, 0=solvent, all values in range [0,1])", "");
	force_mask = parser.checkOption("--force_mask", "Use the mask even when the masked resolution is worse than the unmasked resolution");

	int sharp_section = parser.addSection("Sharpening options");
	fn_mtf = parser.getOption("--mtf", "User-provided STAR-file with the MTF-curve of the detector", "");
	do_auto_bfac = parser.checkOption("--auto_bfac", "Perform automated B-factor determination (Rosenthal and Henderson, 2003)");
	fit_minres = textToFloat(parser.getOption("--autob_lowres", "Lowest resolution (in A) to include in fitting of the B-factor", "10."));
	fit_maxres = textToFloat(parser.getOption("--autob_highres", "Highest resolution (in A) to include in fitting of the B-factor", "0."));
	adhoc_bfac =  textToFloat(parser.getOption("--adhoc_bfac", "User-provided B-factor (in A^2) for map sharpening, e.g. -400", "0."));

	int filter_section = parser.addSection("Filtering options");
	do_fsc_weighting = !parser.checkOption("--skip_fsc_weighting", "Do not use FSC-weighting (Rosenthal and Henderson, 2003) in the sharpening process");
	// include low-pass filter option in the program? This could be useful for structurally heterogeneous reconstructions (instead of FSC-weighting)
	low_pass_freq = textToFloat(parser.getOption("--low_pass", "Resolution (in Angstroms) at which to low-pass filter the final map (0: disable, negative: resolution at FSC=0.143)", "0"));

	int locres_section = parser.addSection("Local-resolution options");
	do_locres = parser.checkOption("--locres", "Perform local resolution estimation");
	locres_sampling = textToFloat(parser.getOption("--locres_sampling", "Sampling rate (in Angstroms) with which to sample the local-resolution map", "25."));
	locres_maskrad = textToFloat(parser.getOption("--locres_maskrad", "Radius (in A) of spherical mask for local-resolution map (default = 0.5*sampling)", "-1"));
	locres_edgwidth = textToFloat(parser.getOption("--locres_edgwidth", "Width of soft edge (in A) on masks for local-resolution map (default = sampling)", "-1"));
	locres_randomize_fsc = textToFloat(parser.getOption("--locres_randomize_at", "Randomize phases from this resolution (in A)", "25."));
	locres_minres = textToFloat(parser.getOption("--locres_minres", "Lowest local resolution allowed (in A)", "50."));

	int expert_section = parser.addSection("Expert options");
	do_ampl_corr = parser.checkOption("--ampl_corr", "Perform amplitude correlation and DPR, also re-normalize amplitudes for non-uniform angular distributions");
	randomize_fsc_at = textToFloat(parser.getOption("--randomize_at_fsc", "Randomize phases from the resolution where FSC drops below this value", "0.8"));
	randomize_at_A  = textToFloat(parser.getOption("--randomize_at_A", "Randomize phases from this resolution (in A) onwards (if positive)", "-1"));
	filter_edge_width = textToInteger(parser.getOption("--filter_edge_width", "Width of the raised cosine on the low-pass filter edge (in resolution shells)", "2"));
	do_interpolate = parser.checkOption("--interpolate", "Interpolate the FSC to obtain an additional, more precise resolution estimate");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));
	int random_seed = textToInteger(parser.getOption("--random_seed", "Seed for random number generator (negative value for truly random)", "0"));

	if (random_seed >= 0)
	{
		init_random_generator(random_seed);
	}

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void Postprocessing::usage()
{
	parser.writeUsage(std::cout);
}

void Postprocessing::clear()
{
	fn_I1 = fn_I2 = fn_OS = "";
	fn_out="postprocess";
	angpix = 1.;
	mtf_angpix = 1.;
	do_auto_mask = false;
	ini_mask_density_threshold = 0.02;
	width_soft_mask_edge = 6.;
	fn_mask = "";
	fn_mtf = "";
	do_auto_bfac = false;
	fit_minres = 10.;
	fit_maxres = 0.;
	adhoc_bfac = 0.;
	do_fsc_weighting = true;
	low_pass_freq = 0.;
	randomize_fsc_at = 0.8;
	randomize_at_A = -1.;
	filter_edge_width = 2.;
	verb = 1;
	do_ampl_corr = false;
}

void Postprocessing::initialise()
{
	// Check if input is a tomo optimisation set file
	if (fn_OS != "")
	{
		optimisationSet.read(fn_OS);

		if (fn_I1 == "")
		{
			if (!optimisationSet.getValue(EMDL_TOMO_REFERENCE_MAP_1_FILE_NAME, fn_I1) ||
				!optimisationSet.getValue(EMDL_TOMO_REFERENCE_MAP_2_FILE_NAME, fn_I2))
				REPORT_ERROR("No halfmap filenames were found in file " + fn_OS);
		}
		else
		{
			if (fn_I2 == "" && !fn_I1.getTheOtherHalf(fn_I2))
				REPORT_ERROR("The input filename does not contain 'half1' or 'half2'");

			optimisationSet.setValue(EMDL_TOMO_REFERENCE_MAP_1_FILE_NAME, fn_I1);
			optimisationSet.setValue(EMDL_TOMO_REFERENCE_MAP_2_FILE_NAME, fn_I2);
		}
	}
	else if (fn_I1 != "") // Read in the input maps
	{
		if (fn_I2 == "" && !fn_I1.getTheOtherHalf(fn_I2))
			REPORT_ERROR("The input filename does not contain 'half1' or 'half2'");
	}
	else
	{
		REPORT_ERROR("Input file is missing. Provide, at least, --i or --ios input.");
	}

	if (verb > 0)
	{
		std::cout <<"== Reading input half-reconstructions: " <<std::endl;
		std::cout.width(35); std::cout << std::left <<"  + half1-map: "; std::cout << fn_I1 << std::endl;
		std::cout.width(35); std::cout << std::left <<"  + half2-map: "; std::cout << fn_I2 << std::endl;
	}

	I1.read(fn_I1);
	I2.read(fn_I2);
	I1().setXmippOrigin();
	I2().setXmippOrigin();

	if (angpix <= 0)
	{
		angpix = I1.samplingRateX();
		std::cerr << "WARNING: You did not specify --angpix. The pixel size in the image header, " << angpix << " A/px, is used." << std::endl;
	}

	if (mtf_angpix < 0.)
	{
		if (verb > 0) std::cout << " + --mtf_angpix was not provided, assuming pixel size in raw micrographs is the same as in particles, " << angpix << " A/px." << std::endl;
		mtf_angpix = angpix;
	}

	// Calculate what fraction of voxels in the box is protein according to the expected ordered molecular weight
	// Protein density is 1.35 g/cm^3, Nav=6.022 E+23, so protein volume = (MW /0.81 Da) A^3
	// 47.6% is volume of sphere relative to box
	if (molweight > 0.)
	{
		frac_molweight = 0.476 * std::pow(XSIZE(I1())*angpix, 3) * 0.81 / ( molweight*1000) ;
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left   << "  + ordered molecular weight (kDa): "; std::cout  << molweight <<std::endl;
			std::cout.width(35); std::cout << std::left   << "  + fraction f (molweight based): "; std::cout  << frac_molweight <<std::endl;
		}
	}

	if (!I1().sameShape(I2()))
	{
		std::cerr << " Size of half1 map: "; I1().printShape(std::cerr); std::cerr << std::endl;
		std::cerr << " Size of half2 map: "; I2().printShape(std::cerr); std::cerr << std::endl;
		REPORT_ERROR("Postprocessing::initialise ERROR: The two half reconstructions are not of the same size!");
	}

	if (do_locres)
	{
		if (locres_maskrad < 0.0)
			locres_maskrad = 0.5*locres_sampling;
		if (locres_edgwidth < 0.0)
			locres_edgwidth = locres_sampling;

		if (fn_mask != "" && verb > 0)
			std::cerr << " WARNING: --mask is used only to make a histogram of local resolutions; it is not used for local resolution calculation itself." << std::endl;
		if (do_auto_bfac)
			REPORT_ERROR("Postprocessing::initialise ERROR: for --locres, you cannot do --auto_bfac, use --adhoc_bfac instead!");
	}

	if (do_auto_mask)
		REPORT_ERROR("Postprocessing:: --auto_mask has been removed. Please make a mask with relion_mask_create beforehand.");

	if (do_auto_bfac && ABS(adhoc_bfac) > 0.)
		REPORT_ERROR("Postprocessing::initialise ERROR: provide either --auto_bfac OR --adhoc_bfac, but not both!");
}

bool Postprocessing::getMask()
{
	if (fn_mask != "")
	{
		if (verb > 0)
		{
			std::cout << "== Using a user-provided mask ... " <<std::endl;
			std::cout.width(35); std::cout << std::left   << "  + input mask: "; std::cout  << fn_mask <<std::endl;
		}

		// Read the mask in memory
		Im.read(fn_mask);
		Im().setXmippOrigin();

		// Check values are between 0 and 1
		RFLOAT avg, stddev, minval, maxval;
		Im().computeStats(avg, stddev, minval, maxval);

		long summask = 0;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Im())
		{
			if (DIRECT_MULTIDIM_ELEM(Im(), n) > 0.5) summask++;
		}
		avg = (RFLOAT)summask / (RFLOAT)NZYXSIZE(Im());
		frac_solvent_mask = 0.476 /avg;
		molweight_frommask = avg * std::pow(XSIZE(Im()) * angpix, 3) * 0.81;

		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left   << "  + fraction f (solvent mask based): "; std::cout  << frac_solvent_mask <<std::endl;
			std::cout.width(35); std::cout << std::left   << "  + molecular weight inside protein mask: "; std::cout  << molweight_frommask <<std::endl;

		}

		if (minval < -1e-6 || maxval - 1. > 1.e-6)
		{
			std::cerr << " minval= " << minval << " maxval= " << maxval << std::endl;
			REPORT_ERROR("Postprocessing::mask ERROR: mask values not in range [0,1]!");
		}

		// Also check the mask is the same size as the input maps
		if (!Im().sameShape(I2()))
		{
			std::cerr << " Size of input mask: "; Im().printShape(std::cerr); std::cerr<< std::endl;
			std::cerr << " Size of input maps: "; I1().printShape(std::cerr); std::cerr<< std::endl;
			REPORT_ERROR("Postprocessing::mask ERROR: mask and input maps do not have the same size!");
		}

	}
	else
	{
		if (verb > 0)
		{
			std::cout << "== Not performing any masking ... " << std::endl;
			frac_solvent_mask = 0.;
		}
		return false;
	}

	return true;
}

void Postprocessing::divideByMtf(MultidimArray<Complex > &FT)
{
	if (fn_mtf != "")
	{
		if (verb > 0)
		{
			std::cout << "== Dividing map by the MTF of the detector ..." << std::endl;
			std::cout.width(35); std::cout << std::left <<"  + mtf STAR-file: "; std::cout << fn_mtf << std::endl;
		}

		MetaDataTable MDmtf;

		if (!fn_mtf.isStarFile())
			REPORT_ERROR("Postprocessing::divideByMtf ERROR: input MTF file is not a STAR file.");

		MDmtf.read(fn_mtf);
		MultidimArray<RFLOAT> mtf_resol, mtf_value;
		mtf_resol.resize(MDmtf.numberOfObjects());
		mtf_value.resize(mtf_resol);

		RFLOAT resol_inv_pixel;
		int i =0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDmtf)
		{
			MDmtf.getValue(EMDL_RESOLUTION_INVPIXEL, resol_inv_pixel);
			DIRECT_A1D_ELEM(mtf_resol, i) = resol_inv_pixel/mtf_angpix; // resolution needs to be given in 1/Ang
			MDmtf.getValue(EMDL_POSTPROCESS_MTF_VALUE, DIRECT_A1D_ELEM(mtf_value, i) );
			if (DIRECT_A1D_ELEM(mtf_value, i) < 1e-10)
			{
				std::cerr << " i= " << i <<  " mtf_value[i]= " << DIRECT_A1D_ELEM(mtf_value, i) << std::endl;
				REPORT_ERROR("Postprocessing::sharpenMap ERROR: zero or negative values encountered in MTF curve!");
			}
			i++;
		}

		// Calculate slope of resolution (in 1/A) per element in the MTF array, in order to interpolate below
		RFLOAT res_per_elem = (DIRECT_A1D_ELEM(mtf_resol, i-1) - DIRECT_A1D_ELEM(mtf_resol, 0)) / (RFLOAT)(i);
		if (res_per_elem < 1e-10) REPORT_ERROR(" ERROR: the resolution in the MTF star file does not go up....");

		RFLOAT xsize_ang = angpix * XSIZE(I1());
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
		{
			int r2 = kp * kp + ip * ip + jp * jp;
			RFLOAT res = sqrt((RFLOAT)r2)/xsize_ang; // get resolution in 1/Ang
			if (res < 1./(2.*angpix) )
			{
				int i_0 = FLOOR(res / res_per_elem);
				RFLOAT mtf;
				// check boundaries of the array
				if (i_0 >= MULTIDIM_SIZE(mtf_value) - 1)
					mtf = DIRECT_A1D_ELEM(mtf_value,  MULTIDIM_SIZE(mtf_value) - 1);
				else if (i_0 <= 0)
					mtf = DIRECT_A1D_ELEM(mtf_value, 0);
				else
				{
					// linear interpolation:
					RFLOAT x_0 = DIRECT_A1D_ELEM(mtf_resol, i_0);
					RFLOAT y_0 = DIRECT_A1D_ELEM(mtf_value, i_0);
					RFLOAT x_1 = DIRECT_A1D_ELEM(mtf_resol, i_0 + 1);
					RFLOAT y_1 = DIRECT_A1D_ELEM(mtf_value, i_0 + 1);
					mtf = y_0 + (y_1 - y_0)*(res - x_0)/(x_1 - x_0);
				}

				// Divide Fourier component by the MTF
				DIRECT_A3D_ELEM(FT, k, i, j) /= mtf;
			}
		}
	}
}

bool Postprocessing::findSurfacePixel(int idx, int kp, int ip, int jp,
		int &best_kpp, int &best_ipp, int &best_jpp,
		int myradius_count, int search)
{
	// bring kp, ip, jp onto the sphere
	RFLOAT frac = (RFLOAT)(myradius_count)/(RFLOAT)idx;
	int kpp = ROUND(frac*kp);
	int ipp = ROUND(frac*ip);
	int jpp = ROUND(frac*jp);

	// Search +/- 2 pixels in all directions and choose voxel closest to the circle
	int best_dist = 999;
	best_kpp=kpp;
	best_ipp=ipp;
	best_jpp=jpp;
	bool found = false;
	for (int kppp = kpp-search; kppp <= kpp+search; kppp++)
	{
		for (int ippp = ipp-search; ippp <= ipp+search; ippp++)
		{
			for (int jppp = jpp-search; jppp <= jpp+search; jppp++)
	 		{
				// Distance to surface on the sphere
				int dist = ABS(ROUND(sqrt((RFLOAT)(kppp*kppp + ippp*ippp + jppp*jppp))) - myradius_count);
				int reldist2 = (kppp-kpp)*(kppp-kpp) + (ippp-ipp)*(ippp-ipp) + (jppp-jpp)*(jppp-jpp);
				if (dist < 0.5 && reldist2 < best_dist)
				{
					best_kpp=kppp;
					best_ipp=ippp;
					best_jpp=jppp;
					best_dist=reldist2;
					found=true;
				}
			}
    		}
	}

	return found;
}

void Postprocessing::correctRadialAmplitudeDistribution(MultidimArray<RFLOAT > &I)
{
	MultidimArray<Complex > FT;
	FourierTransformer transformer;
	transformer.FourierTransform(I, FT, false);

	// First calculate radial average, to normalize the power spectrum
	int myradius = XSIZE(FT);
	MultidimArray< int > radial_count(myradius);
	MultidimArray<RFLOAT> num, ravg;
	num.initZeros(myradius);
	ravg.initZeros(myradius);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
	{
		int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
		if (idx >= myradius)
	    		continue;
		ravg(idx)+= norm(DIRECT_A3D_ELEM(FT, k, i, j));
		radial_count(idx)++;
	}
	FOR_ALL_ELEMENTS_IN_ARRAY1D(ravg)
	{
		if (radial_count(i) > 0)
		{
    			ravg(i) /= radial_count(i);
	    	}
	}

	// Apply correction only beyond low-res fitting of B-factors
	int minr = FLOOR(XSIZE(FT) * angpix / fit_minres);
	int myradius_count = minr;
	MultidimArray<RFLOAT> sum3d;
	MultidimArray<int> count3d;
	sum3d.resize(2*myradius_count+4, 2*myradius_count+4, 2*myradius_count+4);
	sum3d.setXmippOrigin();
	count3d.resize(sum3d);
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
	{
		int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
		// only correct from fit_minres to Nyquist
		if (idx < minr || idx >= myradius)
	     		continue;

		int best_kpp, best_ipp, best_jpp;
		if (!findSurfacePixel(idx, kp, ip, jp, best_kpp, best_ipp, best_jpp, myradius_count, 2))
		{
			std::cerr << "Postprocessing::correctRadialAmplitudeDistribution ERROR! kp= " << kp << " ip= " << ip << " jp= " << jp << std::endl;
		}

		// Apply correction on the spectrum-corrected values!
		RFLOAT aux = norm(DIRECT_A3D_ELEM(FT, k, i, j)) / ravg(idx);
		A3D_ELEM(sum3d, best_kpp, best_ipp, best_jpp) += aux;
		A3D_ELEM(count3d, best_kpp, best_ipp, best_jpp) += 1;
	}

	// Average
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(sum3d)
	{
		if (DIRECT_MULTIDIM_ELEM(count3d, n) > 0)
		{
			DIRECT_MULTIDIM_ELEM(sum3d, n) /= DIRECT_MULTIDIM_ELEM(count3d, n);
		}
	}

	// Now divide all elements by the normalized correction term
	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
	{
		int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
		// only correct from fit_minres to Nyquist
		if (idx < minr || idx >= myradius)
		     	continue;

		int best_kpp, best_ipp, best_jpp;
		if (!findSurfacePixel(idx, kp, ip, jp, best_kpp, best_ipp, best_jpp, myradius_count, 2))
		{
			std::cerr << "Postprocessing::correctRadialAmplitudeDistribution ERROR!  kp= " << kp << " ip= " << ip << " jp= " << jp << std::endl;
		}

		// Apply correction on the spectrum-corrected values!
		RFLOAT aux = sqrt(A3D_ELEM(sum3d, best_kpp, best_ipp, best_jpp));
		DIRECT_A3D_ELEM(FT, k, i, j) /= aux;
	}

	transformer.inverseFourierTransform(FT, I);
}

RFLOAT Postprocessing::sharpenMap()
{
	MultidimArray<Complex > FT;
	FourierTransformer transformer;
	transformer.FourierTransform(I1(), FT, true);

	makeGuinierPlot(FT, guinierin);

	// A. If MTF curve is given, first divide by the MTF
	divideByMtf(FT);
	makeGuinierPlot(FT, guinierinvmtf);

	// B. Then perform B-factor sharpening
	if (do_fsc_weighting)
	{
		if (verb > 0)
		{
			std::cout <<"== Applying sqrt(2*FSC/(FSC+1)) weighting (as in Rosenthal & Henderson, 2003) ..." <<std::endl;
		}
		applyFscWeighting(FT, fsc_true);
	}
	makeGuinierPlot(FT, guinierweighted);

	global_bfactor = 0.;
	if (do_auto_bfac)
	{
		if (verb > 0)
		{
			std::cout <<"== Fitting straight line through Guinier plot to find B-factor ..." <<std::endl;
			std::cout.width(35); std::cout << std::left <<"  + fit from resolution: "; std::cout << fit_minres << std::endl;
			std::cout.width(35); std::cout << std::left  <<"  + fit until resolution: "; std::cout << fit_maxres << std::endl;
		}

		fitStraightLine(guinierweighted, global_slope, global_intercept, global_corr_coeff);
		global_bfactor = 4. * global_slope;
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left  <<"  + slope of fit: "; std::cout << global_slope << std::endl;
			std::cout.width(35); std::cout << std::left  <<"  + intercept of fit: "; std::cout << global_intercept << std::endl;
			std::cout.width(35); std::cout << std::left  <<"  + correlation of fit: "; std::cout << global_corr_coeff << std::endl;
		}
	}
	else if (ABS(adhoc_bfac) > 0.)
	{
		if (verb > 0)
		{
			std::cout <<"== Using a user-provided (ad-hoc) B-factor ..." <<std::endl;
		}
		if (adhoc_bfac > 0.)
			std::cout <<" WARNING: using a positive B-factor. This will effectively dampen your map. Use negative value to sharpen it!" << std::endl;
		global_bfactor = adhoc_bfac;
	}

	// Now apply the B-factor
	if (ABS(global_bfactor) > 0.)
	{
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left  <<"  + apply b-factor of: "; std::cout << global_bfactor << std::endl;
		}
		applyBFactorToMap(FT, XSIZE(I1()), global_bfactor, angpix);
	}
        else
        {
            // Write a warning that the map will not be sharpened
            std::cerr << " WARNING: map will not be sharpened. This map may be used as reference for subsequent refinements, but you would need to run a postprocessing job with B-factor sharpening to better visualise high-resolution features." << std::endl;
        }

	makeGuinierPlot(FT, guiniersharpen);

	RFLOAT applied_filter = low_pass_freq;
	if (low_pass_freq != 0)
	{
		if (low_pass_freq < 0)
			applied_filter = global_resol;

		if (verb > 0)
		{
			std::cout << "== Low-pass filtering final map ... " << std::endl;
			std::cout.width(35); std::cout << std::left  <<"  + filter frequency: "; std::cout << applied_filter << std::endl;
		}

		lowPassFilterMap(FT, XSIZE(I1()), applied_filter, angpix, filter_edge_width);
	}

	transformer.inverseFourierTransform(FT, I1());

	return applied_filter;
}

void Postprocessing::calculateFSCtrue(MultidimArray<RFLOAT> &fsc_true, MultidimArray<RFLOAT> &fsc_unmasked,
		MultidimArray<RFLOAT> &fsc_masked, MultidimArray<RFLOAT> &fsc_random_masked, int randomize_at)
{
	// Now that we have fsc_masked and fsc_random_masked, calculate fsc_true according to Richard's formula
	// FSC_true = FSC_t - FSC_n / ( )

	// Sometimes FSc at origin becomes -1!
	if (DIRECT_A1D_ELEM(fsc_masked, 0) <= 0.)
		DIRECT_A1D_ELEM(fsc_masked, 0) = 1.;
	if (DIRECT_A1D_ELEM(fsc_random_masked, 0) <= 0.)
		DIRECT_A1D_ELEM(fsc_random_masked, 0) = 1.;


	fsc_true.resize(fsc_masked);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
	{
		// 29jan2015: let's move this 2 shells upwards, because of small artefacts near the resolution of randomisation!
		if (i < randomize_at + 2)
		{
			DIRECT_A1D_ELEM(fsc_true, i) = DIRECT_A1D_ELEM(fsc_masked, i);
		}
		else
		{
			RFLOAT fsct = DIRECT_A1D_ELEM(fsc_masked, i);
			RFLOAT fscn = DIRECT_A1D_ELEM(fsc_random_masked, i);
			DIRECT_A1D_ELEM(fsc_true, i) = (fsct - fscn) / (1. - fscn);
		}
	}
}

void Postprocessing::calculateFSCpart(const MultidimArray<RFLOAT> fsc_unmasked, RFLOAT fraction, MultidimArray<RFLOAT> &fsc_part)
{
	// Now that we have fsc_masked and fsc_random_masked, calculate fsc_true according to Richard's formula
	// FSC_true = FSC_t - FSC_n / ( )

	// Sometimes FSc at origin becomes -1!
	if (DIRECT_A1D_ELEM(fsc_masked, 0) <= 0.)
		DIRECT_A1D_ELEM(fsc_masked, 0) = 1.;
	if (DIRECT_A1D_ELEM(fsc_random_masked, 0) <= 0.)
		DIRECT_A1D_ELEM(fsc_random_masked, 0) = 1.;


	fsc_part.resize(fsc_unmasked);
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_part)
	{
		DIRECT_A1D_ELEM(fsc_part, i) = fraction * DIRECT_A1D_ELEM(fsc_unmasked, i) / (1. + (fraction-1) * DIRECT_A1D_ELEM(fsc_unmasked, i));
	}

}

void Postprocessing::applyFscWeighting(MultidimArray<Complex > &FT, MultidimArray<RFLOAT> my_fsc)
{
	// Find resolution where fsc_true drops below zero for the first time
	// Set all weights to zero beyond that resolution
	int ires_max = 0 ;

	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(my_fsc)
	{
		if (DIRECT_A1D_ELEM(my_fsc, i) < 0.0001)
			break;
		ires_max = i;
	}

	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
	{
		int ires = ROUND(sqrt((RFLOAT)kp * kp + ip * ip + jp * jp));
		if (ires <= ires_max)
		{
			RFLOAT fsc = DIRECT_A1D_ELEM(my_fsc, ires);
			if (fsc > 0.)
				DIRECT_A3D_ELEM(FT, k, i, j) *= sqrt((2 * fsc) / (1 + fsc));
			else
				DIRECT_A3D_ELEM(FT, k, i, j) *= 0.;
		}
		else
		{
			DIRECT_A3D_ELEM(FT, k, i, j) = 0.;
		}
	}
}

void Postprocessing::makeGuinierPlot(MultidimArray<Complex > &FT, std::vector<fit_point2D> &guinier)
{
	MultidimArray<int> radial_count(XSIZE(FT));
	MultidimArray<RFLOAT> lnF(XSIZE(FT));
	fit_point2D      onepoint;

	FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
	{
		int r2 = kp * kp + ip * ip + jp * jp;
		int ires = ROUND(sqrt((RFLOAT)r2));
		if (ires < XSIZE(radial_count))
		{
		        lnF(ires) += abs(DIRECT_A3D_ELEM(FT, k, i, j));
		        radial_count(ires)++;
		}
	}

	RFLOAT xsize = XSIZE(I1());
	guinier.clear();
	FOR_ALL_ELEMENTS_IN_ARRAY1D(radial_count)
	{

		RFLOAT res = (xsize * angpix)/(RFLOAT)i; // resolution in Angstrom
		if (res >= angpix * 2.) // Apply B-factor sharpening until Nyquist, then low-pass filter later on (with a soft edge)
		{
			onepoint.x = 1. / (res * res);
			if (DIRECT_A1D_ELEM(lnF, i) > 0.)
			{
				onepoint.y = log ( DIRECT_A1D_ELEM(lnF, i) / DIRECT_A1D_ELEM(radial_count, i) );
				if (res <= fit_minres && res >= fit_maxres)
				{
					onepoint.w = 1.;
				}
				else
				{
					onepoint.w = 0.;
				}
			}
			else
			{
				onepoint.y = -99.;
				onepoint.w = 0.;
			}
			//std::cerr << " onepoint.x= " << onepoint.x << " onepoint.y= " << onepoint.y << " onepoint.w= " << onepoint.w << std::endl;
			guinier.push_back(onepoint);
		}
	}
}

void Postprocessing::writeOutput()
{
	FileName fn_tmp;

	if (verb > 0)
	{
		std::cout <<"== Writing output files ..." <<std::endl;
	}

	writeMaps(fn_out);

	// Write an output STAR file with FSC curves, Guinier plots etc
	std::ofstream  fh;
	fn_tmp = fn_out + ".star";
	if (verb > 0)
	{
		std::cout.width(35); std::cout << std::left <<"  + Metadata file: "; std::cout << fn_tmp<< std::endl;
	}

	fh.open((fn_tmp).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)"MlOptimiser::write: Cannot write file: " + fn_tmp);

	// Write the command line as a comment in the header
	fh << "# RELION postprocess; version " << g_RELION_VERSION << std::endl;
	fh << "# ";
	parser.writeCommandLine(fh);

	MetaDataTable MDlist, MDfsc, MDguinier;

	MDlist.setIsList(true);
	MDlist.setName("general");
	MDlist.addObject();
	MDlist.setValue(EMDL_POSTPROCESS_FINAL_RESOLUTION, global_resol);
	MDlist.setValue(EMDL_POSTPROCESS_BFACTOR, global_bfactor );
	MDlist.setValue(EMDL_POSTPROCESS_UNFIL_HALFMAP1, fn_I1);
	MDlist.setValue(EMDL_POSTPROCESS_UNFIL_HALFMAP2, fn_I2);
	MDlist.setValue(EMDL_POSTPROCESSED_MAP, fn_out + ".mrc");
	if (molweight > 0.)
	{
		MDlist.setValue(EMDL_POSTPROCESS_MOLWEIGHT, molweight);
		MDlist.setValue(EMDL_POSTPROCESS_FRACTION_MOLWEIGHT, frac_molweight);
	}
	if (do_mask)
	{
		MDlist.setValue(EMDL_POSTPROCESS_FRACTION_SOLVENT_MASK, frac_solvent_mask);
		RFLOAT randomize_at_Ang = XSIZE(I1())* angpix / randomize_at;
		MDlist.setValue(EMDL_MASK_NAME, fn_mask);
		MDlist.setValue(EMDL_POSTPROCESS_RANDOMISE_FROM, randomize_at_Ang);
		MDlist.setValue(EMDL_POSTPROCESSED_MAP_MASKED, fn_out + "_masked.mrc");
	}
	if (do_auto_bfac)
	{
		MDlist.setValue(EMDL_POSTPROCESS_GUINIER_FIT_SLOPE, global_slope);
		MDlist.setValue(EMDL_POSTPROCESS_GUINIER_FIT_INTERCEPT, global_intercept);
		MDlist.setValue(EMDL_POSTPROCESS_GUINIER_FIT_CORRELATION, global_corr_coeff);
	}
	MDlist.write(fh);

	// If input optimisation set is provided, also crete it as output
	if (!optimisationSet.isEmpty())
	{
		optimisationSet.setValue(EMDL_TOMO_REFERENCE_FSC_FILE_NAME, fn_tmp);
		optimisationSet.write(fn_out + "_optimisation_set.star");
	}

	MDfsc.setName("fsc");
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
	{
		MDfsc.addObject();
		RFLOAT res = (i > 0) ? (XSIZE(I1()) * angpix / (RFLOAT)i) : 999.;
		MDfsc.setValue(EMDL_SPECTRAL_IDX, (int)i);
		MDfsc.setValue(EMDL_RESOLUTION, 1./res);
		MDfsc.setValue(EMDL_RESOLUTION_ANGSTROM, res);
		if (do_mask)
		{
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_TRUE, DIRECT_A1D_ELEM(fsc_true, i) );
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_PART_FRACMASK, DIRECT_A1D_ELEM(fsc_part_fracmask, i) );
			if (molweight > 0.)
				MDfsc.setValue(EMDL_POSTPROCESS_FSC_PART_MOLWEIGHT, DIRECT_A1D_ELEM(fsc_part_molweight, i) );
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_UNMASKED, DIRECT_A1D_ELEM(fsc_unmasked, i) );
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_MASKED, DIRECT_A1D_ELEM(fsc_masked, i) );
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_RANDOM_MASKED, DIRECT_A1D_ELEM(fsc_random_masked, i) );
			if (do_ampl_corr)
			{
				MDfsc.setValue(EMDL_POSTPROCESS_AMPLCORR_UNMASKED, DIRECT_A1D_ELEM(acorr_unmasked, i) );
				MDfsc.setValue(EMDL_POSTPROCESS_AMPLCORR_MASKED, DIRECT_A1D_ELEM(acorr_masked, i) );
				MDfsc.setValue(EMDL_POSTPROCESS_DPR_UNMASKED, DIRECT_A1D_ELEM(dpr_unmasked, i) );
				MDfsc.setValue(EMDL_POSTPROCESS_DPR_MASKED, DIRECT_A1D_ELEM(dpr_masked, i) );
			}
		}
		else
		{
			MDfsc.setValue(EMDL_POSTPROCESS_FSC_UNMASKED, DIRECT_A1D_ELEM(fsc_true, i) );
			if (molweight > 0.)
				MDfsc.setValue(EMDL_POSTPROCESS_FSC_PART_MOLWEIGHT, DIRECT_A1D_ELEM(fsc_part_molweight, i) );
			if (do_ampl_corr)
			{
				MDfsc.setValue(EMDL_POSTPROCESS_AMPLCORR_UNMASKED, DIRECT_A1D_ELEM(acorr_unmasked, i) );
				MDfsc.setValue(EMDL_POSTPROCESS_DPR_UNMASKED, DIRECT_A1D_ELEM(dpr_unmasked, i) );
			}
		}
	}
	MDfsc.write(fh);

	// Write a plot with the FSC curves
	std::string title= "Final resolution = " + floatToString(global_resol, 5, 2) + " Angstroms";
	CPlot2D *plot2D = new CPlot2D(title);
	plot2D->SetXAxisSize(600);
	plot2D->SetYAxisSize(400);
	MDfsc.addToCPlot2D(plot2D, EMDL_RESOLUTION, EMDL_POSTPROCESS_FSC_TRUE, 0., 0., 0., 2.);
	MDfsc.addToCPlot2D(plot2D, EMDL_RESOLUTION, EMDL_POSTPROCESS_FSC_UNMASKED, 0., 1., 0.);
	MDfsc.addToCPlot2D(plot2D, EMDL_RESOLUTION, EMDL_POSTPROCESS_FSC_MASKED, 0., 0., 1.);
	MDfsc.addToCPlot2D(plot2D, EMDL_RESOLUTION, EMDL_POSTPROCESS_FSC_RANDOM_MASKED, 1., 0., 0.);
	plot2D->SetXAxisTitle("resolution (1/A)");
	plot2D->SetYAxisTitle("Fourier Shell Correlation");
	plot2D->OutputPostScriptPlot(fn_out + "_fsc.eps");
	delete plot2D;

//#define CISTEMFSC
#ifdef CISTEMFSC
	// Write a plot with the FSC curves
	std::string title2= "RELION/cisTEM FSC comparison; MW_mask = " +  floatToString(molweight_frommask/1000., 8,2) + " kDa";
	CPlot2D *plot2Db = new CPlot2D(title2);
	plot2Db->SetXAxisSize(600);
	plot2Db->SetYAxisSize(400);
	MDfsc.addToCPlot2D(plot2Db, EMDL_RESOLUTION, EMDL_POSTPROCESS_FSC_TRUE, 0., 0., 0., 2.);
	MDfsc.addToCPlot2D(plot2Db, EMDL_RESOLUTION, EMDL_POSTPROCESS_FSC_PART_FRACMASK, 1., 0.66, 0., 1.);
	if (molweight > 0.)
		MDfsc.addToCPlot2D(plot2Db, EMDL_RESOLUTION, EMDL_POSTPROCESS_FSC_PART_MOLWEIGHT, 0., 1., 1., 2.);
	plot2Db->SetXAxisTitle("resolution (1/A)");
	plot2Db->SetYAxisTitle("Fourier Shell Correlation");
	plot2Db->OutputPostScriptPlot(fn_out + "_fsc_part.eps");
	delete plot2Db;
#endif

	// Also write XML file with FSC_true curve for EMDB submission
	writeFscXml(MDfsc);

	// Write a DAT file for easier plotting in xmgrace (JZ, June 2nd, 2020)
	writeFscDat(MDfsc);

	MDguinier.setName("guinier");
	MetaDataTable MDextra1, MDextra2; // for postscript plot
	for (int i = 0; i < guinierin.size(); i++)
	{
		MDguinier.addObject();
		MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, guinierin[i].x);
		MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_IN, guinierin[i].y);
		if (fn_mtf != "")
			MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_INVMTF, guinierinvmtf[i].y);
		if (do_fsc_weighting)
		{
			MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_WEIGHTED, guinierweighted[i].y);
			if (guinierweighted[i].y > -99.)
			{
				MDextra1.addObject();
				MDextra1.setValue(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, guinierin[i].x);
				MDextra1.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_WEIGHTED, guinierweighted[i].y);
			}
		}
		if (do_auto_bfac || ABS(adhoc_bfac) > 0.)
		{
			MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_SHARPENED, guiniersharpen[i].y);
			if (guiniersharpen[i].y > -99.)
			{
				MDextra2.addObject();
				MDextra2.setValue(EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, guinierin[i].x);
				MDextra2.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_SHARPENED, guiniersharpen[i].y);
			}
		}
		if (do_auto_bfac)
			MDguinier.setValue(EMDL_POSTPROCESS_GUINIER_VALUE_INTERCEPT, global_intercept);
	}
	MDguinier.write(fh);
	fh.close();

	CPlot2D *plot2Dc = new CPlot2D("Guinier plots");
	plot2Dc->SetXAxisSize(600);
	plot2Dc->SetYAxisSize(400);
	MDguinier.addToCPlot2D(plot2Dc, EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, EMDL_POSTPROCESS_GUINIER_VALUE_IN, 0., 0., 0.);
	if (fn_mtf != "")
		MDguinier.addToCPlot2D(plot2Dc, EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, EMDL_POSTPROCESS_GUINIER_VALUE_INVMTF, 0., 1., 0.);
	if (do_fsc_weighting)
	{
		MDextra1.addToCPlot2D(plot2Dc, EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, EMDL_POSTPROCESS_GUINIER_VALUE_WEIGHTED, 0., 0., 1.);
	}
	if (do_auto_bfac || ABS(adhoc_bfac) > 0.)
	{
		MDextra2.addToCPlot2D(plot2Dc, EMDL_POSTPROCESS_GUINIER_RESOL_SQUARED, EMDL_POSTPROCESS_GUINIER_VALUE_SHARPENED, 1., 0., 0.);
	}
	plot2Dc->SetXAxisTitle("resolution^2 (1/A^2)");
	plot2Dc->SetYAxisTitle("ln(amplitudes)");
	plot2Dc->OutputPostScriptPlot(fn_out + "_guinier.eps");
	delete plot2Dc;

	FileName fn_log = fn_out.beforeLastOf("/") + "/logfile.pdf";
	if (!exists(fn_log))
	{
		std::vector<FileName> fn_eps;
		fn_eps.push_back(fn_out + "_fsc.eps");
		fn_eps.push_back(fn_out + "_fsc_part.eps");
		fn_eps.push_back(fn_out + "_guinier.eps");
		joinMultipleEPSIntoSinglePDF(fn_log, fn_eps);
	}

	if (verb > 0)
	{
		if (do_interpolate)
		{
			std::cout.width(35);
			std::cout << std::left   <<"  + FINAL RESOLUTION: ";
			std::cout << global_resol << " (" << fract_resol << ')' << std::endl;
		}
		else
		{
			std::cout.width(35); std::cout << std::left   <<"  + FINAL RESOLUTION: "; std::cout << global_resol<< std::endl;
		}
	}
}

// This masks I1!
void Postprocessing::writeMaps(FileName fn_root) {
	FileName fn_tmp = fn_root + ".mrc";

	I1.setStatisticsInHeader();
	I1.setSamplingRateInHeader(angpix);
	I1.write(fn_tmp);
	if (verb > 0)
	{
		std::cout.width(35); std::cout << std::left   <<"  + Processed map: "; std::cout << fn_tmp<< std::endl;
	}

	// Also write the masked postprocessed map
	if (fn_mask != "")
	{
		fn_tmp = fn_root + "_masked.mrc";
		I1() *= Im();
		I1.setStatisticsInHeader();
		I1.write(fn_tmp);
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left   <<"  + Processed masked map: "; std::cout << fn_tmp<< std::endl;
		}
	}
}

void Postprocessing::writeFscXml(MetaDataTable &MDfsc)
{
	FileName fn_fsc = fn_out + "_fsc.xml";
	std::ofstream  fh;
	fh.open((fn_fsc).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)"MetaDataTable::write Cannot write to file: " + fn_fsc);

	fh << "<fsc title=\"RELION masked-corrected FSC\" xaxis=\"Resolution (A-1)\" yaxis=\"Correlation Coefficient\">"<<std::endl;

	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDfsc)
	{
		RFLOAT xx, yy;
		MDfsc.getValue(EMDL_RESOLUTION, xx);
		MDfsc.getValue(EMDL_POSTPROCESS_FSC_TRUE, yy);
		fh << "  <coordinate>" << std::endl;
		fh << "    <x>" << xx << "</x>" << std::endl;
		fh << "    <y>" << yy << "</y>" << std::endl;
		fh << "  </coordinate>" << std::endl;
	}
	fh << "</fsc>" << std::endl;
	fh.close();
}

void Postprocessing::writeFscDat(MetaDataTable &MDfsc)
{
	const std::string fn_fsc = fn_out + "_fsc.dat";

	std::ofstream fh(fn_fsc);

	if (!fh)
	{
		REPORT_ERROR( (std::string)"MetaDataTable::write Cannot write to file: " + fn_fsc);
	}

	for (int i = 0; i < MDfsc.numberOfObjects(); i++)
	{
		const RFLOAT xx = MDfsc.getRfloat(EMDL_RESOLUTION, i);
		const RFLOAT yy = MDfsc.getRfloat(EMDL_POSTPROCESS_FSC_TRUE, i);

		fh << xx << ' ' << yy << '\n';
	}
}

void Postprocessing::run_locres(int rank, int size)
{
	// Read input maps and perform some checks
	initialise();

	// Also read the user-provided mask
	//getMask();

	MultidimArray<RFLOAT> I1m, I2m, I1p, I2p, Isum, locmask, Ilocres, Ifil, Isumw;

	// Get sum of two half-maps and sharpen according to estimated or ad-hoc B-factor
	Isum.resize(I1());
	I1m.resize(I1());
	I2m.resize(I1());
	I1p.resize(I1());
	I2p.resize(I1());
	locmask.resize(I1());
	// Initialise local-resolution maps, weights etc
	Ifil.initZeros(I1());
	Ilocres.initZeros(I1());
	Isumw.initZeros(I1());
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I1())
	{
		DIRECT_MULTIDIM_ELEM(Isum, n) = DIRECT_MULTIDIM_ELEM(I1(), n) + DIRECT_MULTIDIM_ELEM(I2(), n);
		DIRECT_MULTIDIM_ELEM(I1p, n) = DIRECT_MULTIDIM_ELEM(I1(), n);
		DIRECT_MULTIDIM_ELEM(I2p, n) = DIRECT_MULTIDIM_ELEM(I2(), n);
	}

	// Pre-sharpen the sum of the two half-maps with the provided MTF curve and adhoc B-factor
	do_fsc_weighting = false;
	MultidimArray<Complex > FTsum;
	FourierTransformer transformer;
	transformer.FourierTransform(Isum, FTsum, true);
	divideByMtf(FTsum);
	applyBFactorToMap(FTsum, XSIZE(Isum), adhoc_bfac, angpix);

	// Step size of locres-sampling in pixels
	int step_size = ROUND(locres_sampling / angpix);
	int maskrad_pix = ROUND(locres_maskrad / angpix);
	int edgewidth_pix = ROUND(locres_edgwidth / angpix);

	// Get the unmasked FSC curve
	getFSC(I1(), I2(), fsc_unmasked);

	// Randomize phases of unmasked maps from user-provided resolution
	int randomize_at = XSIZE(I1())* angpix / locres_randomize_fsc;
	if (verb > 0)
	{
		std::cout.width(35); std::cout << std::left << "  + randomize phases beyond: "; std::cout << XSIZE(I1())* angpix / randomize_at << " Angstroms" << std::endl;
	}
	// Randomize phases
	randomizePhasesBeyond(I1p, randomize_at);
	randomizePhasesBeyond(I2p, randomize_at);

	// Write an output STAR file with FSC curves, Guinier plots etc
	FileName fn_tmp = fn_out + "_locres_fscs.star";
	std::ofstream  fh;
	if (rank == 0)
	{
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left <<"  + Metadata output file: "; std::cout << fn_tmp<< std::endl;
		}

		fh.open((fn_tmp).c_str(), std::ios::out);
		if (!fh)
			REPORT_ERROR( (std::string)"MlOptimiser::write: Cannot write file: " + fn_tmp);
	}

	// Sample the entire volume (within the provided mask)

	int myrad = XSIZE(I1())/2 - maskrad_pix;
	float myradf = (float)myrad/(float)step_size;
	long int nr_samplings = ROUND((4.* PI / 3.) * (myradf*myradf*myradf));
	if (verb > 0)
	{
		std::cout << " Calculating local resolution in " << nr_samplings << " sampling points ..." << std::endl;
		init_progress_bar(nr_samplings);
	}

	long int nn = 0;
	for (long int kk=((I1()).zinit); kk<=((I1()).zinit + (I1()).zdim - 1); kk+= step_size)
	{
		for (long int ii=((I1()).yinit); ii<=((I1()).yinit + (I1()).ydim - 1); ii+= step_size)
		{
			for (long int jj=((I1()).xinit); jj<=((I1()).xinit + (I1()).xdim - 1); jj+= step_size)
			{
				// Abort through the pipeline_control system, TODO: check how this goes with MPI....
				if (pipeline_control_check_abort_job())
					exit(RELION_EXIT_ABORTED);

				// Only calculate local-resolution inside a spherical mask with radius less than half-box-size minus maskrad_pix
				float rad = sqrt(kk*kk + ii*ii + jj*jj);
				if (rad < myrad)
				{
					if (nn%size == rank)
					{
						// Make a spherical mask around (k,i,j), diameter is step_size pixels, soft-edge width is edgewidth_pix
						raisedCosineMask(locmask, maskrad_pix, maskrad_pix + edgewidth_pix, kk, ii, jj);

						// FSC of masked maps
						I1m = I1() * locmask;
						I2m = I2() * locmask;
						getFSC(I1m, I2m, fsc_masked);

						// FSC of masked randomized-phase map
						I1m = I1p * locmask;
						I2m = I2p * locmask;
						getFSC(I1m, I2m, fsc_random_masked);

						// Now that we have fsc_masked and fsc_random_masked, calculate fsc_true according to Richard's formula
						// FSC_true = FSC_t - FSC_n / ( )
						calculateFSCtrue(fsc_true, fsc_unmasked, fsc_masked, fsc_random_masked, randomize_at);

						if (rank == 0)
						{
							MetaDataTable MDfsc;
							FileName fn_name = "fsc_"+integerToString(kk, 5)+"_"+integerToString(ii, 5)+"_"+integerToString(jj, 5);
							MDfsc.setName(fn_name);
							FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
							{
								MDfsc.addObject();
								RFLOAT res = (i > 0) ? (XSIZE(I1()) * angpix / (RFLOAT)i) : 999.;
								MDfsc.setValue(EMDL_SPECTRAL_IDX, (int)i);
								MDfsc.setValue(EMDL_RESOLUTION, 1./res);
								MDfsc.setValue(EMDL_RESOLUTION_ANGSTROM, res);
								MDfsc.setValue(EMDL_POSTPROCESS_FSC_TRUE, DIRECT_A1D_ELEM(fsc_true, i) );
								MDfsc.setValue(EMDL_POSTPROCESS_FSC_UNMASKED, DIRECT_A1D_ELEM(fsc_unmasked, i) );
								MDfsc.setValue(EMDL_POSTPROCESS_FSC_MASKED, DIRECT_A1D_ELEM(fsc_masked, i) );
								MDfsc.setValue(EMDL_POSTPROCESS_FSC_RANDOM_MASKED, DIRECT_A1D_ELEM(fsc_random_masked, i) );
							}
							MDfsc.write(fh);
						}

						float local_resol = 999.;
						// See where corrected FSC drops below 0.143
						FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
						{
							if ( DIRECT_A1D_ELEM(fsc_true, i) < 0.143)
								break;
							local_resol = (i > 0) ? XSIZE(I1())*angpix/(RFLOAT)i : 999.;
						}
						local_resol = XMIPP_MIN(locres_minres, local_resol);
						if (rank == 0)
							fh << " kk= " << kk << " ii= " << ii << " jj= " << jj << " local resolution= " << local_resol << std::endl;

						// Now low-pass filter Isum to the estimated resolution
						MultidimArray<Complex > FT = FTsum;
						applyFscWeighting(FT, fsc_true);
						lowPassFilterMap(FT, XSIZE(I1()), local_resol, angpix, filter_edge_width);

						// Re-use I1m to save some memory
						transformer.inverseFourierTransform(FT, I1m);

						// Store weighted sum of local resolution and filtered map
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I1m)
						{
							DIRECT_MULTIDIM_ELEM(Ifil, n) +=  DIRECT_MULTIDIM_ELEM(locmask, n) * DIRECT_MULTIDIM_ELEM(I1m, n);
							DIRECT_MULTIDIM_ELEM(Ilocres, n) +=  DIRECT_MULTIDIM_ELEM(locmask, n) / local_resol;
							DIRECT_MULTIDIM_ELEM(Isumw, n) +=  DIRECT_MULTIDIM_ELEM(locmask, n);
						}
					}

					nn++;
					if (verb > 0 && nn <= nr_samplings)
						progress_bar(nn);
				}
			}
		}
	}

	fh.close();
	if (verb > 0)
		init_progress_bar(nr_samplings);


	if (size > 1)
	{
		I1m.initZeros();
		MPI_Allreduce(MULTIDIM_ARRAY(Ifil), MULTIDIM_ARRAY(I1m), MULTIDIM_SIZE(Ifil), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		Ifil = I1m;
		I1m.initZeros();
		MPI_Allreduce(MULTIDIM_ARRAY(Ilocres), MULTIDIM_ARRAY(I1m), MULTIDIM_SIZE(Ilocres), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		Ilocres = I1m;
		I1m.initZeros();
		MPI_Allreduce(MULTIDIM_ARRAY(Isumw), MULTIDIM_ARRAY(I1m), MULTIDIM_SIZE(Isumw), MY_MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
		Isumw = I1m;
	}

	if (rank == 0)
	{
		// Now write out the local-resolution map and
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I1m)
		{
			if (DIRECT_MULTIDIM_ELEM(Isumw, n ) > 0.)
			{
				DIRECT_MULTIDIM_ELEM(I1(), n) = 1. / (DIRECT_MULTIDIM_ELEM(Ilocres, n) / DIRECT_MULTIDIM_ELEM(Isumw, n));
				DIRECT_MULTIDIM_ELEM(I2(), n) = DIRECT_MULTIDIM_ELEM(Ifil, n) / DIRECT_MULTIDIM_ELEM(Isumw, n);
			}
			else
			{
				DIRECT_MULTIDIM_ELEM(I1(), n) = 0.;
				DIRECT_MULTIDIM_ELEM(I2(), n) = 0.;
			}
		}

		fn_tmp = fn_out + "_locres.mrc";
		I1.setSamplingRateInHeader(angpix);
		I1.write(fn_tmp);
		fn_tmp = fn_out + "_locres_filtered.mrc";
		I2.setSamplingRateInHeader(angpix);
		I2.write(fn_tmp);

#ifdef DEBUG
		I1() = Isumw;
		fn_tmp = fn_out + "_locres_sumw.mrc";
		I1.write(fn_tmp);
#endif

		if (fn_mask != "")
		{
			std::cout << "Calculating a histogram of local resolutions within the mask." << std::endl;

			std::vector<RFLOAT> values;
			Image<RFLOAT> Imask;
			Imask.read(fn_mask);

			if (!I1().sameShape(Imask(), true))
				REPORT_ERROR("The shape of the input half maps and the mask is not the same.");

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imask())
				if (DIRECT_MULTIDIM_ELEM(Imask(), n) > 0.5)
					values.push_back(DIRECT_MULTIDIM_ELEM(I1(), n));

			std::vector <RFLOAT> histX, histY;
			CPlot2D *plot2D=new CPlot2D("");
                	FileName fn_eps = fn_out + "_histogram.eps";
			MetaDataTable::histogram(values, histX, histY, verb, "local resolution", plot2D);
			plot2D->OutputPostScriptPlot(fn_eps);
			FileName fn_log = fn_out.beforeLastOf("/") + "/histogram.pdf";
			std::vector<FileName> to_convert;
			to_convert.push_back(fn_eps);
			joinMultipleEPSIntoSinglePDF(fn_log, to_convert);
			std::cout << "Written the histogram to " << fn_log << std::endl;
		}
	}

	if (verb > 0)
		std::cout << " done! " << std::endl;

	if (size > 1)
		MPI_Barrier(MPI_COMM_WORLD);
}

void Postprocessing::run()
{
	// Read input maps and perform some checks
	initialise();

	// For amplitude correlation curves: first do radial amplitude correction for non-uniform angular distributions
	if (do_ampl_corr)
	{
		correctRadialAmplitudeDistribution(I1());
		correctRadialAmplitudeDistribution(I2());
	}

	// Calculate FSC of the unmask maps
	getFSC(I1(), I2(), fsc_unmasked);
	if (do_ampl_corr)
		getAmplitudeCorrelationAndDifferentialPhaseResidual(I1(), I2(), acorr_unmasked, dpr_unmasked);

	// Check whether we'll do masking
	do_mask = getMask();
	if (do_mask)
	{
		if (verb > 0)
		{
			std::cout <<"== Masking input maps ..." <<std::endl;
		}
		// Mask I1 and I2 and calculated fsc_masked
		I1() *= Im();
		I2() *= Im();
		getFSC(I1(), I2(), fsc_masked);
		if (do_ampl_corr)
			getAmplitudeCorrelationAndDifferentialPhaseResidual(I1(), I2(), acorr_masked, dpr_masked);

		// To save memory re-read the same input maps again and randomize phases before masking
		I1.read(fn_I1);
		I2.read(fn_I2);
		I1().setXmippOrigin();
		I2().setXmippOrigin();

		// Check at which resolution shell the FSC drops below randomize_fsc_at
		randomize_at = -1;

		if (randomize_at_A > 0.)
		{
			randomize_at = (int)(angpix * XSIZE(I1()) / randomize_at_A );
		}
		else
		{
			// Check when FSC drops below randomize_fsc_at
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_unmasked)
			{
				if (i > 0 && DIRECT_A1D_ELEM(fsc_unmasked, i) < randomize_fsc_at)
				{
					randomize_at = i;
					break;
				}
			}
		}
		if (verb > 0)
		{
			std::cout.width(35); std::cout << std::left << "  + randomize phases beyond: "; std::cout << XSIZE(I1())* angpix / randomize_at << " Angstroms" << std::endl;
		}
		if (randomize_at > 0)
		{
			randomizePhasesBeyond(I1(), randomize_at);
			randomizePhasesBeyond(I2(), randomize_at);
			// Mask randomized phases maps and calculated fsc_random_masked
			I1() *= Im();
			I2() *= Im();
			getFSC(I1(), I2(), fsc_random_masked);
		}
		else
			REPORT_ERROR("Postprocessing::run ERROR: FSC curve never drops below randomize_fsc_at.");

		// Now that we have fsc_masked and fsc_random_masked, calculate fsc_true according to Richard's formula
		// FSC_true = FSC_t - FSC_n / ( )
		calculateFSCtrue(fsc_true, fsc_unmasked, fsc_masked, fsc_random_masked, randomize_at);

		// Also calculate cisTEM-like corrected part_FSC based on expected ordered molecular weight
		calculateFSCpart(fsc_unmasked, frac_molweight, fsc_part_molweight);
		// and based on fraction of white voxels in the solvent mask used for RELION-correction
		calculateFSCpart(fsc_unmasked, frac_solvent_mask, fsc_part_fracmask);

		// Now re-read the original maps yet again into memory
		I1.read(fn_I1);
		I2.read(fn_I2);
		I1().setXmippOrigin();
		I2().setXmippOrigin();

	}
	else
	{
		fsc_true = fsc_unmasked;
	}

	global_resol = 999.;
	// See where corrected FSC drops below 0.143
	int global_resol_i = 0;
	FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_true)
	{
		if ( DIRECT_A1D_ELEM(fsc_true, i) < 0.143)
			break;
		global_resol = XSIZE(I1())*angpix/(RFLOAT)i;
		global_resol_i = i;
	}

	// Determine the fractional resolution (needed for development purposes)
	if (do_interpolate && global_resol_i < fsc_true.xdim - 1)
	{
		const double v0 = DIRECT_A1D_ELEM(fsc_true, global_resol_i) - 0.143;
		const double v1 = DIRECT_A1D_ELEM(fsc_true, global_resol_i + 1) - 0.143;
		const double r = v0 / (v0 - v1);

		fract_resol = XSIZE(I1()) * angpix / (global_resol_i + r);
	}
	else
	{
		fract_resol = global_resol;
	}

	// Perform some checks on phase-randomisation..
	if (do_mask)
	{
		int unmasked_resol_i = 0;
		FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc_unmasked)
		{
			if ( DIRECT_A1D_ELEM(fsc_unmasked, i) < 0.143)
				break;
			unmasked_resol_i = i;
		}
		// Check whether global_resol is worse than the unmasked one
		if (unmasked_resol_i > global_resol_i)
		{
			if (force_mask)
			{
				std::cerr << " WARNING: the unmasked FSC extends beyond the solvent-corrected FSC." << std::endl;
			}
			else
			{
				std::cerr << " WARNING: the unmasked FSC extends beyond the solvent-corrected FSC. Skip masking for now, but you may want to adjust you mask!" << std::endl;
				std::cerr << "          You can force the mask by the '--force_mask' option." << std::endl;
				fsc_true = fsc_unmasked;
				global_resol = XSIZE(I1())*angpix/(RFLOAT)unmasked_resol_i;
			}
		}
		// Check whether the phase-randomised FSC is less than 5% at the resolution estimate, otherwise warn the user
		else if (DIRECT_A1D_ELEM(fsc_random_masked, global_resol_i) > 0.1)
		{
			std::cerr << " WARNING: The phase-randomised FSC is larger than 0.10 at the estimated resolution!" << std::endl;
			std::cerr << " WARNING: This may result in an incorrect resolution estimation. Provide a softer mask with less features to get lower phase-randomised FSCs." << std::endl;
		}
	}

	// Add the two half-maps together for subsequent sharpening
	I1() += I2();
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(I1()) {
		DIRECT_MULTIDIM_ELEM(I1(), n) *= 0.5;
	}

	// Divide by MTF and perform FSC-weighted B-factor sharpening, as in Rosenthal and Henderson, 2003
	// also low-pass filters...
	RFLOAT applied_filter = sharpenMap();

	// Write original and corrected FSC curve, Guinier plot, etc.
	writeOutput();

	if (write_halfmaps) {
		std::cout << "== Writing half maps after applying same sharpening and low-pass filter" << std::endl;
		// Force the filtering as merged map
		do_auto_bfac = false;
		adhoc_bfac = global_bfactor;
		low_pass_freq = applied_filter;
		verb = 0;

		std::cout << "  + half1" << std::endl;
		I1.read(fn_I1);
		I1().setXmippOrigin();
		sharpenMap();
		writeMaps(fn_out + "_half1");

		std::cout << "  + half2" << std::endl;
		I1.read(fn_I2);
		I1().setXmippOrigin();
		sharpenMap();
		writeMaps(fn_out + "_half2");
	}
}

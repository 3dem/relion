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
#include "src/particle_sorter.h"
//#define DEBUG

void ParticleSorter::read(int argc, char **argv)
{

	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "Input STAR file ");
	fn_ref = parser.getOption("--ref", "STAR file with the reference names, or an MRC stack with all references");
	fn_out = parser.getOption("--o", "Output rootname (if empty the input file will be overwritten)", "");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms", "1"));
	angpix_ref = textToFloat(parser.getOption("--angpix_ref", "Pixel size of the references in Angstroms (default is same as micrographs)", "-1"));
	particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circ. mask for the experimental images (in Angstroms, default=automatic)", "-1"));
	lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter in Angstroms for the references (prevent Einstein-from-noise!)","-1"));
	do_invert = parser.checkOption("--invert", "Density in particles is inverted w.r.t. density in template");
	do_ctf = parser.checkOption("--ctf", "Perform CTF correction on the references?");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
	min_z = textToFloat(parser.getOption("--min_z", "Minimum Z-value to count in the sorting of outliers","0"));

	int expert_section = parser.addSection("Expert options");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void ParticleSorter::usage()
{
	parser.writeUsage(std::cout);
}

void ParticleSorter::initialise()
{

	// Read in input metadata file
	MDin.read(fn_in);

	// Check the pixel size
	if (MDin.containsLabel(EMDL_CTF_MAGNIFICATION) && MDin.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
	{
		RFLOAT mag, dstep;
		MDin.goToObject(0);
		MDin.getValue(EMDL_CTF_MAGNIFICATION, mag);
		MDin.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
		angpix = 10000. * dstep / mag;
		if (verb > 0)
			std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
	}
    else if (verb > 0)
    {
    	std::cout << " Warning: input (particle) STAR file does not contain information about pixel size!" << std::endl;
    	std::cout << " Warning: use --angpix to provide the correct value. Now using " << angpix << " Angstroms" << std::endl;
    }

	if (fn_out == "")
		fn_out = fn_in.withoutExtension();
	else
		fn_out = fn_out.withoutExtension();

	// Read in the references
	Mrefs.clear();
	if (fn_ref.isStarFile())
	{

		MetaDataTable MDref;
		if (fn_ref.contains("_model.star"))
		{
			MDref.read(fn_ref, "model_classes");
		}
		else
		{
			// Just read normal STAR file with user-provided references
			MDref.read(fn_ref);
		}

		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDref)
		{
			// Get all reference images and their names
			Image<RFLOAT> Iref;

			FileName fn_img;
			if (!MDref.getValue(EMDL_MLMODEL_REF_IMAGE, fn_img))
			{
				if (!MDref.getValue(EMDL_IMAGE_NAME, fn_img)) // perhaps the reference are called rlnImageName?
					REPORT_ERROR("ParticleSorter::initialise ERROR: either provide rlnReferenceImage or rlnImageName in the reference STAR file!");
			}
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

	// Automated determination of bg_radius (same code as in autopicker.cpp!)
	if (particle_diameter < 0.)
	{
		RFLOAT sumr=0.;
		for (int iref = 0; iref < Mrefs.size(); iref++)
		{
			RFLOAT cornerval = DIRECT_MULTIDIM_ELEM(Mrefs[iref], 0);
			// Look on the central X-axis, which first and last values are NOT equal to the corner value
			int last_corner=99999, first_corner=99999;
			for (long int j=STARTINGX(Mrefs[iref]); j<=FINISHINGX(Mrefs[iref]); j++)
			{
				if (first_corner == 99999)
				{
					if (A3D_ELEM(Mrefs[iref], 0,0,j) != cornerval)
						first_corner = j;
				}
				else if (last_corner  == 99999)
				{
					if (A3D_ELEM(Mrefs[iref], 0,0,j) == cornerval)
						last_corner = j - 1;
				}
			}
			sumr += (last_corner - first_corner);
		}
		sumr /= 2. * Mrefs.size(); // factor 2 to go from diameter to radius; Mref.size() for averaging
		particle_radius2 = sumr*sumr;
		std::cout << " Automatically set the background radius to " << sumr << " pixels in the references" << std::endl;
		std::cout << " You can override this by providing --particle_diameter (in Angstroms)" << std::endl;
	}
	else
	{
		// Get the squared particle radius (in integer pixels)
		particle_radius2 = ROUND(particle_diameter/(2. * angpix));
		particle_radius2*= particle_radius2;
	}


	// Re-scale references if necessary
	if (angpix_ref < 0)
		angpix_ref = angpix;
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

		// Also adjust particle_radius2
		particle_radius2 = ROUND(particle_radius2 * (angpix_ref/angpix) * angpix_ref/angpix);

	}

	particle_size = XSIZE(Mrefs[0]);



	// Invert references if necessary (do this AFTER automasking them!)
	if (do_invert)
	{
		for (int iref = 0; iref < Mrefs.size(); iref++)
		{
			Mrefs[iref] *= -1.;
		}
	}


	// Calculate the size of the FTs given the low-pass filter
	if (lowpass < 0.)
	{
		lowpass = 2. * angpix;
		current_size = particle_size;
	}
	else
	{
		current_size = 2 * ROUND(particle_size * angpix / lowpass);
	}

	// Calculate (downsized) Fourier transforms of the references
	PPref.clear();
	MultidimArray<RFLOAT> dummy;
	Projector projector(particle_size);
	for (int iref = 0; iref < Mrefs.size(); iref++)
	{
		projector.computeFourierTransformMap(Mrefs[iref], dummy, current_size);
		PPref.push_back(projector);
	}

#ifdef DEBUG
	std::cerr << "Finishing initialise" << std::endl;
#endif
}

void ParticleSorter::run()
{


	int nr_parts = MDin.numberOfObjects();
	features.resize(nr_parts, NR_FEATURES);

	int barstep;
	if (verb > 0)
	{
		std::cout << "Calculating sorting features for all input particles..." << std::endl;
		init_progress_bar(nr_parts);
		barstep = XMIPP_MAX(1, nr_parts / 60);
	}

	long int ipart = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
	{

		if (verb > 0 && ipart % barstep == 0)
			progress_bar(ipart);

		calculateFeaturesOneParticle(ipart);

	    ipart++;
	}

	if (verb > 0)
		progress_bar(nr_parts);

	normaliseFeatures();

	writeFeatures();


}

void ParticleSorter::calculateFeaturesOneParticle(long int ipart)
{

	FourierTransformer transformer;
	Image<RFLOAT> img;
	MultidimArray<RFLOAT> Mref_rot;
	MultidimArray<Complex > Fimg, Fref;

	// Get the image
	FileName fn_img;
	MDin.getValue(EMDL_IMAGE_NAME, fn_img, ipart);
	img.read(fn_img);

	if (XSIZE(img()) != particle_size)
	{
		std::cerr << " fn_img= " << fn_img << " XSIZE(img())= " << XSIZE(img()) << " reference size= " << particle_size << std::endl;
		REPORT_ERROR("ParticleSorter::calculateFeaturesOneParticle ERROR: References and images do not have the same size!");
	}
#ifdef DEBUG
	std::cerr << " fn_img= " << fn_img << std::endl;
	Image<RFLOAT> tt;
	tt = img;
	tt.write("debug_ori_image.spi");
#endif

	// If the particle has non-zero shifts, then apply the rounded shifts here!
	Matrix1D<RFLOAT> offset(2);
	MDin.getValue(EMDL_ORIENT_ORIGIN_X, XX(offset), ipart);
	MDin.getValue(EMDL_ORIENT_ORIGIN_Y, YY(offset), ipart);
	offset.selfROUND();
	if ( ABS(XX(offset)) > 0. || ABS(YY(offset)) > 0.)
		selfTranslate(img(), offset, DONT_WRAP);

#ifdef DEBUG
	tt = img;
	tt.write("debug_trans_image.spi");
#endif

	// Low-pass filter the image if necessary
	if (lowpass > 0.)
	{
		RFLOAT radius = XSIZE(img()) * angpix / lowpass;
		radius -= WIDTH_FMASK_EDGEB / 2.;
		RFLOAT radius_p = radius + WIDTH_FMASK_EDGEB;
		transformer.FourierTransform(img(), Fimg);
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fimg)
		{
			RFLOAT r = sqrt((RFLOAT)(kp*kp + ip*ip + jp*jp));
			if (r < radius)
				continue;
			else if (r > radius_p)
				DIRECT_A3D_ELEM(Fimg, k, i, j) = 0.;
			else
			{
				DIRECT_A3D_ELEM(Fimg, k, i, j) *= 0.5 - 0.5 * cos(PI * (radius_p - r) / WIDTH_FMASK_EDGEB);
			}
		}
		transformer.inverseFourierTransform(Fimg, img());
	}

#ifdef DEBUG
	tt = img;
	tt.write("debug_lowpass_image.spi");
#endif

	// Get the transformation parameters and which reference this is
	RFLOAT psi = 0.;
	RFLOAT rot = 0.;
	RFLOAT tilt = 0.;
	int iref;
	MDin.getValue(EMDL_ORIENT_ROT, rot, ipart);
	MDin.getValue(EMDL_ORIENT_TILT, tilt, ipart);
	MDin.getValue(EMDL_ORIENT_PSI, psi, ipart);
	if (!MDin.getValue(EMDL_PARTICLE_CLASS, iref, ipart) && PPref.size() > 1)
		REPORT_ERROR("ParticleSorter::calculateFeatures: cannot find class number in input STAR file");

	// Manually added particles were given a class number of -999 in the relion_display program
	// These particles should be OK, and therefore will be at the top of the sorted list
	if (iref == -999)
	{
		DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_AVG) = 0.;
		DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_SIG) = 0.;
		DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_SKW) = 0.;
		DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_KRT) = 0.;
		DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_QUADSIG) = 0.;
		DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_ROTFOURCORR) = 0.;

		return;
	}

	iref -= 1; // start counting at 0 instead of 1
	if (iref >= PPref.size())
	{
		MetaDataTable MDt;
		MDt.addObject(MDin.getObject(ipart));
		MDt.write(std::cerr);
		REPORT_ERROR("Too large class number for the given number of references: " + integerToString(iref));
	}

	// Get the rotational matrix
	Matrix2D<RFLOAT> A;
	Euler_angles2matrix(rot, tilt, psi, A);

	// Get the reference image in the right orientation
	Fref.resize(current_size, current_size/2 + 1);
	PPref[iref].get2DFourierTransform(Fref, A, IS_NOT_INV);
	if (do_ctf)
	{
		CTF ctf;
		MultidimArray<RFLOAT> Fctf;
		ctf.read(MDin, MDin, ipart);
		Fctf.resize(Fref);
		ctf.getFftwImage(Fctf, XSIZE(img()), YSIZE(img()), angpix, false, false, intact_ctf_first_peak, true);
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fref)
		{
			DIRECT_MULTIDIM_ELEM(Fref, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
		}
	}
	if (lowpass > 0.)
	{
		// Set the Fourier transform back to its original size...
		MultidimArray<Complex > Frefp;
		windowFourierTransform(Fref, Frefp, XSIZE(img()));
		Fref = Frefp;
	}
	Mref_rot.resize(img());
	transformer.inverseFourierTransform(Fref, Mref_rot);
	CenterFFT(Mref_rot, true);

#ifdef DEBUG
	tt() = Mref_rot;
	tt.write("debug_rotated_ctf_ref.spi");
#endif

	/*
	// Calculate the optimal scale between the image and the reference:
	RFLOAT sumxa = 0.;
	RFLOAT suma2 = 0.;
	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mref_rot)
	{
		if (DIRECT_MULTIDIM_ELEM(iMsk, n) > 0.5)
		{
			sumxa += DIRECT_MULTIDIM_ELEM(Mref_rot, n) * DIRECT_MULTIDIM_ELEM(img(), n);
			suma2 += DIRECT_MULTIDIM_ELEM(Mref_rot, n) * DIRECT_MULTIDIM_ELEM(Mref_rot, n);
		}
	}
	if (suma2 < 1e-10)
		REPORT_ERROR("ERROR: Empty reference encountered!!");
	RFLOAT scale = sumxa/suma2;
	// Calculate the difference between the particle and the oriented reference
	img() -= scale * Mref_rot;
	*/

	// Calculate the difference between the particle and the oriented reference
	img() -= Mref_rot;

	// Let's mask it already here (also useful for RotationalSymmetryInFourierTransform)
	img().setXmippOrigin();
	FOR_ALL_ELEMENTS_IN_ARRAY2D(img())
	{
		if (i*i + j*j >= particle_radius2)
			A2D_ELEM(img(), i, j) = 0.;
	}

#ifdef DEBUG
	tt = img;
	tt.write("debug_differene_image.spi");
#endif

	// Calculate real-space statistics on the difference image (only within the masked protein area)
	//iMsk.initConstant(1.);
    RFLOAT mean = 0., stddev = 0., skew = 0., kurt = 0., quadrant_stddev = 0., rotsym = 0.;
	calculateStatsOneImage(img(), mean, stddev, skew, kurt, quadrant_stddev);
	DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_AVG) = mean;
	DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_SIG) = stddev;
	DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_SKW) = skew;
	DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_KRT) = kurt;
	DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_QUADSIG) = quadrant_stddev;

	// Calculate the rotational Symmetry of the Fourier Transform of the masked difference image
	transformer.FourierTransform(img(), Fimg);
	rotsym = rotationalSymmetryFourierTransform(Fimg);
	DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_ROTFOURCORR) = rotsym;


#ifdef DEBUG
	std::cerr << "Difference image protein area" << std::endl;
	std::cerr << " mean= " << mean << " stddev= " << stddev << " skew= " << skew << " kurt= " << kurt << " quadrant_stddev= " <<  quadrant_stddev << " rotSym= " << rotsym << std::endl;
	std::cerr << " Press any key to continue... "  << std::endl;
	char c;
	std::cin >> c;
#endif

}

void ParticleSorter::normaliseFeatures()
{

	// Ignore manually added particles, which will have all features set to 0...
	std::vector< std::vector<RFLOAT> > normfeatures(XSIZE(features));
	MultidimArray<RFLOAT> onerow;
	for (int ipart = 0; ipart < YSIZE(features); ipart++)
	{
		features.getRow(ipart, onerow);
		RFLOAT abssum = 0.;
		for (int ifeat = 0; ifeat < XSIZE(features); ifeat++)
		{
			abssum += ABS(DIRECT_A1D_ELEM(onerow, ifeat));
		}
		if (abssum > 0.)
		{
			for (int ifeat = 0; ifeat < XSIZE(features); ifeat++)
			{
				normfeatures[ifeat].push_back(DIRECT_A1D_ELEM(onerow, ifeat));
			}
		}
	}

	// Now calculate average and stddev for all particles which did not have all features set to 0
	std::vector<RFLOAT> avg(XSIZE(features), 0.), stddev(XSIZE(features), 0.);
	for (int ifeat = 0; ifeat < XSIZE(features); ifeat++)
	{
		// average
		for (long int i =0; i < (normfeatures[ifeat]).size(); i++)
			avg[ifeat] += normfeatures[ifeat][i];
		avg[ifeat] /= (normfeatures[ifeat]).size();

		// stddev
		for (long int i =0; i < (normfeatures[ifeat]).size(); i++)
			stddev[ifeat] += (normfeatures[ifeat][i] - avg[ifeat]) * (normfeatures[ifeat][i] - avg[ifeat]);
		stddev[ifeat] = sqrt(stddev[ifeat] / (normfeatures[ifeat]).size());
	}

	// Calculate Z-scores for all particles, except the ones that were added manually
	for (int ipart = 0; ipart < YSIZE(features); ipart++)
	{
		features.getRow(ipart, onerow);
		RFLOAT abssum = 0.;
		for (int ifeat = 0; ifeat < XSIZE(features); ifeat++)
		{
			abssum += ABS(DIRECT_A1D_ELEM(onerow, ifeat));
		}
		if (abssum > 0.)
		{
			for (int ifeat = 0; ifeat < XSIZE(features); ifeat++)
			{
				DIRECT_A1D_ELEM(onerow, ifeat) =  (DIRECT_A1D_ELEM(onerow, ifeat) - avg[ifeat]) / stddev[ifeat];
				if (DIRECT_A1D_ELEM(onerow, ifeat) < min_z)
					DIRECT_A1D_ELEM(onerow, ifeat) = 0.;
			}
		}
		else
		{
			// Set negative Z-scores (-1) for manually picked particles
			for (int ifeat = 0; ifeat < XSIZE(features); ifeat++)
			{
				DIRECT_A1D_ELEM(onerow, ifeat) = - 1.;
			}
		}
		features.setRow(ipart, onerow);
	}

}

void ParticleSorter::writeFeatures()
{

	FileName fn_img;
	FileName fn_tmp = fn_out +"_features.libsvm";
#ifdef WRITE_LIBSVM
	std::ofstream  fh, fh2, fh3;
	fh.open((fn_tmp).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)" ParticleSorter::writeFeatures: Cannot write file: " + fn_tmp);
	fn_tmp = fn_out +"_features.dat";
	fh2.open((fn_tmp).c_str(), std::ios::out);
	if (!fh2)
		REPORT_ERROR( (std::string)" ParticleSorter::writeFeatures: Cannot write file: " + fn_tmp);
#endif

	MultidimArray<RFLOAT> sumfeatures(YSIZE(features));
	sumfeatures.initZeros();
	for (int ipart = 0; ipart < YSIZE(features); ipart++)
	{
#ifdef WRITE_LIBSVM
		fh2 << ipart << " ";
#endif
		for (int ifeat = 0; ifeat < XSIZE(features); ifeat++)
		{
#ifdef WRITE_LIBSVM
			fh << ifeat <<":"<<DIRECT_A2D_ELEM(features, ipart, ifeat)<<" ";
			fh2 << DIRECT_A2D_ELEM(features, ipart, ifeat)<<" ";
#endif
			DIRECT_A1D_ELEM(sumfeatures, ipart) += DIRECT_A2D_ELEM(features, ipart, ifeat);
		}
		DIRECT_A1D_ELEM(sumfeatures, ipart) /= XSIZE(features);
#ifdef WRITE_LIBSVM
		MultidimArray<RFLOAT> row(1);
		features.getRow(ipart, row);
		fh2<< " sum= " << DIRECT_A1D_ELEM(sumfeatures, ipart);
		MDin.getValue(EMDL_IMAGE_NAME, fn_img, ipart);
		fh2<< " "<< fn_img;
		fh << std::endl;
		fh2 << std::endl;
#endif
	}
#ifdef WRITE_LIBSVM
	fh.close();
	fh2.close();
#endif

	int ipart = 0;
	FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDin)
	{
#ifdef DEBUG
		MDin.setValue(EMDL_IMAGE_STATS_AVG, DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_AVG));
		MDin.setValue(EMDL_IMAGE_STATS_STDDEV, DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_SIG));
		MDin.setValue(EMDL_IMAGE_STATS_SKEW, DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_SKW));
		MDin.setValue(EMDL_IMAGE_STATS_KURT, DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_KRT));
		MDin.setValue(EMDL_PARTICLE_DLL, DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_QUADSIG));
		MDin.setValue(EMDL_POSTPROCESS_GUINIER_FIT_CORRELATION, DIRECT_A2D_ELEM(features, ipart, FEATURE_DF_ROTFOURCORR));
#endif
		MDin.setValue(EMDL_SELECT_PARTICLES_ZSCORE, DIRECT_A1D_ELEM(sumfeatures, ipart));
		ipart++;
	}
	fn_tmp = fn_out + ".star";
	MDin.write(fn_tmp);


	if (verb>0)
		std::cout <<" Written out "<< fn_tmp << std::endl;
}

RFLOAT ParticleSorter::rotationalSymmetryFourierTransform(MultidimArray<Complex > &Fimg)
{
	RFLOAT result;
	// Calculate average and stddev per ring

	int radius = ROUND(YSIZE(Fimg) * angpix / lowpass);
	MultidimArray<RFLOAT> count, power;

    power.initZeros(radius+1);
    count.initZeros(radius+1);
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Fimg)
    {
    	long int idx = ROUND(sqrt(ip*ip + jp*jp));
    	if (idx < radius)
    	{
			power(idx) += norm(DIRECT_A2D_ELEM(Fimg, i, j));
			count(idx) += 1.;
    	}
    }

    for (long int i = 0; i < radius+1; i++)
    {        if (count(i) > 0.)
            power(i) /= count(i);
    }

    RFLOAT sum = 0.;
    RFLOAT sum2 = 0.;
    RFLOAT cc  = 0.;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Fimg)
    {
    	long int idx = ROUND(sqrt(ip*ip + jp*jp));
    	if (idx < radius)
    	{
    		RFLOAT aux = norm(DIRECT_A2D_ELEM(Fimg, i, j)) - power(idx);
    		sum += aux;
    		sum2 += aux*aux;
    		cc++;
    	}
    }

    sum  /= cc;
	sum2 /= cc;
	result = sum2 - sum * sum;
	result *= YXSIZE(Fimg)*YXSIZE(Fimg);
#ifdef DEBUG
	std::cerr <<" rotSym="<<result<<std::endl;
#endif
	return result;
}

void ParticleSorter::calculateStatsOneImage(MultidimArray<RFLOAT> &img,
							RFLOAT &mean, RFLOAT &stddev, RFLOAT &skew, RFLOAT &kurt, RFLOAT &quadrant_stddev)
{

	// First pass: calculate mean
	int na = 0;
	int n1 = 0;
	int n2 = 0;
	int n3 = 0;
	int n4 = 0;
	RFLOAT mean_q1 = 0.;
	RFLOAT mean_q2 = 0.;
	RFLOAT mean_q3 = 0.;
	RFLOAT mean_q4 = 0.;
	RFLOAT stddev_q1 = 0.;
	RFLOAT stddev_q2 = 0.;
	RFLOAT stddev_q3 = 0.;
	RFLOAT stddev_q4 = 0.;
	mean = 0.;
	stddev = 0;
	skew   = 0.;
	kurt   = 0.;

	img.setXmippOrigin();

	FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
	{
		if (i*i + j*j < particle_radius2)
		{
			mean += A2D_ELEM(img, i, j);
			na++;

			if (i > 0 && j > 0)
			{
				mean_q1 += A2D_ELEM(img, i, j);
				n1++;
			}
			else if (i > 0 && j <= 0)
			{
				mean_q2 += A2D_ELEM(img, i, j);
				n2++;
			}
			else if (i <= 0 && j > 0)
			{
				mean_q3 += A2D_ELEM(img, i, j);
				n3++;
			}
			else if (i <= 0 && j <= 0)
			{
				mean_q4 += A2D_ELEM(img, i, j);
				n4++;
			}
		}
	}

	// Second pass calculate higher moments
	mean /= na;
	mean_q1 /= n1;
	mean_q2 /= n2;
	mean_q3 /= n3;
	mean_q4 /= n4;

	FOR_ALL_ELEMENTS_IN_ARRAY2D(img)
	{
		if (i*i + j*j < particle_radius2)
		{
			RFLOAT aux = A2D_ELEM(img, i, j) - mean;
			RFLOAT aux2 = aux * aux;
			stddev += aux2;
			skew += aux2 * aux;
			kurt += aux2 * aux2;
			if (i > 0 && j > 0)
			{
				aux = A2D_ELEM(img, i, j) - mean_q1;
				aux2 = aux * aux;
				stddev_q1 += aux2;
			}
			else if (i > 0 && j <= 0)
			{
				aux = A2D_ELEM(img, i, j) - mean_q2;
				aux2 = aux * aux;
				stddev_q2 += aux2;
			}
			else if (i <= 0 && j > 0)
			{
				aux = A2D_ELEM(img, i, j) - mean_q3;
				aux2 = aux * aux;
				stddev_q3 += aux2;
			}
			else if (i <= 0 && j <= 0)
			{
				aux = A2D_ELEM(img, i, j) - mean_q4;
				aux2 = aux * aux;
				stddev_q4 += aux2;
			}
		}
	}
	RFLOAT var = stddev/na;
	stddev = sqrt(var);
	skew = (skew/na)/(var*stddev);
	kurt = (kurt/na)/(var*var) - 3.;

	var = stddev_q1/n1;
	stddev_q1 = sqrt(var);

	var = stddev_q2/n2;
	stddev_q2 = sqrt(var);

	var = stddev_q3/n3;
	stddev_q3 = sqrt(var);

	var = stddev_q4/n4;
	stddev_q4 = sqrt(var);

	RFLOAT mean_stddevs = 0.;
	RFLOAT stddev_stddevs = 0.;
	mean_stddevs = (stddev_q1 + stddev_q2 + stddev_q3 + stddev_q4) / 4.;
	stddev_stddevs += (stddev_q1-mean_stddevs)*(stddev_q1-mean_stddevs) + (stddev_q2-mean_stddevs)*(stddev_q2-mean_stddevs) + (stddev_q3-mean_stddevs)*(stddev_q3-mean_stddevs) + (stddev_q4-mean_stddevs)*(stddev_q4-mean_stddevs);
	quadrant_stddev = sqrt(stddev_stddevs / 4.);

}



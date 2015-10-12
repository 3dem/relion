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

void AutoPicker::read(int argc, char **argv)
{

	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_in = parser.getOption("--i", "Micrograph STAR file OR filenames from which to autopick particles, e.g. \"Micrographs/*.mrc\"");
	fn_out = parser.getOption("--o", "Output rootname", "autopick");
	angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms"));
	particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms)"));
	do_write_fom_maps = parser.checkOption("--write_fom_maps", "Write calculated probability-ratio maps to disc (for re-reading in subsequent runs)");
	do_read_fom_maps = parser.checkOption("--read_fom_maps", "Skip probability calculations, re-read precalculated maps from disc");

	int ref_section = parser.addSection("References options");
	fn_ref = parser.getOption("--ref", "STAR file with the reference names, or an MRC stack with all references");
	do_invert = parser.checkOption("--invert", "Density in micrograph is inverted w.r.t. density in template");
	psi_sampling = textToFloat(parser.getOption("--ang", "Angular sampling (in degrees); use 360 for no rotations", "10"));
	lowpass = textToFloat(parser.getOption("--lowpass", "Lowpass filter in Angstroms for the references (prevent Einstein-from-noise!)","-1"));
	do_ctf = parser.checkOption("--ctf", "Perform CTF correction on the references?");
	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");

	int peak_section = parser.addSection("Peak-search options");
	min_fraction_expected_Pratio = textToFloat(parser.getOption("--threshold", "Fraction of expected probability ratio in order to consider peaks?", "0.25"));
	min_particle_distance = textToFloat(parser.getOption("--min_distance", "Minimum distance (in A) between any two particles (default is half the box size)","-1"));
	autopick_skip_side = textToInteger(parser.getOption("--skip_side", "Keep this many extra pixels (apart from particle_size/2) away from the edge of the micrograph ","0"));

	int expert_section = parser.addSection("Expert options");
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

}

void AutoPicker::usage()
{
	parser.writeUsage(std::cerr);
}


void AutoPicker::initialise()
{

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
	}
	else
	{

		if (do_ctf)
			REPORT_ERROR("AutoPicker::initialise ERROR: use an input STAR file with the CTF information when using --ctf");

		fn_in.globFiles(fn_micrographs);
		if (fn_micrographs.size() == 0)
			REPORT_ERROR("Cannot find any micrograph called: "+fns_autopick);
	}


	if (verb > 0)
	{
		std::cout << " Run autopicking on the following micrographs: " << std::endl;
		for(unsigned  int  i = 0; i < fn_micrographs.size(); ++i)
			std::cout << "  * " << fn_micrographs[i] << std::endl;
	}

	// Read in the references
	Mrefs.clear();
	if (fn_ref.isStarFile())
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
			Iref().printShape();
			Mrefs.push_back(Iref());
		}
	}

	particle_size = XSIZE(Mrefs[0]);

	// Get the squared particle radius (in integer pixels)
	particle_radius2 = ROUND(particle_diameter/(2. * angpix));
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

	if (min_particle_distance < 0)
	{
		min_particle_distance = particle_size * angpix / 2.;
	}

	// Pre-calculate and store Projectors for all references at the right size
	if (!do_read_fom_maps)
	{
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

		// Now set the mask in the large square and store its FFT
		Maux.initZeros();
		FOR_ALL_ELEMENTS_IN_ARRAY2D(Mcirc_mask)
		{
			A2D_ELEM(Maux, i, j ) = A2D_ELEM(Mcirc_mask, i, j);
		}
		CenterFFT(Maux, true);
		transformer.FourierTransform(Maux, Fmsk);
		transformer.cleanup(); // somehow I get segfaults without this cleanup....

		PPref.clear();
		Projector PP(micrograph_size);
		MultidimArray<RFLOAT> dummy;

		// TODO!!! sum_ref etc should be CTF-dependent!!!
		//sum_ref_under_circ_mask.clear();
		//sum_ref2_under_circ_mask.clear();
		for (int iref = 0; iref < Mrefs.size(); iref++)
		{

			// (Re-)apply the mask to the references
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Mrefs[iref])
			{
				DIRECT_MULTIDIM_ELEM(Mrefs[iref], n) *= DIRECT_MULTIDIM_ELEM(Mcirc_mask, n);
			}

			// Set reference in the large box of the micrograph
			Maux.initZeros();
			FOR_ALL_ELEMENTS_IN_ARRAY2D(Mrefs[iref])
			{
				A2D_ELEM(Maux, i, j) = A2D_ELEM(Mrefs[iref], i, j);
			}

			// And compute its Fourier Transform inside the Projector
			PP.computeFourierTransformMap(Maux, dummy, downsize_mic, 1, false);
			PPref.push_back(PP);

		}
	}

#ifdef DEBUG
	std::cerr << "Finishing initialise" << std::endl;
#endif

}

void AutoPicker::run()
{

	int barstep;
	if (verb > 0)
	{
		std::cout << " Autopicking ..." << std::endl;
		init_progress_bar(fn_micrographs.size());
		barstep = XMIPP_MAX(1, fn_micrographs.size() / 60);
	}


	for (long int imic = 0; imic < fn_micrographs.size(); imic++)
	{

		if (verb > 0 && imic % barstep == 0)
			progress_bar(imic);

		autoPickOneMicrograph(fn_micrographs[imic]);

	}

	if (verb > 0)
		progress_bar(fn_micrographs.size());

}

void AutoPicker::autoPickOneMicrograph(FileName &fn_mic)
{
	Image<RFLOAT> Imic;
	MultidimArray<Complex > Faux, Faux2, Fmic;
	MultidimArray<RFLOAT> Maux, Mstddev, Mmean, Mdiff2, MsumX2, Mccf_best, Mpsi_best, Fctf;
	FourierTransformer transformer;
	RFLOAT sum_ref_under_circ_mask, sum_ref2_under_circ_mask;
	int my_skip_side = autopick_skip_side + particle_size/2;
	CTF ctf;

	int min_distance_pix = ROUND(min_particle_distance / angpix);

#ifdef DEBUG
	Image<RFLOAT> tt;
	tt().resize(micrograph_size, micrograph_size);
	std::cerr << " fn_mic= " << fn_mic << std::endl;
#endif
	// Read in the micrograph
	Imic.read(fn_mic);
	Imic().setXmippOrigin();

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

	// Set mean to zero and stddev to 1 to prevent numerical problems with one-sweep stddev calculations....
	Imic().statisticsAdjust(0., 1.);

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

	Mccf_best.resize(micrograph_size, micrograph_size);
	Mpsi_best.resize(micrograph_size, micrograph_size);

	RFLOAT normfft = (RFLOAT)(micrograph_size * micrograph_size) / (RFLOAT)nr_pixels_circular_mask;;
	if (!do_read_fom_maps)
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

		// Also calculate the FFT of the squared micrograph
		Maux.resize(Imic());
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Maux)
		{
			DIRECT_MULTIDIM_ELEM(Maux, n) = DIRECT_MULTIDIM_ELEM(Imic(), n) * DIRECT_MULTIDIM_ELEM(Imic(), n);
		}
		MultidimArray<Complex > Fmic2;
		transformer.FourierTransform(Maux, Fmic2);

#ifdef DEBUG
		std::cerr << " nr_pixels_circular_invmask= " << nr_pixels_circular_invmask << std::endl;
		std::cerr << " nr_pixels_circular_mask= " << nr_pixels_circular_mask << std::endl;
		windowFourierTransform(Finvmsk, Faux2, micrograph_size);
		transformer.inverseFourierTransform(Faux2, tt());
		CenterFFT(tt(), false);
		tt.write("Minvmask.spi");
		windowFourierTransform(Fmsk, Faux2, micrograph_size);
		transformer.inverseFourierTransform(Faux2, tt());
		CenterFFT(tt(), false);
		tt.write("Mmask.spi");
#endif

		// The following calculate mu and sig under the solvent area at every position in the micrograph
		calculateStddevAndMeanUnderMask(Fmic, Fmic2, Finvmsk, nr_pixels_circular_invmask, Mstddev, Mmean);

		// From now on use downsized Fmic, as the cross-correlation with the references can be done at lower resolution
		windowFourierTransform(Fmic, Faux, downsize_mic);
		Fmic = Faux;

	}// end if do_read_fom_maps

	// Now start looking for the peaks of all references
	// Clear the output vector with all peaks
	std::vector<Peak> peaks;
	peaks.clear();
	for (int iref = 0; iref < Mrefs.size(); iref++)
	{
		RFLOAT expected_Pratio; // the expectedFOM for this (ctf-corrected) reference
		if (do_read_fom_maps)
		{
			FileName fn_tmp;
			Image<RFLOAT> It;
			fn_tmp.compose(fn_mic.withoutExtension()+"_"+fn_out+"_ref", iref,"_bestCCF.spi");
			It.read(fn_tmp);
			Mccf_best = It();
			// Retrieve expected_Pratio from the header of the image..
			It.MDMainHeader.getValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);
			fn_tmp.compose(fn_mic.withoutExtension()+"_"+fn_out+"_ref", iref,"_bestPSI.spi");
			It.read(fn_tmp);
			Mpsi_best = It();

		} //end else if do_read_fom_maps
		else
		{
			Mccf_best.initConstant(-99.e99);
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
				transformer.inverseFourierTransform(Faux2, tt());
				CenterFFT(tt(), false);
				tt.write("Mref_rot.spi");

				windowFourierTransform(Fmic, Faux2, micrograph_size);
				transformer.inverseFourierTransform(Faux2, tt());
				CenterFFT(tt(), false);
				tt.write("Mmic.spi");

#endif

				// Apply the CTF on-the-fly (so same PPref can be used for many different micrographs)
				if (do_ctf)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
					{
						DIRECT_MULTIDIM_ELEM(Faux, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
#ifdef DEBUG
				windowFourierTransform(Faux, Faux2, micrograph_size);
				transformer.inverseFourierTransform(Faux2, Maux);
				CenterFFT(Maux, false);
				Maux.setXmippOrigin();
				tt().resize(particle_size, particle_size);
				tt().setXmippOrigin();
				FOR_ALL_ELEMENTS_IN_ARRAY2D(tt())
				{
					A2D_ELEM(tt(), i, j) = A2D_ELEM(Maux, i, j);
				}
				tt.write("Mref_rot_ctf.spi");
#endif
				}

				if (is_first_psi)
				{
					// Calculate the expected ratio of probabilities for this CTF-corrected reference
					// and the sum_ref_under_circ_mask and sum_ref_under_circ_mask2
					// Do this also if we're not recalculating the fom maps...

					windowFourierTransform(Faux, Faux2, micrograph_size);
					transformer.inverseFourierTransform(Faux2, Maux);
					CenterFFT(Maux, false);
					Maux.setXmippOrigin();
					// TODO: check whether I need CenterFFT(Maux, false)

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
#endif
				}

				// Now multiply template and micrograph to calculate the cross-correlation
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Faux)
				{
					DIRECT_MULTIDIM_ELEM(Faux, n) = conj(DIRECT_MULTIDIM_ELEM(Faux, n)) * DIRECT_MULTIDIM_ELEM(Fmic, n);
				}
				windowFourierTransform(Faux, Faux2, micrograph_size);
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
#ifdef DEBUG
						/*
						if (diff2 < 0. || n==28800 || n==0)
						{
							std::cerr << " n= "<<n<< "diff2= " << diff2 << " old Mdiff2=" <<DIRECT_MULTIDIM_ELEM(Mdiff2, n)
									<< " -2AX/sig " << - 2. * normfft * DIRECT_MULTIDIM_ELEM(Maux, n) / DIRECT_MULTIDIM_ELEM(Mstddev, n)
									<< " 2Amu/sig= " << 2. * DIRECT_MULTIDIM_ELEM(Mmean, n) * sum_ref_under_circ_mask[iref] / DIRECT_MULTIDIM_ELEM(Mstddev, n)
									<< " A2=" <<  sum_ref2_under_circ_mask[iref]
									<< " stddev= " <<  DIRECT_MULTIDIM_ELEM(Mstddev, n) << " avg= "<< DIRECT_MULTIDIM_ELEM(Mmean, n)
									<< std::endl;
						}
						*/
#endif
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
			} // end for psi


			if (do_write_fom_maps)
			{
				// TMP output
				FileName fn_tmp;
				Image<RFLOAT> It;
				It() = Mccf_best;
				// Store expected_Pratio in the header of the image..
				It.MDMainHeader.setValue(EMDL_IMAGE_STATS_MAX, expected_Pratio);;
				fn_tmp.compose(fn_mic.withoutExtension()+"_"+fn_out+"_ref", iref,"_bestCCF.spi");
				It.write(fn_tmp);

				It() = Mpsi_best;
				fn_tmp.compose(fn_mic.withoutExtension()+"_"+fn_out+"_ref", iref,"_bestPSI.spi");
				It.write(fn_tmp);
			} // end if do_write_fom_maps

		} // end if do_read_fom_maps

		// Now that we have Mccf_best and Mpsi_best, get the peaks
		std::vector<Peak> my_ref_peaks;
		Mccf_best.setXmippOrigin();
		Mpsi_best.setXmippOrigin();
		peakSearch(Mccf_best, Mpsi_best, iref, my_skip_side, my_ref_peaks);

		prunePeakClusters(my_ref_peaks, min_distance_pix);

		// append the peaks of this reference to all the other peaks
		peaks.insert(peaks.end(), my_ref_peaks.begin(), my_ref_peaks.end());

	} // end for iref


	//Now that we have done all references, prune the list again...
	prunePeakClusters(peaks, min_distance_pix);

	// And remove all too close neighbours
	removeTooCloselyNeighbouringPeaks(peaks, min_distance_pix);

	// Write out a STAR file with the coordinates
	MetaDataTable MDout;
	for (int ipeak =0; ipeak < peaks.size(); ipeak++)
	{
		MDout.addObject();
		MDout.setValue(EMDL_IMAGE_COORD_X, (RFLOAT)(peaks[ipeak].x));
		MDout.setValue(EMDL_IMAGE_COORD_Y, (RFLOAT)(peaks[ipeak].y));
		MDout.setValue(EMDL_ORIENT_PSI, peaks[ipeak].psi);
		MDout.setValue(EMDL_PARTICLE_CLASS, peaks[ipeak].ref + 1); // start counting at 1
		MDout.setValue(EMDL_PARTICLE_AUTOPICK_FOM, peaks[ipeak].fom);
	}
	FileName fn_tmp = fn_mic.withoutExtension() + "_" + fn_out + ".star";
	MDout.write(fn_tmp);

}

void AutoPicker::calculateStddevAndMeanUnderMask(const MultidimArray<Complex > &_Fmic, const MultidimArray<Complex > &_Fmic2,
		MultidimArray<Complex > &_Fmsk, int nr_nonzero_pixels_mask, MultidimArray<RFLOAT> &_Mstddev, MultidimArray<RFLOAT> &_Mmean)
{

	MultidimArray<Complex > Faux, Faux2;
	MultidimArray<RFLOAT> Maux(micrograph_size, micrograph_size);
	FourierTransformer transformer;

	_Mstddev.initZeros(micrograph_size, micrograph_size);
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
	windowFourierTransform(Faux, Faux2, micrograph_size);
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
	windowFourierTransform(Faux, Faux2, micrograph_size);
	transformer.inverseFourierTransform(Faux2, Maux);

	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(_Mstddev)
	{
		// we already stored minus average-squared in _Mstddev
		DIRECT_MULTIDIM_ELEM(_Mstddev, n) += normfft * DIRECT_MULTIDIM_ELEM(Maux, n);
		if (DIRECT_MULTIDIM_ELEM(_Mstddev, n) > 0.)
			DIRECT_MULTIDIM_ELEM(_Mstddev, n) = sqrt(DIRECT_MULTIDIM_ELEM(_Mstddev, n) );
		else
			DIRECT_MULTIDIM_ELEM(_Mstddev, n) = 0.;
	}

	CenterFFT(_Mstddev, false);

#ifdef DEBUG
	tt()=_Mstddev;
	tt.write("Msig_mic.spi");
#endif


}

void AutoPicker::peakSearch(const MultidimArray<RFLOAT> &Mfom, const MultidimArray<RFLOAT> &Mpsi, int iref,
		int skip_side, std::vector<Peak> &peaks)
{

	peaks.clear();
	Peak peak;
	peak.ref = iref;

	// Skip the pixels along the side of the micrograph!
	// At least 1, so dont have to check for the borders!
	skip_side = XMIPP_MAX(1, skip_side);
	for (int i = FIRST_XMIPP_INDEX(micrograph_ysize) + skip_side; i <= LAST_XMIPP_INDEX(micrograph_ysize) - skip_side; i++)
	{
		for (int j = FIRST_XMIPP_INDEX(micrograph_xsize) + skip_side; j <= LAST_XMIPP_INDEX(micrograph_xsize) - skip_side; j++)
		{

			RFLOAT myval = A2D_ELEM(Mfom, i, j);
			// check if this element is above the threshold
			if (myval  >= min_fraction_expected_Pratio)
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
				peak.x = j - FIRST_XMIPP_INDEX(micrograph_xsize);
				peak.y = i - FIRST_XMIPP_INDEX(micrograph_ysize);
				peak.psi = A2D_ELEM(Mpsi, i, j);
				peak.fom = A2D_ELEM(Mfom, i, j);
				peak.relative_fom = myval;
				peaks.push_back(peak);
			}
		}
	}

}

void AutoPicker::prunePeakClusters(std::vector<Peak> &peaks, int min_distance)
{
	int mind2 = min_distance*min_distance;
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
				if (dx*dx + dy*dy < particle_radius2)
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

void AutoPicker::removeTooCloselyNeighbouringPeaks(std::vector<Peak> &peaks, int min_distance)
{
	// Now only keep those peaks that are at least min_particle_distance number of pixels from any other peak
	std::vector<Peak> pruned_peaks;
	int mind2 = min_distance*min_distance;
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
				if ( d2 < my_mind2 )
					my_mind2 = d2;
			}
		}
		if (my_mind2 > mind2)
			pruned_peaks.push_back(peaks[ipeak]);
	}

	// Set the pruned peaks back into the input vector
	peaks = pruned_peaks;

}



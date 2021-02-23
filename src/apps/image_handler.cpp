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

#include <src/image.h>
#include <src/funcs.h>
#include <src/args.h>
#include <src/fftw.h>
#include <src/time.h>
#include <src/symmetries.h>
#include <src/jaz/single_particle/obs_model.h>
#ifdef HAVE_PNG
#include <src/jaz/gravis/tImage.h>
#endif

#include <map>

class image_handler_parameters
{
	public:
   	FileName fn_in, fn_out, fn_sel, fn_img, fn_sym, fn_sub, fn_mult, fn_div, fn_add, fn_subtract, fn_mask, fn_fsc, fn_adjust_power, fn_correct_ampl, fn_fourfilter, fn_cosDPhi;
	int bin_avg, avg_first, avg_last, edge_x0, edge_xF, edge_y0, edge_yF, filter_edge_width, new_box, minr_ampl_corr, my_new_box_size;
	bool do_add_edge, do_invert_hand, do_flipXY, do_flipmXY, do_flipZ, do_flipX, do_flipY, do_shiftCOM, do_stats, do_calc_com, do_avg_ampl, do_avg_ampl2, do_avg_ampl2_ali, do_average, do_remove_nan, do_average_all_frames, do_power, do_ignore_optics, do_optimise_scale_subtract, write_float16;
	RFLOAT multiply_constant, divide_constant, add_constant, subtract_constant, threshold_above, threshold_below, angpix, requested_angpix, real_angpix, force_header_angpix, lowpass, highpass, logfilter, bfactor, shift_x, shift_y, shift_z, replace_nan, randomize_at, optimise_bfactor_subtract;
	// PNG options
	RFLOAT minval, maxval, sigma_contrast;
	int color_scheme; // There is a global variable called colour_scheme in displayer.h!

	std::string directional;
   	int verb;
	// I/O Parser
	IOParser parser;
	ObservationModel obsModel;

	Image<RFLOAT> Iout;
	Image<RFLOAT> Iop;
	Image<RFLOAT> Imask;
	MultidimArray<RFLOAT> avg_ampl;
	MetaDataTable MD;
	FourierTransformer transformer;
	std::map<FileName, long int> n_images;

	// Image size
	int xdim, ydim, zdim;
	long int ndim;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
		fn_in = parser.getOption("--i", "Input STAR file, image (.mrc) or movie/stack (.mrcs)");
		fn_out = parser.getOption("--o", "Output name (for STAR-input: insert this string before each image's extension)", "");
		write_float16  = parser.checkOption("--float16", "Write in half-precision 16 bit floating point numbers (MRC mode 12), instead of 32 bit (MRC mode 0).");

		int cst_section = parser.addSection("image-by-constant operations");
		multiply_constant = textToFloat(parser.getOption("--multiply_constant", "Multiply the image(s) pixel values by this constant", "1"));
		divide_constant = textToFloat(parser.getOption("--divide_constant", "Divide the image(s) pixel values by this constant", "1"));
		add_constant = textToFloat(parser.getOption("--add_constant", "Add this constant to the image(s) pixel values", "0."));
		subtract_constant = textToFloat(parser.getOption("--subtract_constant", "Subtract this constant from the image(s) pixel values", "0."));
		threshold_above = textToFloat(parser.getOption("--threshold_above", "Set all values higher than this value to this value", "999."));
		threshold_below = textToFloat(parser.getOption("--threshold_below", "Set all values lower than this value to this value", "-999."));

		int img_section = parser.addSection("image-by-image operations");
		fn_mult = parser.getOption("--multiply", "Multiply input image(s) by the pixel values in this image", "");
		fn_div = parser.getOption("--divide", "Divide input image(s) by the pixel values in this image", "");
		fn_add = parser.getOption("--add", "Add the pixel values in this image to the input image(s) ", "");
		fn_subtract = parser.getOption("--subtract", "Subtract the pixel values in this image to the input image(s) ", "");
		fn_fsc = parser.getOption("--fsc", "Calculate FSC curve of the input image with this image", "");
		do_power = parser.checkOption("--power", "Calculate power spectrum (|F|^2) of the input image");
		fn_adjust_power = parser.getOption("--adjust_power", "Adjust the power spectrum of the input image to be the same as this image ", "");
		fn_fourfilter = parser.getOption("--fourier_filter", "Multiply the Fourier transform of the input image(s) with this one image ", "");

		int subtract_section = parser.addSection("additional subtract options");
		do_optimise_scale_subtract = parser.checkOption("--optimise_scale_subtract", "Optimise scale between maps before subtraction?");
		optimise_bfactor_subtract = textToFloat(parser.getOption("--optimise_bfactor_subtract", "Search range for relative B-factor for subtraction (in A^2)", "0."));
		fn_mask = parser.getOption("--mask_optimise_subtract", "Use only voxels in this mask to optimise scale for subtraction", "");

		int four_section = parser.addSection("per-image operations");
		do_stats = parser.checkOption("--stats", "Calculate per-image statistics?");
		do_calc_com = parser.checkOption("--com", "Calculate center of mass?");
		bfactor = textToFloat(parser.getOption("--bfactor", "Apply a B-factor (in A^2)", "0."));
		lowpass = textToFloat(parser.getOption("--lowpass", "Low-pass filter frequency (in A)", "-1."));
		highpass = textToFloat(parser.getOption("--highpass", "High-pass filter frequency (in A)", "-1."));
		directional = parser.getOption("--directional", "Directionality of low-pass filter frequency ('X', 'Y' or 'Z', default non-directional)", "");
		logfilter = textToFloat(parser.getOption("--LoG", "Diameter for optimal response of Laplacian of Gaussian filter (in A)", "-1."));
		angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in A)", "-1"));
		requested_angpix = textToFloat(parser.getOption("--rescale_angpix", "Scale input image(s) to this new pixel size (in A)", "-1."));
		real_angpix = -1;
		force_header_angpix = textToFloat(parser.getOption("--force_header_angpix", "Change the pixel size in the header (in A). Without --rescale_angpix, the image is not scaled.", "-1."));
		new_box = textToInteger(parser.getOption("--new_box", "Resize the image(s) to this new box size (in pixel) ", "-1"));
		filter_edge_width = textToInteger(parser.getOption("--filter_edge_width", "Width of the raised cosine on the low/high-pass filter edge (in resolution shells)", "2"));
		do_flipX = parser.checkOption("--flipX", "Flip (mirror) a 2D image or 3D map in the X-direction?");
		do_flipY = parser.checkOption("--flipY", "Flip (mirror) a 2D image or 3D map in the Y-direction?");
		do_flipZ = parser.checkOption("--flipZ", "Flip (mirror) a 3D map in the Z-direction?");
		do_invert_hand = parser.checkOption("--invert_hand", "Invert hand by flipping X? Similar to flipX, but preserves the symmetry origin. Edge pixels are wrapped around.");
		do_shiftCOM = parser.checkOption("--shift_com", "Shift image(s) to their center-of-mass (only on positive pixel values)");
		shift_x = textToFloat(parser.getOption("--shift_x", "Shift images this many pixels in the X-direction", "0."));
		shift_y = textToFloat(parser.getOption("--shift_y", "Shift images this many pixels in the Y-direction", "0."));
		shift_z = textToFloat(parser.getOption("--shift_z", "Shift images this many pixels in the Z-direction", "0."));
		do_avg_ampl = parser.checkOption("--avg_ampl", "Calculate average amplitude spectrum for all images?");
		do_avg_ampl2 = parser.checkOption("--avg_ampl2", "Calculate average amplitude spectrum for all images?");
		do_avg_ampl2_ali = parser.checkOption("--avg_ampl2_ali", "Calculate average amplitude spectrum for all aligned images?");
		do_average = parser.checkOption("--average", "Calculate average of all images (without alignment)");
		fn_correct_ampl = parser.getOption("--correct_avg_ampl", "Correct all images with this average amplitude spectrum", "");
		minr_ampl_corr = textToInteger(parser.getOption("--minr_ampl_corr", "Minimum radius (in Fourier pixels) to apply average amplitudes", "0"));
		do_remove_nan = parser.checkOption("--remove_nan", "Replace non-numerical values (NaN, inf, etc) in the image(s)");
		replace_nan = textToFloat(parser.getOption("--replace_nan", "Replace non-numerical values (NaN, inf, etc) with this value", "0"));
		randomize_at = textToFloat(parser.getOption("--phase_randomise", "Randomise phases beyond this resolution (in Angstroms)", "-1"));

		int three_d_section = parser.addSection("3D operations");
		fn_sym = parser.getOption("--sym", "Symmetrise 3D map with this point group (e.g. D6)", "");

		int preprocess_section = parser.addSection("2D-micrograph (or movie) operations");
		do_flipXY = parser.checkOption("--flipXY", "Flip the image(s) in the XY direction?");
		do_flipmXY = parser.checkOption("--flipmXY", "Flip the image(s) in the -XY direction?");
		do_add_edge = parser.checkOption("--add_edge", "Add a barcode-like edge to the micrograph/movie frames?");
		edge_x0 = textToInteger(parser.getOption("--edge_x0", "Pixel column to be used for the left edge", "0"));
		edge_y0 = textToInteger(parser.getOption("--edge_y0", "Pixel row to be used for the top edge", "0"));
		edge_xF = textToInteger(parser.getOption("--edge_xF", "Pixel column to be used for the right edge", "4095"));
		edge_yF = textToInteger(parser.getOption("--edge_yF", "Pixel row to be used for the bottom edge", "4095"));

		int avg_section = parser.addSection("Movie-frame averaging options");
		bin_avg = textToInteger(parser.getOption("--avg_bin", "Width (in frames) for binning average, i.e. of every so-many frames", "-1"));
		avg_first = textToInteger(parser.getOption("--avg_first", "First frame to include in averaging", "-1"));
		avg_last = textToInteger(parser.getOption("--avg_last", "Last frame to include in averaging", "-1"));
		do_average_all_frames = parser.checkOption("--average_all_movie_frames", "Average all movie frames of all movies in the input STAR file.");

		int png_section = parser.addSection("PNG options");
		minval = textToFloat(parser.getOption("--black", "Pixel value for black (default is auto-contrast)", "0"));
		maxval = textToFloat(parser.getOption("--white", "Pixel value for white (default is auto-contrast)", "0"));
		sigma_contrast  = textToFloat(parser.getOption("--sigma_contrast", "Set white and black pixel values this many times the image stddev from the mean", "0"));
		if (parser.checkOption("--colour_fire", "Show images in black-grey-white-red colour scheme (highlight high signal)?")) color_scheme = BLACKGREYREDSCALE;
		else if (parser.checkOption("--colour_ice", "Show images in blue-black-grey-white colour scheme (highlight low signal)?")) color_scheme = BLUEGREYWHITESCALE;
		else if (parser.checkOption("--colour_fire-n-ice", "Show images in blue-grey-red colour scheme (highlight high&low signal)?")) color_scheme = BLUEGREYREDSCALE;
		else if (parser.checkOption("--colour_rainbow", "Show images in cyan-blue-black-red-yellow colour scheme?")) color_scheme = RAINBOWSCALE;
		else if (parser.checkOption("--colour_difference", "Show images in cyan-blue-black-red-yellow colour scheme (for difference images)?")) color_scheme = CYANBLACKYELLOWSCALE;
		else color_scheme = GREYSCALE;

		// Hidden
		fn_cosDPhi = getParameter(argc, argv, "--cos_dphi", "");
		// Check for errors in the command-line option
		if (parser.checkForErrors())
			REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

		verb = (do_stats || do_calc_com || fn_fsc !="" || fn_cosDPhi != "" | do_power) ? 0 : 1;

		if (fn_out == "" && verb == 1)
			REPORT_ERROR("Please specify the output file name with --o.");
	}

	void perImageOperations(Image<RFLOAT> &Iin, FileName &my_fn_out, RFLOAT psi = 0.)
	{
		Image<RFLOAT> Iout;
		Iout().resize(Iin());

		bool isPNG = FileName(my_fn_out.getExtension()).toLowercase() == "png";
		if (isPNG && (ZSIZE(Iout()) > 1 || NSIZE(Iout()) > 1))
			REPORT_ERROR("You can only write a 2D image to a PNG file.");

		if (angpix < 0 && (requested_angpix > 0 || fn_fsc != "" || randomize_at > 0 ||
		                   do_power || fn_cosDPhi != "" || fn_correct_ampl != "" ||
		                   fabs(bfactor) > 0 || logfilter > 0 || lowpass > 0 || highpass > 0 || fabs(optimise_bfactor_subtract) > 0))
		{
			angpix = Iin.samplingRateX();
			std::cerr << "WARNING: You did not specify --angpix. The pixel size in the image header, " << angpix << " A/px, is used." << std::endl;
		}

		if (do_add_edge)
		{
			// Treat X-boundaries
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
			{
				if (j < edge_x0)
					DIRECT_A2D_ELEM(Iin(), i, j) = DIRECT_A2D_ELEM(Iin(), i, edge_x0);
				else if (j > edge_xF)
					DIRECT_A2D_ELEM(Iin(), i, j) = DIRECT_A2D_ELEM(Iin(), i, edge_xF);
			}
			// Treat Y-boundaries
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
			{
				if (i < edge_y0)
					DIRECT_A2D_ELEM(Iin(), i, j) = DIRECT_A2D_ELEM(Iin(), edge_y0, j);
				else if (i > edge_yF)
					DIRECT_A2D_ELEM(Iin(), i, j) = DIRECT_A2D_ELEM(Iin(), edge_yF, j);
			}
		}

		// Flipping: this needs to be done from Iin to Iout (i.e. can't be done on-line on Iout only!)
		if (do_flipXY)
		{
			// Flip X/Y
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
			{
				DIRECT_A2D_ELEM(Iout(), i, j) = DIRECT_A2D_ELEM(Iin(), j, i);

			}
		}
		else if (do_flipmXY)
		{
			// Flip mX/Y
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
			{
				DIRECT_A2D_ELEM(Iout(), i, j) = DIRECT_A2D_ELEM(Iin(), XSIZE(Iin()) - 1 - j, YSIZE(Iin()) - 1 - i);
			}
		}
		else
		{
			Iout = Iin;
		}

		// From here on also 3D options
		if (do_remove_nan)
		{
			Iout().setXmippOrigin();
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iout())
			{
				if (std::isnan(DIRECT_A3D_ELEM(Iout(), k, i, j)) || std::isinf(DIRECT_A3D_ELEM(Iout(), k, i, j)))
					DIRECT_A3D_ELEM(Iout(), k, i, j) = replace_nan;
			}
		}

		if (randomize_at > 0.)
		{
			int iran = XSIZE(Iin())* angpix / randomize_at;
			Iout = Iin;
			randomizePhasesBeyond(Iout(), iran);
		}
		if (fabs(multiply_constant - 1.) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) *= multiply_constant;
			}
		}
		else if (fabs(divide_constant - 1.) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) /= divide_constant;
			}
		}
		else if (fabs(add_constant) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) += add_constant;
			}
		}
		else if (fabs(subtract_constant) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) -= subtract_constant;
			}
		}
		else if (fn_mult != "")
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) *= DIRECT_A3D_ELEM(Iop(), k, i, j);
			}
		}
		else if (fn_div != "")
		{
			bool is_first = true;
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				if (ABS(DIRECT_A3D_ELEM(Iop(), k, i, j)) < 1e-10)
				{
					if (is_first)
					{
						std::cout << "Warning: ignore very small pixel values in divide image..." << std::endl;
						is_first = false;
					}
					DIRECT_A3D_ELEM(Iout(), k, i, j) = 0.;
				}
				else
					DIRECT_A3D_ELEM(Iout(), k, i, j) /= DIRECT_A3D_ELEM(Iop(), k, i, j);
			}
		}
		else if (fn_add != "")
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) += DIRECT_A3D_ELEM(Iop(), k, i, j);
			}
		}
		else if (fn_subtract != "")
		{
			RFLOAT my_scale = 1., best_diff2 ;
			if (do_optimise_scale_subtract)
			{
				if (fn_mask == "")
				{
					Imask(). resize(Iop());
					Imask().initConstant(1.);
				}

				if (optimise_bfactor_subtract > 0.)
				{
					MultidimArray< Complex > FTop, FTop_bfac;
					FourierTransformer transformer;
					MultidimArray<RFLOAT> Isharp(Iop());
					transformer.FourierTransform(Iop(), FTop);

					RFLOAT my_bfac, smallest_diff2=99.e99;
					for (RFLOAT bfac = -optimise_bfactor_subtract; bfac <= optimise_bfactor_subtract; bfac+= 10.)
					{
						FTop_bfac = FTop;
						applyBFactorToMap(FTop_bfac, XSIZE(Iop()), bfac, angpix);
						transformer.inverseFourierTransform(FTop_bfac, Isharp);
						RFLOAT scale, diff2;

						RFLOAT sum_aa = 0., sum_xa = 0., sum_xx = 0.;
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iin())
						{
							RFLOAT w = DIRECT_MULTIDIM_ELEM(Imask(), n) * DIRECT_MULTIDIM_ELEM(Imask(), n);
							RFLOAT x = DIRECT_MULTIDIM_ELEM(Iin(), n);
							RFLOAT a = DIRECT_MULTIDIM_ELEM(Isharp, n);
							sum_aa += w*a*a;
							sum_xa += w*x*a;
							sum_xx += w*x*x;
						}
						scale = sum_xa/sum_aa;

						diff2 = 0.;
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iin())
						{
							RFLOAT w = DIRECT_MULTIDIM_ELEM(Imask(), n);
							RFLOAT x = DIRECT_MULTIDIM_ELEM(Iin(), n);
							RFLOAT a = DIRECT_MULTIDIM_ELEM(Isharp, n);
							diff2 += w * w * (x - scale * a) * (x - scale * a);
						}
						if (diff2 < smallest_diff2)
						{
							smallest_diff2 = diff2;
							my_bfac = bfac;
							my_scale = scale;
						}
					}
					std::cout << " Optimised bfactor = " << my_bfac << "; optimised scale = " << my_scale << std::endl;
					applyBFactorToMap(FTop, XSIZE(Iop()), my_bfac, angpix);
					transformer.inverseFourierTransform(FTop, Iop());

				}
				else
				{
					RFLOAT sum_aa = 0., sum_xa = 0.;
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iin())
					{
						RFLOAT w = DIRECT_MULTIDIM_ELEM(Imask(), n);
						RFLOAT x = DIRECT_MULTIDIM_ELEM(Iin(), n);
						RFLOAT a = DIRECT_MULTIDIM_ELEM(Iop(), n);
						sum_aa += w*w*a*a;
						sum_xa += w*w*x*a;
					}
					my_scale = sum_xa/sum_aa;
					std::cout << " Optimised scale = " << my_scale << std::endl;

				}
			}

			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) -= my_scale * DIRECT_A3D_ELEM(Iop(), k, i, j);
			}
		}
		else if (fn_fsc != "")
		{
			MultidimArray<RFLOAT> fsc;
			MetaDataTable MDfsc;
			getFSC(Iout(), Iop(), fsc);
			MDfsc.setName("fsc");
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fsc)
			{
				MDfsc.addObject();
				RFLOAT res = (i > 0) ? (XSIZE(Iout()) * angpix / (RFLOAT)i) : 999.;
				MDfsc.setValue(EMDL_SPECTRAL_IDX, (int)i);
				MDfsc.setValue(EMDL_RESOLUTION, 1./res);
				MDfsc.setValue(EMDL_RESOLUTION_ANGSTROM, res);
				MDfsc.setValue(EMDL_POSTPROCESS_FSC_GENERAL, DIRECT_A1D_ELEM(fsc, i));
			}
			MDfsc.write(std::cout);
		}
		else if (do_power)
		{
			MultidimArray<RFLOAT> spectrum;
			getSpectrum(Iout(), spectrum, POWER_SPECTRUM);
			MetaDataTable MDpower;
			MDpower.setName("power");
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(spectrum)
			{
				if (i > XSIZE(Iout()) / 2 + 1) break; // getSpectrum returns beyond Nyquist!!

				MDpower.addObject();
				RFLOAT res = (i > 0) ? (XSIZE(Iout()) * angpix / (RFLOAT)i) : 999.;
				MDpower.setValue(EMDL_SPECTRAL_IDX, (int)i);
				MDpower.setValue(EMDL_RESOLUTION, 1./res);
				MDpower.setValue(EMDL_RESOLUTION_ANGSTROM, res);
				MDpower.setValue(EMDL_MLMODEL_POWER_REF, DIRECT_A1D_ELEM(spectrum, i));
			}
			MDpower.write(std::cout);
		}
		else if (fn_adjust_power != "")
		{
			MultidimArray<RFLOAT> spectrum;
			getSpectrum(Iop(), spectrum, AMPLITUDE_SPECTRUM);
			adaptSpectrum(Iin(), Iout(), spectrum, AMPLITUDE_SPECTRUM);
		}
		else if (fn_cosDPhi != "")
		{
			MultidimArray<RFLOAT> cosDPhi;
			MetaDataTable MDcos;

			MultidimArray< Complex > FT1, FT2;
			FourierTransformer transformer;
			transformer.FourierTransform(Iout(), FT1);
			transformer.FourierTransform(Iop(), FT2);

			getCosDeltaPhase(FT1, FT2, cosDPhi);
			MDcos.setName("cos");
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(cosDPhi)
			{
				MDcos.addObject();
				RFLOAT res = (i > 0) ? (XSIZE(Iout()) * angpix / (RFLOAT)i) : 999.;
				MDcos.setValue(EMDL_SPECTRAL_IDX, (int)i);
				MDcos.setValue(EMDL_RESOLUTION, 1./res);
				MDcos.setValue(EMDL_RESOLUTION_ANGSTROM, res);
				MDcos.setValue(EMDL_POSTPROCESS_FSC_GENERAL, DIRECT_A1D_ELEM(cosDPhi, i));
			}
			MDcos.write(std::cout);
		}
		else if (fn_correct_ampl != "")
		{
			MultidimArray<Complex> FT;
			transformer.FourierTransform(Iin(), FT, false);
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FT)
			{
				DIRECT_MULTIDIM_ELEM(FT, n) /=  DIRECT_MULTIDIM_ELEM(avg_ampl, n);
			}
			transformer.inverseFourierTransform();
			Iout = Iin;
		}
		else if (fn_fourfilter != "")
		{
			MultidimArray<Complex> FT;
			transformer.FourierTransform(Iin(), FT, false);

			// Note: only 2D rotations are done! 3D application assumes zero rot and tilt!
			Matrix2D<RFLOAT> A;
			rotation2DMatrix(psi, A);

			Iop().setXmippOrigin();
			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT)
			{
				int jpp = ROUND(jp * A(0, 0) + ip * A(0, 1));
				int ipp = ROUND(jp * A(1, 0) + ip * A(1, 1));
				int kpp = kp;
				RFLOAT fil;
				if (jpp >= STARTINGX(Iop()) && jpp <= FINISHINGX(Iop()) && ipp >= STARTINGY(Iop()) && ipp <= FINISHINGY(Iop()))
					fil = A3D_ELEM(Iop(), kpp, ipp, jpp);
				else
					fil = 0.;
				DIRECT_A3D_ELEM(FT, k, i, j) *=  fil;
			}
			transformer.inverseFourierTransform();
			Iout = Iin;
		}

		if (fabs(bfactor) > 0.)
			applyBFactorToMap(Iout(), bfactor, angpix);

		if (logfilter > 0.)
		{
			LoGFilterMap(Iout(), logfilter, angpix);
			RFLOAT avg, stddev, minval, maxval;
			//Iout().statisticsAdjust(0,1);
		}

		if (lowpass > 0.)
		{
			if (directional != "")
				directionalFilterMap(Iout(), lowpass, angpix, directional, filter_edge_width);
			else
				lowPassFilterMap(Iout(), lowpass, angpix, filter_edge_width);
		}

		if (highpass > 0.)
			highPassFilterMap(Iout(), highpass, angpix, filter_edge_width);

		if (do_flipX)
		{
			// For input:  0, 1, 2, 3, 4, 5 (XSIZE = 6)
			// This gives: 5, 4, 3, 2, 1, 0
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) = A3D_ELEM(Iin(), k, i, XSIZE(Iin()) - 1 - j);
			}
		}
		else if (do_flipY)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) = A3D_ELEM(Iin(), k, YSIZE(Iin()) - 1 - i, j);
			}
		}
		else if (do_flipZ)
		{
			if (ZSIZE(Iout()) < 2)
				REPORT_ERROR("ERROR: this is not a 3D map, so cannot be flipped in Z");

			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				DIRECT_A3D_ELEM(Iout(), k, i, j) = A3D_ELEM(Iin(), ZSIZE(Iin()) - 1 - k, i, j);
			}
		}
		else if (do_invert_hand)
		{
			// For input:  0, 1, 2, 3, 4, 5 (XSIZE = 6)
			// This gives: 0, 5, 4, 3, 2, 1
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iin())
			{
				long int dest_x = (j == 0) ? 0 : (XSIZE(Iin()) - j);
				DIRECT_A3D_ELEM(Iout(), k, i, j) = A3D_ELEM(Iin(), k, i, dest_x);
			}
		}

		// Shifting
		if (do_shiftCOM)
			selfTranslateCenterOfMassToCenter(Iout(), DONT_WRAP, true); // verbose=true!
		else if (fabs(shift_x) > 0. || fabs(shift_y) > 0. || fabs(shift_z) > 0.)
		{
			Matrix1D<RFLOAT> shift(2);
			XX(shift) = shift_x;
			YY(shift) = shift_y;
			if (zdim > 1)
			{
				shift.resize(3);
				ZZ(shift) = shift_z;
			}
			selfTranslate(Iout(), shift, DONT_WRAP);
		}

		// Re-scale
		if (requested_angpix > 0.)
		{
			int oldxsize = XSIZE(Iout());
			int oldysize = YSIZE(Iout());
			int oldsize = oldxsize;
			if (oldxsize != oldysize && Iout().getDim() == 2)
			{
				oldsize = XMIPP_MAX( oldxsize, oldysize );
				Iout().setXmippOrigin();
				Iout().window(FIRST_XMIPP_INDEX(oldsize), FIRST_XMIPP_INDEX(oldsize),
			 	              LAST_XMIPP_INDEX(oldsize),  LAST_XMIPP_INDEX(oldsize));
			}

			int newsize = ROUND(oldsize * (angpix / requested_angpix));
			newsize -= newsize % 2; //make even in case it is not already

			real_angpix = oldsize * angpix / newsize;
			if (fabs(real_angpix - requested_angpix) / requested_angpix > 0.001)
				std::cerr << "WARNING: Although the requested pixel size (--rescale_angpix) is " << requested_angpix << " A/px, the actual pixel size will be " << real_angpix << " A/px due to rounding of the box size to an even number. The latter value is set to the image header. You can overwrite the header pixel size by --force_header_angpix." << std::endl;

			resizeMap(Iout(), newsize);
			my_new_box_size = newsize;

			if (oldxsize != oldysize && Iout().getDim() == 2)
			{
				int newxsize = ROUND(oldxsize * (angpix / real_angpix));
				int newysize = ROUND(oldysize * (angpix / real_angpix));;
				newxsize -= newxsize%2; //make even in case it is not already
				newysize -= newysize%2; //make even in case it is not already
				Iout().setXmippOrigin();
				Iout().window(FIRST_XMIPP_INDEX(newysize), FIRST_XMIPP_INDEX(newxsize),
				              LAST_XMIPP_INDEX(newysize),  LAST_XMIPP_INDEX(newxsize));
			}

			// Also reset the sampling rate in the header
			Iout.setSamplingRateInHeader(real_angpix);
		}

		// Re-window
		if (new_box > 0 && XSIZE(Iout()) != new_box)
		{
			Iout().setXmippOrigin();
			if (Iout().getDim() == 2)
			{
				Iout().window(FIRST_XMIPP_INDEX(new_box), FIRST_XMIPP_INDEX(new_box),
						   LAST_XMIPP_INDEX(new_box),  LAST_XMIPP_INDEX(new_box));
			}
			else if (Iout().getDim() == 3)
			{
				Iout().window(FIRST_XMIPP_INDEX(new_box), FIRST_XMIPP_INDEX(new_box), FIRST_XMIPP_INDEX(new_box),
						   LAST_XMIPP_INDEX(new_box),  LAST_XMIPP_INDEX(new_box),  LAST_XMIPP_INDEX(new_box));
			}
			my_new_box_size = new_box;
		}

		if (fn_sym != "")
			symmetriseMap(Iout(), fn_sym);

		// Thresholding (can be done after any other operation)
		if (fabs(threshold_above - 999.) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iout())
			{
				if (DIRECT_A3D_ELEM(Iout(), k, i, j) > threshold_above)
					DIRECT_A3D_ELEM(Iout(), k, i, j) = threshold_above;
			}
		}
		if (fabs(threshold_below + 999.) > 0.)
		{
			FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Iout())
			{
				if (DIRECT_A3D_ELEM(Iout(), k, i, j) < threshold_below)
					DIRECT_A3D_ELEM(Iout(), k, i, j) = threshold_below;
			}
		}

		if (force_header_angpix > 0)
		{
			Iout.setSamplingRateInHeader(force_header_angpix);
			std::cout << "As requested by --force_header_angpix, the pixel size in the image header is set to " << force_header_angpix << " A/px." << std::endl;
		}

		// Write out the result
		// Check whether fn_out has an "@": if so REPLACE the corresponding frame in the output stack!
		long int n;
		FileName fn_tmp;
		my_fn_out.decompose(n, fn_tmp);
		n--;
		if (!isPNG)
		{
			if (n >= 0) // This is a stack...
			{

				// The following assumes the images in the stack come ordered...
				if (n == 0)
					Iout.write(fn_tmp, n, true, WRITE_OVERWRITE, write_float16 ? Float16: Float); // make a new stack
				else
					Iout.write(fn_tmp, n, true, WRITE_APPEND, write_float16 ? Float16: Float);
			}
			else
				Iout.write(my_fn_out, -1, false, WRITE_OVERWRITE, write_float16 ? Float16: Float);
		}
		else
		{
#ifdef HAVE_PNG
			RFLOAT this_minval = minval, this_maxval = maxval; // User setting
			getImageContrast(Iout(), this_minval, this_maxval, sigma_contrast); // Update if neecssary
			const RFLOAT range = this_maxval - this_minval;
			const RFLOAT step = range / 255;

			gravis::tImage<gravis::bRGB> pngOut(XSIZE(Iout()), YSIZE(Iout()));
			pngOut.fill(gravis::bRGB(0));

			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iout())
			{
				const unsigned char val = FLOOR((DIRECT_MULTIDIM_ELEM(Iout(), n) - this_minval) / step);
				unsigned char r, g, b;
				greyToRGB(color_scheme, val, r, g, b);
				pngOut[n] = gravis::bRGB(r, g, b);
			}
			pngOut.writePNG(my_fn_out);
#else
			REPORT_ERROR("You cannot write PNG images because libPNG was not linked during compilation.");
#endif
		}
	}

	void run()
	{
		my_new_box_size = -1;

		long int slice_id;
		std::string fn_stem;
		fn_in.decompose(slice_id, fn_stem);
		bool input_is_stack = (fn_in.getExtension() == "mrcs" || fn_in.getExtension() == "tif" || fn_in.getExtension() == "tiff") && (slice_id == -1);
		bool input_is_star = (fn_in.getExtension() == "star");
		// By default: write single output images

		// Get a MetaDataTable
		if (input_is_star)
		{
			do_ignore_optics = false;
			ObservationModel::loadSafely(fn_in, obsModel, MD, "discover", verb, false); // false means don't die upon failure
			if (obsModel.opticsMdt.numberOfObjects() == 0)
			{
				do_ignore_optics = true;
				std::cout << " + WARNING: reading input STAR file without optics groups ..." << std::endl;
				MD.read(fn_in);
			}
			if (fn_out.getExtension() != "mrcs")
				std::cout << "NOTE: the input (--i) is a STAR file but the output (--o) does not have .mrcs extension. The output is treated as a suffix, not a path." << std::endl;
			FileName fn_img;
			MD.getValue(EMDL_IMAGE_NAME, fn_img, 0);
			fn_img.decompose(slice_id, fn_stem);
			input_is_stack = (fn_in.getExtension() == "mrcs" || fn_in.getExtension() == "tif" || fn_in.getExtension() == "tiff") && (slice_id == -1);
		}
		else if (input_is_stack)
		{
			if (bin_avg > 0 || (avg_first >= 0 && avg_last >= 0))
			{
				MD.addObject();
				MD.setValue(EMDL_IMAGE_NAME, fn_in);
			}
			else
			{
				// Read the header to get the number of images inside the stack and generate that many lines in the MD
				Image<RFLOAT> tmp;
				FileName fn_tmp;
				tmp.read(fn_in, false); //false means do not read image now, only header
				for (int i = 1; i <= NSIZE(tmp()); i++)
				{
					MD.addObject();
					fn_tmp.compose(i, fn_in);
					MD.setValue(EMDL_IMAGE_NAME, fn_tmp);
				}
			}
		}
		else
		{
			// Just individual image input
			MD.addObject();
			MD.setValue(EMDL_IMAGE_NAME, fn_in);
		}

		int i_img = 0;
		time_config();
   		if (verb > 0)
   			init_progress_bar(MD.numberOfObjects());

		bool do_md_out = false;
   		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MD)
		{
			FileName fn_img;
			if (do_average_all_frames)
			{
				MD.getValue(EMDL_MICROGRAPH_MOVIE_NAME, fn_img);
			}
			else
			{
				MD.getValue(EMDL_IMAGE_NAME, fn_img);
			}

			// For fourfilter...
			RFLOAT psi;
			if (!MD.getValue(EMDL_ORIENT_PSI, psi))
				psi =0.;

			Image<RFLOAT> Iin;
			// Initialise for the first image
			if (i_img == 0)
			{
				Image<RFLOAT> Ihead;
				Ihead.read(fn_img, false);
				Ihead.getDimensions(xdim, ydim, zdim, ndim);

				if (zdim > 1 && (do_add_edge || do_flipXY || do_flipmXY))
					REPORT_ERROR("ERROR: you cannot perform 2D operations like --add_edge, --flipXY or --flipmXY on 3D maps. If you intended to operate on a movie, use .mrcs extensions for stacks!");

				if (zdim > 1 && (bin_avg > 0 || (avg_first >= 0 && avg_last >= 0)))
					REPORT_ERROR("ERROR: you cannot perform movie-averaging operations on 3D maps. If you intended to operate on a movie, use .mrcs extensions for stacks!");

				if (fn_mult != "")
					Iop.read(fn_mult);
				else if (fn_div != "")
					Iop.read(fn_div);
				else if (fn_add != "")
					Iop.read(fn_add);
				else if (fn_subtract != "")
				{
					Iop.read(fn_subtract);
					if (do_optimise_scale_subtract && fn_mask != "") Imask.read(fn_mask);
				}
				else if (fn_fsc != "")
					Iop.read(fn_fsc);
				else if (fn_cosDPhi != "")
					Iop.read(fn_cosDPhi);
				else if (fn_adjust_power != "")
					Iop.read(fn_adjust_power);
				else if (fn_fourfilter != "")
					Iop.read(fn_fourfilter);
				else if (fn_correct_ampl != "")
				{
					Iop.read(fn_correct_ampl);

					// Calculate by the radial average in the Fourier domain
					MultidimArray<RFLOAT> spectrum, count;
					spectrum.initZeros(YSIZE(Iop()));
					count.initZeros(YSIZE(Iop()));
					FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Iop())
					{
						long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
						spectrum(idx) += dAkij(Iop(), k, i, j);
						count(idx) += 1.;
					}
					FOR_ALL_ELEMENTS_IN_ARRAY1D(spectrum)
					{
						if (A1D_ELEM(count, i) > 0.)
							A1D_ELEM(spectrum, i) /= A1D_ELEM(count, i);
					}

					FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Iop())
					{
				    		long int idx = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
					    	if (idx > minr_ampl_corr)
					    		dAkij(Iop(), k, i, j) /= spectrum(idx);
					    	else
					    		dAkij(Iop(), k, i, j) = 1.;
					}
					avg_ampl = Iop();
				}

				if (fn_mult != "" || fn_div != "" || fn_add != "" || fn_subtract != "" || fn_fsc != "" || fn_adjust_power != "" ||fn_fourfilter != "")
					if (XSIZE(Iop()) != xdim || YSIZE(Iop()) != ydim || ZSIZE(Iop()) != zdim)
						REPORT_ERROR("Error: operate-image is not of the correct size");

				if (do_avg_ampl || do_avg_ampl2 || do_avg_ampl2_ali)
				{
					avg_ampl.initZeros(zdim, ydim, xdim/2+1);
				}
				else if (do_average || do_average_all_frames)
				{
					avg_ampl.initZeros(zdim, ydim, xdim);
				}

			}

			if (do_stats) // only write statistics to screen
			{
				Iin.read(fn_img);
				RFLOAT avg, stddev, minval, maxval, header_angpix;
				Iin().computeStats(avg, stddev, minval, maxval);
				header_angpix = Iin.samplingRateX();
				std::cout << fn_img << " : (x,y,z,n)= " << XSIZE(Iin()) << " x "<< YSIZE(Iin()) << " x "<< ZSIZE(Iin()) << " x "<< NSIZE(Iin()) << " ; avg= " << avg << " stddev= " << stddev << " minval= " <<minval << " maxval= " << maxval << "; angpix = " << header_angpix << std::endl;
			}
			else if (do_calc_com)
			{
				Matrix1D <RFLOAT> com(3);
				Iin.read(fn_img);
				Iin().setXmippOrigin();
				Iin().centerOfMass(com);
				std::cout << fn_img << " : center of mass (relative to XmippOrigin) x " << com(0);
				if (VEC_XSIZE(com) > 1) std::cout << " y " << YY(com);
				if (VEC_XSIZE(com) > 2) std::cout << " z " << ZZ(com);
				std::cout << std::endl;
			}
			else if (do_avg_ampl || do_avg_ampl2 || do_avg_ampl2_ali)
			{
				Iin.read(fn_img);

				if (do_avg_ampl2_ali)
				{
					RFLOAT xoff = 0.;
					RFLOAT yoff = 0.;
					RFLOAT psi = 0.;
					MD.getValue(EMDL_ORIENT_ORIGIN_X, xoff);
					MD.getValue(EMDL_ORIENT_ORIGIN_Y, yoff);
					MD.getValue(EMDL_ORIENT_PSI, psi);
					// Apply the actual transformation
					Matrix2D<RFLOAT> A;
					rotation2DMatrix(psi, A);
					MAT_ELEM(A,0, 2) = xoff;
					MAT_ELEM(A,1, 2) = yoff;
					selfApplyGeometry(Iin(), A, IS_NOT_INV, DONT_WRAP);
				}

				MultidimArray<Complex> FT;
				transformer.FourierTransform(Iin(), FT);

				if (do_avg_ampl)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FT)
					{
						DIRECT_MULTIDIM_ELEM(avg_ampl, n) +=  abs(DIRECT_MULTIDIM_ELEM(FT, n));
					}
				}
				else if (do_avg_ampl2 || do_avg_ampl2_ali)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(FT)
					{
						DIRECT_MULTIDIM_ELEM(avg_ampl, n) +=  norm(DIRECT_MULTIDIM_ELEM(FT, n));
					}
				}
			}
			else if (do_average)
			{
				Iin.read(fn_img);
				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iin())
				{
					DIRECT_MULTIDIM_ELEM(avg_ampl, n) +=  DIRECT_MULTIDIM_ELEM(Iin(), n);
				}
			}
			else if (do_average_all_frames)
			{
				Iin.read(fn_img);
				for (int n = 0; n < ndim; n++)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(avg_ampl)
					{
						DIRECT_A3D_ELEM(avg_ampl, k, i, j) +=  DIRECT_NZYX_ELEM(Iin(), n, k, i, j);
					}
				}
			}
			else if (bin_avg > 0 || (avg_first >= 0 && avg_last >= 0))
			{
				// movie-frame averaging operations
				int avgndim = 1;
				if (bin_avg > 0)
				{
					avgndim = ndim / bin_avg;
				}
				Image<RFLOAT> Iavg(xdim, ydim, zdim, avgndim);

				if (ndim == 1)
					REPORT_ERROR("ERROR: you are trying to perform movie-averaging options on a single image/volume");

				FileName fn_ext = fn_out.getExtension();
				if (NSIZE(Iavg()) > 1 && ( fn_ext.contains("mrc") && !fn_ext.contains("mrcs") ) )
					REPORT_ERROR("ERROR: trying to write a stack into an MRC image. Use .mrcs extensions for stacks!");

				for (long int nn = 0; nn < ndim; nn++)
				{
					Iin.read(fn_img, true, nn);
					if (bin_avg > 0)
					{
						int myframe = nn / bin_avg;
						if (myframe < avgndim)
						{
							FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(Iin())
							{
								DIRECT_NZYX_ELEM(Iavg(), myframe, 0, i, j) += DIRECT_A2D_ELEM(Iin(), i, j); // just store sum
							}
						}
					}
					else if (avg_first >= 0 && avg_last >= 0 && nn+1 >= avg_first && nn+1 <= avg_last) // add one to start counting at 1
					{
						FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iin())
						{
							DIRECT_MULTIDIM_ELEM(Iavg(), n) += DIRECT_MULTIDIM_ELEM(Iin(), n); // just store sum
						}
					}
				}
				Iavg.write(fn_out, -1, true, WRITE_OVERWRITE, write_float16 ? Float16: Float);
			}
			else
			{
				Iin.read(fn_img);
				FileName my_fn_out;

				if (fn_out.getExtension() == "mrcs" && !fn_out.contains("@"))
				{
					// current_object starts counting from 0, thus needs to be incremented.
					my_fn_out.compose(current_object + 1, fn_out);
				}
				else
				{
					if (input_is_stack)
					{
						my_fn_out = fn_img.insertBeforeExtension("_" + fn_out);
						long int dummy;
						FileName fn_tmp;
						my_fn_out.decompose(dummy, fn_tmp);
						n_images[fn_tmp]++; // this is safe. see https://stackoverflow.com/questions/16177596/stdmapstring-int-default-initialization-of-value.
						my_fn_out.compose(n_images[fn_tmp], fn_tmp);
					}
					else if (input_is_star)
					{
						my_fn_out = fn_img.insertBeforeExtension("_" + fn_out);
					}
					else
					{
						my_fn_out = fn_out;
					}
				}
				perImageOperations(Iin, my_fn_out, psi);
				do_md_out = true;
				MD.setValue(EMDL_IMAGE_NAME, my_fn_out);
			}

			i_img+=ndim;
			if (verb > 0)
				progress_bar(i_img/ndim);
		}


		if (do_avg_ampl || do_avg_ampl2 || do_avg_ampl2_ali || do_average || do_average_all_frames)
		{
			avg_ampl /= (RFLOAT)i_img;
			Iout() = avg_ampl;
			Iout.write(fn_out, -1, false, WRITE_OVERWRITE, write_float16 ? Float16: Float);
		}

		if (verb > 0)
			progress_bar(MD.numberOfObjects());

		if (do_md_out && fn_in.getExtension() == "star")
		{
			FileName fn_md_out = fn_in.insertBeforeExtension("_" + fn_out);

			if (do_ignore_optics)
			{
				MD.write(fn_md_out);
			}
			else
			{
				if (my_new_box_size > 0)
				{
					FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
					{
						obsModel.opticsMdt.setValue(EMDL_IMAGE_SIZE, my_new_box_size);
					}
				}
				if (real_angpix > 0)
				{
					FOR_ALL_OBJECTS_IN_METADATA_TABLE(obsModel.opticsMdt)
					{
						obsModel.opticsMdt.setValue(EMDL_IMAGE_PIXEL_SIZE, real_angpix);
					}
				}
				obsModel.save(MD, fn_md_out);
			}

			std::cout << " Written out new STAR file: " << fn_md_out << std::endl;
		}
	}
};

int main(int argc, char *argv[])
{
	image_handler_parameters prm;

	try
	{
		prm.read(argc, argv);

		prm.run();

	}
	catch (RelionError XE)
	{
        	//prm.usage();
	        std::cerr << XE;
        	return RELION_EXIT_FAILURE;
	}
	return RELION_EXIT_SUCCESS;
}

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

#include <src/backprojector.h>
#include <src/funcs.h>
#include <src/ctf.h>
#include <src/args.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/time.h>

class reconstruct_parameters
{
	public:
   	FileName fn_out, fn_sel, fn_img, fn_sym, fn_sub, fn_fsc, fn_debug;
	MetaDataTable DF;
	int r_max, r_min_nn, blob_order, ref_dim, interpolator, iter, nr_threads, debug_ori_size, debug_size, ctf_dim, nr_helical_asu;
	RFLOAT blob_radius, blob_alpha, angular_error, shift_error, angpix, maxres, beamtilt_x, beamtilt_y, helical_rise, helical_twist;
	bool do_ctf, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, do_fom_weighting, do_3d_rot, do_reconstruct_ctf, do_beamtilt;
	bool skip_gridding, do_reconstruct_ctf2, do_reconstruct_meas;
	float padding_factor;
	// I/O Parser
	IOParser parser;

	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{

		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("General options");
		fn_sel = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
	    fn_out = parser.getOption("--o", "Name for output reconstruction","relion.mrc");
	    fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
       	angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "1"));
       	maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
       	padding_factor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
    	nr_threads = textToInteger(parser.getOption("--j", "Number of threads to use for FFTs", "1"));

	    int ctf_section = parser.addSection("CTF options");
       	do_ctf = parser.checkOption("--ctf", "Apply CTF correction");
    	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
    	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
    	only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");
    	beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in the X-direction (in mrad)", "0."));
    	beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in the Y-direction (in mrad)", "0."));
    	do_beamtilt = (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.);

    	int helical_section = parser.addSection("Helical options");
    	nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
    	helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
    	helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

       	int expert_section = parser.addSection("Expert options");
       	fn_sub = parser.getOption("--subtract","Subtract projections of this map from the images used for reconstruction", "");
        if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation before gridding correction"))
       		interpolator = NEAREST_NEIGHBOUR;
       	else
       		interpolator = TRILINEAR;
        blob_radius   = textToFloat(parser.getOption("--blob_r", "Radius of blob for gridding interpolation", "1.9"));
        blob_order    = textToInteger(parser.getOption("--blob_m", "Order of blob for gridding interpolation", "0"));
        blob_alpha    = textToFloat(parser.getOption("--blob_a", "Alpha-value of blob for gridding interpolation", "15"));
       	iter = textToInteger(parser.getOption("--iter", "Number of gridding-correction iterations", "10"));
       	ref_dim = textToInteger(parser.getOption("--refdim", "Dimension of the reconstruction (2D or 3D)", "3"));
    	angular_error = textToFloat(parser.getOption("--angular_error", "Apply random deviations with this standard deviation (in degrees) to each of the 3 Euler angles", "0."));
    	shift_error = textToFloat(parser.getOption("--shift_error", "Apply random deviations with this standard deviation (in pixels) to each of the 2 translations", "0."));
    	do_fom_weighting = parser.checkOption("--fom_weighting", "Weight particles according to their figure-of-merit (_rlnParticleFigureOfMerit)");
    	fn_fsc = parser.getOption("--fsc", "FSC-curve for regularized reconstruction", "");
    	do_3d_rot = parser.checkOption("--3d_rot", "Perform 3D rotations instead of backprojections from 2D images");
    	ctf_dim  = textToInteger(parser.getOption("--reconstruct_ctf", "Perform a 3D reconstruction from 2D CTF-images, with the given size in pixels", "-1"));
    	do_reconstruct_ctf2 = parser.checkOption("--ctf2", "Reconstruct CTF^2 and then take the sqrt of that");
    	do_reconstruct_meas = parser.checkOption("--measured", "Fill Hermitian half of the CTF reconstruction with how often each voxel was measured.");
    	skip_gridding = parser.checkOption("--skip_gridding", "Skip gridding part of the reconstruction");
    	do_reconstruct_ctf = (ctf_dim > 0);
    	if (do_reconstruct_ctf)
    		do_ctf = false;
		fn_debug = parser.getOption("--debug", "Rootname for debug reconstruction files", "");
		debug_ori_size =  textToInteger(parser.getOption("--debug_ori_size", "Rootname for debug reconstruction files", "1"));
		debug_size =  textToInteger(parser.getOption("--debug_size", "Rootname for debug reconstruction files", "1"));


    	// Hidden
       	r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

    	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

    	// Read MetaData file, which should have the image names and their angles!
    	if (fn_debug == "")
    		DF.read(fn_sel);

     	randomize_random_generator();

     	if (do_beamtilt && ! do_ctf)
     		REPORT_ERROR("ERROR: one can only correct for beamtilt in combination with CTF correction!");

	}

	void reconstruct()
	{

		int data_dim = (do_3d_rot) ? 3 : 2;
		if (fn_debug != "")
		{
			BackProjector backprojector(debug_ori_size, 3, fn_sym, interpolator, padding_factor, r_min_nn, blob_order, blob_radius, blob_alpha, data_dim, skip_gridding);
			backprojector.initialiseDataAndWeight(debug_size);
			backprojector.data.printShape();
			backprojector.weight.printShape();
			Image<RFLOAT> It;
			It.read(fn_debug+"_data_real.mrc");
			It().setXmippOrigin();
			It().xinit=0;

			It().printShape();
			FOR_ALL_ELEMENTS_IN_ARRAY3D(It())
			{
				A3D_ELEM(backprojector.data, k, i, j).real = A3D_ELEM(It(), k, i, j);
			}
			It.read(fn_debug+"_data_imag.mrc");
			It().setXmippOrigin();
			It().xinit=0;
			FOR_ALL_ELEMENTS_IN_ARRAY3D(It())
			{
				A3D_ELEM(backprojector.data, k, i, j).imag = A3D_ELEM(It(), k, i, j);
			}
			It.read(fn_debug+"_weight.mrc");
			It().setXmippOrigin();
			It().xinit=0;
			FOR_ALL_ELEMENTS_IN_ARRAY3D(It())
			{
				A3D_ELEM(backprojector.weight, k, i, j) = A3D_ELEM(It(), k, i, j);
			}

			MultidimArray<RFLOAT> dummy;
			backprojector.reconstruct(It(), iter, false, 1., dummy, dummy, dummy, dummy, dummy, 1., false, true, nr_threads, -1);
	    	It.write(fn_out);
	    	std::cerr<<" Done writing map in "<<fn_out<<std::endl;
                exit(1);

		}
		else
		{

		RFLOAT rot, tilt, psi, fom;
		Matrix2D<RFLOAT> A3D;
		MultidimArray<Complex > Faux, F2D, Fsub;
		MultidimArray<RFLOAT> Fweight, Fctf, dummy;
		Image<RFLOAT> vol, img, sub;
		FourierTransformer transformer;
		Matrix1D< RFLOAT > trans(2);
		Projector proj;
		int mysize;
//#define DEBUG_WW
#ifdef DEBUG_WW

   		// Get dimension of the images
   		(DF).firstObject();
   		DF.getValue(EMDL_IMAGE_NAME, fn_img);
   		img.read(fn_img);
   		mysize=(int)XSIZE(img());
		BackProjector backprojectort(mysize, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn, blob_order, blob_radius, blob_alpha);
		backprojectort.initZeros(2 * r_max);

		Image<RFLOAT> Imagn, Iphas, Iw, tvol;
		Imagn.read("FEW_it24_rank2_data_magn.spi");
		Iphas.read("FEW_it24_rank2_data_phas.spi");
		Iw.read("FEW_it24_rank2_weight.spi");
        Iw().setXmippOrigin();
        Iw().xinit=0;

		// Write out invw
		Image<RFLOAT> oo;
		oo=Iw;
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Iw())
		{
			if (DIRECT_MULTIDIM_ELEM(Iw(), n) > 1e-2)
				DIRECT_MULTIDIM_ELEM(oo(), n) = 1./ DIRECT_MULTIDIM_ELEM(Iw(), n);
		}
		oo.write("invw.spi");

		Imagn().printShape();
		backprojectort.data.printShape();
		backprojectort.weight.printShape();
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Imagn())
		{
			RFLOAT realval = sin(DIRECT_MULTIDIM_ELEM(Iphas(), n)) * DIRECT_MULTIDIM_ELEM(Imagn(), n);
			RFLOAT imagval = cos(DIRECT_MULTIDIM_ELEM(Iphas(), n)) * DIRECT_MULTIDIM_ELEM(Imagn(), n);
			DIRECT_MULTIDIM_ELEM(backprojectort.data, n) = (Complex)(realval, imagval);
		}
		backprojectort.weight = Iw();
  		std::cerr << "Starting the reconstruction ..." << std::endl;
   		backprojectort.reconstruct(tvol(), iter, false, 1., dummy, dummy, dummy, dummy, dummy, 1., false, false, nr_threads, -1);
    	tvol.write(fn_out);
    	std::cerr<<" Done writing TMPPPPPPPPPPPPPPPPP debugging!!!c map in "<<fn_out<<std::endl;
		exit(0);
#endif


   		// Get dimension of the images
   		if (do_reconstruct_ctf)
   		{
   			mysize = ctf_dim;
   			img().resize(ctf_dim, ctf_dim);
   			img().setXmippOrigin();
   		}
   		else
   		{
			(DF).firstObject();
			DF.getValue(EMDL_IMAGE_NAME, fn_img);
			img.read(fn_img);
			mysize=(int)XSIZE(img());
   		}


   		if (DF.containsLabel(EMDL_CTF_MAGNIFICATION) && DF.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
    	{
    		RFLOAT mag, dstep;
   			DF.getValue(EMDL_CTF_MAGNIFICATION, mag);
   			DF.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
   			angpix = 10000. * dstep / mag;
   			std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
    	}

   		BackProjector backprojector(mysize, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn, blob_order, blob_radius, blob_alpha, data_dim, skip_gridding);
   		// This one is only needed for do_reconstruct_meas, but non-empty constructor does not exist...
   		BackProjector backprojector2(mysize, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn, blob_order, blob_radius, blob_alpha, data_dim, skip_gridding);

   		if (maxres < 0.)
   			r_max = -1;
   		else
   			r_max = CEIL(mysize * angpix / maxres);

   		backprojector.initZeros(2 * r_max);
   		Projector projector(mysize, interpolator, padding_factor, r_min_nn);

   		if (fn_sub != "")
   		{
   			sub.read(fn_sub);
   			projector.computeFourierTransformMap(sub(), dummy, 2 * r_max);
   		}

		// Check for beam-tilt parameters in the input star file
   		if (do_beamtilt || ( DF.containsLabel(EMDL_IMAGE_BEAMTILT_X) || DF.containsLabel(EMDL_IMAGE_BEAMTILT_Y) ) )
   				std::cout << " + Using the beamtilt parameters in the input STAR file" << std::endl;

		std::cerr << "Back-projecting all images ..." << std::endl;
   		int imgno = 0;
		time_config();
   		init_progress_bar(DF.size());
   		FOR_ALL_OBJECTS_IN_METADATA_TABLE(DF)
   		{

   			if (!do_reconstruct_ctf)
   			{
   				DF.getValue(EMDL_IMAGE_NAME, fn_img);
   				img.read(fn_img);
   				img().setXmippOrigin();
   			}

			// Rotations
			if (ref_dim==2)
   			{
   				rot = tilt = 0.;
   			}
   			else
   			{
   				DF.getValue( EMDL_ORIENT_ROT, rot);
   				DF.getValue( EMDL_ORIENT_TILT, tilt);
   			}
  			psi = 0.;
  			DF.getValue( EMDL_ORIENT_PSI, psi);
  			if (angular_error > 0.)
  			{
  				rot += rnd_gaus(0., angular_error);
  				tilt += rnd_gaus(0., angular_error);
  				psi += rnd_gaus(0., angular_error);
  				//std::cout << rnd_gaus(0., angular_error) << std::endl;
  			}

   			Euler_angles2matrix(rot, tilt, psi, A3D);

   			// Translations (either through phase-shifts or in real space
   			trans.initZeros();
   			DF.getValue( EMDL_ORIENT_ORIGIN_X, XX(trans));
   			DF.getValue( EMDL_ORIENT_ORIGIN_Y, YY(trans));
  			if (shift_error > 0.)
   			{
   				XX(trans) += rnd_gaus(0., shift_error);
   				YY(trans) += rnd_gaus(0., shift_error);
   			}
  			if (do_3d_rot)
   			{
   				trans.resize(3);
   				DF.getValue( EMDL_ORIENT_ORIGIN_Z, ZZ(trans));
   	  			if (shift_error > 0.)
   	   				ZZ(trans) += rnd_gaus(0., shift_error);
   			}

   			if (do_fom_weighting)
   				DF.getValue( EMDL_PARTICLE_FOM, fom);


   			// Use either selfTranslate OR shiftImageInFourierTransform!!
   			//selfTranslate(img(), trans, WRAP);
   			CenterFFT(img(), true);
   			transformer.FourierTransform(img(), F2D);
   			if (ABS(XX(trans)) > 0. || ABS(YY(trans)) > 0.)
   			{
   				if (do_3d_rot)
   					shiftImageInFourierTransform(F2D, F2D, XSIZE(img()), XX(trans), YY(trans), ZZ(trans));
   				else
   					shiftImageInFourierTransform(F2D, F2D, XSIZE(img()), XX(trans), YY(trans));
   			}

			Fctf.resize(F2D);
			Fctf.initConstant(1.);
			// Apply CTF if necessary
			if (do_ctf || do_reconstruct_ctf)
			{
				CTF ctf;
				ctf.read(DF, DF);
				ctf.getFftwImage(Fctf, mysize, mysize, angpix, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true);

				if (do_beamtilt || (DF.containsLabel(EMDL_IMAGE_BEAMTILT_X) && DF.containsLabel(EMDL_IMAGE_BEAMTILT_Y) ))
				{
					if (DF.containsLabel(EMDL_IMAGE_BEAMTILT_X))
						DF.getValue(EMDL_IMAGE_BEAMTILT_X, beamtilt_x);
					if (DF.containsLabel(EMDL_IMAGE_BEAMTILT_Y))
						DF.getValue(EMDL_IMAGE_BEAMTILT_Y, beamtilt_y);
					selfApplyBeamTilt(F2D, beamtilt_x, beamtilt_y, ctf.lambda, ctf.Cs, angpix, mysize);
				}

			}

			// Subtract reference projection
			if (fn_sub != "")
			{
				Fsub.resize(F2D);
				projector.get2DFourierTransform(Fsub, A3D, IS_NOT_INV);

				// Apply CTF if necessary
				if (do_ctf)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsub)
					{
						DIRECT_MULTIDIM_ELEM(Fsub, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}

				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fsub)
				{
					DIRECT_MULTIDIM_ELEM(F2D, n) -= DIRECT_MULTIDIM_ELEM(Fsub, n);
				}
				// Back-project difference image
				backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV);
			}
			else
			{
				// "Normal" reconstruction, multiply X by CTF, and W by CTF^2
				if (do_reconstruct_ctf)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(F2D, n)  = DIRECT_MULTIDIM_ELEM(Fctf, n);
						if (do_reconstruct_ctf2)
							DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						DIRECT_MULTIDIM_ELEM(Fctf, n) = 1.;
					}
				}
				else if (do_ctf)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
						DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
					}
				}

				// Do the following after squaring the CTFs!
				if (do_fom_weighting)
				{
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(F2D, n)  *= fom;
						DIRECT_MULTIDIM_ELEM(Fctf, n)  *= fom;
					}
				}

//#define DEBUG_RECONSTRUCT_ONLY
#ifdef DEBUG_RECONSTRUCT_ONLY
					if (fn_img == "/lmb/home/scheres/data/betaGal_rh_withnoise/betaGal_2010_all_p_2x2_unflipped/img00001.win100")
					//if (part_id == my_first_particle_id)
					{
						std::cerr << " fn_img= " << fn_img << std::endl;
						//std::cerr << " myscale= " << myscale << std::endl;
						//std::cerr << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << " normcorr= " << normcorr << std::endl;
						//std::cerr << " sigma2_fudge= " << sigma2_fudge << " mymodel.tau2_fudge_factor= " << mymodel.tau2_fudge_factor<< std::endl;
						//std::cerr << " A3D= " << A3D << std::endl;
						std::cerr << " A3D= " << A3D << std::endl;
						//std::cerr << " exp_R_mic= " << exp_R_mic << std::endl;
						std::cerr << " rot= " << rot << " tilt= " << tilt << " psi= " << psi << " xoff= "<< XX(trans)<< " yoff= "<<YY(trans)<<std::endl;
						//std::cerr << "mic_id= "<<mic_id<<" mymodel.sigma2_noise[mic_id]= " << mymodel.sigma2_noise[mic_id] << std::endl;
						Image<RFLOAT> It;
						It()=Fctf;
						It.write("reconstruct_Fctf.spi");
						It().resize(mysize, mysize);
						MultidimArray<Complex > Faux = F2D;
						FourierTransformer transformer;
						transformer.inverseFourierTransform(Faux, It());
						CenterFFT(It(), false);
						It.write("reconstruct_Mimg.spi");
					}
#endif


				backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);


				if (do_reconstruct_meas)
				{
					// Reconstruct from all-one F2Ds
					FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
					{
						DIRECT_MULTIDIM_ELEM(F2D, n) = 1.;
					}
					backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
				}
			}



   			if (imgno++%60==0)
   				progress_bar(imgno);
		}
   		progress_bar(DF.size());


   		bool do_map = false;
   		bool do_use_fsc = false;
   		MultidimArray<RFLOAT> fsc;
   		fsc.resize(mysize/2+1);
   		if (fn_fsc != "")
   		{
   			do_map = true;
   			do_use_fsc =true;
   			MetaDataTable MDfsc;
   			MDfsc.read(fn_fsc);
   			FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDfsc)
   			{
   				int idx;
   				RFLOAT val;
   				MDfsc.getValue(EMDL_SPECTRAL_IDX, idx);
   				MDfsc.getValue(EMDL_MLMODEL_FSC_HALVES_REF, val);
   				fsc(idx) =  val;
   			}
   		}
   		std::cerr << "Starting the reconstruction ..." << std::endl;
   		backprojector.symmetrise(nr_helical_asu, helical_twist, helical_rise);
   		backprojector.reconstruct(vol(), iter, do_map, 1., dummy, dummy, dummy, dummy, fsc, 1., do_use_fsc, true, nr_threads, -1);

   		if (do_reconstruct_ctf)
   		{

   			F2D.clear();
   			transformer.FourierTransform(vol(), F2D);

   			// CenterOriginFFT: Set the center of the FFT in the FFTW origin
			Matrix1D<RFLOAT> shift(3);
			XX(shift)=-(RFLOAT)(int)(ctf_dim / 2);
			YY(shift)=-(RFLOAT)(int)(ctf_dim / 2);
			ZZ(shift)=-(RFLOAT)(int)(ctf_dim / 2);
			shiftImageInFourierTransform(F2D, F2D, (RFLOAT)ctf_dim, XX(shift), YY(shift), ZZ(shift));
			vol().setXmippOrigin();
   			vol().initZeros();
   			FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
   			{
   				// Take care of kp==dim/2, as XmippOrigin lies just right off center of image...
   				if ( kp > FINISHINGZ(vol()) || ip > FINISHINGY(vol()) || jp > FINISHINGX(vol()))
   					continue;
   				A3D_ELEM(vol(), kp, ip, jp)    = FFTW_ELEM(F2D, kp, ip, jp).real;
   				A3D_ELEM(vol(), -kp, -ip, -jp) = FFTW_ELEM(F2D, kp, ip, jp).real;
   			}
   			vol() *= (RFLOAT)ctf_dim;

   			// Take sqrt(CTF^2)
   			if (do_reconstruct_ctf2)
   			{
   				FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(vol())
				{
   					if (DIRECT_MULTIDIM_ELEM(vol(), n) > 0.)
   						DIRECT_MULTIDIM_ELEM(vol(), n) = sqrt(DIRECT_MULTIDIM_ELEM(vol(), n));
   					else
   						DIRECT_MULTIDIM_ELEM(vol(), n) = 0.;
				}
   			}


   			if (do_reconstruct_meas)
   			{
   				backprojector2.reconstruct(vol(), iter, do_map, 1., dummy, dummy, dummy, dummy, fsc, 1., do_use_fsc, true, nr_threads, -1);
   	   			F2D.clear();
   	   			transformer.FourierTransform(vol(), F2D);

   	   			// CenterOriginFFT: Set the center of the FFT in the FFTW origin
   				Matrix1D<RFLOAT> shift(3);
   				XX(shift)=-(RFLOAT)(int)(ctf_dim / 2);
   				YY(shift)=-(RFLOAT)(int)(ctf_dim / 2);
   				ZZ(shift)=-(RFLOAT)(int)(ctf_dim / 2);
   				shiftImageInFourierTransform(F2D, F2D, (RFLOAT)ctf_dim, XX(shift), YY(shift), ZZ(shift));

   				// Divide each value by the average in each resolution shell.
   				// This is the multiplicative correction that will later be applied to sigma2_noise
   				MultidimArray<RFLOAT> sum, count;
   				sum.initZeros(ctf_dim / 2);
   				count.initZeros(ctf_dim / 2);
   				FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
   	   			{
					int ires = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
					DIRECT_A1D_ELEM(sum, ires) += DIRECT_A3D_ELEM(F2D, k, i, j).real;
					DIRECT_A1D_ELEM(count, ires) += 1.;
   	   			}
   				FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sum)
   				{
   					if (DIRECT_A1D_ELEM(count, i) > 0.)
   						DIRECT_A1D_ELEM(sum, i) /= DIRECT_A1D_ELEM(count, i);
   				}
  				FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
   	   			{
					int ires = ROUND(sqrt(kp*kp + ip*ip + jp*jp));
   					if (DIRECT_A1D_ELEM(count, ires) > 0.)
   						DIRECT_A3D_ELEM(F2D, k, i, j) /=  DIRECT_A1D_ELEM(sum, ires);
   	   			}

   				FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
   	   			{
   	   				// Take care of kp==dim/2, as XmippOrigin lies just right off center of image...
   	   				if ( kp > FINISHINGZ(vol()) || ip > FINISHINGY(vol()) || jp > FINISHINGX(vol()))
   	   					continue;
   	   				if (kp == 0 && ip == 0 && jp == 0)
   	   					continue;
   	   				//A3D_ELEM(vol(), kp, ip, jp)    = FFTW_ELEM(F2D, kp, ip, jp).real;
   	   				//NOW only set one Hermitian half to the number of measured components!
   	   				A3D_ELEM(vol(), -kp, -ip, -jp) = FFTW_ELEM(F2D, kp, ip, jp).real;
   	   			}
   	   			vol() *= (RFLOAT)ctf_dim;
   			}

   		}

   		vol.write(fn_out);
    	std::cerr<<" Done writing map in "<<fn_out<<std::endl;

	}
	}


};


int main(int argc, char *argv[])
{
	reconstruct_parameters prm;

	try
    {

		prm.read(argc, argv);

		prm.reconstruct();

    }
    catch (RelionError XE)
    {
        //prm.usage();
        std::cerr << XE;
        exit(1);
    }
    return 0;
}



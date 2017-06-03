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

#include <src/projector.h>
#include <src/backprojector.h>
#include <src/fftw.h>
#include <src/args.h>
#include <src/ctf.h>
#include <src/strings.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/metadata_table.h>
#include <src/ml_model.h>
#include <src/exp_model.h>
#include <src/healpix_sampling.h>
class project_parameters
{
public:

	FileName fn_map, fn_ang, fn_out, fn_img, fn_model, fn_sym, fn_mask;
	RFLOAT rot, tilt, psi, xoff, yoff, zoff, angpix, maxres, stddev_white_noise, particle_diameter, ana_prob_range, ana_prob_step;
	int padding_factor;
	int r_max, r_min_nn, interpolator;
    bool do_only_one, do_ctf, do_ctf2, ctf_phase_flipped, do_ctf_intact_1st_peak, do_timing, do_add_noise, do_subtract_exp, do_ignore_particle_name, do_3d_rot;
	// I/O Parser
	IOParser parser;
	MlModel model;


	void usage()
	{
		parser.writeUsage(std::cerr);
	}

	void read(int argc, char **argv)
	{
		parser.setCommandLine(argc, argv);

		int general_section = parser.addSection("Options");
		fn_map = parser.getOption("--i", "Input map to be projected");
		fn_out = parser.getOption("--o", "Rootname for output projections", "proj");
       	do_ctf = parser.checkOption("--ctf", "Apply CTF to reference projections");
       	ctf_phase_flipped = parser.checkOption("--ctf_phase_flip", "Flip phases of the CTF in the output projections");
       	do_ctf_intact_1st_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
       	angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "1"));
		fn_mask = parser.getOption("--mask", "Mask that will be applied to the input map prior to making projections", "");
       	fn_ang = parser.getOption("--ang", "STAR file with orientations for multiple projections (if None, assume single projection)","None");
       	rot = textToFloat(parser.getOption("--rot", "First Euler angle (for a single projection)", "0"));
       	tilt = textToFloat(parser.getOption("--tilt", "Second Euler angle (for a single projection)", "0"));
       	psi = textToFloat(parser.getOption("--psi", "Third Euler angle (for a single projection)", "0"));
       	xoff = textToFloat(parser.getOption("--xoff", "Origin X-offsets (in pixels) (for a single projection)", "0"));
       	yoff = textToFloat(parser.getOption("--yoff", "Origin Y-offsets (in pixels) (for a single projection)", "0"));
       	zoff = textToFloat(parser.getOption("--zoff", "Origin Z-offsets (in pixels) (for a single 3D rotation)", "0"));
       	do_add_noise = parser.checkOption("--add_noise", "Add noise to the output projections (only with --ang)");
       	stddev_white_noise = textToFloat(parser.getOption("--white_noise", "Standard deviation of added white Gaussian noise", "0"));
       	fn_model = parser.getOption("--model_noise", "Model STAR file with power spectra for coloured Gaussian noise", "");
       	do_subtract_exp = parser.checkOption("--subtract_exp", "Subtract projections from experimental images (in --ang)");
       	do_ignore_particle_name = parser.checkOption("--ignore_particle_name", "Ignore the rlnParticleName column (in --ang)");
       	do_only_one = (fn_ang == "None");
       	do_3d_rot = parser.checkOption("--3d_rot", "Perform 3D rotations instead of projection into 2D images");

       	maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
       	padding_factor = textToInteger(parser.getOption("--pad", "Padding factor", "2"));
    	do_ctf2 = parser.checkOption("--ctf2", "Apply CTF*CTF to reference projections");
       	if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation"))
       		interpolator = NEAREST_NEIGHBOUR;
       	else
       		interpolator = TRILINEAR;

       	// Hidden
       	r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

       	// Check for errors in the command-line option
    	if (parser.checkForErrors())
    		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

	}

	void project()
	{

		MetaDataTable DFo, MDang;
    	Matrix2D<RFLOAT> A3D;
    	FileName fn_expimg;

    	MultidimArray<Complex > F3D, F2D, Fexpimg;
    	MultidimArray<RFLOAT> Fctf, dummy;
    	Image<RFLOAT> vol, img, expimg;
    	FourierTransformer transformer, transformer_expimg;

		std::cout << " Reading map: " << fn_map << std::endl;
    	vol.read(fn_map);
    	std::cout << " Done reading map!" << std::endl;

    	if (fn_mask != "")
    	{
    		Image<RFLOAT> msk;
    		msk.read(fn_mask);
    		if (!msk().sameShape(vol()))
    			REPORT_ERROR("project ERROR: mask and map have different sizes!");
    		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(vol())
    		DIRECT_MULTIDIM_ELEM(vol(), n) *= DIRECT_MULTIDIM_ELEM(msk(), n);
    	}

    	if (!do_only_one)
    	{
    		std::cout << " Reading STAR file with all angles " << fn_ang << std::endl;
    		MDang.read(fn_ang);
    		std::cout << " Done reading STAR file!" << std::endl;
    	}


    	if (angpix < 0.)
    	{
    		if (do_only_one)
    			REPORT_ERROR("project ERROR: please provide pixel size in Angstroms through --angpix");
    		if (MDang.containsLabel(EMDL_CTF_MAGNIFICATION) && MDang.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
            {
                MDang.goToObject(0);
				RFLOAT mag, dstep;
				MDang.getValue(EMDL_CTF_MAGNIFICATION, mag);
				MDang.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
				angpix = 10000. * dstep / mag;
				std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
            }
            else
            {
            	REPORT_ERROR("project ERROR: please provide pixel size in Angstroms through --angpix");
            }
    	}

    	// Now that we have the size of the volume, check r_max
   		if (maxres < 0.)
   			r_max = XSIZE(vol());
   		else
   			r_max = CEIL(XSIZE(vol()) * angpix / maxres);

    	// Set right size of F2D and initialize to zero
    	if (do_3d_rot)
    		img().resize(ZSIZE(vol()), YSIZE(vol()), XSIZE(vol()));
    	else
    		img().resize(YSIZE(vol()), XSIZE(vol()));
    	transformer.setReal(img());
    	transformer.getFourierAlias(F2D);

    	// Set up the projector
    	int data_dim = (do_3d_rot) ? 3 : 2;
    	Projector projector((int)XSIZE(vol()), interpolator, padding_factor, r_min_nn, data_dim);
    	projector.computeFourierTransformMap(vol(), dummy, 2* r_max);

    	if (do_only_one)
    	{
    		Euler_rotation3DMatrix(rot, tilt, psi, A3D);
    		F2D.initZeros();
    		projector.get2DFourierTransform(F2D, A3D, IS_NOT_INV);
            if (ABS(xoff) > 0.001 || ABS(yoff) > 0.001 || (do_3d_rot && ABS(zoff) > 0.001) )
            {
            	Matrix1D<RFLOAT> shift(2);
            	XX(shift) = -xoff;
            	YY(shift) = -yoff;
            	if (do_3d_rot)
            	{
            		shift.resize(3);
            		ZZ(shift) = -zoff;
            		shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift), ZZ(shift));
            	}
            	else
            		shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift));
            }

            // Feb 01,2017 - Shaoda, add white noise to 2D / 3D single images
            if (do_add_noise)
            {
            	if ( (!(stddev_white_noise > 0.)) || (fn_model != "") )
            		REPORT_ERROR("ERROR: Only add --white_noise to a single image!");
            	// fftw normalization and factor sqrt(2) for two-dimensionality of complex plane
            	// TODO: sqrt(2) ??? Why ???
            	stddev_white_noise /= (data_dim == 3) ? (XSIZE(vol()) * XSIZE(vol())) : (XSIZE(vol()) * sqrt(2));
                // Add white noise
                FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
                {
                    DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0., stddev_white_noise);
                    DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0., stddev_white_noise);
                }
            }

        	transformer.inverseFourierTransform();
        	// Shift the image back to the center...
        	CenterFFT(img(), false);
        	img.write(fn_out);
        	std::cout<<" Done writing "<<fn_out<<std::endl;
    	}
    	else
    	{
            init_progress_bar(MDang.numberOfObjects());
            DFo.clear();
            rot = tilt = psi = xoff = yoff = zoff = 0.;

            // Can only add noise to multiple images
            // Feb 01,2017 - Shaoda, now we can add white noise to 2D / 3D single images
            if (do_add_noise)
            {
            	if (fn_model != "")
            		model.read(fn_model);
            	else if (stddev_white_noise > 0.)
            		stddev_white_noise /= XSIZE(vol()) * sqrt(2); // fftw normalization and factor sqrt(2) for two-dimensionality of complex plane
            	else
            		REPORT_ERROR("ERROR: When adding noise provide either --model_noise or --white_noise");
            }


            long int imgno = 0;
            FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDang)
            {

            	MDang.getValue(EMDL_ORIENT_ROT, rot);
                MDang.getValue(EMDL_ORIENT_TILT, tilt);
                MDang.getValue(EMDL_ORIENT_PSI, psi);
                MDang.getValue(EMDL_ORIENT_ORIGIN_X, xoff);
                MDang.getValue(EMDL_ORIENT_ORIGIN_Y, yoff);
                if (do_3d_rot)
                    MDang.getValue(EMDL_ORIENT_ORIGIN_Z, zoff);

                Euler_rotation3DMatrix(rot, tilt, psi, A3D);
                F2D.initZeros();
                projector.get2DFourierTransform(F2D, A3D, IS_NOT_INV);

                if (ABS(xoff) > 0.001 || ABS(yoff) > 0.001 || (do_3d_rot && ABS(zoff) > 0.001) )
                {
                    Matrix1D<RFLOAT> shift(2);
                    XX(shift) = -xoff;
                    YY(shift) = -yoff;

                    if (do_3d_rot)
                    {
                        shift.resize(3);
                        ZZ(shift) = -zoff;
                        shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift), ZZ(shift) );
                    }
                    else
                    	shiftImageInFourierTransform(F2D, F2D, XSIZE(vol()), XX(shift), YY(shift) );
                }

                // Apply CTF if necessary
                CTF ctf;
                if (do_ctf || do_ctf2)
                {
                	if (do_3d_rot)
                	{
                		Image<RFLOAT> Ictf;
                		FileName fn_ctf;
                		MDang.getValue(EMDL_CTF_IMAGE, fn_ctf);
                		Ictf.read(fn_ctf);
        				Ictf().setXmippOrigin();
        				// Set the CTF-image in Fctf
        				FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf)
        				{
        					// Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
        					DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
        				}
                	}
                	else
                	{
						ctf.read(MDang, MDang);
						Fctf.resize(F2D);
						ctf.getFftwImage(Fctf, XSIZE(vol()), XSIZE(vol()), angpix, ctf_phase_flipped, false,  do_ctf_intact_1st_peak, true);
                	}
                	FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
                    {
                        DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                        if (do_ctf2)
                        	DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                    }
                }

                // Apply Gaussian noise
                if (do_add_noise)
                {
                    if (fn_model !="")
                    {
                        //// 23MAY2014: for preparation of 1.3 release: removed reading a exp_model, replaced by just reading MDang
                        // This does however mean that I no longer know mic_id of this image: replace by 0....
                        FileName fn_group;
                        if (MDang.containsLabel(EMDL_MLMODEL_GROUP_NAME))
                        {
                            MDang.getValue(EMDL_MLMODEL_GROUP_NAME, fn_group);
                        }
                        else
                        {
                            if (MDang.containsLabel(EMDL_MICROGRAPH_NAME))
                            {
                                MDang.getValue(EMDL_MICROGRAPH_NAME, fn_group);
                            }
                            else
                            {
                                REPORT_ERROR("ERROR: cannot find rlnGroupName or rlnMicrographName in the input --ang file...");
                            }
                        }
                        int my_mic_id = -1;
                        for (int mic_id = 0; mic_id < model.group_names.size(); mic_id++)
                        {
                            if (fn_group == model.group_names[mic_id])
                            {
                                my_mic_id = mic_id;
                                break;
                            }
                        }
                        if (my_mic_id < 0)
                            REPORT_ERROR("ERROR: cannot find " + fn_group + " in the input model file...");

                        RFLOAT normcorr = 1.;
                        if (MDang.containsLabel(EMDL_IMAGE_NORM_CORRECTION))
                        {
                            MDang.getValue(EMDL_IMAGE_NORM_CORRECTION, normcorr);
                        }

                        // Add coloured noise
                        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
                        {
                            int ires = ROUND( sqrt( (RFLOAT)(kp*kp + ip*ip + jp*jp) ) );
                            ires = XMIPP_MIN(ires, model.ori_size/2); // at freqs higher than Nyquist: use last sigma2 value

                            RFLOAT sigma = sqrt(DIRECT_A1D_ELEM(model.sigma2_noise[my_mic_id], ires));
                            DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0., sigma);
                            DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0., sigma);
                        }
                    }
                    else
                    {
                        // Add white noise
                        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D)
                        {
                            DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0., stddev_white_noise);
                            DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0., stddev_white_noise);
                        }
                    }
                }

                transformer.inverseFourierTransform();
                // Shift the image back to the center...
                CenterFFT(img(), false);

                // Subtract the projection from the corresponding experimental image
                if (do_subtract_exp)
                {
                    MDang.getValue(EMDL_IMAGE_NAME, fn_expimg);
                    MDang.setValue(EMDL_IMAGE_ORI_NAME, fn_expimg); // Store fn_expimg in rlnOriginalParticleName
                    expimg.read(fn_expimg);
                    img() = expimg() - img();
                }

                if (do_3d_rot)
                {
                    fn_img.compose(fn_out, imgno+1,"mrc");
                    img.write(fn_img);
                }
                else
                {
                	// Write this particle to the stack on disc
                    // First particle: write stack in overwrite mode, from then on just append to it
                    fn_img.compose(imgno+1,fn_out+".mrcs");
                    if (imgno == 0)
                        img.write(fn_img, -1, false, WRITE_OVERWRITE);
                    else
                        img.write(fn_img, -1, false, WRITE_APPEND);
                }

                // Set the image name to the output STAR file
                DFo.addObject();
                DFo.setObject(MDang.getObject());
                DFo.setValue(EMDL_IMAGE_NAME,fn_img);

                if (imgno%60==0) progress_bar(imgno);
                imgno++;
            }
            progress_bar(MDang.numberOfObjects());

            // Write out STAR file with all information
            fn_img = fn_out + ".star";
            DFo.write(fn_img);
            std::cout<<" Done writing "<<imgno<<" images in "<<fn_img<<std::endl;

    	} // end else do_only_one

	}// end project function

};

int main(int argc, char *argv[])
{
	time_config();
	project_parameters prm;

	try
    {
		prm.read(argc, argv);

		prm.project();
    }

    catch (RelionError XE)
    {
        //prm.usage();
        std::cerr << XE;
        exit(1);
    }

    return 0;

}

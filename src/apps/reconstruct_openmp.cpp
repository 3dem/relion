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
#include <omp.h>

#include <src/jaz/vtk_helper.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/image_op.h>

class reconstruct_parameters
{
    public:

    FileName fn_out, fn_sel, fn_img, fn_sym, fn_sub, fn_fsc, fn_debug, image_path;

    MetaDataTable DF;

    int r_max, r_min_nn, blob_order, ref_dim, interpolator, iter,
        nr_omp_threads, debug_ori_size, debug_size,
        ctf_dim, nr_helical_asu, newbox, width_mask_edge, nr_sectors, subset;

    RFLOAT blob_radius, blob_alpha, angular_error, shift_error, angpix, maxres,
        beamtilt_x, beamtilt_y,
        beamtilt_xx, beamtilt_xy, beamtilt_yy,
        helical_rise, helical_twist;

    bool do_ctf, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak,
        do_fom_weighting, do_3d_rot, do_reconstruct_ctf, do_beamtilt, cl_beamtilt, anisoTilt, do_ewald;

    bool skip_gridding, do_reconstruct_ctf2, do_reconstruct_meas, is_positive, read_weights;

    float padding_factor, mask_diameter;

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
        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of open-mp threads to use. Memory footprint is multiplied by this value.", "16"));
        image_path = parser.getOption("--img", "Image path", "");
        subset = textToInteger(parser.getOption("--subset", "Subset of images to consider (0: even; 1: odd; other: all)", "-1"));

        int ctf_section = parser.addSection("CTF options");
        do_ctf = parser.checkOption("--ctf", "Apply CTF correction");
        intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
        ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
        only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");
        beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in the X-direction (in mrad)", "0."));
        beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in the Y-direction (in mrad)", "0."));
        cl_beamtilt = (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.);

        beamtilt_xx = textToFloat(parser.getOption("--beamtilt_xx", "Anisotropic beamtilt, XX-coefficient", "1."));
        beamtilt_xy = textToFloat(parser.getOption("--beamtilt_xy", "Anisotropic beamtilt, XY-coefficient", "0."));
        beamtilt_yy = textToFloat(parser.getOption("--beamtilt_yy", "Anisotropic beamtilt, YY-coefficient", "1."));

        anisoTilt = beamtilt_xx != 1.0 || beamtilt_xy != 0.0 || beamtilt_yy != 1.0;

        read_weights = parser.checkOption("--read_weights", "Read freq. weight files");
        do_ewald = parser.checkOption("--ewald", "Correct for Ewald-sphere curvature (developmental)");
        mask_diameter  = textToFloat(parser.getOption("--mask_diameter", "Diameter (in A) of mask for Ewald-sphere curvature correction", "-1."));
        width_mask_edge = textToInteger(parser.getOption("--width_mask_edge", "Width (in pixels) of the soft edge on the mask", "3"));
        is_positive = !parser.checkOption("--reverse_curvature", "Try curvature the other way around");
        newbox = textToInteger(parser.getOption("--newbox", "Box size of reconstruction after Ewald sphere correction", "-1"));
        nr_sectors = textToInteger(parser.getOption("--sectors", "Number of sectors for Ewald sphere correction", "2"));

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

        if (cl_beamtilt || do_ewald)
            do_ctf = true;

    }

    void applyCTFPandCTFQ(MultidimArray<Complex> &Fin, CTF &ctf, FourierTransformer &transformer,
            MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ)
    {
        //FourierTransformer transformer;
        outP.resize(Fin);
        outQ.resize(Fin);
        float angle_step = 180./nr_sectors;
        for (float angle = 0.; angle < 180.;  angle +=angle_step)
        {
            MultidimArray<Complex> CTFP(Fin), Fapp(Fin);
            MultidimArray<RFLOAT> Iapp(YSIZE(Fin), YSIZE(Fin));
            // Two passes: one for CTFP, one for CTFQ
            for (int ipass = 0; ipass < 2; ipass++)
            {
                bool is_my_positive = (ipass == 1) ? is_positive : !is_positive;

                // Get CTFP and multiply the Fapp with it
                ctf.getCTFPImage(CTFP, YSIZE(Fin), YSIZE(Fin), angpix, is_my_positive, angle);

                Fapp = Fin * CTFP; // element-wise complex multiplication!

                // inverse transform and mask out the particle....
                transformer.inverseFourierTransform(Fapp, Iapp);
                CenterFFT(Iapp, false);

                softMaskOutsideMap(Iapp, ROUND(mask_diameter/(angpix*2.)), (RFLOAT)width_mask_edge);

                // Re-box to a smaller size if necessary....
                if (newbox > 0 && newbox < YSIZE(Fin))
                {
                    Iapp.setXmippOrigin();
                    Iapp.window(FIRST_XMIPP_INDEX(newbox), FIRST_XMIPP_INDEX(newbox),
                                   LAST_XMIPP_INDEX(newbox),  LAST_XMIPP_INDEX(newbox));

                }

                // Back into Fourier-space
                CenterFFT(Iapp, true);
                transformer.FourierTransform(Iapp, Fapp, false); // false means: leave Fapp in the transformer

                // First time round: resize the output arrays
                if (ipass == 0 && fabs(angle) < XMIPP_EQUAL_ACCURACY)
                {
                    outP.resize(Fapp);
                    outQ.resize(Fapp);
                }

                // Now set back the right parts into outP (first pass) or outQ (second pass)
                float anglemin = angle + 90. - (0.5*angle_step);
                float anglemax = angle + 90. + (0.5*angle_step);

                // angles larger than 180
                bool is_reverse = false;
                if (anglemin >= 180.)
                {
                    anglemin -= 180.;
                    anglemax -= 180.;
                    is_reverse = true;
                }
                MultidimArray<Complex> *myCTFPorQ, *myCTFPorQb;
                if (is_reverse)
                {
                    myCTFPorQ  = (ipass == 0) ? &outQ : &outP;
                    myCTFPorQb = (ipass == 0) ? &outP : &outQ;
                }
                else
                {
                    myCTFPorQ  = (ipass == 0) ? &outP : &outQ;
                    myCTFPorQb = (ipass == 0) ? &outQ : &outP;
                }

                // Deal with sectors with the Y-axis in the middle of the sector...
                bool do_wrap_max = false;
                if (anglemin < 180. && anglemax > 180.)
                {
                    anglemax -= 180.;
                    do_wrap_max = true;
                }

                // use radians instead of degrees
                anglemin = DEG2RAD(anglemin);
                anglemax = DEG2RAD(anglemax);
                FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(CTFP)
                {
                    RFLOAT x = (RFLOAT)jp;
                    RFLOAT y = (RFLOAT)ip;
                    RFLOAT myangle = (x*x+y*y > 0) ? acos(y/sqrt(x*x+y*y)) : 0; // dot-product with Y-axis: (0,1)
                    // Only take the relevant sector now...
                    if (do_wrap_max)
                    {
                        if (myangle >= anglemin)
                            DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
                        else if (myangle < anglemax)
                            DIRECT_A2D_ELEM(*myCTFPorQb, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
                    }
                    else
                    {
                        if (myangle >= anglemin && myangle < anglemax)
                            DIRECT_A2D_ELEM(*myCTFPorQ, i, j) = DIRECT_A2D_ELEM(Fapp, i, j);
                    }
                }
            }
        }
    }


    void reconstruct()
    {
        int data_dim = (do_3d_rot) ? 3 : 2;

        MultidimArray<RFLOAT> dummy;
        Matrix1D<RFLOAT> trans(2);
        Image<RFLOAT> vol, img0, sub;
        Projector proj;
        int mysize;

        // Get dimension of the images
        if (do_reconstruct_ctf)
        {
            mysize = ctf_dim;
            img0().resize(ctf_dim, ctf_dim);
            img0().setXmippOrigin();
        }
        else
        {
            (DF).firstObject();
            DF.getValue(EMDL_IMAGE_NAME, fn_img);
            img0.read(fn_img);
            mysize=(int)XSIZE(img0());
            // When doing Ewald-curvature correction: allow reconstructing smaller box than the input images (which should have large boxes!!)
            if (do_ewald && newbox > 0)
                mysize = newbox;
        }


        if (DF.containsLabel(EMDL_CTF_MAGNIFICATION) && DF.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
        {
            RFLOAT mag, dstep;
            DF.getValue(EMDL_CTF_MAGNIFICATION, mag);
            DF.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep);
            angpix = 10000. * dstep / mag;
            std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
        }

        Projector projector(mysize, interpolator, padding_factor, r_min_nn);

        if (maxres < 0.)
            r_max = -1;
        else
            r_max = CEIL(mysize * angpix / maxres);

        if (fn_sub != "")
        {
            sub.read(fn_sub);
            projector.computeFourierTransformMap(sub(), dummy, 2 * r_max);
        }

        // Check for beam-tilt parameters in the input star file
        if (cl_beamtilt)
        {
            std::cout << " + Using the beamtilt parameters from the command line" << std::endl;
            do_beamtilt = true;
        }
        else if ( DF.containsLabel(EMDL_IMAGE_BEAMTILT_X) || DF.containsLabel(EMDL_IMAGE_BEAMTILT_Y) )
        {
            std::cout << " + Using the beamtilt parameters in the input STAR file" << std::endl;
            do_beamtilt = true;
        }
        else
        {
            std::cout << " + Assuming zero beamtilt" << std::endl;
            do_beamtilt = false;
        }

        if (anisoTilt)
        {
            std::cout << " + Assuming anisotropic coma model" << std::endl;
        }
        else
        {
            std::cout << " + Assuming isotropic coma model" << std::endl;
        }

        std::vector<BackProjector> backprojectors(nr_omp_threads);

        for (int i = 0; i < nr_omp_threads; i++)
        {
            backprojectors[i] = BackProjector(
                        mysize, ref_dim, fn_sym, interpolator,
                        padding_factor, r_min_nn, blob_order,
                        blob_radius, blob_alpha, data_dim, skip_gridding);
        }

        std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&DF);
        const long gc = mdts.size();

        std::cout << "Back-projecting all images ..." << std::endl;

        time_config();
        init_progress_bar(gc);

        #pragma omp parallel num_threads(nr_omp_threads)
        {
            int threadnum = omp_get_thread_num();

            BackProjector& backprojector = backprojectors[threadnum];
            backprojector.initZeros(2 * r_max);

            RFLOAT rot, tilt, psi, fom, r_ewald_sphere;
            Matrix2D<RFLOAT> A3D;
            MultidimArray<RFLOAT> Fctf;
            Matrix1D<RFLOAT> trans(2);
            FourierTransformer transformer;
            Image<RFLOAT> img;

            #pragma omp for
            for (int g = 0; g < gc; g++)
            {
                std::vector<Image<RFLOAT> > obsR = StackHelper::loadStack(&mdts[g], image_path);

                const long pc = obsR.size();

                for (int p = 0; p < pc; p++)
                {
                    int randSubset;
                    mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
                    randSubset = randSubset - 1;

                    if (subset >= 0 && randSubset != subset) continue;

                    /*if (!do_reconstruct_ctf)
                    {
                        mdts[g].getValue(EMDL_IMAGE_NAME, fn_img, p);
                        img.read(fn_img);
                        img().setXmippOrigin();
                    }*/

                    // Rotations
                    if (ref_dim == 2)
                    {
                        rot = tilt = 0.;
                    }
                    else
                    {
                        mdts[g].getValue(EMDL_ORIENT_ROT, rot, p);
                        mdts[g].getValue(EMDL_ORIENT_TILT, tilt, p);
                    }

                    psi = 0.;
                    mdts[g].getValue(EMDL_ORIENT_PSI, psi, p);

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
                    mdts[g].getValue( EMDL_ORIENT_ORIGIN_X, XX(trans), p);
                    mdts[g].getValue( EMDL_ORIENT_ORIGIN_Y, YY(trans), p);

                    if (shift_error > 0.)
                    {
                        XX(trans) += rnd_gaus(0., shift_error);
                        YY(trans) += rnd_gaus(0., shift_error);
                    }

                    if (do_3d_rot)
                    {
                        trans.resize(3);
                        mdts[g].getValue( EMDL_ORIENT_ORIGIN_Z, ZZ(trans), p);

                        if (shift_error > 0.)
                        {
                            ZZ(trans) += rnd_gaus(0., shift_error);
                        }
                    }

                    if (do_fom_weighting)
                    {
                        mdts[g].getValue( EMDL_PARTICLE_FOM, fom, p);
                    }

                    // Use either selfTranslate OR shiftImageInFourierTransform!!
                    //selfTranslate(img(), trans, WRAP);

                    MultidimArray<Complex> Fsub, F2D, F2DP, F2DQ;
                    CenterFFT(obsR[p](), true);

                    transformer.FourierTransform(obsR[p](), F2D);

                    if (ABS(XX(trans)) > 0. || ABS(YY(trans)) > 0.)
                    {
                        if (do_3d_rot)
                        {
                            shiftImageInFourierTransform(F2D, F2D,
                                XSIZE(obsR[p]()), XX(trans), YY(trans), ZZ(trans));
                        }
                        else
                        {
                            shiftImageInFourierTransform(F2D, F2D,
                                XSIZE(obsR[p]()), XX(trans), YY(trans));
                        }
                    }

                    Fctf.resize(F2D);
                    Fctf.initConstant(1.);

                    // Apply CTF if necessary
                    if (do_ctf || do_reconstruct_ctf)
                    {
                        CTF ctf;
                        ctf.read(mdts[g], mdts[g], p);
                        ctf.getFftwImage(Fctf, mysize, mysize, angpix,
                                         ctf_phase_flipped, only_flip_phases,
                                         intact_ctf_first_peak, true);

                        if (do_beamtilt)
                        {
                            if (!cl_beamtilt)
                            {
                                if (mdts[g].containsLabel(EMDL_IMAGE_BEAMTILT_X))
                                {
                                    mdts[g].getValue(EMDL_IMAGE_BEAMTILT_X, beamtilt_x, p);
                                }

                                if (mdts[g].containsLabel(EMDL_IMAGE_BEAMTILT_Y))
                                {
                                    mdts[g].getValue(EMDL_IMAGE_BEAMTILT_Y, beamtilt_y, p);
                                }
                            }

                            if (anisoTilt)
                            {
                                selfApplyBeamTilt(
                                    F2D, beamtilt_x, beamtilt_y,
                                    beamtilt_xx, beamtilt_xy, beamtilt_yy,
                                    ctf.lambda, ctf.Cs, angpix, mysize);
                            }
                            else
                            {
                                selfApplyBeamTilt(
                                    F2D, beamtilt_x, beamtilt_y,
                                    ctf.lambda, ctf.Cs, angpix, mysize);
                            }
                        }

                        // Ewald-sphere curvature correction
                        if (do_ewald)
                        {
                            applyCTFPandCTFQ(F2D, ctf, transformer, F2DP, F2DQ);

                            // Also calculate W, store again in Fctf
                            //std::cerr << " temporarily using very large diameter for weight for debugging...." << std::endl;
                            //ctf.applyWeightEwaldSphereCurvature(Fctf, mysize, mysize, angpix, 100000.*mask_diameter);
                            ctf.applyWeightEwaldSphereCurvature(Fctf, mysize, mysize, angpix, mask_diameter);

                            // Also calculate the radius of the Ewald sphere (in pixels)
                            //std::cerr << " temporarily switching off Ewald sphere curvature for debugging...." << std::endl;
                            //r_ewald_sphere = -1.;
                            r_ewald_sphere = mysize * angpix / ctf.lambda;
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
                        else if (do_ewald)
                        {
                            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
                            {
                                DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                            }
                        }
                        // "Normal" reconstruction, multiply X by CTF, and W by CTF^2
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
                                DIRECT_MULTIDIM_ELEM(Fctf, n) *= fom;
                            }
                        }

                        if (read_weights)
                        {
                            std::string name, fullName;

                            mdts[g].getValue(EMDL_IMAGE_NAME, fullName, 0);
                            name = fullName.substr(fullName.find("@")+1);

                            if (image_path != "")
                            {
                                name = image_path + "/" + name.substr(name.find_last_of("/")+1);
                            }

                            std::string wghName = name;
                            wghName = wghName.substr(0, wghName.find_last_of('.')) + "_weight.mrc";

                            Image<RFLOAT> wgh;
                            wgh.read(wghName);

                            if (   Fctf.ndim != wgh().ndim
                                || Fctf.zdim != wgh().zdim
                                || Fctf.ydim != wgh().ydim
                                || Fctf.xdim != wgh().xdim)
                            {
                                REPORT_ERROR(wghName + " and " + name + " are of unequal size.\n");
                            }

                            for (long int n = 0; n < Fctf.ndim; n++)
                            for (long int z = 0; z < Fctf.zdim; z++)
                            for (long int y = 0; y < Fctf.ydim; y++)
                            for (long int x = 0; x < Fctf.xdim; x++)
                            {
                                DIRECT_NZYX_ELEM(Fctf, n, z, y, x)
                                        *= DIRECT_NZYX_ELEM(wgh(), n, z, y, x);
                            }
                        }

                        DIRECT_A2D_ELEM(F2D, 0, 0) = 0.0;

                        if (do_ewald)
                        {
                            backprojector.set2DFourierTransform(F2DP, A3D, IS_NOT_INV, &Fctf, r_ewald_sphere, true);
                            backprojector.set2DFourierTransform(F2DQ, A3D, IS_NOT_INV, &Fctf, r_ewald_sphere, false);
                        }
                        else
                        {
                            backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
                        }
                    }

                    if (threadnum == 0)
                    {
                        progress_bar(g);
                    }
                }
            }
        }

        progress_bar(gc);

        std::cerr << "\nMerging volumes..." << std::endl;

        BackProjector& backprojector = backprojectors[0];

        for (int bpi = 1; bpi < nr_omp_threads; bpi++)
        {
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(backprojector.data)
            {
                DIRECT_MULTIDIM_ELEM(backprojector.data, n)
                        += DIRECT_MULTIDIM_ELEM(backprojectors[bpi].data, n);
            }

            backprojectors[bpi].data.clear();

            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(backprojector.weight)
            {
                DIRECT_MULTIDIM_ELEM(backprojector.weight, n)
                        += DIRECT_MULTIDIM_ELEM(backprojectors[bpi].weight, n);
            }

            backprojectors[bpi].weight.clear();
        }

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
        backprojector.symmetrise(nr_helical_asu, helical_twist, helical_rise/angpix);
        backprojector.reconstruct(vol(), iter, do_map, 1., dummy, dummy, dummy, dummy,
                                  fsc, 1., do_use_fsc, true, nr_omp_threads, -1, false);

        MultidimArray<Complex> F2D;

        if (do_reconstruct_ctf)
        {
            FourierTransformer transformer;

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
        }

        vol.write(fn_out);
        std::cerr<<" Done writing map in "<<fn_out<<std::endl;
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



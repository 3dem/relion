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

    FileName fn_out, fn_sel, fn_sym, fn_sub, fn_fsc, fn_debug, fn_vol;
    std::string image_path;

    MetaDataTable DF0;

    int r_max, r_min_nn, blob_order, ref_dim, interpolator, iter, nr_fft_threads,
    nr_omp_threads, debug_ori_size, debug_size, ctf_dim, nr_helical_asu, subset;

    RFLOAT blob_radius, blob_alpha, angular_error, shift_error, angpix, maxres,
    beamtilt_x, beamtilt_y, helical_rise, helical_twist;

    bool do_ctf, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak,
    do_fom_weighting, do_3d_rot, do_reconstruct_ctf, do_beamtilt, skip_gridding,
    read_weights, ctf_nyq_wgh;

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

        parser.addSection("General options");

        fn_sel = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
        fn_out = parser.getOption("--o", "Name for output reconstruction","relion.mrc");
        fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
        fn_vol = parser.getOption("--ov", "Name for output volume prior to gridding","");
       	angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "1"));
       	maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
       	padding_factor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
        nr_fft_threads = textToInteger(parser.getOption("--jfft", "Number of threads to use for FFTs", "1"));
        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of open-mp threads to use. Memory footprint is multiplied by this value.", "16"));
        image_path = parser.getOption("--ipath", "Image path", "");
        subset = textToInteger(parser.getOption("--subset", "Subset of images to consider (0: even; 1: odd; other: all)", "-1"));

        if (subset > 1)
        {
            subset = -1;
        }

        parser.addSection("CTF options");

       	do_ctf = parser.checkOption("--ctf", "Apply CTF correction");
    	intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
    	ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
    	only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");
    	beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in the X-direction (in mrad)", "0."));
    	beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in the Y-direction (in mrad)", "0."));
    	do_beamtilt = (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.);
        read_weights = parser.checkOption("--read_weights", "Read freq. weight files");
        ctf_nyq_wgh = parser.checkOption("--ctf_nyq", "Downweight pixels where the CTF-freq. apporaches nyquist");

        parser.addSection("Helical options");

    	nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
    	helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
    	helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

        parser.addSection("Expert options");

       	fn_sub = parser.getOption("--subtract","Subtract projections of this map from the images used for reconstruction", "");

        if (parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation before gridding correction"))
        {
            interpolator = NEAREST_NEIGHBOUR;
        }
       	else
        {
            interpolator = TRILINEAR;
        }

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
    	skip_gridding = parser.checkOption("--skip_gridding", "Skip gridding part of the reconstruction");
    	do_reconstruct_ctf = (ctf_dim > 0);

    	if (do_reconstruct_ctf)
        {
            do_ctf = false;
        }

    	// Hidden
       	r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

    	// Check for errors in the command-line option
    	if (parser.checkForErrors())
        {
            REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
        }

    	// Read MetaData file, which should have the image names and their angles!
    	if (fn_debug == "")
        {
            DF0.read(fn_sel);
        }

     	randomize_random_generator();

     	if (do_beamtilt && ! do_ctf)
        {
            REPORT_ERROR("ERROR: one can only correct for beamtilt in combination with CTF correction!");
        }

    }

    void reconstruct()
    {

        int data_dim = (do_3d_rot) ? 3 : 2;

        Image<RFLOAT> vol, img0, sub;
        int mysize;
        //#define DEBUG_WW

        // Get dimension of the images
        if (do_reconstruct_ctf)
        {
            mysize = ctf_dim;
            img0().resize(ctf_dim, ctf_dim);
            img0().setXmippOrigin();
        }
        else
        {
            FileName fn_img;
            DF0.getValue(EMDL_IMAGE_NAME, fn_img, 0);
            img0.read(fn_img);
            mysize=(int)XSIZE(img0());
        }

        if (DF0.containsLabel(EMDL_CTF_MAGNIFICATION) && DF0.containsLabel(EMDL_CTF_DETECTOR_PIXEL_SIZE))
        {
            RFLOAT mag, dstep;
            DF0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
            DF0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
            angpix = 10000. * dstep / mag;
            std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
        }

        Projector projector0(mysize, interpolator, padding_factor, r_min_nn);
        MultidimArray<RFLOAT> dummy;

        if (maxres < 0.)
            r_max = -1;
        else
            r_max = CEIL(mysize * angpix / maxres);

        if (fn_sub != "")
        {
            sub.read(fn_sub);
            projector0.computeFourierTransformMap(sub(), dummy, 2 * r_max);
        }

        std::vector<BackProjector> backprojectors(nr_omp_threads);

        backprojectors[0] = BackProjector(
                mysize, ref_dim, fn_sym, interpolator,
                padding_factor, r_min_nn, blob_order,
                blob_radius, blob_alpha, data_dim, skip_gridding);

        for (int i = 1; i < nr_omp_threads; i++)
        {
            backprojectors[i] = backprojectors[0];
        }

        std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&DF0);
        const long gc = mdts.size();

        std::cerr << "Back-projecting all images ..." << std::endl;

        time_config();
        init_progress_bar(gc/nr_omp_threads);

        #pragma omp parallel num_threads(nr_omp_threads)
        {
            int threadnum = omp_get_thread_num();

            BackProjector& backprojector = backprojectors[threadnum];
            backprojector.initZeros(2 * r_max);

            RFLOAT rot, tilt, psi, fom;
            Matrix2D<RFLOAT> A3D;
            MultidimArray<Complex> Fsub;
            MultidimArray<RFLOAT> Fctf;
            Matrix1D<RFLOAT> trans(2);

            Projector projector(mysize, interpolator, padding_factor, r_min_nn);

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
                    mdts[g].getValue( EMDL_ORIENT_PSI, psi, p);

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
                    mdts[g].getValue(EMDL_ORIENT_ORIGIN_X, XX(trans), p);
                    mdts[g].getValue(EMDL_ORIENT_ORIGIN_Y, YY(trans), p);

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
                    MultidimArray<Complex> F2D;
                    CenterFFT(obsR[p](), true);

                    FourierTransformer transformer;
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

                        if (do_beamtilt || (   mdts[g].containsLabel(EMDL_IMAGE_BEAMTILT_X)
                                            && mdts[g].containsLabel(EMDL_IMAGE_BEAMTILT_Y) ))
                        {
                            if (mdts[g].containsLabel(EMDL_IMAGE_BEAMTILT_X))
                            {
                                mdts[g].getValue(EMDL_IMAGE_BEAMTILT_X, beamtilt_x, p);
                            }

                            if (mdts[g].containsLabel(EMDL_IMAGE_BEAMTILT_Y))
                            {
                                mdts[g].getValue(EMDL_IMAGE_BEAMTILT_Y, beamtilt_y, p);
                            }

                            selfApplyBeamTilt(F2D, beamtilt_x, beamtilt_y,
                                              ctf.lambda, ctf.Cs, angpix, mysize);
                        }

                        if (ctf_nyq_wgh)
                        {
                            RFLOAT as = (RFLOAT)Fctf.ydim * angpix;

                            for (long int y = 0; y < Fctf.ydim; y++)
                            for (long int x = 0; x < Fctf.xdim; x++)
                            {
                                RFLOAT cf = ctf.getCtfFreq(x/as, y/as) / (as * PI);

                                if (cf > 0.5)
                                {
                                    RFLOAT t = 2.0 * (1.0 - cf);
                                    if (t < 0.0) t = 0.0;

                                    DIRECT_NZYX_ELEM(Fctf, 0, 0, y, x) *= t;
                                    DIRECT_NZYX_ELEM(F2D, 0, 0, y, x) *= t;
                                }
                            }
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

                        backprojector.set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
                    }

                    if (threadnum == 0)
                    {
                        progress_bar(g);
                    }
                }
            }
        }

        //progress_bar(DF0.size());

        std::cerr << "\nMerging volumes..." << std::endl;
        MultidimArray<Complex > F2D0;

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

        backprojector.symmetrise(nr_helical_asu, helical_twist, helical_rise);

        if (fn_vol != "")
        {
            VtkHelper::writeVTK_Complex(backprojector.data, fn_vol+"_data.vtk");
            VtkHelper::writeVTK(backprojector.weight, fn_vol+"_weight.vtk");

            ComplexIO::write(backprojector.data, fn_vol+"_data", ".mrc");
            Image<RFLOAT> wt;
            wt.data = backprojector.weight;
            wt.write(fn_vol+"_weight.mrc");
        }

        std::cerr << "Starting the reconstruction ..." << std::endl;
        backprojector.reconstruct(vol(), iter, do_map, 1.,
            dummy, dummy, dummy, fsc, 1., do_use_fsc, true, nr_fft_threads, -1);

        if (do_reconstruct_ctf)
        {

            F2D0.clear();
            FourierTransformer transformer;
            transformer.FourierTransform(vol(), F2D0);

            // CenterOriginFFT: Set the center of the FFT in the FFTW origin
            Matrix1D<RFLOAT> shift(3);
            XX(shift)=-(RFLOAT)(int)(ctf_dim / 2);
            YY(shift)=-(RFLOAT)(int)(ctf_dim / 2);
            ZZ(shift)=-(RFLOAT)(int)(ctf_dim / 2);
            shiftImageInFourierTransform(F2D0, F2D0, (RFLOAT)ctf_dim,
                                         XX(shift), YY(shift), ZZ(shift));
            vol().setXmippOrigin();
            vol().initZeros();
            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D0)
            {
                // Take care of kp==dim/2, as XmippOrigin lies just right off center of image...
                if ( kp > FINISHINGZ(vol()) || ip > FINISHINGY(vol()) || jp > FINISHINGX(vol()))
                    continue;
                A3D_ELEM(vol(), kp, ip, jp)    = FFTW_ELEM(F2D0, kp, ip, jp).real;
                A3D_ELEM(vol(), -kp, -ip, -jp) = FFTW_ELEM(F2D0, kp, ip, jp).real;
            }
            vol() *= (RFLOAT)ctf_dim;
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
        prm.usage();
        std::cerr << XE;
        exit(1);
    }
    return 0;
}



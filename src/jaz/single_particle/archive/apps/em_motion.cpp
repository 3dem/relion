
#include <unistd.h>
#include <string.h>
#include <fstream>

#include <src/args.h>
#include <src/image.h>
#include <src/fftw.h>
#include <src/complex.h>
#include <src/metadata_table.h>
#include <src/backprojector.h>
#include <src/euler.h>
#include <src/jaz/single_particle/image_log.h>
#include <src/jaz/single_particle/slice_helper.h>
#include <src/jaz/single_particle/spectral_helper.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/backprojection_helper.h>
#include <src/jaz/single_particle/volume_converter.h>
#include <src/jaz/single_particle/complex_io.h>
#include <src/jaz/single_particle/fftw_helper.h>
#include <src/jaz/single_particle/resampling_helper.h>
#include <src/jaz/single_particle/ctf_helper.h>
#include <src/jaz/single_particle/defocus_refinement.h>
#include <src/jaz/single_particle/magnification_refinement.h>
#include <src/jaz/single_particle/refinement_helper.h>
#include <src/jaz/single_particle/stack_helper.h>
#include <src/jaz/single_particle/tilt_refinement.h>
#include <src/jaz/single_particle/motion/motion_refinement.h>
#include <src/jaz/single_particle/img_proc/image_op.h>
#include <src/jaz/single_particle/Fourier_helper.h>
#include <src/jaz/single_particle/fsc_helper.h>
#include <src/jaz/single_particle/damage_helper.h>
#include <src/jaz/single_particle/interpolation.h>
#include <src/jaz/single_particle/distribution_helper.h>
#include <src/jaz/single_particle/noise_helper.h>
#include <src/jaz/single_particle/convolution_helper.h>
#include <src/jaz/single_particle/motion_em.h>
#include <src/jaz/single_particle/local_motion_fit.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/single_particle/parallel_ft.h>

#include <omp.h>

using namespace gravis;

int main(int argc, char *argv[])
{
    std::string starFn, reconFn0, reconFn1, maskFn, moviePath, outPath, fn_sym;

    RFLOAT beamtilt_x, beamtilt_y, paddingFactor, dmga, dmgb, dmgc, totalDose,
            sig_pos, sig_cutoff, k_cutoff,
            sig_vel0, sig_vel1, sig_vel2, sig_vel,
            sig_div0, sig_div1, sig_div2, sig_div,
            sig_acc,
            helical_rise, helical_twist;

    bool applyTilt, evaluate;

    bool rigid = true;
    bool rounded = false;
    bool localOpt = true;
    bool rigidLocalOpt = true;
    bool rigidLocalOptOnly = true;

    long maxMG = -1, minMG = 0;
    int it_number, nr_omp_threads, nr_helical_asu, evalFrames;

    Image<RFLOAT> map0, map1, dummy;

    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        starFn = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
        reconFn0 = parser.getOption("--m0", "Reference, half 1", "");
        reconFn1 = parser.getOption("--m1", "Reference, half 2", "");
        maskFn = parser.getOption("--mask", "Reference mask", "");

        moviePath = parser.getOption("--movies", "Input path to movie files", "");

        outPath = parser.getOption("--out", "Output path", "tracks");

        it_number = textToInteger(parser.getOption("--iters", "Number of iterations", "5"));

        paddingFactor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));

        beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in the X-direction (in mrad)", "0."));
        beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in the Y-direction (in mrad)", "0."));
        applyTilt = (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.);

        dmga = textToFloat(parser.getOption("--dmg_a", "Damage model, parameter a", " 3.40406"));
        dmgb = textToFloat(parser.getOption("--dmg_b", "                        b", "-1.06027"));
        dmgc = textToFloat(parser.getOption("--dmg_c", "                        c", "-0.540896"));

        totalDose = textToFloat(parser.getOption("--dose", "Total electron dose (in e^-/A^2)", "1"));

        sig_pos = textToFloat(parser.getOption("--s_pos", "Position sigma", "30.0"));

        sig_vel0 = textToFloat(parser.getOption("--s_vel_0", "Velocity sigma, frame 1", "6.0"));
        sig_vel1 = textToFloat(parser.getOption("--s_vel_1", "Velocity sigma, frame 2", "4.0"));
        sig_vel2 = textToFloat(parser.getOption("--s_vel_2", "Velocity sigma, frame 3", "3.0"));
        sig_vel = textToFloat(parser.getOption("--s_vel", "Velocity sigma, other frames", "2.0"));

        sig_div0 = textToFloat(parser.getOption("--s_div_0", "Divergence sigma, frame 1", "0.15"));
        sig_div1 = textToFloat(parser.getOption("--s_div_1", "Divergence sigma, frame 2", "0.1"));
        sig_div2 = textToFloat(parser.getOption("--s_div_2", "Divergence sigma, frame 3", "0.05"));
        sig_div = textToFloat(parser.getOption("--s_div", "Divergence sigma, other frames", "0.01"));

        sig_acc = textToFloat(parser.getOption("--s_acc", "Acceleration sigma", "-1.0"));

        sig_cutoff = textToFloat(parser.getOption("--s_cut", "Crop range (in sigma)", "3.0"));
        k_cutoff = textToFloat(parser.getOption("--k_cut", "Freq. cutoff (in pixels)", "-1.0"));

        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of OMP threads", "1"));

        maxMG = textToInteger(getParameter(argc, argv, "--max_MG", "-1"));
        minMG = textToInteger(getParameter(argc, argv, "--min_MG", "0"));

        rigid = parser.checkOption("--rigid", "Rigid alignment instead of EM algorithm");
        localOpt = parser.checkOption("--local", "Refine tracks locally");
        rigidLocalOpt = parser.checkOption("--rigid-local", "Refine rigid track locally");
        rigidLocalOptOnly = parser.checkOption("--rigid-local-only", "Refine only the rigid track locally");

        evalFrames = textToInteger(parser.getOption("--eval", "Measure FSC for this many initial frames", "0"));
        fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
        nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
        helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
        helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

        evaluate = evalFrames > 0;

        try
        {
            map0.read(reconFn0);
        }
        catch (RelionError XE)
        {
            std::cout << "Unable to read map: " << reconFn0 << "\n";
            exit(1);
        }

        try
        {
            map1.read(reconFn1);
        }
        catch (RelionError XE)
        {
            std::cout << "Unable to read map: " << reconFn1 << "\n";
            exit(1);
        }
    }
    catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        exit(1);
    }

    if (map0.data.xdim != map0.data.ydim || map0.data.ydim != map0.data.zdim)
    {
        REPORT_ERROR(reconFn0 + " is not cubical.\n");
    }

    if (map1.data.xdim != map1.data.ydim || map1.data.ydim != map1.data.zdim)
    {
        REPORT_ERROR(reconFn1 + " is not cubical.\n");
    }

    if (   map0.data.xdim != map1.data.xdim
        || map0.data.ydim != map1.data.ydim
        || map0.data.zdim != map1.data.zdim)
    {
        REPORT_ERROR(reconFn0 + " and " + reconFn1 + " are of unequal size.\n");
    }

    const int s = map0.data.xdim;
    const int sh = map0.data.xdim/2 + 1;

    if (maskFn != "")
    {
        std::cout << "masking references...\n";
        Image<RFLOAT> mask, maskedRef;

        try
        {
            mask.read(maskFn);
        }
        catch (RelionError XE)
        {
            std::cout << "Unable to read mask: " << maskFn << "\n";
            exit(1);
        }

        mask.read(maskFn);

        ImageOp::multiply(mask, map0, maskedRef);
        map0 = maskedRef;

        ImageOp::multiply(mask, map1, maskedRef);
        map1 = maskedRef;
    }

    std::vector<double> sigma2ref(sh, 0.0);

    const bool refVar = false;

    if (refVar)
    {
        std::cout << "measuring reference variance...\n";
        Image<Complex> ref0, ref1, deltaRef(sh,s,s);

        FourierTransformer ft;
        ft.FourierTransform(map0(), ref0());
        ft.FourierTransform(map1(), ref1());

        ImageOp::linearCombination(ref0, ref1, 0.5, -0.5, deltaRef);
        sigma2ref = FscHelper::powerSpectrum3D(deltaRef);

        std::ofstream sigR_out(outPath + "_sig2_ref.dat");

        for (int i = 0; i < sigma2ref.size(); i++)
        {
            sigR_out << i << " " << sigma2ref[i] << "\n";
        }
    }

    MetaDataTable mdt0;
    mdt0.read(starFn);

    RFLOAT Cs, lambda, kV;

    mdt0.getValue(EMDL_CTF_CS, Cs, 0);
    mdt0.getValue(EMDL_CTF_VOLTAGE, kV, 0);

    RFLOAT V = kV * 1e3;
    lambda = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));

    std::cout << "transforming references...\n";

    Projector projector0(s, TRILINEAR, paddingFactor, 10, 2);
    projector0.computeFourierTransformMap(map0.data, dummy.data, s);

    Projector projector1(s, TRILINEAR, paddingFactor, 10, 2);
    projector1.computeFourierTransformMap(map1.data, dummy.data, s);

    std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&mdt0);

    RFLOAT mag, dstep;
    mdts[0].getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
    mdts[0].getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
    RFLOAT angpix = 10000 * dstep / mag;

    ObservationModel obsModel(angpix);

    if (applyTilt)
    {
        obsModel = ObservationModel(angpix, Cs, kV * 1e3, beamtilt_x, beamtilt_y);
    }

    const long gc = maxMG >= 0? maxMG : mdts.size()-1;
    const long g0 = minMG;


    std::string name, fullName, movieName;
    mdts[0].getValue(EMDL_IMAGE_NAME, fullName, 0);
    mdts[0].getValue(EMDL_MICROGRAPH_NAME, movieName, 0);
    name = fullName.substr(fullName.find("@")+1);

    std::string finName;

    if (moviePath == "")
    {
        finName = name;
    }
    else
    {
        finName = moviePath + "/" + movieName.substr(movieName.find_last_of("/")+1);
    }

    Image<RFLOAT> stack0;
    stack0.read(finName, false);

    const int pc0 = mdts[0].numberOfObjects();
    const bool zstack = stack0.data.zdim > 1;
    const int stackSize = zstack? stack0.data.zdim : stack0.data.ndim;
    const int fc = stackSize / pc0;


    std::vector<double> sig_vel_vec(4);
    sig_vel_vec[0] = sig_vel0;
    sig_vel_vec[1] = sig_vel1;
    sig_vel_vec[2] = sig_vel2;
    sig_vel_vec[3] = sig_vel;

    std::vector<double> sig_div_vec(4);
    sig_div_vec[0] = sig_div0;
    sig_div_vec[1] = sig_div1;
    sig_div_vec[2] = sig_div2;
    sig_div_vec[3] = sig_div;

    std::cout << "pc0 = " << pc0 << "\n";
    std::cout << "fc = " << fc << "\n";

    std::vector<Image<RFLOAT> > dmgWeight = DamageHelper::damageWeights(
        s, angpix, 0, fc, totalDose, dmga, dmgb, dmgc);


    int k_out = k_cutoff + 21;

    for (int f = 0; f < fc; f++)
    {
        dmgWeight[f].data.xinit = 0;
        dmgWeight[f].data.yinit = 0;

        if (k_cutoff > 0.0)
        {
            std::stringstream stsf;
            stsf << f;
            dmgWeight[f] = FilterHelper::ButterworthEnvFreq2D(dmgWeight[f], k_cutoff-1, k_cutoff+1);
        }
    }

    std::cout << "mg range: " << g0 << ".." << gc << "\n";

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    double t0 = omp_get_wtime();

    int pctot = 0;

    const bool writeDebugImages = false;
    const bool measureFCC = true;

    std::vector<Image<RFLOAT> > tables(nr_omp_threads), tablesV(nr_omp_threads), tablesVW(nr_omp_threads);
    std::vector<Image<RFLOAT> > weights0(nr_omp_threads), weights0V(nr_omp_threads), weights0VW(nr_omp_threads);
    std::vector<Image<RFLOAT> > weights1(nr_omp_threads), weights1V(nr_omp_threads), weights1VW(nr_omp_threads);

    if (evaluate)
    {
        if (measureFCC)
        {
            for (int i = 0; i < nr_omp_threads; i++)
            {
                FscHelper::initFscTable(sh, fc, tables[i], weights0[i], weights1[i]);
                FscHelper::initFscTable(sh, fc, tablesV[i], weights0V[i], weights1V[i]);
                FscHelper::initFscTable(sh, fc, tablesVW[i], weights0VW[i], weights1VW[i]);
            }
        }
    }

    for (long g = g0; g <= gc; g++)
    {
        std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

        std::stringstream stsg;
        stsg << g;

        const int pc = mdts[g].numberOfObjects();
        pctot += pc;

        std::vector<std::vector<Image<Complex> > > movie;

        try
        {
            movie = StackHelper::loadMovieStackFS(
                    &mdts[g], moviePath, false, nr_omp_threads, &fts);
        }
        catch (RelionError XE)
        {
            continue;
        }

        std::vector<double> sigma2 = StackHelper::powerSpectrum(movie);

        if (refVar)
        {
            std::ofstream sigD_out(outPath + "_mg" + stsg.str() + "_sig2_data.dat");

            for (int i = 0; i < sigma2.size(); i++)
            {
                sigD_out << i << " " << sigma2[i] << "\n";

                sigma2[i] += sigma2ref[i];
            }
        }

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        for (int f = 0; f < fc; f++)
        {
            MotionRefinement::noiseNormalize(movie[p][f], sigma2, movie[p][f]);
        }

        std::vector<gravis::d2Vector> positions(pc);
        std::vector<double> defoci(pc);

        for (int p = 0; p < pc; p++)
        {
            mdts[g].getValue(EMDL_IMAGE_COORD_X, positions[p].x, p);
            mdts[g].getValue(EMDL_IMAGE_COORD_Y, positions[p].y, p);

            double du, dv;
            mdts[g].getValue(EMDL_CTF_DEFOCUSU, du, p);
            mdts[g].getValue(EMDL_CTF_DEFOCUSV, dv, p);

            defoci[p] = 0.5*(du + dv)/angpix;
        }

        bool allDefociEqual = true;

        for (int p = 1; p < pc; p++)
        {
            if (defoci[p] != defoci[0])
            {
                allDefociEqual = false;
            }
        }

        if (allDefociEqual)
        {
            std::cout << "WARNING: The defoci are identical for all particles!\n";
            std::cout << "         You may want to determine per-particle defoci first.\n";
        }

        MotionEM motEM(
            projector0, projector1, obsModel, mdts[g], movie, positions,
            sigma2, dmgWeight, sig_pos, sig_vel_vec, sig_div_vec, sig_cutoff, nr_omp_threads);

        if (it_number > 0)
        {
            std::cout << "    computing initial correlations...\n";
            motEM.computeInitial();
        }

        if (writeDebugImages)
        {
            for (int p = 0; p < pc; p++)
            {
                std::stringstream stsp;
                stsp << p;

                ImageLog::write(motEM.posProb[p], "ppDebug/i"+stsg.str()+"_p"+stsp.str()+"_pos_initial", CenterXY);
            }
        }

        if (!rigid)
        for (int it = 0; it < it_number; it++)
        {
            std::stringstream stsit;
            stsit << it;

            std::cout << "    iteration " << it << "\n";
            motEM.iterate();

            if (writeDebugImages)
            {
                for (int p = 0; p < pc; p++)
                {
                    std::stringstream stsp;
                    stsp << p;

                    ImageLog::write(motEM.velProb[p], "ppDebug/i"+stsg.str()+"_p"+stsp.str()+"_vel_"+stsit.str(), CenterXY);
                    ImageLog::write(motEM.posProb[p], "ppDebug/i"+stsg.str()+"_p"+stsp.str()+"_pos_"+stsit.str(), CenterXY);
                }
            }
        }

        std::vector<std::vector<gravis::d2Vector>> tracks(pc);

        if (rigid && it_number > 0)
        {
            std::vector<gravis::d2Vector> globTrack = motEM.getGlobalTrack();

            if (rounded)
            {
                for (int f = 0; f < fc; f++)
                {
                    globTrack[f].x = std::round(globTrack[f].x);
                    globTrack[f].y = std::round(globTrack[f].y);
                }
            }

            for (int p = 0; p < pc; p++)
            {
                tracks[p] = globTrack;
            }
        }
        else
        {
            for (int p = 0; p < pc; p++)
            {
                tracks[p] = motEM.getTrack(p);
            }
        }

        if (localOpt)
        {
            if (rigidLocalOpt)
            {
                std::vector<double> velWgh(fc-1);
                std::vector<double> accWgh(fc-1, sig_acc > 0.0? 1.0/(sig_acc*sig_acc) : 0.0);

                for (int f = 0; f < fc-1; f++)
                {
                    double sv;

                    if (f < sig_vel_vec.size())
                    {
                        sv = sig_vel_vec[f];
                    }
                    else
                    {
                        sv = sig_vel_vec[sig_vel_vec.size()-1];
                    }

                    velWgh[f] = 1.0 / (sv*sv);
                }

                std::vector<std::vector<std::vector<double>>> divWgh(0);

                std::vector<gravis::d2Vector> globTrack = motEM.getGlobalTrack();
                std::vector<std::vector<Image<RFLOAT>>> ccSum(1);
                ccSum[0] = motEM.e_sum;

                for (int f = 0; f < fc-1; f++)
                {
                    velWgh[f] *= pc;
                }

                LocalMotionFit lmf(ccSum, velWgh, accWgh, divWgh,
                    std::vector<d2Vector>(fc, d2Vector(0,0)), nr_omp_threads);

                std::vector<double> initial(2*fc);

                for (int f = 0; f < fc; f++)
                {
                    initial[2*f]     = globTrack[f].x;
                    initial[2*f + 1] = globTrack[f].y;
                }

                std::vector<double> grad0(2*fc*pc);

                lmf.grad(initial, grad0, 0);

                double gl = 0.0;

                for (int i = 0; i < grad0.size(); i++)
                {
                    double gi = grad0[i];
                    gl += gi*gi;
                }

                gl = sqrt(gl);

                std::cout << "gl = " << gl << "\n";

                std::cout << "    optimizing rigid path locally...\n";

                std::vector<double> optPos = GradientDescent::optimize(
                    initial, lmf, 0.05/gl, 1e-9/gl, 1e-9, 10000, 0.0, true);

                for (int p = 0; p < pc; p++)
                {
                    for (int f = 0; f < fc; f++)
                    {
                        tracks[p][f].x = optPos[2*f];
                        tracks[p][f].y = optPos[2*f + 1];
                    }
                }
            }

            if (!rigidLocalOptOnly)
            {
                std::vector<double> velWgh(fc-1);
                std::vector<double> accWgh(fc-1, sig_acc > 0.0? 0.5/(sig_acc*sig_acc) : 0.0);

                for (int f = 0; f < fc-1; f++)
                {
                    double sv;

                    if (f < sig_vel_vec.size())
                    {
                        sv = sig_vel_vec[f];
                    }
                    else
                    {
                        sv = sig_vel_vec[sig_vel_vec.size()-1];
                    }

                    velWgh[f] = 0.5 / (sv*sv);
                }

                std::vector<std::vector<std::vector<double>>> divWgh(fc-1);

                for (int f = 0; f < fc-1; f++)
                {
                    divWgh[f] = std::vector<std::vector<double>>(pc);

                    for (int p = 0; p < pc; p++)
                    {
                        divWgh[f][p] = std::vector<double>(pc);

                        for (int q = 0; q < pc; q++)
                        {
                            d2Vector dp = positions[p] - positions[q];
                            double dd = defoci[p] - defoci[q];

                            double dist = sqrt(dp.x*dp.x + dp.y*dp.y + dd*dd);

                            double sd;

                            if (f < sig_div_vec.size())
                            {
                                sd = sig_div_vec[f];
                            }
                            else
                            {
                                sd = sig_div_vec[sig_div_vec.size()-1];
                            }

                            divWgh[f][p][q] = 0.5 / (sd * sd * dist);
                        }
                    }
                }

                LocalMotionFit lmf(motEM.initialCC, velWgh, accWgh, divWgh,
                            std::vector<d2Vector>(fc, d2Vector(0,0)), nr_omp_threads);

                std::vector<double> initial(2*pc*fc);

                for (int p = 0; p < pc; p++)
                {
                    for (int f = 0; f < fc; f++)
                    {
                        initial[2*(p*fc + f)]     = tracks[p][f].x;
                        initial[2*(p*fc + f) + 1] = tracks[p][f].y;
                    }
                }

                std::vector<double> grad0(2*fc*pc);

                lmf.grad(initial, grad0, 0);

                double gl = 0.0;

                for (int i = 0; i < grad0.size(); i++)
                {
                    double gi = grad0[i];
                    gl += gi*gi;
                }

                gl = sqrt(gl);

                std::cout << "gl = " << gl << "\n";

                std::cout << "    optimizing locally...\n";

                std::vector<double> optPos = GradientDescent::optimize(
                    initial, lmf, 0.05/gl, 1e-9/gl, 1e-9, 10000, 0.0, false);

                for (int p = 0; p < pc; p++)
                {
                    for (int f = 0; f < fc; f++)
                    {
                        tracks[p][f].x = optPos[2*(p*fc + f)];
                        tracks[p][f].y = optPos[2*(p*fc + f) + 1];
                    }
                }
            }
        }

        if (evaluate)
        {
            if (measureFCC)
            {
                #pragma omp parallel for num_threads(nr_omp_threads)
                for (int p = 0; p < pc; p++)
                {
                    int threadnum = omp_get_thread_num();

                    Image<Complex> pred;
                    std::vector<Image<Complex>> obs = movie[p];

                    if (it_number > 0)
                    {
                        for (int f = 0; f < fc; f++)
                        {
                            shiftImageInFourierTransform(obs[f](), obs[f](), s, -tracks[p][f].x, -tracks[p][f].y);
                        }
                    }

                    int randSubset;
                    mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
                    randSubset -= 1;

                    if (randSubset == 0)
                    {
                        pred = obsModel.predictObservation(projector1, mdts[g], p, true, true);
                    }
                    else
                    {
                        pred = obsModel.predictObservation(projector0, mdts[g], p, true, true);
                    }

                    FscHelper::updateFscTable(obs, pred, tables[threadnum],
                                              weights0[threadnum], weights1[threadnum]);

                    std::vector<d2Vector> vel(fc);

                    vel[0] = tracks[p][1] - tracks[p][0];

                    for (int f = 1; f < fc-1; f++)
                    {
                        vel[f] = 0.5*(tracks[p][f+1] - tracks[p][f-1]);
                    }

                    vel[fc-1] = tracks[p][fc-1] - tracks[p][fc-2];

                    FscHelper::updateFscTableVelWgh(obs, vel, pred, tablesVW[threadnum],
                                              weights0VW[threadnum], weights1VW[threadnum]);


                    FscHelper::updateVelFscTable(
                        obs, vel, pred, tablesV[threadnum],
                        weights0V[threadnum], weights1V[threadnum], k_cutoff, k_out);
                }
            }
        }

        if (it_number > 0)
        {
            std::vector<std::vector<gravis::d2Vector>>
                    centTracks(pc), visTracks(pc), centVisTracks(pc);

            double visScale = 30.0;

            for (int p = 0; p < pc; p++)
            {
                centTracks[p] = std::vector<gravis::d2Vector>(fc);
                visTracks[p] = std::vector<gravis::d2Vector>(fc);
                centVisTracks[p] = std::vector<gravis::d2Vector>(fc);
            }

            std::vector<gravis::d2Vector> globalTrack(fc);

            for (int f = 0; f < fc; f++)
            {
                globalTrack[f] = d2Vector(0,0);

                for (int p = 0; p < pc; p++)
                {
                    globalTrack[f] += tracks[p][f];
                }

                globalTrack[f] /= pc;
                for (int p = 0; p < pc; p++)
                {
                    centTracks[p][f] = tracks[p][f] - globalTrack[f];
                    visTracks[p][f] = positions[p] + visScale * tracks[p][f];
                    centVisTracks[p][f] = positions[p] + visScale * centTracks[p][f];
                }
            }

            std::ofstream rawOut(outPath + "_mg" + stsg.str() + "_tracks.dat");
            std::ofstream visOut(outPath + "_mg" + stsg.str() + "_visTracks.dat");
            //std::ofstream centVisOut(outPath + "_mg" + stsg.str() + "_centVisTracks.dat");
            std::ofstream visOut15(outPath + "_mg" + stsg.str() + "_visTracks_first15.dat");

            for (int p = 0; p < pc; p++)
            {
                rawOut << "#particle " << p << "\n";
                visOut << "#particle " << p << "\n";
                //centVisOut << "#particle " << p << "\n";
                visOut15 << "#particle " << p << "\n";

                for (int f = 0; f < fc; f++)
                {
                    rawOut << tracks[p][f].x << " " << tracks[p][f].y << "\n";
                    visOut << visTracks[p][f].x << " " << visTracks[p][f].y << "\n";
                    //centVisOut << centVisTracks[p][f].x << " " << centVisTracks[p][f].y << "\n";

                    if (f < 15) visOut15 << visTracks[p][f].x << " " << visTracks[p][f].y << "\n";
                }

                rawOut << "\n";
                visOut << "\n";
                //centVisOut << "\n";
                visOut15 << "\n";
            }

            std::ofstream glbOut(outPath + "_mg" + stsg.str() + "_globTrack.dat");

            for (int f = 0; f < fc; f++)
            {
                glbOut << globalTrack[f].x << " " << globalTrack[f].y << "\n";
            }
        }

    } // micrographs

    if (evaluate)
    {
        if (measureFCC)
        {
            Image<RFLOAT> table, weight;
            Image<RFLOAT> tableV, weightV;
            Image<RFLOAT> tableVW, weightVW;

            FscHelper::mergeFscTables(tables, weights0, weights1, table, weight);
            FscHelper::mergeFscTables(tablesV, weights0V, weights1V, tableV, weightV);
            FscHelper::mergeFscTables(tablesVW, weights0VW, weights1VW, tableVW, weightVW);

            ImageLog::write(tableV, outPath + "_FCC_data_V");

            const int ccvBins = 10;
            std::vector<double> ccByV(ccvBins, 0.0);
            std::vector<double> ccByV_w(ccvBins, 0.0);

            for (int f = 0; f < 4; f++)
            for (int k = 0; k < sh; k++)
            {
                int b = k*ccvBins/sh;
                ccByV[b] += tableV(f,k)*weightV(f,k);
                ccByV_w[b] += weightV(f,k);
            }

            std::ofstream ccvOut(outPath + "_CC_by_V.dat");

            for (int b = 0; b < ccvBins; b++)
            {
                if (ccByV_w[b] > 0.0)
                {
                    ccByV[b] /= ccByV_w[b];
                }

                ccvOut << (b+0.5)*sh/ccvBins << " " << ccByV[b] << "\n";
            }

            int f_max = fc;
            double total = 0.0;
            double totalVW = 0.0;

            std::ofstream fccOut(outPath + "_FCC_perFrame.dat");
            std::ofstream fccOutVW(outPath + "_FCC_perFrame_velWgh_Gauss.dat");

            for (int y = 0; y < f_max; y++)
            {
                double avg = 0.0;
                double avgVW = 0.0;

                for (int k = k_cutoff+2; k < k_out; k++)
                {
                    avg += table(y,k);
                    avgVW += tableVW(y,k);
                }

                avg /= k_out - k_cutoff - 1;
                avgVW /= k_out - k_cutoff - 1;

                fccOut << y << " " << avg << "\n";
                fccOutVW << y << " " << avgVW << "\n";

                total += avg;
                totalVW += avgVW;
            }

            total /= f_max;
            totalVW /= f_max;

            std::cout << "total: " << total << " (" << totalVW <<")\n";

            table.write(outPath + "_FCC_data.mrc");
            /*weight.write(outPath + "_FCC_weight.mrc");*/

            ImageLog::write(table, outPath + "_FCC_data");
        }
    }

    double t1 = omp_get_wtime();
    double diff = t1 - t0;
    std::cout << "elapsed (total): " << diff << " sec\n";

    return 0;
}


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
#include <src/jaz/image_log.h>
#include <src/jaz/slice_helper.h>
#include <src/jaz/spectral_helper.h>
#include <src/jaz/filter_helper.h>
#include <src/jaz/backprojection_helper.h>
#include <src/jaz/volume_converter.h>
#include <src/jaz/complex_io.h>
#include <src/jaz/fftw_helper.h>
#include <src/jaz/resampling_helper.h>
#include <src/jaz/ctf_helper.h>
#include <src/jaz/defocus_refinement.h>
#include <src/jaz/magnification_refinement.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/stack_helper.h>
#include <src/jaz/tilt_refinement.h>
#include <src/jaz/image_op.h>
#include <src/jaz/Fourier_helper.h>
#include <src/jaz/fsc_helper.h>
#include <src/jaz/damage_helper.h>

#include <omp.h>

using namespace gravis;

#define BILLION 1000000000L

int main(int argc, char *argv[])
{
    std::string frcPathIn, frcPathOut, fn_sel, moviePath, fn_sym, reconFn;
    bool frcReady, applyTilt, computeFcc;
    int nr_omp_threads, nr_helical_asu, k_max;
    double beamtilt_x, beamtilt_y, helical_rise, helical_twist, totalDose, angpix;
    Image<RFLOAT> map;

    IOParser parser;

    try
    {
        parser.setCommandLine(argc, argv);

        parser.addSection("General options");

        frcPathIn = parser.getOption("--fsc_in", "Precomputed FSC/FCC table", "");
        totalDose = textToFloat(parser.getOption("--dose", "Total electron dose (in e^-/A^2)", "1"));
        angpix = textToFloat(parser.getOption("--angpix", "Pixel resolution (in A/pix)", "1"));

        frcReady = frcPathIn != "";


        parser.addSection("FSC/FCC computation");

        fn_sel = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
        moviePath = parser.getOption("--mov", "Path to movies", "");
        k_max = textToInteger(parser.getOption("--k_max", "Maximum freq. used for fit", "100"));

        nr_omp_threads = textToInteger(parser.getOption("--jomp", "Number of open-mp threads to use. Memory footprint is multiplied by this value.", "16"));

        frcPathOut = parser.getOption("--fsc_out", "Output path for FSC/FCC table", "");
        if (!frcReady && frcPathOut == "")
        {
            std::cout << "An output path is required for the FSC/FCC (--fsc_out).\n";
            return -1;
        }

        parser.addSection("FCC");

        reconFn = parser.getOption("--map", "Input MRC file of an initial reconstruction (when present, an FCC is computed istead of the FSC)", "");

        beamtilt_x = textToFloat(parser.getOption("--beamtilt_x", "Beamtilt in the X-direction (in mrad)", "0."));
        beamtilt_y = textToFloat(parser.getOption("--beamtilt_y", "Beamtilt in the Y-direction (in mrad)", "0."));
        applyTilt = (ABS(beamtilt_x) > 0. || ABS(beamtilt_y) > 0.);

        parser.addSection("FSC");

        fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
        nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
        helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
        helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

        computeFcc = reconFn != "" && !frcReady;

        if (computeFcc)
        {
            try
            {
                map.read(reconFn);
            }
            catch (RelionError XE)
            {
                std::cout << "Unable to read map: " << reconFn << "\n";
                exit(1);
            }
        }

    }
    catch (RelionError XE)
    {
        parser.writeUsage(std::cout);
        std::cerr << XE;
        exit(1);
    }

    //const bool individual = false;

    Image<RFLOAT> frcD, frcW;
    const bool fcc = true;

    if (frcReady)
    {
        frcD.read(frcPathIn+"_data.mrc");
        //frcW.read(frcPathIn+"_weight.mrc");
    }
    else
    {

        MetaDataTable DF0;
        DF0.read(fn_sel);

        double Cs, kV;

        DF0.getValue(EMDL_CTF_CS, Cs, 0);
        DF0.getValue(EMDL_CTF_VOLTAGE, kV, 0);

        double V = kV * 1e3;
        double lambda = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));

        std::vector<MetaDataTable> mdts = StackHelper::splitByStack(&DF0);
        int stack_count = mdts.size();

        std::vector<std::vector<Image<Complex> > > movie = StackHelper::loadMovieStackFS(&mdts[0], moviePath);
        int frame_count = movie[0].size();
        int sh = movie[0][0]().xdim;
        int s = movie[0][0]().ydim;

        std::cout << frame_count << " " << sh << "x" << s << " frames\n";

        RFLOAT mag, dstep;
        DF0.getValue(EMDL_CTF_MAGNIFICATION, mag, 0);
        DF0.getValue(EMDL_CTF_DETECTOR_PIXEL_SIZE, dstep, 0);
        double angpix = 10000. * dstep / mag;

        std::cout << "angpix: " << angpix << std::endl;

        const int data_dim = 2;
        const int ref_dim = 3;
        const int interpolator = TRILINEAR;
        double padding_factor = 1;
        const int r_min_nn = 10;

        MultidimArray<RFLOAT> dummy;
        Projector projector(s, interpolator, padding_factor, r_min_nn, data_dim);
        projector.computeFourierTransformMap(map.data, dummy, 2*s);

        if (fcc)
        {
            std::vector<Image<RFLOAT> > tables(nr_omp_threads);
            std::vector<Image<RFLOAT> > weights0(nr_omp_threads);
            std::vector<Image<RFLOAT> > weights1(nr_omp_threads);
            std::vector<bool> ready(nr_omp_threads,false);

            #pragma omp parallel for num_threads(nr_omp_threads)
            for (long g = 0; g < stack_count; g++)
            {
                int threadnum = omp_get_thread_num();

                std::cout << "micrograph " << g << " / " << mdts.size() <<"\n";

                std::vector<Image<Complex> > pred;

                if (nr_omp_threads > 1)
                {
                    pred = StackHelper::projectStackPar(&projector, &mdts[g], nr_omp_threads);
                }
                else
                {
                    pred = StackHelper::projectStack(&projector, &mdts[g]);
                }

                const int pc = pred.size();

                for (int p = 0; p < pc; p++)
                {
                    CTF ctf0;
                    ctf0.read(mdts[g], mdts[g], p);

                    FilterHelper::modulate(pred[p](), ctf0, angpix);

                    if (applyTilt)
                    {
                        selfApplyBeamTilt(pred[p].data, -beamtilt_x, -beamtilt_y, lambda, Cs, angpix, s);
                    }
                }

                std::vector<std::vector<Image<Complex> > > movie = StackHelper::loadMovieStackFS(&mdts[g], moviePath);

                if (!ready[threadnum])
                {
                    FscHelper::initFscTable(sh, movie[0].size(), tables[threadnum], weights0[threadnum], weights1[threadnum]);
                    ready[threadnum] = true;
                }

                FscHelper::updateFscTable(movie, pred, tables[threadnum], weights0[threadnum], weights1[threadnum]);
            }

            Image<RFLOAT> table, weight;

            FscHelper::mergeFscTables(tables, weights0, weights1, table, weight);

            table.write(frcPathOut+"_data.mrc");
            weight.write(frcPathOut+"_weight.mrc");

            ImageLog::write(table, frcPathOut+"_data");
            ImageLog::write(weight, frcPathOut+"_weight");
        }
        else
        {
            double blob_radius = 1.9;
            int blob_order = 0;
            double blob_alpha = 15.0;

            const bool skip_gridding = true;

            int r_max = -1;

            std::vector<std::vector<BackProjector> > backprojs(2);
            backprojs[0] = std::vector<BackProjector>(frame_count);
            backprojs[1] = std::vector<BackProjector>(frame_count);

            for (int j = 0; j < 2; j++)
            for (int i = 0; i < frame_count; i++)
            {
                backprojs[j][i] = BackProjector(
                            s, ref_dim, fn_sym, interpolator,
                            padding_factor, r_min_nn, blob_order,
                            blob_radius, blob_alpha, data_dim, skip_gridding);

                backprojs[j][i].initZeros(2 * r_max);
            }

            std::cerr << "Back-projecting all images...\n";

            //time_config();
            //init_progress_bar(gc);

            for (int g = 0; g < stack_count; g++)
            {
                std::cout << "   " << (g+1) << "/" << stack_count << "\n";

                if (g > 0)
                {
                    movie = StackHelper::loadMovieStackFS(&mdts[g], moviePath, true, nr_omp_threads);
                }

                const int particle_count = movie.size();

                for (int p = 0; p < particle_count; p++)
                {
                    std::cout << "      " << (p+1) << "/" << particle_count << "\n";

                    RFLOAT rot, tilt, psi;
                    Matrix2D<RFLOAT> A3D;
                    Matrix1D<RFLOAT> trans(2);

                    int randSubset;
                    mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
                    randSubset = randSubset - 1;

                    mdts[g].getValue(EMDL_ORIENT_ROT, rot, p);
                    mdts[g].getValue(EMDL_ORIENT_TILT, tilt, p);
                    mdts[g].getValue(EMDL_ORIENT_PSI, psi, p);

                    Euler_angles2matrix(rot, tilt, psi, A3D);

                    trans.initZeros();
                    mdts[g].getValue(EMDL_ORIENT_ORIGIN_X, XX(trans), p);
                    mdts[g].getValue(EMDL_ORIENT_ORIGIN_Y, YY(trans), p);

                    #pragma omp parallel for num_threads(nr_omp_threads)
                    for (int f = 0; f < frame_count; f++)
                    {
                        MultidimArray<Complex> F2D = movie[p][f]();

                        shiftImageInFourierTransform(F2D, F2D, s, XX(trans), YY(trans));

                        MultidimArray<RFLOAT> Fctf;
                        Fctf.resize(F2D);
                        Fctf.initConstant(1.);

                        CTF ctf;
                        ctf.read(mdts[g], mdts[g], p);
                        ctf.getFftwImage(Fctf, s, s, angpix, false, false, false, true);

                        if (applyTilt || (   mdts[g].containsLabel(EMDL_IMAGE_BEAMTILT_X)
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

                            selfApplyBeamTilt(F2D, beamtilt_x, beamtilt_y, ctf.lambda, ctf.Cs, angpix, sh);
                        }

                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D)
                        {
                            DIRECT_MULTIDIM_ELEM(F2D, n)  *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                            DIRECT_MULTIDIM_ELEM(Fctf, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                        }

                        backprojs[randSubset][f].set2DFourierTransform(F2D, A3D, IS_NOT_INV, &Fctf);
                    }
                }

            }

            frcD = Image<RFLOAT>(sh, frame_count);
            frcW = Image<RFLOAT>(sh, frame_count);

            std::cout << "computing fsc...\n";

            for (int f = 0; f < frame_count; f++)
            {
                for (int c = 0; c < 2; c++)
                {
                    backprojs[c][f].symmetrise(nr_helical_asu, helical_twist, helical_rise);
                    MultidimArray<Complex> dcopy = backprojs[c][f].data;
                    MultidimArray<RFLOAT> wcopy = backprojs[c][f].weight;

                    backprojs[c][f].decenter(dcopy, backprojs[c][f].data, sh*sh);
                    backprojs[c][f].decenter(wcopy, backprojs[c][f].weight, sh*sh);
                }

                FscHelper::computeFscRow(backprojs[0][f].data, backprojs[1][f].data, f, frcD, frcW);
            }

            frcD.write(frcPathOut+"_data.mrc");
            frcW.write(frcPathOut+"_weight.mrc");

            ImageLog::write(frcD, frcPathOut+"_data");
            ImageLog::write(frcW, frcPathOut+"_weight");
        }
    }

    const int kc = frcD.data.xdim;
    const int tc = frcD.data.ydim;

    const double dosePerFrame = totalDose/(double)tc;

    const std::string tag = "fcc";
    const double fsc_thresh = 0.001;// 0.143
    const double t0 = 10;
    const double k0 = 20;
    const double k1 = kc;
    const bool L1 = false, root = false;

    const double rho = 1.0 / (2.0 * (kc-1) * angpix);

    frcW = Image<RFLOAT>(kc, tc);

    for (int t = 0; t < tc; t++)
    for (int k = 0; k < kc; k++)
    {
        DIRECT_A2D_ELEM(frcW.data, t, k) = k < k_max? 1.0 : 0.0;

        /*if (DIRECT_A2D_ELEM(frcD.data, t, k) > fsc_thresh)
        {
            DIRECT_A2D_ELEM(frcW.data, t, k) = 1.0;
        }
        else
        {
            DIRECT_A2D_ELEM(frcW.data, t, k) = 0.0;
        }*/
    }

    ImageLog::write(frcW, "damage/frcMask_"+tag);
    ImageLog::write(frcD, "damage/frcTest_"+tag);


    std::vector<double> amp;
    double a, b, c;

    DamageHelper::fitGlobalDamage(frcD, frcW, amp, &a, &b, &c, k0, k1, t0, angpix, dosePerFrame, L1);

    std::ofstream os("damage/damage_glob_"+tag+".dat");

    for (int k = 1; k < kc; k++)
    {
        double dd = a*pow(rho*k,b)+c;
        double dout = dd > k_max? k_max : dd;
        os << rho*k << " " << dout << "\n";
    }

    os.close();

    ImageLog::write(frcD, "damage/glob_frcObs_"+tag);

    Image<RFLOAT> snrSynth = DamageHelper::plotGlobalDamage(a, b, c, amp, kc, tc, angpix, dosePerFrame, false);
    Image<RFLOAT> snrSynthN = DamageHelper::plotGlobalDamage(a, b, c, kc, tc, angpix, dosePerFrame, false);

    ImageLog::write(snrSynth, "damage/glob_snrSynth_"+tag);
    ImageLog::write(snrSynthN, "damage/glob_snrSynthN_"+tag);



    return 0;
}

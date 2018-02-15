
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
#include <src/jaz/vtk_helper.h>
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
#include <src/jaz/motion_refinement.h>
#include <src/jaz/image_op.h>
#include <src/jaz/refinement_program.h>
#include <src/jaz/damage_helper.h>
#include <src/jaz/fsc_helper.h>
#include <src/jaz/local_motion_fit.h>
#include <src/jaz/gradient_descent.h>
#include <src/jaz/distribution_helper.h>
#include <src/jaz/parallel_ft.h>

#include <src/jaz/motion_em.h>

#include <omp.h>

using namespace gravis;

class MotionFitProg : public RefinementProgram
{
    public:

        MotionFitProg();

            int evalFrames;
            RFLOAT dmga, dmgb, dmgc, totalDose,
                sig_vel, sig_div, sig_acc,
                k_cutoff;
            bool evaluate;

        int readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();
};

int main(int argc, char *argv[])
{
    MotionFitProg mt;

    int rc0 = mt.init(argc, argv);
    if (rc0 != 0) return rc0;

    int rc1 = mt.run();
    if (rc1 != 0) return rc1;
}

MotionFitProg::MotionFitProg()
:   RefinementProgram(false, true)
{
}

int MotionFitProg::readMoreOptions(IOParser& parser, int argc, char *argv[])
{
    dmga = textToFloat(parser.getOption("--dmg_a", "Damage model, parameter a", " 3.40"));
    dmgb = textToFloat(parser.getOption("--dmg_b", "                        b", "-1.06"));
    dmgc = textToFloat(parser.getOption("--dmg_c", "                        c", "-0.54"));

    totalDose = textToFloat(parser.getOption("--dose", "Total electron dose (in e^-/A^2)", "1"));

    sig_vel = textToFloat(parser.getOption("--s_vel", "Velocity sigma, other frames", "2.0"));
    sig_div = textToFloat(parser.getOption("--s_div", "Divergence sigma, other frames", "0.01"));
    sig_acc = textToFloat(parser.getOption("--s_acc", "Acceleration sigma", "-1.0"));

    k_cutoff = textToFloat(parser.getOption("--k_cut", "Freq. cutoff (in pixels)", "-1.0"));

    evalFrames = textToInteger(parser.getOption("--eval", "Measure FSC for this many initial frames", "0"));
    evaluate = evalFrames > 0;

    return 0;
}

int MotionFitProg::_init()
{
    return 0;
}

int MotionFitProg::_run()
{
    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    if (preextracted)
    {
        std::string name, fullName, movieName;
        mdts[0].getValue(EMDL_IMAGE_NAME, fullName, 0);
        mdts[0].getValue(EMDL_MICROGRAPH_NAME, movieName, 0);
        name = fullName.substr(fullName.find("@")+1);

        std::string finName;

        if (imgPath == "")
        {
            finName = name;
        }
        else
        {
            finName = imgPath + "/" + movieName.substr(movieName.find_last_of("/")+1);
        }

        Image<RFLOAT> stack0;
        stack0.read(finName, false);

        const int pc0 = mdts[0].numberOfObjects();
        const bool zstack = stack0.data.zdim > 1;
        const int stackSize = zstack? stack0.data.zdim : stack0.data.ndim;

        fc = stackSize / pc0;
    }
    else
    {
        std::vector<std::vector<Image<Complex>>> movie = StackHelper::extractMovieStackFS(
            &mdts[0], meta_path, imgPath, bin, coords_bin, movie_bin, s,
            nr_omp_threads, !nogain, binType, false, hotCutoff, debug);

        fc = movie[0].size();
    }

    std::vector<Image<RFLOAT> > dmgWeight = DamageHelper::damageWeights(
                s, angpix, fc, totalDose, dmga, dmgb, dmgc);

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


    double t0 = omp_get_wtime();

    int pctot = 0;

    std::vector<Image<RFLOAT> > tables(nr_omp_threads);
    std::vector<Image<RFLOAT> > weights0(nr_omp_threads);
    std::vector<Image<RFLOAT> > weights1(nr_omp_threads);

    if (evaluate)
    {
        for (int i = 0; i < nr_omp_threads; i++)
        {
            FscHelper::initFscTable(sh, fc, tables[i], weights0[i], weights1[i]);
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
            if (preextracted)
            {
                movie = StackHelper::loadMovieStackFS(
                    &mdts[g], imgPath, false, nr_omp_threads, &fts);
            }
            else
            {
                movie = StackHelper::extractMovieStackFS(
                    &mdts[g], meta_path, imgPath,
                    bin, coords_bin, movie_bin, s,
                    nr_omp_threads, !nogain, binType,
                    true, hotCutoff, debug);

                #pragma omp parallel for num_threads(nr_omp_threads)
                for (int p = 0; p < pc; p++)
                {
                    StackHelper::varianceNormalize(movie[p], false);
                }
            }
        }
        catch (RelionError XE)
        {
            std::cerr << "warning: unable to load micrograph #" << (g+1) << "\n";
            continue;
        }

        std::vector<double> sigma2 = StackHelper::powerSpectrum(movie);

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

        //@TODO - refactor:

        std::vector<RFLOAT> sig_vel_vec(1, sig_vel);
        std::vector<RFLOAT> sig_div_vec(1, sig_div);

        std::cout << "    computing initial correlations...\n";

        std::vector<std::vector<Image<RFLOAT>>> movieCC = MotionRefinement::movieCC(
                projectors[0], projectors[1], obsModel, mdts[g], movie,
                sigma2, dmgWeight, fts, nr_omp_threads);

        std::vector<std::vector<gravis::d2Vector>> tracks(pc);

        std::vector<gravis::d2Vector> globTrack = MotionRefinement::getGlobalTrack(movieCC);

        for (int p = 0; p < pc; p++)
        {
            tracks[p] = globTrack;
        }

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

        LocalMotionFit lmf(movieCC, velWgh, accWgh, divWgh, nr_omp_threads);

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

        lmf.grad(initial, grad0);

        double gl = 0.0;

        for (int i = 0; i < grad0.size(); i++)
        {
            double gi = grad0[i];
            gl += gi*gi;
        }

        gl = sqrt(gl);

        std::cout << "gl = " << gl << "\n";

        std::cout << "    optimizing...\n";

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

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            int threadnum = omp_get_thread_num();

            Image<Complex> pred;
            std::vector<Image<Complex>> obs = movie[p];

            for (int f = 0; f < fc; f++)
            {
                shiftImageInFourierTransform(obs[f](), obs[f](), s, -tracks[p][f].x, -tracks[p][f].y);
            }

            int randSubset;
            mdts[g].getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
            randSubset -= 1;

            if (randSubset == 0)
            {
                pred = obsModel.predictObservation(projectors[1], mdts[g], p, true, true);
            }
            else
            {
                pred = obsModel.predictObservation(projectors[0], mdts[g], p, true, true);
            }

            FscHelper::updateFscTable(obs, pred, tables[threadnum],
                                      weights0[threadnum], weights1[threadnum]);
        }

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
            std::ofstream visOut15(outPath + "_mg" + stsg.str() + "_visTracks_first15.dat");

            for (int p = 0; p < pc; p++)
            {
                rawOut << "#particle " << p << "\n";
                visOut << "#particle " << p << "\n";
                visOut15 << "#particle " << p << "\n";

                for (int f = 0; f < fc; f++)
                {
                    rawOut << tracks[p][f].x << " " << tracks[p][f].y << "\n";
                    visOut << visTracks[p][f].x << " " << visTracks[p][f].y << "\n";

                    if (f < 15) visOut15 << visTracks[p][f].x << " " << visTracks[p][f].y << "\n";
                }

                rawOut << "\n";
                visOut << "\n";
                visOut15 << "\n";
            }

            std::ofstream glbOut(outPath + "_mg" + stsg.str() + "_globTrack.dat");

            for (int f = 0; f < fc; f++)
            {
                glbOut << globalTrack[f].x << " " << globalTrack[f].y << "\n";
            }
        }

    } // micrographs

    Image<RFLOAT> table, weight;

    FscHelper::mergeFscTables(tables, weights0, weights1, table, weight);
    table.write(outPath + "_FCC_data.mrc");
    VtkHelper::writeVTK(table, outPath + "_FCC_data.vtk");

    if (evaluate)
    {
        int f_max = fc;
        double total = 0.0;

        std::ofstream fccOut(outPath + "_FCC_perFrame.dat");

        for (int y = 0; y < f_max; y++)
        {
            double avg = 0.0;

            for (int k = k_cutoff+2; k < k_out; k++)
            {
                avg += table(y,k);
            }

            avg /= k_out - k_cutoff - 1;

            fccOut << y << " " << avg << "\n";

            total += avg;
        }

        total /= f_max;

        std::cout << "total: " << total << "\n";
    }

    double t1 = omp_get_wtime();
    double diff = t1 - t0;
    std::cout << "elapsed (total): " << diff << " sec\n";

    return 0;
}

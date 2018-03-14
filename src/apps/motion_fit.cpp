
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

            bool unregGlob, noGlobOff;
            int maxIters;
            double dmga, dmgb, dmgc, dosePerFrame,
                sig_vel, sig_div, sig_acc,
                k_cutoff, maxStep, maxDistDiv;

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

    dosePerFrame = textToFloat(parser.getOption("--fdose", "Electron dose per frame (in e^-/A^2)", "1"));

    sig_vel = textToFloat(parser.getOption("--s_vel", "Velocity sigma [Angst/dose]", "1.6"));
    sig_div = textToFloat(parser.getOption("--s_div", "Divergence sigma [Angst/(sqrt(Angst)*dose)]", "0.067"));
    sig_acc = textToFloat(parser.getOption("--s_acc", "Acceleration sigma [Angst]", "-1.0"));

    maxDistDiv = textToFloat(parser.getOption("--max_dist_div", "Ignore neighbors if they are further away than this [pixels]", "-1.0"));

    k_cutoff = textToFloat(parser.getOption("--k_cut", "Freq. cutoff (in pixels)", "-1.0"));
    maxIters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "10000"));
    maxStep = textToFloat(parser.getOption("--max_step", "Maximum step size", "0.05"));

    unregGlob = parser.checkOption("--unreg_glob", "Do not regularize global component of motion");
    noGlobOff = parser.checkOption("--no_glob_off", "Do not compute initial per-particle offsets");

    return 0;
}

int MotionFitProg::_init()
{
    return 0;
}

int MotionFitProg::_run()
{
    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    loadInitialMovieValues();

    // @TODO: warn if angpix < coords_anpix or movie_angpix

    std::vector<Image<RFLOAT> > dmgWeight = DamageHelper::damageWeights(
        s, angpix, firstFrame, fc, dosePerFrame, dmga, dmgb, dmgc);

    const double sig_vel_nrm = dosePerFrame * sig_vel / angpix;
    const double sig_acc_nrm = dosePerFrame * sig_acc / angpix;
    const double sig_div_nrm = dosePerFrame * sqrt(coords_angpix) * sig_div / angpix;

    // @TODO: replace k_out by .143 res
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

    for (int i = 0; i < nr_omp_threads; i++)
    {
        FscHelper::initFscTable(sh, fc, tables[i], weights0[i], weights1[i]);
    }

    for (long g = g0; g <= gc; g++)
    {
        std::cout << "micrograph " << (g+1) << " / " << mdts.size() <<"\n";

        std::stringstream stsg;
        stsg << g;

        const int pc = mdts[g].numberOfObjects();
        pctot += pc;

        std::vector<std::vector<Image<Complex>>> movie;

        try
        {
            movie = loadMovie(g, pc, fts);
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

        std::cout << "    computing initial correlations...\n";

        std::vector<std::vector<Image<RFLOAT>>> movieCC = MotionRefinement::movieCC(
                projectors[0], projectors[1], obsModel, mdts[g], movie,
                sigma2, dmgWeight, fts, nr_omp_threads);

        // @TODO: make global offsets optional

        std::vector<Image<RFLOAT>> ccSum = MotionRefinement::addCCs(movieCC);
        std::vector<gravis::d2Vector> globTrack = MotionRefinement::getGlobalTrack(ccSum);
        std::vector<gravis::d2Vector> globOffsets;

        if (noGlobOff)
        {
            globOffsets = std::vector<d2Vector>(pc, d2Vector(0,0));
        }
        else
        {
            globOffsets = MotionRefinement::getGlobalOffsets(
                    movieCC, globTrack, 0.25*s, nr_omp_threads);
        }

        if (debug)
        {
            ImageLog::write(ccSum, outPath + "_CCsum_mg"+stsg.str(), CenterXY);
        }

        std::vector<std::vector<gravis::d2Vector>> tracks(pc);

        for (int p = 0; p < pc; p++)
        {
            tracks[p] = std::vector<d2Vector>(fc);

            for (int f = 0; f < fc; f++)
            {
                if (unregGlob)
                {
                    tracks[p][f] = globOffsets[p];
                }
                else
                {
                    tracks[p][f] = globTrack[f] + globOffsets[p];
                }
            }
        }

        std::vector<double> velWgh(fc-1, 0.5/(sig_vel_nrm*sig_vel_nrm));
        std::vector<double> accWgh(fc-1, sig_acc > 0.0? 0.5/(sig_acc_nrm*sig_acc_nrm) : 0.0);
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
                    double dist = sqrt(dp.x*dp.x + dp.y*dp.y);

                    if (q == p || (maxDistDiv >= 0.0 && dist > maxDistDiv))
                    {
                        divWgh[f][p][q] = 0.0;
                    }
                    else
                    {
                        divWgh[f][p][q] = 0.5 / (sig_div_nrm * sig_div_nrm * dist);
                    }
                }
            }
        }

        std::vector<d2Vector> globComp = unregGlob? globTrack : std::vector<d2Vector>(fc, d2Vector(0,0));

        LocalMotionFit lmf(movieCC, velWgh, accWgh, divWgh, globComp, nr_omp_threads);

        std::vector<double> initial(2*pc*fc);

        for (int p = 0; p < pc; p++)
        {
            for (int f = 0; f < fc; f++)
            {
                initial[2*(p*fc + f)]     = tracks[p][f].x;
                initial[2*(p*fc + f) + 1] = tracks[p][f].y;
            }
        }

        std::cout << "    optimizing...\n";

        std::vector<double> optPos = GradientDescent::optimize(
            initial, lmf, maxStep, 1e-9, 1e-9, maxIters, 0.0, debug);

        for (int p = 0; p < pc; p++)
        {
            for (int f = 0; f < fc; f++)
            {
                if (unregGlob)
                {
                    tracks[p][f].x = optPos[2*(p*fc + f)] + globTrack[f].x;
                    tracks[p][f].y = optPos[2*(p*fc + f) + 1] + globTrack[f].y;
                }
                else
                {
                    tracks[p][f].x = optPos[2*(p*fc + f)];
                    tracks[p][f].y = optPos[2*(p*fc + f) + 1];
                }
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
    ImageLog::write(table, outPath+"_FCC_data");


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

    double t1 = omp_get_wtime();
    double diff = t1 - t0;
    std::cout << "elapsed (total): " << diff << " sec\n";

    return 0;
}

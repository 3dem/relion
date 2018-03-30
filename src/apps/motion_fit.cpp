
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
#include <src/jaz/gp_motion_fit.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/distribution_helper.h>
#include <src/jaz/parallel_ft.h>
#include <src/jaz/d3x3/dsyev2.h>
#include <src/jaz/alignment_set.h>

#include <src/jaz/motion_em.h>

#include <omp.h>

using namespace gravis;

#define TIMING 1

#ifdef TIMING
    #define RCTIC(timer,label) (timer.tic(label))
    #define RCTOC(timer,label) (timer.toc(label))
#else
    #define RCTIC(timer,label)
    #define RCTOC(timer,label)
#endif


class MotionFitProg : public RefinementProgram
{
    public:

        MotionFitProg();

            bool unregGlob, noGlobOff,
                paramEstim2, paramEstim3,
                debugOpt, diag, expKer, global_init,
                useAlignmentSet, useLbfgs;

            int maxIters, paramEstimIters, paramEstimSteps, maxEDs, maxRange;

            double dmga, dmgb, dmgc, dosePerFrame,
                sig_vel, sig_div, sig_acc,
                k_cutoff, maxStep, maxDistDiv,
                param_rV, param_rD, param_rA,
                optEps;

            AlignmentSet alignmentSet;

            #ifdef TIMING
                Timer motionTimer;
                int timeInit, timeSetup0, timeSetup1, timeOpt, timeEval;
            #endif

        int readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();

        void prepAlignment(
                int k_out,
                const std::vector<Image<RFLOAT>>& dmgWeight);

        void prepMicrograph(
                int g,
                std::vector<ParFourierTransformer>& fts,
                const std::vector<Image<RFLOAT>>& dmgWeight,
                std::vector<std::vector<Image<Complex>>>& movie,
                std::vector<std::vector<Image<RFLOAT>>>& movieCC,
                std::vector<gravis::d2Vector>& positions,
                std::vector<std::vector<gravis::d2Vector>>& initialTracks,
                std::vector<d2Vector>& globComp);

        void estimateMotion(
                std::vector<ParFourierTransformer>& fts,
                const std::vector<Image<RFLOAT>>& dmgWeight);

        d2Vector estimateTwoParams(
                std::vector<ParFourierTransformer>& fts,
                const std::vector<Image<RFLOAT>>& dmgWeight,
                int k_out, double sig_v_0, double sig_d_0,
                double sig_v_step, double sig_d_step,
                int maxIters, int recDepth);

        d3Vector estimateThreeParams(
                std::vector<ParFourierTransformer>& fts,
                const std::vector<Image<RFLOAT>>& dmgWeight,
                int k_out, double sig_v_0, double sig_d_0, double sig_a_0,
                double sig_v_step, double sig_d_step, double sig_a_step,
                int maxIters, int recDepth);

        void evaluateParams(
                std::vector<ParFourierTransformer>& fts,
                const std::vector<Image<RFLOAT>>& dmgWeight,
                int k_out,
                const std::vector<d3Vector>& sig_vals,
                std::vector<double>& TSCs);

        void computeWeights(
                double sig_vel_nrm, double sig_acc_nrm, double sig_div_nrm,
                const std::vector<d2Vector>& positions, int fc,
                std::vector<double>& velWgh,
                std::vector<double>& accWgh,
                std::vector<std::vector<std::vector<double>>>& divWgh);

        std::vector<std::vector<d2Vector>> optimize(
                const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
                const std::vector<std::vector<d2Vector>>& inTracks,
                double sig_vel_px,
                double sig_acc_px,
                double sig_div_px,
                const std::vector<d2Vector>& positions,
                const std::vector<d2Vector>& globComp,
                double step, double minStep, double minDiff,
                long maxIters, double inertia);

        void updateFCC(const std::vector<std::vector<Image<Complex>>>& movie,
                const std::vector<std::vector<d2Vector>>& tracks,
                const MetaDataTable& mdt,
                std::vector<Image<RFLOAT>>& tables,
                std::vector<Image<RFLOAT>>& weights0,
                std::vector<Image<RFLOAT>>& weights1);

        void writeOutput(
                const std::vector<std::vector<d2Vector>>& tracks,
                const std::vector<Image<RFLOAT>>& fccData,
                const std::vector<Image<RFLOAT>>& fccWeight0,
                const std::vector<Image<RFLOAT>>& fccWeight1,
                const std::vector<d2Vector>& positions,
                std::string outPath, int mg,
                double visScale);
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

    sig_vel = textToFloat(parser.getOption("--s_vel", "Velocity sigma [Angst/dose]", "6"));
    sig_div = textToFloat(parser.getOption("--s_div", "Divergence sigma [Angst]", "720.0"));
    sig_acc = textToFloat(parser.getOption("--s_acc", "Acceleration sigma [Angst/dose]", "-1.0"));

    global_init = parser.checkOption("--gi", "Initialize with global trajectories instead of loading them from metadata file");

    expKer = parser.checkOption("--exp_k", "Use exponential kernel instead of sq. exponential");
    maxEDs = textToInteger(parser.getOption("--max_ed", "Maximum number of eigendeformations", "-1"));

    k_cutoff = textToFloat(parser.getOption("--k_cut", "Freq. cutoff (in pixels)", "-1.0"));
    maxIters = textToInteger(parser.getOption("--max_iters", "Maximum number of iterations", "1000"));

    unregGlob = parser.checkOption("--unreg_glob", "Do not regularize global component of motion");
    noGlobOff = parser.checkOption("--no_glob_off", "Do not compute initial per-particle offsets");

    debugOpt = parser.checkOption("--debug_opt", "Write optimization debugging info");
    diag = parser.checkOption("--diag", "Write out diagnostic data");

    parser.addSection("Parameter estimation");

    paramEstim2 = parser.checkOption("--params2", "Estimate 2 parameters instead of motion");
    paramEstim3 = parser.checkOption("--params3", "Estimate 3 parameters instead of motion");
    param_rV = textToFloat(parser.getOption("--r_vel", "Test s_vel +/- r_vel * s_vel", "0.5"));
    param_rD = textToFloat(parser.getOption("--r_div", "Test s_div +/- r_div * s_div", "0.5"));
    param_rA = textToFloat(parser.getOption("--r_acc", "Test s_acc +/- r_acc * s_acc", "0.5"));
    paramEstimIters = textToInteger(parser.getOption("--par_iters", "Parameter estimation is iterated this many times, each time halving the search range", "3"));
    paramEstimSteps = textToInteger(parser.getOption("--par_steps", "Parameter estimation takes max. this many steps before halving the range", "10"));
    useAlignmentSet = !parser.checkOption("--exact_params", "Use slower but more precise parameter evaluation");
    maxRange = textToInteger(parser.getOption("--mot_range", "Limit allowed motion range for parameter estimation [Px]", "50"));

    parser.addSection("Development options");

    optEps = textToFloat(parser.getOption("--eps", "Abort optimization once ||grad|| < eps * max(1, ||x||)", "1e-4"));
    useLbfgs = !parser.checkOption("--grad_desc", "Use gradient descent instead of LBFGS for optimization (don't!)");
    maxStep = textToFloat(parser.getOption("--max_step", "Maximum step size for gradient descent", "0.05"));

    if ((paramEstim2 || paramEstim3) && k_cutoff < 0)
    {
        std::cerr << "\nParameter estimation requires a freq. cutoff (--k_cut).\n";
        return 1138;
    }

    if (!global_init && corrMicFn == "")
    {
        std::cerr << "\nWarning: in the absence of a corrected_micrographs.star file (--corr_mic), global paths are used for initialization.\n";
        global_init = true;
    }

    if (paramEstim2 && paramEstim3)
    {
        std::cerr << "\nOnly 2 or 3 parameters can be estimated (--params2 or --params3), not both.\n";
        return 1139;
    }

    return 0;
}

int MotionFitProg::_init()
{
    return 0;
}

int MotionFitProg::_run()
{
    #ifdef TIMING
        timeInit = motionTimer.setNew(" time_Init ");
        timeSetup0 = motionTimer.setNew(" time_Setup0 ");
        timeSetup1 = motionTimer.setNew(" time_Setup1 ");
        timeOpt = motionTimer.setNew(" time_Opt ");
        timeEval = motionTimer.setNew(" time_Eval ");
    #endif

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    loadInitialMovieValues();

    std::vector<Image<RFLOAT>> dmgWeight = DamageHelper::damageWeights(
        s, angpix, firstFrame, fc, dosePerFrame, dmga, dmgb, dmgc);

    int k_out = sh;

    for (int i = 1; i < sh; i++)
    {
        if (freqWeight1D[i] <= 0.0)
        {
            k_out = i;
            break;
        }
    }

    std::cout << "max freq. = " << k_out << " px\n";

    for (int f = 0; f < fc; f++)
    {
        dmgWeight[f].data.xinit = 0;
        dmgWeight[f].data.yinit = 0;

        if (k_cutoff > 0.0)
        {
            std::stringstream stsf;
            stsf << f;
            dmgWeight[f] = FilterHelper::ButterworthEnvFreq2D(dmgWeight[f], k_cutoff-1, k_cutoff+1);

            ImageOp::multiplyBy(dmgWeight[f], freqWeight);
        }
    }

    if (useAlignmentSet && (paramEstim2 || paramEstim3))
    {
        prepAlignment(k_out, dmgWeight);
    }

    double t0 = omp_get_wtime();

    if (paramEstim2)
    {
        d2Vector par = estimateTwoParams(
            fts, dmgWeight, k_out, sig_vel, sig_div,
            std::abs(sig_vel*param_rV),
            std::abs(sig_div*param_rD),
            paramEstimSteps, paramEstimIters);

        std::cout << "opt = " << par << "\n";
    }
    else if (paramEstim3)
    {
        d3Vector par = estimateThreeParams(
            fts, dmgWeight, k_out,
            sig_vel, sig_div, sig_acc,
            std::abs(sig_vel*param_rV),
            std::abs(sig_div*param_rD),
            std::abs(sig_acc*param_rA),
            paramEstimSteps, paramEstimIters);

        std::cout << "opt = " << par << "\n";
    }
    else
    {
        estimateMotion(fts, dmgWeight);
    }


    double t1 = omp_get_wtime();
    double diff = t1 - t0;
    std::cout << "elapsed (total): " << diff << " sec\n\n";

    #ifdef TIMING
        motionTimer.printTimes(true);
    #endif

    std::cout << "\n";

    return 0;
}

void MotionFitProg::prepAlignment(int k_out, const std::vector<Image<RFLOAT>>& dmgWeight)
{
    RCTIC(motionTimer,timeInit);

    std::cout << "preparing alignment data...\n";

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    std::vector<MetaDataTable> mdtsAct(gc - g0 + 1);

    for (int m = g0; m <= gc; m++)
    {
        mdtsAct[m-g0] = mdts[m];
    }

    alignmentSet = AlignmentSet(mdtsAct, fc, s, k_cutoff+2, k_out);

    int pctot = 0;

    for (long g = g0; g <= gc; g++)
    {
        const int pc = mdts[g].numberOfObjects();

        if (pc < 2) continue;

        pctot += pc;

        std::cout << "    micrograph " << (g+1) << " / " << mdts.size() << ": "
            << pc << " particles [" << pctot << " total]\n";

        std::vector<std::vector<Image<Complex>>> movie;
        std::vector<std::vector<Image<RFLOAT>>> movieCC;

        try
        {
            prepMicrograph(
                g, fts, dmgWeight,
                movie, movieCC,
                alignmentSet.positions[g],
                alignmentSet.initialTracks[g],
                alignmentSet.globComp[g]);
        }
        catch (RelionError XE)
        {
            std::cerr << "warning: unable to load micrograph #" << (g+1) << "\n";
            continue;
        }

        Image<float> CCfloat(s,s);

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            for (int f = 0; f < fc; f++)
            {
                if (maxRange > 0)
                {
                    movieCC[p][f] = FilterHelper::cropCorner2D(movieCC[p][f], 2*maxRange, 2*maxRange);
                }

                for (int y = 0; y < movieCC[p][f]().ydim; y++)
                for (int x = 0; x < movieCC[p][f]().xdim; x++)
                {
                    CCfloat(y,x) = movieCC[p][f](y,x);
                }

                alignmentSet.CCs[g][p][f] = movieCC[p][f];
                //alignmentSet.CCs[g][p][f] = CCfloat;
                alignmentSet.obs[g][p][f] = alignmentSet.accelerate(movie[p][f]);

                Image<Complex> pred;
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

                alignmentSet.pred[g][p] = alignmentSet.accelerate(pred);
            }
        }
    }

    std::cout << "done\n";

    RCTOC(motionTimer,timeInit);
}

void MotionFitProg::prepMicrograph(
        int g, std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<RFLOAT>>& dmgWeight,
        std::vector<std::vector<Image<Complex>>>& movie,
        std::vector<std::vector<Image<RFLOAT>>>& movieCC,
        std::vector<d2Vector>& positions,
        std::vector<std::vector<d2Vector>>& initialTracks,
        std::vector<d2Vector>& globComp)
{
    const int pc = mdts[g].numberOfObjects();

    movie = loadMovie(g, pc, fts); // throws exceptions

    std::vector<double> sigma2 = StackHelper::powerSpectrum(movie);

    #pragma omp parallel for num_threads(nr_omp_threads)
    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++)
    {
        MotionRefinement::noiseNormalize(movie[p][f], sigma2, movie[p][f]);
    }

    positions = std::vector<gravis::d2Vector>(pc);

    for (int p = 0; p < pc; p++)
    {
        mdts[g].getValue(EMDL_IMAGE_COORD_X, positions[p].x, p);
        mdts[g].getValue(EMDL_IMAGE_COORD_Y, positions[p].y, p);
    }

    if (!paramEstim2 && !paramEstim3)
    {
        std::cout << "    computing initial correlations...\n";
    }

    movieCC = MotionRefinement::movieCC(
            projectors[0], projectors[1], obsModel, mdts[g], movie,
            sigma2, dmgWeight, fts, nr_omp_threads);

    initialTracks = std::vector<std::vector<d2Vector>>(pc);

    if (global_init)
    {
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

        if (diag)
        {
            std::string tag = outPath + "/" + getMicrographTag(g);
            std::string path = tag.substr(0, tag.find_last_of('/'));
            mktree(path);

            ImageLog::write(ccSum, tag + "_CCsum", CenterXY);
        }

        for (int p = 0; p < pc; p++)
        {
            initialTracks[p] = std::vector<d2Vector>(fc);

            for (int f = 0; f < fc; f++)
            {
                if (unregGlob)
                {
                    initialTracks[p][f] = globOffsets[p];
                }
                else
                {
                    initialTracks[p][f] = globTrack[f] + globOffsets[p];
                }
            }
        }

        globComp = unregGlob? globTrack : std::vector<d2Vector>(fc, d2Vector(0,0));
    }
    else
    {
        const d2Vector inputScale(
                coords_angpix / (movie_angpix * micrograph.getWidth()),
                coords_angpix / (movie_angpix * micrograph.getHeight()));

        const double outputScale = movie_angpix / angpix;

        globComp = std::vector<d2Vector>(fc, d2Vector(0,0));

        if (unregGlob)
        {
            for (int f = 0; f < fc; f++)
            {
                RFLOAT sx, sy;
                micrograph.getShiftAt(f+1, 0, 0, sx, sy, false);

                globComp[f] = -outputScale * d2Vector(sx, sy);
            }
        }

        for (int p = 0; p < pc; p++)
        {
            initialTracks[p] = std::vector<d2Vector>(fc);

            for (int f = 0; f < fc; f++)
            {
                d2Vector in(inputScale.x * positions[p].x - 0.5,
                            inputScale.y * positions[p].y - 0.5);

                RFLOAT sx, sy;

                micrograph.getShiftAt(f+1, in.x, in.y, sx, sy, true);

                initialTracks[p][f] = -outputScale * d2Vector(sx,sy) - globComp[f];
            }
        }
    }
}

void MotionFitProg::estimateMotion(
        std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<RFLOAT>>& dmgWeight)
{
    std::vector<Image<RFLOAT>>
            tables(nr_omp_threads),
            weights0(nr_omp_threads),
            weights1(nr_omp_threads);

    for (int i = 0; i < nr_omp_threads; i++)
    {
        FscHelper::initFscTable(sh, fc, tables[i], weights0[i], weights1[i]);
    }

    const double sig_vel_nrm = dosePerFrame * sig_vel / angpix;
    const double sig_acc_nrm = dosePerFrame * sig_acc / angpix;
    const double sig_div_nrm = dosePerFrame * sig_div / coords_angpix;

    int pctot = 0;

    // initialize parameter-estimation:

    for (long g = g0; g <= gc; g++)
    {        
        const int pc = mdts[g].numberOfObjects();

        std::cout << "micrograph " << (g+1) << " / " << mdts.size()
                  << ": " << pc << " particles\n";

        if (pc < 2) continue;

        std::stringstream stsg;
        stsg << g;

        std::vector<std::vector<Image<Complex>>> movie;
        std::vector<std::vector<Image<RFLOAT>>> movieCC;
        std::vector<d2Vector> positions;
        std::vector<std::vector<d2Vector>> initialTracks;
        std::vector<d2Vector> globComp;

        try
        {
            prepMicrograph(
                g, fts, dmgWeight,
                movie, movieCC, positions, initialTracks, globComp);
        }
        catch (RelionError XE)
        {
            std::cerr << "warning: unable to load micrograph #" << (g+1) << "\n";
            continue;
        }

        pctot += pc;

        std::cout << "    optimizing...\n";

        std::vector<std::vector<gravis::d2Vector>> tracks = optimize(
                movieCC, initialTracks,
                sig_vel_nrm, sig_acc_nrm, sig_div_nrm,
                positions, globComp,
                maxStep, 1e-9, 1e-9, maxIters, 0.0);

        updateFCC(movie, tracks, mdts[g], tables, weights0, weights1);

        writeOutput(tracks, tables, weights0, weights1, positions, outPath, g, 30.0);

        for (int i = 0; i < nr_omp_threads; i++)
        {
            tables[i].data.initZeros();
            weights0[i].data.initZeros();
            weights1[i].data.initZeros();
        }

    } // micrographs
}

d2Vector MotionFitProg::estimateTwoParams(
        std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<RFLOAT>>& dmgWeight,
        int k_out, double sig_v_0, double sig_d_0,
        double sig_v_step, double sig_d_step,
        int maxIters, int recDepth)
{
    std::cout << "estimating sig_v and sig_d\n";
    std::cout << "sig_v = " << sig_v_0 << " +/- " << sig_v_step << "\n";
    std::cout << "sig_d = " << sig_d_0 << " +/- " << sig_d_step << "\n";

    int tot_vi = sig_v_0 <= 0.0? (int)(-sig_v_0/sig_v_step + 1) : 0;
    int tot_di = sig_d_0 <= 0.0? (int)(-sig_d_0/sig_d_step + 1) : 0;

    std::vector<d3Vector> all_sig_vals(9);

    for (int p = 0; p < 9; p++)
    {
        int vi = (p%3) - 1;
        int di = (p/3) - 1;

        all_sig_vals[p][0] = sig_v_0 + (vi+tot_vi) * sig_v_step;
        all_sig_vals[p][1] = sig_d_0 + (di+tot_di) * sig_d_step;
        all_sig_vals[p][2] = sig_acc;
    }

    std::vector<double> all_TSCs(9);

    std::vector<d3Vector> unknown_sig_vals = all_sig_vals;
    std::vector<double> unknown_TSCs(9);

    std::vector<int> unknown_ind(9);

    for (int p = 0; p < 9; p++)
    {
        unknown_ind[p] = p;
    }

    bool centerBest = false;
    int iters = 0;

    d3Vector bestParams = all_sig_vals[4];

    const bool verbose = true;

    while (!centerBest && iters < maxIters)
    {
        if (verbose)
        {
            std::cout << "\nevaluating:\n";

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    std::cout << all_sig_vals[3*i + j] << " \t ";
                }

                std::cout << "\n";
            }

            std::cout << "\n";
        }

        evaluateParams(fts, dmgWeight, k_out, unknown_sig_vals, unknown_TSCs);

        for (int p = 0; p < unknown_ind.size(); p++)
        {
            all_TSCs[unknown_ind[p]] = unknown_TSCs[p];
        }

        if (verbose)
        {
            std::cout << "\nresult:\n";

            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    std::cout << all_TSCs[3*i + j] << " \t ";
                }

                std::cout << "\n";
            }

            std::cout << "\n";
        }

        int bestIndex = 0;
        double bestTSC = all_TSCs[0];

        for (int p = 0; p < 9; p++)
        {
            if (all_TSCs[p] > bestTSC)
            {
                bestTSC = all_TSCs[p];
                bestIndex = p;
            }
        }

        bestParams = all_sig_vals[bestIndex];

        int shift_v = bestIndex % 3 - 1;
        int shift_d = bestIndex / 3 - 1;

        if (verbose)
        {
            std::cout << "optimum: " << all_sig_vals[bestIndex] << " ";
            std::cout << "(" << all_sig_vals[bestIndex] << ")\n";
            std::cout << "grid_0: [" << shift_v << ", " << shift_d << "]\n";
        }

        if (shift_v < 0 && all_sig_vals[3][0] <= 0.0)
        {
            shift_v = 0;
        }

        if (shift_d < 0 && all_sig_vals[1][1] <= 0.0)
        {
            shift_d = 0;
        }

        if (verbose)
        {
            std::cout << "grid_1: [" << shift_v << ", " << shift_d << "]\n";
        }

        if (shift_v == 0 && shift_d == 0)
        {
            if (recDepth <= 0)
            {
                return d2Vector(all_sig_vals[bestIndex][0], all_sig_vals[bestIndex][1]);
            }
            else
            {
                std::cout << "\nrepeating at half scan range.\n\n";

                estimateTwoParams(
                    fts, dmgWeight, k_out,
                    all_sig_vals[bestIndex][0], all_sig_vals[bestIndex][1],
                    sig_v_step/2.0, sig_d_step/2.0,
                    maxIters, recDepth-1);
            }
        }

        tot_vi += shift_v;
        tot_di += shift_d;

        std::vector<d3Vector> next_sig_vals(9, d3Vector(0,0,0));
        std::vector<double> next_TSCs(9);
        std::vector<bool> known(9, false);

        for (int p = 0; p < 9; p++)
        {
            int vi = (p%3) - 1;
            int di = (p/3) - 1;

            int vi_next = vi - shift_v;
            int di_next = di - shift_d;

            if (vi_next >= -1 && vi_next <= 1
                && di_next >= -1 && di_next <= 1)
            {
                int p_next = (vi_next + 1) + 3*(di_next + 1);

                next_sig_vals[p_next] = all_sig_vals[p];
                next_TSCs[p_next] = all_TSCs[p];
                known[p_next] = true;
            }
        }

        all_sig_vals = next_sig_vals;
        all_TSCs = next_TSCs;

        unknown_sig_vals.clear();
        unknown_ind.clear();

        for (int p = 0; p < 9; p++)
        {
            if (!known[p])
            {
                int vi = (p%3) - 1 + tot_vi;
                int di = (p/3) - 1 + tot_di;

                all_sig_vals[p][0] = sig_v_0 + vi * sig_v_step;
                all_sig_vals[p][1] = sig_d_0 + di * sig_d_step;
                all_sig_vals[p][2] = sig_acc;

                unknown_sig_vals.push_back(all_sig_vals[p]);
                unknown_ind.push_back(p);
            }
        }

        unknown_TSCs.resize(unknown_ind.size());

        iters++;
    }

    return d2Vector(bestParams[0], bestParams[1]);
}


d3Vector MotionFitProg::estimateThreeParams(
        std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<RFLOAT>>& dmgWeight,
        int k_out, double sig_v_0, double sig_d_0, double sig_a_0,
        double sig_v_step, double sig_d_step, double sig_a_step,
        int maxIters, int recDepth)
{
    int tot_vi = sig_v_0 <= 0.0? (int)(-sig_v_0/sig_v_step + 1) : 0;
    int tot_di = sig_d_0 <= 0.0? (int)(-sig_d_0/sig_d_step + 1) : 0;
    int tot_ai = sig_a_0 <= 0.0? (int)(-sig_a_0/sig_a_step + 1) : 0;

    std::vector<d3Vector> all_sig_vals(27);

    for (int p = 0; p < 27; p++)
    {
        int vi = (p%3) - 1;
        int di = ((p/3)%3) - 1;
        int ai = p/9 - 1;

        all_sig_vals[p][0] = sig_v_0 + (vi+tot_vi) * sig_v_step;
        all_sig_vals[p][1] = sig_d_0 + (di+tot_di) * sig_d_step;
        all_sig_vals[p][2] = sig_a_0 + (ai+tot_ai) * sig_a_step;
    }

    std::vector<double> all_TSCs(27);

    std::vector<d3Vector> unknown_sig_vals = all_sig_vals;
    std::vector<double> unknown_TSCs(27);

    std::vector<int> unknown_ind(27);

    for (int p = 0; p < 27; p++)
    {
        unknown_ind[p] = p;
    }

    bool centerBest = false;
    int iters = 0;

    d3Vector bestParams = all_sig_vals[13];

    const bool verbose = true;

    while (!centerBest && iters < maxIters)
    {
        if (verbose)
        {
            std::cout << "evaluating:\n";

            for (int k = 0; k < 3; k++)
            {
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        std::cout << all_sig_vals[9*k + 3*i + j] << " \t ";
                    }

                    std::cout << "\n";
                }

                std::cout << "\n";
            }

            std::cout << "\n";
        }

        evaluateParams(fts, dmgWeight, k_out, unknown_sig_vals, unknown_TSCs);

        for (int p = 0; p < unknown_ind.size(); p++)
        {
            all_TSCs[unknown_ind[p]] = unknown_TSCs[p];
        }

        if (verbose)
        {
            std::cout << "result:\n";

            for (int k = 0; k < 3; k++)
            {
                for (int i = 0; i < 3; i++)
                {
                    for (int j = 0; j < 3; j++)
                    {
                        std::cout << all_TSCs[9*k + 3*i + j] << " \t ";
                    }

                    std::cout << "\n";
                }

                std::cout << "\n";
            }

            std::cout << "\n";
        }

        int bestIndex = 0;
        double bestTSC = all_TSCs[0];

        for (int p = 0; p < 27; p++)
        {
            if (all_TSCs[p] > bestTSC)
            {
                bestTSC = all_TSCs[p];
                bestIndex = p;
            }
        }

        bestParams = all_sig_vals[bestIndex];

        int shift_v = bestIndex % 3 - 1;
        int shift_d = ((bestIndex / 3)%3) - 1;
        int shift_a = bestIndex / 9 - 1;

        if (shift_v < 0 && all_sig_vals[12][0] <= 0.0)
        {
            shift_v = 0;
        }

        if (shift_d < 0 && all_sig_vals[10][1] <= 0.0)
        {
            shift_d = 0;
        }

        if (shift_a < 0 && all_sig_vals[4][2] <= 0.0)
        {
            shift_a = 0;
        }

        if (shift_v == 0 && shift_d == 0 && shift_a == 0)
        {
            if (recDepth <= 0)
            {
                return all_sig_vals[bestIndex];
            }
            else
            {
                std::cout << "current optimum: " << all_sig_vals[bestIndex] << "\n";
                std::cout << "repeating at half scan range.\n";

                estimateThreeParams(
                    fts, dmgWeight, k_out,
                    all_sig_vals[bestIndex][0], all_sig_vals[bestIndex][1], all_sig_vals[bestIndex][2],
                    sig_v_step/2.0, sig_d_step/2.0, sig_a_step/2.0,
                    maxIters, recDepth-1);
            }
        }

        tot_vi += shift_v;
        tot_di += shift_d;
        tot_ai += shift_a;

        std::vector<d3Vector> next_sig_vals(27, d3Vector(0,0,0));
        std::vector<double> next_TSCs(27);
        std::vector<bool> known(27, false);

        for (int p = 0; p < 27; p++)
        {
            int vi = (p%3) - 1;
            int di = ((p/3)%3) - 1;
            int ai = p/9 - 1;

            int vi_next = vi - shift_v;
            int di_next = di - shift_d;
            int ai_next = ai - shift_a;

            if (vi_next >= -1 && vi_next <= 1
                && di_next >= -1 && di_next <= 1
                && ai_next >= -1 && ai_next <= 1)
            {
                int p_next = (vi_next + 1) + 3*(di_next + 1) + 9*(ai_next + 1);

                next_sig_vals[p_next] = all_sig_vals[p];
                next_TSCs[p_next] = all_TSCs[p];
                known[p_next] = true;
            }
        }

        all_sig_vals = next_sig_vals;
        all_TSCs = next_TSCs;

        unknown_sig_vals.clear();
        unknown_ind.clear();

        for (int p = 0; p < 27; p++)
        {
            if (!known[p])
            {
                int vi = (p%3) - 1 + tot_vi;
                int di = ((p/3)%3) - 1 + tot_di;
                int ai = (p/9) - 1 + tot_ai;

                all_sig_vals[p][0] = sig_v_0 + vi * sig_v_step;
                all_sig_vals[p][1] = sig_d_0 + di * sig_d_step;
                all_sig_vals[p][2] = sig_a_0 + ai * sig_a_step;

                unknown_sig_vals.push_back(all_sig_vals[p]);
                unknown_ind.push_back(p);
            }
        }

        unknown_TSCs.resize(unknown_ind.size());

        iters++;
    }

    return bestParams;
}

void MotionFitProg::evaluateParams(
        std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<RFLOAT>>& dmgWeight,
        int k_out,
        const std::vector<d3Vector>& sig_vals,
        std::vector<double>& TSCs)
{    
    RCTIC(motionTimer,timeSetup0);

    const int paramCount = sig_vals.size();
    TSCs.resize(paramCount);

    std::vector<double> sig_v_vals_nrm(paramCount);
    std::vector<double> sig_d_vals_nrm(paramCount);
    std::vector<double> sig_a_vals_nrm(paramCount);

    for (int i = 0; i < paramCount; i++)
    {
        sig_v_vals_nrm[i] = dosePerFrame * sig_vals[i][0] / angpix;
        sig_d_vals_nrm[i] = dosePerFrame * sig_vals[i][1] / coords_angpix;
        sig_a_vals_nrm[i] = dosePerFrame * sig_vals[i][2] / angpix;
    }

    std::vector<std::vector<Image<RFLOAT>>>
        paramTables(paramCount), paramWeights0(paramCount), paramWeights1(paramCount);

    int pctot = 0;

    std::vector<d3Vector> tscsAs(paramCount, d3Vector(0.0, 0.0, 0.0));

    if (!useAlignmentSet)
    {
        for (int i = 0; i < paramCount; i++)
        {
            paramTables[i] = std::vector<Image<RFLOAT>>(nr_omp_threads);
            paramWeights0[i] = std::vector<Image<RFLOAT>>(nr_omp_threads);
            paramWeights1[i] = std::vector<Image<RFLOAT>>(nr_omp_threads);

            for (int j = 0; j < nr_omp_threads; j++)
            {
                FscHelper::initFscTable(sh, fc, paramTables[i][j],
                    paramWeights0[i][j], paramWeights1[i][j]);
            }
        }
    }

    RCTOC(motionTimer,timeSetup0);

    for (long g = g0; g <= gc; g++)
    {
        const int pc = mdts[g].numberOfObjects();

        if (pc < 2) continue;

        pctot += pc;

        std::cout << "    micrograph " << (g+1) << " / " << mdts.size() << ": "
            << pc << " particles [" << pctot << " total]\n";

        std::vector<std::vector<Image<Complex>>> movie;
        std::vector<std::vector<Image<RFLOAT>>> movieCC;
        std::vector<d2Vector> positions;
        std::vector<std::vector<d2Vector>> initialTracks;
        std::vector<d2Vector> globComp;

        if (!useAlignmentSet)
        {
            try
            {
                prepMicrograph(
                    g, fts, dmgWeight,
                    movie, movieCC, positions, initialTracks, globComp);
            }
            catch (RelionError XE)
            {
                std::cerr << "warning: unable to load micrograph #" << (g+1) << "\n";
                continue;
            }
        }

        for (int i = 0; i < paramCount; i++)
        {
            if (debug)
            {
                std::cout << "        evaluating: " << sig_vals[i] << "\n";
            }

            if (useAlignmentSet)
            {
                RCTIC(motionTimer,timeSetup1);

                /*std::vector<std::vector<Image<RFLOAT>>> CCs(pc, std::vector<Image<RFLOAT>>(fc));

                #pragma omp parallel for num_threads(nr_omp_threads)
                for (int p = 0; p < pc; p++)
                for (int f = 0; f < fc; f++)
                {
                    const int w = alignmentSet.CCs[g][p][f].data.xdim;
                    const int h = alignmentSet.CCs[g][p][f].data.ydim;

                    CCs[p][f] = Image<RFLOAT>(w,h);

                    for (int y = 0; y < h; y++)
                    for (int x = 0; x < w; x++)
                    {
                        CCs[p][f](y,x) = alignmentSet.CCs[g][p][f](y,x);
                    }
                }*/

                RCTOC(motionTimer,timeSetup1);

                std::vector<std::vector<gravis::d2Vector>> tracks = optimize(
                        alignmentSet.CCs[g],
                        //CCs,
                        alignmentSet.initialTracks[g],
                        sig_v_vals_nrm[i], sig_a_vals_nrm[i], sig_d_vals_nrm[i],
                        alignmentSet.positions[g], alignmentSet.globComp[g],
                        maxStep, 1e-9, 1e-9, maxIters, 0.0);

                RCTIC(motionTimer,timeEval);

                tscsAs[i] += alignmentSet.updateTsc(tracks, g, nr_omp_threads);
            }
            else
            {
                RCTIC(motionTimer,timeSetup1);
                RCTOC(motionTimer,timeSetup1);

                std::vector<std::vector<gravis::d2Vector>> tracks = optimize(
                        movieCC, initialTracks,
                        sig_v_vals_nrm[i], sig_a_vals_nrm[i], sig_d_vals_nrm[i],
                        positions, globComp, maxStep, 1e-9, 1e-9, maxIters, 0.0);

                RCTIC(motionTimer,timeEval);

                updateFCC(movie, tracks, mdts[g], paramTables[i], paramWeights0[i], paramWeights1[i]);
            }
        }

    } // micrographs

    if (useAlignmentSet)
    {
        for (int i = 0; i < paramCount; i++)
        {
            double wg = tscsAs[i][1] * tscsAs[i][2];

            if (wg > 0.0)
            {
                TSCs[i] = tscsAs[i][0] / sqrt(wg);
            }
        }
    }
    else
    {
        for (int i = 0; i < paramCount; i++)
        {
            TSCs[i] = FscHelper::computeTsc(
                paramTables[i], paramWeights0[i], paramWeights1[i], k_cutoff+2, k_out);
        }
    }


    RCTOC(motionTimer,timeEval);
}

std::vector<std::vector<d2Vector>> MotionFitProg::optimize(
        const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
        const std::vector<std::vector<d2Vector>>& inTracks,
        double sig_vel_px, double sig_acc_px, double sig_div_px,
        const std::vector<d2Vector>& positions,
        const std::vector<d2Vector>& globComp,
        double step, double minStep, double minDiff,
        long maxIters, double inertia)
{
    const double eps = 1e-20;

    if (sig_vel_px < eps)
    {
        if (debug) std::cerr << "Warning: sig_vel < " << eps << " px. Setting to " << eps << ".\n";
        sig_vel_px = eps;
    }

    if (sig_div_px < eps)
    {
        if (debug) std::cerr << "Warning: sig_div < " << eps << " px. Setting to " << eps << ".\n";
        sig_div_px = eps;
    }

    const int pc = inTracks.size();

    if (pc == 0) return std::vector<std::vector<d2Vector>>(0);

    const int fc = inTracks[0].size();


    RCTIC(motionTimer,timeOpt);

    GpMotionFit gpmf(movieCC, sig_vel_px, sig_div_px, sig_acc_px,
                     maxEDs, positions,
                     globComp, nr_omp_threads, expKer);


    std::vector<double> initialCoeffs(2*(pc + pc*(fc-1)));

    gpmf.posToParams(inTracks, initialCoeffs);

    std::vector<double> optCoeffs;

    if (useLbfgs)
    {
        optCoeffs = LBFGS::optimize(
            initialCoeffs, gpmf, debugOpt, maxIters, optEps);
    }
    else
    {
        optCoeffs = GradientDescent::optimize(
            initialCoeffs, gpmf, step, minStep, minDiff, maxIters, inertia, debugOpt);
    }

    std::vector<std::vector<d2Vector>> out(pc, std::vector<d2Vector>(fc));
    gpmf.paramsToPos(optCoeffs, out);

    RCTOC(motionTimer,timeOpt);

    return out;
}

void MotionFitProg::updateFCC(
        const std::vector<std::vector<Image<Complex>>>& movie,
        const std::vector<std::vector<d2Vector>>& tracks,
        const MetaDataTable& mdt,
        std::vector<Image<RFLOAT>>& tables,
        std::vector<Image<RFLOAT>>& weights0,
        std::vector<Image<RFLOAT>>& weights1)
{
    const int pc = mdt.numberOfObjects();

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
        mdt.getValue(EMDL_PARTICLE_RANDOM_SUBSET, randSubset, p);
        randSubset -= 1;

        if (randSubset == 0)
        {
            pred = obsModel.predictObservation(projectors[1], mdt, p, true, true);
        }
        else
        {
            pred = obsModel.predictObservation(projectors[0], mdt, p, true, true);
        }

        FscHelper::updateFscTable(obs, pred, tables[threadnum],
                                  weights0[threadnum], weights1[threadnum]);
    }
}

void MotionFitProg::writeOutput(
        const std::vector<std::vector<d2Vector>>& tracks,
        const std::vector<Image<RFLOAT>>& fccData,
        const std::vector<Image<RFLOAT>>& fccWeight0,
        const std::vector<Image<RFLOAT>>& fccWeight1,
        const std::vector<d2Vector>& positions,
        std::string outPath, int mg,
        double visScale)
{
    const int pc = tracks.size();

    if (pc == 0) return;

    const int fc = tracks[0].size();

    std::string tag = getMicrographTag(mg);
    MotionRefinement::writeTracks(tracks, outPath + "/" + tag + "_tracks.star");

    Image<RFLOAT> fccDataSum(sh,fc), fccWeight0Sum(sh,fc), fccWeight1Sum(sh,fc);
    fccDataSum.data.initZeros();
    fccWeight0Sum.data.initZeros();
    fccWeight1Sum.data.initZeros();

    for (int i = 0; i < fccData.size(); i++)
    {
        for (int y = 0; y < fc; y++)
        for (int x = 0; x < sh; x++)
        {
            fccDataSum(y,x) += fccData[i](y,x);
            fccWeight0Sum(y,x) += fccWeight0[i](y,x);
            fccWeight1Sum(y,x) += fccWeight1[i](y,x);
        }
    }

    fccDataSum.write(outPath + "/" + tag + "_FCC_cc.mrc");
    fccWeight0Sum.write(outPath + "/" + tag + "_FCC_w0.mrc");
    fccWeight1Sum.write(outPath + "/" + tag + "_FCC_w1.mrc");

    if (!diag) return;

    std::stringstream sts;
    sts << mg;

    mktree(outPath + "/diag");

    std::string diagPath = outPath + "/diag/mg" + sts.str();

    // plot graphs here:

    std::vector<std::vector<gravis::d2Vector>>
            centTracks(pc), visTracks(pc), centVisTracks(pc);

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

    std::ofstream rawOut(diagPath + "_tracks.dat");
    std::ofstream visOut(diagPath + "_visTracks.dat");
    std::ofstream visOut15(diagPath + "_visTracks_first15.dat");

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

    std::ofstream glbOut(diagPath + "_globTrack.dat");

    for (int f = 0; f < fc; f++)
    {
        glbOut << globalTrack[f].x << " " << globalTrack[f].y << "\n";
    }
}


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
#include <src/jaz/d3x3/dsyev2.h>

#include <src/jaz/motion_em.h>

#include <omp.h>

using namespace gravis;

class MotionFitProg : public RefinementProgram
{
    public:

        MotionFitProg();

            bool unregGlob, noGlobOff, paramEstim, debugOpt, diag;
            int maxIters;
            double dmga, dmgb, dmgc, dosePerFrame,
                sig_vel, sig_div, sig_acc,
                k_cutoff, maxStep, maxDistDiv,
                param_rV, param_rD, param_thresh;

        int readMoreOptions(IOParser& parser, int argc, char *argv[]);
        int _init();
        int _run();

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
                const std::vector<Image<RFLOAT>>& dmgWeight,
                int k_out);

        d2Vector estimateParams(
                std::vector<ParFourierTransformer>& fts,
                const std::vector<Image<double>>& dmgWeight,
                int k_out, double sig_v_0, double sig_d_0,
                double sig_v_step, double sig_d_step,
                int maxIters);

        void evaluateParams(
                std::vector<ParFourierTransformer>& fts,
                const std::vector<Image<double>>& dmgWeight,
                int k_out,
                const std::vector<d2Vector>& sig_vals,
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
                const std::vector<double>& velWgh,
                const std::vector<double>& accWgh,
                const std::vector<std::vector<std::vector<double>>>& divWgh,
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
                const std::vector<d2Vector>& positions,
                std::string outPath, std::string mgIndex,
                double visScale);

        /*void writeParamEstOutput(
                const std::vector<std::vector<double>>& paramTsc,
                const std::vector<double>& bestV,
                const std::vector<double>& bestD,
                const std::vector<int>& cumulPartCount,
                const std::vector<std::string>& paramTags,
                std::string outPath, std::string mgIndex);*/

        /*d2Vector optimalParams(
                const std::vector<double>& paramTsc,
                const std::vector<double>& sig_v_vals,
                const std::vector<double>& sig_d_vals);*/

        d2Vector interpolateMax(
                const std::vector<d2Vector>& all_sig_vals,
                const std::vector<double>& all_TSCs);
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

    debugOpt = parser.checkOption("--debug_opt", "Write optimization debugging info");
    diag = parser.checkOption("--diag", "Write out diagnostic data");

    parser.addSection("Parameter estimation");

    paramEstim = parser.checkOption("--params", "Estimate parameters instead of motion");
    param_rV = textToFloat(parser.getOption("--r_vel", "Test s_vel +/- r_vel * s_vel", "0.5"));
    param_rD = textToFloat(parser.getOption("--r_div", "Test s_div +/- r_div * s_div", "0.2"));
    param_thresh = textToFloat(parser.getOption("--pthresh", "Abort when relative TSC change is smaller than this", "0.1"));

    if (paramEstim && k_cutoff < 0)
    {
        std::cerr << "Parameter estimation requires a freq. cutoff (--k_cut).\n";
        return 1138;
    }

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
            // @TODO: test!
        }
    }

    double t0 = omp_get_wtime();

    if (paramEstim)
    {
        d2Vector par = estimateParams(
            fts, dmgWeight, k_out, sig_vel, sig_div,
            sig_vel*param_rV, sig_div*param_rD, 10);

        std::cout << "opt = " << par << "\n";
    }
    else
    {
        estimateMotion(fts, dmgWeight, k_out);
    }


    double t1 = omp_get_wtime();
    double diff = t1 - t0;
    std::cout << "elapsed (total): " << diff << " sec\n";

    return 0;
}

void MotionFitProg::prepMicrograph(
        int g, std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<double>>& dmgWeight,
        std::vector<std::vector<Image<Complex>>>& movie,
        std::vector<std::vector<Image<double>>>& movieCC,
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

    std::cout << "    computing initial correlations...\n";

    movieCC = MotionRefinement::movieCC(
            projectors[0], projectors[1], obsModel, mdts[g], movie,
            sigma2, dmgWeight, fts, nr_omp_threads);

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
        std::stringstream stsg;
        stsg << g;

        ImageLog::write(ccSum, outPath + "_CCsum_mg"+stsg.str(), CenterXY);
    }

    initialTracks = std::vector<std::vector<d2Vector>>(pc);

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

void MotionFitProg::estimateMotion(
        std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<double>>& dmgWeight,
        int k_out)
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
    const double sig_div_nrm = dosePerFrame * sqrt(coords_angpix) * sig_div / angpix;

    int pctot = 0;

    // initialize parameter-estimation:

    for (long g = g0; g <= gc; g++)
    {
        std::cout << "micrograph " << (g+1) << " / " << mdts.size() <<"\n";

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

        const int pc = movie.size();
        pctot += pc;

        std::vector<double> velWgh, accWgh;
        std::vector<std::vector<std::vector<double>>> divWgh;

        std::cout << "    optimizing...\n";

        computeWeights(sig_vel_nrm, sig_acc_nrm, sig_div_nrm,
                       positions, fc, velWgh, accWgh, divWgh);

        std::vector<std::vector<gravis::d2Vector>> tracks = optimize(
                movieCC, initialTracks, velWgh, accWgh, divWgh, globComp,
                maxStep, 1e-9, 1e-9, maxIters, 0.0);

        updateFCC(movie, tracks, mdts[g], tables, weights0, weights1);

        writeOutput(tracks, positions, outPath, stsg.str(), 30.0);

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

}

/*void MotionFitProg::estimateParams_old(
        std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<double>>& dmgWeight,
        int k_out)
{
    const double sig_vel_nrm = dosePerFrame * sig_vel / angpix;
    const double sig_acc_nrm = dosePerFrame * sig_acc / angpix;
    const double sig_div_nrm = dosePerFrame * sqrt(coords_angpix) * sig_div / angpix;

    // initialize parameter-estimation:
    const int paramCount = 9;
    int pctot = 0;

    std::vector<double> sig_v_vals(paramCount), sig_d_vals(paramCount);
    std::vector<double> sig_v_vals_nrm(paramCount), sig_d_vals_nrm(paramCount);
    std::vector<std::string> paramTags(paramCount);

    std::vector<std::vector<double>> allParamTsc(paramCount);
    std::vector<double> paramTsc(paramCount);
    std::vector<double> bestV(0), bestD(0);

    for (int i = 0; i < paramCount; i++)
    {
        const double dv = ((i%3)-1) * param_rV;
        const double dd = ((i/3)-1) * param_rD;

        sig_v_vals[i] = sig_vel * (1.0 + dv);
        sig_d_vals[i] = sig_div * (1.0 + dd);
        sig_v_vals_nrm[i] = sig_vel_nrm * (1.0 + dv);
        sig_d_vals_nrm[i] = sig_div_nrm * (1.0 + dd);

        std::stringstream sts;
        sts << "v" << sig_v_vals[i] << "-d" << sig_d_vals[i];

        paramTags[i] = sts.str();

        allParamTsc[i] = std::vector<double>(0);
    }

    std::vector<std::vector<Image<RFLOAT>>>
        paramTables(paramCount), paramWeights0(paramCount), paramWeights1(paramCount);

    std::vector<int> cumulPartCount(0);

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

    for (long g = g0; g <= gc; g++)
    {
        std::cout << "micrograph " << (g+1) << " / " << mdts.size() <<"\n";

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

        const int pc = movie.size();
        pctot += pc;

        std::vector<double> velWgh, accWgh;
        std::vector<std::vector<std::vector<double>>> divWgh;

        std::cout << "    optimizing...\n";

        for (int i = 0; i < paramCount; i++)
        {
            if (debug)
            {
                std::cout << "        optimizing: sig_vel = " << sig_v_vals[i] << " px\n";
                std::cout << "                    sig_div = " << sig_d_vals[i] << " px^(1/2)\n";
            }

            computeWeights(sig_v_vals_nrm[i], sig_acc_nrm, sig_d_vals_nrm[i],
                           positions, fc, velWgh, accWgh, divWgh);

            std::vector<std::vector<gravis::d2Vector>> tracks = optimize(
                    movieCC, initialTracks, velWgh, accWgh, divWgh, globComp,
                    maxStep, 1e-9, 1e-9, maxIters, 0.0);

            updateFCC(movie, tracks, mdts[g], paramTables[i], paramWeights0[i], paramWeights1[i]);

            double tsc = FscHelper::computeTsc(
                paramTables[i], paramWeights0[i], paramWeights1[i], k_cutoff+2, k_out);

            paramTsc[i] = tsc;
            allParamTsc[i].push_back(tsc);
        }

        cumulPartCount.push_back(pctot);

        d2Vector optPar = optimalParams(paramTsc, sig_v_vals, sig_d_vals);

        bestV.push_back(optPar[0]);
        bestD.push_back(optPar[1]);

        writeParamEstOutput(allParamTsc, bestV, bestD,
                cumulPartCount, paramTags, outPath, stsg.str());

    } // micrographs
}
*/

d2Vector MotionFitProg::estimateParams(
        std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<double>>& dmgWeight,
        int k_out, double sig_v_0, double sig_d_0,
        double sig_v_step, double sig_d_step,
        int maxIters)
{
    std::vector<d2Vector> all_sig_vals(9);

    for (int p = 0; p < 9; p++)
    {
        int vi = (p%3) - 1;
        int di = (p/3) - 1;

        all_sig_vals[p][0] = sig_v_0 + vi * sig_v_step;
        all_sig_vals[p][1] = sig_d_0 + di * sig_d_step;
    }

    std::vector<double> all_TSCs(9);

    std::vector<d2Vector> unknown_sig_vals = all_sig_vals;
    std::vector<double> unknown_TSCs(9);

    std::vector<int> unknown_ind(9);

    for (int p = 0; p < 9; p++)
    {
        unknown_ind[p] = p;
    }

    bool centerBest = false;
    int iters = 0;

    int tot_vi = 0, tot_di = 0;

    while (!centerBest && iters < maxIters)
    {
        if (debug)
        {
            std::cout << "evaluating:\n";

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

        if (debug)
        {
            std::cout << "result:\n";

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

        // say something

        if (bestTSC == 4)
        {
            return interpolateMax(all_sig_vals, all_TSCs);
        }

        int shift_v = bestIndex % 3 - 1;
        int shift_d = bestIndex / 3 - 1;

        tot_vi += shift_v;
        tot_di += shift_d;

        std::vector<d2Vector> next_sig_vals(9, d2Vector(0,0));
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

                unknown_sig_vals.push_back(all_sig_vals[p]);
                unknown_ind.push_back(p);
            }
        }

        unknown_TSCs.resize(unknown_ind.size());

        iters++;
    }
}

void MotionFitProg::evaluateParams(
        std::vector<ParFourierTransformer>& fts,
        const std::vector<Image<double>>& dmgWeight,
        int k_out,
        const std::vector<d2Vector>& sig_vals,
        std::vector<double>& TSCs)
{
    const int paramCount = sig_vals.size();
    TSCs.resize(paramCount);

    std::vector<double> sig_v_vals_nrm(paramCount);
    std::vector<double> sig_d_vals_nrm(paramCount);

    for (int i = 0; i < paramCount; i++)
    {
        sig_v_vals_nrm[i] = dosePerFrame * sig_vals[i][0] / angpix;
        sig_d_vals_nrm[i] = dosePerFrame * sqrt(coords_angpix) * sig_vals[i][1] / angpix;
    }

    double sig_acc_nrm = dosePerFrame * sig_acc / angpix;

    std::vector<std::vector<Image<RFLOAT>>>
        paramTables(paramCount), paramWeights0(paramCount), paramWeights1(paramCount);

    int pctot = 0;

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

    for (long g = g0; g <= gc; g++)
    {
        std::cout << "    micrograph " << (g+1) << " / " << mdts.size() << ": ";

        std::cout.flush();

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

        const int pc = movie.size();
        pctot += pc;

        std::cout << pc << " particles \t [" << pctot << "]\n";

        std::vector<double> velWgh, accWgh;
        std::vector<std::vector<std::vector<double>>> divWgh;

        for (int i = 0; i < paramCount; i++)
        {
            computeWeights(sig_v_vals_nrm[i], sig_acc_nrm, sig_d_vals_nrm[i],
                           positions, fc, velWgh, accWgh, divWgh);

            std::vector<std::vector<gravis::d2Vector>> tracks = optimize(
                    movieCC, initialTracks, velWgh, accWgh, divWgh, globComp,
                    maxStep, 1e-9, 1e-9, maxIters, 0.0);

            updateFCC(movie, tracks, mdts[g], paramTables[i], paramWeights0[i], paramWeights1[i]);
        }

    } // micrographs

    for (int i = 0; i < paramCount; i++)
    {
        TSCs[i] = FscHelper::computeTsc(
            paramTables[i], paramWeights0[i], paramWeights1[i], k_cutoff+2, k_out);
    }
}

void MotionFitProg::computeWeights(
        double sig_vel_nrm, double sig_acc_nrm, double sig_div_nrm,
        const std::vector<d2Vector> &positions, int fc,
        std::vector<double> &velWgh,
        std::vector<double> &accWgh,
        std::vector<std::vector<std::vector<double> > > &divWgh)
{
    const int pc = positions.size();

    velWgh = std::vector<double>(fc-1, 0.5/(sig_vel_nrm*sig_vel_nrm));
    accWgh = std::vector<double>(fc-1, sig_acc > 0.0? 0.5/(sig_acc_nrm*sig_acc_nrm) : 0.0);
    divWgh = std::vector<std::vector<std::vector<double>>>(fc-1);

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
}

std::vector<std::vector<d2Vector>> MotionFitProg::optimize(
        const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
        const std::vector<std::vector<d2Vector>>& inTracks,
        const std::vector<double>& velWgh,
        const std::vector<double>& accWgh,
        const std::vector<std::vector<std::vector<double>>>& divWgh,
        const std::vector<d2Vector>& globComp,
        double step, double minStep, double minDiff,
        long maxIters, double inertia)
{
    const int pc = inTracks.size();

    if (pc == 0) return std::vector<std::vector<d2Vector>>(0);

    const int fc = inTracks[0].size();

    LocalMotionFit lmf(movieCC, velWgh, accWgh, divWgh, globComp, nr_omp_threads);

    std::vector<double> initial(2*pc*fc);

    for (int p = 0; p < pc; p++)
    {
        for (int f = 0; f < fc; f++)
        {
            initial[2*(p*fc + f)]     = inTracks[p][f].x;
            initial[2*(p*fc + f) + 1] = inTracks[p][f].y;
        }
    }

    std::vector<double> optPos = GradientDescent::optimize(
        initial, lmf, step, minStep, minDiff, maxIters, inertia, debugOpt);

    std::vector<std::vector<d2Vector>> out(pc);

    for (int p = 0; p < pc; p++)
    {
        out[p] = std::vector<d2Vector>(fc);

        for (int f = 0; f < fc; f++)
        {
            out[p][f].x = optPos[2*(p*fc + f)] + globComp[f].x;
            out[p][f].y = optPos[2*(p*fc + f) + 1] + globComp[f].y;
        }
    }

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
        const std::vector<d2Vector>& positions,
        std::string outPath, std::string mgIndex,
        double visScale)
{
    const int pc = tracks.size();

    if (pc == 0) return;

    const int fc = tracks[0].size();

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

    std::ofstream rawOut(outPath + "_mg" + mgIndex + "_tracks.dat");
    std::ofstream visOut(outPath + "_mg" + mgIndex + "_visTracks.dat");
    std::ofstream visOut15(outPath + "_mg" + mgIndex + "_visTracks_first15.dat");

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

    std::ofstream glbOut(outPath + "_mg" + mgIndex + "_globTrack.dat");

    for (int f = 0; f < fc; f++)
    {
        glbOut << globalTrack[f].x << " " << globalTrack[f].y << "\n";
    }
}

d2Vector MotionFitProg::interpolateMax(
    const std::vector<d2Vector> &all_sig_vals,
    const std::vector<double> &all_TSCs)
{
    const int parCt = all_sig_vals.size();

    Matrix2D<RFLOAT> A(parCt,6);
    Matrix1D<RFLOAT> b(parCt);

    int bestP = 0;
    double bestTsc = all_TSCs[0];

    for (int p = 0; p < parCt; p++)
    {
        if (all_TSCs[p] > bestTsc)
        {
            bestTsc = all_TSCs[p];
            bestP = p;
        }
    }

    if (bestP != 4)
    {
        std::cerr << "Waring: value not maximal at the center.\n";
        return all_sig_vals[bestP];
    }

    for (int p = 0; p < parCt; p++)
    {
        const double v = all_sig_vals[p][0];
        const double d = all_sig_vals[p][1];

        A(p,0) = v*v;
        A(p,1) = 2.0*v*d;
        A(p,2) = 2.0*v;
        A(p,3) = d*d;
        A(p,4) = 2.0*d;
        A(p,5) = 1.0;

        b(p) = all_TSCs[p];
    }

    const double tol = 1e-20;
    Matrix1D<RFLOAT> x(6);
    solve(A, b, x, tol);

    d2Matrix C2(x(0), x(1),
                x(1), x(3));

    d2Vector l(x(2), x(4));

    d2Matrix C2i = C2;
    C2i.invert();

    d2Vector min = -(C2i * l);

    return min;
}

/*void MotionFitProg::writeParamEstOutput(
    const std::vector<std::vector<double>>& allParamTsc,
    const std::vector<double>& bestV,
    const std::vector<double>& bestD,
    const std::vector<int>& cumulPartCount,
    const std::vector<std::string>& paramTags,
    std::string outPath, std::string mgIndex)
{
    const int parCt = allParamTsc.size();
    const int sc = allParamTsc[0].size();

    for (int p = 0; p < parCt; p++)
    {
        std::ofstream outTsc(outPath + "_mg" + mgIndex + "_TSC_" + paramTags[p] + ".dat");
        std::ofstream outTscRel(outPath + "_mg" + mgIndex + "_TSC-rel_" + paramTags[p] + ".dat");

        for (int s = 0; s < sc; s++)
        {
            outTsc << cumulPartCount[s] << " " << allParamTsc[p][s] << "\n";

            if (allParamTsc[4][s] > 0.0)
            {
                outTscRel << cumulPartCount[s] << " " << allParamTsc[p][s]/allParamTsc[4][s] << "\n";
            }
        }
    }

    std::ofstream outV(outPath + "_mg" + mgIndex + "_best_sigma_v.dat");
    std::ofstream outD(outPath + "_mg" + mgIndex + "_best_sigma_d.dat");

    for (int s = 0; s < sc; s++)
    {
        outV << cumulPartCount[s] << " " << bestV[s] << "\n";
        outD << cumulPartCount[s] << " " << bestD[s] << "\n";
    }
}*/

/*d2Vector MotionFitProg::optimalParams(
    const std::vector<double>& paramTsc,
    const std::vector<double>& sig_v_vals,
    const std::vector<double>& sig_d_vals)
{
    const int parCt = paramTsc.size();

    Matrix2D<RFLOAT> A(parCt,6);
    Matrix1D<RFLOAT> b(parCt);

    int bestP = 0;
    double bestTsc = paramTsc[0];

    for (int p = 0; p < parCt; p++)
    {
        if (paramTsc[p] > bestTsc)
        {
            bestTsc = paramTsc[p];
            bestP = p;
        }
    }

    if (bestP != 4)
    {

        return d2Vector(sig_v_vals[bestP], sig_d_vals[bestP]);
    }

    for (int p = 0; p < parCt; p++)
    {
        const double v = sig_v_vals[p];
        const double d = sig_d_vals[p];

        A(p,0) = v*v;
        A(p,1) = 2.0*v*d;
        A(p,2) = 2.0*v;
        A(p,3) = d*d;
        A(p,4) = 2.0*d;
        A(p,5) = 1.0;

        b(p) = paramTsc[p];
    }

    const double tol = 1e-20;
    Matrix1D<RFLOAT> x(6);
    solve(A, b, x, tol);

    d2Matrix C2(x(0), x(1),
                x(1), x(3));

    d2Vector l(x(2), x(4));

    double Cd[2][2];

    Cd[0][0] = (double)x(0);
    Cd[0][1] = (double)x(1);
    Cd[1][0] = (double)x(1);
    Cd[1][1] = (double)x(3);

    double eigVal0, eigVal1;
    d2Vector eigVec0;

    dsyev2(x(0), x(1), x(3), &eigVal0, &eigVal1, &eigVec0.x, &eigVec0.y);

    d2Vector eigVec1(eigVec0.y, -eigVec0.x);

    if (eigVal0 < 0.0 && eigVal1 < 0.0)
    {
        if (debug)
        {
            std::cout << "negative curvature.\n";
        }

        d2Matrix C2i = C2;
        C2i.invert();

        d2Vector min = -(C2i * l);

        return min;
    }
    else if (eigVal0 > 0.0 && eigVal1 <= 0.0)
    {
        if (debug)
        {
            std::cout << "mixed curvature.\n";
        }

        double m0 = l.dot(eigVec0);
    }
}*/

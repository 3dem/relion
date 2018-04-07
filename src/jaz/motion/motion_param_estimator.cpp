#include "motion_param_estimator.h"
#include "motion_refiner.h"
#include "two_hyperparameter_fit.h"
#include "three_hyperparameter_fit.h"

#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/index_sort.h>
#include <src/jaz/filter_helper.h>


using namespace gravis;

const double MotionParamEstimator::velScale = 1000.0;
const double MotionParamEstimator::divScale = 1.0;
const double MotionParamEstimator::accScale = 1000.0;

MotionParamEstimator::MotionParamEstimator()
:   paramsRead(false), ready(false)
{
}

int MotionParamEstimator::read(IOParser& parser, int argc, char *argv[])
{
    parser.addSection("Parameter estimation");

    estim2 = parser.checkOption("--params2", "Estimate 2 parameters instead of motion");
    estim3 = parser.checkOption("--params3", "Estimate 3 parameters instead of motion");
    k_cutoff = textToFloat(parser.getOption("--k_cut", "Freq. cutoff for parameter estimation [Pixels]", "-1.0"));
    k_cutoff_Angst = textToFloat(parser.getOption("--k_cut_A", "Freq. cutoff for parameter estimation [Angstrom]", "-1.0"));
    k_eval = textToFloat(parser.getOption("--k_eval", "Threshold freq. for parameter evaluation [Pixels]", "-1.0"));
    k_eval_Angst = textToFloat(parser.getOption("--k_eval_A", "Threshold freq. for parameter evaluation [Angstrom]", "-1.0"));

    minParticles = textToInteger(parser.getOption("--min_p", "Minimum number of particles on which to estimate the parameters", "1000"));
    sV = textToFloat(parser.getOption("--s_vel_0", "Initial s_vel", "0.6"));
    sD = textToFloat(parser.getOption("--s_div_0", "Initial s_div", "3000"));
    sA = textToFloat(parser.getOption("--s_acc_0", "Initial s_acc", "-1"));
    iniStep = textToFloat(parser.getOption("--in_step", "Initial step size in s_div", "500"));
    conv = textToFloat(parser.getOption("--conv", "Abort when simplex diameter falls below this", "10"));
    maxIters = textToInteger(parser.getOption("--par_iters", "Max. number of iterations", "300"));
    maxRange = textToInteger(parser.getOption("--mot_range", "Limit allowed motion range [Px]", "50"));
    seed = textToInteger(parser.getOption("--seed", "Random seed for micrograph selection", "23"));

    paramsRead = true;
}

void MotionParamEstimator::init(
    int verb, int nr_omp_threads, bool debug,
    int s, int fc,
    const std::vector<MetaDataTable>& allMdts,
    MotionEstimator* motionEstimator,
    ReferenceMap* reference,
    ObservationModel* obsModel)
{
    if (!paramsRead)
    {
        REPORT_ERROR("ERROR: MotionParamEstimator::init: MotionParamEstimator has not read its cmd-line parameters.");
    }

    this->verb = verb;
    this->nr_omp_threads = nr_omp_threads;
    this->debug = debug;
    this->s = s;
    this->fc = fc;
    this->motionEstimator = motionEstimator;
    this->obsModel = obsModel;
    this->reference = reference;

    if (!motionEstimator->isReady())
    {
        REPORT_ERROR("ERROR: MotionParamEstimator initialized before MotionEstimator.");
    }

    if (k_cutoff_Angst > 0.0 && k_cutoff > 0.0)
    {
        REPORT_ERROR("ERROR: Cutoff frequency can only be provided in pixels (--k_cut) or Angstrom (--k_eval_A), not both.");
    }

    if (k_eval_Angst > 0.0 && k_eval > 0.0)
    {
        REPORT_ERROR("ERROR: Evaluation frequency can only be provided in pixels (--k_eval) or Angstrom (--k_cut_A), not both.");
    }

    if (k_cutoff_Angst > 0.0 && k_cutoff < 0)
    {
        k_cutoff = obsModel->angToPix(k_cutoff_Angst, s);
    }

    else if (k_cutoff > 0 && k_cutoff_Angst < 0.0)
    {
        k_cutoff_Angst = obsModel->angToPix(k_cutoff, s);
    }

    if ((estim2 || estim3) && k_cutoff < 0)
    {
        REPORT_ERROR("ERROR: Parameter estimation requires a freq. cutoff (--k_cut or --k_cut_A).");
    }

    if (estim2 && estim3)
    {
        REPORT_ERROR("ERROR: Only 2 or 3 parameters can be estimated (--params2 or --params3), not both.");
    }

    if (k_eval < 0 && k_eval_Angst > 0.0)
    {
        k_eval = obsModel->angToPix(k_eval_Angst, s);
    }
    else if (k_eval > 0 && k_eval_Angst < 0.0)
    {
        k_eval_Angst = obsModel->angToPix(k_eval, s);
    }
    else
    {
        k_eval = k_cutoff;
        k_eval_Angst = k_cutoff_Angst;
    }

    if (verb > 0)
    {
        std::cout << " + maximum frequency to consider for alignment: "
            << k_cutoff_Angst << " A (" << k_cutoff << " px)\n";

        std::cout << " + frequency range to consider for evaluation:  "
            << k_eval_Angst << " - " << obsModel->pixToAng(reference->k_out,s) << " A ("
            << k_eval << " - " << reference->k_out << " px)\n";

    }

    const long mc = allMdts.size();

    srand(seed);

    std::vector<double> randNums(mc);

    for (int m = 0; m < mc; m++)
    {
        randNums[m] = rand() / (double)RAND_MAX;
    }

    std::vector<int> order = IndexSort<RFLOAT>::sortIndices(randNums);

    int pc = 0;
    mdts.clear();

    if (verb > 0)
    {
        std::cout << " + micrographs randomly selected for parameter optimization:\n";
    }

    for (int i = 0; i < order.size(); i++)
    {
        const int m = order[i];

        const int pcm = allMdts[m].numberOfObjects();

        // motion estimation does not work on one single particle
        if (pcm < 2) continue;

        mdts.push_back(allMdts[m]);
        pc += pcm;

        if (verb > 0)
        {
            std::string mn;
            allMdts[m].getValue(EMDL_MICROGRAPH_NAME, mn, 0);

            std::cout << "        " << m << ": " << mn << "\n";
        }

        if (pc >= minParticles)
        {
            if (verb > 0)
            {
                std::cout << "\n + " << pc << " particles found in "
                          << mdts.size() << " micrographs\n";
            }

            break;
        }
    }

    if (verb > 0 && pc < minParticles)
    {
        std::cout << "\n   - Warning: this dataset does not contain " << minParticles
                  << " particles (--min_p) in micrographs with at least 2 particles\n";
    }

    k_out = reference->k_out;

    ready = true;
}

void MotionParamEstimator::run()
{
    #ifdef TIMING
        timeSetup = paramTimer.setNew(" time_Setup ");
        timeOpt = paramTimer.setNew(" time_Opt ");
        timeEval = paramTimer.setNew(" time_Eval ");
    #endif

    if (!ready)
    {
        REPORT_ERROR("ERROR: MotionParamEstimator::run: MotionParamEstimator not initialized.");
    }

    if (!estim2 && !estim3) return;

    RCTIC(paramTimer, timeSetup);

    prepAlignment();

    RCTOC(paramTimer, timeSetup);

    d4Vector opt;

    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(6);

    if (estim2)
    {
        //estimateTwoParamsRec(sV, sD, sA, rV*sV, rD*sD, steps, recursions);
        opt = estimateTwoParamsNM(sV, sD, sA, iniStep, conv, maxIters);
    }

    if (estim3)
    {
        //estimateThreeParamsRec(sV, sD, sA, rV*sV, rD*sD, rA*sA, steps, recursions);
        opt = estimateThreeParamsNM(sV, sD, sA, iniStep, conv, maxIters);
    }

    std::cout.setf(std::ios::floatfield);

    d3Vector nrm(
        opt[0] * velScale,
        opt[1] * divScale,
        opt[2] * accScale);

    // round result to conv / 2 (the min. radius of the optimization simplex)
    d3Vector rnd(
        conv * 0.5 * ((int)(2.0*nrm[0]/conv + 0.5)) / velScale,
        conv * 0.5 * ((int)(2.0*nrm[1]/conv + 0.5)) / divScale,
        conv * 0.5 * ((int)(2.0*nrm[2]/conv + 0.5)) / accScale);

    if (estim2)
    {
        rnd[2] = sA;
    }

    if (opt[2] <= 0.0)
    {
        rnd[2] = -1;
    }

    std::cout << "\ngood parameters:"
              << " --s_vel " << rnd[0]
              << " --s_div " << rnd[1]
              << " --s_acc " << rnd[2] << "\n\n";

    #ifdef TIMING
        paramTimer.printTimes(true);
    #endif
}

bool MotionParamEstimator::anythingToDo()
{
    return estim2 || estim3;
}

d4Vector MotionParamEstimator::estimateTwoParamsRec(
        double sig_v_0, double sig_d_0, double sig_acc,
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
    double bestTSC = 0.0;

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

        evaluateParams(unknown_sig_vals, unknown_TSCs);

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
        bestTSC = all_TSCs[0];

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
            std::cout << "(" << all_TSCs[bestIndex] << ")\n";
        }

        if (shift_v < 0 && all_sig_vals[3][0] <= 0.0)
        {
            shift_v = 0;
        }

        if (shift_d < 0 && all_sig_vals[1][1] <= 0.0)
        {
            shift_d = 0;
        }

        if (shift_v == 0 && shift_d == 0)
        {
            if (recDepth <= 0)
            {
                return d4Vector(bestParams[0],
                                bestParams[1],
                                bestParams[2],
                                bestTSC);
            }
            else
            {
                std::cout << "\nrepeating at half scan range.\n\n";

                return estimateTwoParamsRec(
                    all_sig_vals[bestIndex][0], all_sig_vals[bestIndex][1], sig_acc,
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

    std::cout << "\nsearch aborted after " << maxIters << " steps\n";

    return d4Vector(bestParams[0], bestParams[1], sig_acc, bestTSC);
}

d4Vector MotionParamEstimator::estimateThreeParamsRec(
    double sig_v_0, double sig_d_0, double sig_a_0,
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
    double bestTSC = 0.0;

    const bool verbose = true;

    while (!centerBest && iters < maxIters)
    {
        if (verbose)
        {
            std::cout << "\nevaluating:\n";

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

        evaluateParams(unknown_sig_vals, unknown_TSCs);

        for (int p = 0; p < unknown_ind.size(); p++)
        {
            all_TSCs[unknown_ind[p]] = unknown_TSCs[p];
        }

        if (verbose)
        {
            std::cout << "\nresult:\n";

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
        bestTSC = all_TSCs[0];

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

        if (verbose)
        {
            std::cout << "optimum: " << all_sig_vals[bestIndex] << " ";
            std::cout << "(" << all_TSCs[bestIndex] << ")\n";
        }

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
                return d4Vector(bestParams[0],
                                bestParams[1],
                                bestParams[2],
                                bestTSC);
            }
            else
            {
                std::cout << "\nrepeating at half scan range.\n\n";

                return estimateThreeParamsRec(
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

    std::cout << "\nsearch aborted after " << maxIters << " steps\n";

    return d4Vector(bestParams[0],
                    bestParams[1],
                    bestParams[2],
                    bestTSC);
}

d4Vector MotionParamEstimator::estimateTwoParamsNM(
        double sig_v_0, double sig_d_0, double sig_acc, double inStep, double conv, int maxIters)
{
    std::cout << "\nit: \t s_vel: \t s_div: \t s_acc: \t fsc:\n\n";

    TwoHyperParameterProblem thpp(*this, sig_acc);

    std::vector<double> initial = TwoHyperParameterProblem::motionToProblem(
            d2Vector(sig_v_0, sig_d_0));

    double minTsc;

    std::vector<double> final = NelderMead::optimize(
            initial, thpp, inStep, conv, maxIters,
            1.0, 2.0, 0.5, 0.5, false, &minTsc);

    d2Vector vd = TwoHyperParameterProblem::problemToMotion(final);

    return d4Vector(vd[0], vd[1], sig_acc, -minTsc);
}

d4Vector MotionParamEstimator::estimateThreeParamsNM(
        double sig_v_0, double sig_d_0, double sig_a_0, double inStep, double conv, int maxIters)
{
    std::cout << "\nit: \t s_vel: \t s_div: \t s_acc: \t fsc:\n\n";

    ThreeHyperParameterProblem thpp(*this);

    std::vector<double> initial = ThreeHyperParameterProblem::motionToProblem(
            d3Vector(sig_v_0, sig_d_0, sig_a_0));

    double minTsc;

    std::vector<double> final = NelderMead::optimize(
            initial, thpp, inStep, conv, maxIters,
            1.0, 2.0, 0.5, 0.5, false, &minTsc);

    d3Vector vd = ThreeHyperParameterProblem::problemToMotion(final);

    return d4Vector(vd[0], vd[1], vd[2], -minTsc);
}

void MotionParamEstimator::evaluateParams(
    const std::vector<d3Vector>& sig_vals,
    std::vector<double>& TSCs)
{
    const int paramCount = sig_vals.size();
    TSCs.resize(paramCount);

    std::vector<double> sig_v_vals_px(paramCount);
    std::vector<double> sig_d_vals_px(paramCount);
    std::vector<double> sig_a_vals_px(paramCount);

    for (int i = 0; i < paramCount; i++)
    {
        sig_v_vals_px[i] = motionEstimator->normalizeSigVel(sig_vals[i][0]);
        sig_d_vals_px[i] = motionEstimator->normalizeSigDiv(sig_vals[i][1]);
        sig_a_vals_px[i] = motionEstimator->normalizeSigAcc(sig_vals[i][2]);
    }

    int pctot = 0;
    std::vector<d3Vector> tscsAs(paramCount, d3Vector(0.0, 0.0, 0.0));

    const int gc = mdts.size();

    for (long g = 0; g < gc; g++)
    {
        const int pc = mdts[g].numberOfObjects();

        if (pc < 2) continue; // not really needed, mdts are pre-screened

        pctot += pc;

        if (debug)
        {
            std::cout << "    micrograph " << (g+1) << " / " << mdts.size() << ": "
                << pc << " particles [" << pctot << " total]\n";
        }

        for (int i = 0; i < paramCount; i++)
        {
            if (debug)
            {
                std::cout << "        evaluating: " << sig_vals[i] << "\n";
            }

            RCTIC(paramTimer,timeOpt);

            std::vector<std::vector<gravis::d2Vector>> tracks =
                motionEstimator->optimize(
                    alignmentSet.CCs[g],
                    alignmentSet.initialTracks[g],
                    sig_v_vals_px[i], sig_a_vals_px[i], sig_d_vals_px[i],
                    alignmentSet.positions[g], alignmentSet.globComp[g]);

            RCTOC(paramTimer,timeOpt);

            RCTIC(paramTimer,timeEval);

            tscsAs[i] += alignmentSet.updateTsc(tracks, g, nr_omp_threads);

            RCTOC(paramTimer,timeEval);
        }

    } // micrographs

    if (debug)
    {
        std::cout << "\n";
    }

    RCTIC(paramTimer,timeEval);

    // compute final TSC
    for (int i = 0; i < paramCount; i++)
    {
        double wg = tscsAs[i][1] * tscsAs[i][2];

        if (wg > 0.0)
        {
            TSCs[i] = tscsAs[i][0] / sqrt(wg);
        }
    }

    RCTOC(paramTimer,timeEval);
}

void MotionParamEstimator::prepAlignment()
{
    std::cout << " + preparing alignment data... \n";

    const std::vector<Image<RFLOAT>>& dmgWgh = motionEstimator->getDamageWeights();
    std::vector<Image<RFLOAT>> alignDmgWgh(fc);

    for (int f = 0; f < fc; f++)
    {
        alignDmgWgh[f] = FilterHelper::ButterworthEnvFreq2D(dmgWgh[f], k_cutoff-1, k_cutoff+1);
    }

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    alignmentSet = AlignmentSet(mdts, fc, s, k_eval+2, k_out);

    for (int f = 0; f < fc; f++)
    {
        alignmentSet.damage[f] = alignmentSet.accelerate(dmgWgh[f]);
    }

    const int gc = mdts.size();

    int pctot = 0;

    for (long g = 0; g < gc; g++)
    {
        const int pc = mdts[g].numberOfObjects();

        if (pc < 2) continue;

        pctot += pc;

        std::cout << "        micrograph " << (g+1) << " / " << gc << ": "
            << pc << " particles [" << pctot << " total]\n";

        std::vector<std::vector<Image<Complex>>> movie;
        std::vector<std::vector<Image<RFLOAT>>> movieCC;

        try
        {
            motionEstimator->prepMicrograph(
                mdts[g], fts, alignDmgWgh,
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

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
            for (int f = 0; f < fc; f++)
            {
                if (maxRange > 0)
                {
                    movieCC[p][f] = FilterHelper::cropCorner2D(movieCC[p][f], 2*maxRange, 2*maxRange);
                }

                alignmentSet.CCs[g][p][f] = movieCC[p][f];
                alignmentSet.obs[g][p][f] = alignmentSet.accelerate(movie[p][f]);

                Image<Complex> pred = reference->predict(
                    mdts[g], p, *obsModel, ReferenceMap::Opposite);

                alignmentSet.pred[g][p] = alignmentSet.accelerate(pred);
            }
        }

        alignmentSet.initialTracks[g] = motionEstimator->optimize(
                alignmentSet.CCs[g],
                alignmentSet.initialTracks[g],
                motionEstimator->normalizeSigVel(sV),
                motionEstimator->normalizeSigAcc(sA),
                motionEstimator->normalizeSigDiv(sD),
                alignmentSet.positions[g],
                alignmentSet.globComp[g]);
    }

    std::cout << "   done\n";
}


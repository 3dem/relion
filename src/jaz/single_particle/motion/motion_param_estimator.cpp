/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#include "motion_param_estimator.h"
#include "motion_refiner.h"
#include "two_hyperparameter_fit.h"
#include "three_hyperparameter_fit.h"

#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/single_particle/index_sort.h>
#include <src/jaz/single_particle/img_proc/filter_helper.h>
#include <src/jaz/single_particle/vtk_helper.h>

#include <src/jaz/util/zio.h>

using namespace gravis;


const double MotionParamEstimator::velScale = 10000.0;
const double MotionParamEstimator::divScale = 1.0;
const double MotionParamEstimator::accScale = 1000.0;

MotionParamEstimator::MotionParamEstimator()
:   paramsRead(false), ready(false)
{
}

void MotionParamEstimator::read(IOParser& parser, int argc, char *argv[])
{
    parser.addSection("Parameter estimation");

    estim2 = parser.checkOption("--params2", "Estimate 2 parameters instead of motion");
    estim3 = parser.checkOption("--params3", "Estimate 3 parameters instead of motion");

    align_frac = textToDouble(parser.getOption("--align_frac", "Fraction of pixels to be used for alignment", "0.5"));
    eval_frac = textToDouble(parser.getOption("--eval_frac", "Fraction of pixels to be used for evaluation", "0.5"));

    minParticles = textToInteger(parser.getOption("--min_p", "Minimum number of particles on which to estimate the parameters", "1000"));
    group = textToInteger(parser.getOption("--par_group", "Estimate parameters for this optics group only (negative means all)", "-1")) - 1;
    sV = textToDouble(parser.getOption("--s_vel_0", "Initial s_vel", "0.6"));
    sD = textToDouble(parser.getOption("--s_div_0", "Initial s_div", "10000"));
    sA = textToDouble(parser.getOption("--s_acc_0", "Initial s_acc", "3"));
    iniStep = textToDouble(parser.getOption("--in_step", "Initial step size in s_div", "3000"));
    conv = textToDouble(parser.getOption("--conv", "Abort when simplex diameter falls below this", "30"));
    maxIters = textToInteger(parser.getOption("--par_iters", "Max. number of iterations", "100"));
    maxRange = textToInteger(parser.getOption("--mot_range", "Limit allowed motion range [Px]", "50"));
    seed = textToInteger(parser.getOption("--seed", "Random seed for micrograph selection", "23"));

    paramsRead = true;
}

void MotionParamEstimator::init(
    int verb, int nr_omp_threads, bool debug,
    std::string outPath, int fc,
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
	this->outPath = outPath;
    this->fc = fc;
    this->motionEstimator = motionEstimator;
    this->obsModel = obsModel;
    this->reference = reference;

	this->s_ref = reference->s;

	std::vector<int> allS, allSh;
	obsModel->getBoxSizes(allS, allSh);

	if (group < 0)
	{
		if (!obsModel->allPixelSizesIdentical()
		 || !obsModel->allBoxSizesIdentical())
		{
			REPORT_ERROR_STR("MotionParamEstimator::init: unable to estimate motion parameters for all \n"
						 << "optics groups simultaneously due to varying pixel and box sizes.\n"
						 << "Please estimate them separately for each optics group (--par_group).");
		}

		s = allS[0];
		sh = allSh[0];

		group = 0;

		allGroups = true;

		if (verb > 0)
		{
			std::cout << " + estimating motion parameters for all optics groups simultaneously ..."
					  << std::endl;
		}
	}
	else
	{
		s = allS[group];
		sh = allSh[group];

		allGroups = false;

		if (verb > 0)
		{
			std::cout << " + estimating motion parameters for optics group "
					  << obsModel->getGroupName(group) << " ..." << std::endl;
		}
	}

    if (!motionEstimator->isReady())
    {
        REPORT_ERROR("ERROR: MotionParamEstimator initialized before MotionEstimator.");
    }

    if (eval_frac + align_frac > 1.000001)
    {
        REPORT_ERROR_STR("ERROR: Alignment and evaluation sets are intersecting. "
			<< "Please make sure that --align_frac and --eval_frac do not add up to more than 1.");
    }

    if (estim2 && estim3)
    {
        REPORT_ERROR("ERROR: Only 2 or 3 parameters can be estimated (--params2 or --params3), not both.");
    }

	k_out = reference->k_out;

	k_cutoff = (int)(k_out * sqrt(align_frac) + 0.5);
	k_eval = (int)(k_out * sqrt(1.0 - eval_frac) + 0.5);

    if (verb > 0)
    {
		double k_cutoff_Angst = reference->angToPix(k_cutoff);
		double k_eval_Angst = reference->angToPix(k_eval);

        std::cout << " + maximum frequency to consider for alignment: "
            << k_cutoff_Angst << " A (" << k_cutoff << " ref. px)" << std::endl;

        std::cout << " + frequency range to consider for evaluation:  "
            << k_eval_Angst << " - " << reference->pixToAng(reference->k_out) << " A ("
            << k_eval << " - " << reference->k_out << " ref. px)" << std::endl;
    }

    const long mc = allMdts.size();

    srand(seed);

    std::vector<double> randNums(mc);

    for (int m = 0; m < mc; m++)
    {
        randNums[m] = rand() / (double)RAND_MAX;
    }

    std::vector<int> order = IndexSort<double>::sortIndices(randNums);

    int pc = 0;
    mdts.clear();

    if (verb > 0)
    {
        std::cout << " + micrographs randomly selected for parameter optimization:" << std::endl;
    }

    for (int i = 0; i < order.size(); i++)
    {
        const int m = order[i];

        const int pcm = allMdts[m].numberOfObjects();

        // motion estimation does not work on one single particle
        if (pcm < 2) continue;

		if (allGroups)
		{
			mdts.push_back(allMdts[m]);
			pc += pcm;
		}
		else
		{
			MetaDataTable rightGroup;

			for (int p = 0; p < pcm; p++)
			{
				if (obsModel->getOpticsGroup(allMdts[m], p) == group)
				{
					rightGroup.addObject(allMdts[m].getObject(p));
				}
			}

			mdts.push_back(rightGroup);
			pc += rightGroup.numberOfObjects();
		}

        if (verb > 0)
        {
            std::string mn;
            allMdts[m].getValue(EMDL_MICROGRAPH_NAME, mn, 0);

            std::cout << "        " << m << ": " << mn << std::endl;
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
	std::cout.precision(5);

    if (estim2)
    {
        opt = estimateTwoParamsNM(sV, sD, sA, iniStep, conv, maxIters);
    }

    if (estim3)
    {
        opt = estimateThreeParamsNM(sV, sD, sA, iniStep, conv, maxIters);
	}

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

	FileName newdir = FileName(outPath).beforeLastOf("/");

	ZIO::makeDir(newdir);

	std::string paramFn;

	if (allGroups)
	{
		paramFn = "opt_params_all_groups.txt";
	}
	else
	{
		paramFn = "opt_params_group_" + obsModel->getGroupName(group) + ".txt";
	}

	std::ofstream ofs(outPath+paramFn);
	ofs << rnd[0] << " ";
	ofs << rnd[1] << " ";
	ofs << rnd[2] << std::endl;
	ofs.close();

	std::cout << "written to " << (outPath+paramFn) << std::endl;

    #ifdef TIMING
        paramTimer.printTimes(true);
    #endif
}

bool MotionParamEstimator::anythingToDo()
{
    return estim2 || estim3;
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
			1.0, 2.0, 0.5, 0.5, true, &minTsc);

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
			1.0, 2.0, 0.5, 0.5, true, &minTsc);

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
        sig_v_vals_px[i] = motionEstimator->normalizeSigVel(sig_vals[i][0], reference->angpix);
        sig_d_vals_px[i] = motionEstimator->normalizeSigDiv(sig_vals[i][1], reference->angpix);
        sig_a_vals_px[i] = motionEstimator->normalizeSigAcc(sig_vals[i][2], reference->angpix);
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
                << pc << " particles [" << pctot << " total]" << std::endl;
        }

        for (int i = 0; i < paramCount; i++)
        {
            if (debug)
            {
                std::cout << "        evaluating: " << sig_vals[i] << std::endl;
            }

            RCTIC(paramTimer,timeOpt);

            std::vector<std::vector<gravis::d2Vector>> tracks =
                motionEstimator->optimize(
                    alignmentSet.CCs[g],
                    alignmentSet.initialTracks[g],
                    sig_v_vals_px[i], sig_a_vals_px[i], sig_d_vals_px[i],
                    alignmentSet.positions[g], alignmentSet.globComp[g]);

			if (debug)
            {
				std::stringstream sts;
				sts << "debug-track_" << sig_vals[i][0] << "_" << sig_vals[i][1] << "_" << sig_vals[i][2] << ".dat";

				std::ofstream debugStr(sts.str());

				for (int p = 0; p < pc; p++)
				{
					for (int f = 0; f < fc; f++)
					{
						debugStr << tracks[p][f] << std::endl;
					}

					debugStr << std::endl;
				}

				debugStr.close();
			}

            RCTOC(paramTimer,timeOpt);

            RCTIC(paramTimer,timeEval);

            tscsAs[i] += alignmentSet.updateTsc(tracks, g, nr_omp_threads);

            RCTOC(paramTimer,timeEval);
        }

    } // micrographs

    if (debug)
    {
        std::cout << std::endl;
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
    std::cout << " + preparing alignment data... " << std::endl;

    const std::vector<Image<RFLOAT>>& dmgWgh = motionEstimator->computeDamageWeights(group);
    std::vector<Image<RFLOAT>> alignDmgWgh(fc);

    for (int f = 0; f < fc; f++)
    {
        alignDmgWgh[f] = FilterHelper::raisedCosEnvFreq2D(dmgWgh[f], k_cutoff-1, k_cutoff+1);
    }

    alignmentSet = AlignmentSet<float>(mdts, fc, s, k_eval+2, k_out);

    for (int f = 0; f < fc; f++)
    {
        alignmentSet.accelerate(dmgWgh[f], alignmentSet.damage[f]);
    }

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    const int gc = mdts.size();

    int pctot = 0;

    for (long g = 0; g < gc; g++)
    {
        const int pc = mdts[g].numberOfObjects();

        if (pc < 2) continue;

        pctot += pc;

        std::cout << "        micrograph " << (g+1) << " / " << gc << ": "
            << pc << " particles [" << pctot << " total]" << std::endl;

        std::vector<std::vector<Image<Complex>>> movie;
        std::vector<std::vector<Image<RFLOAT>>> movieCC;

        try
        {
            motionEstimator->prepMicrograph(
                mdts[g], fts, alignDmgWgh,
				0,
                movie, movieCC,
                alignmentSet.positions[g],
                alignmentSet.initialTracks[g],
                alignmentSet.globComp[g]);
        }
        catch (RelionError XE)
        {
            std::cerr << "warning: unable to load micrograph #" << (g+1) << std::endl;
            continue;
        }

		const int maxRangeP = 2 * motionEstimator->getCCPad() * maxRange;

        #pragma omp parallel for num_threads(nr_omp_threads)
        for (int p = 0; p < pc; p++)
        {
			for (int f = 0; f < fc; f++)
            {
				if (maxRange > 0)
                {
                    movieCC[p][f] = FilterHelper::cropCorner2D(movieCC[p][f], maxRangeP, maxRangeP);
                }

                alignmentSet.copyCC(g, p, f, movieCC[p][f]);

                Image<Complex> pred = reference->predict(
                    mdts[g], p, *obsModel, ReferenceMap::Opposite);

                alignmentSet.accelerate(movie[p][f], alignmentSet.obs[g][p][f]);
                alignmentSet.accelerate(pred, alignmentSet.pred[g][p]);
            }
        }
    }

    // release all unneeded heap space back to the OS
    // (this can free tens of Gb)
#if !defined(__APPLE__)
    malloc_trim(0);
#endif

    std::cout << "   done\n";
}


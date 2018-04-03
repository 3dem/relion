#include "motion_param_estimator.h"
#include "motion_refiner.h"

#include <src/jaz/index_sort.h>
#include <src/jaz/filter_helper.h>

using namespace gravis;

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

    minParticles = textToInteger(parser.getOption("--min_p", "Minimum number of particles on which to estimate the parameters", "1000"));
    rV = textToFloat(parser.getOption("--r_vel", "Test s_vel +/- r_vel * s_vel", "0.5"));
    rD = textToFloat(parser.getOption("--r_div", "Test s_div +/- r_div * s_div", "0.5"));
    rA = textToFloat(parser.getOption("--r_acc", "Test s_acc +/- r_acc * s_acc", "0.5"));
    recursions = textToInteger(parser.getOption("--par_recs", "Parameter estimation is iterated recursively this many times, each time halving the search range", "3"));
    steps = textToInteger(parser.getOption("--par_steps", "Parameter estimation takes max. this many steps before halving the range", "10"));
    maxRange = textToInteger(parser.getOption("--mot_range", "Limit allowed motion range [Px]", "50"));

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
        REPORT_ERROR("ERROR: Cutoff frequency can only be provided in pixels (--k_cut) or Angstrom (--k_cut_A), not both.");
    }

    if (k_cutoff_Angst > 0.0 && k_cutoff < 0.0)
    {
        k_cutoff = obsModel->angToPix(k_cutoff_Angst, s);
    }

    else if (k_cutoff > 0.0 && k_cutoff_Angst < 0.0)
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

    const long mc = allMdts.size();

    std::vector<double> randNums(mc);

    for (int m = 0; m < mc; m++)
    {
        randNums[m] = rand() / (double)RAND_MAX;
    }

    std::vector<int> order = IndexSort<RFLOAT>::sortIndices(randNums);

    int pc = 0;
    mdts.clear();

    for (int i = 0; i < order.size(); i++)
    {
        const int m = order[i];

        const int pcm = allMdts[m].numberOfObjects();

        // motion estimation does not work on one single particle
        if (pcm < 2) continue;

        mdts.push_back(allMdts[m]);
        pc += pcm;

        if (pc > minParticles)
        {
            if (verb > 0)
            {
                std::cout << " + " << pc << " particles found in "
                          << mdts.size() << " micrographs\n";
            }

            break;
        }
    }

    if (verb > 0 && pc < minParticles)
    {
        std::cout << " - Warning: this dataset does not contain " << minParticles
                  << " particles (--min_p) in micrographs with at least 2 particles\n";
    }

    k_out = reference->k_out;

    ready = true;
}

void MotionParamEstimator::run()
{
    if (!ready)
    {
        REPORT_ERROR("ERROR: MotionParamEstimator::run: MotionParamEstimator not initialized.");
    }

    if (!estim2 && !estim3) return;

    prepAlignment();

    if (estim2)
    {
        estimateTwoParamsRec();
    }
}

bool MotionParamEstimator::anythingToDo()
{
    return estim2 || estim3;
}

d4Vector MotionParamEstimator::estimateTwoParamsRec()
{

}

void MotionParamEstimator::prepAlignment()
{
    std::cout << "preparing alignment data...\n";

    const std::vector<Image<RFLOAT>>& dmgWgh = motionEstimator->getDamageWeights();
    std::vector<Image<RFLOAT>> alignDmgWgh(fc);

    for (int f = 0; f < fc; f++)
    {
        alignDmgWgh[f] = FilterHelper::ButterworthEnvFreq2D(dmgWgh[f], k_out-1, k_out+1);
    }

    std::vector<ParFourierTransformer> fts(nr_omp_threads);

    alignmentSet = AlignmentSet(mdts, fc, s, k_cutoff+2, k_out);

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

        std::cout << "    micrograph " << (g+1) << " / " << gc << ": "
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
    }

    std::cout << "done\n";
}


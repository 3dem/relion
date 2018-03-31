#include "motion_param_estimator.h"
#include "motion_refiner.h"

MotionParamEstimator::MotionParamEstimator(MotionRefiner& motionRefiner)
:   motionRefiner(motionRefiner)
{

}

int MotionParamEstimator::read(IOParser &parser, int argc, char *argv[])
{
    parser.addSection("Parameter estimation");

    estim2 = parser.checkOption("--params2", "Estimate 2 parameters instead of motion");
    estim3 = parser.checkOption("--params3", "Estimate 3 parameters instead of motion");
    rV = textToFloat(parser.getOption("--r_vel", "Test s_vel +/- r_vel * s_vel", "0.5"));
    rD = textToFloat(parser.getOption("--r_div", "Test s_div +/- r_div * s_div", "0.5"));
    rA = textToFloat(parser.getOption("--r_acc", "Test s_acc +/- r_acc * s_acc", "0.5"));
    recursions = textToInteger(parser.getOption("--par_recs", "Parameter estimation is iterated recursively this many times, each time halving the search range", "3"));
    steps = textToInteger(parser.getOption("--par_steps", "Parameter estimation takes max. this many steps before halving the range", "10"));
    maxRange = textToInteger(parser.getOption("--mot_range", "Limit allowed motion range [Px]", "50"));

}

void MotionParamEstimator::prepare()
{
    if ((estim2 || estim3) && motionRefiner.k_cutoff < 0)
    {
        REPORT_ERROR("ERROR: Parameter estimation requires a freq. cutoff (--k_cut).");
    }

    if (estim2 && estim3)
    {
        REPORT_ERROR("ERROR: Only 2 or 3 parameters can be estimated (--params2 or --params3), not both.");
    }


}

void MotionParamEstimator::run()
{

}

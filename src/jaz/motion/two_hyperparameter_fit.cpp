#include "two_hyperparameter_fit.h"
#include "motion_param_estimator.h"

using namespace gravis;

double TwoHyperParameterProblem::velScale = 1000.0;
double TwoHyperParameterProblem::divScale = 1.0;
double TwoHyperParameterProblem::accScale = 1000.0;


TwoHyperParameterProblem::TwoHyperParameterProblem(
        MotionParamEstimator& motionParamEstimator, double s_acc)
:   motionParamEstimator(motionParamEstimator),
    s_acc(s_acc)
{
}

double TwoHyperParameterProblem::f(const std::vector<double>& x, void *tempStorage) const
{
    d2Vector vd = problemToMotion(x);

    std::vector<d3Vector> vda(1);
    vda[0] = d3Vector(vd[0], vd[1], s_acc);

    std::vector<double> tsc(1);

    motionParamEstimator.evaluateParams(vda, tsc);

    return -tsc[0];
}

void TwoHyperParameterProblem::report(int iteration, double cost, const std::vector<double>& x) const
{
    d2Vector vd = problemToMotion(x);

    std::cout << iteration << ": " << vd[0] << ", " << vd[1]
              << ", " << s_acc << " @ " << -cost << "\n";
}

d2Vector TwoHyperParameterProblem::problemToMotion(const std::vector<double>& x)
{
    return d2Vector(x[0]/velScale, x[1]/divScale);
}

std::vector<double> TwoHyperParameterProblem::motionToProblem(d2Vector vd)
{
    return std::vector<double>{vd[0]*velScale, vd[1]*divScale};
}

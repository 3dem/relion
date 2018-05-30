#include "two_hyperparameter_fit.h"
#include "motion_param_estimator.h"

using namespace gravis;

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
	
	std::cout.precision(5);
	
    std::cout << iteration << ": \t "
              << vd[0] << "  \t " << vd[1] << "  \t " << s_acc
              << "  \t ";
	
	std::cout.precision(12);
	
	std::cout << -cost << std::endl;
	
	std::cout.precision(5);
}

d2Vector TwoHyperParameterProblem::problemToMotion(const std::vector<double>& x)
{
    return d2Vector(
        x[0] / MotionParamEstimator::velScale,
        x[1] / MotionParamEstimator::divScale);
}

std::vector<double> TwoHyperParameterProblem::motionToProblem(d2Vector vd)
{
    return std::vector<double>{
        vd[0] * MotionParamEstimator::velScale,
        vd[1] * MotionParamEstimator::divScale};
}

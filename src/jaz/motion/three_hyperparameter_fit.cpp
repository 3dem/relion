#include "three_hyperparameter_fit.h"
#include "motion_param_estimator.h"

using namespace gravis;


double ThreeHyperParameterProblem::accThresh = 500.0;
double ThreeHyperParameterProblem::accEps = 1e-5;


ThreeHyperParameterProblem::ThreeHyperParameterProblem(
        MotionParamEstimator& motionParamEstimator)
:   motionParamEstimator(motionParamEstimator)
{
}

double ThreeHyperParameterProblem::f(const std::vector<double>& x, void *tempStorage) const
{
    d3Vector vda = problemToMotion(x);
    std::vector<double> tsc(1);

    std::vector<d3Vector> vdav{vda};

    motionParamEstimator.evaluateParams(vdav, tsc);
	
	/*std::cout.precision(16);
	std::cout << "  " << vda << " => " << tsc[0] << std::endl;
	std::cout.precision(5);*/

    return -tsc[0];
}

void ThreeHyperParameterProblem::report(int iteration, double cost, const std::vector<double>& x) const
{
	d3Vector vda = problemToMotion(x);
	
	std::cout.precision(5);
	
    std::cout << iteration << ": \t "
              << vda[0] << "  \t " << vda[1] << "  \t " << vda[2]
              << "  \t ";
	
	std::cout.precision(12);
	
	std::cout << -cost << std::endl;
	
	std::cout.precision(5);	
}

d3Vector ThreeHyperParameterProblem::problemToMotion(const std::vector<double>& x)
{
    const double w = x[2] / MotionParamEstimator::accScale;

    return d3Vector(
        x[0] / MotionParamEstimator::velScale,
        x[1] / MotionParamEstimator::divScale,
        //(w > 1.0/accThresh ? (w < 1.0/accEps? 1.0/w : 1.0/accEps) : -1.0));
		(x[2] > accEps)? x[2] / MotionParamEstimator::accScale : accEps / MotionParamEstimator::accScale);
}

std::vector<double> ThreeHyperParameterProblem::motionToProblem(d3Vector vd)
{
    const double s = vd[2];
    const double w = s < accThresh? (s > accEps? 1.0/s : 1.0/accEps): 0.0;

    return std::vector<double>{
        vd[0] * MotionParamEstimator::velScale,
        vd[1] * MotionParamEstimator::divScale,
        //w * MotionParamEstimator::accScale};
		vd[2] * MotionParamEstimator::accScale};
}

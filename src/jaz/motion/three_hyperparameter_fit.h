#ifndef THREE_HYPERPARAMETER_PROBLEM_H
#define THREE_HYPERPARAMETER_PROBLEM_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/gravis/t3Vector.h>

class MotionParamEstimator;

class ThreeHyperParameterProblem : public Optimization
{
    public:

        ThreeHyperParameterProblem(
            MotionParamEstimator& motionParamEstimator);

        double f(const std::vector<double>& x, void* tempStorage) const;
        void report(int iteration, double cost, const std::vector<double>& x) const;

        static gravis::d3Vector problemToMotion(const std::vector<double>& x);
        static std::vector<double> motionToProblem(gravis::d3Vector vd);

    protected:

        MotionParamEstimator& motionParamEstimator;

        static double accThresh, accEps;
};

#endif

#ifndef TWO_HYPERPARAMETER_PROBLEM_H
#define TWO_HYPERPARAMETER_PROBLEM_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/gravis/t2Vector.h>

class MotionParamEstimator;

class TwoHyperParameterProblem : public Optimization
{
    public:

        TwoHyperParameterProblem(
            MotionParamEstimator& motionParamEstimator,
            double s_acc);

        double f(const std::vector<double>& x, void* tempStorage) const;
        void report(int iteration, double cost, const std::vector<double>& x) const;

        static gravis::d2Vector problemToMotion(const std::vector<double>& x);
        static std::vector<double> motionToProblem(gravis::d2Vector vd);

    protected:

        MotionParamEstimator& motionParamEstimator;
        double s_acc;

        static double velScale, divScale, accScale;
};

#endif

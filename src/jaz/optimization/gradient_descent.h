#ifndef GRADIENT_DESCENT_H
#define GRADIENT_DESCENT_H

#include <vector>
#include "optimization.h"

class GradientDescent
{
    public:

        static std::vector<double> optimize(const std::vector<double>& initial,
                const DifferentiableOptimization& opt,
                double step, double minStep, double minDiff, long maxIters,
                double inertia = 0.0, bool verbose = false);
};

#endif

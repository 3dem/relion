#ifndef NELDER_MEAD_H
#define NELDER_MEAD_H

#include <vector>
#include <src/jaz/optimization.h>

class NelderMead
{
    public:

        static std::vector<double> optimize(const std::vector<double>& initial,
                const Optimization& opt,
                double initialStep, double tolerance, long maxIters,
                double alpha = 1.0, double gamma = 2.0,
                double rho = 0.5, double sigma = 0.5, bool verbose = false);

        static void test();
};

#endif

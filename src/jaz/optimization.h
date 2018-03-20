#ifndef JAZ_OPTIMIZATION_H
#define JAZ_OPTIMIZATION_H

#include <vector>

// abstract class for optimization problems
class Optimization
{
    public:

    virtual double f(const std::vector<double>& x) const = 0;
};

class DifferentiableOptimization : public Optimization
{
    public:

    virtual double f(const std::vector<double>& x) const = 0;
    virtual void grad(const std::vector<double>& x, std::vector<double>& gradDest) const = 0;
};

class RosenbrockBanana : public Optimization
{
    public:

    double f(const std::vector<double>& x) const;
};

#endif

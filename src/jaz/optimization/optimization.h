#ifndef JAZ_OPTIMIZATION_H
#define JAZ_OPTIMIZATION_H

#include <vector>

// abstract class for optimization problems
class Optimization
{
    public:

    virtual double f(const std::vector<double>& x, void* tempStorage) const = 0;
    virtual void* allocateTempStorage() const {return 0;}
    virtual void deallocateTempStorage(void* ts) const {}
    virtual void report(int iteration, double cost, const std::vector<double>& x) const {}
};

class DifferentiableOptimization : public Optimization
{
    public:

    virtual void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const = 0;
};

class RosenbrockBanana : public DifferentiableOptimization
{
    public:

    double f(const std::vector<double>& x, void* tempStorage) const;
    void grad(const std::vector<double>& x, std::vector<double>& gradDest, void* tempStorage) const;
};

#endif

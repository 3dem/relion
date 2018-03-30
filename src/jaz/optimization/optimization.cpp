#include "optimization.h"

double RosenbrockBanana::f(const std::vector<double> &x, void *tempStorage) const
{
    const double a = 1.0;
    const double b = 100.0;
    const double& xx = x[0];
    const double& yy = x[1];

    return (a - xx) * (a - xx) + b * (yy - xx * xx) * (yy - xx * xx);
}

void RosenbrockBanana::grad(
    const std::vector<double> &x, std::vector<double> &gradDest, void *tempStorage) const
{
    const double a = 1.0;
    const double b = 100.0;
    const double& xx = x[0];
    const double& yy = x[1];

    gradDest[0] = -2.0 * (a - xx) - 4.0 * b * (yy - xx * xx) * xx;
    gradDest[1] = 2.0 * b * (yy - xx * xx);
}

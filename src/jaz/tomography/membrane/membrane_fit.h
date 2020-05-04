
#ifndef MEMBRANE_FIT
#define MEMBRANE_FIT

#include <string>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/cutting.h>
#include <src/jaz/image/resampling.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/error.h>

class MembraneFit : public Optimization
{
    public:

        MembraneFit(const BufferedImage<float>& box, double sigma);

            const BufferedImage<float>& box;
            double sigma;

        double f(const std::vector<double>& x, void* tempStorage) const;
        void report(int iteration, double cost, const std::vector<double>& x) const;

        std::vector<double> average(const std::vector<double> &x) const;
        BufferedImage<float> expand(const std::vector<double> &x) const;
};

#endif

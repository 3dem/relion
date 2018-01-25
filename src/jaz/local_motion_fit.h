#ifndef LOCAL_MOTION_FIT
#define LOCAL_MOTION_FIT

#include <src/image.h>
#include <src/jaz/optimization.h>
#include <vector>

class LocalMotionFit : public DifferentiableOptimization
{
    public:

        LocalMotionFit(
                const std::vector<std::vector<Image<RFLOAT>>>& correlation,
                const std::vector<double>& velWgh,
                const std::vector<double>& accWgh,
                const std::vector<std::vector<std::vector<double>>>& divWgh,
                int threads);

        double f(const std::vector<double>& x) const;
        void grad(const std::vector<double>& x, std::vector<double>& gradDest) const;

    private:

        int pc, fc, threads;
        const std::vector<std::vector<Image<RFLOAT>>>& correlation;
        const std::vector<double>& velWgh;
        const std::vector<double>& accWgh;
        const std::vector<std::vector<std::vector<double>>>& divWgh;
};

#endif

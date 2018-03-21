#ifndef GP_MOTION_FIT
#define GP_MOTION_FIT

#include <src/image.h>
#include <src/jaz/optimization.h>
#include <src/jaz/gravis/t2Vector.h>
#include <vector>

class GpMotionFit : public DifferentiableOptimization
{
    public:

        GpMotionFit(
                const std::vector<std::vector<Image<RFLOAT>>>& correlation,
                double sig_vel_px, double sig_div_px, double sig_acc_px,
                int maxDims,
                const std::vector<gravis::d2Vector>& positions,
                const std::vector<gravis::d2Vector>& perFrameOffsets,
                int threads);

        double f(const std::vector<double>& x) const;
        void grad(const std::vector<double>& x, std::vector<double>& gradDest) const;

        void paramsToPos(const std::vector<double>& x,
                         std::vector<std::vector<gravis::d2Vector>>& pos) const;

        void posToParams(const std::vector<std::vector<gravis::d2Vector>>& pos,
                         std::vector<double>& x) const;

    private:

        int pc, fc, dc, threads;
        double sig_vel_px, sig_div_px, sig_acc_px;

        Matrix2D<RFLOAT> basis;
        std::vector<double> eigenVals;

        const std::vector<std::vector<Image<RFLOAT>>>& correlation;
        const std::vector<gravis::d2Vector>& positions;
        const std::vector<gravis::d2Vector>& perFrameOffsets;
};

#endif

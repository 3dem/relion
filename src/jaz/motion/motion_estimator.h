#ifndef MOTION_ESTIMATOR_H
#define MOTION_ESTIMATOR_H

#include <src/image.h>
#include <src/jaz/gravis/t2Vector.h>
#include <vector>

class MotionRefiner;
class IOParser;
class ParFourierTransformer;

class MotionEstimator
{
    public:

        MotionEstimator(MotionRefiner& motionRefiner);

            MotionRefiner& motionRefiner;

            bool ready;

            int s, sh, fc;
            int maxEDs, maxIters;

            bool unregGlob, noGlobOff, cutoffOut,
                diag, expKer, global_init, debugOpt;

            double dmga, dmgb, dmgc, dosePerFrame,
                sig_vel, sig_div, sig_acc,
                k_cutoff, k_cutoff_Angst, optEps;

            std::vector<Image<RFLOAT>> dmgWeight, dmgWeightEval;


        void read(IOParser& parser, int argc, char *argv[]);
        void init();
        void process(const std::vector<MetaDataTable> &mdts, long g_start, long g_end);

        // load micrograph g and compute all data required for the optimization;
        // also used by MotionParamEstimator
        void prepMicrograph(
            // in:
            const MetaDataTable& mdt, std::vector<ParFourierTransformer>& fts,
            const std::vector<Image<RFLOAT>>& dmgWeight,
            // out:
            std::vector<std::vector<Image<Complex>>>& movie,
            std::vector<std::vector<Image<RFLOAT>>>& movieCC,
            std::vector<gravis::d2Vector>& positions,
            std::vector<std::vector<gravis::d2Vector>>& initialTracks,
            std::vector<gravis::d2Vector>& globComp);

        // perform the actual optimization (also used by MotionParamEstimator)
        std::vector<std::vector<gravis::d2Vector>> optimize(
            const std::vector<std::vector<Image<RFLOAT>>>& movieCC,
            const std::vector<std::vector<gravis::d2Vector>>& inTracks,
            double sig_vel_px, double sig_acc_px, double sig_div_px,
            const std::vector<gravis::d2Vector>& positions,
            const std::vector<gravis::d2Vector>& globComp);

        void updateFCC(
            const std::vector<std::vector<Image<Complex>>>& movie,
            const std::vector<std::vector<gravis::d2Vector>>& tracks,
            const MetaDataTable& mdt,
            std::vector<Image<RFLOAT>>& tables,
            std::vector<Image<RFLOAT>>& weights0,
            std::vector<Image<RFLOAT>>& weights1);

        void writeOutput(
            const std::vector<std::vector<gravis::d2Vector>>& tracks,
            const std::vector<Image<RFLOAT>>& fccData,
            const std::vector<Image<RFLOAT>>& fccWeight0,
            const std::vector<Image<RFLOAT>>& fccWeight1,
            const std::vector<gravis::d2Vector>& positions,
            std::string fn_root, double visScale);

        static bool isFinished(std::string filenameRoot);

};

#endif

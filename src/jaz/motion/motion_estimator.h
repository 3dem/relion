#ifndef MOTION_ESTIMATOR_H
#define MOTION_ESTIMATOR_H

#include <src/image.h>
#include <src/jaz/gravis/t2Vector.h>
#include <vector>

class IOParser;
class ParFourierTransformer;
class ReferenceMap;
class ObservationModel;
class MicrographHandler;

class MotionEstimator
{
    public:

        MotionEstimator();


        void read(IOParser& parser, int argc, char *argv[]);

        void init(int verb, int s, int fc, int nr_omp_threads,
                  bool debug, std::string outPath,
                  ReferenceMap* reference,
                  ObservationModel* obsModel,
                  MicrographHandler* micrographHandler);

        void process(const std::vector<MetaDataTable> &mdts, long g_start, long g_end);


        // load micrograph from mdt and compute all data required for the optimization;
        // positions, initialTracks and globComp need to have the right sizes already (pc, pc*fc, fc)
        // (also used by MotionParamEstimator)
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
            const std::vector<gravis::d2Vector>& globComp) const;

        const std::vector<Image<RFLOAT>>& getDamageWeights();

        bool isReady();

        double getDosePerFrame();
        void proposeDosePerFrame(double dpf, std::string metaFn, int verb);


        static std::vector<MetaDataTable> findUnfinishedJobs(
                const std::vector<MetaDataTable>& mdts, std::string path);

        // translates the given parameters (in A or A/dose) into pixels
        // done in one place to ensure consistency
        double normalizeSigVel(double sig_vel);
        double normalizeSigDiv(double sig_div);
        double normalizeSigAcc(double sig_acc);


    protected:

            bool paramsRead, ready;

            // read from cmd line
            int maxEDs, maxIters;

            bool unregGlob, noGlobOff, cutoffOut,
                diag, expKer, global_init, debugOpt;

            double dmga, dmgb, dmgc, dosePerFrame,
                sig_vel, sig_div, sig_acc, optEps;

            // set at init
            int s, sh, fc, verb, nr_omp_threads;
            double angpix;
            bool debug;

            std::vector<Image<RFLOAT>> dmgWeight;

            std::string outPath;

            ReferenceMap* reference;
            ObservationModel* obsModel;
            MicrographHandler* micrographHandler;


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


        static bool isJobFinished(std::string filenameRoot);

};

#endif

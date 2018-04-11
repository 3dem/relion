#ifndef MOTION_PARAM_ESTIMATOR_H
#define MOTION_PARAM_ESTIMATOR_H

#include "alignment_set.h"

#include <src/image.h>
#include <src/time.h>

#include <src/jaz/gravis/t4Vector.h>


//#define TIMING 1

#ifdef TIMING
    #define RCTIC(timer,label) (timer.tic(label))
    #define RCTOC(timer,label) (timer.toc(label))
#else
    #define RCTIC(timer,label)
    #define RCTOC(timer,label)
#endif

class MotionEstimator;
class ReferenceMap;
class ObservationModel;
class ParFourierTransformer;

class MotionParamEstimator
{
    public:

        MotionParamEstimator();


        int read(IOParser& parser, int argc, char *argv[]);

        void init(int verb, int nr_omp_threads, bool debug,
                  int s, int fc,
                  const std::vector<MetaDataTable>& allMdts,
                  MotionEstimator* motionEstimator,
                  ReferenceMap* reference,
                  ObservationModel* obsModel);

        void run();

        bool anythingToDo();

        // to be used by instances of OptimizationProblem
        void evaluateParams(const std::vector<gravis::d3Vector>& sig_vals,
                            std::vector<double>& TSCs);

        static const double velScale, divScale, accScale;


    protected:

            bool paramsRead, ready;

            AlignmentSet<float> alignmentSet;

            // read from cmd. line:
            bool estim2, estim3;
            int minParticles, maxRange, maxIters, seed;
            double sV, sD, sA;
            double iniStep, conv;
            double k_cutoff, k_cutoff_Angst;
            double k_eval, k_eval_Angst;

            // set at init:
            std::vector<MetaDataTable> mdts;

            MotionEstimator* motionEstimator;
            ObservationModel* obsModel;
            ReferenceMap* reference;

            int fc, s, k_out, verb, nr_omp_threads;
            bool debug;

            #ifdef TIMING
                Timer paramTimer;
                int timeSetup, timeOpt, timeEval;
            #endif


        gravis::d4Vector estimateTwoParamsNM(
                double sig_v_0, double sig_d_0, double sig_acc,
                double inStep, double conv, int maxIters);

        gravis::d4Vector estimateThreeParamsNM(
                double sig_v_0, double sig_d_0, double sig_a_0,
                double inStep, double conv, int maxIters);


        void prepAlignment();
};

#endif

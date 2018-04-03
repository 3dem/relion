#ifndef MOTION_PARAM_ESTIMATOR_H
#define MOTION_PARAM_ESTIMATOR_H

#include "alignment_set.h"

#include <src/image.h>
#include <src/jaz/gravis/t4Vector.h>

class MotionEstimator;
class ReferenceMap;
class ObservationModel;

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


    protected:

            bool paramsRead, ready;

            AlignmentSet alignmentSet;

            // read from cmd. line:
            bool estim2, estim3;
            int minParticles, maxRange, recursions, steps;
            double rV, rD, rA;
            double k_cutoff, k_cutoff_Angst;

            // set at init:
            std::vector<MetaDataTable> mdts;

            MotionEstimator* motionEstimator;
            ObservationModel* obsModel;
            ReferenceMap* reference;

            int fc, s, k_out, verb, nr_omp_threads;
            bool debug;


        gravis::d4Vector estimateTwoParamsRec();

        void prepAlignment();
};

#endif

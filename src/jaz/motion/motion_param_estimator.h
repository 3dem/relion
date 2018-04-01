#ifndef MOTION_PARAM_ESTIMATOR_H
#define MOTION_PARAM_ESTIMATOR_H

#include "alignment_set.h"
#include <src/jaz/gravis/t4Vector.h>

class MotionRefiner;

class MotionParamEstimator
{
    public:

        MotionParamEstimator(MotionRefiner& motionRefiner);

            MotionRefiner& motionRefiner;
            AlignmentSet alignmentSet;

            bool ready, estim2, estim3;
            int maxRange, recursions, steps;
            double rV, rD, rA;


        int read(IOParser& parser, int argc, char *argv[]);

        void init();
        void run();

        gravis::d4Vector estimateTwoParamsRec();

        void prepAlignment(int k_out,
                const std::vector<Image<RFLOAT>>& dmgWeight0,
                const std::vector<Image<RFLOAT>>& dmgWeight);
};

#endif

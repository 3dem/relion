#ifndef MOTION_PARAM_ESTIMATOR_H
#define MOTION_PARAM_ESTIMATOR_H

#include <src/jaz/alignment_set.h>

class MotionRefiner;

class MotionParamEstimator
{
    public:

        MotionParamEstimator(MotionRefiner& motionRefiner);

            MotionRefiner& motionRefiner;
            AlignmentSet alignmentSet;

            bool estim2, estim3;
            int maxRange, recursions, steps;
            double rV, rD, rA;


        int read(IOParser& parser, int argc, char *argv[]);

        void prepare();
        void run();

        gravis::d2Vector estimateTwoParamsRec();
};

#endif

#ifndef FRAME_RECOMBINER_H
#define FRAME_RECOMBINER_H

#include <src/image.h>
#include <vector>
#include <string>

class MotionRefiner;
class IOParser;

class FrameRecombiner
{
    public:

        FrameRecombiner(MotionRefiner& motionRefiner);

            MotionRefiner& motionRefiner;

            int s, sh, fc;

            bool doCombineFrames, bfac_diag;
            int k0, k1;
            double k0a, k1a;
            std::string trackFn, bfacFn;

        void read(IOParser& parser, int argc, char *argv[]);
        void init();
        void process(long g_start, long g_end);


    protected:

        std::vector<Image<RFLOAT>> weightsFromFCC();
        std::vector<Image<RFLOAT>> weightsFromBfacs();
};

#endif

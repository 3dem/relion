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
            std::vector<Image<RFLOAT>> freqWeights;

        void read(IOParser& parser, int argc, char *argv[]);
        void init(const std::vector<MetaDataTable>& allMdts);
        void process(const std::vector<MetaDataTable>& mdts, long g_start, long g_end);

        static bool isFinished(std::string filenameRoot);


    protected:

        std::vector<Image<RFLOAT>> weightsFromFCC(const std::vector<MetaDataTable>& allMdts);
        std::vector<Image<RFLOAT>> weightsFromBfacs();
};

#endif

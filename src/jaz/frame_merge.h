#ifndef FRAME_MERGE_H
#define FRAME_MERGE_H

#include <src/image.h>

class FrameMerge
{
    public:

        static void mergeAvg(Image<RFLOAT>& stack, Image<RFLOAT>& tgt);
        static void valueHistogram(Image<RFLOAT>& stack, Image<RFLOAT>& tgt);
};

#endif

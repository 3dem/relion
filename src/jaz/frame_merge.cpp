#include <src/jaz/frame_merge.h>

void FrameMerge :: mergeAvg(Image<RFLOAT>& stack, Image<RFLOAT>& tgt)
{
    const int bs = stack.data.xdim / tgt.data.xdim;
    const int fc = stack.data.zdim;

    for (int y = 0; y < tgt.data.ydim; y++)
    for (int x = 0; x < tgt.data.xdim; x++)
    {
        double sum = 0.0;

        for (int yb = 0; yb < bs; yb++)
        for (int xb = 0; xb < bs; xb++)
        for (int n = 0; n < fc; n++)
        {
            sum += DIRECT_NZYX_ELEM(stack.data, 0, n, y*bs + yb, x*bs + xb);
        }

        DIRECT_A2D_ELEM(tgt.data, y, x) = sum / (double)(bs*bs*fc);
    }
}

void FrameMerge :: valueHistogram(Image<RFLOAT>& stack, Image<RFLOAT>& tgt)
{
    const int mv = tgt.data.xdim;
    std::vector<int> bins(mv);

    for (int i = 0; i < mv; i++)
    {
        bins[i] = 0;
    }

    for (int z = 0; z < tgt.data.zdim; z++)
    for (int y = 0; y < tgt.data.ydim; y++)
    for (int x = 0; x < tgt.data.xdim; x++)
    {
        double v = DIRECT_NZYX_ELEM(stack.data, 0, z, y, x);

        int vb = (int)v;

        if (vb < 0) vb = 0;
        else if (vb > mv) vb = mv;

        if (!(vb == vb)) continue;

        bins[vb]++;
    }

    double bmax = 0;
    for (int i = 1; i < mv; i++)
    {
        if (bins[i] > bmax) bmax = bins[i];
    }

    bmax += 2.0;

    for (int x = 0; x < tgt.data.xdim; x++)
    {
        double bv = (double)tgt.data.ydim * (double)bins[x] / bmax;

        for (int y = 0; y < tgt.data.ydim; y++)
        {
            DIRECT_A2D_ELEM(tgt.data, y, x) = y >= bv? 1.0 : 0.0;
        }
    }
}

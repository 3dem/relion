#include <src/jaz/fftw_helper.h>

void FftwHelper::decenterUnflip2D(const MultidimArray<double> &src, MultidimArray<double> &dest)
{
    const long int w = src.xdim;
    const long int h = src.ydim;

    dest.reshape(h, 2*(w - 1));

    const long int yc = dest.ydim/2;

    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        long int xs = x - w;

        if (xs < 0)
        {
            long int ys = (y + yc - 1) % dest.ydim;
            DIRECT_A2D_ELEM(dest, y, x) = -DIRECT_A2D_ELEM(src, h-ys-1, -xs-1);
        }
        else
        {
            long int ys = (y + yc) % dest.ydim;
            DIRECT_A2D_ELEM(dest, y, x) =  DIRECT_A2D_ELEM(src, ys, xs);
        }
    }
}

void FftwHelper::decenterDouble2D(const MultidimArray<double> &src, MultidimArray<double> &dest)
{
    const long int w = src.xdim;
    const long int h = src.ydim;

    dest.reshape(h, 2*(w - 1));

    const long int yc = dest.ydim/2;

    for (long int y = 0; y < dest.ydim; y++)
    for (long int x = 0; x < dest.xdim; x++)
    {
        long int xs = x - w;

        if (xs < 0)
        {
            long int ys = (y + yc - 1) % dest.ydim;
            DIRECT_A2D_ELEM(dest, y, x) = DIRECT_A2D_ELEM(src, h-ys-1, -xs-1);
        }
        else
        {
            long int ys = (y + yc) % dest.ydim;
            DIRECT_A2D_ELEM(dest, y, x) =  DIRECT_A2D_ELEM(src, ys, xs);
        }
    }
}

#include "membrane_fit.h"

using namespace gravis;

MembraneFit::MembraneFit(const BufferedImage<float> &box, double sigma)
:   box(box),
    sigma(sigma)
{
}

double MembraneFit::f(const std::vector<double> &x, void *tempStorage) const
{
    /*
     * [0 1 2 3]
     * [1 4 5 6]
     * [2 5 7 8]
     * [3 6 8 9]
     * */

    const size_t w = box.xdim;
    const size_t h = box.ydim;
    const size_t d = box.zdim;

    std::vector<double> avg = average(x);

    const size_t s = avg.size();
    const int orig = s/2;

    double out(0.0);

        /*std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
        std::cout << x[1] << " " << x[4] << " " << x[5] << " " << x[6] << std::endl;
        std::cout << x[2] << " " << x[5] << " " << x[7] << " " << x[8] << std::endl;
        std::cout << x[3] << " " << x[6] << " " << x[8] << " " << 0 << std::endl << std::endl;*/

    for (size_t zi = 0; zi < d; zi++)
    for (size_t yi = 0; yi < h; yi++)
    for (size_t xi = 0; xi < w; xi++)
    {
        const double xx = xi - w/2.0;
        const double yy = yi - h/2.0;
        const double zz = zi - d/2.0;

        const double f = x[0] * xx * xx  +  2 * x[1] * xx * yy  +  2 * x[2] * xx * zz  +  2 * x[3] * xx
                       + x[4] * yy * yy  +  2 * x[5] * yy * zz  +  2 * x[6] * yy
                       + x[7] * zz * zz  +  2 * x[8] * zz  +  orig;

        const int f0 = (int) f;
        const int f1 = f0 + 1;
        const double df = f - f0;

        const float val = box(xi,yi,zi);

        double su(0.0), wg(0.0);

        if (f0 >= 0 && f0 < s)
        {
            su += (1.0 - df) * avg[f0];
            wg += (1.0 - df);
        }

        if (f1 >= 0 && f1 < s)
        {
            su += df * avg[f1];
            wg += df;
        }

        const double ival = wg > 0.0? su / wg : 0.0;

        const double err = val - ival;

        out += err * err;
    }

    return out;
}

void MembraneFit::report(int iteration, double cost, const std::vector<double> &x) const
{
    std::cout << iteration << ": " << cost << std::endl;
}

std::vector<double> MembraneFit::average(const std::vector<double> &x) const
{
    const size_t w = box.xdim;
    const size_t h = box.ydim;
    const size_t d = box.zdim;

    const double diag = sqrt(w*w + h*h + d*d);
    const size_t s = (int)(2*diag);

    const int orig = s/2;

    std::vector<double> sum(s, 0.0);
    std::vector<double> wgh(s, 0.0);

    for (size_t zi = 0; zi < d; zi++)
    for (size_t yi = 0; yi < h; yi++)
    for (size_t xi = 0; xi < w; xi++)
    {
        const double xx = xi - w/2.0;
        const double yy = yi - h/2.0;
        const double zz = zi - d/2.0;

        const double f = x[0] * xx * xx  +  2 * x[1] * xx * yy  +  2 * x[2] * xx * zz  +  2 * x[3] * xx
                       + x[4] * yy * yy  +  2 * x[5] * yy * zz  +  2 * x[6] * yy
                       + x[7] * zz * zz  +  2 * x[8] * zz  +  orig;

            /*if (xi%10 == 0 && yi%10 == 0 && zi%10 == 0)
            {
                std::cout << "AVG: " << xi << ", " << yi << ", " << zi << ": " << f << std::endl;
            }*/

        const int f0 = (int) f;
        const int f1 = f0 + 1;
        const double df = f - f0;

        const float val = box(xi,yi,zi);

        if (f0 >= 0 && f0 < s)
        {
            sum[f0] += (1.0 - df) * val;
            wgh[f0] += (1.0 - df);
        }

        if (f1 >= 0 && f1 < s)
        {
            sum[f1] += df * val;
            wgh[f1] += df;
        }

            /*if (xi%10 == 0 && yi%10 == 0 && zi%10 == 0)
            {
                std::cout << "AVG: " << f0 << ", " << f1 << ", " << std::endl;
            }*/


    }

    for (int i = 0; i < s; i++)
    {
        if (wgh[i] > 0.0)
        {
            sum[i] /= wgh[i];
        }
    }

    return sum;
}

BufferedImage<float> MembraneFit::expand(const std::vector<double> &x) const
{
    const size_t w = box.xdim;
    const size_t h = box.ydim;
    const size_t d = box.zdim;

    std::vector<double> avg = average(x);

    const size_t s = avg.size();
    const int orig = s/2;

    BufferedImage<float> out(w,h,d);

    for (size_t zi = 0; zi < d; zi++)
    for (size_t yi = 0; yi < h; yi++)
    for (size_t xi = 0; xi < w; xi++)
    {
        const double xx = xi - w/2.0;
        const double yy = yi - h/2.0;
        const double zz = zi - d/2.0;

        const double f = x[0] * xx * xx  +  2 * x[1] * xx * yy  +  2 * x[2] * xx * zz  +  2 * x[3] * xx
                       + x[4] * yy * yy  +  2 * x[5] * yy * zz  +  2 * x[6] * yy
                       + x[7] * zz * zz  +  2 * x[8] * zz  +  orig;

                /*if (xi%10 == 0 && yi%10 == 0 && zi%10 == 0)
                {

                    std::cout << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
                    std::cout << x[1] << " " << x[4] << " " << x[5] << " " << x[6] << std::endl;
                    std::cout << x[2] << " " << x[5] << " " << x[7] << " " << x[8] << std::endl;
                    std::cout << x[3] << " " << x[6] << " " << x[8] << " " << 0 << std::endl << std::endl;

                    std::cout << "EXP: " << xi << ", " << yi << ", " << zi << ": " << f << std::endl;

                    std::cout << "test1: " << (x[0] * xx * xx  +  2 * x[1] * xx * yy  +  2 * x[2] * xx * zz  +  2 * x[3] * xx
                            + x[4] * yy * yy  +  2 * x[5] * yy * zz  +  2 * x[6] * yy
                            + x[7] * zz * zz  +  2 * x[8] * zz  +  orig) << std::endl;

                    std::cout << "test2: " << (x[0] * xx * xx  +  2 * x[1] * xx * yy  +  2 * x[2] * xx * zz  +  2 * x[3] * xx
                            + x[4] * yy * yy  +  2 * x[5] * yy * zz  +  2 * x[6] * yy
                            + x[7] * zz * zz  +  2 * x[8] * zz) << std::endl;

                    std::cout << "xx: " << xx << std::endl;
                    std::cout << "yy: " << yy << std::endl;
                    std::cout << "zz: " << zz << std::endl;
                }*/

        const int f0 = (int) f;
        const int f1 = f0 + 1;
        const double df = f - f0;

        double su(0.0), wg(0.0);

        if (f0 >= 0 && f0 < s)
        {
            su += (1.0 - df) * avg[f0];
            wg += (1.0 - df);
        }

        if (f1 >= 0 && f1 < s)
        {
            su += df * avg[f1];
            wg += df;
        }

                /*if (xi%10 == 0 && yi%10 == 0 && zi%10 == 0)
                {
                    std::cout << "EXP: " << f0 << ", " << f1 << ", " << su << std::endl;
                }*/

        const double ival = wg > 0.0? su / wg : 0.0;

        out(xi,yi,zi) = ival;
    }

    return out;
}

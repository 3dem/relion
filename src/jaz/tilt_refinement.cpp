#include <src/jaz/tilt_refinement.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/nelder_mead.h>

using namespace gravis;

void TiltRefinement::updateTiltShift(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        CTF &ctf, RFLOAT angpix,
        Image<Complex>& xyDest,
        Image<RFLOAT>& wDest)
{

    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    const RFLOAT as = (RFLOAT)h * angpix;

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        const double xf = x;
        const double yf = y < w? y : y - h;

        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);

        RFLOAT c = ctf.getCTF(xf/as, yf/as);

        DIRECT_A2D_ELEM(xyDest.data, y, x) += c * vx.conj() * vy;
        DIRECT_A2D_ELEM(wDest.data, y, x) += c * c * vx.norm();
    }
}

void TiltRefinement::updateTiltShiftPar(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        CTF &ctf, RFLOAT angpix,
        Image<Complex>& xyDest,
        Image<RFLOAT>& wDest)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    const RFLOAT as = (RFLOAT)h * angpix;

    #pragma omp parallel for
    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++)
    {
        const double xf = x;
        const double yf = y < w? y : y - h;

        Complex vx = DIRECT_A2D_ELEM(prediction.data, y, x);
        Complex vy = DIRECT_A2D_ELEM(observation.data, y, x);

        RFLOAT c = ctf.getCTF(xf/as, yf/as);

        DIRECT_A2D_ELEM(xyDest.data, y, x) += c * vx.conj() * vy;
        DIRECT_A2D_ELEM(wDest.data, y, x) += c * c * vx.norm();
    }
}

void TiltRefinement::fitTiltShift(const Image<RFLOAT>& phase,
                             const Image<RFLOAT>& weight,
                             RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
                             RFLOAT* shift_x, RFLOAT* shift_y,
                             RFLOAT* tilt_x, RFLOAT* tilt_y,
                             Image<RFLOAT>* fit,
                             d2Matrix magCorr)
{
    const long w = phase.data.xdim;
    const long h = phase.data.ydim;

    double axx = 0.0, axy = 0.0, axz = 0.0,//axw == ayz,
                      ayy = 0.0, ayz = 0.0, ayw = 0.0,
                                 azz = 0.0, azw = 0.0,
                                            aww = 0.0;

    double bx = 0.0, by = 0.0, bz = 0.0, bw = 0.0;

    const RFLOAT as = (RFLOAT)h * angpix;

    for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++)
    {
        double x = xi;
        double y = yi < w? yi : ((yi-h));

        d2Vector p = magCorr * d2Vector(x,y);
        x = p.x/as;
        y = p.y/as;

        double q = x*x + y*y;

        double v = DIRECT_A2D_ELEM(phase.data, yi, xi);
        double g = DIRECT_A2D_ELEM(weight.data, yi, xi);

        axx += g     * x * x;
        axy += g     * x * y;
        axz += g * q * x * x;

        ayy += g     * y * y;
        ayz += g * q * x * y;
        ayw += g * q * y * y;

        azz += g * q * q * x * x;
        azw += g * q * q * x * y;

        aww += g * q * q * y * y;

        bx += g * x * v;
        by += g * y * v;
        bz += g * q * x * v;
        bw += g * q * y * v;
    }

    gravis::d4Matrix A;
    gravis::d4Vector b(bx, by, bz, bw);

    A(0,0) = axx;
    A(0,1) = axy;
    A(0,2) = axz;
    A(0,3) = ayz;

    A(1,0) = axy;
    A(1,1) = ayy;
    A(1,2) = ayz;
    A(1,3) = ayw;

    A(2,0) = axz;
    A(2,1) = ayz;
    A(2,2) = azz;
    A(2,3) = azw;

    A(3,0) = ayz;
    A(3,1) = ayw;
    A(3,2) = azw;
    A(3,3) = aww;

    gravis::d4Matrix Ainv = A;
    Ainv.invert();

    gravis::d4Vector opt = Ainv * b;

    //std::cout << opt[0] << ", " << opt[1] << ", " << opt[2] << ", " << opt[3] << "\n";

    *shift_x = opt[0];
    *shift_y = opt[1];

    *tilt_x = -opt[2]*180.0/(0.360 * Cs * 10000000 * lambda * lambda * 3.141592654);
    *tilt_y = -opt[3]*180.0/(0.360 * Cs * 10000000 * lambda * lambda * 3.141592654);

    /*destA = gravis::d2Vector(opt[0], opt[1]).length();
    destPhiA = (180.0/3.1416)*std::atan2(opt[1], opt[0]);

    destB = gravis::d2Vector(opt[2], opt[3]).length();
    destPhiB = (180.0/3.1416)*std::atan2(opt[3], opt[2]);

    std::cout << "linear: " << destA << " @ " << destPhiA << "°\n";
    std::cout << "cubic:  " << destB << " @ " << destPhiB << "°\n";
    std::cout << "    =  -" << destB << " @ " << (destPhiB + 180.0) << "°\n";*/

    //std::cout << "tilt_x = " << *tilt_x << "\n";
    //std::cout << "tilt_y = " << *tilt_y << "\n";

    if (fit != 0)
    {
        *fit = Image<RFLOAT>(w,h);
        drawPhaseShift(opt[0], opt[1], opt[2], opt[3], w, h, as, magCorr, fit);
    }
}

void TiltRefinement::optimizeTilt(
    const Image<Complex> &xy,
    const Image<RFLOAT> &weight,
    RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
    bool L1,
    RFLOAT shift0_x, RFLOAT shift0_y,
    RFLOAT tilt0_x, RFLOAT tilt0_y,
    RFLOAT *shift_x, RFLOAT *shift_y,
    RFLOAT *tilt_x, RFLOAT *tilt_y,
    Image<RFLOAT> *fit)
{
    TiltOptimization prob(xy, weight, angpix, L1, false);

    double scale = 180.0/(0.360 * Cs * 10000000 * lambda * lambda * 3.141592654);

    std::vector<double> initial{shift0_x, shift0_y, -tilt0_x/scale, -tilt0_y/scale};

    std::vector<double> opt = NelderMead::optimize(
                initial, prob, 0.01, 0.000001, 100000,
                1.0, 2.0, 0.5, 0.5, false);

    *shift_x = opt[0];
    *shift_y = opt[1];

    *tilt_x = -opt[2]*scale;
    *tilt_y = -opt[3]*scale;

    //std::cout << opt[0] << ", " << opt[1] << ", " << opt[2] << ", " << opt[3] << "\n";
    //std::cout << "tilt_x = " << *tilt_x << "\n";
    //std::cout << "tilt_y = " << *tilt_y << "\n";

    if (fit != 0)
    {
        const int w = xy.data.xdim;
        const int h = xy.data.ydim;
        const RFLOAT as = (RFLOAT)h * angpix;

        *fit = Image<RFLOAT>(w,h);
        drawPhaseShift(opt[0], opt[1], opt[2], opt[3], w, h, as, d2Matrix(), fit);
    }
}

void TiltRefinement::optimizeAnisoTilt(
    const Image<Complex> &xy,
    const Image<RFLOAT> &weight,
    RFLOAT Cs, RFLOAT lambda, RFLOAT angpix,
    bool L1,
    RFLOAT shift0_x, RFLOAT shift0_y,
    RFLOAT tilt0_x, RFLOAT tilt0_y,
    RFLOAT *shift_x, RFLOAT *shift_y,
    RFLOAT *tilt_x, RFLOAT *tilt_y,
    RFLOAT* tilt_xx, RFLOAT* tilt_xy, RFLOAT* tilt_yy,
    Image<RFLOAT> *fit)
{
    TiltOptimization prob(xy, weight, angpix, L1, true);

    double scale = 180.0/(0.360 * Cs * 10000000 * lambda * lambda * 3.141592654);

    std::vector<double> initial{shift0_x, shift0_y, -tilt0_x/scale, -tilt0_y/scale, 0.0, 1.0};

    std::vector<double> opt = NelderMead::optimize(
                initial, prob, 0.01, 0.000001, 100000,
                1.0, 2.0, 0.5, 0.5, false);

    *shift_x = opt[0];
    *shift_y = opt[1];

    *tilt_x = -opt[2]*scale;
    *tilt_y = -opt[3]*scale;

    *tilt_xx = 1.0;
    *tilt_xy = opt[4];
    *tilt_yy = opt[5];

    std::cout << opt[0] << ", " << opt[1] << ", " << opt[2]
        << ", " << opt[3] << ", " << opt[4] << ", " << opt[5] << "\n";

    std::cout << "tilt_x = " << *tilt_x << "\n";
    std::cout << "tilt_y = " << *tilt_y << "\n";
    std::cout << "tilt_xx = " << *tilt_xx << "\n";
    std::cout << "tilt_xy = " << *tilt_xy << "\n";
    std::cout << "tilt_yy = " << *tilt_yy << "\n";

    if (fit != 0)
    {
        const int w = xy.data.xdim;
        const int h = xy.data.ydim;
        const RFLOAT as = (RFLOAT)h * angpix;

        *fit = Image<RFLOAT>(w,h);
        drawPhaseShift(opt[0], opt[1], opt[2], opt[3], 1.0, opt[4], opt[5], w, h, as, d2Matrix(), fit);
    }
}

void TiltRefinement::drawPhaseShift(
        double shift_x, double shift_y,
        double tilt_x, double tilt_y,
        int w, int h, double as,
        gravis::d2Matrix magCorr,
        Image<double> *tgt)
{
    for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++)
    {
        double x = xi;
        double y = yi < w? yi : yi-h;

        d2Vector p = magCorr * d2Vector(x,y);
        x = p.x/as;
        y = p.y/as;

        double q = x*x + y*y;

        DIRECT_A2D_ELEM(tgt->data, yi, xi) = x * shift_x + y * shift_y + q * x * tilt_x + q * y * tilt_y;
    }
}

void TiltRefinement::drawPhaseShift(
        double shift_x, double shift_y,
        double tilt_x, double tilt_y,
        double tilt_xx, double tilt_xy, double tilt_yy,
        int w, int h, double as,
        gravis::d2Matrix magCorr,
        Image<double> *tgt)
{
    for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++)
    {
        double x = xi;
        double y = yi < w? yi : yi-h;

        d2Vector p = magCorr * d2Vector(x,y);
        x = p.x/as;
        y = p.y/as;

        double q = tilt_xx * x * x + 2.0 * tilt_xy * x * y + tilt_yy * y * y;

        DIRECT_A2D_ELEM(tgt->data, yi, xi) = x * shift_x + y * shift_y + q * x * tilt_x + q * y * tilt_y;
    }
}

TiltOptimization::TiltOptimization(const Image<Complex> &xy,
    const Image<double> &weight, double angpix,
    bool L1,
    bool anisotropic)
    :   xy(xy),
        weight(weight),
        angpix(angpix),
        L1(L1),
        anisotropic(anisotropic)
{
}

double TiltOptimization::f(const std::vector<double> &x) const
{
    double out = 0.0;

    const int w = xy.data.xdim;
    const int h = xy.data.ydim;

    const double as = (double)h * angpix;

    for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++)
    {
        const double xd = xi/as;
        const double yd = yi < w? yi/as : (yi-h)/as;
        double rr;

        if (anisotropic)
        {
            rr = xd*xd + 2.0*x[4]*xd*yd + x[5]*yd*yd;
        }
        else
        {
            rr = xd*xd + yd*yd;
        }

        const double phi = x[0]*xd + x[1]*yd + rr*x[2]*xd + rr*x[3]*yd;
        const double e = (xy(yi,xi) - Complex(cos(phi), sin(phi))).norm();

        if (L1)
        {
            out += weight(yi,xi) * sqrt(e);
        }
        else
        {
            out += weight(yi,xi) * e;
        }

        /*double d = phi - atan2(xy(yi,xi).imag, xy(yi,xi).real);
        out += weight(yi,xi) * d*d;*/
    }

    return out;
}


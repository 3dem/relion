#include <src/jaz/convolution_helper.h>

Image<double> ConvolutionHelper::convolve2D(Image<double> &img0, Image<double> &img1)
{
    FourierTransformer ft;
    return convolve2D(img0, img1, ft);
}

Image<RFLOAT> ConvolutionHelper::convolve2D(Image<RFLOAT> &img0, Image<RFLOAT> &img1, FourierTransformer &ft)
{
    int w = img0().xdim;
    int wf = w/2 + 1;
    int h = img0().ydim;

    Image<Complex> I0, I1, P(wf,h);

    ft.FourierTransform(img0(), I0());
    ft.FourierTransform(img1(), I1());

    double sc = w*h;

    for (int y = 0; y < h; y++)
    for (int x = 0; x < wf; x++)
    {
        DIRECT_A2D_ELEM(P.data, y, x) = sc * DIRECT_A2D_ELEM(I0.data, y, x) * DIRECT_A2D_ELEM(I1.data, y, x).conj();
    }

    Image<RFLOAT> out(w,h);
    ft.inverseFourierTransform(P(), out());

    return out;
}

Image<RFLOAT> ConvolutionHelper::gaussianKernel2D(double sigma, int w, int h, bool normalize, bool centered, bool half)
{
    Image<RFLOAT> out(w,h);

    const double s2 = sigma*sigma;

    double sum = 0.0;

    for (int yy = 0; yy < h; yy++)
    for (int xx = 0; xx < w; xx++)
    {
        double x, y;

        if (centered)
        {
            if (half)
            {
                x = xx;
                y = yy - h/2 - 1;
            }
            else
            {
                x = xx - w/2 - 1;
                y = yy - h/2 - 1;
            }
        }
        else
        {
            if (half)
            {
                x = xx;
                y = yy <= h/2 + 1? yy : yy - h;
            }
            else
            {
                x = xx <= w/2 + 1? xx : xx - w;
                y = yy <= h/2 + 1? yy : yy - h;
            }
        }

        out(yy,xx) = exp(-0.5*(x*x+y*y)/s2);
        sum += out(yy,xx);
    }

    if (normalize)
    {
        for (int yy = 0; yy < h; yy++)
        for (int xx = 0; xx < w; xx++)
        {
            out(yy,xx) /= sum;
        }
    }

    return out;
}

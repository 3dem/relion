#include "new_ft.h"

#include "src/macros.h"
#include "src/fftw.h"
#include "src/args.h"
#include <string.h>
#include <math.h>

pthread_mutex_t NewFFT::fftw_plan_mutex_new = PTHREAD_MUTEX_INITIALIZER;


void NewFFT::FourierTransform(
    MultidimArray<double>& src,
    MultidimArray<dComplex>& dest,
    const NewFFT::DoublePlan& plan,
    Normalization normalization)
{
    if (!plan.isCompatible(src))
    {
        REPORT_ERROR("ERROR: plan incompatible with input array\n");
    }

    if (!plan.isCompatible(dest))
    {
        dest.resizeNoCp(src.ndim, src.zdim, src.ydim, src.xdim);
    }

    fftw_execute_dft_r2c(plan.getForward(),
        MULTIDIM_ARRAY(src), (fftw_complex*) MULTIDIM_ARRAY(dest));

    if (normalization == FwdOnly)
    {
        const double scale = MULTIDIM_SIZE(src);

        for (long int i = 0; i < NZYXSIZE(dest); i++)
        {
            dest.data[i] /= scale;
        }
    }
    else if (normalization == Both)
    {
        const double scale = sqrt(MULTIDIM_SIZE(src));

        for (long int i = 0; i < NZYXSIZE(dest); i++)
        {
            dest.data[i] /= scale;
        }
    }
}

void NewFFT::inverseFourierTransform(
    MultidimArray<dComplex>& src,
    MultidimArray<double>& dest,
    const NewFFT::DoublePlan& plan,
    Normalization normalization,
    bool preserveInput)
{
    if (!plan.isCompatible(src))
    {
        REPORT_ERROR("ERROR: plan incompatible with input array\n");
    }

    if (!plan.isCompatible(dest))
    {
        dest.resizeNoCp(src.ndim, src.zdim, src.ydim, src.xdim);
    }

    fftw_complex* in;
    MultidimArray<dComplex> src2;

    if (preserveInput)
    {
        src2 = src;
        in = (fftw_complex*) MULTIDIM_ARRAY(src2);
    }
    else
    {
        in = (fftw_complex*) MULTIDIM_ARRAY(src);
    }

    fftw_execute_dft_c2r(plan.getBackward(),
        in, MULTIDIM_ARRAY(dest));

    if (normalization == Both)
    {
        const double scale = sqrt(MULTIDIM_SIZE(dest));

        for (long int i = 0; i < NZYXSIZE(dest); i++)
        {
            dest.data[i] /= scale;
        }
    }
}

void NewFFT::FourierTransform(
    MultidimArray<float>& src,
    MultidimArray<fComplex>& dest,
    const NewFFT::FloatPlan& plan,
    Normalization normalization)
{
    if (!plan.isCompatible(src))
    {
        REPORT_ERROR("ERROR: plan incompatible with input array\n");
    }

    if (!plan.isCompatible(dest))
    {
        dest.resizeNoCp(src.ndim, src.zdim, src.ydim, src.xdim);
    }

    fftwf_execute_dft_r2c(plan.getForward(),
        MULTIDIM_ARRAY(src), (fftwf_complex*) MULTIDIM_ARRAY(dest));

    if (normalization == FwdOnly)
    {
        const float scale = MULTIDIM_SIZE(src);

        for (long int i = 0; i < NZYXSIZE(dest); i++)
        {
            dest.data[i] /= scale;
        }
    }
    else if (normalization == Both)
    {
        const float scale = sqrt(MULTIDIM_SIZE(src));

        for (long int i = 0; i < NZYXSIZE(dest); i++)
        {
            dest.data[i] /= scale;
        }
    }
}

void NewFFT::inverseFourierTransform(
    MultidimArray<fComplex>& src,
    MultidimArray<float>& dest,
    const NewFFT::FloatPlan& plan,
    Normalization normalization,
    bool preserveInput)
{
    if (!plan.isCompatible(src))
    {
        REPORT_ERROR("ERROR: plan incompatible with input array\n");
    }

    if (!plan.isCompatible(dest))
    {
        dest.resizeNoCp(src.ndim, src.zdim, src.ydim, src.xdim);
    }

    fftwf_complex* in;
    MultidimArray<fComplex> src2;

    if (preserveInput)
    {
        src2 = src;
        in = (fftwf_complex*) MULTIDIM_ARRAY(src2);
    }
    else
    {
        in = (fftwf_complex*) MULTIDIM_ARRAY(src);
    }

    fftwf_execute_dft_c2r(plan.getBackward(),
        in, MULTIDIM_ARRAY(dest));

    if (normalization == Both)
    {
        const float scale = sqrt(MULTIDIM_SIZE(dest));

        for (long int i = 0; i < NZYXSIZE(dest); i++)
        {
            dest.data[i] /= scale;
        }
    }
}



void NewFFT::FourierTransform(
    MultidimArray<double>& src,
    MultidimArray<dComplex>& dest,
    Normalization normalization)
{
    DoublePlan p(src, dest);
    FourierTransform(src, dest, p, normalization);
}

void NewFFT::inverseFourierTransform(
    MultidimArray<dComplex>& src,
    MultidimArray<double>& dest,
    Normalization normalization,
    bool preserveInput)
{
    if (preserveInput)
    {
        MultidimArray<dComplex> src2 = src;
        DoublePlan p(dest, src2);
        inverseFourierTransform(src2, dest, p, normalization, false);
    }
    else
    {
        DoublePlan p(dest, src);
        inverseFourierTransform(src, dest, p, normalization, preserveInput);
    }
}

void NewFFT::FourierTransform(
    MultidimArray<float>& src,
    MultidimArray<fComplex>& dest,
    Normalization normalization)
{
    FloatPlan p(src, dest);
    FourierTransform(src, dest, p, normalization);
}

void NewFFT::inverseFourierTransform(
    MultidimArray<fComplex>& src,
    MultidimArray<float>& dest,
    Normalization normalization,
    bool preserveInput)
{
    if (preserveInput)
    {
        MultidimArray<fComplex> src2 = src;
        FloatPlan p(dest, src2);
        inverseFourierTransform(src2, dest, p, normalization, false);
    }
    else
    {
        FloatPlan p(dest, src);
        inverseFourierTransform(src, dest, p, normalization, preserveInput);
    }
}




NewFFT::DoublePlan::DoublePlan(int w, int h, int d,
                               unsigned int flags)
:   w(w), h(h), d(d)
{
    MultidimArray<double> realDummy(d,h,w);
    MultidimArray<dComplex> complexDummy(d,h,w/2+1);

    std::vector<int> N(1,w);
    if (h > 1) N.push_back(h);
    if (d > 1) N.push_back(d);

    const int ndim = N.size();

    pthread_mutex_lock(&fftw_plan_mutex_new);

    fftw_plan planForward = fftw_plan_dft_r2c(
                ndim, &N[0],
                MULTIDIM_ARRAY(realDummy),
                (fftw_complex*) MULTIDIM_ARRAY(complexDummy),
                FFTW_UNALIGNED | flags);

    fftw_plan planBackward = fftw_plan_dft_c2r(
                ndim, &N[0],
                (fftw_complex*) MULTIDIM_ARRAY(complexDummy),
                MULTIDIM_ARRAY(realDummy),
                FFTW_UNALIGNED | flags);

    pthread_mutex_unlock(&fftw_plan_mutex_new);

    if (planForward == NULL || planBackward == NULL)
    {
        REPORT_ERROR("FFTW plans cannot be created");
    }

    plan = std::shared_ptr<Plan>(new Plan(planForward, planBackward));
}

NewFFT::DoublePlan::DoublePlan(
    MultidimArray<double>& real,
    MultidimArray<dComplex>& complex,
    unsigned int flags)
:   w(real.xdim), h(real.ydim), d(real.zdim)
{
    std::vector<int> N(1,w);
    if (h > 1) N.push_back(h);
    if (d > 1) N.push_back(d);

    const int ndim = N.size();

    pthread_mutex_lock(&fftw_plan_mutex_new);

    fftw_plan planForward = fftw_plan_dft_r2c(
                ndim, &N[0],
                MULTIDIM_ARRAY(real),
                (fftw_complex*) MULTIDIM_ARRAY(complex),
                flags);

    fftw_plan planBackward = fftw_plan_dft_c2r(
                ndim, &N[0],
                (fftw_complex*) MULTIDIM_ARRAY(complex),
                MULTIDIM_ARRAY(real),
                flags);

    pthread_mutex_unlock(&fftw_plan_mutex_new);

    if (planForward == NULL || planBackward == NULL)
    {
        REPORT_ERROR("FFTW plans cannot be created");
    }

    plan = std::shared_ptr<Plan>(new Plan(planForward, planBackward));
}

NewFFT::FloatPlan::FloatPlan(int w, int h, int d,
                             unsigned int flags)
:   w(w), h(h), d(d)
{
    MultidimArray<float> realDummy(d,h,w);
    MultidimArray<fComplex> complexDummy(d,h,w/2+1);

    std::vector<int> N(1,w);
    if (h > 1) N.push_back(h);
    if (d > 1) N.push_back(d);

    const int ndim = N.size();

    pthread_mutex_lock(&fftw_plan_mutex_new);

    fftwf_plan planForward = fftwf_plan_dft_r2c(
                ndim, &N[0],
                MULTIDIM_ARRAY(realDummy),
                (fftwf_complex*) MULTIDIM_ARRAY(complexDummy),
                FFTW_UNALIGNED | flags);

    fftwf_plan planBackward = fftwf_plan_dft_c2r(
                ndim, &N[0],
                (fftwf_complex*) MULTIDIM_ARRAY(complexDummy),
                MULTIDIM_ARRAY(realDummy),
                FFTW_UNALIGNED | flags);

    pthread_mutex_unlock(&fftw_plan_mutex_new);

    if (planForward == NULL || planBackward == NULL)
    {
        REPORT_ERROR("FFTW plans cannot be created");
    }

    plan = std::shared_ptr<Plan>(new Plan(planForward, planBackward));
}

NewFFT::FloatPlan::FloatPlan(
        MultidimArray<float>& real,
        MultidimArray<fComplex>& complex,
        unsigned int flags)
:   w(real.xdim), h(real.ydim), d(real.zdim)
{
    std::vector<int> N(1,w);
    if (h > 1) N.push_back(h);
    if (d > 1) N.push_back(d);

    const int ndim = N.size();

    pthread_mutex_lock(&fftw_plan_mutex_new);

    fftwf_plan planForward = fftwf_plan_dft_r2c(
                ndim, &N[0],
                MULTIDIM_ARRAY(real),
                (fftwf_complex*) MULTIDIM_ARRAY(complex),
                flags);

    fftwf_plan planBackward = fftwf_plan_dft_c2r(
                ndim, &N[0],
                (fftwf_complex*) MULTIDIM_ARRAY(complex),
                MULTIDIM_ARRAY(real),
                flags);

    pthread_mutex_unlock(&fftw_plan_mutex_new);

    if (planForward == NULL || planBackward == NULL)
    {
        REPORT_ERROR("FFTW plans cannot be created");
    }

    plan = std::shared_ptr<Plan>(new Plan(planForward, planBackward));
}

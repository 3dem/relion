/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/


#include "src/macros.h"
#include "parallel_ft.h"
#include "src/fftw.h"
#include "src/args.h"
#include <string.h>
#include <math.h>

//#define DEBUG_PLANS

// Constructors and destructors --------------------------------------------
ParFourierTransformer::ParFourierTransformer():
        plans_are_set(false)
{
    init();

#ifdef DEBUG_PLANS
    std::cerr << "INIT this= "<<this<< std::endl;
#endif
}

ParFourierTransformer::~ParFourierTransformer()
{
    clear();
#ifdef DEBUG_PLANS
    std::cerr << "CLEARED this= "<<this<< std::endl;
#endif
}

ParFourierTransformer::ParFourierTransformer(const ParFourierTransformer& op) :
        plans_are_set(false)
{
    // Clear current object
    clear();
    // New object is an extact copy of op
    *this = op;
}

void ParFourierTransformer::init()
{
    fReal            = NULL;
    fComplex         = NULL;
    fPlanForward     = NULL;
    fPlanBackward    = NULL;
    dataPtr          = NULL;
    complexDataPtr   = NULL;
}

void ParFourierTransformer::clear()
{
    fFourier.clear();
    // Clean-up all other FFTW-allocated memory
    destroyPlans();
    // Initialise all pointers to NULL
    init();

}

void ParFourierTransformer::cleanup()
{
    // First clear object and destroy plans
    clear();
    // Then clean up all the junk fftw keeps lying around
    // SOMEHOW THE FOLLOWING IS NOT ALLOWED WHEN USING MULTPLE TRANSFORMER OBJECTS....
#ifdef RELION_SINGLE_PRECISION
    fftwf_cleanup();
#else
    fftw_cleanup();
#endif

#ifdef DEBUG_PLANS
    std::cerr << "CLEANED-UP this= "<<this<< std::endl;
#endif

}

void ParFourierTransformer::destroyPlans()
{
    // Anything to do with plans has to be protected for threads!
    #pragma omp critical(FourierTransformer_fftw_plan)
    {
        if (plans_are_set)
        {
#ifdef RELION_SINGLE_PRECISION
            fftwf_destroy_plan(fPlanForward);
            fftwf_destroy_plan(fPlanBackward);
#else
            fftw_destroy_plan(fPlanForward);
            fftw_destroy_plan(fPlanBackward);
#endif
            plans_are_set = false;
        }
    }
}

// Initialization ----------------------------------------------------------
const MultidimArray<RFLOAT> &ParFourierTransformer::getReal() const
{
    return (*fReal);
}

const MultidimArray<Complex > &ParFourierTransformer::getComplex() const
{
    return (*fComplex);
}


void ParFourierTransformer::setReal(MultidimArray<RFLOAT> &input)
{
    bool recomputePlan=false;
    if (fReal==NULL)
        recomputePlan=true;
    /*else if (dataPtr!=MULTIDIM_ARRAY(input))
        recomputePlan=true;*/
    else
        recomputePlan=!(fReal->sameShape(input));

    fFourier.reshape(ZSIZE(input),YSIZE(input),XSIZE(input)/2+1);
    fReal=&input;

    if (recomputePlan)
    {
        int ndim=3;
        if (ZSIZE(input)==1)
        {
            ndim=2;
            if (YSIZE(input)==1)
                ndim=1;
        }
        int *N = new int[ndim];
        switch (ndim)
        {
        case 1:
            N[0]=XSIZE(input);
            break;
        case 2:
            N[0]=YSIZE(input);
            N[1]=XSIZE(input);
            break;
        case 3:
            N[0]=ZSIZE(input);
            N[1]=YSIZE(input);
            N[2]=XSIZE(input);
            break;
        }

        // Destroy both forward and backward plans if they already exist
        destroyPlans();

        // Make new plans
        plans_are_set = true;

        
        #pragma omp critical(FourierTransformer_fftw_plan)
        {
#ifdef RELION_SINGLE_PRECISION
            fPlanForward = fftwf_plan_dft_r2c(ndim, N, MULTIDIM_ARRAY(*fReal),
                                             (fftwf_complex*) MULTIDIM_ARRAY(fFourier), FFTW_ESTIMATE);
            fPlanBackward = fftwf_plan_dft_c2r(ndim, N,
                                              (fftwf_complex*) MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(*fReal),
                                              FFTW_ESTIMATE);
#else
            fPlanForward = fftw_plan_dft_r2c(ndim, N, MULTIDIM_ARRAY(*fReal),
                                             (fftw_complex*) MULTIDIM_ARRAY(fFourier), FFTW_ESTIMATE);
            fPlanBackward = fftw_plan_dft_c2r(ndim, N,
                                              (fftw_complex*) MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(*fReal),
                                          FFTW_ESTIMATE);
#endif
        }

        if (fPlanForward == NULL || fPlanBackward == NULL)
            REPORT_ERROR("FFTW plans cannot be created");

#ifdef DEBUG_PLANS
        std::cerr << " SETREAL fPlanForward= " << fPlanForward << " fPlanBackward= " << fPlanBackward  <<" this= "<<this<< std::endl;
#endif

        delete [] N;
        dataPtr=MULTIDIM_ARRAY(*fReal);
    }
}

void ParFourierTransformer::setReal(MultidimArray<Complex > &input)
{
    bool recomputePlan=false;
    if (fComplex==NULL)
        recomputePlan=true;
    /*else if (complexDataPtr!=MULTIDIM_ARRAY(input))
        recomputePlan=true;*/
    else
        recomputePlan=!(fComplex->sameShape(input));

    fFourier.resize(input);
    fComplex=&input;

    if (recomputePlan)
    {
        int ndim=3;
        if (ZSIZE(input)==1)
        {
            ndim=2;
            if (YSIZE(input)==1)
                ndim=1;
        }
        int *N = new int[ndim];
        switch (ndim)
        {
        case 1:
            N[0]=XSIZE(input);
            break;
        case 2:
            N[0]=YSIZE(input);
            N[1]=XSIZE(input);
            break;
        case 3:
            N[0]=ZSIZE(input);
            N[1]=YSIZE(input);
            N[2]=XSIZE(input);
            break;
        }

        // Destroy both forward and backward plans if they already exist
        destroyPlans();

        plans_are_set = true;

        #pragma omp critical(FourierTransformer_fftw_plan)
        {
#ifdef RELION_SINGLE_PRECISION
            fPlanForward = fftwf_plan_dft(ndim, N, (fftwf_complex*) MULTIDIM_ARRAY(*fComplex),
                                          (fftwf_complex*) MULTIDIM_ARRAY(fFourier), FFTW_FORWARD, FFTW_ESTIMATE);
            fPlanBackward = fftwf_plan_dft(ndim, N, (fftwf_complex*) MULTIDIM_ARRAY(fFourier),
                                          (fftwf_complex*) MULTIDIM_ARRAY(*fComplex), FFTW_BACKWARD, FFTW_ESTIMATE);
#else
            fPlanForward = fftw_plan_dft(ndim, N, (fftw_complex*) MULTIDIM_ARRAY(*fComplex),
                                         (fftw_complex*) MULTIDIM_ARRAY(fFourier), FFTW_FORWARD, FFTW_ESTIMATE);
            fPlanBackward = fftw_plan_dft(ndim, N, (fftw_complex*) MULTIDIM_ARRAY(fFourier),
                                          (fftw_complex*) MULTIDIM_ARRAY(*fComplex), FFTW_BACKWARD, FFTW_ESTIMATE);
#endif
        } 

        if (fPlanForward == NULL || fPlanBackward == NULL)
            REPORT_ERROR("FFTW plans cannot be created");

        delete [] N;
        complexDataPtr=MULTIDIM_ARRAY(*fComplex);
    }
}

void ParFourierTransformer::setFourier(const MultidimArray<Complex> &inputFourier)
{
    memcpy(MULTIDIM_ARRAY(fFourier),MULTIDIM_ARRAY(inputFourier),
           MULTIDIM_SIZE(inputFourier)*2*sizeof(RFLOAT));
}

// Transform ---------------------------------------------------------------
void ParFourierTransformer::Transform(int sign)
{
    if (sign == FFTW_FORWARD)
    {
#ifdef RELION_SINGLE_PRECISION
        fftwf_execute_dft_r2c(fPlanForward,MULTIDIM_ARRAY(*fReal),
                (fftwf_complex*) MULTIDIM_ARRAY(fFourier));
#else
        fftw_execute_dft_r2c(fPlanForward,MULTIDIM_ARRAY(*fReal),
                (fftw_complex*) MULTIDIM_ARRAY(fFourier));
#endif
        // Normalisation of the transform
        unsigned long int size=0;
        if(fReal!=NULL)
            size = MULTIDIM_SIZE(*fReal);
        else if (fComplex!= NULL)
            size = MULTIDIM_SIZE(*fComplex);
        else
            REPORT_ERROR("No complex nor real data defined");

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fFourier)
            DIRECT_MULTIDIM_ELEM(fFourier,n) /= size;
    }
    else if (sign == FFTW_BACKWARD)
    {
#ifdef RELION_SINGLE_PRECISION
        fftwf_execute_dft_c2r(fPlanBackward,
                (fftwf_complex*) MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(*fReal));
#else
        fftw_execute_dft_c2r(fPlanBackward,
                (fftw_complex*) MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(*fReal));
#endif
    }
}

void ParFourierTransformer::FourierTransform()
{
    Transform(FFTW_FORWARD);
}

void ParFourierTransformer::inverseFourierTransform()
{
    Transform(FFTW_BACKWARD);
}

// Enforce Hermitian symmetry ---------------------------------------------
void ParFourierTransformer::enforceHermitianSymmetry()
{
    int ndim=3;
    if (ZSIZE(*fReal)==1)
    {
        ndim=2;
        if (YSIZE(*fReal)==1)
            ndim=1;
    }
    long int yHalf=YSIZE(*fReal)/2;
    if (YSIZE(*fReal)%2==0)
        yHalf--;
    long int zHalf=ZSIZE(*fReal)/2;
    if (ZSIZE(*fReal)%2==0)
        zHalf--;
    switch (ndim)
    {
    case 2:
        for (long int i=1; i<=yHalf; i++)
        {
            long int isym=intWRAP(-i,0,YSIZE(*fReal)-1);
            Complex mean=0.5*(
                                          DIRECT_A2D_ELEM(fFourier,i,0)+
                                          conj(DIRECT_A2D_ELEM(fFourier,isym,0)));
            DIRECT_A2D_ELEM(fFourier,i,0)=mean;
            DIRECT_A2D_ELEM(fFourier,isym,0)=conj(mean);
        }
        break;
    case 3:
        for (long int k=0; k<ZSIZE(*fReal); k++)
        {
            long int ksym=intWRAP(-k,0,ZSIZE(*fReal)-1);
            for (long int i=1; i<=yHalf; i++)
            {
                long int isym=intWRAP(-i,0,YSIZE(*fReal)-1);
                Complex mean=0.5*(
                                              DIRECT_A3D_ELEM(fFourier,k,i,0)+
                                              conj(DIRECT_A3D_ELEM(fFourier,ksym,isym,0)));
                DIRECT_A3D_ELEM(fFourier,k,i,0)=mean;
                DIRECT_A3D_ELEM(fFourier,ksym,isym,0)=conj(mean);
            }
        }
        for (long int k=1; k<=zHalf; k++)
        {
            long int ksym=intWRAP(-k,0,ZSIZE(*fReal)-1);
            Complex mean=0.5*(
                                          DIRECT_A3D_ELEM(fFourier,k,0,0)+
                                          conj(DIRECT_A3D_ELEM(fFourier,ksym,0,0)));
            DIRECT_A3D_ELEM(fFourier,k,0,0)=mean;
            DIRECT_A3D_ELEM(fFourier,ksym,0,0)=conj(mean);
        }
        break;
    }
}

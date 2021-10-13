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

#ifndef NEW_FFTW_H
#define NEW_FFTW_H

#include <fftw3.h>
#include <memory>
#include <src/multidim_array.h>
#include "t_complex.h"

/*
    Usage patterns:

    #1: one-time transform

        - e.g.: call NewFFT::FourierTransform(U,V);
                for MultidimArray<double> U and MultidimArray<dComplex> V

    #2: reusable, non-array-specific plans

        (does not allow for SIMD, since the arrays that are actually transformed
         might not be memory aligned)

        - create a plan using only the size of the arrays:
          e.g.: NewFFT::DoublePlan p(512, 512);

        - plans can be copied safely:
          e.g.: NewFFT::DoublePlan p2 = p;

        - transform arbitrary arrays of that size:
          e.g.: NewFFT::FourierTransform(U,V,p2);
          (in parallel, if desired)


    #3: reusable, array-specific plans

        (allows for SIMD if the arrays are memory-aligned)

        - create a plan specific to a pair of arrays:
          e.g.: NewFFT::DoublePlan p(U,V);

        - transform those arrays using their plan:
          e.g.: NewFFT::FourierTransform(U,V,p);

*/
class NewFFT
{
    public:

        class DoublePlan;
        class FloatPlan;

        typedef enum
        {
            None,
            FwdOnly, // old default: divide by size when forward-transforming
            Both     // divide by sqrt(size) when going both ways (Parseval's theorem holds)
        }
        Normalization;

        /* Four transforms using reusable plans
           (forward and backward for doubles and floats).

           NOTE: reusable plans created with the Double/FloatPlan(w,h,d)
           constructor always assume that the arrays
           are not memory-aligned - SIMD is thus disabled.

           The output array will be rescaled if and only if
           it does not have the correct size already.
        */
        static void FourierTransform(
                MultidimArray<double>& src,
                MultidimArray<dComplex>& dest,
                const DoublePlan& plan,
                Normalization normalization = FwdOnly);

        /* FFTW always destroys the input Complex-array in an inverse-FFT.
           If preserveInput is set, a copy will be made prior to transforming. */
        static void inverseFourierTransform(
                MultidimArray<dComplex>& src,
                MultidimArray<double>& dest,
                const DoublePlan& plan,
                Normalization normalization = FwdOnly,
                bool preserveInput = true);


        static void FourierTransform(
                MultidimArray<float>& src,
                MultidimArray<fComplex>& dest,
                const FloatPlan& plan,
                Normalization normalization = FwdOnly);

        static void inverseFourierTransform(
                MultidimArray<fComplex>& src,
                MultidimArray<float>& dest,
                const FloatPlan& plan,
                Normalization normalization = FwdOnly,
                bool preserveInput = true);


        // Four transforms using array-specific plans.
        // The plan will be created ad-hoc.
        // If the two arrays are memory-aligned,
        // SIMD will be used (if available)
        static void FourierTransform(
                MultidimArray<double>& src,
                MultidimArray<dComplex>& dest,
                Normalization normalization = FwdOnly);

        static void inverseFourierTransform(
                MultidimArray<dComplex>& src,
                MultidimArray<double>& dest,
                Normalization normalization = FwdOnly,
                bool preserveInput = true);


        static void FourierTransform(
                MultidimArray<float>& src,
                MultidimArray<fComplex>& dest,
                Normalization normalization = FwdOnly);

        static void inverseFourierTransform(
                MultidimArray<fComplex>& src,
                MultidimArray<float>& dest,
                Normalization normalization = FwdOnly,
                bool preserveInput = true);
	
		
		template<class T>
		static bool areSizesCompatible(
				const MultidimArray<T>& real,
				const MultidimArray<tComplex<T> >& complex)
		{
			return real.xdim == 2 * (complex.xdim - 1)
			    && real.ydim == complex.ydim
			    && real.zdim == complex.zdim
			    && real.ndim == complex.ndim;
		}
		
		template<class T>
		static void resizeRealToMatch(
				MultidimArray<T>& real,
				const MultidimArray<tComplex<T> >& complex)
		{
			real.resizeNoCp(complex.ndim, complex.zdim, complex.ydim, 2 * (complex.xdim - 1));
		}
		
		template<class T>
		static void resizeComplexToMatch(
				const MultidimArray<T>& real,
				MultidimArray<tComplex<T> >& complex)
		{
			complex.resizeNoCp(real.ndim, real.zdim, real.ydim, real.xdim/2 + 1);
		}
	
		
	private:
		
		/* Private static methods that perform the actual transforms.
		   The arrays are guaranteed to have the correct sizes 
		   and a compatible plan when these are called. 
		   Also, the complex array is always destroyed in the 
		   inverse transform - a copy has been made before.*/
		
		static void _FourierTransform(
                MultidimArray<double>& src,
                MultidimArray<dComplex>& dest,
                const DoublePlan& plan,
                Normalization normalization);

        static void _inverseFourierTransform(
                MultidimArray<dComplex>& src,
                MultidimArray<double>& dest,
                const DoublePlan& plan,
                Normalization normalization);

        static void _FourierTransform(
                MultidimArray<float>& src,
                MultidimArray<fComplex>& dest,
                const FloatPlan& plan,
                Normalization normalization);

        static void _inverseFourierTransform(
                MultidimArray<fComplex>& src,
                MultidimArray<float>& dest,
                const FloatPlan& plan,
                Normalization normalization);
		
	public:

        /* These plan classes can be copied freely.
           The corresponding pairs of fftw_plan instances
           will be automatically deallocated using fftw_destroy_plan()
           when no instance of Double/FloatPlan points to them.
           (note the std::shared_ptr<Plan> 'plan') */
        class DoublePlan
        {
            public:

                /* Constructor for reusable plans.

                   'w', 'h' and 'd' are the dimensions of the *real*-array
                       to be transformed. The corresponding complex array has
                       the dimensions (w/2+1)*h*d

                   'flags' allows for controlling planning rigor and
                       setting algorithmic restrictions.
                       (cf. http://www.fftw.org/fftw3_doc/Planner-Flags.html)
                */
                DoublePlan(int w, int h = 1, int d = 1,
                           unsigned int flags = FFTW_ESTIMATE);

                // constructor for array-specific plans
                DoublePlan( MultidimArray<double>& real,
                            MultidimArray<dComplex>& complex,
                            unsigned int flags = FFTW_ESTIMATE);

                fftw_plan getForward() const
                {
                    return plan.get()->forward;
                }

                fftw_plan getBackward() const
                {
                    return plan.get()->backward;
                }

                bool isCompatible(const MultidimArray<double>& real) const
                {
                    return real.xdim == w && real.ydim == h && real.zdim == d
							&& (reusable || realPtr == MULTIDIM_ARRAY(real));
                }

                bool isCompatible(const MultidimArray<dComplex>& complex) const
                {
                    return complex.xdim == w/2+1 && complex.ydim == h && complex.zdim == d							
							&& (reusable || complexPtr == (double*)MULTIDIM_ARRAY(complex));
                }
				
		bool isReusable() const
		{
			return reusable;
		}

            private:

                class Plan
                {
                    public:

                        Plan(fftw_plan forward, fftw_plan backward)
                        :   forward(forward), backward(backward)
                        {}

                        ~Plan()
                        {
                            #pragma omp critical(FourierTransformer_fftw_plan) 
                            {
                                fftw_destroy_plan(forward);
                                fftw_destroy_plan(backward);
                            }
                        }

                        fftw_plan forward, backward;
                };

				bool reusable;
                int w, h, d;
				double *realPtr, *complexPtr;
                std::shared_ptr<Plan> plan;
        };

        class FloatPlan
        {
            public:

                FloatPlan(int w, int h = 1, int d = 1,
                          unsigned int flags = FFTW_ESTIMATE);

                FloatPlan(MultidimArray<float>& real,
                          MultidimArray<fComplex>& complex,
                          unsigned int flags = FFTW_ESTIMATE);

                fftwf_plan getForward() const
                {
                    return plan.get()->forward;
                }

                fftwf_plan getBackward() const
                {
                    return plan.get()->backward;
                }

                bool isCompatible(const MultidimArray<float>& real) const
                {
                    return (real.xdim == w && real.ydim == h && real.zdim == d)
							&& (reusable || realPtr == MULTIDIM_ARRAY(real));
                }

                bool isCompatible(const MultidimArray<fComplex>& complex) const
                {
                    return (complex.xdim == w/2+1 && complex.ydim == h && complex.zdim == d)
							&& (reusable || complexPtr == (float*)MULTIDIM_ARRAY(complex));
                }
				
		bool isReusable() const
		{
			return reusable;
		}

            private:

                class Plan
                {
                	public:

                        Plan(fftwf_plan forward, fftwf_plan backward)
                        :   forward(forward), backward(backward)
                        {}

                        ~Plan()
                        {
                            #pragma omp critical(FourierTransformer_fftw_plan)
                            {
                                fftwf_destroy_plan(forward);
                                fftwf_destroy_plan(backward);
                            }
                        }

                        fftwf_plan forward, backward;
                };
				
				bool reusable;
                int w, h, d;
				float *realPtr, *complexPtr;
                std::shared_ptr<Plan> plan;
        };

};

// This is to get NewFFTPlan::Plan<RFLOAT>
template <typename T> struct NewFFTPlan {};

template <> struct NewFFTPlan<double>
{
	typedef NewFFT::DoublePlan type;
};

template <> struct NewFFTPlan<float>
{
	typedef NewFFT::FloatPlan type;
};

#endif

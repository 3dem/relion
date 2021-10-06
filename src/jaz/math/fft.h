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

#ifndef VERY_NEW_FFTW_H
#define VERY_NEW_FFTW_H

#include <fftw3.h>
#include <omp.h>
#include <memory>
#include <src/jaz/image/buffered_image.h>
#include "t_complex.h"

/*
	Usage patterns:

	#1: one-time transform

		- e.g.: call FFT::FourierTransform(U,V);
				for Image<double> U and Image<dComplex> V

	#2: reusable, non-array-specific plans

		(does not allow for SIMD, since the arrays that are actually transformed
		 might not be memory aligned)

		- create a plan using only the size of the arrays:
		  e.g.: FFT::DoublePlan p(512, 512);

		- plans can be copied safely:
		  e.g.: FFT::DoublePlan p2 = p;

		- transform arbitrary arrays of that size:
		  e.g.: FFT::FourierTransform(U,V,p2);
		  (in parallel, if desired)


	#3: reusable, array-specific plans

		(allows for SIMD if the arrays are memory-aligned)

		- create a plan specific to a pair of arrays:
		  e.g.: FFT::DoublePlan p(U,V);

		- transform those arrays using their plan:
		  e.g.: FFT::FourierTransform(U,V,p);

*/
class FFT
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
				BufferedImage<double>& src,
				BufferedImage<dComplex>& dest,
				const DoublePlan& plan,
				Normalization normalization = FwdOnly);

		/* FFTW always destroys the input Complex-array in an inverse-FFT.
		   If preserveInput is set, a copy will be made prior to transforming. */
		static void inverseFourierTransform(
				BufferedImage<dComplex>& src,
				BufferedImage<double>& dest,
				const DoublePlan& plan,
				Normalization normalization = FwdOnly,
				bool preserveInput = true);


		static void FourierTransform(
				BufferedImage<float>& src,
				BufferedImage<fComplex>& dest,
				const FloatPlan& plan,
				Normalization normalization = FwdOnly);

		static void inverseFourierTransform(
				BufferedImage<fComplex>& src,
				BufferedImage<float>& dest,
				const FloatPlan& plan,
				Normalization normalization = FwdOnly,
				bool preserveInput = true);


		// Four transforms using array-specific plans.
		// The plan will be created ad-hoc.
		// If the two arrays are memory-aligned,
		// SIMD will be used (if available)
		static void FourierTransform(
				BufferedImage<double>& src,
				BufferedImage<dComplex>& dest,
				Normalization normalization = FwdOnly);

		static void inverseFourierTransform(
				BufferedImage<dComplex>& src,
				BufferedImage<double>& dest,
				Normalization normalization = FwdOnly,
				bool preserveInput = true);


		static void FourierTransform(
				BufferedImage<float>& src,
				BufferedImage<fComplex>& dest,
				Normalization normalization = FwdOnly);

		static void inverseFourierTransform(
				BufferedImage<fComplex>& src,
				BufferedImage<float>& dest,
				Normalization normalization = FwdOnly,
				bool preserveInput = true);

		
		template<class T>
		static bool areSizesCompatible(
				const BufferedImage<T>& real,
				const BufferedImage<tComplex<T>>& complex)
		{
			return real.xdim == 2 * (complex.xdim - 1)
					&& real.ydim == complex.ydim
					&& real.zdim == complex.zdim;
		}
		
		template<class T>
		static void resizeRealToMatch(
				BufferedImage<T>& real,
				const BufferedImage<tComplex<T>>& complex)
		{
			real.resize(2 * (complex.xdim - 1), complex.ydim, complex.zdim);
		}
		
		template<class T>
		static void resizeComplexToMatch(
				const BufferedImage<T>& real,
				BufferedImage<tComplex<T>>& complex)
		{
			complex.resize(real.xdim/2 + 1, real.ydim, real.zdim);
		}
		
		template<class T>
		static BufferedImage<tComplex<T>> toComplex(const BufferedImage<T>& FourierSpaceReal)
		{
			BufferedImage<tComplex<T>> out(
						FourierSpaceReal.xdim,
						FourierSpaceReal.ydim,
						FourierSpaceReal.zdim);

			const long long int size =
					FourierSpaceReal.xdim *
					FourierSpaceReal.ydim *
					FourierSpaceReal.zdim;
			
			for (long long int i = 0; i < size; i++)
			{
				out[i] = tComplex<T>(FourierSpaceReal[i], T(0));
			}
			
			return out;
		}
		
		template<class T>
		static BufferedImage<T> toReal(const BufferedImage<tComplex<T>>& FourierSpaceComplex)
		{
			BufferedImage<T> out(
						FourierSpaceComplex.xdim,
						FourierSpaceComplex.ydim,
						FourierSpaceComplex.zdim);

			const long long int size =
					FourierSpaceComplex.xdim *
					FourierSpaceComplex.ydim *
					FourierSpaceComplex.zdim;
			
			for (long long int i = 0; i < size; i++)
			{
				out[i] = FourierSpaceComplex[i].real;
			}
			
			return out;
		}
		
		
		// syntactic sugar (the RawImages are copied to BufferedImages internally)
		
		static void FourierTransform(
				const RawImage<double>& src,
				BufferedImage<dComplex>& dest,
				Normalization normalization = FwdOnly);

		static void inverseFourierTransform(
				const RawImage<dComplex>& src,
				BufferedImage<double>& dest,
				Normalization normalization = FwdOnly);


		static void FourierTransform(
				const RawImage<float>& src,
				BufferedImage<fComplex>& dest,
				Normalization normalization = FwdOnly);

		static void inverseFourierTransform(
				const RawImage<fComplex>& src,
				BufferedImage<float>& dest,
				Normalization normalization = FwdOnly);

		
	private:
		
		/* Private static methods that perform the actual transforms.
		   The arrays are guaranteed to have the correct sizes
		   and a compatible plan when these are called.
		   Also, the complex array is always destroyed in the
		   inverse transform - a copy has been made before.*/
		
		static void _FourierTransform(
				BufferedImage<double>& src,
				BufferedImage<dComplex>& dest,
				const DoublePlan& plan,
				Normalization normalization);

		static void _inverseFourierTransform(
				BufferedImage<dComplex>& src,
				BufferedImage<double>& dest,
				const DoublePlan& plan,
				Normalization normalization);

		static void _FourierTransform(
				BufferedImage<float>& src,
				BufferedImage<fComplex>& dest,
				const FloatPlan& plan,
				Normalization normalization);

		static void _inverseFourierTransform(
				BufferedImage<fComplex>& src,
				BufferedImage<float>& dest,
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
				DoublePlan( BufferedImage<double>& real,
							BufferedImage<dComplex>& complex,
							unsigned int flags = FFTW_ESTIMATE);

				DoublePlan(){}


				fftw_plan getForward() const
				{
					return plan.get()->forward;
				}

				fftw_plan getBackward() const
				{
					return plan.get()->backward;
				}

				bool isCompatible(const BufferedImage<double>& real) const
				{
					return real.xdim == w && real.ydim == h && real.zdim == d
							&& (reusable || realPtr == real.getData());
				}

				bool isCompatible(const BufferedImage<dComplex>& complex) const
				{
					return complex.xdim == w/2+1 && complex.ydim == h && complex.zdim == d
							&& (reusable || complexPtr == (double*)complex.getData());
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

				FloatPlan(BufferedImage<float>& real,
						  BufferedImage<fComplex>& complex,
						  unsigned int flags = FFTW_ESTIMATE);

				FloatPlan(){}


				fftwf_plan getForward() const
				{
					return plan.get()->forward;
				}

				fftwf_plan getBackward() const
				{
					return plan.get()->backward;
				}

				bool isCompatible(const BufferedImage<float>& real) const
				{
					return (real.xdim == w && real.ydim == h && real.zdim == d)
							&& (reusable || realPtr == real.getData());
				}

				bool isCompatible(const BufferedImage<fComplex>& complex) const
				{
					return (complex.xdim == w/2+1 && complex.ydim == h && complex.zdim == d)
							&& (reusable || complexPtr == (float*)complex.getData());
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

#endif

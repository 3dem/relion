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

#include "fft.h"

#include <src/macros.h>
//#include "fftw.h"
//#include "args.h"
#include <string.h>
#include <math.h>

void FFT::FourierTransform(
		BufferedImage<double>& src,
		BufferedImage<dComplex>& dest,
		const FFT::DoublePlan& plan,
		Normalization normalization)
{
	if (!plan.isCompatible(src))
	{
		REPORT_ERROR("ERROR: plan incompatible with input array\n");
	}
	
	if (!plan.isCompatible(dest))
	{
		if (plan.isReusable())
		{
			dest.resize(src.xdim/2 + 1, src.ydim, src.zdim);
		}
		else
		{
			REPORT_ERROR("FFT::FourierTransform: plan incompatible with output array\n");
		}
	}
	
	_FourierTransform(src, dest, plan, normalization);
}

void FFT::inverseFourierTransform(
		BufferedImage<dComplex>& src,
		BufferedImage<double>& dest,
		const FFT::DoublePlan& plan,
		Normalization normalization,
		bool preserveInput)
{
	if (preserveInput && !plan.isReusable())
	{
		REPORT_ERROR("FFT::inverseFourierTransform: preserveInput is only supported for reusable plans\n");
	}
	
	if (!plan.isCompatible(src))
	{
		REPORT_ERROR("FFT::inverseFourierTransform: plan incompatible with input array\n");
	}
	
	if (!plan.isCompatible(dest))
	{
		if (plan.isReusable())
		{
			dest.resize(2*(src.xdim-1), src.ydim, src.zdim);
		}
		else
		{
			REPORT_ERROR("FFT::inverseFourierTransform: plan incompatible with output array\n");
		}
	}
	
	BufferedImage<dComplex> src2;
	
	if (preserveInput)
	{
		src2 = src;
		_inverseFourierTransform(src2, dest, plan, normalization);
	}
	else
	{
		_inverseFourierTransform(src, dest, plan, normalization);
	}
}

void FFT::FourierTransform(
		BufferedImage<float>& src,
		BufferedImage<fComplex>& dest,
		const FFT::FloatPlan& plan,
		Normalization normalization)
{
	if (!plan.isCompatible(src))
	{
		REPORT_ERROR("FFT::FourierTransform: plan incompatible with input array\n");
	}
	
	if (!plan.isCompatible(dest))
	{
		if (plan.isReusable())
		{
			dest.resize(src.xdim/2 + 1, src.ydim, src.zdim);
		}
		else
		{
			REPORT_ERROR("FFT::FourierTransform: plan incompatible with output array\n");
		}
	}
	
	_FourierTransform(src, dest, plan, normalization);
}

void FFT::inverseFourierTransform(
		BufferedImage<fComplex>& src,
		BufferedImage<float>& dest,
		const FFT::FloatPlan& plan,
		Normalization normalization,
		bool preserveInput)
{
	if (preserveInput && !plan.isReusable())
	{
		REPORT_ERROR("FFT::inverseFourierTransform: preserveInput is only supported for reusable plans\n");
	}
	
	if (!plan.isCompatible(src))
	{
		REPORT_ERROR("FFT::inverseFourierTransform: plan incompatible with input array\n");
	}
	
	if (!plan.isCompatible(dest))
	{
		if (plan.isReusable())
		{
			dest.resize(2*(src.xdim-1), src.ydim, src.zdim);
		}
		else
		{
			REPORT_ERROR("FFT::inverseFourierTransform: plan incompatible with output array\n");
		}
	}
	
	BufferedImage<fComplex> src2;
	
	if (preserveInput)
	{
		src2 = src;
		_inverseFourierTransform(src2, dest, plan, normalization);
	}
	else
	{
		_inverseFourierTransform(src, dest, plan, normalization);
	}
}



void FFT::FourierTransform(
		BufferedImage<double>& src,
		BufferedImage<dComplex>& dest,
		Normalization normalization)
{
	if (!areSizesCompatible(src, dest))
	{
		resizeComplexToMatch(src, dest);
	}
			
	DoublePlan p(src, dest);
	_FourierTransform(src, dest, p, normalization);
}

void FFT::inverseFourierTransform(
		BufferedImage<dComplex>& src,
		BufferedImage<double>& dest,
		Normalization normalization,
		bool preserveInput)
{
	if (!areSizesCompatible(dest, src))
	{
		resizeRealToMatch(dest, src);
	}
			
	if (preserveInput)
	{
		BufferedImage<dComplex> src2 = src;
		DoublePlan p(dest, src2);
		_inverseFourierTransform(src2, dest, p, normalization);
	}
	else
	{
		DoublePlan p(dest, src);
		_inverseFourierTransform(src, dest, p, normalization);
	}
}

void FFT::FourierTransform(
		BufferedImage<float>& src,
		BufferedImage<fComplex>& dest,
		Normalization normalization)
{
	if (!areSizesCompatible(src, dest))
	{
		resizeComplexToMatch(src, dest);
	}
			
	FloatPlan p(src, dest);
	_FourierTransform(src, dest, p, normalization);
}

void FFT::inverseFourierTransform(
		BufferedImage<fComplex>& src,
		BufferedImage<float>& dest,
		Normalization normalization,
		bool preserveInput)
{
	if (!areSizesCompatible(dest, src))
	{
		resizeRealToMatch(dest, src);
	}
			
	if (preserveInput)
	{
		BufferedImage<fComplex> src2 = src;
		FloatPlan p(dest, src2);
		_inverseFourierTransform(src2, dest, p, normalization);
	}
	else
	{
		FloatPlan p(dest, src);
		_inverseFourierTransform(src, dest, p, normalization);
	}
}

void FFT::FourierTransform(
        const RawImage<double>& src,
        BufferedImage<dComplex>& dest,
        Normalization normalization)
{
	BufferedImage<double> copy(src);
	FourierTransform(copy, dest, normalization);
}

void FFT::inverseFourierTransform(
        const RawImage<dComplex>& src,
        BufferedImage<double>& dest,
        Normalization normalization)
{
	BufferedImage<dComplex> copy(src);
	inverseFourierTransform(copy, dest, normalization, false);
}

void FFT::FourierTransform(
        const RawImage<float>& src,
        BufferedImage<fComplex>& dest,
        Normalization normalization)
{
	BufferedImage<float> copy(src);
	FourierTransform(copy, dest, normalization);
}

void FFT::inverseFourierTransform(
        const RawImage<fComplex>& src,
        BufferedImage<float>& dest,
        Normalization normalization)
{
	BufferedImage<fComplex> copy(src);
	inverseFourierTransform(copy, dest, normalization, false);
}

void FFT::_FourierTransform(
		BufferedImage<double>& src,
		BufferedImage<dComplex>& dest,
		const FFT::DoublePlan& plan,
		Normalization normalization)
{
	fftw_execute_dft_r2c(
		plan.getForward(),
		src.getData(), 
		(fftw_complex*) dest.getData());
	
	if (normalization == FwdOnly)
	{
		const double scale = src.getSize();
		
		for (long int i = 0; i < dest.getSize(); i++)
		{
			dest.data[i] /= scale;
		}
	}
	else if (normalization == Both)
	{
		const double scale = sqrt(src.getSize());
		
		for (long int i = 0; i < dest.getSize(); i++)
		{
			dest.data[i] /= scale;
		}
	}
}

void FFT::_inverseFourierTransform(
		BufferedImage<dComplex>& src,
		BufferedImage<double>& dest,
		const FFT::DoublePlan& plan,
		Normalization normalization)
{
	fftw_complex* in = (fftw_complex*) src.getData();	
	fftw_execute_dft_c2r(plan.getBackward(), in, dest.getData());
	
	if (normalization == Both)
	{
		const double scale = sqrt(dest.getSize());
		
		for (long int i = 0; i < dest.getSize(); i++)
		{
			dest[i] /= scale;
		}
	}
}

void FFT::_FourierTransform(
		BufferedImage<float>& src,
		BufferedImage<fComplex>& dest,
		const FFT::FloatPlan& plan,
		Normalization normalization)
{
	fftwf_execute_dft_r2c(
		plan.getForward(),
		src.getData(), 
		(fftwf_complex*) dest.getData());
	
	if (normalization == FwdOnly)
	{
		const float scale = src.getSize();
		
		for (long int i = 0; i < dest.getSize(); i++)
		{
			dest.data[i] /= scale;
		}
	}
	else if (normalization == Both)
	{
		const float scale = sqrt(src.getSize());
		
		for (long int i = 0; i < dest.getSize(); i++)
		{
			dest.data[i] /= scale;
		}
	}
}

void FFT::_inverseFourierTransform(
		BufferedImage<fComplex>& src,
		BufferedImage<float>& dest,
		const FFT::FloatPlan& plan,
		Normalization normalization)
{	
	fftwf_complex* in = (fftwf_complex*) src.getData();	
	fftwf_execute_dft_c2r(plan.getBackward(), in, dest.getData());
	
	if (normalization == Both)
	{
		const float scale = sqrt(dest.getSize());
		
		for (long int i = 0; i < dest.getSize(); i++)
		{
			dest.data[i] /= scale;
		}
	}
}


FFT::DoublePlan::DoublePlan(int w, int h, int d, unsigned int flags)
:   
	reusable(true), 
	w(w), h(h), d(d),
	realPtr(0), 
	complexPtr(0)
{
	BufferedImage<double> realDummy(w,h,d);
	BufferedImage<dComplex> complexDummy(w/2+1,h,d);
	
	std::vector<int> N(0);
	if (d > 1) N.push_back(d);
	if (h > 1) N.push_back(h);
	           N.push_back(w);
	
	const int ndim = N.size();

	fftw_plan planForward, planBackward;	
	#pragma omp critical(FourierTransformer_fftw_plan)
	{
		planForward = fftw_plan_dft_r2c(
				ndim, &N[0],
				realDummy.getData(),
				(fftw_complex*) complexDummy.getData(),
				FFTW_UNALIGNED | flags);
		
		planBackward = fftw_plan_dft_c2r(
				ndim, &N[0],
				(fftw_complex*) complexDummy.getData(),
				realDummy.getData(),
				FFTW_UNALIGNED | flags);
	}
	
	if (planForward == NULL || planBackward == NULL)
	{
		REPORT_ERROR("FFTW plans cannot be created");
	}
	
	plan = std::shared_ptr<Plan>(new Plan(planForward, planBackward));
}

FFT::DoublePlan::DoublePlan(
		BufferedImage<double>& real,
		BufferedImage<dComplex>& complex,
		unsigned int flags)
:   
	reusable(flags & FFTW_UNALIGNED), 
	w(real.xdim), h(real.ydim), d(real.zdim),
	realPtr(real.getData()), 
	complexPtr((double*)complex.getData())
{
	std::vector<int> N(0);
	if (d > 1) N.push_back(d);
	if (h > 1) N.push_back(h);
	           N.push_back(w);
	
	const int ndim = N.size();

	fftw_plan planForward, planBackward;	
	#pragma omp critical(FourierTransformer_fftw_plan)
	{	
		planForward = fftw_plan_dft_r2c(
				ndim, &N[0],
				real.getData(),
				(fftw_complex*) complex.getData(),
				flags);
		
		planBackward = fftw_plan_dft_c2r(
				ndim, &N[0],
				(fftw_complex*) complex.getData(),
				real.getData(),
				flags);
	}
	
	if (planForward == NULL || planBackward == NULL)
	{
		REPORT_ERROR("FFTW plans cannot be created");
	}
	
	plan = std::shared_ptr<Plan>(new Plan(planForward, planBackward));
}

FFT::FloatPlan::FloatPlan(int w, int h, int d, unsigned int flags)
:	reusable(true), 
	w(w), h(h), d(d),
	realPtr(0), 
	complexPtr(0)
{
	BufferedImage<float> realDummy(w,h,d);
	BufferedImage<fComplex> complexDummy(w/2+1,h,d);
	
	std::vector<int> N(0);
	if (d > 1) N.push_back(d);
	if (h > 1) N.push_back(h);
	           N.push_back(w);
	
	const int ndim = N.size();

	fftwf_plan planForward, planBackward;	
	#pragma omp critical(FourierTransformer_fftw_plan)
	{	
		planForward = fftwf_plan_dft_r2c(
				ndim, &N[0],
				realDummy.getData(),
				(fftwf_complex*) complexDummy.getData(),
				FFTW_UNALIGNED | flags);
		
		planBackward = fftwf_plan_dft_c2r(
				ndim, &N[0],
				(fftwf_complex*) complexDummy.getData(),
				realDummy.getData(),
				FFTW_UNALIGNED | flags);
	}
	
	if (planForward == NULL || planBackward == NULL)
	{
		REPORT_ERROR("FFTW plans cannot be created");
	}
	
	plan = std::shared_ptr<Plan>(new Plan(planForward, planBackward));
}

FFT::FloatPlan::FloatPlan(
		BufferedImage<float>& real,
		BufferedImage<fComplex>& complex,
		unsigned int flags)
:
	reusable(flags & FFTW_UNALIGNED), 
	w(real.xdim), h(real.ydim), d(real.zdim),
	realPtr(real.getData()), 
	complexPtr((float*)complex.getData())
{
	std::vector<int> N(0);
	if (d > 1) N.push_back(d);
	if (h > 1) N.push_back(h);
	           N.push_back(w);
	
	const int ndim = N.size();

	fftwf_plan planForward, planBackward;	
	#pragma omp critical(FourierTransformer_fftw_plan)
	{	
		planForward = fftwf_plan_dft_r2c(
				ndim, &N[0],
				real.getData(),
				(fftwf_complex*) complex.getData(),
				flags);
		
		planBackward = fftwf_plan_dft_c2r(
				ndim, &N[0],
				(fftwf_complex*) complex.getData(),
				real.getData(),
				flags);
	}	
	
	if (planForward == NULL || planBackward == NULL)
	{
		REPORT_ERROR("FFTW plans cannot be created");
	}
	
	plan = std::shared_ptr<Plan>(new Plan(planForward, planBackward));
}

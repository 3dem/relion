#ifndef NEW_STACK_HELPER_H
#define NEW_STACK_HELPER_H

#include "buffered_image.h"
#include <src/jaz/math/t_complex.h>
#include <src/jaz/math/fft.h>
#include <vector>
#include <omp.h>


class NewStackHelper
{
	public:
		
		template <typename T>
		static BufferedImage<T> extractSliceZ(const RawImage<T>& img, int z);
		
		template <typename SrcType, typename DestType>
		static void insertSliceZ(const RawImage<SrcType>& slice, RawImage<DestType>& dest, int z);
		
		template <typename T>
		static void writeAsStack(const std::vector<RawImage<T>>& vec, std::string fn);
		
		template <typename T>
		static void writeAsStack(const std::vector<BufferedImage<T>>& vec, std::string fn);
		

		template <typename T>
		static void FourierTransformStack(
				const RawImage<T>& stack,
				RawImage<tComplex<T>>& stackOut,
				bool center = true,
				int num_threads = 1);

		static void FourierTransformStack_fast(
				const RawImage<float>& stack,
				RawImage<fComplex>& stackOut,
				bool center = true,
				int num_threads = 1);
		
		template <typename T>
		static BufferedImage<tComplex<T>> FourierTransformStack(
				const RawImage<T>& stack, 
				bool center = true,
				int num_threads = 1);
		
		template <typename T>
		static void inverseFourierTransformStack(
				const RawImage<tComplex<T>>& stack, 
				RawImage<T>& stackOut, 
				bool decenter = true,
				int num_threads = 1);
		
		template <typename T>
		static BufferedImage<T> inverseFourierTransformStack(
				const RawImage<tComplex<T>>& stack,  
				bool decenter = true,
				int num_threads = 1);
		
		template <typename T>
		static void shiftStack(
				const RawImage<tComplex<T>>& stack,
				const std::vector<gravis::d2Vector>& pos,
				RawImage<tComplex<T>>& stackOut,
				bool centered = true,
				int num_threads = 1); 
};

template<typename T>
BufferedImage<T> NewStackHelper::extractSliceZ(const RawImage<T>& img, int z)
{
	const int w = img.xdim;
	const int h = img.ydim;
	const int d = img.zdim;
	
	if (z < 0 || z >= d)
	{
		REPORT_ERROR_STR("NewStackHelper::extractSliceZ: unable to extract slice index " << z
						 << " from " << d << " slices.\n");
	}
	
	BufferedImage<T> out(w,h);
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		out(x,y) = img(x,y,z);
	}
	
	return out;
}

template<typename SrcType, typename DestType>
void NewStackHelper::insertSliceZ(const RawImage<SrcType> &slice, RawImage<DestType> &dest, int z)
{
	const int w = dest.xdim;
	const int h = dest.ydim;
	const int d = dest.zdim;
	
	if (z < 0 || z >= d)
	{
		REPORT_ERROR_STR("NewStackHelper::insertSliceZ: unable to insert slice index " << z
						 << " into " << d << " slices.\n");
	}
	
	if (slice.xdim != w || slice.ydim != h)
	{
		REPORT_ERROR_STR("NewStackHelper::insertSliceZ: image size mismatch - slice: "
						 << slice.xdim << "x" << slice.ydim << ", dest: " 
						 << w << "x" << h << "x" << d << ".\n");
	}
	
	for (int y = 0; y < h; y++)
	for (int x = 0; x < w; x++)
	{
		dest(x,y,z) = (DestType) slice(x,y);
	}
}

template<typename T>
void NewStackHelper::writeAsStack(const std::vector<RawImage<T>>& vec, std::string fn)
{
	const int fc = vec.size();
	
	if (fc < 1) 
	{
		REPORT_ERROR_STR("NewStackHelper::writeAsStack: unable to write " << fn << " - empty stack.");
	}
	
	const size_t w = vec[0].xdim;
	const size_t h = vec[0].ydim;
	
	BufferedImage<T> out(w,h,fc);
	
	for (int f = 0; f < fc; f++)
	{
		insertSliceZ(vec[f], out, f);
	}
	
	out.write(fn);
}

template<typename T>
void NewStackHelper::writeAsStack(const std::vector<BufferedImage<T>>& vec, std::string fn)
{
	const int fc = vec.size();
	
	if (fc < 1) 
	{
		REPORT_ERROR_STR("NewStackHelper::writeAsStack: unable to write " << fn << " - empty stack.");
	}
	
	const size_t w = vec[0].xdim;
	const size_t h = vec[0].ydim;
	
	BufferedImage<T> out(w,h,fc);
	
	for (int f = 0; f < fc; f++)
	{
		insertSliceZ(vec[f], out, f);
	}
	
	out.write(fn);
}

template <typename T>
void NewStackHelper::FourierTransformStack(
		const RawImage<T>& stack,
		RawImage<tComplex<T>>& stackOut,
		bool center,
		int num_threads)
{
	const int w = stack.xdim;
	const int h = stack.ydim;
	const int fc = stack.zdim;

	const int wh = w/2 + 1;

	std::vector<BufferedImage<T>> tempRS(num_threads);
	std::vector<BufferedImage<tComplex<T>>> tempFS(num_threads);

	for (int t = 0; t < num_threads; t++)
	{
		tempRS[t] = BufferedImage<T>(w,h);
		tempFS[t] = BufferedImage<tComplex<T>>(wh,h);
	}

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();

		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < w; x++)
		{
			tempRS[t](x,y) = center? stack((x+w/2)%w, (y+h/2)%h, f) : stack(x,y,f);
		}

		FFT::FourierTransform(tempRS[t], tempFS[t], FFT::Both);

		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < wh; x++)
		{
			stackOut(x,y,f) = tempFS[t](x,y);
		}
	}
}

template <typename T>
BufferedImage<tComplex<T>> NewStackHelper::FourierTransformStack(
		const RawImage<T>& stack, 
		bool center,
		int num_threads)
{
	BufferedImage<tComplex<T>> out(stack.xdim/2 + 1, stack.ydim, stack.zdim);
	FourierTransformStack(stack, out, center, num_threads);
	return out;
}

template <typename T>
void NewStackHelper::inverseFourierTransformStack(
		const RawImage<tComplex<T>>& stack, 
		RawImage<T>& stackOut, 
		bool decenter,
		int num_threads)
{
	const int wh = stack.xdim;
	const int h = stack.ydim;
	const int fc = stack.zdim;
	
	const int w = (wh - 1) * 2;
	
	if (stackOut.xdim != w || stackOut.ydim != h || stackOut.zdim != fc)
	{
		REPORT_ERROR_STR("Extraction::inverseFourierTransformStack: wrong output size: "
						 << stackOut.xdim << "x" << stackOut.ydim << "x" << stackOut.zdim
						 << " (expected: " 
						 << w << "x" << h << "x" << fc << ")");
	}
		
	std::vector<BufferedImage<T>> tempRS(num_threads);
	std::vector<BufferedImage<tComplex<T>>> tempFS(num_threads);
	
	for (int t = 0; t < num_threads; t++)
	{
		tempRS[t] = BufferedImage<T>(w,h);
		tempFS[t] = BufferedImage<tComplex<T>>(wh,h);
	}
	
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		const int t = omp_get_thread_num();
		
		for (long int y = 0; y < h;  y++)
		for (long int x = 0; x < wh; x++)
		{
			const int mod = decenter? ((1 - 2 * (x % 2)) * (1 - 2 * (y % 2))) : 1;
			tempFS[t](x,y) = mod * stack(x,y,f);
		}
				
		FFT::inverseFourierTransform(tempFS[t], tempRS[t], FFT::Both);
		
		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < w; x++)
		{
			stackOut(x,y,f) = tempRS[t](x,y);
		}
	}
}

template <typename T>
BufferedImage<T> NewStackHelper::inverseFourierTransformStack(
		const RawImage<tComplex<T>>& stack, 
		bool center,
		int num_threads)
{
	BufferedImage<T> out(2*(stack.xdim - 1), stack.ydim, stack.zdim);
	inverseFourierTransformStack(stack, out, center, num_threads);
	return out;
}

template <typename T>
void NewStackHelper::shiftStack(
		const RawImage<tComplex<T>>& stack,
		const std::vector<gravis::d2Vector> &pos,
		RawImage<tComplex<T>>& stackOut,
		bool centered,
		int num_threads)
{
	const int wh = stack.xdim;
	const int w = (wh-1)*2;
	const int h = stack.ydim;
	const int fc = stack.zdim;
	
	#pragma omp parallel for num_threads(num_threads)	
	for (int f = 0; f < fc; f++)
	{
		const double tx = centered? pos[f].x - w/2 : pos[f].x;
		const double ty = centered? pos[f].y - h/2 : pos[f].y;
		
		const double mx = 2 * PI * tx / (double) w;
		const double my = 2 * PI * ty / (double) h;
		
		for (int y = 0; y < h; y++)
		for (int x = 0; x < wh; x++)
		{
			const double xx = x;
			const double yy = y < h/2? y : y - h;
			
			const double phi = mx * xx + my * yy;
			
			tComplex<T> z(cos(phi), sin(phi));
			
			stackOut(x,y,f) = z * stack(x,y,f);
		}
	}
}

#endif

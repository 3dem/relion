#ifndef TOMO_EXTRACTION_HELPER_H
#define TOMO_EXTRACTION_HELPER_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <vector>
#include <string>
#include <omp.h>

#include <src/jaz/math/fft.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/image/interpolation.h>
#include <src/jaz/image/resampling.h>

#define EDGE_FALLOFF 5

class TomoExtraction
{
	public:
		
		template <typename T>
		static void extractFrameAt3D_Fourier(
				const RawImage<T>& stack, int f, int s, double bin,
				const gravis::d4Matrix& projIn,
				gravis::d3Vector center,
				RawImage<tComplex<T>>& out,
				gravis::d4Matrix& projOut,
				int num_threads = 1, 
				bool no_subpixel_shift = false,
				bool circle_crop = true);
		
		template <typename T>
		static void extractAt3D_Fourier(
				const RawImage<T>& stack, int s, double bin,
				const std::vector<gravis::d4Matrix>& projIn,
				const std::vector<gravis::d3Vector>& trajectory,
				RawImage<tComplex<T>>& out,
				std::vector<gravis::d4Matrix>& projOut,
				int num_threads = 1, 
				bool no_subpixel_shift = false,
				bool circle_crop = true);
		
		template <typename T>
		static void extractAt2D_Fourier(
				const RawImage<T>& stack, int s, double bin,
				const std::vector<gravis::d4Matrix>& projIn,
				const std::vector<gravis::d2Vector>& centers,
				RawImage<tComplex<T>>& out,
				std::vector<gravis::d4Matrix>& projOut,
				int num_threads = 1, 
				bool no_subpixel_shift = false,
				bool circle_crop = true);
		
		template <typename T>
		static void extractAt3D_real(
				const RawImage<T>& stack, int s, double bin,
				const std::vector<gravis::d4Matrix>& projIn,
				gravis::d3Vector center,
				RawImage<T>& out,
				std::vector<gravis::d4Matrix>& projOut,
				int num_threads = 1, 
				bool no_subpixel_shift = false,
				bool circle_crop = true);
		
		template <typename T>
		static void extractAt2D_real(
				const RawImage<T>& stack, int s, double bin,
				const std::vector<gravis::d4Matrix>& projIn,
				const std::vector<gravis::d2Vector>& centers,
				RawImage<T>& out,
				std::vector<gravis::d4Matrix>& projOut,
				int num_threads = 1, 
				bool no_subpixel_shift = false,
				bool circle_crop = true);
		
		template <typename T>
		static void extractSquares(
				const RawImage<T>& stack, 
				int w, int h, 
				const std::vector<gravis::d2Vector>& origins,
				RawImage<T>& out,
				bool center,
				int num_threads = 1);

		template <typename T>
		static void cropCircle(
				RawImage<T>& stack,
				double boundary,
				double falloff,
				int num_threads = 1);

		template <typename T>
		static void griddingPreCorrect(
				RawImage<T>& stack,
				double boundary,
				int num_threads = 1);
		
		/*template <typename T>
		static void extractSquares(
				const RawImage<T>& stack, 
				int w2D, int h2D, 
				const std::vector<gravis::d4Matrix>& projIn,
				int w3D, int h3D, int d3D,
				RawImage<T>& out,
				std::vector<gravis::d4Matrix>& projOut,
				bool center,
				int num_threads = 1);*/
		
				
};

template <typename T>
void TomoExtraction::extractFrameAt3D_Fourier(
		const RawImage<T>& stack, int f, int s, double bin,
		const gravis::d4Matrix& projIn,
		gravis::d3Vector center,
		RawImage<tComplex<T>>& out,
		gravis::d4Matrix& projOut,
		int num_threads, 
		bool no_subpixel_shift,
		bool circle_crop)
{
	std::vector<gravis::d4Matrix> projVec;
	
	extractAt3D_Fourier(
		stack.getConstSliceRef(f), s, bin, {projIn}, {center}, out, projVec, 
		num_threads, no_subpixel_shift, circle_crop);
	
	projOut = projVec[0];
}

template <typename T>
void TomoExtraction::extractAt3D_Fourier(
        const RawImage<T>& stack, int s, double bin,
		const std::vector<gravis::d4Matrix>& projIn,
		const std::vector<gravis::d3Vector>& trajectory,
		RawImage<tComplex<T>>& out,
		std::vector<gravis::d4Matrix>& projOut,
		int num_threads, 
		bool no_subpixel_shift,
		bool circle_crop)
{
	const int fc = projIn.size();
	std::vector<gravis::d2Vector> centers(fc);
	
	for (int f = 0; f < fc; f++)
	{
		const gravis::d4Vector q = projIn[f] * gravis::d4Vector(trajectory[f]);
		centers[f] = gravis::d2Vector(q.x, q.y);
	}
	
	extractAt2D_Fourier(stack, s, bin, projIn, centers, out, projOut, num_threads,
						no_subpixel_shift, circle_crop);
}

template <typename T>
void TomoExtraction::extractAt2D_Fourier(
        const RawImage<T>& stack, int s, double bin,
		const std::vector<gravis::d4Matrix>& projIn,
		const std::vector<gravis::d2Vector>& centers,
		RawImage<tComplex<T>>& out,
		std::vector<gravis::d4Matrix>& projOut,
		int num_threads, 
		bool no_subpixel_shift,
		bool circle_crop)
{
	const int sh = s/2 + 1;
	const int fc = stack.zdim;

    const int sb = (int)(s / bin + 0.5);
    const int sbh = sb/2 + 1;
	
    BufferedImage<T> smallStack(s,s,fc);
	projOut.resize(fc);
			
	std::vector<gravis::d2Vector> integralShift(fc);
			
	for (int f = 0; f < fc; f++)
	{
		integralShift[f] = gravis::d2Vector(
                round(centers[f].x) - s/2,
                round(centers[f].y) - s/2);
	}
	
	extractSquares(stack, s, s, integralShift, smallStack, false, num_threads);
	
	if (circle_crop) 
	{
		cropCircle(smallStack, 0, EDGE_FALLOFF, num_threads);
	}
	
	std::vector<gravis::d2Vector> posInNewImg(fc);
	
	for (int f = 0; f < fc; f++)
	{
		projOut[f] = projIn[f];
		
		if (no_subpixel_shift)
		{
			projOut[f](0,3) += sb/2 - round(centers[f].x);
			projOut[f](1,3) += sb/2 - round(centers[f].y);
		}
		else
		{
			projOut[f](0,3) += sb/2 - centers[f].x;
			projOut[f](1,3) += sb/2 - centers[f].y;
		}
		
        posInNewImg[f] = (centers[f] - integralShift[f]) / bin;
	}	
	
	BufferedImage<tComplex<T>> smallStackFS(sh,s,fc);
	
	if (no_subpixel_shift && bin == 1.0)
	{
		NewStackHelper::FourierTransformStack(smallStack, smallStackFS, true, num_threads);
		out.copyFrom(smallStackFS);
	}
	else
	{	
		NewStackHelper::FourierTransformStack(smallStack, smallStackFS, true, num_threads);
	
		if (bin != 1.0)
		{
			smallStackFS = Resampling::FourierCrop_fftwHalfStack(
					smallStackFS, bin, num_threads);
		}

		if (no_subpixel_shift)
		{
			//NewStackHelper::inverseFourierTransformStack(smallStackFS, smallStack, true, num_threads);
			out.copyFrom(smallStackFS);
		}
		else
		{
			NewStackHelper::shiftStack(smallStackFS, posInNewImg, out, true, num_threads);
		}
	}
}

template <typename T>
void TomoExtraction::extractAt3D_real(
        const RawImage<T>& stack, int s, double bin,
		const std::vector<gravis::d4Matrix>& projIn,
		gravis::d3Vector center,
		RawImage<T>& out,
		std::vector<gravis::d4Matrix>& projOut,
		int num_threads, 
		bool no_subpixel_shift,
		bool circle_crop)
{
	const int fc = projIn.size();
	std::vector<gravis::d2Vector> centers(fc);
	
	const gravis::d4Vector p(center.x, center.y, center.z, 1.0);
	
	for (int f = 0; f < fc; f++)
	{
		const gravis::d4Vector q = projIn[f] * p;
		centers[f] = gravis::d2Vector(q.x, q.y);
	}
	
    extractAt2D_real(stack, s, bin, projIn, centers, out, projOut, num_threads, 
					 no_subpixel_shift, circle_crop);
}

template <typename T>
void TomoExtraction::extractAt2D_real(
        const RawImage<T>& stack, int s, double bin,
		const std::vector<gravis::d4Matrix>& projIn,
		const std::vector<gravis::d2Vector>& centers,
		RawImage<T>& out,
		std::vector<gravis::d4Matrix>& projOut,
		int num_threads, 
		bool no_subpixel_shift,
		bool circle_crop)
{
	const int sh = s/2 + 1;
	const int fc = stack.zdim;

    const int sb = (int)(s / bin + 0.5);
    const int sbh = sb/2 + 1;
	
    BufferedImage<T> smallStack(s,s,fc);
	projOut.resize(fc);
			
	std::vector<gravis::d2Vector> integralShift(fc);
			
	for (int f = 0; f < fc; f++)
	{
		integralShift[f] = gravis::d2Vector(
                round(centers[f].x) - s/2,
                round(centers[f].y) - s/2);
	}
	
	extractSquares(stack, s, s, integralShift, smallStack, false, num_threads);
	
	if (circle_crop) 
	{
		cropCircle(smallStack, 0, EDGE_FALLOFF, num_threads);
	}
	
	std::vector<gravis::d2Vector> posInNewImg(fc);
	
	for (int f = 0; f < fc; f++)
	{
		projOut[f] = projIn[f];
		
		if (no_subpixel_shift)
		{
			projOut[f](0,3) += sb/2 - round(centers[f].x);
			projOut[f](1,3) += sb/2 - round(centers[f].y);
		}
		else
		{
			projOut[f](0,3) += sb/2 - centers[f].x;
			projOut[f](1,3) += sb/2 - centers[f].y;
		}
		
        posInNewImg[f] = (centers[f] - integralShift[f]) / bin;
	}	
	
	if (no_subpixel_shift && bin == 1.0)
	{
		out.copyFrom(smallStack);
	}
	else
	{	
		BufferedImage<tComplex<T>> smallStackFS(sh,s,fc);
		NewStackHelper::FourierTransformStack(smallStack, smallStackFS, true, num_threads);
	
		if (bin != 1.0)
		{
			smallStackFS = Resampling::FourierCrop_fftwHalfStack(
					smallStackFS, bin, num_threads);
			
			smallStack = BufferedImage<T>(sb,sb,fc);
		}

		if (no_subpixel_shift)
		{
			NewStackHelper::inverseFourierTransformStack(smallStackFS, smallStack, true, num_threads);
			out.copyFrom(smallStack);
		}
		else
		{
			NewStackHelper::shiftStack(smallStackFS, posInNewImg, smallStackFS, true, num_threads);
			NewStackHelper::inverseFourierTransformStack(smallStackFS, smallStack, true, num_threads);
			out.copyFrom(smallStack);
		}
	}
}

template <typename T>
void TomoExtraction::extractSquares(
		const RawImage<T>& stack, 
		int w, int h, 
		const std::vector<gravis::d2Vector>& origins,
		RawImage<T>& out,
		bool center,
		int num_threads)
{
	const int w0 = stack.xdim;
	const int h0 = stack.ydim;
	const int fc = stack.zdim;
		
	#pragma omp parallel for num_threads(num_threads)		
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			int xx = (center? (x + w/2) % w : x) + origins[f].x;
			int yy = (center? (y + h/2) % h : y) + origins[f].y;
			
			if (xx < 0) xx = 0;
			else if (xx >= w0) xx = w0 - 1;
			
			if (yy < 0) yy = 0;
			else if (yy >= h0) yy = h0 - 1;
			
			out(x,y,f) = stack(xx,yy,f);
		}
	}
}

template <typename T>
void TomoExtraction::cropCircle(
		RawImage<T>& stack,
		double boundary,
		double falloff,
		int num_threads)
{
	const int  w = stack.xdim;
	const int  h = stack.ydim;
	const int fc = stack.zdim;
	
	const double mx = w / 2;
	const double my = h / 2;

	const double crop_rad = mx - boundary;
		
	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		double meanInside = 0.0;
		double sumInside = 0.0;
		
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			double xx = x - mx;
			double yy = y - my;
			double  r = sqrt(xx*xx + yy*yy);
			
			double c;
			
			if (r < crop_rad - falloff) c = 1.0;
			else if (r < crop_rad) c = 0.5 - 0.5 * cos(PI * (crop_rad - r) / falloff);
			else c = 0.0;
			
			meanInside += c * stack(x,y,f);
			sumInside  += c;
		}
		
		if (sumInside > 0.0) meanInside /= sumInside;
		
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			double xx = x - mx;
			double yy = y - my;
			double  r = sqrt(xx*xx + yy*yy);
			
			double c;
			
			if (r < crop_rad - falloff) c = 1.0;
			else if (r < crop_rad) c = 0.5 - 0.5 * cos(PI * (crop_rad - r) / falloff);
			else c = 0.0;
						
			stack(x,y,f) = c * (stack(x,y,f) - meanInside);
		}
	}
}

template <typename T>
void TomoExtraction::griddingPreCorrect(
		RawImage<T>& stack,
		double boundary,
		int num_threads)
{
	const int  w = stack.xdim;
	const int  h = stack.ydim;
	const int fc = stack.zdim;

	const double mx = w / 2;
	const double my = h / 2;

	const double crop_rad = mx - boundary;
	const double eps = 1.0;

	#pragma omp parallel for num_threads(num_threads)
	for (int f = 0; f < fc; f++)
	{
		for (int y = 0; y < h; y++)
		for (int x = 0; x < w; x++)
		{
			double xx = x - mx;
			double yy = y - my;
			double  r = sqrt(xx*xx + yy*yy);

			double c;

			if (r < eps) c = 1.0;
			else if (r < crop_rad - eps) c = r / sin(PI * r / crop_rad);
			else c = 0.0;

			stack(x,y,f) *= c;
		}
	}
}


#endif

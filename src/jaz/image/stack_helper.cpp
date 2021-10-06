#include "stack_helper.h"


void NewStackHelper::FourierTransformStack_fast(
		const RawImage<float>& stack,
		RawImage<fComplex>& stackOut,
		bool center,
		int num_threads)
{
	const int w = stack.xdim;
	const int h = stack.ydim;
	const int fc = stack.zdim;

	const int wh = w/2 + 1;


	std::vector<BufferedImage<float>> tempRS(num_threads);
	std::vector<BufferedImage<fComplex>> tempFS(num_threads);
	std::vector<FFT::FloatPlan> plans(num_threads);

	for (int t = 0; t < num_threads; t++)
	{
		tempRS[t] = BufferedImage<float>(w,h);
		tempFS[t] = BufferedImage<fComplex>(wh,h);
		plans[t] = FFT::FloatPlan(tempRS[t], tempFS[t], FFTW_ESTIMATE);
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

		FFT::FourierTransform(tempRS[t], tempFS[t], plans[t], FFT::Both);

		for (long int y = 0; y < h; y++)
		for (long int x = 0; x < wh; x++)
		{
			stackOut(x,y,f) = tempFS[t](x,y);
		}
	}
}


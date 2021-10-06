#include "fft_buffer.h"

FloatFFTBuffer::FloatFFTBuffer()
{
}

FloatFFTBuffer::FloatFFTBuffer(int w, int h, int d)
:	real(w,h,d),
	complex(w/2 + 1, h, d)
{
	plan = FFT::FloatPlan(real, complex, FFTW_ESTIMATE);
}

void FloatFFTBuffer::FourierTransform()
{
	FFT::FourierTransform(real, complex, plan, FFT::Both);
}

void FloatFFTBuffer::inverseFourierTransform()
{
	FFT::inverseFourierTransform(complex, real, plan, FFT::Both);
}


FloatFFTBufferArray::FloatFFTBufferArray(int w, int h, int d, int number)
:	buffers(number)
{
	for (int i = 0; i < number; i++)
	{
		buffers[i] = FloatFFTBuffer(w,h,d);
	}
}

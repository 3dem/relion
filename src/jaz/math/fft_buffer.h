#ifndef FFT_BUFFER_H
#define FFT_BUFFER_H

#include "fft.h"

class FloatFFTBuffer
{
	public:

		FloatFFTBuffer();
		FloatFFTBuffer(int w, int h, int d);

			BufferedImage<float> real;
			BufferedImage<fComplex> complex;
			FFT::FloatPlan plan;

		void FourierTransform();
		void inverseFourierTransform();
};

class FloatFFTBufferArray
{
	public:

		FloatFFTBufferArray(int w, int h, int d, int number);

			std::vector<FloatFFTBuffer> buffers;
};

#endif

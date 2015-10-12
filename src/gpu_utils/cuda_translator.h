#ifndef CUDA_TRANSLATOR_H_
#define CUDA_TRANSLATOR_H_

#include "src/multidim_array.h"
#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <cuda_runtime.h>

class CudaTranslator
{
	unsigned img_size;
	XFLOAT *d_shifter_real, *d_shifter_imag;

public:
	class Plan
	{
		friend class CudaTranslator;
		unsigned long itrans_min,
			itrans_max;
		XFLOAT scale_correction;

		cudaStream_t stream;

		CudaGlobalPtr<XFLOAT> ctf,
			imgs_real,
			imgs_imag;

		bool do_ctf, do_scale;

	public:
		Plan(
			Complex *h_Fimgs,
			unsigned image_size,
			unsigned long itrans_min,
			unsigned long itrans_max,
			CudaCustomAllocator *allocator,
			cudaStream_t stream,
			RFLOAT scale_correction = 1,
			RFLOAT *h_ctf = NULL):
				itrans_min(itrans_min),
				itrans_max(itrans_max),
				scale_correction(scale_correction),
				stream(stream),
				ctf(stream, allocator),
				imgs_real(image_size, stream, allocator),
				imgs_imag(image_size, stream, allocator),
				do_ctf(h_ctf != NULL),
				do_scale(scale_correction != 1)
		{
			for (int i = 0; i < image_size; i ++) //TODO Make better copy than this!
			{
				imgs_real[i] = (XFLOAT) h_Fimgs[i].real;
				imgs_imag[i] = (XFLOAT) h_Fimgs[i].imag;
			}

			imgs_real.put_on_device();
			imgs_imag.put_on_device();

			if (do_ctf)
			{
				ctf.size = image_size;
				ctf.host_alloc();

				for (int i = 0; i < image_size; i ++) //TODO Make better copy than this!
					ctf[i] = (XFLOAT) h_ctf[i];

				ctf.put_on_device();
			}
		};
	};

	CudaTranslator():
		img_size(0), d_shifter_real(0), d_shifter_imag(0)
	{};

	CudaTranslator(std::vector<MultidimArray<Complex> >  &shifter):
		img_size(0), d_shifter_real(0), d_shifter_imag(0)
	{ setShifters(shifter); };

	void setShifters(XFLOAT *h_real, XFLOAT *h_imag, unsigned count, unsigned img_size);
	void setShifters(std::vector<MultidimArray<Complex> >  &shifter);

	void translate(Plan &plan, XFLOAT *d_imgs_shifted_real, XFLOAT *d_imgs_shifted_imag);

	void clear();
	~CudaTranslator();

};

#endif

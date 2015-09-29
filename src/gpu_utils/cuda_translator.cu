#include "src/gpu_utils/cuda_translator.h"

#define TRANS_BLOCK_SIZE 128

void CudaTranslator::setShifters(XFLOAT *h_real, XFLOAT *h_imag, unsigned count, unsigned img_sz)
{
	img_size = img_sz;

	if (d_shifter_real != NULL)
		DEBUG_HANDLE_ERROR(cudaFree(d_shifter_real));

	HANDLE_ERROR(cudaMalloc( (void**) &d_shifter_real, img_sz * count * sizeof(XFLOAT)));
	DEBUG_HANDLE_ERROR(cudaMemcpy( d_shifter_real, h_real, img_sz * count * sizeof(XFLOAT), cudaMemcpyHostToDevice));

	if (d_shifter_imag != NULL)
		DEBUG_HANDLE_ERROR(cudaFree(d_shifter_imag));

	HANDLE_ERROR(cudaMalloc( (void**) &d_shifter_imag, img_sz * count * sizeof(XFLOAT)));
	DEBUG_HANDLE_ERROR(cudaMemcpy( d_shifter_imag, h_imag, img_sz * count * sizeof(XFLOAT), cudaMemcpyHostToDevice));
}

void CudaTranslator::setShifters(std::vector<MultidimArray<Complex> >  &shifter)
{
	size_t img_sz = shifter[0].nzyxdim;

	XFLOAT *real = new XFLOAT[shifter.size() * img_sz];
	XFLOAT *imag = new XFLOAT[shifter.size() * img_sz];

	for (int i = 0; i < shifter.size(); i ++) //TODO Make better copy than this!
	{
		for (int j = 0; j < img_sz; j ++)
		{
			real[i * img_sz + j] = (XFLOAT) shifter[i].data[j].real;
			imag[i * img_sz + j] = (XFLOAT) shifter[i].data[j].imag;
		}
	}

	setShifters(real, imag, shifter.size(), img_sz);

	delete [] real;
	delete [] imag;
}

template<bool do_ctf_correction, bool do_scale_correction>
__global__ void cuda_translation_kernel(
		XFLOAT *g_ctfs,
		XFLOAT *g_Fimgs_real,
		XFLOAT *g_Fimgs_imag,
		XFLOAT *g_Fimgs_shifted_real,
		XFLOAT *g_Fimgs_shifted_imag,
		XFLOAT *g_fftshifts_real,
		XFLOAT *g_fftshifts_imag,
		unsigned long image_size,
		unsigned long image_idx_offset,
		double scale_correction)
{
	int img_idx = blockIdx.x + image_idx_offset;
	int tid = threadIdx.x;

	XFLOAT real, imag;

	unsigned pixel_pass_num( ceilf( (float)image_size / (float)TRANS_BLOCK_SIZE ) );

	for (unsigned pass = 0; pass < pixel_pass_num; pass++)
	{
		unsigned local_pixel = (pass * TRANS_BLOCK_SIZE) + tid;

		if(local_pixel >= image_size)
			break;

		unsigned global_pixel = img_idx * image_size + local_pixel;

		real = g_fftshifts_real[global_pixel] * g_Fimgs_real[local_pixel]
		     - g_fftshifts_imag[global_pixel] * g_Fimgs_imag[local_pixel];

		imag = g_fftshifts_real[global_pixel] * g_Fimgs_imag[local_pixel]
		     + g_fftshifts_imag[global_pixel] * g_Fimgs_real[local_pixel];

		if (do_scale_correction)
		{
			real *= scale_correction;
			imag *= scale_correction;
		}

		if (do_ctf_correction)
		{
			real /= g_ctfs[local_pixel];
			imag /= g_ctfs[local_pixel];
		}

		g_Fimgs_shifted_real[global_pixel] = real;
		g_Fimgs_shifted_imag[global_pixel] = imag;
	}
}

void CudaTranslator::translate(Plan &plan, XFLOAT *d_imgs_shifted_real, XFLOAT *d_imgs_shifted_imag)
{
#ifdef DEBUG_CUDA
	if (d_shifter_real == NULL || d_shifter_imag == NULL)
	{
		printf("DEBUG_ERROR: Shifters not set in call to CudaTranslator::translate.");
		raise(SIGSEGV);
	}
#endif

	unsigned img_count = plan.itrans_max - plan.itrans_min;

	if (plan.do_ctf && plan.do_scale)
			cuda_translation_kernel<true,true><<<img_count,TRANS_BLOCK_SIZE,0,plan.stream>>>(
					~plan.ctf,
					~plan.imgs_real,
					~plan.imgs_imag,
					d_imgs_shifted_real,
					d_imgs_shifted_imag,
					d_shifter_real,
					d_shifter_imag,
					img_size,
					plan.itrans_min,
					1/plan.scale_correction
					);
		else if (plan.do_ctf)
			cuda_translation_kernel<true,false><<<img_count,TRANS_BLOCK_SIZE,0,plan.stream>>>(
					~plan.ctf,
					~plan.imgs_real,
					~plan.imgs_imag,
					d_imgs_shifted_real,
					d_imgs_shifted_imag,
					d_shifter_real,
					d_shifter_imag,
					img_size,
					plan.itrans_min,
					1
					);
		else if (plan.do_scale)
			cuda_translation_kernel<false,true><<<img_count,TRANS_BLOCK_SIZE,0,plan.stream>>>(
					NULL,
					~plan.imgs_real,
					~plan.imgs_imag,
					d_imgs_shifted_real,
					d_imgs_shifted_imag,
					d_shifter_real,
					d_shifter_imag,
					img_size,
					plan.itrans_min,
					1/plan.scale_correction
					);
		else
			cuda_translation_kernel<false,false><<<img_count,TRANS_BLOCK_SIZE,0,plan.stream>>>(
					NULL,
					~plan.imgs_real,
					~plan.imgs_imag,
					d_imgs_shifted_real,
					d_imgs_shifted_imag,
					d_shifter_real,
					d_shifter_imag,
					img_size,
					plan.itrans_min,
					1
					);
}

void CudaTranslator::clear()
{
	if (d_shifter_real != NULL)
		DEBUG_HANDLE_ERROR(cudaFree(d_shifter_real));
	if (d_shifter_imag != NULL)
		DEBUG_HANDLE_ERROR(cudaFree(d_shifter_imag));

	d_shifter_real = d_shifter_imag = NULL;
}

CudaTranslator::~CudaTranslator()
{
	clear();
}

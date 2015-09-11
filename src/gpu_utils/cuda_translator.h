#ifndef CUDA_TRANSLATOR_H_
#define CUDA_TRANSLATOR_H_

#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <cuda_runtime.h>

class CudaTranslator
{
	friend class CudaProjectorKernel;

	unsigned coarse_size, fine_size;

	XFLOAT *ab1_coarse, *ab1_current, *ab2_coarse, *ab2_current;

public:
	CudaTranslator():
		coarse_size(0), fine_size(0),
		ab1_coarse(0), ab1_current(0),
		ab2_coarse(0), ab2_current(0)
	{};

	void setFftShiftsCoarse(XFLOAT ab, XFLOAT ab2, unsigned image_size);
	void setFftShiftsFine(XFLOAT ab, XFLOAT ab2, unsigned image_size);

	~CudaTranslator();

};

#endif

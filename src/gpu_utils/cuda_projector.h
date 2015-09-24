#ifndef CUDA_PROJECTOR_H_
#define CUDA_PROJECTOR_H_

#include "src/complex.h"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_mem_utils.h"
#include <cuda_runtime.h>

class CudaProjector
{
	friend class CudaProjectorKernel;

	int mdlX, mdlY, mdlZ, mdlXYZ, mdlMaxR,
	    mdlInitY, mdlInitZ,
	    padding_factor;

#ifndef CUDA_DOUBLE_PRECISION
	float *texArrayReal2D, *texArrayImag2D;
	cudaArray_t *texArrayReal, *texArrayImag;
	cudaTextureObject_t *mdlReal, *mdlImag;
	size_t pitch2D;
#else
	double *mdlReal, *mdlImag;
#endif

public:
	CudaProjector():
			mdlX(0), mdlY(0), mdlZ(0),
			mdlXYZ(0), mdlMaxR(0),
			mdlInitY(0), mdlInitZ(0),
			padding_factor(0)
	{
#ifndef CUDA_DOUBLE_PRECISION
		texArrayReal2D = 0;
		texArrayImag2D = 0;
		texArrayReal = 0;
		texArrayImag = 0;
		mdlReal = 0;
		mdlImag = 0;
		pitch2D = 0;
#else
		mdlReal = 0;
		mdlImag = 0;
#endif
	}

	void setMdlDim(
			int xdim, int ydim, int zdim,
			int inity, int initz,
			int maxr, int paddingFactor);

	void setMdlData(XFLOAT *real, XFLOAT *imag);
	void setMdlData(Complex *data);

	void clear();

	~CudaProjector()
	{
		clear();
	};

};

#endif
